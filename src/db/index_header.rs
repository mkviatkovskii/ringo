const MAGIC: [u8; 8] = *b"RINGOFP\0";
const VERSION: u16 = 2;

/// A multiple of 8, so the records after it stay 8-byte aligned.
pub const INDEX_HEADER_SIZE: usize = 32;

/// Start of every `.fp` file: magic (8 bytes), format version (`u16`), ECFP radius
/// (`u16`), fingerprint bits (`u32`), SMILES file length (`u64`), record count
/// (`u64`), all little-endian. Search builds the query fingerprint from radius and
/// bits, and refuses an index whose SMILES file or record count no longer matches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexHeader {
    pub radius: usize,
    pub bits: usize,
    pub smiles_len: u64,
    pub count: u64,
}

impl IndexHeader {
    pub fn new(radius: usize, bits: usize) -> Result<IndexHeader, String> {
        if radius > u16::MAX as usize {
            return Err(format!("radius {radius} is too large"));
        }
        // whole 64-bit blocks keep the byte layout the same on 32- and 64-bit targets
        if bits == 0 || !bits.is_multiple_of(64) || bits > u32::MAX as usize {
            return Err(format!(
                "fingerprint length {bits} must be a positive multiple of 64"
            ));
        }
        Ok(IndexHeader {
            radius,
            bits,
            smiles_len: 0,
            count: 0,
        })
    }

    #[expect(
        clippy::cast_possible_truncation,
        reason = "new() keeps radius within u16 and bits within u32"
    )]
    pub fn to_bytes(&self) -> [u8; INDEX_HEADER_SIZE] {
        let mut bytes = [0u8; INDEX_HEADER_SIZE];
        bytes[..8].copy_from_slice(&MAGIC);
        bytes[8..10].copy_from_slice(&VERSION.to_le_bytes());
        bytes[10..12].copy_from_slice(&(self.radius as u16).to_le_bytes());
        bytes[12..16].copy_from_slice(&(self.bits as u32).to_le_bytes());
        bytes[16..24].copy_from_slice(&self.smiles_len.to_le_bytes());
        bytes[24..32].copy_from_slice(&self.count.to_le_bytes());
        bytes
    }

    pub fn from_bytes(bytes: &[u8; INDEX_HEADER_SIZE]) -> Result<IndexHeader, String> {
        if bytes[..8] != MAGIC {
            return Err("not a ringo index, or written by an older version: re-index".into());
        }
        let version = u16::from_le_bytes([bytes[8], bytes[9]]);
        if version != VERSION {
            return Err(format!(
                "index format version {version} is not supported, expected {VERSION}: re-index"
            ));
        }
        let radius = u16::from_le_bytes([bytes[10], bytes[11]]);
        let bits = u32::from_le_bytes([bytes[12], bytes[13], bytes[14], bytes[15]]);
        let u64_at = |i: usize| {
            let mut word = [0u8; 8];
            word.copy_from_slice(&bytes[i..i + 8]);
            u64::from_le_bytes(word)
        };
        Ok(IndexHeader {
            smiles_len: u64_at(16),
            count: u64_at(24),
            ..IndexHeader::new(radius as usize, bits as usize)?
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::db::index_header::{IndexHeader, INDEX_HEADER_SIZE};

    #[test]
    fn test_header_round_trip() {
        let header = IndexHeader {
            smiles_len: u64::MAX - 1,
            count: 12345,
            ..IndexHeader::new(3, 1024).unwrap()
        };
        assert_eq!(IndexHeader::from_bytes(&header.to_bytes()), Ok(header));
    }

    #[test]
    fn test_header_rejects_bad_input() {
        assert!(IndexHeader::new(2, 0).is_err());
        assert!(IndexHeader::new(2, 100).is_err());
        assert!(IndexHeader::new(70_000, 512).is_err());
        // a pre-header index starts with a line offset, not the magic
        assert!(IndexHeader::from_bytes(&[0u8; INDEX_HEADER_SIZE]).is_err());

        let mut bytes = IndexHeader::new(2, 512).unwrap().to_bytes();
        bytes[8] = 1;
        assert!(IndexHeader::from_bytes(&bytes)
            .unwrap_err()
            .contains("version 1"));
    }
}
