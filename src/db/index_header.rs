const MAGIC: [u8; 8] = *b"RINGOFP\0";
// Bump once per release when index compatibility changes.
const VERSION: u16 = 1;

/// A multiple of 8, so the group table and the records after it stay 8-byte aligned.
pub const INDEX_HEADER_SIZE: usize = 40;

/// Start of every `.fp` file: magic (8 bytes), format version (`u16`), ECFP radius
/// (`u16`), fingerprint bits (`u32`), SMILES file length (`u64`), record count
/// (`u64`), smallest and largest bit count of a record (`u32` each), all little-endian.
/// Search builds the query fingerprint from radius and bits, and refuses an index whose
/// SMILES file or record count no longer matches.
///
/// The group table follows: a `u64` record count for each bit count from `min_ones` to
/// `max_ones`. Then the records, by bit count and in line order within a bit count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexHeader {
    pub radius: usize,
    pub bits: usize,
    pub smiles_len: u64,
    pub count: u64,
    pub min_ones: usize,
    pub max_ones: usize,
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
            min_ones: 0,
            max_ones: 0,
        })
    }

    /// Entries in the group table
    pub fn groups(&self) -> usize {
        self.max_ones - self.min_ones + 1
    }

    /// Byte offset of the first record, after the header and the group table
    pub fn records_start(&self) -> usize {
        INDEX_HEADER_SIZE + 8 * self.groups()
    }

    #[expect(
        clippy::cast_possible_truncation,
        reason = "new() keeps radius within u16 and bits within u32, and bit counts are at most bits"
    )]
    pub fn to_bytes(&self) -> [u8; INDEX_HEADER_SIZE] {
        let mut bytes = [0u8; INDEX_HEADER_SIZE];
        bytes[..8].copy_from_slice(&MAGIC);
        bytes[8..10].copy_from_slice(&VERSION.to_le_bytes());
        bytes[10..12].copy_from_slice(&(self.radius as u16).to_le_bytes());
        bytes[12..16].copy_from_slice(&(self.bits as u32).to_le_bytes());
        bytes[16..24].copy_from_slice(&self.smiles_len.to_le_bytes());
        bytes[24..32].copy_from_slice(&self.count.to_le_bytes());
        bytes[32..36].copy_from_slice(&(self.min_ones as u32).to_le_bytes());
        bytes[36..40].copy_from_slice(&(self.max_ones as u32).to_le_bytes());
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
        let u32_at =
            |i: usize| u32::from_le_bytes([bytes[i], bytes[i + 1], bytes[i + 2], bytes[i + 3]]);
        let bits = u32_at(12);
        let u64_at = |i: usize| {
            let mut word = [0u8; 8];
            word.copy_from_slice(&bytes[i..i + 8]);
            u64::from_le_bytes(word)
        };
        let (min_ones, max_ones) = (u32_at(32), u32_at(36));
        if min_ones > max_ones || max_ones > bits {
            return Err(format!(
                "bit counts {min_ones} to {max_ones} do not fit {bits}-bit fingerprints: re-index"
            ));
        }
        Ok(IndexHeader {
            smiles_len: u64_at(16),
            count: u64_at(24),
            min_ones: min_ones as usize,
            max_ones: max_ones as usize,
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
            min_ones: 7,
            max_ones: 1024,
            ..IndexHeader::new(3, 1024).unwrap()
        };
        assert_eq!(&header.to_bytes()[8..10], &1u16.to_le_bytes());
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
        bytes[8..10].copy_from_slice(&2u16.to_le_bytes());
        assert!(IndexHeader::from_bytes(&bytes)
            .unwrap_err()
            .contains("version 2"));

        // bit counts out of order, or more than the fingerprint holds
        for (min_ones, max_ones) in [(5, 4), (0, 513)] {
            let header = IndexHeader {
                min_ones,
                max_ones,
                ..IndexHeader::new(2, 512).unwrap()
            };
            assert!(IndexHeader::from_bytes(&header.to_bytes())
                .unwrap_err()
                .contains("bit counts"));
        }
    }
}
