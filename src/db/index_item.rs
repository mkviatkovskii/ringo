use crate::model::fingerprint::Fingerprint;
use std::io;

/// On-disk record: `u64` little-endian line offset, then `bits / 8` fingerprint bytes.
pub fn index_item_size(bits: usize) -> usize {
    8 + bits / 8
}

#[derive(Debug)]
pub struct IndexItem {
    pub position: usize,
    pub fingerprint: Fingerprint,
}

impl IndexItem {
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = (self.position as u64).to_le_bytes().to_vec();
        bytes.extend(self.fingerprint.to_bytes());
        bytes
    }

    pub fn from_bytes(bytes: &[u8]) -> io::Result<IndexItem> {
        let invalid = |msg: &str| io::Error::new(io::ErrorKind::InvalidData, msg);
        let (position, fingerprint) = bytes
            .split_first_chunk::<8>()
            .ok_or_else(|| invalid("index record is too short"))?;
        Ok(IndexItem {
            position: usize::try_from(u64::from_le_bytes(*position))
                .map_err(|_| invalid("offset does not fit in usize"))?,
            fingerprint: Fingerprint::from_bytes(fingerprint),
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::db::index_item::{index_item_size, IndexItem};
    use crate::model::fingerprint::Fingerprint;
    use fixedbitset::FixedBitSet;

    #[test]
    fn test_index_item_encode_decode() {
        let fp = Fingerprint(FixedBitSet::with_capacity(512));
        let mut ii = IndexItem {
            position: 0,
            fingerprint: fp,
        };
        ii.position = 0;
        ii.fingerprint.0.set(1, true);
        ii.fingerprint.0.set(17, true);

        let decoded = IndexItem::from_bytes(&ii.to_bytes()).unwrap();
        assert_eq!(decoded.position, 0);
        assert_eq!(
            decoded.fingerprint.0.ones().collect::<Vec<usize>>(),
            vec![1, 17]
        );
    }

    #[test]
    fn test_index_item_dense_large_offset() {
        // previously overflowed the fixed 72-byte slot with bincode varint encoding
        let mut fp = FixedBitSet::with_capacity(512);
        fp.insert_range(..);
        let ii = IndexItem {
            position: u32::MAX as usize,
            fingerprint: Fingerprint(fp),
        };
        let bytes = ii.to_bytes();
        assert_eq!(bytes.len(), index_item_size(512));
        assert_eq!(index_item_size(512), 72);

        let decoded = IndexItem::from_bytes(&bytes).unwrap();
        assert_eq!(decoded.position, u32::MAX as usize);
        assert_eq!(decoded.fingerprint.0.count_ones(..), 512);
    }
}
