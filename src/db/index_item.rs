use crate::model::fingerprint::{Fingerprint, FINGERPRINT_BYTES};

/// On-disk record: `u64` little-endian line offset, then fingerprint bytes.
pub const INDEX_ITEM_SIZE: usize = 8 + FINGERPRINT_BYTES;

#[derive(Debug)]
pub struct IndexItem {
    pub position: usize,
    pub fingerprint: Fingerprint,
}

impl IndexItem {
    pub fn to_bytes(&self) -> [u8; INDEX_ITEM_SIZE] {
        let mut bytes = [0u8; INDEX_ITEM_SIZE];
        bytes[..8].copy_from_slice(&(self.position as u64).to_le_bytes());
        bytes[8..].copy_from_slice(&self.fingerprint.to_bytes());
        bytes
    }

    pub fn from_bytes(bytes: &[u8; INDEX_ITEM_SIZE]) -> IndexItem {
        let position = u64::from_le_bytes(bytes[..8].try_into().unwrap());
        IndexItem {
            position: usize::try_from(position).expect("offset does not fit in usize"),
            fingerprint: Fingerprint::from_bytes(bytes[8..].try_into().unwrap()),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::db::index_item::{IndexItem, INDEX_ITEM_SIZE};
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

        let decoded = IndexItem::from_bytes(&ii.to_bytes());
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
        assert_eq!(bytes.len(), INDEX_ITEM_SIZE);
        assert_eq!(INDEX_ITEM_SIZE, 72);

        let decoded = IndexItem::from_bytes(&bytes);
        assert_eq!(decoded.position, u32::MAX as usize);
        assert_eq!(decoded.fingerprint.0.count_ones(..), 512);
    }
}
