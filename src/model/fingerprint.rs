use fixedbitset::{Block, FixedBitSet};

pub const FINGERPRINT_SIZE: usize = 512;
pub const FINGERPRINT_BYTES: usize = FINGERPRINT_SIZE / 8;

#[derive(Debug)]
pub struct Fingerprint(pub FixedBitSet);

impl Fingerprint {
    /// Bit `n` goes to byte `n / 8` (little-endian blocks), so the layout is the same
    /// on 32- and 64-bit targets.
    pub fn to_bytes(&self) -> [u8; FINGERPRINT_BYTES] {
        assert_eq!(
            self.0.len(),
            FINGERPRINT_SIZE,
            "unexpected fingerprint size"
        );
        let mut bytes = [0u8; FINGERPRINT_BYTES];
        for (chunk, block) in bytes
            .chunks_exact_mut(size_of::<Block>())
            .zip(self.0.as_slice())
        {
            chunk.copy_from_slice(&block.to_le_bytes());
        }
        bytes
    }

    pub fn from_bytes(bytes: &[u8; FINGERPRINT_BYTES]) -> Fingerprint {
        let blocks = bytes
            .chunks_exact(size_of::<Block>())
            .map(|chunk| Block::from_le_bytes(chunk.try_into().unwrap()));
        Fingerprint(FixedBitSet::with_capacity_and_blocks(
            FINGERPRINT_SIZE,
            blocks,
        ))
    }
}

#[cfg(test)]
mod tests {
    use crate::model::fingerprint::{Fingerprint, FINGERPRINT_SIZE};
    use fixedbitset::FixedBitSet;

    #[test]
    fn test_fingerprint_encode_decode() {
        let mut fp = Fingerprint(FixedBitSet::with_capacity(FINGERPRINT_SIZE));
        fp.0.set(1, true);
        fp.0.set(17, true);

        let decoded = Fingerprint::from_bytes(&fp.to_bytes());
        assert_eq!(decoded.0.ones().collect::<Vec<usize>>(), vec![1, 17]);
    }

    #[test]
    fn test_fingerprint_byte_layout() {
        let mut fp = Fingerprint(FixedBitSet::with_capacity(FINGERPRINT_SIZE));
        fp.0.insert(9);
        fp.0.insert(511);
        let bytes = fp.to_bytes();
        assert_eq!(bytes[1], 0b10);
        assert_eq!(bytes[63], 0b1000_0000);
        assert_eq!(bytes.iter().filter(|b| **b != 0).count(), 2);
    }

    #[test]
    fn test_fingerprint_dense_round_trip() {
        let mut fp = Fingerprint(FixedBitSet::with_capacity(FINGERPRINT_SIZE));
        fp.0.insert_range(..);
        let decoded = Fingerprint::from_bytes(&fp.to_bytes());
        assert_eq!(decoded.0.count_ones(..), FINGERPRINT_SIZE);
    }
}
