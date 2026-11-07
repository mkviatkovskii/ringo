use fixedbitset::{Block, FixedBitSet};

#[derive(Debug)]
pub struct Fingerprint(pub FixedBitSet);

impl Fingerprint {
    /// Bit `n` goes to byte `n / 8` (little-endian blocks), so the layout is the same
    /// on 32- and 64-bit targets. The length must be a multiple of 64 bits.
    pub fn to_bytes(&self) -> Vec<u8> {
        assert_eq!(
            self.0.len() % 64,
            0,
            "fingerprint length must be a multiple of 64 bits"
        );
        self.0
            .as_slice()
            .iter()
            .flat_map(|block| block.to_le_bytes())
            .collect()
    }

    /// The fingerprint gets `bytes.len() * 8` bits.
    pub fn from_bytes(bytes: &[u8]) -> Fingerprint {
        let (chunks, rest) = bytes.as_chunks::<{ size_of::<Block>() }>();
        assert!(rest.is_empty(), "fingerprint bytes must be whole blocks");
        let blocks = chunks.iter().map(|chunk| Block::from_le_bytes(*chunk));
        Fingerprint(FixedBitSet::with_capacity_and_blocks(
            bytes.len() * 8,
            blocks,
        ))
    }
}

#[cfg(test)]
mod tests {
    use crate::model::fingerprint::Fingerprint;
    use fixedbitset::FixedBitSet;

    #[test]
    fn test_fingerprint_encode_decode() {
        let mut fp = Fingerprint(FixedBitSet::with_capacity(512));
        fp.0.set(1, true);
        fp.0.set(17, true);

        let decoded = Fingerprint::from_bytes(&fp.to_bytes());
        assert_eq!(decoded.0.ones().collect::<Vec<usize>>(), vec![1, 17]);
    }

    #[test]
    fn test_fingerprint_byte_layout() {
        let mut fp = Fingerprint(FixedBitSet::with_capacity(512));
        fp.0.insert(9);
        fp.0.insert(511);
        let bytes = fp.to_bytes();
        assert_eq!(bytes.len(), 64);
        assert_eq!(bytes[1], 0b10);
        assert_eq!(bytes[63], 0b1000_0000);
        assert_eq!(bytes.iter().filter(|b| **b != 0).count(), 2);
    }

    #[test]
    fn test_fingerprint_dense_round_trip() {
        for bits in [64, 512, 1024] {
            let mut fp = Fingerprint(FixedBitSet::with_capacity(bits));
            fp.0.insert_range(..);
            let decoded = Fingerprint::from_bytes(&fp.to_bytes());
            assert_eq!(decoded.0.len(), bits);
            assert_eq!(decoded.0.count_ones(..), bits);
        }
    }
}
