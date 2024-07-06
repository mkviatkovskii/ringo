use bincode::de::BorrowDecoder;
use bincode::error::{DecodeError, EncodeError};
use fixedbitset::{Block, FixedBitSet};

pub const FINGERPRINT_SIZE: usize = 512;

#[derive(Debug)]
pub struct Fingerprint(pub FixedBitSet);

impl bincode::Encode for Fingerprint {
    fn encode<E: bincode::enc::Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        self.0.as_slice().encode(encoder)?;
        Ok(())
    }
}

impl bincode::Decode for Fingerprint {
    fn decode<D: bincode::de::Decoder>(decoder: &mut D) -> Result<Self, DecodeError> {
        let slice = Vec::<Block>::decode(decoder)?;
        let fp = FixedBitSet::with_capacity_and_blocks(FINGERPRINT_SIZE, slice);
        Ok(Fingerprint(fp))
    }
}

impl<'de> bincode::BorrowDecode<'de> for Fingerprint {
    fn borrow_decode<D: BorrowDecoder<'de>>(decoder: &mut D) -> Result<Self, DecodeError> {
        let slice = Vec::<Block>::borrow_decode(decoder)?;
        let fp = FixedBitSet::with_capacity_and_blocks(FINGERPRINT_SIZE, slice);
        Ok(Fingerprint(fp))
    }
}

#[cfg(test)]
mod tests {
    use crate::ringo::ringo::fingerprint::{Fingerprint, FINGERPRINT_SIZE};
    use fixedbitset::FixedBitSet;

    #[test]
    fn test_fingerprint_encode_decode() {
        let mut fp = Fingerprint(FixedBitSet::with_capacity(FINGERPRINT_SIZE));
        fp.0.set(1, true);
        fp.0.set(17, true);

        let mut buf = vec![0u8; FINGERPRINT_SIZE / 8];
        let encoded = bincode::encode_into_slice(&fp, buf.as_mut_slice(), bincode::config::standard()).unwrap();

        let decoded: Fingerprint =
            bincode::decode_from_slice(&buf, bincode::config::standard())
                .unwrap()
                .0;
        assert_eq!(decoded.0.ones().collect::<Vec<usize>>(), vec![1, 17]);
    }
}
