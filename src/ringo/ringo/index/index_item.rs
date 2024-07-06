use crate::ringo::fingerprint::fingerprint::Fingerprint;
use bincode::{Decode, Encode};

#[derive(Debug, Encode, Decode)]
pub struct IndexItem {
    pub position: usize,
    pub fingerprint: Fingerprint,
}

#[cfg(test)]
mod tests {
    use crate::ringo::fingerprint::fingerprint::Fingerprint;
    use crate::ringo::ringo::index::index_item::IndexItem;
    use bincode::config::standard;
    use bincode::{decode_from_slice, encode_to_vec};
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

        let encoded = encode_to_vec(&ii, standard()).unwrap();
        let decoded: IndexItem = decode_from_slice(&encoded, standard()).unwrap().0;
        assert_eq!(
            decoded.fingerprint.0.ones().collect::<Vec<usize>>(),
            vec![1, 17]
        );
    }
}
