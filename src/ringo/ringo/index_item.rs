use bit_vec::BitVec;

pub struct IndexItem {
    pub position: usize,
    pub fingerprint: Vec<u8>
}

impl IndexItem {
    pub fn new(position: usize, fingerprint: BitVec) -> IndexItem {
        IndexItem {
            position,
            fingerprint: fingerprint.to_bytes()
        }
    }
}
