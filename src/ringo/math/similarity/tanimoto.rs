use bit_vec::BitVec;

pub fn tanimoto_vec(a: &BitVec, b: &BitVec) -> f32 {
    let mut and_ = a.clone();
    let mut or_ = a.clone();
    and_.and(b);
    or_.or(b);

    let mut dividend: u32 = 0;
    for b in and_.blocks() {
        dividend += b.count_ones();
    }
    let mut divisor: u32 = 0;
    for b in or_.blocks() {
        divisor += b.count_ones();
    }

    return dividend as f32 / divisor as f32;
}

pub unsafe fn tanimoto_array(a: &[u64; 4], b: &[u64; 4]) -> f32 {
    let mut dividend: u32 = 0;
    let mut divisor: u32 = 0;

    for i in 0..4 {
        dividend += ((a[i] & b[i]) as i64).count_ones();
        divisor += ((a[i] | b[i]) as i64).count_ones();
    }
    return dividend as f32 / divisor as f32;
}

#[cfg(test)]
mod tests {
    use bit_vec::BitVec;
    use crate::ringo::math::similarity::tanimoto::{tanimoto_array, tanimoto_vec};

    #[test]
    fn test_tanimoto_vec_033() {
        let a: BitVec = BitVec::from_bytes(&[0b00000101]);
        let b = BitVec::from_bytes(&[0b00000011]);

        assert_eq!(tanimoto_vec(&a, &b), 0.33333334);
    }

    #[test]
    fn test_tanimoto_vec_05() {
        let a: BitVec = BitVec::from_bytes(&[0b0000001]);
        let b = BitVec::from_bytes(&[0b00000011]);

        assert_eq!(tanimoto_vec(&a, &b), 0.5);
    }

    #[test]
    fn test_tanimoto_array_033() {
        let a: [u64; 4] = [0b00000101, 0, 0, 0];
        let b = [0b00000011, 0, 0, 0];

        unsafe {
            assert_eq!(tanimoto_array(&a, &b), 0.33333334);
        }
    }

    #[test]
    fn test_tanimoto_array_05() {
        let a: [u64; 4] = [0b00000001, 0, 0, 0];
        let b = [0b00000011, 0, 0, 0];

        unsafe {
            assert_eq!(tanimoto_array(&a, &b), 0.5);
        }
    }
}