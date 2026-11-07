use fixedbitset::FixedBitSet;

/// Returns 0.0 when both sets are empty, so empty records never match.
pub fn tanimoto_bitset(a: &FixedBitSet, b: &FixedBitSet) -> f32 {
    let union = a.union_count(b);
    if union == 0 {
        return 0.0;
    }
    a.intersection_count(b) as f32 / union as f32
}

#[cfg(test)]
mod tests {
    use crate::math::similarity::tanimoto::tanimoto_bitset;
    use fixedbitset::FixedBitSet;

    #[test]
    fn test_tanimoto_bitset_033() {
        let mut a = FixedBitSet::with_capacity(8);
        a.insert(0);
        a.insert(2);
        let mut b = FixedBitSet::with_capacity(8);
        b.insert(0);
        b.insert(1);
        assert_eq!(tanimoto_bitset(&a, &b), 0.33333334);
    }

    #[test]
    fn test_tanimoto_bitset_05() {
        let mut a = FixedBitSet::with_capacity(8);
        a.insert(0);
        let mut b = FixedBitSet::with_capacity(8);
        b.insert(0);
        b.insert(1);
        assert_eq!(tanimoto_bitset(&a, &b), 0.5);
    }

    #[test]
    fn test_tanimoto_bitset_empty() {
        // used to be 0/0 = NaN
        let a = FixedBitSet::with_capacity(8);
        let b = FixedBitSet::with_capacity(8);
        assert_eq!(tanimoto_bitset(&a, &b), 0.0);
    }
}
