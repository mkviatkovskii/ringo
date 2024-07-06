use fixedbitset::FixedBitSet;

pub fn tanimoto_bitset(a: &FixedBitSet, b: &FixedBitSet) -> f32 {
    let mut and_ = a.clone();
    and_.intersect_with(b);
    return and_.count_ones(..) as f32
        / (a.count_ones(..) + b.count_ones(..) - and_.count_ones(..)) as f32;
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
}
