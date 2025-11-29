use crate::counters::Counter;
use std::collections::hash_map::RandomState;
use std::hash::BuildHasher;

pub struct LinearCounter<S = RandomState> {
    bit_array: Vec<u8>,
    size: usize,
    hasher: S,
}

impl<S: BuildHasher + Default> Counter for LinearCounter<S> {
    fn new(size: usize) -> Self {
        LinearCounter {
            bit_array: vec![0; size.div_ceil(8)],
            size,
            hasher: S::default(),
        }
    }

    fn add(&mut self, item: &[u8]) {
        let hash = self.hasher.hash_one(item);

        let index = (hash % self.size as u64) as usize;
        self.bit_array[index / 8] |= 1 << (index % 8);
    }

    fn estimate(&self) -> f64 {
        let num_unset_bits = std::cmp::max(
            1,
            self.bit_array
                .iter()
                .map(|byte| byte.count_zeros() as usize)
                .sum::<usize>(),
        );

        self.size as f64 * (self.size as f64 / num_unset_bits as f64).ln()
    }
}
