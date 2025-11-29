use crate::counters::Counter;
use std::collections::hash_map::RandomState;
use std::hash::BuildHasher;

const PHI: f64 = 0.77351;

pub struct FMCounter<S = RandomState> {
    size: usize,
    bitset: Vec<u8>,
    hasher: S,
}

impl<S: BuildHasher + Default> Counter for FMCounter<S> {
    fn new(size: usize) -> Self {
        FMCounter {
            size,
            bitset: vec![0; size.div_ceil(8)],
            hasher: S::default(),
        }
    }

    fn add(&mut self, item: &[u8]) {
        let hash = self.hasher.hash_one(item);

        let num_trailing_zeros = hash.trailing_zeros() as usize;
        let index = std::cmp::min(num_trailing_zeros, self.size - 1) as usize;
        self.bitset[index / 8] |= 1 << (index % 8);
    }

    fn estimate(&self) -> f64 {
        let first_zero_bit = self
            .bitset
            .iter()
            .enumerate()
            .filter_map(|(idx, &byte)| {
                if byte == u8::MAX {
                    None
                } else {
                    Some(idx * 8 + byte.trailing_ones() as usize)
                }
            })
            .next()
            .unwrap_or(self.size - 1);

        (1_usize << first_zero_bit) as f64 / PHI
    }
}
