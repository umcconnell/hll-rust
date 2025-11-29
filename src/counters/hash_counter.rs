use crate::counters::Counter;
use std::collections::HashSet;
use std::collections::hash_map::RandomState;
use std::hash::BuildHasher;

pub struct HashCounter<S: BuildHasher + Default = RandomState> {
    hasher: S,
    counter: HashSet<u64>,
}

impl<S: BuildHasher + Default> Counter for HashCounter<S> {
    fn new(_size: usize) -> Self {
        HashCounter {
            hasher: S::default(),
            counter: HashSet::new(),
        }
    }

    fn add(&mut self, item: &[u8]) {
        let hash = self.hasher.hash_one(item);
        self.counter.insert(hash);
    }

    fn estimate(&self) -> f64 {
        self.counter.len() as f64
    }
}
