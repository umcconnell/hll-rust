use crate::counters::Counter;
use std::collections::hash_map::RandomState;
use std::hash::BuildHasher;

const AM_4: f64 = 0.673;
const AM_5: f64 = 0.697;
const AM_6: f64 = 0.709;

pub struct HLLCounter<S = RandomState> {
    size: usize,
    am: f64,
    registers: Vec<u8>,
    hasher: S,
}

impl<S: BuildHasher + Default> Counter for HLLCounter<S> {
    fn new(size: usize) -> Self {
        let num_registers = 1 << size;
        let am = match size {
            0..=4 => AM_4,
            5 => AM_5,
            6 => AM_6,
            _ => 0.7213 / (1.0 + 1.079 / num_registers as f64),
        };
        HLLCounter {
            size,
            am,
            registers: vec![u8::MIN; num_registers],
            hasher: S::default(),
        }
    }

    fn add(&mut self, item: &[u8]) {
        let hash = self.hasher.hash_one(item);
        self.add_hash(hash);
    }

    fn estimate(&self) -> f64 {
        let num_registers = (1 << self.size) as f64;

        let numerator = self.am * num_registers * num_registers;

        let denominator: f64 = self
            .registers
            .iter()
            .map(|&reg| 2f64.powi(-(reg as i32)))
            .sum();

        let mut estimate = numerator / denominator;

        // Small range correction
        if estimate <= 2.5 * num_registers {
            let zeros = self.registers.iter().filter(|&&reg| reg == 0).count();
            if zeros > 0 {
                estimate = num_registers * (num_registers / zeros as f64).ln();
            }
        } else if estimate > (2f64.powi(64) / 30f64) {
            estimate = -2f64.powi(64) * (1f64 - estimate * 2f64.powi(-64)).ln()
        }

        estimate
    }
}

impl<S: BuildHasher + Default> HLLCounter<S> {
    // Some specialized high-performance methods
    #[inline(always)]
    pub fn add_u64(&mut self, item: u64) {
        let hash = self.hasher.hash_one(item);
        self.add_hash(hash);
    }

    #[inline(always)]
    fn add_hash(&mut self, hash: u64) {
        let index = (hash & ((1u64 << self.size) - 1)) as usize;
        let remainder = hash >> self.size;
        // trailing_zeros() will usually be compiled to a single instruction
        // like BSF on x86 architectures
        // see this example: https://godbolt.org/z/eGejof3Kz
        let rho = std::cmp::min(remainder.trailing_zeros() + 1, 64 - self.size as u32) as u8;

        self.registers[index] = std::cmp::max(self.registers[index], rho);
    }

    pub fn merge(&mut self, other: &HLLCounter<S>) {
        assert_eq!(self.size, other.size);
        for (reg_self, reg_other) in self.registers.iter_mut().zip(other.registers.iter()) {
            *reg_self = std::cmp::max(*reg_self, *reg_other);
        }
    }
}
