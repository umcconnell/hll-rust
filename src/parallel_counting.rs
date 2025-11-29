use crate::Counter;
use crate::HLLCounter;
use crate::fasta::FastaReader;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader};

// A=00, C=01, G=10, T=11
const ENCODING: [u8; 256] = {
    let mut table = [0xFF; 256];
    table[b'A' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'T' as usize] = 3;
    // Handle lowercase as well, though we uppercase the sequence
    table[b'a' as usize] = 0;
    table[b'c' as usize] = 1;
    table[b'g' as usize] = 2;
    table[b't' as usize] = 3;
    table
};

const K_MER_LENGTH: usize = 31;
const K_MER_MASK: u64 = (1u64 << 2 * K_MER_LENGTH) - 1; // Mask for 31-mer (62 bits)

#[inline(always)]
fn get_canonical_u64(kmer: u64) -> u64 {
    // Reverse complement for 2-bit encoding (A=00, C=01, G=10, T=11)
    // 1. Reverse bits
    // 2. Shift right by 2 (since we use 62 bits for 31-mer)
    // 3. Swap adjacent bits (to fix 2-bit chunk order)
    // 4. XOR with mask (to complement)

    let mut r = kmer.reverse_bits();
    r >>= 64 - 2 * K_MER_LENGTH; // Align to LSB (64 - 2*K_MER_LENGTH)

    // Swap adjacent bits: (r >> 1) & 0x55... | (r & 0x55...) << 1
    // 0x5555... is the mask 0101..., allowing us to select every 2nd bit
    r = ((r >> 1) & 0x5555555555555555) | ((r & 0x5555555555555555) << 1);

    // Complement: XOR with 11...11 (62 bits)
    // 11 binary is 3 decimal. We want to XOR each 2-bit pair with 11.
    // So we XOR with all ones (masked to 62 bits).
    r ^= (1u64 << 2 * K_MER_LENGTH) - 1;

    if kmer < r { kmer } else { r }
}

pub fn run_parallel_fasta_analysis<S: std::hash::BuildHasher + Default + Send + Sync>(
    path: &str,
) -> io::Result<(u64, HLLCounter<S>)> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut fasta_reader = FastaReader::new(reader);

    let sequences = std::iter::from_fn(move || match fasta_reader.next_record() {
        Ok(true) => match fasta_reader.read_sequence() {
            Ok(seq) => Some(Ok(seq)),
            Err(e) => Some(Err(e)),
        },
        Ok(false) => None,
        Err(e) => Some(Err(e)),
    });

    let final_counter = sequences
        .par_bridge()
        .map(|res| {
            let seq = res.expect("Error reading sequence");
            let mut counter = HLLCounter::<S>::new(16);
            let mut kmers_seen: u64 = 0;

            // Fast path using u64 for 31-mers
            // We use a rolling window with 2-bit encoding
            let mut kmer_u64: u64 = 0;
            let mut valid_len = 0;

            for &byte in seq.iter() {
                let code = ENCODING[byte as usize];
                if code == 0xFF {
                    // Skip unknown characters
                    valid_len = 0;
                    kmer_u64 = 0;
                } else {
                    kmer_u64 = ((kmer_u64 << 2) & K_MER_MASK) | (code as u64);
                    valid_len += 1;

                    if valid_len >= K_MER_LENGTH {
                        let canonical = get_canonical_u64(kmer_u64);
                        counter.add_u64(canonical);
                        kmers_seen += 1;
                    }
                }
            }

            (kmers_seen, counter)
        })
        .reduce(
            || (0, HLLCounter::<S>::new(16)),
            |(count_a, mut a), (count_b, b)| {
                a.merge(&b);
                (count_a + count_b, a)
            },
        );

    Ok(final_counter)
}
