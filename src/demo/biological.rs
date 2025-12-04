use hll_rust::fasta::FastaReader;
use hll_rust::parallel_counting;
use hll_rust::{Counter, FMCounter, HLLCounter, HashCounter, LinearCounter};
use std::fs::File;
use std::io::{self, BufReader};

pub fn run_sequential<S: std::hash::BuildHasher + Default>(
    dataset: &[(&str, &str)],
    verbose: bool,
) -> io::Result<()> {
    // Store (length)
    let mut hll_estimated_complexity: Vec<f64> = Vec::new();

    for (name, path) in dataset.iter() {
        println!("Processing dataset: {}", name);

        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut fasta_reader = FastaReader::new(reader);

        let mut linear_counter: LinearCounter<S> = LinearCounter::new(1_000_000);
        let mut hash_counter: HashCounter<S> = HashCounter::new(0);
        let mut fm_counter: FMCounter<S> = FMCounter::new(32);
        let mut hll_counter: HLLCounter<S> = HLLCounter::new(14);

        let mut total_kmers_seen: u64 = 0;

        while fasta_reader.next_record()? {
            if let Some(id) = &fasta_reader.id
                && verbose
            {
                println!(">{}", String::from_utf8_lossy(id));
            }

            for kmer_result in fasta_reader.canonical_kmers(31) {
                let kmer = kmer_result?;
                if verbose {
                    linear_counter.add(&kmer);
                    fm_counter.add(&kmer);
                    hash_counter.add(&kmer);
                }
                hll_counter.add(&kmer);

                total_kmers_seen += 1;
            }
        }

        let hll_estimate = hll_counter.estimate();
        hll_estimated_complexity.push(hll_estimate / total_kmers_seen as f64);

        if verbose {
            let linear_estimate = linear_counter.estimate();
            let fm_estimate = fm_counter.estimate();
            let hash_estimate = hash_counter.estimate();

            println!(
                "\n{:<15} {:<15} {:<20} {:<15}",
                "Counter", "Complexity", "Estimate", "Rel Error (%)"
            );
            println!("{:-<65}", "");
            println!(
                "{:<15} {:<15.6} {:<20.0} {:<15.4}",
                "Linear",
                linear_estimate / total_kmers_seen as f64,
                linear_estimate.round(),
                100f64 * (linear_estimate - hash_estimate) / hash_estimate
            );
            println!(
                "{:<15} {:<15.6} {:<20.0} {:<15.4}",
                "FM",
                fm_estimate / total_kmers_seen as f64,
                fm_estimate.round(),
                100f64 * (fm_estimate - hash_estimate) / hash_estimate
            );
            println!(
                "{:<15} {:<15.6} {:<20.0} {:<15.4}",
                "HLL",
                hll_estimate / total_kmers_seen as f64,
                hll_estimate.round(),
                100f64 * (hll_estimate - hash_estimate) / hash_estimate
            );
            println!(
                "{:<15} {:<15.6} {:<20.0} {:<15.4}",
                "True (Hash)",
                hash_estimate / total_kmers_seen as f64,
                hash_estimate.round(),
                0.0
            );
            println!();
        }
    }

    for (est, (name, _)) in hll_estimated_complexity.iter().zip(dataset.iter()) {
        println!("Dataset: {}, HLL estimated complexity: {:.6}", name, est);
    }

    Ok(())
}

pub fn run_parallel<S: std::hash::BuildHasher + Default + Send + Sync>(
    dataset: &[(&str, &str)],
    _verbose: bool,
) -> io::Result<()> {
    println!(
        "\n{:<20} | {:<15} | {:<15} | {:<15}",
        "Dataset", "Complexity", "Total K-mers", "Time"
    );
    println!("{:-<80}", "");

    for (name, path) in dataset.iter() {
        // println!("Processing dataset: {}", name);
        let start = std::time::Instant::now();
        let (total_count, counter) = parallel_counting::run_parallel_fasta_analysis::<S>(path)?;
        let duration = start.elapsed();

        let unique_count_estimate = counter.estimate();
        let complexity_estimate =
            f64::clamp(unique_count_estimate / (total_count as f64), 0.0, 1.0);

        println!(
            "{:<20} | {:<15.4} | {:<15} | {:?}",
            name, complexity_estimate, total_count, duration
        );
    }
    println!();
    Ok(())
}
