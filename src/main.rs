mod demo;

use xxhash_rust::xxh64::Xxh64Builder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sample_dataset = [
        ("SARS-CoV-2", "data/SARS-CoV-2/NC_045512v2.fa"),
        ("Thale Cress", "data/ThaleCress/TAIR9_chr_all.fa"),
        ("Zebrafish", "data/Zebrafish/danRer11.fa"),
        ("Human", "data/Human/hs1.fa"),
    ];

    // The following examples use the Xxh64 hash function.
    // You can use a different hash function by providing a different BuildHasher.
    //
    // For example, to use the default RandomState hasher:
    // let mut counter = FMCounter::<RandomState>::new(32);

    // Generate the comparison plot
    println!("Synthetic data plot");
    println!("===================");
    demo::synthetic::plot_comparison::<Xxh64Builder>(true)?;

    println!();
    println!("Real biological data (parallel)");
    println!("===============================");
    // Optionally run single-threaded analysis
    // println!("Real biological data");
    // demo::biological::run_sequential::<Xxh64Builder>(&sample_dataset, false)?;
    demo::biological::run_parallel::<Xxh64Builder>(&sample_dataset, false)?;

    Ok(())
}
