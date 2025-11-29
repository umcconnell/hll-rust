use hll_rust::Counter;
use hll_rust::fasta::FastaReader;
use hll_rust::{FMCounter, HLLCounter, HashCounter, LinearCounter};

use plotters::prelude::*;
use std::fs::File;
use std::hash::RandomState;
use std::io::{self, BufReader};

use xxhash_rust::xxh64::Xxh64Builder;

fn test_linear<S: std::hash::BuildHasher + Default>(n: u64, seed: u64) -> f64 {
    let mut counter: LinearCounter<S> = LinearCounter::new(1 << 20);
    for i in 0..n {
        // XOR with seed to get different sequences for each seed
        let value = i ^ seed;
        counter.add(&value.to_le_bytes());
    }
    counter.estimate()
}

fn test_fm<S: std::hash::BuildHasher + Default>(n: u64, seed: u64) -> f64 {
    let mut counter: FMCounter<S> = FMCounter::new(32);

    for i in 0..n {
        let value = i ^ seed;
        counter.add(&value.to_le_bytes());
    }
    counter.estimate()
}

fn test_hll<S: std::hash::BuildHasher + Default>(n: u64, seed: u64) -> f64 {
    let mut counter: HLLCounter<S> = HLLCounter::new(20);
    for i in 0..n {
        let value = i ^ seed;
        counter.add(&value.to_le_bytes());
    }
    counter.estimate()
}

/// Collect test data for all counter types
/// Returns (linear_data, fm_data, hll_data) where each is Vec<(seed, Vec<(n, estimate)>)>
fn collect_test_data<S: std::hash::BuildHasher + Default>() -> (
    Vec<(u64, Vec<(f64, f64)>)>,
    Vec<(u64, Vec<(f64, f64)>)>,
    Vec<(u64, Vec<(f64, f64)>)>,
) {
    let seeds: Vec<u64> = (1..=9).collect();
    let ns: Vec<u64> = (0..25).map(|i| 1u64 << i).collect();

    let mut linear_data = Vec::new();

    let mut fm_data = Vec::new();
    let mut hll_data = Vec::new();

    for &seed in &seeds {
        let mut linear_points = Vec::new();
        let mut fm_points = Vec::new();
        let mut hll_points = Vec::new();

        for &n in &ns {
            linear_points.push((n as f64, test_linear::<S>(n, seed)));
            fm_points.push((n as f64, test_fm::<S>(n, seed)));
            hll_points.push((n as f64, test_hll::<S>(n, seed)));
        }

        linear_data.push((seed, linear_points));
        fm_data.push((seed, fm_points));
        hll_data.push((seed, hll_points));
    }

    (linear_data, fm_data, hll_data)
}

/// Plot comparison charts using plotters
fn plot_comparison<S: std::hash::BuildHasher + Default>() -> Result<(), Box<dyn std::error::Error>>
{
    println!("Collecting test data...");
    let (linear_data, fm_data, hll_data) = collect_test_data::<S>();

    // Find the max value across all data for consistent scaling
    let max_val = [&linear_data, &fm_data, &hll_data]
        .iter()
        .flat_map(|data| {
            data.iter()
                .flat_map(|(_, points)| points.iter().map(|(_, y)| *y))
        })
        .fold(0.0f64, f64::max);

    let max_n = 16777216.0f64;

    // Define colors for each seed (matching matplotlib default colors)
    let colors = [
        RGBColor(31, 119, 180),  // blue
        RGBColor(255, 127, 14),  // orange
        RGBColor(44, 160, 44),   // green
        RGBColor(214, 39, 40),   // red
        RGBColor(148, 103, 189), // purple
        RGBColor(140, 86, 75),   // brown
        RGBColor(227, 119, 194), // pink
        RGBColor(127, 127, 127), // gray
        RGBColor(188, 189, 34),  // olive
    ];

    // Create the plot with higher resolution
    let root = BitMapBackend::new("counter_comparison.png", (2400, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let areas = root.split_evenly((1, 3));

    let datasets = [
        ("LinearCounting", &linear_data),
        ("FM", &fm_data),
        ("HLL", &hll_data),
    ];

    for (idx, (area, (title, data))) in areas.iter().zip(datasets.iter()).enumerate() {
        let mut chart = ChartBuilder::on(area)
            .caption(*title, ("sans-serif", 32).into_font())
            .margin(15)
            .x_label_area_size(50)
            .y_label_area_size(80)
            .build_cartesian_2d(
                (1.0f64..max_n).log_scale(),
                (1.0f64..max_val * 1.5).log_scale(),
            )?;

        chart
            .configure_mesh()
            .x_desc("n")
            .y_desc("estimate")
            .label_style(("sans-serif", 18))
            .draw()?;

        // Draw the perfect counter line (y = x)
        chart.draw_series(LineSeries::new(
            vec![(1.0, 1.0), (max_n, max_n)],
            ShapeStyle::from(&BLACK).stroke_width(2),
        ))?;

        // Draw each seed's data
        for (i, (seed, points)) in data.iter().enumerate() {
            let color = colors[i % colors.len()];

            let series = LineSeries::new(points.clone(), color.stroke_width(3));

            // Only add legend for the last (rightmost) chart
            if idx == 2 {
                chart
                    .draw_series(series)?
                    .label(format!("seed {}", seed))
                    .legend(move |(x, y)| {
                        PathElement::new(vec![(x, y), (x + 30, y)], color.stroke_width(3))
                    });
            } else {
                chart.draw_series(series)?;
            }
        }

        // Configure legend only for the last chart
        if idx == 2 {
            chart
                .configure_series_labels()
                .position(SeriesLabelPosition::UpperLeft)
                .label_font(("sans-serif", 18))
                .border_style(BLACK)
                .background_style(WHITE.mix(0.8))
                .draw()?;
        }
    }

    root.present()?;
    println!("Plot saved to counter_comparison.png");

    Ok(())
}

fn run_fasta_analysis<S: std::hash::BuildHasher + Default>(
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

            println!("Estimated complexity / Estimated unique 10-mers (rel. error)");
            println!(
                "LinearCounter:\t{}\t{:.0}\t({:.4}%)",
                linear_estimate / total_kmers_seen as f64,
                linear_estimate.round(),
                100f64 * (linear_estimate - hash_estimate) / hash_estimate
            );
            println!(
                "FMCounter:\t{}\t{:.0}\t({:.4}%)",
                fm_estimate / total_kmers_seen as f64,
                fm_estimate.round(),
                100f64 * (fm_estimate - hash_estimate) / hash_estimate
            );
            println!(
                "HLL Counter:\t{}\t{:.0}\t({:.4}%)",
                hll_estimate / total_kmers_seen as f64,
                hll_estimate.round(),
                100f64 * (hll_estimate - hash_estimate) / hash_estimate
            );
            println!(
                "true value:\t{}\t{:.0}\t({:.4}%)",
                hash_estimate / total_kmers_seen as f64,
                hash_estimate.round(),
                0
            );
        }
    }

    for (est, (name, _)) in hll_estimated_complexity.iter().zip(dataset.iter()) {
        println!("Dataset: {}, HLL estimated complexity: {:.6}", name, est);
    }

    Ok(())
}

fn run_parallel_analysis<S: std::hash::BuildHasher + Default + Send + Sync>(
    dataset: &[(&str, &str)],
    _verbose: bool,
) -> io::Result<()> {
    for (name, path) in dataset.iter() {
        println!("Processing dataset: {}", name);
        let start = std::time::Instant::now();
        let (total_count, counter) =
            hll_rust::parallel_counting::run_parallel_fasta_analysis::<S>(path)?;
        let duration = start.elapsed();

        let unique_count_estimate = counter.estimate();
        let complexity_estimate =
            f64::clamp(unique_count_estimate / (total_count as f64), 0.0, 1.0);

        println!(
            "Dataset: {},\nComplexity estimate: {:.4} (out of {} kmers in total)\nTime: {:?}\n=================================",
            name, complexity_estimate, total_count, duration
        );
    }
    Ok(())
}

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
    println!("Synthetic data");
    plot_comparison::<Xxh64Builder>()?;

    // Optionally run single-threaded analysis
    // println!("Real biological data");
    // run_fasta_analysis::<Xxh64Builder>(&sample_dataset, false)?;

    println!("Real biological data (parallel)");
    run_parallel_analysis::<Xxh64Builder>(&sample_dataset, false)?;

    Ok(())
}
