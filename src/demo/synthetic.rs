use hll_rust::{Counter, FMCounter, HLLCounter, LinearCounter};
use plotters::prelude::*;
use rayon::prelude::*;

pub type SeedData = (u64, Vec<(f64, f64)>);

fn process_seed<S: std::hash::BuildHasher + Default>(
    seed: u64,
    ns: &[u64],
) -> (SeedData, SeedData, SeedData) {
    let mut linear_points = Vec::new();
    let mut fm_points = Vec::new();
    let mut hll_points = Vec::new();

    let mut linear_counter: LinearCounter<S> = LinearCounter::new(1 << 20);
    let mut fm_counter: FMCounter<S> = FMCounter::new(32);
    let mut hll_counter: HLLCounter<S> = HLLCounter::new(20);

    let mut last_n = 0;
    for &n in ns {
        for i in last_n..n {
            let value = i ^ seed;
            let bytes = value.to_le_bytes();
            linear_counter.add(&bytes);
            fm_counter.add(&bytes);
            hll_counter.add(&bytes);
        }
        last_n = n;

        linear_points.push((n as f64, linear_counter.estimate()));
        fm_points.push((n as f64, fm_counter.estimate()));
        hll_points.push((n as f64, hll_counter.estimate()));
    }

    ((seed, linear_points), (seed, fm_points), (seed, hll_points))
}

pub fn collect_test_data_sequential<S: std::hash::BuildHasher + Default>()
-> (Vec<SeedData>, Vec<SeedData>, Vec<SeedData>) {
    let seeds: Vec<u64> = (1..=9).collect();
    let ns: Vec<u64> = (0..25).map(|i| 1u64 << i).collect();

    let mut linear_data = Vec::new();
    let mut fm_data = Vec::new();
    let mut hll_data = Vec::new();

    for &seed in &seeds {
        let (l, f, h) = process_seed::<S>(seed, &ns);
        linear_data.push(l);
        fm_data.push(f);
        hll_data.push(h);
    }

    (linear_data, fm_data, hll_data)
}

pub fn collect_test_data_parallel<S: std::hash::BuildHasher + Default + Send + Sync>()
-> (Vec<SeedData>, Vec<SeedData>, Vec<SeedData>) {
    let seeds: Vec<u64> = (1..=9).collect();
    let ns: Vec<u64> = (0..25).map(|i| 1u64 << i).collect();

    let results: Vec<_> = seeds
        .par_iter()
        .map(|&seed| process_seed::<S>(seed, &ns))
        .collect();

    let mut linear_data = Vec::new();
    let mut fm_data = Vec::new();
    let mut hll_data = Vec::new();

    for (l, f, h) in results {
        linear_data.push(l);
        fm_data.push(f);
        hll_data.push(h);
    }

    (linear_data, fm_data, hll_data)
}

pub fn plot_comparison<S: std::hash::BuildHasher + Default + Send + Sync>(
    parallel: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Collecting test data (parallel={})...", parallel);
    let (linear_data, fm_data, hll_data) = if parallel {
        collect_test_data_parallel::<S>()
    } else {
        collect_test_data_sequential::<S>()
    };

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
