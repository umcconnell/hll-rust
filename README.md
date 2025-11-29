# hll-rust

A Rust implementation of probabilistic data structures for k-mer counting in biological sequences. This project explores and compares different cardinality estimation algorithms against exact counting methods.

## Overview

`hll-rust` implements several probabilistic counters to estimate the number of distinct elements (cardinality) in large datasets, specifically focusing on k-mers in DNA sequences (FASTA files). For HyperLogLog in particular, you can read more on Wikipedia (https://en.wikipedia.org/wiki/HyperLogLog) and in the original publication by Flajolet et al. [[1]].

## Features

*   **Probabilistic Counters**:
    *   **HyperLogLog (HLL)**: State-of-the-art cardinality estimation with low memory footprint.
    *   **Flajolet-Martin**: A classic probabilistic counting algorithm.
    *   **Linear Counter**: Efficient for smaller cardinalities.
*   **Exact Counting**:
    *   **Hash Counter**: Baseline exact counter for validation.
*   **High Performance**:
    *   Written in **Rust**.
    *   Uses **xxHash** (`xxh64`) for fast, high-quality hashing.
    *   **Parallel Processing**: Multi-threaded processing of FASTA files using rayon for maximum throughput.
*   **Visualization**: Automatically generates performance comparison plots (`counter_comparison.png`).

## Project Structure

*   `src/counters/`: Implementations of the various counting algorithms.
*   `src/parallel_counting.rs`: Logic for parallel processing of datasets.
*   `src/fasta.rs`: FASTA file parsing utilities.

## Usage

Ensure you have Rust and Cargo installed.

To run the analysis (synthetic benchmarks and biological data processing):

```bash
cargo run --release
```

This command will:
1.  Run synthetic benchmarks comparing Linear, FM, and HLL counters.
2.  Generate a plot `counter_comparison.png`.
3.  Process the configured biological datasets (FASTA files) in parallel.

### Custom Hash Function

The counters in this library are generic over the hash function. By default, the examples use `xxHash` (`Xxh64Builder`) for performance, but you can easily swap it for any other hasher that implements `std::hash::BuildHasher`.

For example, to use the standard library's `RandomState`:

```rust
use std::hash::RandomState;
use hll_rust::FMCounter;

// Initialize a Flajolet-Martin counter with RandomState
let mut counter = FMCounter::<RandomState>::new(20);
```

Or to use the default `Xxh64Builder` as in the main example:

```rust
use xxhash_rust::xxh64::Xxh64Builder;
use hll_rust::HLLCounter;

// Initialize a HyperLogLog counter with xxHash
let mut counter = HLLCounter::<Xxh64Builder>::new(12);
```

## Results

The following table compares the complexity estimates (distinct k-mers / total k-mers) obtained using our HyperLogLog implementation against the ground truth calculated by [Jellyfish](https://github.com/gmarcais/Jellyfish).

| Dataset     | Ground truth complexity | HLL complexity | Relative error |
|-------------|-------------------------|----------------|----------------|
|  SARS-CoV-2 |                  0.9999 |         1.0000 |         +0.01% | 
| Thale Cress |                  0.9368 |         0.9381 |         +0.14% |
|   Zebrafish |                  0.6694 |         0.6706 |         +0.18% |
|       Human |                  0.8060 |         0.8051 |         -0.11% |

*Note: Ground truth complexity computed using Jellyfish. Sample run using `xxh3`.*

### Generating Ground Truth

To generate the ground truth data using Jellyfish, we provide a script that runs the analysis in a Docker container. As prerequisite, Docker must be installed and running.

To run the Jellyfish analysis with default settings (looking for data in `data/` and outputting to `jellyfish_results/`):

```bash
./run_jellyfish.sh
```

You can also specify custom data and results directories:

```bash
./run_jellyfish.sh <path_to_data> <path_to_results>
```

This script will:
1.  Build a Docker image containing Jellyfish.
2.  Process all `.fa` files found in the specified data directory.
3.  Output the results (counts and stats) to the specified results directory.

## Data Availability

The biological datasets used for this project are internal and cannot be shared publicly. However, the code is designed to work with standard FASTA files, so you can easily adapt it to run on your own data.

## References

[1] P. Flajolet, É. Fusy, O. Gandouet, and F. Meunier, “HyperLogLog: the analysis of a near-
optimal cardinality estimation algorithm,” in Proceedings of the International Conference on
Analysis of Algorithms (AOFA), pp. 127–146, 2007.
https://algo.inria.fr/flajolet/Publications/FlFuGaMe07.pdf

## License

This project is released under the [MIT License](LICENSE.md).