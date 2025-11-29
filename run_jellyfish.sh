#!/bin/bash
set -e

DATA_DIR=${1:-data}
RESULTS_DIR=${2:-jellyfish_results}

# Build the docker image
docker build -t jellyfish -f JellyfishDocker .

# Create results directory
mkdir -p "$RESULTS_DIR"

# Find all .fa files in data/
find "$DATA_DIR" -name "*.fa" | while read file; do
    echo "Processing $file..."
    
    filename=$(basename "$file")
    dirname=$(dirname "$file")
    # Get relative path (e.g., Human)
    rel_dir=${dirname#$DATA_DIR/}
    
    # Create output directory in results
    mkdir -p "$RESULTS_DIR/$rel_dir"
    
    # Define paths
    # Inside container, we map $(pwd)/$(DATA_DIR) to /data
    # So data/Human/hs1.fa becomes /data/Human/hs1.fa
    input_path="/data/$rel_dir/$filename"
    
    # Output path inside container
    # We map $(pwd)/results to /results
    output_jf="/results/$rel_dir/${filename%.*}.jf"
    
    # Host path for stats file (since we redirect stdout)
    stats_file="$RESULTS_DIR/$rel_dir/${filename%.*}.stats"
    
    echo "Counting kmers for $filename..."
    docker run --rm \
        -v "$(pwd)/$DATA_DIR:/data" \
        -v "$(pwd)/$RESULTS_DIR:/results" \
        jellyfish jellyfish count \
        -m 31 \
        -s 1G \
        -t 6 \
        -C \
        -o "$output_jf" \
        "$input_path"
        
    echo "Computing stats for $filename..."
    docker run --rm \
        -v "$(pwd)/$RESULTS_DIR:/results" \
        jellyfish jellyfish stats \
        "$output_jf" > "$stats_file"
        
    echo "Stats saved to $stats_file"
done
