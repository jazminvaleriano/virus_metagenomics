#!/bin/bash
#SBATCH --job-name=spike_subsamples
#SBATCH --output=logs/spike_subsamples_%A.out
#SBATCH --error=logs/spike_subsamples_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate kaiju_env

# -------- CONFIG --------
MOCK_DIRECTORY=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/00_mock_control_reads
BACKGROUND=$MOCK_DIRECTORY/Zymo/LONG/03_Zymo-200K_length_adjusted.fastq
VIRAL=$MOCK_DIRECTORY/viral_reads_sim/viral_sim_sample0_aligned_reads.fastq
OUTDIR=$MOCK_DIRECTORY/viral_spiked
mkdir -p "$OUTDIR"
# ------------------------

# Format: NAME  VIRAL_COUNT  BACKGROUND_COUNT
declare -a DATASETS=(
    "spike_5p_20K     1000    19000"
    "spike_2p_20K      400    19600"
    "spike_1p_20K      200    19800"
    "spike_0.5p_20K    100    19900"
    "spike_0.1p_20K     20    19980"
)

for entry in "${DATASETS[@]}"; do
    read NAME VCOUNT BCOUNT <<< "$entry"
    echo "Generating $NAME with $VCOUNT viral reads and $BCOUNT background reads..."

    # Subsample reads
    seqtk sample -s100 "$VIRAL" "$VCOUNT" > "$OUTDIR/${NAME}_viral.fastq"
    seqtk sample -s100 "$BACKGROUND" "$BCOUNT" > "$OUTDIR/${NAME}_bg.fastq"

    # Combine
    cat "$OUTDIR/${NAME}_bg.fastq" "$OUTDIR/${NAME}_viral.fastq" > "$OUTDIR/${NAME}.fastq"

    # Shuffle reads
    seqkit shuffle "$OUTDIR/${NAME}.fastq" -o "$OUTDIR/${NAME}_shuffled.fastq"

    # Clean up intermediates
    #rm "$OUTDIR/${NAME}_bg.fastq" "$OUTDIR/${NAME}_viral.fastq" "$OUTDIR/${NAME}.fastq"

    echo "✅ $NAME created: $OUTDIR/${NAME}_shuffled.fastq"
done

echo "All spike-in datasets generated in $OUTDIR"
