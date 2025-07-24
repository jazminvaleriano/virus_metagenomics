#!/bin/bash
#SBATCH --job-name=metrics_map
#SBATCH --error=logs/metrics_map_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

module load BEDTools/2.30.0-GCC-10.3.0
module load SAMtools/1.13-GCC-10.3.0

#----Generate per-base depth files-------------

BAM_DIR=02_remapping_results
OUTDIR=03_coverage_metrics/depth_files

mkdir -p "$OUTDIR"

module load BEDTools/2.30.0-GCC-10.3.0

for BAM in "$BAM_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM" .bam)
    bedtools genomecov -ibam "$BAM" -dz > "$OUTDIR/${SAMPLE}_depth.txt"
done

# -------- Generate mapping per accession --------
BAM_DIR=02_remapping_results
OUTDIR=03_coverage_metrics/accession_mapping_stats
mkdir -p "$OUTDIR"

module load SAMtools/1.13-GCC-10.3.0

for BAM in "$BAM_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM" .bam)
    OUTFILE="${OUTDIR}/${SAMPLE}_accession_stats.tsv"

    echo "Processing $SAMPLE..."

    # Estimate average read length from first 1000 reads
    AVG_READ_LEN=$(samtools view "$BAM" | head -n 1000 | \
        awk '{sum += length($10)} END {if (NR > 0) print sum / NR; else print 0}')

    # Extract per-accession stats from BAM
    samtools idxstats "$BAM" | awk -v avglen="$AVG_READ_LEN" 'BEGIN {
        OFS = "\t";
        print "Accession", "GenomeLength", "MappedReads", "AvgReadLength";
    }
    {
        print $1, $2, $3, avglen;
    }' > "$OUTFILE"

    echo "  → Saved to $OUTFILE"
done