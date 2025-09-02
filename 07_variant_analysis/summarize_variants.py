import os
import glob
import csv

# Input/output paths
medaka_dir = "05_variant_calling/medaka"
output_csv = "05_variant_calling/variant_summary_medaka.csv"

with open(output_csv, "w", newline="") as out_csv:
    writer = csv.writer(out_csv)
    writer.writerow(["Sample", "Virus", "Chrom", "Pos", "Ref", "Alt", "Type", "DP", "GQ"])

    for vcf_path in glob.glob(f"{medaka_dir}/*/medaka.annotated.vcf"):
        sample_virus = os.path.basename(os.path.dirname(vcf_path))
        if "_" not in sample_virus:
            continue
        sample, virus = sample_virus.split("_", 1)

        with open(vcf_path) as vcf:
            for line in vcf:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                chrom = fields[0]
                pos = fields[1]
                ref = fields[3]
                alt = fields[4]
                info = fields[7]
                format_field = fields[8]
                sample_field = fields[9]

                # Determine variant type
                if len(ref) == 1 and len(alt) == 1:
                    var_type = "SNP"
                elif len(ref) < len(alt):
                    var_type = "Insertion"
                elif len(ref) > len(alt):
                    var_type = "Deletion"
                else:
                    var_type = "Complex"

                # Get depth (DP) from INFO
                dp = "-"
                for item in info.split(";"):
                    if item.startswith("DP="):
                        dp = item.split("=")[1]

                # Get GQ (genotype quality) from FORMAT and SAMPLE columns
                gq = "-"
                format_keys = format_field.split(":")
                sample_values = sample_field.split(":")
                if "GQ" in format_keys:
                    gq_index = format_keys.index("GQ")
                    if gq_index < len(sample_values):
                        gq = sample_values[gq_index]

                writer.writerow([sample, virus, chrom, pos, ref, alt, var_type, dp, gq])
