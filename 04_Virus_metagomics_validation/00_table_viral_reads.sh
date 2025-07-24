#!/bin/bash

# Input FASTQ file
INPUT_FASTQ="00_simulated_reads/full_sim_om_200000.fastq"

# Output TSV file (final result)
OUTPUT_TSV="200K_virus_reads_chimeric.tsv"

# Extract and annotate virus-containing reads
grep '^@' "$INPUT_FASTQ" | awk '
BEGIN {
    OFS="\t";
    # Mapping virus tags to pretty names
    map["Betaarterivirus"] = "Betaarterivirus suid 1";
    map["Betacoronavirus"] = "Betacoronavirus 1";
    map["betaherpesvirus"] = "Suid betaherpesvirus";
    map["parvovirus-NC-001718"] = "Porcine parvovirus";
    map["Virus"] = "Influenza A virus";
}
{
    header = $0;
    sub(/^@/, "", header);                        # Remove '@'
    read_id = header;                             # Full read ID
    is_chimeric = "not_chimeric";
    if (tolower(header) ~ /chimeric/) {
        is_chimeric = "chimeric";
    }

    # Extract virus tag
    n = split(header, parts, /[;_]/);
    for (i = 1; i <= n; i++) {
        if (tolower(parts[i]) ~ /virus/) {
            tag = parts[i];
            pretty = (tag in map) ? map[tag] : "Unknown";
            print read_id, tag, is_chimeric, pretty;
        }
    }
}' > "$OUTPUT_TSV"

echo "Done. Results saved to $OUTPUT_TSV"
