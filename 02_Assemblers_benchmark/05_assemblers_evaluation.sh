#!/bin/bash

RESULTS_DIR=results/03_assemblies_evaluation
TOOLS="flye raven megahit"
COMBINED_REPORT=$RESULTS_DIR/combined_report.tsv
GENOMES_REPORT=$RESULTS_DIR/combined_gen_recov_report.tsv
NGA50_REPORT=$RESULTS_DIR/combined_nga50_report.tsv
GENOMES_LENGTH=$RESULTS_DIR/genomes_length.tsv

# Write the header 
echo -e "Assembly\t# contigs (>= 0 bp)\t# contigs (>= 1000 bp)\t# contigs (>= 5000 bp)\t# contigs (>= 10000 bp)\t# contigs (>= 25000 bp)\t# contigs (>= 50000 bp)\tTotal length (>= 0 bp)\tTotal length (>= 1000 bp)\tTotal length (>= 5000 bp)\tTotal length (>= 10000 bp)\tTotal length (>= 25000 bp)\tTotal length (>= 50000 bp)\t# contigs\tLargest contig\tTotal length\tReference length\tN50\tN90\tauN\tL50\tL90\t# misassemblies\t# misassembled contigs\tMisassembled contigs length\t# local misassemblies\t# scaffold gap ext. mis.\t# scaffold gap loc. mis.\t# unaligned mis. contigs\t# unaligned contigs\tUnaligned length\tGenome fraction (%)\tDuplication ratio\t# N's per 100 kbp\t# mismatches per 100 kbp\t# indels per 100 kbp\tLargest alignment\tTotal aligned length\tNA50\tNA90\tauNA\tLA50\tLA90\tAssembler" > $COMBINED_REPORT

# Process each assembler's report
for tool in $TOOLS; do 
    REPORT_PATH="$RESULTS_DIR/$tool/combined_reference/transposed_report.tsv"
    
    if [[ -f "$REPORT_PATH" ]]; then
        tail -n +2 "$REPORT_PATH" | while IFS= read -r line; do
            echo -e "${line}\t${tool}" >> $COMBINED_REPORT
        done
    else
        echo "Warning: Report not found for $tool at $REPORT_PATH"
    fi
done

# Genome recovery report

# 1. Detectar todas las columnas posibles
ALL_COLUMNS=()
for tool in $TOOLS; do
    REPORT_PATH="$RESULTS_DIR/$tool/summary/TSV/Genome_fraction.tsv"
    if [[ -f "$REPORT_PATH" ]]; then
        header=$(head -n 1 "$REPORT_PATH")
        IFS=$'\t' read -r -a cols <<< "$header"
        for col in "${cols[@]:1}"; do  # omitimos "Assembly"
            if [[ ! " ${ALL_COLUMNS[*]} " =~ " ${col} " ]]; then
                ALL_COLUMNS+=("$col")
            fi
        done
    fi
done

# 2. Escribir headers completos
{
    echo -ne "Assembly"
    for col in "${ALL_COLUMNS[@]}"; do
        echo -ne "\t${col}"
    done
    echo -e "\tAssembler"
} > "$GENOMES_REPORT"

# 3. Agregar filas por herramienta
for tool in $TOOLS; do
    REPORT_PATH="$RESULTS_DIR/$tool/summary/TSV/Genome_fraction.tsv"
    if [[ -f "$REPORT_PATH" ]]; then
        header=$(head -n 1 "$REPORT_PATH")
        IFS=$'\t' read -r -a tool_cols <<< "$header"

        # Map: column index -> column name
        declare -A col_index_map
        for i in "${!tool_cols[@]}"; do
            col_index_map["${tool_cols[$i]}"]=$i
        done

        tail -n +2 "$REPORT_PATH" | while IFS=$'\t' read -r -a values; do
            assembly="${values[0]}"
            output="$assembly"

            for col in "${ALL_COLUMNS[@]}"; do
                if [[ -n "${col_index_map[$col]}" ]]; then
                    idx=${col_index_map[$col]}
                    output+="\t${values[$idx]}"
                else
                    output+="\t"
                fi
            done
            output+="\t$tool"
            echo -e "$output" >> "$GENOMES_REPORT"
        done

        unset col_index_map
    else
        echo "Warning: Report not found for $tool at $REPORT_PATH"
    fi
done




# NGA50 report

# 1. Detect todas las columnas posibles
ALL_COLUMNS=()
for tool in $TOOLS; do
    REPORT_PATH="$RESULTS_DIR/$tool/summary/TSV/NGA50.tsv"
    if [[ -f "$REPORT_PATH" ]]; then
        header=$(head -n 1 "$REPORT_PATH")
        IFS=$'\t' read -r -a cols <<< "$header"
        for col in "${cols[@]:1}"; do  # omitimos "Assembly"
            if [[ ! " ${ALL_COLUMNS[*]} " =~ " ${col} " ]]; then
                ALL_COLUMNS+=("$col")
            fi
        done
    fi
done

# 2. Escribir headers completos
{
    echo -ne "Assembly"
    for col in "${ALL_COLUMNS[@]}"; do
        echo -ne "\t${col}"
    done
    echo -e "\tAssembler"
} > "$NGA50_REPORT"

# 3. Agregar filas por herramienta
for tool in $TOOLS; do
    REPORT_PATH="$RESULTS_DIR/$tool/summary/TSV/NGA50.tsv"
    if [[ -f "$REPORT_PATH" ]]; then
        header=$(head -n 1 "$REPORT_PATH")
        IFS=$'\t' read -r -a tool_cols <<< "$header"

        # Map: column index -> column name
        declare -A col_index_map
        for i in "${!tool_cols[@]}"; do
            col_index_map["${tool_cols[$i]}"]=$i
        done

        tail -n +2 "$REPORT_PATH" | while IFS=$'\t' read -r -a values; do
            assembly="${values[0]}"
            output="$assembly"

            for col in "${ALL_COLUMNS[@]}"; do
                if [[ -n "${col_index_map[$col]}" ]]; then
                    idx=${col_index_map[$col]}
                    output+="\t${values[$idx]}"
                else
                    output+="\t"
                fi
            done
            output+="\t$tool"
            echo -e "$output" >> "$NGA50_REPORT"
        done

        unset col_index_map
    else
        echo "Warning: Report not found for $tool at $REPORT_PATH"
    fi
done


# Genome lengths

# write headers
echo -e "genome\tlength" > $GENOMES_LENGTH

for genome_dir in $RESULTS_DIR/megahit/runs_per_reference/*; do
    genome_report=$genome_dir/report.tsv
    length=$(awk -F'\t' 'NR==15 {print $2}' "$genome_report")
    genome=$(basename "$genome_dir")
    echo -e "${genome}\t${length}" >> $GENOMES_LENGTH
done
