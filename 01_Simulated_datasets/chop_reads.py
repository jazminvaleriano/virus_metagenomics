import random
import pandas as pd

# ----------------------
# Settings
subset= "LONG"
input_fastq = f"00_mock_control_reads/Zymo/{subset}/02_Zymo-200K_{subset}.fastq"  # Uncompressed input FASTQ
length_distribution_file = f"read_lengths/real_DS_filt.tsv"  # Output from script 00_get_read_lengths.sh for length filterd reads
output_fastq = f"00_mock_control_reads/Zymo/{subset}/length-adj/03_Zymo-200K_len-adj_{subset}.fastq"  # Output
# ----------------------

# Load real read length distribution from TSV
try:
    df_lengths = pd.read_csv(length_distribution_file, sep="\t")
    target_lengths = df_lengths["read_length"].astype(int).tolist()
except Exception as e:
    raise ValueError(f"ERROR: Couldn't load valid lengths from {length_distribution_file}.\n{e}")

if not target_lengths:
    raise ValueError("ERROR: Length distribution file is empty or doesn't contain valid lengths.")

# Stats counters
total_reads = 0
total_written = 0
total_cut = 0

# Process FASTQ
with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
    while True:
        header = infile.readline()
        if not header:
            break  # EOF
        seq = infile.readline().strip()
        plus = infile.readline()
        qual = infile.readline().strip()

        total_reads += 1
        read_len = len(seq)
        target_len = random.choice(target_lengths)

        if read_len <= target_len:
            outfile.write(header)
            outfile.write(seq + '\n')
            outfile.write(plus)
            outfile.write(qual + '\n')
            total_written += 1
        else:
            start = random.randint(0, read_len - target_len)
            end = start + target_len
            new_header = f"{header.strip()}_frag_{start}_{end}\n"
            new_seq = seq[start:end]
            new_qual = qual[start:end]

            outfile.write(new_header)
            outfile.write(new_seq + '\n')
            outfile.write(plus)
            outfile.write(new_qual + '\n')

            total_written += 1
            total_cut += 1

# Report
print(f"Finished processing {total_reads} reads.")
print(f"Total written reads: {total_written}")
print(f"Reads cut and fragmented: {total_cut}")
print(f"Output saved as: {output_fastq}")
