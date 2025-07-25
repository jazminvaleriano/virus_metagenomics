import sys
import re
import pandas as pd

# Path to your log file
LOG_FILE = sys.argv[1]              
OUT_FILE = sys.argv[2]

# Regular expressions to extract information
cmd_pattern = re.compile(r'\[M::main\] CMD: minimap2 .*?/(barcode\d+)_length-filt\.fastq')
mapped_pattern = re.compile(r'\[M::worker_pipeline::[^\]]+\] mapped (\d+) sequences')
processed_pattern = re.compile(r'\[M::bam2fq_mainloop\] processed (\d+) reads')

# Read the log file
with open(LOG_FILE, 'r') as f:
    log = f.read()
# Extract values
samples = cmd_pattern.findall(log)
mapped = [int(x) for x in mapped_pattern.findall(log)]
processed = [int(x) for x in processed_pattern.findall(log)]

# Validate input sizes
assert len(samples) * 2 == len(mapped), "Unexpected number of mapping entries"
assert len(samples) * 2 == len(processed), "Unexpected number of processed entries"

# Summarize
records = []
for i in range(0, len(samples)):
    sample = samples[i]
    total_reads = mapped[i * 2]  # pig mapping
    pig_remaining = processed[i * 2]  # remaining after pig
    human_remaining = processed[i * 2 + 1]  # remaining after human
    
    pig_reads = total_reads - pig_remaining
    human_reads = pig_remaining - human_remaining
    non_host = human_remaining
    host_percentage = pig_reads / total_reads
    human_percentage =  human_reads / total_reads
    non_host_percentage = non_host / total_reads

    records.append({
        "Sample": sample,
        "Total_reads": total_reads,
        "Pig_reads": pig_reads,
        "Human_reads": human_reads,
        "Non_host_reads": non_host,
        "Host_contamination": host_percentage,
        "Human_percentage": human_percentage,
        "non_host_percent": non_host_percentage
    })

# Create and export table
df = pd.DataFrame(records)
df.to_csv(OUT_FILE, index=False)
print(f"Summary saved as {OUT_FILE}")
