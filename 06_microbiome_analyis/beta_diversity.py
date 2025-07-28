#!/usr/bin/env python3

################################################################
#Script adapted from bray_curtis.py by 2019 Jennifer Lu, jlu26@jhmi.edu,
#distributed as part of KrakenTools
################################################################

import os
import sys
import argparse
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="Compute Bray–Curtis dissimilarity")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input files (already filtered by rank)")
    parser.add_argument("-o", "--output", required=False, help="Optional output file to write Bray–Curtis matrix")
    return parser.parse_args()

def parse_abundance_file(filepath):
    abundances = {}
    with open(filepath, 'r') as f:
        for line in f:
            cols = line.strip().split("\t")
            if not cols or len(cols) < 9:
                continue
            taxid = cols[6]  # Column 7: taxon ID
            try:
                abundance = float(cols[-1])  # Last column = relative abundance
            except ValueError:
                continue
            if abundance > 0:
                abundances[taxid] = abundance
    return abundances

def compute_bray_curtis_matrix(samples):
    n = len(samples)
    matrix = np.zeros((n, n))
    totals = [sum(sample.values()) for sample in samples]

    for i in range(n):
        for j in range(i, n):
            shared = sum(min(samples[i].get(tid, 0), samples[j].get(tid, 0))
                         for tid in set(samples[i]) | set(samples[j]))
            denom = totals[i] + totals[j]
            bc = 1 - (2 * shared / denom) if denom > 0 else 0
            matrix[i][j] = matrix[j][i] = bc
    return matrix

def save_matrix(matrix, sample_names, filepath):
    with open(filepath, 'w') as out:
        out.write("x\t" + "\t".join(sample_names) + "\n")
        for i in range(len(sample_names)):
            row = [f"{matrix[i][j]:.3f}" for j in range(len(sample_names))]
            out.write(f"{sample_names[i]}\t" + "\t".join(row) + "\n")

def main():
    args = parse_args()
    filepaths = args.input
    sample_names = [os.path.basename(f).split('_')[0] for f in filepaths]
    sample_abundances = [parse_abundance_file(f) for f in filepaths]
    matrix = compute_bray_curtis_matrix(sample_abundances)

    print("# Sample indices:")
    for i, name in enumerate(sample_names):
        print(f"#{i}\t{name}")
    print("x\t" + "\t".join(str(i) for i in range(len(sample_names))))
    for i in range(len(sample_names)):
        row = [f"{matrix[i][j]:.3f}" if i <= j else "x.xxx" for j in range(len(sample_names))]
        print(f"{i}\t" + "\t".join(row))

    if args.output:
        save_matrix(matrix, sample_names, args.output)
        print(f"\nMatrix written to: {args.output}")

if __name__ == "__main__":
    main()