#!/usr/bin/env python3

import sys

def gtf_to_bed_coords_with_type(gtf_file):
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            chrom = fields[0]
            start = int(fields[3]) - 1  # 0-based start
            end = int(fields[4])        # 1-based inclusive end
            feature_type = fields[2]
            print(f"{chrom}\t{start}\t{end}\t{feature_type}")

if __name__ == '__main__':
    gtf_to_bed_coords_with_type(sys.argv[1])

