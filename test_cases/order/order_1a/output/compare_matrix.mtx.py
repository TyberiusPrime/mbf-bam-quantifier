#!/usr/bin/env python

import sys
from pathlib import Path

should = Path(sys.argv[1])
actual = Path(sys.argv[2])


def read_matrix(fn):
    barcode_file = fn.with_name("barcodes.tsv")
    feature_file = fn.with_name("features.tsv")
    # Read barcodes
    with open(barcode_file) as f:
        barcodes = [line.strip() for line in f]
    with open(feature_file) as f:
        features = [line.strip().split("\t")[0] for line in f]
    print(len(barcodes), "barcodes")
    print(len(features), "features")
    result = {}
    first = True
    with open(fn) as f:
        for line in f:
            if line.startswith("%"):
                continue
            if first:
                # Skip the first line, which row count, column count, entries...
                parts = line.strip().split(" ")
                if int(parts[0]) != len(features):
                    raise ValueError(
                        f"Expected {len(features)} features, got {parts[0]} in {fn}"
                    )
                if int(parts[1]) != len(barcodes):
                    raise ValueError(
                        f"Expected {len(barcodes)} barcodes, got {parts[1]} in {fn}"
                    )
                no_of_expected_entries = int(parts[2])
                first = False
                continue
            parts = line.strip().split(" ")
            feature_idx = int(parts[0]) - 1
            barcode_idx = int(parts[1]) - 1
            value = int(parts[2])
            try:
                barcode = barcodes[barcode_idx]
                feature = features[feature_idx]
            except IndexError:
                print(
                    f"Index error: 'feature_idx'={feature_idx}, barcode_idx={barcode_idx} in {fn}, no features: {len(features)}, # barcodes {len(barcodes)}"
                )
                raise
            result[(feature, barcode)] = value
    if len(result) != no_of_expected_entries:
        raise ValueError(
            f"Expected {no_of_expected_entries} entries, got {len(result)} in {fn}"
        )
    return result


print("loading should")
a = read_matrix(should)
print("loading actual")
b = read_matrix(actual)

if a != b:
    print("Matrix differed")
    print(len(a), len(b))
    union = set(a.keys()).union(set(b.keys()))
    added_count = 0
    added_sum = 0
    lost_count = 0
    differ_count = 0
    for k in sorted(union):
        ...
        if k in a and k in b:
            if a[k] != b[k]:
                print(f"{k}: {a.get(k, 0)} != {b.get(k, 0)}")
                differ_count += 1
        elif k in a:
            print("lost", k, a[k])
            lost_count += 1
        elif k in b:
            print("added", k, b[k])
            added_sum += b[k]
            added_count += 1
        else:
            print("wtf")
    print("added", added_count, 'added_sum', added_sum)
    print("lost", lost_count)
    print("differ", differ_count)

    sys.exit(1)
else:
    print("matrix identical", len(a), len(b))

sys.exit(0)
