---
weight: 101
---
# Cell Barcode Section

```toml
[cell_barcodes]
    extract = # see Extractors section

    separator_char = '_' # if multiple/joined barcodes, they're separated with this.
    max_hamming = 0
    allowlist_files = ['file1.txt', 'file2.txt'] # list of files with barcodes, one per line
```

## 

Activate per cell (barcode) counting for single cell applications.

Optionally, filter to barcodes that are within max_hamming distance to an entry
entry in the allow list.

Supports combinatorial barcodes separated by `separator_char`.
In this case, you need as many allowlist files as there are distinct barcodes sections,
and we accept any combination.

This changes the output into matrix.mtx format.
