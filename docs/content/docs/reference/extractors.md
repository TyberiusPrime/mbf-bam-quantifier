---
weight: 110
---

# Extractors

UMIS and barcodes can be extracted from read names, sequences or tags.


To do so you define an extractor.

Available extractors are


## Regular expression within read names 

Use a regular expression to extract from name

Example

```toml

[umi]
    # ...
    extract = { mode= 'RegexName', regex= 'umi=([ACTG]{8})' }
```

Extracts the first regex group. Follows [rust regex syntax](https://docs.rs/regex/latest/regex/).


## Search in name

Search for a fixed string in in read names,
then take the next N letters, optionally skipping some.

Example:

```toml
[umi]
    extract = {mode= 'SearchinName', search='umi|', length=8, skip=0}
```

(The extract sequence starts after `search`!).


## BAM tag

Read from a BAM tag. Must be a string tag.

```toml
[umi]
    extract = {mode = 'Tag',  tag='XQ'}
```


## ReadRegion
Extract from a fixed region within the read

```toml
[umi]
    extract = {mode = 'ReadRegion',  start=4, stop=8}
```

