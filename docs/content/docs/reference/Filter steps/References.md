---
weight: 10
---

# Filter References

```
[[filter]]
    mode = 'references'
    references = ['chr1', 'chrX'] # at least one
    action = 'remove' # or keep to filter all other bam references
```

Filters reads depending on their BAM reference.

Regions filtered by this filter will not have their BAM records fetched,
making this a great time saver if you only need some references.

