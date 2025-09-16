---
title: Read1 / Read2
weight: 12
---

# Read1 / Read2

```toml
[[filter]]
    mode = 'read1'
    action = 'remove|keep_only' # see 'Filter steps' for details
```

```toml
[[filter]]
    mode = 'read1'
    action = 'remove|keep_only' # see 'Filter steps' for details
```

Filter to reads that have either the Read1 (0x20) or Read2 (0x40) SAM flag.

Note that while these are labeled 'read1/read2' they're defined as 'first segment in the template'
and 'last segment in the template'. 
