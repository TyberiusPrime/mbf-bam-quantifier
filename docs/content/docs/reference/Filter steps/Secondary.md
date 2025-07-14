---
weight: 12
---

# Secondary / NonPrimary

```toml
[[filter]]
    mode = 'secondary'
    action = 'remove|keep_only' # see 'Filter steps' for details
```


Filter to primary (secondary & remove) or secondary (secondary & keep) alignments.

Based on the SAM flag 0x100 (0x100 = 256, 0x100 & SAM flag = 0 means primary alignment).

 
