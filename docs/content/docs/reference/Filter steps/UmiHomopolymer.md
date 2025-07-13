---
weight: 15
title: UMI homopolymer
---


# UMIHomopolymer

```toml
[[filter]]
    mode = 'umi_homopolymer'
    action = 'remove|keep'
```

Filter reads that consist only of one base. E.g. 'AAAAA', 'TTTTT', 'GGGGG', 'CCCCC'.


