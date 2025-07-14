---
weight: 25
---


# UmiAllA

```toml
[[filter]]
    mode = 'umi_all_a'
    action = 'remove|keep_only' # see 'Filter steps' for details
```

Filter reads that have all 'A' bases in their UMI.


## Why is this here?
This is what STARsolo (2.7.11b) does. 

Mistakenly.

The doc and code suggests STAR should filter all homopolymers 
but there's a bug where the UMI length is not read correctly from 
`--soloUMIposition`, and that makes all homopolymer checks but A always fail, 
except if your UMI [length happens to be 10, or you're setting `--soloUMIlen` in addition](https://github.com/alexdobin/STAR/issues/1909#issuecomment-1648204555).

(It's claimed to be fixed in STAR 2.7.11a, but I've observed it even later).

You probably do not want to use this in production., but [UMI holopolymer](../umihomopolymer) instead.


