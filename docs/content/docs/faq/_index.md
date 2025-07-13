+++
title = "FAQ"
BookFlatSection = true
+++


## Why are there so few defaults?

mbf-bam-quantifier is following the python mantra 'explicit is better than implicit'.

It's presumptuous to assume our user's use case, and mismatches between an assumed
and actual use case lead to unwelcome surprises that the user might only discover 
much later if at all.

Defaults also make for difficult upgrade paths - you can't really change them later on
without silently breaking your user's outputs (They'll be different, but it will be 
not immediately clear to the user why).

Instead have a look at the how-to section.


