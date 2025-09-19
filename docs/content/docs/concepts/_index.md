---
weight: 10
bookFlatSection: true
title: "Concepts"
---
TODO



### The problem of same-gene-multi-mappers
(this section needs a home, just writing down thoughts for now)

Multi mappers are reads that are aligned to multiple regions in the genome
(as opposed to one read spanning regions).

They can either hit multiple genes (in which case, not counting them is 
essentially throwing away signal that says 'was from A or B'),
or they can hit the same gene repeatedly (in which case, counting them
is over inflating the expression of that gene).

UBC is a good candidate gene for the later, it's highly repetitive.


Our spiritual ancestor 'mbf-bam' collected read names per gene,
and then returned the number of unique names as the count.

That works around the problem, 
and could be implemented here as well, 
but given that a read may start in chunk X, 
but have a block in chunk X+N, 
we'd need to keep the lists around / shared between chunks,
and I feel that's going to tank performance fiercely.

Experimentation shows that we can use the PerRegion deduplication
and get *almost* the same results. Only difference is, I believe
that a read that's aligned twice, once so it matches a, and once where it matches a&b
would have been counted twice in mbf-bam, but only once here.






