- smart exon tag counting like we do...

- star solo comparison tests

- coverage quantification

- faster interval format instead of gtf parsing?

- can we split our intervals bgfz borders? 

- extract for umi / barcodes

- unify from stranded/unstraded code bases?


- rewrite central engine
    - from the chunked iterator, pass each read to a trait
    impl that calculates two HashMaps, 
    <id, 1>, and <id, 1> for matching in direction and reverse matching.
    Filter then set the counts to 0.
    Unstranded adds both to the first.
    Stranded leaves them as is.
    Funky 'overlap' rules manipulate them in place.
    then, if requested write this as a tag on the read 
    for visibility/debugging purposes.



