import pysam
import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]

ih = pysam.Samfile(input_filename, "rb")
oh = pysam.Samfile(output_filename, "wb", template=ih)


# umi tools has *weird* multi-alignment handling
# sometimes not deduplicating...
# sometimes ommitting primary reads, though the multi alignments are not in the same read pos,
# but actually outputting teh secondary alignment.

# there are some missed-deduplicaton samples in umitools_broken/

# anyway, the weird primary / secondary handling means we have to filter 
# to primary alignments here, before umi-tools sees them.
# and then we get the same numbers.

for read in ih.fetch(until_eof=True):
    umi = read.query_sequence[:6]
    read.set_tag("XX", umi, "Z")  # If no UMI, set XT to NA
    if not read.is_secondary:
        oh.write(read)
oh.close()

pysam.index(output_filename)
