queries = {
    "HWI-C00113:73:HVM2VBCXX:1:1213:4598:72277",
    "HWI-C00113:73:HVM2VBCXX:2:1103:2169:93070",
}

import pysam
import sys

input_filename = "umitools_input.bam"
output_filename = "umitools_broken.bam"

ih = pysam.Samfile(input_filename, "rb")
oh = pysam.Samfile(output_filename, "wb", template=ih)


for read in ih.fetch(until_eof=True):
    if read.query_name not in queries:
        continue
    oh.write(read)
oh.close()

pysam.index(output_filename)
