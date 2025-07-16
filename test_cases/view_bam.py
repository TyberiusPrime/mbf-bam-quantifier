#!/usr/bin/env python
import sys
import pysam


def print_reads(reads):
    max_qname_len = max((len(read.qname) for read in reads), default=0)
    for read in reads:
        try :
            cor_pos = read.get_tag('XP') -1
        except:
            if read.is_reverse:
                cor_pos = read.pos  + read.infer_query_length()
            else:
                cor_pos = read.pos
        if not read.has_tag('XF'):
            print(
        read.qname.ljust(max_qname_len),
        read.flag,
        read.reference_name,
        read.pos + 1,
        cor_pos + 1,
        sep="\t",
    )


try:
    try:
        f = pysam.Samfile(sys.argv[1])
        print_reads(list(f.fetch(sys.argv[2])))
    except ValueError as e:
        if "fetching by region is not available for SAM files" in str(e):
            print_reads(list(f.fetch(until_eof=True)))
except BrokenPipeError:
    pass
