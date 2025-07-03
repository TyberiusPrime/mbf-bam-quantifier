import pysam


a = pysam.AlignmentFile("actual/output/annotated.bam")
total = 0
error = 0
minimal_block_len_in_aggreement = 1000
longest_primary = 0
longest_secondary = 0
for read in a.fetch(until_eof=True):
    total += 1
    try:
        gn = read.get_tag("GN")
        xq = read.get_tag("XQ")
        if "=" in xq:
            xq = xq.split("=")[0]
        if gn != xq:
            error += 1
            # how many bases does the read actually have aligned
            # ll = len(list(read.get_aligned_pairs(matches_only=True)))
            # if (ll <= 36) and (gn == "-"):
            #     continue
            ll = 0
            for block in read.get_blocks():
                bl = block[1] - block[0]
                ll = max(ll, bl)
            if read.is_secondary and ll > longest_secondary:
                longest_secondary = ll
            elif not read.is_secondary and ll > longest_primary:
                longest_primary = ll
            print(read.query_name, read.cigarstring, 'secondary=', read.is_secondary)
            print(ll)
            print(gn)
            print(xq)
            print(read.get_tag("CB"))
            print("----------------------")
        else:
            ll = 0
            for block in read.get_blocks():
                bl = block[1] - block[0]
                ll = max(ll, bl)
            if ll == 12:
                print('ll', read.qname)
                print(read.is_secondary)
            minimal_block_len_in_aggreement = min(minimal_block_len_in_aggreement, ll)

    except KeyError:
        pass

print("error", error)
print("total", total)
print("minimal_block_len_in_aggreement", minimal_block_len_in_aggreement)
print("longest_primary", longest_primary)
print("longest_secondary", longest_secondary)
# Ok, so, where we assign, but starsolo doesn,t have aligned blocks <= 52


#es sind  87 reads... das ist jetz nicht was mir fehlt
