import sys
import gzip
import os
import subprocess

gtf_file = sys.argv[1]
bam_file = sys.argv[2]
query = sys.argv[3]

# Open GTF file (gz or plain text)
def open_file(filename):
    return gzip.open(filename, 'rt') if filename.endswith('.gz') else open(filename, 'r')

# First pass: collect regions matching the query
regions = []
with open_file(gtf_file) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if query in line:
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            regions.append((chrom, start, end))

# Second pass: write overlapping lines to test.gtf
def overlaps(c1, s1, e1, c2, s2, e2):
    return c1 == c2 and not (e1 < s2 or e2 < s1)

with open_file(gtf_file) as f_in, open('test.gtf', 'w') as f_out:
    for line in f_in:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        for rchrom, rstart, rend in regions:
            if overlaps(chrom, start, end, rchrom, rstart, rend):
                f_out.write(line)
                break

# Filter BAM file using samtools for all regions
region_strs = [f"{chrom}:{start}-{end}" for chrom, start, end in regions]
cmd = ['samtools', 'view', '-M', '-b', bam_file] + region_strs
with open('test.unsorted.bam', 'wb') as bam_out:
    subprocess.check_call(cmd, stdout=bam_out)
subprocess.check_call(['samtools', 'sort', 'test.unsorted.bam'], stdout=open("test.sorted.bam", 'wb'))
# Index the sorted BAM file
# move to final name
os.rename('test.sorted.bam', 'test.bam')
subprocess.check_call(['samtools', 'index', 'test.bam'])
os.remove('test.unsorted.bam')

