from pathlib import Path
def iterate_fasta(filename_or_handle, keyFunc=None):
    """An iterator over a fasta file
    Yields tupples of key, sequence on each iteration
    """
    o = open_file(filename_or_handle)
    try:
        key = ""
        sequence = []
        for chunk in chunkify(o, b"\n>"):
            if chunk:
                key = chunk[: chunk.find(b"\n")].strip()
                if key.startswith(b">"):
                    key = key[1:]
                if keyFunc:
                    key = keyFunc(key)
                if chunk.find(b"\n") != -1:
                    seq = (
                        chunk[chunk.find(b"\n") + 1 :]
                        .replace(b"\r", b"")
                        .replace(b"\n", b"")
                    )
                else:
                    seq = ""
                yield (key.decode("utf-8"), seq.decode("utf-8"))
        return

        for line in o:
            if line.startswith(b">"):
                if key != "" and len(sequence) > 0:
                    yield (
                        key,
                        b"".join(sequence).replace(b"\n", b"").replace(b"\r", b""),
                    )
                key = line[1:].strip()
                sequence = []
            else:
                sequence.append(line)
        if key != "" and len(sequence) > 0:
            yield (
                key.decode("utf-8"),
                (b"".join(sequence).replace(b"\n", b"").replace(b"\r", b"")).decode(
                    "utf-8"
                ),
            )
    finally:
        o.close()
def open_file(fileNameOrHandle, mode="rb"):
    """Transparently open compressed or uncompressed files"""
    if hasattr(fileNameOrHandle, "read"):
        return fileNameOrHandle
    elif isinstance(fileNameOrHandle, Path):
        fileNameOrHandle = str(fileNameOrHandle)
    if fileNameOrHandle.endswith(".gz"):
        import gzip

        return gzip.GzipFile(fileNameOrHandle, mode)
    elif fileNameOrHandle.endswith(".bz2"):
        import bz2

        return bz2.BZ2File(fileNameOrHandle, mode)
    else:
        return open(fileNameOrHandle, mode)


def chunkify(handle, separator, block_size=None):
    """take a file handle and split it at separator, reading in efficently in 50 mb blocks or so"""
    if block_size is None:
        block_size = 50 * 1024 * 1024
    chunk = handle.read(block_size)
    chunk = chunk.split(separator)
    while True:
        for k in chunk[:-1]:
            yield k
        next = handle.read(block_size)
        if next:
            chunk = chunk[-1] + next
            chunk = chunk.split(separator)
        else:
            yield chunk[-1]
            break


def pathify(output_filename, default, create_parents=True):
    if output_filename is None:
        res = Path(default)
    else:
        res = Path(output_filename)

    if create_parents:
        res.parent.mkdir(exist_ok=True)
    return res
class GFF3:
    def escape(self, str):
        if str is None:
            return "."
        escape = '"\t\n\r=;'
        for k in escape:
            str = str.replace(k, "%" + "%2x" % ord(k))
        return str

    def format_attributes(self, attributes):
        if attributes is None:
            return "."
        if isinstance(attributes, dict):
            attributes = list(attributes.items())
        # valid_attributes = [
        #     "Name",
        #     "Alias",
        #     "Parent",
        #     "Target",
        #     "Gap",
        #     "Derives_from",
        #     "Note",
        #     "Dbxref",
        # ]
        res = []
        for id, value in attributes:
            #    if not id in valid_attributes and id != id.lower(): #lower case names are not reserved
            #       raise ValueError("Not a valid tag: %s" % id)
            res.append("%s=%s" % (self.escape(id), self.escape(value)))
        return ";".join(res)

    def dump_row(
        self,
        file_handle,
        seqid=None,
        source=None,
        type=None,
        start=None,
        end=None,
        score=None,
        strand=None,
        phase=None,
        attributes=None,
    ):
        file_handle.write(
            "\t".join(
                (
                    self.escape(seqid),
                    self.escape(source),
                    self.escape(type),
                    self.escape(start),
                    self.escape(end),
                    self.escape(score),
                    self.escape(strand),
                    self.escape(phase),
                    self.format_attributes(attributes),
                )
            )
        )
        file_handle.write("\n")
class GTF(GFF3):
    def escape(self, str):
        if str is None:
            return "."
        escape = '"\t\n\r";'
        for k in escape:
            str = str.replace(k, "%" + "%2x" % ord(k))
        return str

    def format_attributes(self, attributes):
        if attributes is None:
            return "."
        if isinstance(attributes, dict):
            attributes = list(attributes.items())
        # valid_attributes = [
        #     "Name",
        #     "Alias",
        #     "Parent",
        #     "Target",
        #     "Gap",
        #     "Derives_from",
        #     "Note",
        #     "Dbxref",
        # ]
        res = []
        for id, value in attributes:
            #    if not id in valid_attributes and id != id.lower(): #lower case names are not reserved
            #       raise ValueError("Not a valid tag: %s" % id)
            res.append('%s "%s"' % (self.escape(id), self.escape(value)))
        return "; ".join(res)

def fasta_to_gtf(
    input_filename,
    output_filename_or_handle,
    strand="+",
    source=None,
    name_mangler=None,
):
    """Mostly for creating a fake tf/gff for FileBasedgenomes"""
    return _fasta_to_gfx(
        GTF, input_filename, output_filename_or_handle, strand, source, name_mangler
    )


def fasta_to_gff(
    input_filename,
    output_filename_or_handle,
    strand="+",
    source=None,
    name_mangler=None,
):
    return _fasta_to_gfx(
        GFF3, input_filename, output_filename_or_handle, strand, source, name_mangler
    )


def _fasta_to_gfx(
    gfx,
    input_filename,
    output_filename_or_handle,
    strand="+",
    source=None,
    name_mangler=None,
):

    output_file_handle = open_file(output_filename_or_handle, "w")
    output_file_handle.write("##gff-version 3\n")
    gff = gfx()
    for org_name, seq in iterate_fasta(input_filename):
        if name_mangler:
            name = name_mangler(org_name)
            output_file_handle.write(f"##sequence-region {name} 1 {len(seq)}\n")

    for org_name, seq in iterate_fasta(input_filename):
        if name_mangler:
            name = name_mangler(org_name)
        else:
            name = org_name
        gff.dump_row(
            output_file_handle,
            name,
            source=source,
            type="gene",
            start=str(1),
            end=str(len(seq)),
            strand=strand,
            phase="0",
            attributes={
                # 'org_name': org_name,
                "gene_id": name,
                "gene_biotype": "protein_coding",
                "gene_name": org_name,
            },
        )
        gff.dump_row(
            output_file_handle,
            name,
            source=source,
            type="transcript",
            start=str(1),
            end=str(len(seq)),
            strand=strand,
            phase="0",
            attributes={
                # 'org_name': org_name,
                "transcript_id": "TR-" + name,
                "gene_id": name,
                "gene_name": org_name,
                "transcript_biotype": "protein_coding",
            },
        )
        gff.dump_row(
            output_file_handle,
            name,
            source=source,
            type="exon",
            start=str(1),
            end=str(len(seq)),
            strand=strand,
            phase="0",
            attributes={
                # 'org_name': org_name,
                "gene_id": name,
                "transcript_id": "TR-" + name,
                "exon_id": "EX-" + name,
            },
        )

    output_file_handle.close()


fasta_to_gff('genome.fasta','genome.gtf')
