"""Collect human hit IDs from the BLAST chunk outputs (snakemake script).

This reproduces q2_quality_control._blast._extract_hits / _perc_coverage
VERBATIM -- same coverage formula, same '*' guard, same >= comparison -- so the
chunked search yields byte-identical decisions to `qiime quality-control
exclude-seqs`. Because _extract_hits is a set union over independent lines, it
does not matter whether the lines come from one BLAST run or from N concatenated
chunk files.

Output is a QIIME2 metadata TSV of the hit IDs, ready for
`qiime feature-table filter-seqs --m-metadata-file` (with and without
--p-exclude-ids).
"""


def _perc_coverage(end, start, query_len):
    # query start is one-based relative to start of sequence, hence the +1.
    return (float(end) - float(start) + 1) / float(query_len)


chunks = list(snakemake.input)                            # noqa: F821
out = snakemake.output[0]                                 # noqa: F821
min_cov = float(snakemake.params.perc_query_aligned)      # noqa: F821
log = open(snakemake.log[0], "w")                         # noqa: F821

hits = set()
lines = 0
for path in chunks:
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            lines += 1
            query_id, subject_id, query_len, start, end = line.split("\t")
            if _perc_coverage(end, start, query_len) >= min_cov \
                    and subject_id != "*":
                hits.add(query_id)

with open(out, "w") as o:
    o.write("id\n")
    for h in sorted(hits):
        o.write(f"{h}\n")

print(f"{len(chunks)} chunk file(s), {lines} alignment(s), "
      f"{len(hits)} sequence(s) above perc_query_aligned={min_cov}", file=log)
log.close()
