"""Split a FASTA into N chunks of whole records (snakemake script).

Why: filter_human decides per SEQUENCE, independently of every other sequence
and of which sample it came from -- the dereplication has already merged the
cohort into one feature set. So the BLAST search parallelises exactly over
chunks of sequences, while the cohort (and with it the shared ASV table) stays
intact. Records are dealt out round-robin so every chunk gets a comparable mix
of easy and hard queries, which keeps the parallel runtimes even.
"""
import sys

src = snakemake.input[0]                      # noqa: F821
outs = list(snakemake.output)                 # noqa: F821
log = open(snakemake.log[0], "w")             # noqa: F821

handles = [open(p, "w") for p in outs]
n = len(handles)
idx = -1
cur = None
with open(src) as fh:
    for line in fh:
        if line.startswith(">"):
            idx += 1
            cur = handles[idx % n]
        if cur is None:
            continue                          # stray leading text, ignore
        cur.write(line)
for h in handles:
    h.close()
print(f"{idx + 1} sequences -> {n} chunks", file=log)
log.close()
