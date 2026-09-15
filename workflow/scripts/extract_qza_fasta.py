"""Pull data/dna-sequences.fasta straight out of a QIIME2 artifact.

A .qza is a ZIP, and reading one in place is what this repo already does in
complete_filter.py, extract_humancount.py, prepare_diversity.py and
rename_qzv.py -- see commit 0c69404, "read qiime2 artifacts in place".

Doing it here too replaces a `qiime tools export` call that had to start the
container and load QIIME2 just to copy one file out. That call failed on
2026-09-15 (job 172048) after three and a half hours of pipeline -- an
expensive way to find out, and a whole class of failure this avoids.
"""
import shutil
import zipfile

src = snakemake.input[0]                      # noqa: F821
dst = snakemake.output[0]                     # noqa: F821

with open(snakemake.log[0], "w") as log:      # noqa: F821
    with zipfile.ZipFile(src) as z:
        names = [n for n in z.namelist()
                 if n.endswith("data/dna-sequences.fasta")]
        if len(names) != 1:
            raise RuntimeError(
                f"expected exactly one data/dna-sequences.fasta in {src}, "
                f"found {len(names)}: {names[:5]}")
        print(f"reading {names[0]} from {src}", file=log)
        with z.open(names[0]) as f, open(dst, "wb") as out:
            shutil.copyfileobj(f, out)

    with open(dst) as fh:
        n = sum(1 for line in fh if line.startswith(">"))
    print(f"{n} sequence(s) -> {dst}", file=log)
