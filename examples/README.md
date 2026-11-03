# Examples

Three data sets, the commands that assemble them, and the results to expect. Every command is run
from the repository root. The runner assembles all three into `examples/*/out/` (ignored by git)
and compares the contigs with the stored `expected/` files:

```bash
node examples/run-examples.js
```

The contigs must come out byte-identical whatever the machine or thread count; the reports and logs
differ only in the elapsed time and the number of threads. Timings below are from a four-core
machine with the default three worker threads.

## 1. Synthetic R10 reads of three amplicons (`synthetic/reads-r10.fastq.gz`)

300 reads: three random 2000 nt templates (`synthetic/templates.fasta`), 100 reads each, with an
R10-like error profile (about 2 % error, insertions and deletions concentrated in homopolymer
runs), read ends jittered by up to 20 nt, either strand at random, and a `barcodeNN` tag in every
header. Made by `synthetic/make-reads.js`:

```bash
node examples/synthetic/make-reads.js reads-r10.fastq --amplicons 3 --depth 100 --length 2000 --model R10 --seed 11 --templates templates.fasta
```

Assembled with the defaults, grouped by the barcode tag:

```bash
node tools/cap-cli.js examples/synthetic/reads-r10.fastq.gz --out examples/synthetic/out/r10 --bin-regex "barcode(\d+)"
```

Expected (`synthetic/expected/r10.bins.tsv`): every barcode gives **one contig of exactly 2000 bp,
all 100 reads placed, no singleton**, and each contig is identical to its template (edit distance
0, checked by `synthetic/check.js`). About 25 s in all, 220 MB peak memory.

```
group  reads_seen  reads_used  contigs  longest_bp  total_consensus_bp  singletons
01     100         100         1        2000        2000                0
02     100         100         1        2000        2000                0
03     100         100         1        2000        2000                0
```

The three reports show the typical anatomy of a deep amplicon assembly: about 4950 candidate pairs
per barcode, all but one accepted, a substitution error rate of about 0.6 % measured from the
mismatches in the run's own overlaps, and no polymorphic site, since the reads come from one template.

## 2. The same amplicons with R9 reads (`synthetic/reads-r9.fastq.gz`)

The same templates and seed with the R9-like profile (about 7 % error). With the default identity
of 90 % every overlap is refused as too divergent; with `--min-identity 85`, the setting the help
recommends for R9 flow cells, every read is placed:

```bash
node tools/cap-cli.js examples/synthetic/reads-r9.fastq.gz --out examples/synthetic/out/r9 --bin-regex "barcode(\d+)" --min-identity 85
```

Expected (`synthetic/expected/r9.bins.tsv`): one contig per barcode, 1995, 1986 and 1999 bp, all
reads placed. The contigs differ from their templates by 5, 14 and 1 edits (0.25 %, 0.70 % and
0.05 %), all of them in homopolymer runs that most reads under-call in the same way, which is the
expected limit of R9 data. About 20–30 s.

To check any contigs file against the templates:

```bash
node examples/synthetic/check.js examples/synthetic/out/r9.02.contigs.fasta
```

## 3. A real Oxford Nanopore run file (`nanopore/FAS90246_pass_barcode79_cfafe70b_765ba3f3_1.fastq.gz`)

One of the files MinKNOW writes per barcode: 4000 reads from an R9.4.1 flow cell, basecalled
with the high-accuracy model, barcode 79 of a pooled amplicon library. The reads are what a real
run gives: 1 to 14,687 nt long (median 381 nt; 3119 of the 4000 are 200 nt or longer), many
different products in one barcode, most of them at low depth.

```bash
node tools/cap-cli.js examples/nanopore/FAS90246_pass_barcode79_cfafe70b_765ba3f3_1.fastq.gz --out examples/nanopore/out/barcode79 --bin-regex "barcode(\d+)" --min-identity 85 --min-length 200 --no-singletons
```

Reads shorter than 200 nt are dropped, the first 1000 of the rest are kept (the default
`--max-reads`), and the singletons are left out of the FASTA to keep it small (they are still listed
in the report). Expected (`nanopore/expected/barcode79.bins.tsv` and the report's run summary):

| | |
|---|---|
| Reads seen / used | 3119 / 1000 |
| Reads placed in contigs | 128 |
| Singletons | 872 |
| Contigs | 47 |
| Longest contig | 11,858 bp |
| N50 | 4062 bp |
| Candidate pairs aligned | 7888 |
| Overlaps accepted / rejected | 103 / 7785 |
| Time | about 10 s on three threads, 260 MB |

This example is included for what it shows about real data rather than for a tidy result. Most
reads stay singletons because the barcode holds many products and the sample of 1000 reads
rarely contains two reads of the same one; the rejected-overlaps table of
`nanopore/expected/barcode79.79.report.txt` lists why each candidate pair was refused, mostly
identity below 85 % between two noisy R9 reads and unmatched overhangs where reads share a primer
or a repeat but not a product. On the same file, lowering `--min-identity` to 75 roughly triples
the number of reads placed, and using all 3119 reads instead of 1000 (`--max-reads 4000`) gives
about 200 contigs with the longest 13–17 kb at 1–6× coverage, in about 36 s. For a run of this
kind the useful step is to group the reads by product (an amplicon tag added by a classifier,
matched with `--bin-regex`) rather than by barcode alone.

## Making your own test data

`synthetic/make-reads.js` writes any number of amplicons, reads per amplicon, read length and
error model, with wrapped records (`--wrap 80`) or a file of a given size (`--bytes`) for load
tests:

```bash
node examples/synthetic/make-reads.js big.fastq --amplicons 3 --depth 1000 --length 2000 --model R10 --bytes 2000000000 --templates big-templates.fasta
```

A 2 GB file made this way (529,235 reads) streamed and assembled, 1000 reads per barcode, in
483 s at 904 MB peak memory, every contig within two edits of its template.
