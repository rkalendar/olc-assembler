# OLC Assembler

An overlap–layout–consensus (OLC) assembler for Sanger fragments and Oxford Nanopore amplicon
reads, written in JavaScript with no dependencies. The engine (`js/cap.js`) follows CAP3
(Huang & Madan, *Genome Research* 9:868, 1999) and runs unchanged in the browser, as the
[OLC Assembler page of DigitalGens](https://digitalgens.org/cap.html), and under Node.js. This
repository holds the engine and its **command-line batch mode**, `tools/cap-cli.js`, which assembles
whole Nanopore sequencing runs: it streams FASTQ files of any size, groups the reads by a barcode or
amplicon tag in the header, keeps a sample of each group and assembles the groups one after another
on worker threads, writing the same contigs, report and log for each group that the web page shows.

```bash
node tools/cap-cli.js fastq_pass/barcode79 --out runs/barcode79 --bin-regex "barcode(\d+)" --min-identity 85
```

## Contents

- [What it does](#what-it-does)
- [Requirements and installation](#requirements-and-installation)
- [Quick start](#quick-start)
- [Input](#input)
- [Grouping and sampling reads](#grouping-and-sampling-reads)
- [Options](#options)
- [Output files](#output-files)
- [Settings for Oxford Nanopore data](#settings-for-oxford-nanopore-data)
- [Performance and reproducibility](#performance-and-reproducibility)
- [Examples](#examples)
- [The browser version](#the-browser-version)
- [Repository layout](#repository-layout)
- [Citation](#citation)
- [Licence](#licence)

## What it does

The engine assembles one or more contigs from a set of fragments that may be repeated, imperfect
and in either orientation: Sanger reads, cloned inserts, NGS contigs, synthetic oligonucleotides, or
Oxford Nanopore R9/R10 reads of PCR products. The stages are those of CAP3:

1. **Seeding.** Every read is indexed by its *k*-mers (*k* = 12). Two reads sharing at least three
   seed words in a consistent diagonal become a candidate pair. Words seen in more reads than a
   repeat threshold derived from the run itself are treated as repeat words and left out of seeding,
   except that in deep amplicon data, where every word is deep, the depth typical of the run is
   recognised as coverage rather than repetition and the deep words are seeded.
2. **Banded overlap alignment** of each candidate pair (match 2, mismatch −3, gap open −5, gap
   extend −2; band ±32, widened to ±256 when the alignment runs out of it). Reverse-complement
   overlaps are searched unless `--no-rc`.
3. **Overlap acceptance.** An alignment is accepted as an overlap when it is at least
   `--min-overlap` nt long, at least `--min-identity` % identical and leaves no more than
   `--max-overhang` nt unmatched at the ends of the reads. The overhang filter is CAP3's defence
   against interspersed repeats: two reads that merely share a repeat align over it and then
   diverge, leaving unmatched tails on both reads, whereas in a genuine join one read runs off each
   end of the alignment. Two further guards refuse overlaps that lie almost entirely inside repeated
   sequence (measured by word counts across the input) and joins across a repeat proven by
   divergence (two reads that align and then go separate ways, each way confirmed by other reads).
4. **Chimera detection.** Three checks. Before seeding, a read whose end is the reverse complement
   of sequence earlier in the same read — a fold-back read, the strand read on into its own
   complement — has the duplicated stretch cut away (the copy must be at least 50 nt long and reach
   an end of the read; a hairpin keeps its longer arm, copies that lie apart keep only what lies
   between them) and is listed in the report; left whole, its tail is laid beyond the end of the
   product and written into the consensus. At overlap acceptance, an alignment whose identity falls
   over any 50 columns to 35 points below `--min-identity` is refused as running across a chimeric
   junction: two reads that each join the same two products agree on either side of the stretch
   between their junctions and not across it, however high the identity of the whole alignment.
   After the overlaps, a read carries two molecules (two products ligated end to end, or two
   molecules read through one adapter) when, at one point inside it, the reads overlapping one part
   stop and those overlapping the other part start, or when the reads that align with it part from
   it there — they agree with it up to that point and then go on with sequence of their own, so the
   overhang filter refused the overlap — and almost none run across; such reads are set aside and
   listed in the report with the pattern seen. `--no-chimera-check` turns all three off.
5. **Greedy layout.** Accepted overlaps are taken in order of score and chained into contigs; an
   overlap that contradicts the layout built so far is recorded as a layout conflict, and a read end
   that overlaps several reads which do not overlap each other stops the contig rather than guessing.
   The overlaps are sorted deterministically, so the layout does not depend on how many threads
   computed them.
6. **Column-wise multiple alignment and consensus.** Each read is aligned under its contig, with a
   cheaper gap than at the overlap stage (gap open −2). At the overlap stage's price a read with one
   base more than the growing consensus in a homopolymer and a deletion a few bases away is aligned
   straight across the two as mismatches, its extra base lands outside the run, and a run that came
   out short stays short whatever the depth. The first base of a gap stays dearer than a mismatch,
   or two neighbouring SNPs would be written as an insertion and a deletion. The column vote on a
   run still favours the length the reads were fitted to, so runs of three or more are then checked
   against the reads themselves: each read's run is counted in its own sequence between anchors just
   outside the run, and the run takes the longest length held by at least three-quarters as many
   reads as the commonest length, and by two reads at least. The check only lengthens, because reads
   lose bases in runs far more often than they gain them. `--set consGapOpen=-5` restores the old
   price, `--set runVoteMin=0` turns the check off, and the report counts the runs lengthened. A
   column is called polymorphic, and given the IUPAC code of its alleles, only when the minority base is
   carried by at least two reads, makes up at least 25 % of a column covered by at least three reads,
   and passes a binomial test against the sequencing error rate measured from the data itself. An
   evenly split column is a tie; it goes to the base that lies further inside its read, because read
   ends carry the poorest base calls and, in PCR products, the primer. Indel polymorphisms are listed
   in the report and the consensus takes the majority. A contig end carried by fewer than three
   reads, or by fewer than 5 % of the contig's median coverage, over more than `--max-overhang`
   bases is cut back to the last well-supported column (`--set minEndCov=0` keeps it); the bases
   cut from each end are given per contig in the report. A contig whose own median coverage is
   under twice that threshold is left alone: at one or two reads of coverage every end is carried
   by one read.

The design is explained stage by stage in the comments of `js/cap.js`, and the parameters are
documented, with the reasoning behind each default, in `CAP_DEFAULTS` at the top of that file.

Three things happen for Nanopore amplicon data that a Sanger assembly never needs: deep coverage is
told apart from repetition by the word counts typical of the run; where neighbouring PCR products
overlap (tiled amplicons) the doubly covered junction is recognised as a tiling junction and the
join is accepted; and reads carrying two molecules are set aside instead of splicing products
together.

## Requirements and installation

- [Node.js](https://nodejs.org) 18 or later. Nothing else: the code has no dependencies and uses
  only the `fs`, `zlib`, `readline` and `worker_threads` modules that ship with Node.
- Any operating system Node runs on. The commands below are the same in bash, PowerShell and cmd.

```bash
git clone https://github.com/rkalendar/olc-assembler.git
cd olc-assembler
node tools/cap-cli.js --help
```

The three files that matter are `tools/cap-cli.js`, `js/cap.js` and `js/cap-report.js`; copying
them with that layout is a complete installation. `npm link` (or `npm install -g .`) in the
repository folder makes the command available as `olc-assembler` from any directory.

## Quick start

One FASTQ file, one group, default parameters (R10 reads):

```bash
node tools/cap-cli.js reads.fastq.gz --out runs/sample1
```

A MinKNOW output folder, reads grouped by the barcode in their headers, R9 chemistry, reads shorter
than 200 nt dropped and a uniform sample of 1000 reads per barcode:

```bash
node tools/cap-cli.js fastq_pass/barcode79 --out runs/barcode79 --bin-regex "barcode(\d+)" --min-identity 85 --min-length 200 --sample 1
```

Several files, with the read-under-consensus alignment and the read placements written too:

```bash
node tools/cap-cli.js part_1.fastq.gz part_2.fastq.gz part_3.fastq.gz --out runs/all --alignment --placements
```

The output folder is created if it does not exist. Progress goes to stderr; the final line names
the `bins.tsv` summary. Run with `--help` for the complete option list and `--version` for the
version.

## Input

- **FASTA or FASTQ**, plain or gzip-compressed (`.gz`), any size: the file is streamed record by
  record and never held in memory as a whole. A 2 GB FASTQ of 529,235 reads streamed in about eight
  minutes on four cores at 904 MB peak memory, assembling 1000 reads of each of its three barcodes.
- **Several files or a folder.** Every positional argument is a file or a folder; a folder stands
  for every `.fastq`, `.fq`, `.fasta`, `.fa`, `.fna` or `.fas` file inside it, plain or `.gz`, in
  natural name order (`part_1`, `part_2`, … `part_10`), which is how a basecaller leaves a run. The
  files are read one after another into the same groups.
- **FASTQ layout.** A record is `@name`, the sequence on one or more lines, a `+` line and a quality
  string of the same length as the sequence. Wrapped and indented records and a byte-order mark are
  read correctly; the quality string may start with `@`, `+` or `>`. A record whose quality outruns
  its sequence is a read error (exit code 2); a last record cut off at the end of a file is dropped
  with a warning. The quality values are not used by the assembly.
- **Read names** are the first word of the header; the whole header line is what `--bin-regex` is
  matched against.
- `--min-length` and `--max-length` drop reads outside a length range before grouping, so that the
  per-group sample is drawn from reads worth assembling.

## Grouping and sampling reads

A sequencing run usually holds several samples or several PCR products, and an all-against-all
assembler should not align reads that cannot belong together. `--bin-regex RE` groups the reads by
a regular expression matched against the header line: the first capture group (or the whole match
when there is none) names the group. `"barcode(\d+)"` groups a MinKNOW run by barcode;
`"^(\S+)_"` groups by the part of the name before the first underscore; any header field a
demultiplexer or an amplicon classifier adds can be used the same way. Reads that the expression
does not match are counted and left out, or assembled as a group of their own with
`--keep-unbinned`. Without `--bin-regex` all reads form one group.

Of each group, `--max-reads` reads are kept (default 1000; the engine accepts up to 20,000): the
first N as they come, or a uniform random sample with `--sample [SEED]` (reservoir sampling, so a
sample of a 2 GB file costs no more memory than the sample itself; the seed makes it reproducible).
The consensus of one PCR product stops improving after a few hundred reads, and the pairwise stage
grows with the square of the number of reads, so a few hundred to a thousand reads per product is
the useful range. Only one group is in memory at a time.

## Options

```
Input and grouping
  --out PREFIX          Output prefix (default: the input file name without its extension,
                        in the current folder). Files: PREFIX[.BIN].contigs.fasta,
                        .report.txt, .log.txt and PREFIX.bins.tsv.
  --bin-regex RE        Group reads by RE matched against the whole header line; the first
                        capture group (or the whole match) names the group.
  --keep-unbinned       Assemble reads the regex does not match as group "unbinned"
                        (default: count them and leave them out).
  --max-reads N         Reads kept per group (default 1000; the engine accepts up to 20000).
  --sample [SEED]       Keep a uniform random sample of N reads per group instead of the
                        first N (reservoir sampling; SEED makes it reproducible, default 1).
  --min-length N        Drop reads shorter than N nt before grouping (default 0).
  --max-length N        Drop reads longer than N nt before grouping (default: none).

Assembly
  --threads K           Worker threads for the overlap stage (default: CPUs minus one, at
                        most 16). 1 runs everything on this thread; the output is identical.
  --min-identity P      Minimum overlap identity, % (default 90; 85 suits R9 reads).
  --min-overlap N       Minimum overlap length, nt (default 40).
  --max-overhang N      Maximum unmatched overhang, nt (default 30).
  --kmer K              Seed word length (default 12).
  --min-read-len N      Reads shorter than this are not used (default 40).
  --end-trim N          Bases clipped from both ends of every read (default 0).
  --end-penalty N       End-gap penalty (default 0).
  --no-rc               Do not consider reverse-complement overlaps.
  --no-chimera-check    Keep reads that look like two molecules joined in one read: no
                        fold-back trimming, no junction test, no chimeras set aside.
  --no-repeat-guard     Do not refuse joins across a repeat proven by divergence.
  --no-singletons       Leave unassembled reads out of the FASTA.
  --set KEY=VALUE       Any other engine parameter by name (repeatable), e.g. --set minSeeds=2.

Output
  --alignment           Also write PREFIX[.BIN].alignment.txt (reads under the consensus).
  --placements          Also write PREFIX[.BIN].placements.tsv.
  --quiet               No progress on stderr, only the final summary.
  --help, --version
```

The assembly defaults are the defaults of the web page, so a group assembled here gives the same
result as the same reads pasted into the page. Every engine parameter not listed above can be set
by its name in `CAP_DEFAULTS` with `--set`.

Exit codes: 0 done; 1 bad arguments; 2 the input could not be read; 3 at least one group failed to
assemble (the others are still written, and `bins.tsv` says which failed).

## Output files

For `--out runs/sample1` and a group named `79`:

| File | Content |
|---|---|
| `runs/sample1.79.contigs.fasta` | The consensus sequences, upper case, one record per contig. The header carries `length`, `reads`, `cov` (mean coverage), `snp` and `indel` (polymorphic sites called) and `tie` (evenly split columns, carried as IUPAC codes). Unassembled reads follow, marked as singletons, unless `--no-singletons`. |
| `runs/sample1.79.report.txt` | The assembly report: a run summary (reads used, placed, singletons, chimeric and fold-back reads, contigs, total consensus, longest contig, N50, mean coverage, candidate pairs, overlaps accepted and rejected, layout conflicts, polymorphic sites, contig ends trimmed, elapsed time and threads), the parameters used, then tab-separated tables: one row per contig (length, reads, coverage, variants, strands, bases trimmed from each end); one per read placement (contig, read, strand, start, end, length, clipped bases); one per polymorphic site or tie (position, coverage, counts of A, C, G, T and gap, the IUPAC code, minor-allele fraction); one per accepted overlap (reads, strand, offset, length, identity, score, mismatches, gaps, containment and whether it joined, confirmed or conflicted with the layout); one per rejected overlap with its reason; one per layout conflict; the singletons with their length, GC content and repeat-word hits; the chimeric reads with their junctions; and the fold-back reads with their two copies and the part kept. The tables paste straight into a spreadsheet. |
| `runs/sample1.79.log.txt` | The short log: input, seeding statistics, the measured error rate, threads used and any notices. |
| `runs/sample1.79.placements.tsv` | With `--placements`: the read placement table alone, for downstream scripts. |
| `runs/sample1.79.alignment.txt` | With `--alignment`: every read under its contig's consensus, a dot where the read agrees, `-` for a gap, the strand after the read name and `*` over polymorphic columns. |
| `runs/sample1.bins.tsv` | One row per group: `group`, `reads_seen`, `reads_used`, `contigs`, `longest_bp`, `total_consensus_bp`, `singletons`, `threads`, `seconds`, `status`. |

Without `--bin-regex` the group segment of the file names is omitted (`runs/sample1.contigs.fasta`).
Contig orientation is arbitrary: with no reference a contig may come out reverse-complemented, and
an `R` in one orientation is a `Y` in the other.

The report is the first thing to read. The counts of rejected overlaps and layout conflicts, with
their reasons, show where the data are ambiguous: a cluster of rejections mentioning repeat content
or unmatched overhang is the boundary of a repeat, and a cluster of overlaps refused as too
divergent means `--min-identity` is above what the reads can reach.

## Settings for Oxford Nanopore data

- **R10 reads** (about 1–3 % error) assemble with the defaults.
- **R9 reads** (5–10 % error) sit close to the identity limit: when two reads each carry 5–8 %
  error their alignment falls below 90 % identity, every overlap is refused as too divergent and
  every read is a singleton. Use `--min-identity 85`. In simulations this placed every R9 read while
  leaving the R10 result unchanged, and going lower gained nothing on single products.
- **Short reads.** A run contains fragments, adapter dimers and broken reads; `--min-length 200`
  (or the length of the shortest product) keeps them out of the sample.
- **Depth.** A few hundred reads per product give the same consensus as thousands; `--max-reads` and
  `--sample` keep the run time proportional to the number of products rather than to the run.
- **What to expect.** On simulated amplicons of 1–2 kb with R10-like and R9-like error profiles,
  5 to 500 reads per product (sets S1 and S2 of `benchmark/`), R10-like reads gave the template
  exactly from 50 reads, and R9-like reads gave at most 0.2 % error from 50 reads and the template
  exactly at 200 and 500. Homopolymer runs of four to seven bases came out at full length (2 of 573
  lost, from 50 reads on). A longer run can stay a base short, because few reads carry the whole of
  it: 28 of 67 runs of eight or more were, and a stress set with 20 inserted runs of 4–10 nt per
  2 kb gave 0.16 % at 200 reads and 0.30 % at 50. Below 20 reads, random errors remain: about 0.3 %
  at 10 reads. Real reads differ from any simulation, so check the alignment at homopolymer runs.
- **Chimeric reads.** A few percent of the reads of a PCR library carry two molecules. Fold-back
  reads have their duplicated stretch cut away, overlaps across a junction are refused and ligation
  chimeras are set aside, so the products come out at their own length; the run summary counts
  each kind. Measured on two 2 kb products at 200 reads with 1, 5 and 10 % chimeric reads (half of
  them fold-backs, three replicates each): on R10-like reads every product came out exact and at
  its own length, on R9-like reads at its own length with 0.05 % error or less, and every fold-back
  and ligation chimera was listed. The engine before these checks had
  extended the contigs by the fold-back tails and built contigs out of the ligation chimeras. One
  caveat: a product whose two ends are inverted repeats of each other over 50 nt or more
  (mirrored primer sites, terminal inverted repeats) looks like a fold-back in every full-length
  read, and each read would be cut to its middle; turn the trimming off for such a product with
  `--set foldbackMin=0` (the report then says "fold-back copies off").
- **Pooled libraries.** When one barcode holds many different products, each at low depth, most
  reads stay singletons whatever the settings, because no second read of the same product is in the
  sample; the number of contigs and the singleton count in `bins.tsv` say so at a glance. Raising
  `--max-reads`, or grouping by an amplicon tag instead of the barcode, is the remedy.

## Performance and reproducibility

Almost all the time goes into aligning candidate read pairs, which is spread across
`worker_threads` with the same work division the web page uses for Web Workers; seeding, layout and
consensus stay on the main thread. Measured on a four-core machine (Node 26):

| Data | Threads | Time | Peak memory |
|---|---|---|---|
| 1000 R10-like reads of one 2 kb amplicon | 1 | 375 s | 356 MB |
| the same | 4 | 128 s | 491 MB |
| 300 R10-like reads, three 2 kb amplicons (the bundled example) | 3 | 25 s | 217 MB |
| 2 GB FASTQ, 529,235 reads, three barcodes, 1000 reads kept per barcode | 3 | 483 s | 904 MB |
| Nanopore R9.4.1 run file, 4000 reads, 1000 of ≥ 200 nt kept | 3 | 9 s | 259 MB |

The output is byte-identical whatever the thread count and whatever form the input takes: the same
reads as plain, wrapped, indented or gzipped FASTQ, or as FASTA, give the same contigs and the same
report (only the elapsed time and the thread count in the summary differ). The default heap of Node
suffices for groups of a few thousand reads; for larger groups raise it with
`node --max-old-space-size=16000 tools/cap-cli.js …`.

## Examples

The `examples/` folder holds two synthetic data sets with the templates they were made from and
the generator that made them, one file of a real Oxford Nanopore R9.4.1 run, the expected contigs of
each, and a runner that assembles them all and compares the result with the expected output:

```bash
node examples/run-examples.js
```

See [`examples/README.md`](examples/README.md) for the commands, what each example shows and the
figures to expect.

## The browser version

The same engine runs in the browser as the [OLC Assembler](https://digitalgens.org/cap.html) page
of DigitalGens, with an interactive alignment view, downloads of every artefact and the overlap
stage spread across Web Workers; its help page describes the parameters and the report in detail.
The page takes what fits in a text box, up to 5000 reads; this command-line mode exists for whole
runs. The two share `js/cap.js` and `js/cap-report.js` unchanged, so results agree exactly.

## Repository layout

```
tools/cap-cli.js        the command-line batch mode (also the worker it spawns)
js/cap.js               the assembly engine, DOM-free, shared with the web page
js/cap-report.js        report, FASTA, placements and alignment text from an assembly result
examples/               bundled data, expected results, generator, runner (see examples/README.md)
package.json            metadata; "bin" makes `olc-assembler` available after npm link
CITATION.cff            citation metadata
LICENSE.txt             GNU General Public License v3
```

## Citation

If this software contributed to your work, please cite the DigitalGens platform paper:

- Kalendar R. (2025) Comprehensive web-based platform for advanced PCR design, genotyping,
  synthetic biology, molecular diagnostics, and sequence analysis. *Molecular Therapy Nucleic
  Acids* **36**(4):102716. [doi:10.1016/j.omtn.2025.102716](https://doi.org/10.1016/j.omtn.2025.102716)

The assembly strategy follows CAP3:

- Huang X., Madan A. (1999) CAP3: a DNA sequence assembly program. *Genome Research* **9**:868–877.
  [doi:10.1101/gr.9.9.868](https://doi.org/10.1101/gr.9.9.868)

`CITATION.cff` carries the same information in machine-readable form.

## Licence

Free and open-source software under the **GNU General Public License version 3** (see
[`LICENSE.txt`](LICENSE.txt)), which permits use, study, modification and redistribution by
anyone, including commercially. The code is an independent re-implementation; CAP3 is acknowledged
in the source header. If you modify or extend this code, please include:

> Based on the work by Ruslan Kalendar — https://github.com/rkalendar/
