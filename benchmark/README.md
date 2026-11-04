# Benchmark scripts

Scripts for comparing the OLC Assembler with other reference-free consensus tools on synthetic
Oxford Nanopore-like amplicon reads. They need Node.js 18 or later and nothing else; the other
tools are run through their own command lines.

| Script | Purpose |
|---|---|
| `generate-sets.js` | makes the sets S1–S5 (single product, homopolymer stress, product mixtures, alleles, chimeras) with `../examples/synthetic/make-reads.js`, each case as `reads/<case>.fastq` and `truth/<case>.templates.fasta` (plus `.variants.tsv` for the allele set), listed in `manifest.tsv` |
| `bench.js` | runs one tool over one set, case by case, recording wall time, CPU time, peak memory and exit code (`time.json`), copies the tool's consensus to `consensus.fasta` and scores it with `evaluate.js` (`score.tsv`, `details.tsv`); `--collect` gathers every `score.tsv` into one table |
| `evaluate.js` | scores any consensus FASTA or FASTQ against the templates: closest template on either strand, global alignment, substitutions, insertions and deletions, homopolymer context, end trimming, ambiguity codes (IUPAC-aware and strict counts), recovered, missed, duplicate and chimeric consensuses, and the handling of a two-allele template (split, coded, major, minor, mixed, lost) |
| `summarise.js` | aggregates the collected table into one row per condition (replicates pooled): recovery, unmatched, extended, duplicate and chimeric consensuses, mean internal error, error types, homopolymer share, end trimming, time, memory and allele classes, as TSV and Markdown |

```bash
node benchmark/generate-sets.js --runs D:/OLC/runs            # pilot design; --full for the complete one
node benchmark/bench.js --set S1 --tool olc --runs D:/OLC/runs --threads 4
node benchmark/bench.js --set S1 --tool amplicon_sorter --runs D:/OLC/runs --filter "_r1$"
node benchmark/bench.js --collect --runs D:/OLC/runs --results D:/OLC/results/summary.tsv
node benchmark/summarise.js --summary D:/OLC/results/summary.tsv     # -> summary-by-condition.tsv and .md
node benchmark/evaluate.js --consensus consensus.fasta --truth truth/case.templates.fasta [--variants truth/case.variants.tsv]
```

`bench.js` expects the other tools under `<runs>/../software/<tool>/` (or `--software DIR`):
amplicon_sorter as a clone of its repository (the single-core script is used on Windows, the
multiprocessing one elsewhere). Linux-only tools (NGSpeciesID, Decona, medaka smolecule,
wf-amplicon) are timed with `/usr/bin/time -v` and are added to `TOOLS` in `bench.js` as they are
installed.

## Scoring rules

- A consensus is assigned to the template (either strand) with the smallest IUPAC-aware edit
  distance; candidate templates are found by shared 12-mers first.
- The alignment has three parts: leading gaps in either sequence, a core aligned globally, and
  trailing gaps. A gap base at either end costs 1 and an edit inside the core costs 4, so a
  consensus shorter than the template or overhanging it is reported as `term_del` and `term_ins`
  (end trimming or overhang, not internal errors) and a chance match inside a trimmed or junk end
  cannot pull it into the core; the price is that an error within the two terminal bases is
  reported as end trimming. `err_internal_*` is (substitutions + insertions + deletions) over the
  covered template length; `err_global_*` is the plain IUPAC-aware edit distance over the template
  length, end gaps included. The `strict` variants count every ambiguity code as an error; the
  default counts a code as correct when it contains the template base and reports such codes as
  `codes_matching`.
- An error is `hp` when the template has a run of three or more identical bases at that position
  (for an insertion, next to it).
- Extended: more than 10 % of the template length (at least 30 nt) beyond the template's ends,
  such as a contig carrying the tail of a foldback chimera.
- Matched: at most 1 % internal error over at least 90 % of the template, and not extended.
  Recovered: a template with at least one matched consensus. Duplicates: matched consensuses beyond the first per
  template (per allele for two-allele templates). Chimeric: the two halves of a consensus fit
  different templates (infix alignment within 5 %).
- Alleles (`--variants`): `split` when the two alleles get separate matched consensuses; otherwise
  the variant positions of the single consensus are read as `coded` (ambiguity code containing
  both bases), `major` (allele 1), `minor` (allele 2) or `other`, and the class is the majority
  (`mixed` when none reaches half); `lost` when nothing matched; `+dup` when extra consensuses of
  the same allele exist.
