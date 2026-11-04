// make-reads.js - synthetic Oxford Nanopore-like amplicon reads for trying the OLC Assembler and
// for the benchmark sets.
//
// usage: node make-reads.js <out.fastq> [--amplicons A] [--depth D[,D2,...]] [--length L] [--model R10|R9]
//                           [--templates out.fasta] [--templates-in in.fasta] [--homopolymers N]
//                           [--allele P] [--allele-ratio R] [--variants out.tsv] [--chimera F]
//                           [--shuffle] [--bytes N] [--wrap W] [--seed S]
//
// A random templates of L nt are drawn from the seed (or read from --templates-in); D reads are
// written per template (a comma list gives one depth per template, for unequal abundance), each a
// copy with its end positions jittered by up to 20 nt, the errors of the chosen model applied
// (R10: about 2 %, R9: about 7 %, both concentrated in homopolymer runs) and the strand chosen
// at random. Each template gets its own "barcodeNN" tag in the read header, so the reads can be
// grouped with --bin-regex "barcode(\d+)". The same seed gives the same templates for both models.
//
// Benchmark options. --homopolymers N inserts N runs of 4-10 nt into every template before any
// reads are made. --allele P makes a second allele of every template with P % substitutions at
// random positions (at least one) and draws the fraction R (--allele-ratio, default 0.5) of the
// template's reads from it; the header carries "allele=1" or "allele=2" and --variants writes the
// positions (template, 1-based position, base of allele 1, base of allele 2). --chimera F turns
// the fraction F of reads into chimeras: a read of one template continued by the end of another
// ("chimera=junction:a,b"), or a read followed by the reverse complement of its own start
// ("chimera=foldback"). --shuffle writes the reads in random order instead of template by
// template. --bytes keeps writing records, cycling over the templates, until the file reaches N
// bytes (for load tests); --wrap W writes sequence and quality over lines of W characters.
'use strict';
const fs = require('fs');
const args = process.argv.slice(2);
const out = args[0];
const opt = (k, d) => { const i = args.indexOf(k); return i < 0 ? d : args[i + 1]; };
const has = k => args.includes(k);
const A = +opt('--amplicons', 1), L = +opt('--length', 2000);
const DEPTHS = String(opt('--depth', '1000')).split(',').map(Number);
const MODEL = opt('--model', 'R10'), BYTES = +opt('--bytes', 0), WRAP = +opt('--wrap', 0), SEED = +opt('--seed', 20261102);
const TEMPL = opt('--templates', ''), TEMPL_IN = opt('--templates-in', '');
const HP = +opt('--homopolymers', 0);
const ALLELE = +opt('--allele', 0), RATIO = +opt('--allele-ratio', 0.5), VARIANTS = opt('--variants', '');
const CHIMERA = +opt('--chimera', 0);
const SHUFFLE = has('--shuffle');
if (!out) { console.error('usage: node make-reads.js <out.fastq> [options]'); process.exit(1); }

function rng(seed) { let s = seed >>> 0; return () => { s = (Math.imul(s, 1664525) + 1013904223) >>> 0; return s / 4294967296; }; }
function randSeq(n, r) { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; }
function rc(s) { const m = { A: 'T', C: 'G', G: 'C', T: 'A' }; let o = ''; for (let i = s.length - 1; i >= 0; i--) o += m[s[i]]; return o; }
const MODELS = {
  R10: { sub: 0.004, ins: 0.005, del: 0.008, hp: 0.03 },   // ~2 % error, homopolymer-biased
  R9:  { sub: 0.020, ins: 0.020, del: 0.030, hp: 0.08 },   // ~7 % error
};
const p = MODELS[MODEL];
if (!p) { console.error('unknown model ' + MODEL + ' (R10 or R9)'); process.exit(1); }
function mutate(t, r) {
  let o = '';
  for (let i = 0; i < t.length; i++) {
    const c = t[i];
    let run = 1, j = i - 1; while (j >= 0 && t[j] === c) { run++; j--; } j = i + 1; while (j < t.length && t[j] === c) { run++; j++; }
    const del = p.del + (run >= 4 ? p.hp : 0);
    let x = r();
    if (x < del) continue;
    x -= del;
    if (x < p.sub) { o += 'ACGT'.replace(c, '')[Math.floor(r() * 3)]; continue; }
    x -= p.sub;
    o += c;
    if (x < p.ins) o += 'ACGT'[Math.floor(r() * 4)];
  }
  return o;
}
// Nanopore-like quality: PHRED 7..40, skewed high, as Sanger characters.
function qual(n, r) { let q = ''; for (let i = 0; i < n; i++) q += String.fromCharCode(33 + 7 + Math.floor(Math.pow(r(), 0.6) * 34)); return q; }
function cut(s, w) { if (!w || s.length <= w) return s; let o = ''; for (let i = 0; i < s.length; i += w) o += (i ? '\n' : '') + s.slice(i, i + w); return o; }
function readFasta(file) {
  const recs = []; let cur = null;
  for (const line of fs.readFileSync(file, 'utf8').split(/\r?\n/)) {
    if (line[0] === '>') { cur = { name: line.slice(1).trim().split(/\s+/)[0], seq: '' }; recs.push(cur); }
    else if (cur) cur.seq += line.trim().toUpperCase();
  }
  return recs;
}

const r = rng(SEED);
const templates = TEMPL_IN ? readFasta(TEMPL_IN)
  : Array.from({ length: A }, (_, i) => ({ name: 'amplicon' + (i + 1), seq: randSeq(L, r) }));
const NT = templates.length;
if (!NT) { console.error('no templates'); process.exit(1); }

// Homopolymer stress: insert N runs of 4-10 nt into every template.
if (HP > 0) for (const t of templates) {
  for (let k = 0; k < HP; k++) {
    const pos = 1 + Math.floor(r() * (t.seq.length - 1));
    const base = 'ACGT'[Math.floor(r() * 4)];
    const len = 4 + Math.floor(r() * 7);
    t.seq = t.seq.slice(0, pos) + base.repeat(len) + t.seq.slice(pos);
  }
}

// Alleles: a second copy of every template with P % substitutions.
const alleles = [], variants = [];
if (ALLELE > 0) for (const t of templates) {
  const n = Math.max(1, Math.round(ALLELE / 100 * t.seq.length));
  const positions = new Set();
  while (positions.size < n) positions.add(Math.floor(r() * t.seq.length));
  const s = t.seq.split('');
  for (const pos of [...positions].sort((a, b) => a - b)) {
    const b1 = s[pos], b2 = 'ACGT'.replace(b1, '')[Math.floor(r() * 3)];
    s[pos] = b2;
    variants.push([t.name, pos + 1, b1, b2]);
  }
  alleles.push(s.join(''));
}

const tag = i => 'barcode' + String(i + 1).padStart(2, '0');
if (TEMPL) {
  let txt = '';
  templates.forEach((t, i) => {
    txt += '>' + t.name + ' ' + tag(i) + (ALLELE ? ' allele=1' : '') + '\n' + t.seq + '\n';
    if (ALLELE) txt += '>' + t.name + '_allele2 ' + tag(i) + ' allele=2\n' + alleles[i] + '\n';
  });
  fs.writeFileSync(TEMPL, txt);
}
if (VARIANTS) fs.writeFileSync(VARIANTS, 'template\tposition\tallele1\tallele2\n' + variants.map(v => v.join('\t') + '\n').join(''));

// Which template and allele each read comes from.
const plan = [];
if (!BYTES) {
  templates.forEach((t, a) => {
    const d = DEPTHS[a % DEPTHS.length];
    const n2 = ALLELE ? Math.round(RATIO * d) : 0;
    for (let i = 0; i < d; i++) plan.push({ a, allele: i < n2 ? 2 : 1 });
  });
  if (SHUFFLE) for (let i = plan.length - 1; i > 0; i--) { const j = Math.floor(r() * (i + 1)); [plan[i], plan[j]] = [plan[j], plan[i]]; }
}

function readOf(seq, r) {
  const from = Math.floor(r() * 21), to = seq.length - Math.floor(r() * 21);
  return mutate(seq.slice(from, to), r);
}

const ws = fs.createWriteStream(out);
let written = 0, n = 0;
const uuid = (k) => (k.toString(16).padStart(8, '0') + '-1234-4abc-8def-' + (k * 2654435761 % 0xffffffffffff).toString(16).padStart(12, '0')).slice(0, 36);
function record(k) {
  const { a, allele } = BYTES ? { a: k % NT, allele: ALLELE && r() < RATIO ? 2 : 1 } : plan[k];
  const t = allele === 2 ? alleles[a] : templates[a].seq;
  let s, extra = '';
  if (CHIMERA > 0 && r() < CHIMERA) {
    if (NT > 1 && r() < 0.5) {
      let b = Math.floor(r() * (NT - 1)); if (b >= a) b++;
      const tb = allele === 2 ? alleles[b] : templates[b].seq;
      const cutA = Math.round((0.3 + 0.4 * r()) * t.length), cutB = Math.round((0.3 + 0.4 * r()) * tb.length);
      s = mutate(t.slice(Math.floor(r() * 21), cutA), r) + mutate(tb.slice(tb.length - cutB, tb.length - Math.floor(r() * 21)), r);
      extra = ' chimera=junction:' + (a + 1) + ',' + (b + 1);
    } else {
      const body = readOf(t, r);
      s = body + rc(body.slice(0, Math.round((0.3 + 0.4 * r()) * body.length)));
      extra = ' chimera=foldback';
    }
  } else {
    s = readOf(t, r);
  }
  if (r() < 0.5) s = rc(s);
  const hdr = '@' + uuid(k + 1) + ' runid=5f3a9c ' + tag(a) + (ALLELE ? ' allele=' + allele : '') + extra + ' read=' + (k + 1) + ' ch=' + (1 + (k % 512)) + ' start_time=2026-11-02T08:00:00Z';
  return hdr + '\n' + cut(s, WRAP) + '\n+\n' + cut(qual(s.length, r), WRAP) + '\n';
}
const total = BYTES ? Infinity : plan.length;
(async () => {
  while (n < total && (!BYTES || written < BYTES)) {
    let chunk = '';
    for (let i = 0; i < 200 && n < total && (!BYTES || written < BYTES); i++) { const rec = record(n++); chunk += rec; written += Buffer.byteLength(rec); }
    if (!ws.write(chunk)) await new Promise(res => ws.once('drain', res));
  }
  ws.end();
  await new Promise(res => ws.on('finish', res));
  const lens = templates.map(t => t.seq.length);
  const lenTxt = Math.min(...lens) === Math.max(...lens) ? lens[0] + ' nt' : Math.min(...lens) + '-' + Math.max(...lens) + ' nt';
  console.log(`${out}: ${n} reads, ${NT} amplicon(s) of ${lenTxt}, model ${MODEL}` +
    (ALLELE ? `, alleles at ${ALLELE} % (${variants.length} variant${variants.length === 1 ? '' : 's'}, ratio ${RATIO})` : '') +
    (CHIMERA ? `, chimera fraction ${CHIMERA}` : '') + (HP ? `, ${HP} homopolymer runs per template` : '') +
    `, ${(written / 1048576).toFixed(1)} MB`);
})();
