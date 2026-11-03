// make-reads.js - synthetic Oxford Nanopore-like amplicon reads for trying the OLC Assembler.
//
// usage: node make-reads.js <out.fastq> [--amplicons A] [--depth D] [--length L] [--model R10|R9]
//                           [--bytes N] [--wrap W] [--seed S] [--templates out.fasta]
//
// A random templates of L nt are drawn from the seed; D reads are written per template, each a
// copy with its end positions jittered by up to 20 nt, the errors of the chosen model applied
// (R10: about 2 %, R9: about 7 %, both concentrated in homopolymer runs) and the strand chosen
// at random. Each template gets its own "barcodeNN" tag in the read header, so the reads can be
// grouped with --bin-regex "barcode(\d+)". --bytes keeps writing records, cycling over the
// templates, until the file reaches N bytes (for load tests); --wrap W writes sequence and
// quality over lines of W characters. The same seed gives the same templates for both models.
'use strict';
const fs = require('fs');
const args = process.argv.slice(2);
const out = args[0];
const opt = (k, d) => { const i = args.indexOf(k); return i < 0 ? d : args[i + 1]; };
const A = +opt('--amplicons', 1), D = +opt('--depth', 1000), L = +opt('--length', 2000);
const MODEL = opt('--model', 'R10'), BYTES = +opt('--bytes', 0), WRAP = +opt('--wrap', 0), SEED = +opt('--seed', 20261102);
const TEMPL = opt('--templates', '');
if (!out) { console.error('usage: node make-reads.js <out.fastq> [options]'); process.exit(1); }

function rng(seed) { let s = seed >>> 0; return () => { s = (Math.imul(s, 1664525) + 1013904223) >>> 0; return s / 4294967296; }; }
function randSeq(n, r) { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; }
function rc(s) { const m = { A: 'T', C: 'G', G: 'C', T: 'A' }; let o = ''; for (let i = s.length - 1; i >= 0; i--) o += m[s[i]]; return o; }
const MODELS = {
  R10: { sub: 0.004, ins: 0.005, del: 0.008, hp: 0.03 },   // ~2 % error, homopolymer-biased
  R9:  { sub: 0.020, ins: 0.020, del: 0.030, hp: 0.08 },   // ~7 % error
};
const p = MODELS[MODEL];
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

const r = rng(SEED);
const templates = Array.from({ length: A }, () => randSeq(L, r));
if (TEMPL) fs.writeFileSync(TEMPL, templates.map((t, i) => '>amplicon' + (i + 1) + ' barcode' + String(i + 1).padStart(2, '0') + '\n' + t + '\n').join(''));

const ws = fs.createWriteStream(out);
let written = 0, n = 0;
const uuid = (k) => (k.toString(16).padStart(8, '0') + '-1234-4abc-8def-' + (k * 2654435761 % 0xffffffffffff).toString(16).padStart(12, '0')).slice(0, 36);
function record(k) {
  const a = BYTES ? k % A : Math.floor(k / D);
  const t = templates[a];
  const from = Math.floor(r() * 21), to = t.length - Math.floor(r() * 21);
  let s = mutate(t.slice(from, to), r);
  if (r() < 0.5) s = rc(s);
  const hdr = '@' + uuid(k + 1) + ' runid=5f3a9c barcode' + String(a + 1).padStart(2, '0') + ' read=' + (k + 1) + ' ch=' + (1 + (k % 512)) + ' start_time=2026-11-02T08:00:00Z';
  return hdr + '\n' + cut(s, WRAP) + '\n+\n' + cut(qual(s.length, r), WRAP) + '\n';
}
const total = BYTES ? Infinity : A * D;
(async () => {
  while (n < total && (!BYTES || written < BYTES)) {
    let chunk = '';
    for (let i = 0; i < 200 && n < total && (!BYTES || written < BYTES); i++) { const rec = record(n++); chunk += rec; written += Buffer.byteLength(rec); }
    if (!ws.write(chunk)) await new Promise(res => ws.once('drain', res));
  }
  ws.end();
  await new Promise(res => ws.on('finish', res));
  console.log(`${out}: ${n} reads, ${A} amplicon(s) of ${L} nt, model ${MODEL}, ${(written / 1048576).toFixed(1)} MB`);
})();
