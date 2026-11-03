// check.js - compare assembled contigs with the templates the reads were made from.
//
// usage: node check.js <contigs.fasta> [templates.fasta]
//
// For every contig the closest template is found (either strand) and the edit distance
// between the two is printed, an IUPAC code in the contig counting as a match when the
// template base is one of its alternatives. Exit code 0 when every contig is within 1 % of
// its template: R10-like reads give an exact match, R9-like reads leave a few homopolymer
// errors (0.1-0.7 % in the bundled example).
'use strict';
const fs = require('fs');
const path = require('path');
const contigsFile = process.argv[2];
const templatesFile = process.argv[3] || path.join(__dirname, 'templates.fasta');
if (!contigsFile) { console.error('usage: node check.js <contigs.fasta> [templates.fasta]'); process.exit(1); }

const BITS = { a: 1, c: 2, g: 4, t: 8, r: 5, y: 10, s: 6, w: 9, k: 12, m: 3, b: 14, d: 13, h: 11, v: 7, n: 15 };
function readFasta(file) {
  const out = [];
  let cur = null;
  for (const line of fs.readFileSync(file, 'utf8').split(/\r?\n/)) {
    if (line[0] === '>') { cur = { name: line.slice(1).split(/\s+/)[0], seq: '' }; out.push(cur); }
    else if (cur) cur.seq += line.trim().toLowerCase();
  }
  return out;
}
function revComp(s) {
  const m = { a: 't', c: 'g', g: 'c', t: 'a' };
  let o = '';
  for (let i = s.length - 1; i >= 0; i--) o += m[s[i]] || 'n';
  return o;
}
// Levenshtein distance; a contig base matches a template base when their 4-bit sets intersect.
function editDistance(a, b) {
  const n = a.length, m = b.length;
  let prev = new Int32Array(m + 1), cur = new Int32Array(m + 1);
  for (let j = 0; j <= m; j++) prev[j] = j;
  for (let i = 1; i <= n; i++) {
    cur[0] = i;
    const ai = BITS[a[i - 1]] || 0;
    for (let j = 1; j <= m; j++) {
      const cost = (ai & (BITS[b[j - 1]] || 0)) ? 0 : 1;
      cur[j] = Math.min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost);
    }
    [prev, cur] = [cur, prev];
  }
  return prev[m];
}

const templates = readFasta(templatesFile);
const contigs = readFasta(contigsFile).filter(c => /^contig/i.test(c.name));
if (!contigs.length) { console.log(`${path.basename(contigsFile)}: no contig`); process.exit(1); }
let ok = true;
for (const c of contigs) {
  let best = null;
  for (const t of templates) {
    const d = Math.min(editDistance(c.seq, t.seq), editDistance(revComp(c.seq), t.seq));
    if (!best || d < best.d) best = { t, d };
  }
  const pct = (100 * best.d / best.t.seq.length).toFixed(2);
  console.log(`${path.basename(contigsFile)}: ${c.name} (${c.seq.length} nt) is ${best.t.name} (${best.t.seq.length} nt), edit distance ${best.d} (${pct} %)`);
  if (best.d > 0.01 * best.t.seq.length) ok = false;
}
process.exit(ok ? 0 : 1);
