// evaluate.js - score the consensus sequences of any tool against the templates the reads were
// made from (examples/synthetic/make-reads.js --templates ... --variants ...).
//
// usage: node evaluate.js --consensus FILE --truth templates.fasta [--variants variants.tsv]
//                         [--pattern REGEX] [--time time.json] [--tag key=value ...]
//                         [--append summary.tsv] [--details details.tsv] [--quiet]
//
// Every record of the consensus file (or those whose name matches --pattern) is assigned to the
// closest template on either strand and aligned to it globally. The distance is counted twice:
// IUPAC-aware (an ambiguity code in the consensus matches when the template base is one of its
// alternatives) and strict (every ambiguity code is an error). Substitutions, insertions and
// deletions are counted separately, and each is "hp" when the template has a homopolymer run of
// three or more at that position (an insertion: next to one). Gap runs at the two ends of the
// alignment are end trimming or overhang, not internal errors. A consensus is "matched" when its
// internal error is at most 1 % over at least 90 % of the template and it is not "extended" (more
// than 10 % of the template length, at least 30 nt, beyond the template's ends); a template is
// "recovered" when at least one consensus matches it; matched consensuses beyond the first per template (or per
// allele) are "duplicates"; a consensus whose two halves best fit different templates is
// "chimeric". With --variants the consensus of a two-allele template is classified: "split" (one
// consensus per allele), "coded" (one consensus, ambiguity codes at the variant positions),
// "major" or "minor" (one consensus carrying one allele), "mixed", or "lost".
//
// --append adds one summary row (the --tag columns first, then the metrics, then the wall time,
// peak memory and exit code read from --time) to a TSV file, writing the header when the file is
// new; --details adds one row per consensus.
'use strict';
const fs = require('fs');
const path = require('path');

const argv = process.argv.slice(2);
function opt(k, d) { const i = argv.indexOf(k); return i < 0 ? d : argv[i + 1]; }
const consensusFile = opt('--consensus', ''), truthFile = opt('--truth', ''), variantsFile = opt('--variants', '');
const patternTxt = opt('--pattern', ''), timeFile = opt('--time', ''), appendFile = opt('--append', ''), detailsFile = opt('--details', '');
const quiet = argv.includes('--quiet');
const tags = [];
for (let i = 0; i < argv.length; i++) if (argv[i] === '--tag') { const kv = argv[++i]; const e = kv.indexOf('='); tags.push([kv.slice(0, e), kv.slice(e + 1)]); }
if (!consensusFile || !truthFile) { console.error('usage: node evaluate.js --consensus FILE --truth templates.fasta [options]'); process.exit(1); }

const BITS = { a: 1, c: 2, g: 4, t: 8, r: 5, y: 10, s: 6, w: 9, k: 12, m: 3, b: 14, d: 13, h: 11, v: 7, n: 15, u: 8 };
const COMP = { a: 't', c: 'g', g: 'c', t: 'a', u: 'a', r: 'y', y: 'r', s: 's', w: 'w', k: 'm', m: 'k', b: 'v', v: 'b', d: 'h', h: 'd', n: 'n' };
const bits = ch => BITS[ch] || 0;
const ambiguous = ch => !'acgtu'.includes(ch);

function readFasta(file) {
  const recs = []; let cur = null;
  for (const line of fs.readFileSync(file, 'utf8').split(/\r?\n/)) {
    if (line[0] === '>') { cur = { name: line.slice(1).trim().split(/\s+/)[0], header: line.slice(1).trim(), seq: '' }; recs.push(cur); }
    else if (cur && line[0] !== '@' && line[0] !== '+') cur.seq += line.trim().toLowerCase().replace(/[^a-z]/g, '');
  }
  return recs;
}
function readFastq(file) {
  const lines = fs.readFileSync(file, 'utf8').split(/\r?\n/);
  const recs = [];
  for (let i = 0; i + 1 < lines.length; i += 4) if (lines[i][0] === '@') recs.push({ name: lines[i].slice(1).trim().split(/\s+/)[0], header: lines[i].slice(1).trim(), seq: lines[i + 1].trim().toLowerCase() });
  return recs;
}
function readSeqs(file) {
  const txt = fs.readFileSync(file, 'utf8');
  const first = txt.trimStart()[0];
  return first === '@' ? readFastq(file) : readFasta(file);
}
function revComp(s) { let o = ''; for (let i = s.length - 1; i >= 0; i--) o += COMP[s[i]] || 'n'; return o; }

// Global edit distance, IUPAC-aware, no traceback.
function editDistance(a, b) {
  const n = a.length, m = b.length;
  let prev = new Int32Array(m + 1), cur = new Int32Array(m + 1);
  const bb = new Uint8Array(m); for (let j = 0; j < m; j++) bb[j] = bits(b[j]);
  for (let j = 0; j <= m; j++) prev[j] = j;
  for (let i = 1; i <= n; i++) {
    cur[0] = i;
    const ai = bits(a[i - 1]);
    for (let j = 1; j <= m; j++) {
      const cost = (ai & bb[j - 1]) ? 0 : 1;
      let v = prev[j - 1] + cost;
      if (prev[j] + 1 < v) v = prev[j] + 1;
      if (cur[j - 1] + 1 < v) v = cur[j - 1] + 1;
      cur[j] = v;
    }
    [prev, cur] = [cur, prev];
  }
  return prev[m];
}
// Edit distance of the whole of q against any substring of t (free leading and trailing text).
function infixDistance(q, t) {
  const n = q.length, m = t.length;
  let prev = new Int32Array(m + 1), cur = new Int32Array(m + 1);
  const tb = new Uint8Array(m); for (let j = 0; j < m; j++) tb[j] = bits(t[j]);
  for (let i = 1; i <= n; i++) {
    cur[0] = i;
    const qi = bits(q[i - 1]);
    for (let j = 1; j <= m; j++) {
      const cost = (qi & tb[j - 1]) ? 0 : 1;
      let v = prev[j - 1] + cost;
      if (prev[j] + 1 < v) v = prev[j] + 1;
      if (cur[j - 1] + 1 < v) v = cur[j - 1] + 1;
      cur[j] = v;
    }
    [prev, cur] = [cur, prev];
  }
  let best = prev[0]; for (let j = 1; j <= m; j++) if (prev[j] < best) best = prev[j];
  return best;
}
// Alignment of a consensus (a) to its template (b) in three parts: a run of leading gaps in
// either sequence, a core aligned globally, and a run of trailing gaps. A gap base at either end
// costs 1 and an edit inside the core costs EDIT, so that trimmed or overhanging ends come out as
// end gaps and a chance match inside such an end cannot pull it into the core; an error within
// the two terminal bases is reported as end trimming instead. Returns the operations as [op, i, j]
// with op M (match), A (ambiguity code containing the template base), X (substitution), I
// (insertion: consensus base not in the template) and D (deletion: template base not in the
// consensus), i and j being 1-based end positions in a and b, or null when the matrix would be
// too large.
const EDIT = 4;
function align(a, b) {
  const n = a.length, m = b.length, W = m + 1;
  if (!n || !m) return { dist: n + m, ops: [...Array(n)].map((_, i) => ['I', i + 1, 0]).concat([...Array(m)].map((_, j) => ['D', n, j + 1])) };
  if ((n + 1) * W > 1e8) return null;
  const INF = 0x3fffffff;
  // Layers: L (leading gaps only, cost i + j, needs no table), C (core), T (trailing gaps only).
  // tbC: 0 diagonal from C, 1 diagonal from L, 2 I from C, 3 D from C.
  // tbT: 0 I from C, 1 I from T, 2 D from C, 3 D from T.
  const tbC = new Uint8Array((n + 1) * W), tbT = new Uint8Array((n + 1) * W);
  let pC = new Int32Array(W).fill(INF), cC = new Int32Array(W).fill(INF);
  let pT = new Int32Array(W).fill(INF), cT = new Int32Array(W).fill(INF);
  const bb = new Uint8Array(m); for (let j = 0; j < m; j++) bb[j] = bits(b[j]);
  for (let i = 1; i <= n; i++) {
    const ai = bits(a[i - 1]), row = i * W;
    cC[0] = INF; cT[0] = INF;
    for (let j = 1; j <= m; j++) {
      const cost = (ai & bb[j - 1]) ? 0 : EDIT;
      let v = pC[j - 1] + cost, t = 0;
      const fromL = i - 1 + j - 1 + cost; if (fromL < v) { v = fromL; t = 1; }
      if (pC[j] + EDIT < v) { v = pC[j] + EDIT; t = 2; }
      if (cC[j - 1] + EDIT < v) { v = cC[j - 1] + EDIT; t = 3; }
      cC[j] = v; tbC[row + j] = t;
      let w = pC[j] + 1, u = 0;
      if (pT[j] + 1 < w) { w = pT[j] + 1; u = 1; }
      if (cC[j - 1] + 1 < w) { w = cC[j - 1] + 1; u = 2; }
      if (cT[j - 1] + 1 < w) { w = cT[j - 1] + 1; u = 3; }
      cT[j] = w; tbT[row + j] = u;
    }
    [pC, cC] = [cC, pC]; [pT, cT] = [cT, pT];
  }
  let best = pC[m], layer = 'C';
  if (pT[m] < best) { best = pT[m]; layer = 'T'; }
  if (n + m < best) { best = n + m; layer = 'L'; }
  const ops = [];
  let i = n, j = m;
  while (i > 0 || j > 0) {
    if (layer === 'L') { if (i > 0) { ops.push(['I', i, j]); i--; } else { ops.push(['D', i, j]); j--; } continue; }
    if (layer === 'T') {
      const u = tbT[i * W + j];
      if (u < 2) { ops.push(['I', i, j]); i--; layer = u === 0 ? 'C' : 'T'; } else { ops.push(['D', i, j]); j--; layer = u === 2 ? 'C' : 'T'; }
      continue;
    }
    const t = tbC[i * W + j];
    if (t < 2) {
      const ca = a[i - 1], cb = b[j - 1];
      ops.push([(bits(ca) & bits(cb)) ? (ambiguous(ca) ? 'A' : 'M') : 'X', i, j]); i--; j--;
      if (t === 1) layer = 'L';
    } else if (t === 2) { ops.push(['I', i, j]); i--; }
    else { ops.push(['D', i, j]); j--; }
  }
  ops.reverse();
  return { dist: best, ops };
}
function runLengths(s) {
  const run = new Int32Array(s.length);
  let i = 0;
  while (i < s.length) { let j = i; while (j < s.length && s[j] === s[i]) j++; for (let k = i; k < j; k++) run[k] = j - i; i = j; }
  return run;
}
const K = 12;
function kmerSet(s) { const set = new Set(); for (let i = 0; i + K <= s.length; i++) set.add(s.substr(i, K)); return set; }
function kmerHits(s, set) { let h = 0; for (let i = 0; i + K <= s.length; i++) if (set.has(s.substr(i, K))) h++; return h; }

// Truth: templates and their allele groups.
const truth = readFasta(truthFile).filter(t => t.seq.length);
for (const t of truth) { t.group = t.name.replace(/_allele\d+$/, ''); t.allele = (/_allele(\d+)$/.exec(t.name) || [, '1'])[1]; t.kmers = kmerSet(t.seq); t.run = runLengths(t.seq); }
const groups = [...new Set(truth.map(t => t.group))];
const variants = {};
if (variantsFile && fs.existsSync(variantsFile)) {
  for (const line of fs.readFileSync(variantsFile, 'utf8').split(/\r?\n/).slice(1)) {
    const f = line.split('\t'); if (f.length < 4) continue;
    (variants[f[0]] = variants[f[0]] || []).push({ pos: +f[1], b1: f[2].toLowerCase(), b2: f[3].toLowerCase() });
  }
}

// Consensus records.
let cons = fs.existsSync(consensusFile) ? readSeqs(consensusFile).filter(c => c.seq.length) : [];
if (patternTxt) { const re = new RegExp(patternTxt, 'i'); cons = cons.filter(c => re.test(c.header)); }

// Assignment of every consensus to its closest template.
function bestTruth(seq) {
  const rcSeq = revComp(seq);
  const scored = truth.map(t => { const f = kmerHits(seq, t.kmers), r = kmerHits(rcSeq, t.kmers); return { t, hits: Math.max(f, r), strand: r > f ? '-' : '+' }; });
  scored.sort((x, y) => y.hits - x.hits);
  const candidates = scored[0].hits > 0 ? scored.slice(0, 3) : scored;
  let best = null;
  for (const c of candidates) {
    for (const strand of (scored[0].hits > 0 ? [c.strand] : ['+', '-'])) {
      const s = strand === '-' ? rcSeq : seq;
      const d = editDistance(s, c.t.seq);
      if (!best || d < best.dist) best = { t: c.t, strand, dist: d, seq: s };
    }
  }
  return best;
}
function halfGroup(half) {
  const rcHalf = revComp(half);
  const scored = truth.map(t => { const f = kmerHits(half, t.kmers), r = kmerHits(rcHalf, t.kmers); return { t, hits: Math.max(f, r), strand: r > f ? '-' : '+' }; });
  scored.sort((x, y) => y.hits - x.hits);
  if (scored[0].hits === 0) return null;
  const c = scored[0];
  const d = infixDistance(c.strand === '-' ? rcHalf : half, c.t.seq);
  return d <= 0.05 * half.length ? c.t.group : null;
}

const rows = [];
for (const c of cons) {
  const b = bestTruth(c.seq);
  const row = { name: c.name, len: c.seq.length, truth: b.t.name, group: b.t.group, allele: b.t.allele, strand: b.strand, distIupac: b.dist, truthLen: b.t.seq.length };
  row.matched = b.dist <= 0.01 * b.t.seq.length;
  row.ambiguous = [...c.seq].filter(ambiguous).length;
  const al = align(b.seq, b.t.seq);
  row.subs = row.ins = row.del = row.hp = row.nonhp = row.codes = row.termDel = row.termIns = 0;
  if (al) {
    // Gaps at the two ends of the alignment are end trimming or overhang, not internal errors.
    const ops = al.ops;
    let lead = 0; while (lead < ops.length && (ops[lead][0] === 'I' || ops[lead][0] === 'D')) lead++;
    let trail = 0; while (trail < ops.length - lead && (ops[ops.length - 1 - trail][0] === 'I' || ops[ops.length - 1 - trail][0] === 'D')) trail++;
    for (let k = 0; k < ops.length; k++) {
      const [op, , j] = ops[k];
      if (op === 'M') continue;
      if (op === 'A') { row.codes++; continue; }
      if (k < lead || k >= ops.length - trail) { if (op === 'D') row.termDel++; else row.termIns++; continue; }
      let hp;
      if (op === 'X') { row.subs++; hp = b.t.run[j - 1] >= 3; }
      else if (op === 'D') { row.del++; hp = b.t.run[j - 1] >= 3; }
      else { row.ins++; hp = (j > 0 && b.t.run[j - 1] >= 3) || (j < b.t.seq.length && b.t.run[j] >= 3); }
      if (hp) row.hp++; else row.nonhp++;
    }
    row.alignment = al; row.alignedSeq = b.seq;
  }
  row.distStrict = row.distIupac + row.codes;
  row.pctIupac = 100 * row.distIupac / b.t.seq.length;
  row.pctStrict = 100 * row.distStrict / b.t.seq.length;
  row.internal = row.subs + row.ins + row.del;
  row.covered = b.t.seq.length - row.termDel;
  row.pctInternal = row.covered > 0 ? 100 * row.internal / row.covered : 100;
  row.pctInternalStrict = row.covered > 0 ? 100 * (row.internal + row.codes) / row.covered : 100;
  // Extended: more than 10 % of the template length (at least 30 nt) beyond the template's ends.
  row.extended = al ? row.termIns > Math.max(30, 0.1 * b.t.seq.length) : false;
  // Matched: at most 1 % internal error over at least 90 % of the template, and not extended.
  row.matched = al ? (row.pctInternal <= 1 && row.covered >= 0.9 * b.t.seq.length && !row.extended) : row.matched;
  row.chimeric = false;
  if (groups.length > 1 && c.seq.length >= 200) {
    const h = c.seq.length >> 1;
    const g1 = halfGroup(c.seq.slice(0, h)), g2 = halfGroup(c.seq.slice(h));
    row.chimeric = !!(g1 && g2 && g1 !== g2);
  }
  rows.push(row);
}

// Recovery per template group.
const perGroup = {};
for (const g of groups) perGroup[g] = { matched: [], alleles: new Set() };
for (const r of rows) if (r.matched) { perGroup[r.group].matched.push(r); perGroup[r.group].alleles.add(r.allele); }
const recovered = groups.filter(g => perGroup[g].matched.length > 0);
const duplicates = groups.reduce((s, g) => s + Math.max(0, perGroup[g].matched.length - perGroup[g].alleles.size), 0);
const matchedRows = rows.filter(r => r.matched);
const unmatched = rows.length - matchedRows.length;
const chimeric = rows.filter(r => r.chimeric).length;
const extended = rows.filter(r => r.extended).length;

// Allele classification (first template with variants; the benchmark sets have one).
let alleleClass = '', varN = 0, varCoded = 0, varMajor = 0, varMinor = 0, varOther = 0;
const varGroups = Object.keys(variants);
if (varGroups.length) {
  const g = varGroups[0], vs = variants[g];
  varN = vs.length;
  const m = perGroup[g] ? perGroup[g].matched : [];
  const base = truth.find(t => t.group === g && t.allele === '1');
  if (!m.length) alleleClass = 'lost';
  else if (perGroup[g].alleles.size >= 2) alleleClass = 'split';
  else {
    // One allele's worth of consensus: look at the variant positions of the best-matching one.
    const r = m.slice().sort((x, y) => x.distIupac - y.distIupac)[0];
    const al = (r.truth === base.name && r.alignment) ? r.alignment : align(r.alignedSeq, base.seq);
    const at = {};
    if (al) for (const [op, i, j] of al.ops) { if (op === 'M' || op === 'A' || op === 'X') at[j] = r.alignedSeq[i - 1]; else if (op === 'D') at[j] = '-'; }
    for (const v of vs) {
      const ch = at[v.pos] || '-';
      const bb = bits(ch), b1 = bits(v.b1), b2 = bits(v.b2);
      if ((bb & b1) && (bb & b2)) varCoded++;
      else if (bb & b1) varMajor++;
      else if (bb & b2) varMinor++;
      else varOther++;
    }
    const top = Math.max(varCoded, varMajor, varMinor, varOther);
    alleleClass = top < varN / 2 ? 'mixed' : top === varCoded ? 'coded' : top === varMajor ? 'major' : top === varMinor ? 'minor' : 'mixed';
    if (m.length > 1) alleleClass += '+dup';
  }
}

const mean = a => a.length ? a.reduce((s, x) => s + x, 0) / a.length : NaN;
const median = a => { if (!a.length) return NaN; const s = a.slice().sort((x, y) => x - y); const h = s.length >> 1; return s.length % 2 ? s[h] : (s[h - 1] + s[h]) / 2; };
const sum = (a, k) => a.reduce((s, r) => s + r[k], 0);
const fmt = x => Number.isNaN(x) ? '' : (Math.round(x * 10000) / 10000).toString();
let time = {};
if (timeFile && fs.existsSync(timeFile)) { try { time = JSON.parse(fs.readFileSync(timeFile, 'utf8')); } catch (e) { /* ignore */ } }

const summary = {
  n_consensus: rows.length, n_truth: groups.length, recovered: recovered.length, missed: groups.length - recovered.length,
  matched: matchedRows.length, unmatched, duplicates, chimeric, extended,
  err_internal_mean_pct: fmt(mean(matchedRows.map(r => r.pctInternal))), err_internal_median_pct: fmt(median(matchedRows.map(r => r.pctInternal))),
  err_internal_max_pct: fmt(matchedRows.length ? Math.max(...matchedRows.map(r => r.pctInternal)) : NaN),
  err_internal_strict_mean_pct: fmt(mean(matchedRows.map(r => r.pctInternalStrict))),
  err_global_mean_pct: fmt(mean(matchedRows.map(r => r.pctIupac))), err_global_strict_mean_pct: fmt(mean(matchedRows.map(r => r.pctStrict))),
  subs: sum(matchedRows, 'subs'), ins: sum(matchedRows, 'ins'), del: sum(matchedRows, 'del'), hp: sum(matchedRows, 'hp'), nonhp: sum(matchedRows, 'nonhp'),
  term_del: sum(matchedRows, 'termDel'), term_ins: sum(matchedRows, 'termIns'),
  codes_matching: sum(matchedRows, 'codes'), ambiguous: sum(matchedRows, 'ambiguous'),
  allele_class: alleleClass, var_n: varN, var_coded: varCoded, var_major: varMajor, var_minor: varMinor, var_other: varOther,
  wall_s: time.wall_s === undefined ? '' : time.wall_s, cpu_s: time.cpu_s === undefined ? '' : time.cpu_s, peak_mb: time.peak_mb === undefined ? '' : time.peak_mb, exit: time.exit === undefined ? '' : time.exit,
};

function appendRow(file, cols, values) {
  const line = values.join('\t') + '\n';
  if (!fs.existsSync(file) || fs.statSync(file).size === 0) { fs.mkdirSync(path.dirname(file), { recursive: true }); fs.writeFileSync(file, cols.join('\t') + '\n' + line); }
  else fs.appendFileSync(file, line);
}
if (appendFile) appendRow(appendFile, [...tags.map(t => t[0]), ...Object.keys(summary)], [...tags.map(t => t[1]), ...Object.values(summary)]);
if (detailsFile) {
  const cols = ['name', 'len', 'truth', 'strand', 'matched', 'chimeric', 'extended', 'covered', 'internal', 'pct_internal', 'pct_internal_strict', 'dist_global', 'dist_global_strict', 'subs', 'ins', 'del', 'hp', 'nonhp', 'term_del', 'term_ins', 'codes_matching', 'ambiguous'];
  for (const r of rows) appendRow(detailsFile, [...tags.map(t => t[0]), ...cols], [...tags.map(t => t[1]), r.name, r.len, r.truth, r.strand, r.matched ? 1 : 0, r.chimeric ? 1 : 0, r.extended ? 1 : 0, r.covered, r.internal, fmt(r.pctInternal), fmt(r.pctInternalStrict), r.distIupac, r.distStrict, r.subs, r.ins, r.del, r.hp, r.nonhp, r.termDel, r.termIns, r.codes, r.ambiguous]);
}
if (!quiet) {
  const label = path.basename(consensusFile);
  console.log(`${label}: ${rows.length} consensus, ${groups.length} template(s): ${recovered.length} recovered, ${groups.length - recovered.length} missed, ${unmatched} unmatched, ${duplicates} duplicate(s), ${chimeric} chimeric, ${extended} extended` +
    (matchedRows.length ? `; internal error ${summary.err_internal_mean_pct} % (strict ${summary.err_internal_strict_mean_pct} %): subs ${summary.subs}, ins ${summary.ins}, del ${summary.del} (${summary.hp} in homopolymers); ends: ${summary.term_del} nt missing, ${summary.term_ins} nt extra` : '') +
    (alleleClass ? `; alleles: ${alleleClass} (${varCoded} coded, ${varMajor} major, ${varMinor} minor, ${varOther} other of ${varN})` : ''));
  for (const r of rows) console.log(`  ${r.name} (${r.len} nt) -> ${r.truth} ${r.strand}: ${r.internal} internal error(s) over ${r.covered} nt${r.termDel ? ', ' + r.termDel + ' nt missing at the ends' : ''}${r.termIns ? ', ' + r.termIns + ' nt extra at the ends' : ''}${r.codes ? ', ' + r.codes + ' matching code(s)' : ''}${r.matched ? '' : ' UNMATCHED'}${r.chimeric ? ' CHIMERIC' : ''}${r.extended ? ' EXTENDED' : ''}`);
}
