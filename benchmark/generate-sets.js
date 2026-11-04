// generate-sets.js - make the synthetic benchmark sets S1-S5 of the comparison plan.
//
// usage: node benchmark/generate-sets.js --runs DIR [--sets S1,S2,S3,S4,S5] [--full] [--force]
//
// Writes DIR/<set>/reads/<case>.fastq, DIR/<set>/truth/<case>.templates.fasta (and
// .variants.tsv for S4) and DIR/<set>/manifest.tsv. The pilot design (default) is a subset of
// the full design (--full) of benchmark-plan.md section 3.1:
//   S1 single product: lengths x models x depths x replicates
//   S2 homopolymer stress: 2000 nt templates with 20 inserted runs of 4-10 nt
//   S3 product mixtures: 3 and 5 products of different lengths, equal and 100:10:1 abundance, a nested pair
//   S4 alleles: one 2000 nt product and a second allele at P % divergence, mixed at ratio R
//   S5 chimeras: two 2000 nt products with a fraction F of chimeric reads
// Seeds are derived from the set, the template design and the replicate, so the templates of
// one replicate are the same across depths and models.
'use strict';
const fs = require('fs');
const path = require('path');
const { spawnSync } = require('child_process');

const argv = process.argv.slice(2);
const opt = (k, d) => { const i = argv.indexOf(k); return i < 0 ? d : argv[i + 1]; };
const RUNS = opt('--runs', ''), FULL = argv.includes('--full'), FORCE = argv.includes('--force');
const SETS = String(opt('--sets', 'S1,S2,S3,S4,S5')).split(',');
if (!RUNS) { console.error('usage: node benchmark/generate-sets.js --runs DIR [--sets S1,...] [--full]'); process.exit(1); }
const maker = path.join(__dirname, '..', 'examples', 'synthetic', 'make-reads.js');

const D = FULL ? {
  s1: { lengths: [500, 1000, 2000, 5000], depths: [5, 10, 20, 50, 100, 200, 500, 1000], reps: 10 },
  s2: { depths: [50, 200], reps: 10 },
  s3: { reps: 5 },
  s4: { divergence: [0.1, 0.5, 1, 2, 3, 5], ratios: [0.5, 0.2, 0.05], depths: [50, 200, 1000], reps: 5 },
  s5: { fractions: [0.01, 0.05, 0.10], reps: 5 },
} : {
  s1: { lengths: [1000, 2000], depths: [5, 10, 20, 50, 100, 200, 500], reps: 3 },
  s2: { depths: [50, 200], reps: 3 },
  s3: { reps: 3 },
  s4: { divergence: [0.1, 0.5, 1, 2, 5], ratios: [0.5, 0.2], depths: [50, 200], reps: 2 },
  s5: { fractions: [0.01, 0.05, 0.10], reps: 3 },
};
const MODELS = ['R10', 'R9'];

function hash(s) { let h = 2166136261; for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619) >>> 0; } return h || 1; }
function rng(seed) { let s = seed >>> 0; return () => { s = (Math.imul(s, 1664525) + 1013904223) >>> 0; return s / 4294967296; }; }
function randSeq(n, r) { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; }

const manifest = {};
function emit(set, c) {
  const dir = path.join(RUNS, set);
  fs.mkdirSync(path.join(dir, 'reads'), { recursive: true }); fs.mkdirSync(path.join(dir, 'truth'), { recursive: true });
  const reads = path.join(dir, 'reads', c.name + '.fastq');
  const truth = path.join(dir, 'truth', c.name + '.templates.fasta');
  const variants = c.allele ? path.join(dir, 'truth', c.name + '.variants.tsv') : '';
  const args = [maker, reads, '--model', c.model, '--seed', String(c.seed), '--templates', truth, '--shuffle', '--depth', String(c.depth)];
  if (c.templatesIn) args.push('--templates-in', c.templatesIn); else args.push('--amplicons', String(c.amplicons || 1), '--length', String(c.length));
  if (c.hp) args.push('--homopolymers', String(c.hp));
  if (c.allele) args.push('--allele', String(c.allele), '--allele-ratio', String(c.ratio), '--variants', variants);
  if (c.chimera) args.push('--chimera', String(c.chimera));
  (manifest[set] = manifest[set] || []).push({ case: c.name, model: c.model, length: c.length, depth: c.depth, rep: c.rep, seed: c.seed, params: c.params || '', reads, truth, variants });
  if (!FORCE && fs.existsSync(reads) && fs.existsSync(truth)) return;
  const r = spawnSync(process.execPath, args, { encoding: 'utf8' });
  if (r.status !== 0) { console.error(`${set}/${c.name} FAILED\n${r.stderr}`); process.exitCode = 1; return; }
  process.stdout.write(`${set}/${c.name}: ${r.stdout.trim().split(':').slice(1).join(':').trim()}\n`);
}

if (SETS.includes('S1')) for (const L of D.s1.lengths) for (let rep = 1; rep <= D.s1.reps; rep++) {
  const seed = hash(`S1 L${L} rep${rep}`);
  for (const model of MODELS) for (const depth of D.s1.depths)
    emit('S1', { name: `L${L}_${model}_d${depth}_r${rep}`, model, length: L, depth, rep, seed, params: 'single product' });
}
if (SETS.includes('S2')) for (let rep = 1; rep <= D.s2.reps; rep++) {
  const seed = hash(`S2 rep${rep}`);
  for (const model of MODELS) for (const depth of D.s2.depths)
    emit('S2', { name: `HP20_${model}_d${depth}_r${rep}`, model, length: 2000, depth, rep, seed, hp: 20, params: '20 homopolymer runs of 4-10 nt' });
}
if (SETS.includes('S3')) {
  const designs = [
    { key: 'mix3eq', lengths: [500, 1200, 3000], depth: '100,100,100', params: '3 products, equal' },
    { key: 'mix3ratio', lengths: [500, 1200, 3000], depth: '200,20,2', params: '3 products, 100:10:1' },
    { key: 'mix5eq', lengths: [400, 800, 1500, 2500, 5000], depth: '100,100,100,100,100', params: '5 products, equal' },
    { key: 'nested2', lengths: [3000, 'nested'], depth: '100,100', params: 'product 2 = bases 601-1800 of product 1' },
  ];
  for (const d of designs) for (let rep = 1; rep <= D.s3.reps; rep++) {
    const seed = hash(`S3 ${d.key} rep${rep}`);
    const r = rng(seed);
    const dir = path.join(RUNS, 'S3', 'truth'); fs.mkdirSync(dir, { recursive: true });
    const tin = path.join(dir, `${d.key}_r${rep}.templates-in.fasta`);
    const seqs = [];
    for (const L of d.lengths) seqs.push(L === 'nested' ? seqs[0].slice(600, 1800) : randSeq(L, r));
    fs.writeFileSync(tin, seqs.map((s, i) => `>amplicon${i + 1}\n${s}\n`).join(''));
    for (const model of MODELS)
      emit('S3', { name: `${d.key}_${model}_r${rep}`, model, length: seqs.map(s => s.length).join(','), depth: d.depth, rep, seed: seed + 1, templatesIn: tin, params: d.params });
  }
}
if (SETS.includes('S4')) for (const P of D.s4.divergence) for (let rep = 1; rep <= D.s4.reps; rep++) {
  const seed = hash(`S4 P${P} rep${rep}`);
  for (const ratio of D.s4.ratios) for (const model of MODELS) for (const depth of D.s4.depths)
    emit('S4', { name: `P${P}_R${ratio}_${model}_d${depth}_r${rep}`, model, length: 2000, depth, rep, seed, allele: P, ratio, params: `alleles ${P} %, minor fraction ${ratio}` });
}
if (SETS.includes('S5')) for (const F of D.s5.fractions) for (let rep = 1; rep <= D.s5.reps; rep++) {
  const seed = hash(`S5 rep${rep}`);
  for (const model of MODELS)
    emit('S5', { name: `C${F}_${model}_d200_r${rep}`, model, length: 2000, depth: 200, rep, seed, amplicons: 2, chimera: F, params: `2 products, chimera fraction ${F}` });
}

for (const set of Object.keys(manifest)) {
  const rows = manifest[set];
  const cols = ['case', 'model', 'length', 'depth', 'rep', 'seed', 'params', 'reads', 'truth', 'variants'];
  fs.writeFileSync(path.join(RUNS, set, 'manifest.tsv'), cols.join('\t') + '\n' + rows.map(r => cols.map(c => r[c]).join('\t') + '\n').join(''));
  console.log(`${set}: ${rows.length} cases -> ${path.join(RUNS, set, 'manifest.tsv')}`);
}
