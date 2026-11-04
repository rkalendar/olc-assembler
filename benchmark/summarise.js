#!/usr/bin/env node
// summarise.js - aggregate the per-run table written by `bench.js --collect` into one row per
// condition (set, tool, condition, case without its replicate suffix): recovery, unmatched,
// extended, duplicate and chimeric consensuses, internal error of the matched consensuses, the
// error types, the homopolymer share, end trimming, time, memory and the allele classes.
//
// usage: node summarise.js [--summary D:/OLC/results/summary.tsv] [--out D:/OLC/results/summary-by-condition]
//        writes <out>.tsv and <out>.md
'use strict';
const fs = require('fs');
const path = require('path');
const args = process.argv.slice(2);
const opt = (k, d) => { const i = args.indexOf(k); return i >= 0 && i + 1 < args.length ? args[i + 1] : d; };
const summaryFile = opt('--summary', 'D:/OLC/results/summary.tsv');
const outBase = opt('--out', path.join(path.dirname(summaryFile), 'summary-by-condition'));

const lines = fs.readFileSync(summaryFile, 'utf8').split(/\r?\n/).filter(l => l.trim());
const header = lines[0].split('\t');
const rows = lines.slice(1).map(l => { const f = l.split('\t'); const o = {}; header.forEach((h, i) => o[h] = f[i] === undefined ? '' : f[i]); return o; });
const num = v => v === '' || v === undefined ? NaN : +v;

const groups = new Map();
for (const r of rows) {
  const cond = r.case.replace(/_r\d+$/, '');
  const key = [r.set, r.tool, r.condition, cond].join('\t');
  if (!groups.has(key)) groups.set(key, { set: r.set, tool: r.tool, condition: r.condition, cond, model: r.model, length: num(r.length), depth: num(r.depth), params: r.params, rows: [] });
  groups.get(key).rows.push(r);
}
const sum = (rs, k) => rs.reduce((s, r) => s + (num(r[k]) || 0), 0);
const mean = a => a.length ? a.reduce((s, x) => s + x, 0) / a.length : NaN;
const fmt = (x, d = 2) => Number.isNaN(x) || x === undefined ? '' : (Math.round(x * 10 ** d) / 10 ** d).toString();
const natural = (a, b) => a.localeCompare(b, 'en', { numeric: true });
const modelRank = m => m === 'R10' ? 0 : m === 'R9' ? 1 : 2;

const out = [];
for (const g of groups.values()) {
  const rs = g.rows, withMatch = rs.filter(r => num(r.matched) > 0);
  const hp = sum(withMatch, 'hp'), nonhp = sum(withMatch, 'nonhp'), matched = sum(rs, 'matched');
  const classes = {};
  for (const r of rs) if (r.allele_class) classes[r.allele_class] = (classes[r.allele_class] || 0) + 1;
  out.push({
    set: g.set, tool: g.tool, condition: g.condition, case_group: g.cond, model: g.model, length: g.length, depth: g.depth, params: g.params,
    reps: rs.length, truth: sum(rs, 'n_truth'), recovered: sum(rs, 'recovered'),
    recovered_pct: fmt(100 * sum(rs, 'recovered') / Math.max(1, sum(rs, 'n_truth')), 1),
    consensus: sum(rs, 'n_consensus'), matched, unmatched: sum(rs, 'unmatched'), extended: sum(rs, 'extended'), duplicates: sum(rs, 'duplicates'), chimeric: sum(rs, 'chimeric'),
    err_internal_mean_pct: fmt(mean(withMatch.map(r => num(r.err_internal_mean_pct))), 4),
    err_internal_max_pct: fmt(withMatch.length ? Math.max(...withMatch.map(r => num(r.err_internal_max_pct))) : NaN, 4),
    err_internal_strict_mean_pct: fmt(mean(withMatch.map(r => num(r.err_internal_strict_mean_pct))), 4),
    subs: sum(withMatch, 'subs'), ins: sum(withMatch, 'ins'), del: sum(withMatch, 'del'),
    hp_share_pct: fmt(hp + nonhp ? 100 * hp / (hp + nonhp) : NaN, 0),
    codes_matching: sum(withMatch, 'codes_matching'), ambiguous: sum(withMatch, 'ambiguous'),
    term_del_per_consensus: fmt(matched ? sum(withMatch, 'term_del') / matched : NaN, 1),
    term_ins_per_consensus: fmt(matched ? sum(withMatch, 'term_ins') / matched : NaN, 1),
    wall_s_mean: fmt(mean(rs.map(r => num(r.wall_s)).filter(x => !Number.isNaN(x))), 1),
    cpu_s_mean: fmt(mean(rs.map(r => num(r.cpu_s)).filter(x => !Number.isNaN(x))), 1),
    peak_mb_max: fmt(Math.max(...rs.map(r => num(r.peak_mb)).filter(x => !Number.isNaN(x)).concat([-Infinity])) === -Infinity ? NaN : Math.max(...rs.map(r => num(r.peak_mb)).filter(x => !Number.isNaN(x))), 0),
    nonzero_exit: rs.filter(r => r.exit !== '' && r.exit !== '0').length,
    allele_classes: Object.keys(classes).sort().map(k => `${k}:${classes[k]}`).join(' '),
  });
}
out.sort((a, b) => natural(a.set, b.set) || natural(a.tool, b.tool) || natural(a.condition, b.condition) || modelRank(a.model) - modelRank(b.model) || natural(a.case_group, b.case_group));

const cols = Object.keys(out[0] || {});
fs.writeFileSync(outBase + '.tsv', cols.join('\t') + '\n' + out.map(o => cols.map(c => o[c]).join('\t')).join('\n') + '\n');

// Markdown: one table per set.
const md = [`# Pilot benchmark by condition`, ``, `Aggregated from \`${path.basename(summaryFile)}\` (${rows.length} runs) by \`benchmark/summarise.js\`. "rec" is recovered templates over templates in the condition's runs; "err %" is the mean internal error of the matched consensuses (ends excluded), "strict" counts ambiguity codes as errors; "ends" is the mean number of template bases missing at the ends per matched consensus; "hp" is the share of internal errors in homopolymer runs.`, ``];
for (const set of [...new Set(out.map(o => o.set))]) {
  md.push(`## ${set}`, ``, `| tool | condition | case | n | rec | unmatched | extended | dup | chimeric | err % | strict % | subs / ins / del | hp | ends | wall s | peak MB | alleles |`, `|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|`);
  for (const o of out.filter(o => o.set === set)) md.push(`| ${o.tool} | ${o.condition} | ${o.case_group} | ${o.reps} | ${o.recovered}/${o.truth} | ${o.unmatched} | ${o.extended} | ${o.duplicates} | ${o.chimeric} | ${o.err_internal_mean_pct} | ${o.err_internal_strict_mean_pct} | ${o.subs} / ${o.ins} / ${o.del} | ${o.hp_share_pct}${o.hp_share_pct === '' ? '' : ' %'} | ${o.term_del_per_consensus} | ${o.wall_s_mean} | ${o.peak_mb_max} | ${o.allele_classes} |`);
  md.push('');
}
fs.writeFileSync(outBase + '.md', md.join('\n'));
console.log(`${out.length} condition(s) from ${rows.length} runs -> ${outBase}.tsv, ${outBase}.md`);
