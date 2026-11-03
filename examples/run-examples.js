// run-examples.js - assemble the bundled examples and compare the contigs with examples/*/expected/.
//
// usage: node examples/run-examples.js [--threads K]
//
// Writes into examples/synthetic/out/ and examples/nanopore/out/ (ignored by git). The contigs
// must come out byte-identical to the stored ones whatever the thread count; the reports differ
// only in the elapsed time and the number of threads.
'use strict';
const { spawnSync } = require('child_process');
const fs = require('fs');
const path = require('path');

const root = path.join(__dirname, '..');
const cli = path.join(root, 'tools', 'cap-cli.js');
const threads = (() => { const i = process.argv.indexOf('--threads'); return i < 0 ? [] : ['--threads', process.argv[i + 1]]; })();

const examples = [
  {
    name: 'synthetic R10 reads, three 2 kb amplicons, defaults',
    dir: 'synthetic', input: 'reads-r10.fastq.gz', prefix: 'r10',
    args: ['--bin-regex', 'barcode(\\d+)'],
    outputs: ['01', '02', '03'],
  },
  {
    name: 'synthetic R9 reads, the same amplicons, --min-identity 85',
    dir: 'synthetic', input: 'reads-r9.fastq.gz', prefix: 'r9',
    args: ['--bin-regex', 'barcode(\\d+)', '--min-identity', '85'],
    outputs: ['01', '02', '03'],
  },
  {
    name: 'Oxford Nanopore R9.4.1 run, barcode 79, first 1000 reads of 200 nt or more',
    dir: 'nanopore', input: 'FAS90246_pass_barcode79_cfafe70b_765ba3f3_1.fastq.gz', prefix: 'barcode79',
    args: ['--bin-regex', 'barcode(\\d+)', '--min-identity', '85', '--min-length', '200', '--no-singletons'],
    outputs: ['79'],
  },
];

let failed = 0;
for (const ex of examples) {
  const dir = path.join(root, 'examples', ex.dir);
  const out = path.join(dir, 'out');
  fs.mkdirSync(out, { recursive: true });
  const args = [cli, path.join(dir, ex.input), '--out', path.join(out, ex.prefix), ...ex.args, ...threads, '--quiet'];
  process.stdout.write(`\n== ${ex.name}\n   node ${args.map(a => /[\s()\\]/.test(a) ? '"' + a + '"' : a).join(' ')}\n`);
  const t0 = Date.now();
  const r = spawnSync(process.execPath, args, { stdio: ['ignore', 'inherit', 'inherit'] });
  if (r.status !== 0) { console.log(`   FAILED with exit code ${r.status}`); failed++; continue; }
  console.log(`   ${((Date.now() - t0) / 1000).toFixed(1)} s`);
  for (const bin of ex.outputs) {
    const file = `${ex.prefix}.${bin}.contigs.fasta`;
    const got = fs.readFileSync(path.join(out, file));
    const want = fs.readFileSync(path.join(dir, 'expected', file));
    const same = got.equals(want);
    console.log(`   ${file}: ${same ? 'identical to expected' : 'DIFFERS from expected'}`);
    if (!same) failed++;
  }
  if (ex.dir === 'synthetic') {
    for (const bin of ex.outputs) {
      const c = spawnSync(process.execPath, [path.join(dir, 'check.js'), path.join(out, `${ex.prefix}.${bin}.contigs.fasta`)], { encoding: 'utf8' });
      process.stdout.write(c.stdout.split('\n').filter(l => l.trim()).map(l => '   ' + l + '\n').join(''));
    }
  }
}
console.log(failed ? `\n${failed} problem(s)` : '\nall examples reproduce the expected contigs');
process.exit(failed ? 1 : 0);
