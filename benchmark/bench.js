// bench.js - run one tool over one benchmark set, time it, and score every case with evaluate.js.
//
// usage: node benchmark/bench.js --set S1 --tool olc|amplicon_sorter --runs DIR
//                                [--condition NAME] [--threads K] [--filter REGEX] [--limit N]
//                                [--force] [--reeval] [--dry-run]
//        node benchmark/bench.js --collect --runs DIR --results FILE
//
// Each case of DIR/<set>/manifest.tsv (made by generate-sets.js) is run in
// DIR/<set>/<tool>/<condition>/<case>/ with the tool's own output, stdout.txt, stderr.txt,
// time.json (wall seconds, CPU seconds, peak resident memory in MB, exit code, command line),
// consensus.fasta (a copy of the tool's final consensus file) and score.tsv / details.tsv from
// evaluate.js. A case with a time.json is skipped unless --force (rerun) or --reeval (score
// again). --collect concatenates every score.tsv under DIR into one table.
//
// Tools and conditions:
//   olc              default: --min-identity 90 (R10 cases) or 85 (R9 cases), --max-reads 20000
//   amplicon_sorter  default: -min 100 -maxr 100000 -ra -a;  r10: the same plus -ssg 94 -sc 98
// Linux tools are timed with /usr/bin/time -v; on Windows the peak working set and processor
// time are sampled every 250 ms through PowerShell.
'use strict';
const fs = require('fs');
const path = require('path');
const { spawn, spawnSync } = require('child_process');

const argv = process.argv.slice(2);
const opt = (k, d) => { const i = argv.indexOf(k); return i < 0 ? d : argv[i + 1]; };
const has = k => argv.includes(k);
const RUNS = opt('--runs', ''), SET = opt('--set', ''), TOOL = opt('--tool', ''), COND = opt('--condition', 'default');
const THREADS = opt('--threads', '4'), FILTER = opt('--filter', '') ? new RegExp(opt('--filter')) : null, LIMIT = +opt('--limit', 0);
const FORCE = has('--force'), REEVAL = has('--reeval'), DRY = has('--dry-run');
const root = path.join(__dirname, '..');
const evaluator = path.join(__dirname, 'evaluate.js');
const SOFTWARE = opt('--software', path.join(RUNS, '..', 'software'));
// Linux side (WSL): conda environments made by D:\OLC\software\setup-wsl.sh, GNU time, the NGSpeciesID binary.
const HOME = process.env.HOME || '';
const ENVS = opt('--conda-envs', path.join(HOME, 'miniforge3', 'envs'));
const TIME_BIN = process.env.TIME_BIN || (fs.existsSync('/usr/bin/time') ? '/usr/bin/time' : (fs.existsSync(path.join(HOME, 'miniforge3', 'bin', 'time')) ? path.join(HOME, 'miniforge3', 'bin', 'time') : 'time'));
const PYTHON = opt('--python', process.platform === 'win32' ? 'python' : (fs.existsSync(path.join(ENVS, 'ampsort', 'bin', 'python')) ? path.join(ENVS, 'ampsort', 'bin', 'python') : 'python3'));
const NGS_BIN = opt('--ngspeciesid', '');
const medakaModel = c => c.model === 'R9' ? 'r941_min_sup_g507' : 'r1041_e82_400bps_sup_v5.2.0';
const withPath = dirs => ({ PATH: dirs.concat(process.env.PATH || '').join(path.delimiter) });
// Manifests are written on Windows; under WSL the same files are reachable through /mnt/<drive>/.
function toLocal(p) {
  if (process.platform === 'win32' || !p) return p;
  const m = /^([A-Za-z]):[\\/](.*)$/.exec(p);
  return m ? '/mnt/' + m[1].toLowerCase() + '/' + m[2].replace(/\\/g, '/') : p;
}

const TOOLS = {
  olc: {
    conditions: { default: c => ['--min-identity', c.model === 'R9' ? '85' : '90'] },
    command(c, dir, cond) {
      return { cmd: process.execPath, args: [path.join(root, 'tools', 'cap-cli.js'), c.reads, '--out', path.join(dir, 'olc'), '--threads', THREADS, '--max-reads', '20000', '--quiet', ...cond(c)] };
    },
    consensus: dir => path.join(dir, 'olc.contigs.fasta'),
    pattern: '^contig',
  },
  amplicon_sorter: {
    conditions: { default: () => [], r10: () => ['-ssg', '94', '-sc', '98'] },
    command(c, dir, cond) {
      const script = path.join(SOFTWARE, 'amplicon_sorter', process.platform === 'win32' ? 'amplicon_sorter_single.py' : 'amplicon_sorter.py');
      const args = [script, '-i', c.reads, '-o', path.join(dir, 'as'), '-min', '100', '-maxr', '100000', '-ra', '-a', ...cond(c)];
      if (process.platform !== 'win32') args.push('-np', THREADS);
      return { cmd: PYTHON, args };
    },
    consensus: dir => findFile(path.join(dir, 'as'), /_consensussequences\.fasta$/) || findFile(path.join(dir, 'as'), /^consensusfile\.fasta$/),
    pattern: '',
  },
  // Linux only (WSL). default = medaka polishing with the model for the chemistry; racon = two racon
  // rounds; draft = the spoa consensus alone.
  ngspeciesid: {
    conditions: { default: c => ['--medaka', '--medaka_model', medakaModel(c)], racon: () => ['--racon'], draft: () => [] },
    command(c, dir, cond) {
      const bin = NGS_BIN || findFile(path.join(HOME, 'tools', 'ngspeciesid'), /^NGSpeciesID(-[\w.]+)?$/);
      return { cmd: bin, args: ['--ont', '--fastq', c.reads, '--outfolder', path.join(dir, 'ngs'), '--t', '1', '--seed', '1', '--consensus', ...cond(c)], env: withPath([path.join(ENVS, 'medaka', 'bin')]) };
    },
    consensus(dir) {
      const out = path.join(dir, 'ngs');
      for (const re of [/^medaka_cl_id_\d+$/, /^racon_cl_id_\d+$/]) {
        const files = fs.existsSync(out) ? fs.readdirSync(out).filter(n => re.test(n)).map(n => path.join(out, n, 'consensus.fasta')).filter(f => fs.existsSync(f)) : [];
        if (files.length) return files;
      }
      return fs.existsSync(out) ? fs.readdirSync(out).filter(n => /^consensus_reference_\d+\.fasta$/.test(n)).map(n => path.join(out, n)) : [];
    },
    pattern: '',
  },
  medaka_smolecule: {
    conditions: { default: () => [] },
    command(c, dir) {
      return { cmd: path.join(ENVS, 'medaka', 'bin', 'medaka'), args: ['smolecule', '--model', medakaModel(c), '--depth', '2', '--length', '100', '--threads', THREADS, path.join(dir, 'smol'), c.reads], env: withPath([path.join(ENVS, 'medaka', 'bin')]) };
    },
    consensus: dir => path.join(dir, 'smol', 'consensus.fasta'),
    pattern: '',
  },
  decona: {
    conditions: { default: () => [], medaka: () => ['-M'] },
    command(c, dir, cond) {
      // Decona takes a folder of FASTQ files and writes into ./data of the working directory.
      const inDir = path.join(dir, 'in');
      fs.rmSync(inDir, { recursive: true, force: true }); fs.rmSync(path.join(dir, 'data'), { recursive: true, force: true });
      fs.mkdirSync(inDir, { recursive: true }); fs.copyFileSync(c.reads, path.join(inDir, 'reads.fastq'));
      return { cmd: 'bash', args: [path.join(ENVS, 'decona', 'bin', 'decona'), '-i', inDir, '-l', '100', '-m', '20000', '-q', '7', '-c', '0.95', '-n', '5', '-k', '15', '-T', THREADS, ...cond(c)], env: withPath([path.join(ENVS, 'decona', 'bin')]) };
    },
    consensus: dir => findFile(path.join(dir, 'data'), /^medaka_fastas_.*\.fasta$/) || findFile(path.join(dir, 'data'), /^polished_fastas_.*\.fasta$/) || findFile(dir, /^(all_medaka_fastas|all_polished_fastas)\.fasta$/),
    pattern: '',
  },
};

function findFile(dir, re) {
  if (!fs.existsSync(dir)) return null;
  for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, e.name);
    if (e.isDirectory()) { const f = findFile(p, re); if (f) return f; }
    else if (re.test(e.name)) return p;
  }
  return null;
}
function readManifest(file) {
  const lines = fs.readFileSync(file, 'utf8').split(/\r?\n/).filter(l => l.trim());
  const cols = lines[0].split('\t');
  return lines.slice(1).map(l => { const f = l.split('\t'); const o = {}; cols.forEach((c, i) => o[c] = f[i] || ''); return o; });
}

function runTimed(cmd, args, dir, extraEnv) {
  return new Promise(resolve => {
    const outFd = fs.openSync(path.join(dir, 'stdout.txt'), 'w'), errFd = fs.openSync(path.join(dir, 'stderr.txt'), 'w');
    const t0 = Date.now();
    const env = Object.assign({}, process.env, extraEnv || {});
    let child, timeTxt = path.join(dir, 'time.txt');
    if (process.platform === 'win32') child = spawn(cmd, args, { cwd: dir, stdio: ['ignore', outFd, errFd], windowsHide: true, env });
    else child = spawn(TIME_BIN, ['-v', '-o', timeTxt, cmd, ...args], { cwd: dir, stdio: ['ignore', outFd, errFd], env });
    let peak = 0, cpu = 0, poller = null, pollerDone = Promise.resolve();
    if (process.platform === 'win32') {
      const script = `while($true){ try { $p = Get-Process -Id ${child.pid} -ErrorAction Stop; [Console]::Out.WriteLine("$($p.PeakWorkingSet64) $($p.TotalProcessorTime.TotalSeconds)") } catch { break }; Start-Sleep -Milliseconds 250 }`;
      poller = spawn('powershell', ['-NoProfile', '-NonInteractive', '-Command', script], { windowsHide: true, stdio: ['ignore', 'pipe', 'ignore'] });
      poller.stdout.on('data', d => { for (const l of d.toString().split(/\r?\n/)) { const [m, c] = l.trim().split(/\s+/); if (+m > peak) peak = +m; if (+c > cpu) cpu = +c; } });
      pollerDone = new Promise(res => { poller.on('exit', res); poller.on('error', res); });
    }
    child.on('error', e => { fs.writeFileSync(path.join(dir, 'stderr.txt'), String(e)); });
    child.on('exit', async code => {
      const wall = (Date.now() - t0) / 1000;
      fs.closeSync(outFd); fs.closeSync(errFd);
      await Promise.race([pollerDone, new Promise(res => setTimeout(res, 3000))]);
      if (poller) try { poller.kill(); } catch (e) { /* gone */ }
      const t = { wall_s: Math.round(wall * 100) / 100, exit: code === null ? -1 : code, cmd: [cmd, ...args].join(' ') };
      if (process.platform === 'win32') {
        t.peak_mb = Math.round(peak / 1048576); t.cpu_s = Math.round(cpu * 100) / 100;
        // A run shorter than the first memory sample: take the figure the OLC CLI prints itself.
        if (!t.peak_mb) { const m = /peak memory (\d+) MB/.exec(fs.readFileSync(path.join(dir, 'stderr.txt'), 'utf8')); if (m) t.peak_mb = +m[1]; }
      }
      else if (fs.existsSync(timeTxt)) {
        const txt = fs.readFileSync(timeTxt, 'utf8');
        const g = re => { const m = re.exec(txt); return m ? +m[1] : undefined; };
        t.peak_mb = Math.round((g(/Maximum resident set size \(kbytes\): (\d+)/) || 0) / 1024);
        t.cpu_s = Math.round(((g(/User time \(seconds\): ([\d.]+)/) || 0) + (g(/System time \(seconds\): ([\d.]+)/) || 0)) * 100) / 100;
      }
      resolve(t);
    });
  });
}

async function main() {
  if (!RUNS) { console.error('usage: node benchmark/bench.js --set S1 --tool olc --runs DIR [options]'); process.exit(1); }
  if (has('--collect')) {
    const out = opt('--results', path.join(RUNS, '..', 'results', 'summary.tsv'));
    const files = [];
    (function walk(d, depth) { if (depth > 5 || !fs.existsSync(d)) return; for (const e of fs.readdirSync(d, { withFileTypes: true })) { const p = path.join(d, e.name); if (e.isDirectory()) walk(p, depth + 1); else if (e.name === 'score.tsv') files.push(p); } })(RUNS, 0);
    let header = null, rows = [];
    for (const f of files.sort()) { const lines = fs.readFileSync(f, 'utf8').split(/\r?\n/).filter(l => l.trim()); if (!header) header = lines[0]; rows.push(...lines.slice(1)); }
    fs.mkdirSync(path.dirname(out), { recursive: true });
    fs.writeFileSync(out, (header || '') + '\n' + rows.join('\n') + '\n');
    console.log(`${rows.length} rows from ${files.length} runs -> ${out}`);
    return;
  }
  const tool = TOOLS[TOOL];
  if (!tool || !SET) { console.error('--set and --tool (' + Object.keys(TOOLS).join(', ') + ') are required'); process.exit(1); }
  const cond = tool.conditions[COND];
  if (!cond) { console.error(`${TOOL} has no condition "${COND}" (${Object.keys(tool.conditions).join(', ')})`); process.exit(1); }
  let cases = readManifest(path.join(RUNS, SET, 'manifest.tsv'));
  for (const c of cases) { c.reads = toLocal(c.reads); c.truth = toLocal(c.truth); c.variants = toLocal(c.variants); }
  if (FILTER) cases = cases.filter(c => FILTER.test(c.case));
  if (LIMIT) cases = cases.slice(0, LIMIT);
  console.log(`${SET}: ${cases.length} case(s), tool ${TOOL}, condition ${COND}`);
  let n = 0, failed = 0;
  for (const c of cases) {
    const dir = path.join(RUNS, SET, TOOL, COND, c.case);
    fs.mkdirSync(dir, { recursive: true });
    const timeFile = path.join(dir, 'time.json');
    const { cmd, args, env } = tool.command(c, dir, cond);
    if (DRY) { console.log(`${c.case}: ${[cmd, ...args].join(' ')}`); continue; }
    if (!cmd) { console.log(`${c.case}: the ${TOOL} executable was not found`); failed++; continue; }
    let t;
    if (fs.existsSync(timeFile) && !FORCE) {
      if (!REEVAL) { continue; }
      t = JSON.parse(fs.readFileSync(timeFile, 'utf8'));
    } else {
      for (const f of ['consensus.fasta', 'score.tsv', 'details.tsv']) if (fs.existsSync(path.join(dir, f))) fs.unlinkSync(path.join(dir, f));
      for (const sub of ['as', 'ngs', 'smol', 'data']) fs.rmSync(path.join(dir, sub), { recursive: true, force: true });
      t = await runTimed(cmd, args, dir, env);
      fs.writeFileSync(timeFile, JSON.stringify(t, null, 1));
    }
    const src = tool.consensus(dir);
    const consensus = path.join(dir, 'consensus.fasta');
    if (Array.isArray(src)) {
      // One file per cluster: join them, naming each record after its folder so names stay unique.
      let txt = '';
      for (const f of src) { const label = path.basename(path.dirname(f)) === path.basename(dir) || /consensus_reference/.test(path.basename(f)) ? path.basename(f, '.fasta') : path.basename(path.dirname(f)); txt += fs.readFileSync(f, 'utf8').replace(/^>/gm, '>' + label + '|'); if (!txt.endsWith('\n')) txt += '\n'; }
      fs.writeFileSync(consensus, txt);
    } else if (src && fs.existsSync(src)) fs.copyFileSync(src, consensus); else if (!fs.existsSync(consensus)) fs.writeFileSync(consensus, '');
    for (const f of ['score.tsv', 'details.tsv']) if (fs.existsSync(path.join(dir, f))) fs.unlinkSync(path.join(dir, f));
    const tags = [['set', SET], ['tool', TOOL], ['condition', COND], ['case', c.case], ['model', c.model], ['length', c.length], ['depth', c.depth], ['rep', c.rep], ['params', c.params]];
    const evArgs = [evaluator, '--consensus', consensus, '--truth', c.truth, '--time', timeFile, '--append', path.join(dir, 'score.tsv'), '--details', path.join(dir, 'details.tsv')];
    if (c.variants) evArgs.push('--variants', c.variants);
    if (tool.pattern) evArgs.push('--pattern', tool.pattern);
    for (const [k, v] of tags) evArgs.push('--tag', `${k}=${v}`);
    const ev = spawnSync(process.execPath, evArgs, { encoding: 'utf8' });
    const first = (ev.stdout || '').split('\n')[0] || ev.stderr;
    console.log(`${c.case}: exit ${t.exit}, ${t.wall_s} s, ${t.peak_mb === undefined ? '?' : t.peak_mb} MB | ${first.replace(/^consensus\.fasta: /, '')}`);
    if (t.exit !== 0) failed++;
    n++;
  }
  console.log(`done: ${n} run(s), ${failed} with a non-zero exit code`);
}
main();
