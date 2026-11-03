#!/usr/bin/env node
'use strict';
/* ─────────────────────────────────────────────────
   tools/cap-cli.js — command-line batch mode for the OLC Assembler (js/cap.js).

   The page assembles what fits in a text box: a few thousand reads at most. A
   Nanopore run is gigabytes of FASTQ, more than one JavaScript string can hold
   (V8 stops at about 512 million characters) and far more than an all-against-
   all assembler should be asked to align. This script is the bridge between the
   two: it STREAMS the file record by record, keeps only as many reads of each
   amplicon as consensus building can use, and runs the unchanged engine on each
   group in turn — hundreds of small assemblies instead of one impossible one.

     node tools/cap-cli.js reads.fastq --out run1 --bin-regex "barcode(\d+)" --threads 8

   What it does, in order:
     1. Reads FASTA or FASTQ (plain or .gz) as a stream — several files one
        after another, or a whole folder of them — applying the engine's
        own FASTQ rule — "@name", sequence on one or more lines, "+", quality up
        to the sequence length — so wrapped or indented files read correctly.
     2. Drops reads outside --min-length/--max-length, groups the rest by the
        capture group of --bin-regex matched against the header (one group when
        no regex is given), and keeps --max-reads of each group: the first N, or
        a uniform reservoir sample with --sample.
     3. Assembles each group with CapJS.assemble(). The pairwise overlap stage,
        which dominates the runtime, is spread over worker_threads with exactly
        the contract of the page's Worker pool (js/cap-pool.js), so the result is
        identical to a single-thread run; only the wall-clock differs.
     4. Writes, per group: <out>[.<bin>].contigs.fasta, .report.txt and .log.txt
        (the same artefacts the page offers), plus <out>.bins.tsv summarising
        every group.

   Memory: each group is held in memory only while it is assembled. If a group
   is large (thousands of long reads) and Node runs out of heap, raise it:
     node --max-old-space-size=16000 tools/cap-cli.js ...
   The engine's behaviour is not changed here in any way; js/cap.js is loaded as
   the page loads it.
   ───────────────────────────────────────────────── */

const fs = require('fs');
const os = require('os');
const path = require('path');
const zlib = require('zlib');
const readline = require('readline');
const { Worker, isMainThread, parentPort, workerData } = require('worker_threads');

const CAP_PATH = path.join(__dirname, '..', 'js', 'cap.js');
const REPORT_PATH = path.join(__dirname, '..', 'js', 'cap-report.js');

/* ───────────────────────── worker role ─────────────────────────
   The same file is the worker: it loads the engine and answers the two
   messages js/cap-worker.js answers — 'init' with the slimmed reads and the
   parameters, then 'shard' with a slice of candidate pairs to align. */
if (!isMainThread && workerData && workerData.role === 'cap-overlap-worker') {
  const CapJS = require(workerData.capPath);
  let reads = null, P = null, repeatThreshold = Infinity;
  parentPort.on('message', (msg) => {
    try {
      if (msg.type === 'init') {
        reads = msg.reads;
        P = CapJS.params(msg.params);
        repeatThreshold = typeof msg.repeatThreshold === 'number' ? msg.repeatThreshold : Infinity;
        parentPort.postMessage({ type: 'ready' });
      } else if (msg.type === 'shard') {
        if (!reads) throw new Error('shard received before init');
        const res = CapJS.computeOverlaps(reads, msg.pairs, P, repeatThreshold);
        parentPort.postMessage({ type: 'shard-done', id: msg.id, overlaps: res.overlaps, rejected: res.rejected });
      }
    } catch (err) {
      parentPort.postMessage({ type: 'error', id: msg.id, message: err && err.message ? err.message : String(err) });
    }
  });
  return;
}

const CapJS = require(CAP_PATH);
const CapReport = require(REPORT_PATH);

/* ───────────────────────── arguments ───────────────────────── */

const USAGE = `OLC Assembler — command-line batch mode (DigitalGens)

Usage
  node tools/cap-cli.js <reads.fastq|reads.fasta>[.gz] [more files or a folder ...] [options]

Input and grouping
  Several input files are read one after another. A folder stands for every FASTQ or FASTA
  file inside it (plain or .gz), in name order: a MinKNOW fastq_pass/barcodeNN folder, say.
  --out PREFIX          Output prefix (default: the input file name without its extension,
                        in the current folder). Files: PREFIX[.BIN].contigs.fasta,
                        .report.txt, .log.txt and PREFIX.bins.tsv.
  --bin-regex RE        Group reads by RE matched against the whole header line; the first
                        capture group (or the whole match) names the group, e.g.
                        "barcode(\\d+)" or "^(\\S+)_". Without it all reads form one group.
  --keep-unbinned       Assemble reads the regex does not match as group "unbinned"
                        (default: count them and leave them out).
  --max-reads N         Reads kept per group (default 1000; the engine accepts up to 20000).
                        Consensus of one amplicon stops improving after a few hundred reads.
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
  --no-chimera-check    Keep reads that look like two molecules joined in one read.
  --no-repeat-guard     Do not refuse joins across a repeat proven by divergence.
  --no-singletons       Leave unassembled reads out of the FASTA.
  --set KEY=VALUE       Any other engine parameter by name (repeatable), e.g. --set minSeeds=2.

Output
  --alignment           Also write PREFIX[.BIN].alignment.txt (reads under the consensus).
  --placements          Also write PREFIX[.BIN].placements.tsv.
  --quiet               No progress on stderr, only the final summary.
  --help, --version

Memory
  Only one group is in memory at a time. For groups of thousands of long reads raise the
  Node heap: node --max-old-space-size=16000 tools/cap-cli.js reads.fastq ...

Exit codes: 0 done; 1 bad arguments; 2 the input could not be read; 3 a group failed to
assemble (the others were still written).
`;

function fail(code, message) {
  process.stderr.write(message.replace(/\s+$/, '') + '\n');
  process.exit(code);
}

function siteVersion() {
  const read = p => { try { return fs.readFileSync(p, 'utf8'); } catch (e) { return ''; } };
  const theme = read(path.join(__dirname, '..', 'assets', 'js', 'theme.js'));
  if (theme) {
    const v = /DIGITALGENS_VERSION\s*=\s*"([^"]+)"/.exec(theme);
    const d = /DIGITALGENS_UPDATED\s*=\s*"([^"]+)"/.exec(theme);
    return 'DigitalGens v' + (v ? v[1] : '?') + (d ? ' (' + d[1] + ')' : '');
  }
  // The stand-alone package has no site around it: name and version of its package.json.
  try {
    const p = JSON.parse(read(path.join(__dirname, '..', 'package.json')));
    if (p.name && p.version) return p.name + ' v' + p.version;
  } catch (e) { /* no package.json either */ }
  return 'OLC Assembler';
}

function parseArgs(argv) {
  const o = {
    inputs: [], out: null, binRegex: null, keepUnbinned: false, maxReads: 1000,
    sample: false, seed: 1, minLength: 0, maxLength: Infinity, threads: null,
    threadsExplicit: false, engine: {}, alignment: false, placements: false, quiet: false
  };
  const num = (name, v, lo, hi) => {
    const n = Number(v);
    if (!Number.isFinite(n) || n < lo || n > hi) fail(1, `${name} needs a number between ${lo} and ${hi}, got "${v}".`);
    return n;
  };
  const takeValue = (i, name) => {
    if (i + 1 >= argv.length) fail(1, `${name} needs a value.`);
    return argv[i + 1];
  };
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    switch (a) {
      case '--help': case '-h': process.stdout.write(USAGE); process.exit(0); break;
      case '--version': process.stdout.write('OLC Assembler CLI, ' + siteVersion() + '\n'); process.exit(0); break;
      case '--out': o.out = takeValue(i, a); i++; break;
      case '--bin-regex':
        try { o.binRegex = new RegExp(takeValue(i, a)); } catch (e) { fail(1, `--bin-regex is not a valid regular expression: ${e.message}`); }
        i++; break;
      case '--keep-unbinned': o.keepUnbinned = true; break;
      case '--max-reads': o.maxReads = Math.round(num(a, takeValue(i, a), 2, 20000)); i++; break;
      case '--sample':
        o.sample = true;
        if (i + 1 < argv.length && /^\d+$/.test(argv[i + 1])) { o.seed = Number(argv[i + 1]); i++; }
        break;
      case '--min-length': o.minLength = num(a, takeValue(i, a), 0, 1e9); i++; break;
      case '--max-length': o.maxLength = num(a, takeValue(i, a), 1, 1e9); i++; break;
      case '--threads': o.threads = Math.round(num(a, takeValue(i, a), 1, 64)); o.threadsExplicit = true; i++; break;
      case '--min-identity': o.engine.minIdentity = num(a, takeValue(i, a), 50, 100); i++; break;
      case '--min-overlap': o.engine.minOverlap = num(a, takeValue(i, a), 10, 100000); i++; break;
      case '--max-overhang': o.engine.maxOverhang = num(a, takeValue(i, a), 0, 100000); i++; break;
      case '--kmer': o.engine.kmer = num(a, takeValue(i, a), 6, 15); i++; break;
      case '--min-read-len': o.engine.minReadLen = num(a, takeValue(i, a), 10, 100000); i++; break;
      case '--end-trim': o.engine.endTrim = num(a, takeValue(i, a), 0, 500); i++; break;
      case '--end-penalty': o.engine.endPenalty = num(a, takeValue(i, a), 0, 200); i++; break;
      case '--no-rc': o.engine.allowRC = false; break;
      case '--no-chimera-check': o.engine.chimeraCheck = false; break;
      case '--no-repeat-guard': o.engine.repeatGuard = false; break;
      case '--no-singletons': o.engine.includeSingletons = false; break;
      case '--set': {
        const kv = takeValue(i, a); i++;
        const m = /^([A-Za-z]\w*)=(.*)$/.exec(kv);
        if (!m) fail(1, `--set needs KEY=VALUE, got "${kv}".`);
        const v = m[2];
        o.engine[m[1]] = v === 'true' ? true : v === 'false' ? false : (v !== '' && Number.isFinite(Number(v)) ? Number(v) : v);
        break;
      }
      case '--alignment': o.alignment = true; break;
      case '--placements': o.placements = true; break;
      case '--quiet': case '-q': o.quiet = true; break;
      default:
        if (a[0] === '-') fail(1, `Unknown option ${a}. Run with --help.`);
        o.inputs.push(a);
    }
  }
  if (!o.inputs.length) fail(1, USAGE);
  if (o.maxLength < o.minLength) fail(1, '--max-length is smaller than --min-length.');
  if (!o.out) o.out = path.basename(path.resolve(o.inputs[0])).replace(/\.gz$/i, '').replace(/\.(fastq|fq|fasta|fa|fna|fas)$/i, '') || 'assembly';
  if (o.threads === null) o.threads = Math.max(1, Math.min(16, os.cpus().length - 1));
  return o;
}

/* ───────────────────────── streaming reader ─────────────────────────
   Yields {header, seq} per record without ever holding the file. The FASTQ
   branch is the engine's rule (capParseFastq in js/cap.js): the quality is read
   up to the length of the sequence and never looked at. A record cut off by the
   end of the file is reported and dropped; a record whose quality outruns its
   sequence is a corrupt file and stops the run. */
async function* readRecords(file, progress) {
  const raw = fs.createReadStream(file);
  const input = /\.gz$/i.test(file) ? raw.pipe(zlib.createGunzip()) : raw;
  const rl = readline.createInterface({ input, crlfDelay: Infinity });
  let mode = null;                 // 'fastq' | 'fasta' | 'plain', from the first non-blank line
  let header = null, seq = [], seqLen = 0, qualLen = 0, inQual = false;
  let lineNo = 0, plainCount = 0;
  const rec = () => ({ header, seq: seq.join('') });

  for await (const line of rl) {
    lineNo++;
    if (progress && (lineNo & 0xFFFF) === 0) progress(raw.bytesRead);
    const t = line.trim();         // also removes a byte-order mark and indentation
    if (mode === null) {
      if (!t) continue;
      mode = t[0] === '@' ? 'fastq' : t[0] === '>' ? 'fasta' : 'plain';
    }
    if (mode === 'fastq') {
      if (header === null) {
        if (!t) continue;
        if (t[0] !== '@') throw new Error(`line ${lineNo}: expected a FASTQ header starting with "@" (the four-line layout broke here)`);
        header = t.slice(1).trim(); seq = []; seqLen = 0; qualLen = 0; inQual = false;
        continue;
      }
      if (!inQual) {
        if (t[0] === '+') {
          inQual = true;
          if (seqLen === 0) { yield rec(); header = null; }     // zero-length read: no quality follows
        } else { seq.push(t); seqLen += t.length; }
        continue;
      }
      qualLen += t.length;
      if (qualLen >= seqLen) {
        if (qualLen > seqLen) throw new Error(`line ${lineNo}: the quality string of "${header.split(/\s+/)[0]}" is longer than its sequence; the file is not valid FASTQ`);
        yield rec(); header = null;
      }
    } else if (mode === 'fasta') {
      if (t[0] === '>') {
        if (header !== null) yield rec();
        header = t.slice(1).trim(); seq = [];
      } else if (t) {
        seq.push(t);
      }
    } else if (t) {
      yield { header: 'Seq' + (++plainCount), seq: t };
    }
  }
  if (header !== null) {
    if (mode === 'fasta') yield rec();
    else process.stderr.write(`warning: the last record "${header.split(/\s+/)[0]}" is cut off at the end of "${file}" and was dropped\n`);
  }
  if (progress) progress(raw.bytesRead, true);
}

/* ───────────────────────── input files ─────────────────────────
   A folder on the command line means every sequence file inside it, in natural
   name order (part_1, part_2, … part_10), which is how a basecaller leaves a
   run: fastq_pass/barcodeNN/ holding a few thousand reads per file. */
const SEQ_FILE = /\.(fastq|fq|fasta|fa|fna|fas)(\.gz)?$/i;

function expandInputs(list) {
  const files = [];
  for (const p of list) {
    let st;
    try { st = fs.statSync(p); } catch (e) { fail(2, `Cannot read "${p}": ${e.message}`); }
    if (st.isDirectory()) {
      const names = fs.readdirSync(p).filter(n => SEQ_FILE.test(n)).sort((a, b) => a.localeCompare(b, 'en', { numeric: true }));
      if (!names.length) fail(2, `No FASTQ or FASTA files in the folder "${p}".`);
      for (const n of names) files.push(path.join(p, n));
    } else {
      files.push(p);
    }
  }
  return files.map(file => ({ file, size: fs.statSync(file).size }));
}

/* ───────────────────────── grouping and sampling ─────────────────────────
   One reservoir per group: the first N reads fill it; after that each new read
   replaces a random slot with probability N/seen, which leaves a uniform sample
   of everything seen so far (Vitter's Algorithm R). Without --sample the
   reservoir simply stops accepting reads once full. */
function mulberry32(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function binKey(opts, header) {
  if (!opts.binRegex) return '';
  const m = opts.binRegex.exec(header);
  if (!m) return null;
  return m[1] !== undefined ? m[1] : m[0];
}

function safeName(bin) {
  return String(bin).replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 80) || 'group';
}

/* ───────────────────────── overlap runner on worker_threads ─────────────────────────
   Mirrors js/cap-pool.js: the reads go to every worker once, stripped to what
   the overlap stage uses (seq, len, mult, pile — the last two are what keep the
   repeat filter and its amplicon exemption identical to the main thread), and
   the candidate pairs are handed out in shards from a queue, so a slow shard
   cannot idle the pool. The layout stage sorts the overlaps into a total order
   before using them, which is why the order in which shards finish does not
   change the assembly. */
const MIN_PAIRS_FOR_POOL = 400;

function makePoolRunner(threads, explicit) {
  let workers = null;
  const minPairs = explicit ? 1 : MIN_PAIRS_FOR_POOL;

  function spawn() {
    const ws = [];
    for (let i = 0; i < threads; i++) {
      ws.push(new Worker(__filename, { workerData: { role: 'cap-overlap-worker', capPath: CAP_PATH } }));
    }
    return ws;
  }

  function runOnPool(reads, pairs, P, ctx, repeatThreshold) {
    const slim = reads.map(r => ({ idx: r.idx, seq: r.seq, len: r.len, mult: r.mult, pile: r.pile }));
    const chunk = Math.max(25, Math.min(2000, Math.ceil(pairs.length / (workers.length * 8))));
    const total = Math.ceil(pairs.length / chunk);
    return new Promise((resolve, reject) => {
      let next = 0, done = 0, settled = false;
      const overlaps = [], rejected = [];
      const finish = (err) => {
        if (settled) return;
        settled = true;
        if (err) reject(err); else resolve({ overlaps, rejected });
      };
      const feed = (w) => {
        if (settled) return;
        if (next >= total) { if (done >= total) finish(null); return; }
        const id = next++;
        const lo = id * chunk;
        w.postMessage({ type: 'shard', id, pairs: pairs.slice(lo, Math.min(pairs.length, lo + chunk)) });
      };
      for (const w of workers) {
        w.removeAllListeners('message');
        w.removeAllListeners('error');
        w.on('message', (m) => {
          if (m.type === 'ready') { feed(w); return; }
          if (m.type === 'error') { finish(new Error('worker: ' + m.message)); return; }
          if (m.type === 'shard-done') {
            for (const o of m.overlaps) overlaps.push(o);
            for (const r of m.rejected) rejected.push(r);
            done++;
            if (ctx && ctx.tick) ctx.tick('overlaps', Math.min(pairs.length, done * chunk), pairs.length).catch(() => {});
            if (done >= total) { finish(null); return; }
            feed(w);
          }
        });
        w.on('error', (err) => finish(err));
        w.postMessage({ type: 'init', reads: slim, params: P, repeatThreshold });
      }
    });
  }

  async function runner(reads, pairs, P, ctx, repeatThreshold) {
    if (!pairs.length) return { overlaps: [], rejected: [] };
    if (threads <= 1 || pairs.length < minPairs) {
      runner.threadsUsed = 1;
      return CapJS.serialOverlapRunner(reads, pairs, P, ctx, repeatThreshold);
    }
    if (!workers) workers = spawn();
    runner.threadsUsed = workers.length;
    try {
      return await runOnPool(reads, pairs, P, ctx, repeatThreshold);
    } catch (err) {
      // A worker that dies mid-run (out of memory) must not lose the group: redo
      // the stage on this thread, as the page does.
      process.stderr.write(`warning: worker pool failed (${err.message}); aligning on the main thread\n`);
      terminate();
      runner.threadsUsed = 1;
      return CapJS.serialOverlapRunner(reads, pairs, P, ctx, repeatThreshold);
    }
  }
  runner.threadsUsed = 1;

  function terminate() {
    if (workers) for (const w of workers) w.terminate();
    workers = null;
  }
  return { runner, terminate };
}

/* ───────────────────────── main ───────────────────────── */

function fmt(n) { return Number(n).toLocaleString('en-GB'); }

async function main() {
  const opts = parseArgs(process.argv.slice(2));
  const t0 = Date.now();
  const say = opts.quiet ? () => {} : (s) => process.stderr.write(s + '\n');

  const inputs = expandInputs(opts.inputs);
  const size = inputs.reduce((s, f) => s + f.size, 0);
  if (!size) fail(2, inputs.length === 1 ? `"${inputs[0].file}" is empty.` : 'The input files are empty.');
  const gz = inputs.some(f => /\.gz$/i.test(f.file));

  // 1. Stream, filter, group, sample.
  const bins = new Map();           // key -> {reads: [{header, seq}], seen, kept}
  const rnd = mulberry32(opts.seed);
  let total = 0, tooShort = 0, tooLong = 0, unbinned = 0, lastPct = -1;
  const progress = (bytes, final) => {
    if (opts.quiet) return;
    const pct = Math.min(100, Math.floor(bytes / size * 100));
    if (final || pct >= lastPct + 5) {
      lastPct = pct;
      process.stderr.write(`reading ${gz ? '(gzip) ' : ''}${pct}%  ${fmt(total)} reads\r`);
    }
  };
  const shown = inputs.length === 1 ? inputs[0].file
              : inputs.length + ' files, ' + inputs[0].file + ' … ' + path.basename(inputs[inputs.length - 1].file);
  say(`${siteVersion()} — OLC Assembler CLI\ninput  ${shown}  (${(size / 1048576).toFixed(1)} MB${gz ? ', gzip' : ''})`);
  let doneBytes = 0, current = '';
  try {
    for (const f of inputs) {
    current = f.file;
    const last = f === inputs[inputs.length - 1];
    for await (const r of readRecords(f.file, (bytes, final) => progress(doneBytes + bytes, final && last))) {
      total++;
      const len = r.seq.length;
      if (len < opts.minLength) { tooShort++; continue; }
      if (len > opts.maxLength) { tooLong++; continue; }
      let key = binKey(opts, r.header);
      if (key === null) {
        if (!opts.keepUnbinned) { unbinned++; continue; }
        key = 'unbinned';
      }
      let b = bins.get(key);
      if (!b) { b = { reads: [], seen: 0 }; bins.set(key, b); }
      b.seen++;
      if (b.reads.length < opts.maxReads) b.reads.push(r);
      else if (opts.sample) {
        const j = Math.floor(rnd() * b.seen);
        if (j < opts.maxReads) b.reads[j] = r;
      }
    }
    doneBytes += f.size;
    }
  } catch (e) {
    fail(2, `\nCannot read "${current}": ${e.message}`);
  }
  if (!opts.quiet) process.stderr.write('\n');
  say(`read ${fmt(total)} records in ${((Date.now() - t0) / 1000).toFixed(1)} s` +
      (tooShort || tooLong ? `; dropped ${fmt(tooShort)} shorter than ${opts.minLength} nt and ${fmt(tooLong)} longer than ${opts.maxLength} nt` : '') +
      (opts.binRegex ? `; ${fmt(bins.size)} group(s) by ${opts.binRegex}` + (unbinned ? `, ${fmt(unbinned)} unmatched read(s) left out` : '') : ''));
  if (!bins.size) fail(2, 'No reads to assemble.');

  // 2. Assemble each group.
  const outDir = path.dirname(path.resolve(opts.out));
  fs.mkdirSync(outDir, { recursive: true });
  const pool = makePoolRunner(opts.threads, opts.threadsExplicit);
  const rows = [['group', 'reads_seen', 'reads_used', 'contigs', 'longest_bp', 'total_consensus_bp', 'singletons', 'threads', 'seconds', 'status']];
  let failed = 0, k = 0;
  const keys = [...bins.keys()].sort((a, b) => a.localeCompare(b, 'en', { numeric: true }));
  for (const key of keys) {
    const b = bins.get(key);
    k++;
    const label = key === '' ? '' : key;
    const prefix = opts.out + (label === '' ? '' : '.' + safeName(label));
    const t1 = Date.now();
    const fasta = b.reads.map(r => '>' + r.header + '\n' + r.seq + '\n').join('');
    b.reads = null;                                   // the FASTA text is all the engine needs
    const params = Object.assign({}, opts.engine, { input: fasta, maxReads: Math.max(2, opts.maxReads) });
    let lastStage = '', lastShown = 0;
    const onProgress = (stage, i, n) => {
      if (opts.quiet) return;
      const now = Date.now();
      if (stage !== lastStage) { lastStage = stage; lastShown = now; process.stderr.write(`  [${k}/${keys.length}] ${label || 'all reads'}: ${stage}${n ? '' : '…'}\r`); return; }
      if (n && now - lastShown > 500) { lastShown = now; process.stderr.write(`  [${k}/${keys.length}] ${label || 'all reads'}: ${stage} ${Math.floor(i / n * 100)}%\r`); }
    };
    let status = 'ok', result = null;
    try {
      result = await CapJS.assemble(params, { overlapRunner: pool.runner, onProgress, threads: pool.runner.threadsUsed });
      if (!result || !result.ok) status = 'failed: ' + ((result && result.message) || 'no result');
    } catch (e) {
      status = 'failed: ' + e.message;
    }
    const secs = (Date.now() - t1) / 1000;
    if (!opts.quiet) process.stderr.write(''.padEnd(78) + '\r');
    if (status !== 'ok') {
      failed++;
      say(`  [${k}/${keys.length}] ${label || 'all reads'}: ${status}`);
      rows.push([label || 'all', b.seen, 0, 0, 0, 0, 0, pool.runner.threadsUsed, secs.toFixed(1), status]);
      continue;
    }
    result.stats.threads = pool.runner.threadsUsed;
    const includeSingletons = opts.engine.includeSingletons !== false;
    result.params.includeSingletons = includeSingletons;
    const contigs = result.contigs || [];
    const longest = contigs.reduce((m, c) => Math.max(m, (c.consensus || '').length), 0);
    const totalBp = contigs.reduce((m, c) => m + (c.consensus || '').length, 0);
    fs.writeFileSync(prefix + '.contigs.fasta', CapReport.fasta(result, { includeSingletons }));
    fs.writeFileSync(prefix + '.report.txt', CapReport.text(result));
    fs.writeFileSync(prefix + '.log.txt', CapReport.log(result).replace('(Web Workers)', '(worker threads)'));
    if (opts.alignment) fs.writeFileSync(prefix + '.alignment.txt', CapReport.alignmentText(result, { width: 100 }));
    if (opts.placements) fs.writeFileSync(prefix + '.placements.tsv', CapReport.placements(result));
    const used = (result.reads || []).length;
    say(`  [${k}/${keys.length}] ${label || 'all reads'}: ${fmt(used)} of ${fmt(b.seen)} reads → ${contigs.length} contig(s), longest ${fmt(longest)} bp, ` +
        `${(result.singletons || []).length} singleton(s), ${secs.toFixed(1)} s on ${pool.runner.threadsUsed} thread(s)`);
    rows.push([label || 'all', b.seen, used, contigs.length, longest, totalBp, (result.singletons || []).length, pool.runner.threadsUsed, secs.toFixed(1), status]);
  }
  pool.terminate();

  fs.writeFileSync(opts.out + '.bins.tsv', rows.map(r => r.join('\t')).join('\n') + '\n');
  const rss = process.resourceUsage().maxRSS / 1024;   // kilobytes -> MB
  process.stderr.write(`done: ${keys.length} group(s), ${failed ? failed + ' failed, ' : ''}${((Date.now() - t0) / 1000).toFixed(1)} s, peak memory ${rss.toFixed(0)} MB → ${opts.out}.bins.tsv\n`);
  process.exit(failed ? 3 : 0);
}

main().catch(e => fail(2, e && e.stack ? e.stack : String(e)));
