'use strict';
/* ════════════════════════════════════════════════════════════════════════════
   cap.js — OLC-JS, the engine of the OLC Assembler: an overlap–layout–consensus DNA sequence assembler.

   Assembles one or more contigs from a set of short, possibly repeated,
   possibly imperfect fragments (Sanger reads, cloned inserts, NGS contigs,
   synthetic oligos, Oxford Nanopore R9/R10 amplicon reads from FASTQ). The design follows CAP3 (Huang & Madan, Genome Res 9:868,
   1999): k-mer seeding → banded overlap alignment → overhang-filtered overlap
   acceptance → greedy layout → column-wise multiple alignment → consensus.

   Two things matter for the intended use and drive most of the design:

     • Reads may be REPEATED. A read pair that shares only an interspersed
       repeat aligns over the repeat and then diverges, leaving a large
       unmatched overhang on both sides. Rejecting overlaps by *overhang*
       (not just by length or identity) is what keeps repeats from collapsing
       into one contig — this is CAP3's key filter and it is implemented in
       acceptOverlap() below.

     • Reads may DIFFER. Divergence has two causes that must be told apart:
       random sequencing error (isolated, unreproducible) and real
       polymorphism (a minority allele seen in several independent reads).
       The consensus caller separates them by minor-allele fraction AND
       minor-allele read count, and emits an IUPAC ambiguity code (R, Y, S,
       W, K, M, B, D, H, V, N) at genuinely polymorphic columns.

   The file is DOM-free on purpose: the same source runs on the main thread,
   inside a Web Worker via importScripts (js/cap-worker.js), and under Node.js
   for testing. Nothing here touches window or document.

   Coordinates are 0-based, half-open [start, end) throughout.
   Sequences are lower-case IUPAC letters as produced by DNA() in js/dna.js.
   ════════════════════════════════════════════════════════════════════════════ */

/* ───────────────────────────── 0. Constants & tables ─────────────────────── */

// IUPAC letter -> 4-bit set {A=1, C=2, G=4, T=8}. Used both for match scoring
// (an R is a partial match to A and to G) and for consensus encoding.
var CAP_BITS = {
    a: 1, c: 2, g: 4, t: 8, u: 8,
    r: 5, y: 10, s: 6, w: 9, k: 12, m: 3,
    b: 14, d: 13, h: 11, v: 7, n: 15
};
// 4-bit set -> IUPAC letter. Index 0 is unused (empty set); it maps to 'n' so a
// degenerate call can never produce an out-of-alphabet character.
var CAP_FROM_BITS = ['n', 'a', 'c', 'm', 'g', 'r', 's', 'v',
                     't', 'w', 'y', 'h', 'k', 'd', 'b', 'n'];

var CAP_NEG = -1073741824;              // -2^30: "impossible" DP cell
var CAP_HALF_NEG = CAP_NEG / 2;         // anything below this is still impossible

// Default parameters. Every entry is overridable from the UI; the comments give
// the reasoning so the defaults can be argued with rather than guessed at.
var CAP_DEFAULTS = {
    // ── input handling ──
    minReadLen: 40,      // reads shorter than this cannot support a minOverlap join
    endTrim: 0,          // fixed bases clipped from both ends of every read
    qualClip: 0,         // PHRED floor for end clipping (0 = off; needs quality input)
    maxReads: 5000,      // refuse absurd inputs rather than hanging the tab

    // ── seeding ──
    kmer: 12,            // P(error-free 12-mer) = 0.99^12 ≈ 0.89 at 1% error, 0.28 at 10%;
                         // a 100 nt overlap still yields ~25 exact seeds at 10% error
    // A pair needs this many seeds on one diagonal band before it is worth
    // aligning. Chance k-mer collisions between unrelated reads are common —
    // a million indexed 12-mers collide tens of thousands of times — and each
    // one costs a full banded alignment, so the threshold is what keeps the
    // overlap stage proportional to the real overlaps rather than to the square
    // of the input. 3 was measured to give results identical to 2 on clean data,
    // on a 300-read repeat-rich set, AND on an 8%-error set, while cutting the
    // work. Lower it to 2 for extremely noisy input, where few words survive.
    minSeeds: 3,
    maxSeedOcc: 100,     // k-mers occurring more often than this are repeat words: skipped
                         // for seeding but counted toward a read's repetitiveness score
    maxCandidates: 400000,   // hard cap on aligned pairs; truncation is reported

    // ── overlap acceptance (CAP3-style) ──
    minOverlap: 40,      // shortest alignment accepted as evidence of a join
    minIdentity: 90,     // percent identity over aligned columns
    maxOverhang: 30,     // KEY REPEAT FILTER: unmatched bases allowed at a join's
                         // outer edges. A repeat-only overlap has large overhangs
                         // on both sides and is rejected here.
    bandWidth: 32,       // half-width of the DP band around the seed diagonal
    maxBandWidth: 256,   // band is doubled and retried when the path hugs the edge

    // ── scoring ──
    matchScore: 2,
    mismatchScore: -3,
    gapOpen: -5,         // one-off cost of opening a gap …
    gapExtend: -2,       // … plus this per gap position

    // ── layout ──
    offsetTolerance: 20, // an overlap disagreeing with the layout by more than this
                         // many bases is a conflict, not a join
    allowRC: true,       // consider reverse-complemented reads
    chimeraCheck: true,  // set aside reads that carry two molecules (capChimeras, capPartingChimeras),
                         // trim fold-back reads (capFoldbacks) and refuse overlaps
                         // that run across a chimeric junction (capOverlapPair)
    // A read whose last (or first) bases are the reverse complement of sequence
    // elsewhere in the same read is a fold-back read. The copy must be at least
    // this long; a shorter one is left to the overhang filter and to the
    // contig-end trim, which between them contain the damage a short tail can do.
    // 0 = off. Measured on two 2000 nt products with 1-10 % chimeric reads: the
    // tails were 600-1400 nt and every one was found; see capFoldbacks.
    foldbackMin: 50,
    // Local identity, percent over CAP_JUNCTION_WINDOW (50) aligned columns,
    // below which an overlap is refused as running across a chimeric junction
    // even when its overall identity passes. 0 = minIdentity − 35, not below
    // 40. Measured: in 32,722 genuine overlaps between R9-level reads (7 %
    // error each, minIdentity 85) the worst 50-column window of an overlap was
    // 54 % once and under 60 % in one overlap per thousand; between R10-level
    // reads (minIdentity 90) it never went under 72 %. Between two chimeras of
    // the same two products the stretch between their junctions aligned at
    // 32-47 %. Real reads carry noisier stretches than the model: on the
    // bundled R9.4.1 run (minIdentity 85) three accepted overlaps of 427-1151
    // nt had a worst window of 52-54 % and one of 46 %, so the margin is 35
    // points rather than the 30 the model would allow.
    junctionIdentity: 0,
    // A contig end carried by fewer reads than this — or than 5 % of the
    // contig's median coverage, whichever is larger — over more than maxOverhang
    // bases is cut back to the last well-supported column: it is one or two
    // reads' unmatched tail, not the product. The maxOverhang tolerance leaves
    // the natural taper of read ends alone (reads of one product end within a
    // few bases of one another, and the last base may be carried by few of
    // them), and a contig whose own median coverage is under twice the
    // threshold is left alone (see capTrimEnds). 0 = off.
    minEndCov: 3,
    repeatGuard: true,   // stop a contig at a repeat boundary rather than guessing
                         // which copy continues (see capBranchGuard)
    // These two defaults are measured, not guessed. On a 30 kb target carrying a
    // 400 bp exact repeat at five loci, sequenced as 300 reads of 600-900 nt at
    // 7.5x with 1.5% error, they take the accepted-overlap set from 152 wrong
    // joins (of 2390) to ZERO — all 2162 accepted overlaps agree with the known
    // truth — and assemble the whole target into a single 29889 bp contig at
    // 99.92% identity. Loosening repeatFrac to 0.85 restores 34 wrong joins;
    // tightening repeatFactor to 1.6 with repeatFrac 0.4 starts splitting
    // contigs that are correct.
    repeatFactor: 1.8,   // word copy number, as a multiple of the median, above which
                         // a base counts as repeat sequence
    repeatFrac: 0.6,     // an overlap whose aligned region is at least this fraction
                         // repeat on BOTH reads is refused as unresolvable

    // ── consensus & polymorphism ──
    polyMinFrac: 25,     // minor allele must reach this % of the column …
    polyMinCount: 2,     // … in at least this many reads …
    polyMinCov: 3,       // … and the column must have at least this coverage.
                         // Together these keep a 1–2% error rate from being called
                         // as polymorphism while a true 50/50 het site is caught.
    polyMaxP: 0,         // significance level for the binomial test; 0 = auto
                         // (1/total columns, i.e. ~1 false call per assembly)
    qualWeight: true,    // weight base votes by PHRED quality when quality is supplied
    endPenalty: 0,       // down-weight this many bases at each read end (0 = off)
    includeSingletons: true
};

/* ───────────────────────────── 1. Small utilities ────────────────────────── */

var CAP_COMP = { a: 't', c: 'g', g: 'c', t: 'a', u: 'a', r: 'y', y: 'r', s: 's',
                 w: 'w', k: 'm', m: 'k', b: 'v', d: 'h', h: 'd', v: 'b', n: 'n' };

// Reverse complement. Deliberately local rather than js/dna.js's ComplementDNA()
// so the engine can be loaded in a Worker or in Node without pulling dna.js in.
function capRC(s) {
    var out = new Array(s.length);
    for (var i = 0, j = s.length - 1; j >= 0; i++, j--) {
        var c = CAP_COMP[s.charAt(j)];
        out[i] = c === undefined ? 'n' : c;
    }
    return out.join('');
}

// Merge user parameters over the defaults, clamping each to a sane range so a
// stray value from the UI cannot produce a pathological run.
function capParams(p) {
    var q = {}, k;
    for (k in CAP_DEFAULTS) if (Object.prototype.hasOwnProperty.call(CAP_DEFAULTS, k)) q[k] = CAP_DEFAULTS[k];
    if (p) for (k in p) if (Object.prototype.hasOwnProperty.call(p, k) && p[k] !== undefined && p[k] !== null && p[k] !== '') q[k] = p[k];

    var num = function (name, lo, hi) {
        var v = typeof q[name] === 'number' ? q[name] : parseFloat(q[name]);
        if (!isFinite(v)) v = CAP_DEFAULTS[name];
        q[name] = Math.min(hi, Math.max(lo, v));
    };
    num('minReadLen', 10, 100000);
    num('endTrim', 0, 500);
    num('qualClip', 0, 93);
    num('maxReads', 2, 20000);
    num('kmer', 6, 15);              // >15 does not fit the 2-bit packing in an int32
    num('minSeeds', 1, 100);
    num('maxSeedOcc', 2, 100000);
    num('maxCandidates', 1000, 20000000);
    num('minOverlap', 10, 100000);
    num('minIdentity', 50, 100);
    num('maxOverhang', 0, 100000);
    num('bandWidth', 4, 4096);
    num('maxBandWidth', 8, 8192);
    num('matchScore', 1, 100);
    num('mismatchScore', -100, -1);
    num('gapOpen', -100, 0);
    num('gapExtend', -100, -1);
    num('offsetTolerance', 0, 10000);
    num('repeatFactor', 1.2, 100);
    num('repeatFrac', 0.1, 1);
    num('polyMinFrac', 1, 100);
    num('polyMinCount', 1, 1000);
    num('polyMinCov', 1, 1000);
    num('endPenalty', 0, 200);
    num('polyMaxP', 0, 1);
    num('foldbackMin', 0, 100000);
    num('junctionIdentity', 0, 100);
    num('minEndCov', 0, 100000);

    q.kmer = Math.round(q.kmer);
    q.minSeeds = Math.round(q.minSeeds);
    if (q.maxBandWidth < q.bandWidth) q.maxBandWidth = q.bandWidth;
    q.allowRC = !!q.allowRC;
    q.chimeraCheck = !!q.chimeraCheck;
    q.repeatGuard = !!q.repeatGuard;
    q.qualWeight = !!q.qualWeight;
    q.includeSingletons = !!q.includeSingletons;
    return q;
}

/* ───────────────────────────── 2. Input parsing ──────────────────────────── */

// FASTQ is read by shape, not by line prefix: a record is "@name", the sequence
// on one or more lines, "+" (optionally repeating the name) and quality lines
// up to the length of the sequence. The length is the only safe delimiter,
// because a quality line may itself begin with "@", "+" or ">" — all legal
// PHRED characters. Indentation and a byte-order mark are ignored (trim()
// removes both). The quality is not read: the Quality tab takes PHRED values
// as numbers, and nothing here decodes the FASTQ encoding. The names are kept.
function capIsFastq(lines) {
    var i = 0, j, c;
    while (i < lines.length && !lines[i].trim()) i++;
    if (i >= lines.length || lines[i].trim().charAt(0) !== '@') return false;
    // The "+" line must follow within the first record; a second "@" (or a ">")
    // before it means this is not FASTQ. Quality lines only come after the "+".
    for (j = i + 1; j < lines.length && j < i + 400; j++) {
        c = lines[j].trim().charAt(0);
        if (c === '+') return true;
        if (c === '@' || c === '>') return false;
    }
    return false;
}

// Returns the records read; `recs.stoppedAt` is set to the 1-based line number
// of the record where the layout broke, so the caller can say what was not read.
function capParseFastq(lines) {
    var recs = [], i = 0, n = lines.length;
    while (i < n) {
        var head = lines[i].trim();
        if (!head) { i++; continue; }
        if (head.charAt(0) !== '@') { recs.stoppedAt = i + 1; break; }
        var start = i, seq = [], seqLen = 0, qualLen = 0, s;
        for (i++; i < n && lines[i].trim().charAt(0) !== '+'; i++) {
            s = lines[i].trim();
            seq.push(s); seqLen += s.length;
        }
        if (i >= n) { recs.stoppedAt = start + 1; break; }          // no "+" line: a truncated record
        for (i++; i < n && qualLen < seqLen; i++) qualLen += lines[i].trim().length;
        if (qualLen > seqLen) { recs.stoppedAt = start + 1; break; } // quality longer than the sequence: not FASTQ
        head = head.substring(1).trim();
        recs.push({ name: head.split(/\s+/)[0] || ('Seq' + (recs.length + 1)),
                    desc: head, body: seq.join('') });
    }
    return recs;
}

// Minimal FASTA/FASTQ/tab/plain reader. The site's ReadingSeq() (js/readseq.js)
// is not used here because the engine must also load inside a Worker and under
// Node, and because we need per-record parsing of the optional quality file with
// the same record splitting. Behaviour matches ReadingSeq() for the formats
// users actually paste: multi-FASTA, "name<TAB>sequence" lines, or a bare
// sequence; FASTQ, four-line or wrapped, is read as well (see capParseFastq).
function capParseFasta(text) {
    var recs = [];
    if (!text) return recs;
    var t = String(text).replace(/\r\n?/g, '\n');
    var fqLines = t.split('\n');
    if (capIsFastq(fqLines)) return capParseFastq(fqLines);
    // Text that opens with "@" but never reaches a "+" line is FASTQ that cannot
    // be read (or no sequence at all). Refuse it rather than make a record of
    // every line, quality strings included; capParseReads explains.
    var f = 0;
    while (f < fqLines.length && !fqLines[f].trim()) f++;
    if (f < fqLines.length && fqLines[f].trim().charAt(0) === '@') {
        recs.stoppedAt = f + 1;
        recs.fastqLike = true;
        return recs;
    }
    if (t.indexOf('>') !== -1) {
        var blocks = t.split('\n>');
        for (var i = 0; i < blocks.length; i++) {
            var b = blocks[i];
            if (i === 0) {
                var p = b.indexOf('>');
                if (p === -1) continue;              // text before the first ">"
                b = b.substring(p + 1);
            }
            var nl = b.indexOf('\n');
            var head = (nl === -1 ? b : b.substring(0, nl)).trim();
            var body = nl === -1 ? '' : b.substring(nl + 1);
            if (!head && !body.trim()) continue;
            recs.push({ name: head.split(/\s+/)[0] || ('Seq' + (recs.length + 1)),
                        desc: head, body: body });
        }
        return recs;
    }
    var lines = t.split('\n');
    var anyTab = t.indexOf('\t') !== -1;
    for (var L = 0; L < lines.length; L++) {
        var line = lines[L].trim();
        if (!line) continue;
        if (anyTab) {
            var tab = line.indexOf('\t');
            if (tab === -1) continue;
            recs.push({ name: line.substring(0, tab).trim() || ('Seq' + (recs.length + 1)),
                        desc: line.substring(0, tab).trim(), body: line.substring(tab + 1) });
        } else {
            recs.push({ name: 'Seq' + (recs.length + 1), desc: '', body: line });
        }
    }
    return recs;
}

// Keep only IUPAC nucleotide letters, lower-cased; u→t and i→g as the rest of
// the suite does. Everything else (digits, gaps, whitespace) is dropped.
function capCleanDNA(s) {
    var out = [], i, c;
    s = String(s).toLowerCase();
    for (i = 0; i < s.length; i++) {
        c = s.charAt(i);
        if (c === 'u') c = 't';
        else if (c === 'i') c = 'g';
        if (CAP_BITS[c] !== undefined) out.push(c);
    }
    return out.join('');
}

/**
 * Parse reads (and optional PHRED qualities) and apply end clipping.
 *
 * @param {string} text     FASTA / tab-delimited / bare sequence input
 * @param {string} qualText optional FASTA-quality text (whitespace-separated ints)
 * @param {object} P        parameters from capParams()
 * @returns {{reads:Array, notes:Array<string>, dropped:Array}}
 */
function capParseReads(text, qualText, P) {
    var recs = capParseFasta(text);
    var notes = [], dropped = [], reads = [];
    if (recs.fastqLike) {
        notes.push('The text starts with "@" like FASTQ, but no record could be read. A record is "@name", ' +
                   'the sequence (one or more lines), a "+" line and a quality string of the same length.');
    } else if (recs.stoppedAt) {
        notes.push('FASTQ input: ' + recs.length + ' record(s) were read; from line ' + recs.stoppedAt +
                   ' the text no longer followed the FASTQ layout ("@name", sequence, "+", quality of the same length) ' +
                   'and was not read.');
    }

    // Optional quality file, matched to reads by name first, by order as a fallback.
    var qualByName = null, qualByOrder = null;
    if (qualText && String(qualText).trim()) {
        var qrecs = capParseFasta(qualText);
        qualByName = {};
        qualByOrder = [];
        for (var qi = 0; qi < qrecs.length; qi++) {
            var nums = String(qrecs[qi].body).trim().split(/[^0-9]+/);
            var arr = [];
            for (var qj = 0; qj < nums.length; qj++) {
                if (nums[qj] === '') continue;
                var v = parseInt(nums[qj], 10);
                arr.push(isFinite(v) ? Math.min(93, Math.max(0, v)) : 0);
            }
            qualByOrder.push(arr);
            // Matching by name is only safe while the names are unique. A repeated
            // name would otherwise silently resolve to whichever record came last
            // and attach one read's quality values to a different read — worse
            // than having no qualities at all, because nothing would look wrong.
            // Duplicated names are dropped from the map so those reads fall back
            // to matching by order.
            var qn = qrecs[qi].name;
            if (qn) {
                if (Object.prototype.hasOwnProperty.call(qualByName, qn)) qualByName[qn] = null;
                else qualByName[qn] = arr;
            }
        }
        if (!qrecs.length) notes.push('The quality input contained no usable records; qualities were ignored.');
    }

    for (var i = 0; i < recs.length; i++) {
        if (reads.length >= P.maxReads) {
            notes.push('Input truncated at ' + P.maxReads + ' reads (maxReads). ' +
                       (recs.length - i) + ' further record(s) were not read.');
            break;
        }
        var raw = capCleanDNA(recs[i].body);
        var name = recs[i].name || ('Seq' + (i + 1));
        // How much of the record survived the nucleotide filter. A protein
        // sequence is not empty — roughly half of its residues happen to be valid
        // IUPAC letters — so it arrives here as a short "read" and would otherwise
        // be reported as merely too short, which sends the user off adjusting a
        // length threshold that was never the problem.
        var letters = (String(recs[i].body).match(/[A-Za-z]/g) || []).length;
        var looksProtein = letters >= 20 && raw.length < letters * 0.7;
        if (!raw.length) {
            dropped.push({ name: name, reason: letters ? 'no nucleotide bases (are these amino acids?)'
                                                       : 'no nucleotide bases' });
            continue;
        }

        var qual = null;
        if (qualByName) {
            var byName = Object.prototype.hasOwnProperty.call(qualByName, name) ? qualByName[name] : undefined;
            var q = byName || qualByOrder[i];      // null = the name was ambiguous
            if (q && q.length) {
                if (q.length !== raw.length) {
                    notes.push('Quality record "' + name + '" has ' + q.length + ' values for ' +
                               raw.length + ' bases; qualities for this read were ignored.');
                } else {
                    qual = q;
                }
            }
        }

        // Clipping: a fixed trim, then (if qualities exist and qualClip > 0) walk
        // in from each end while the quality is below the floor. Sanger reads are
        // unreliable at both ends and an unclipped tail is a common source of
        // spurious mismatches in the consensus.
        var lo = 0, hi = raw.length;
        if (P.endTrim > 0) { lo = Math.min(P.endTrim, raw.length); hi = Math.max(lo, raw.length - P.endTrim); }
        if (qual && P.qualClip > 0) {
            while (lo < hi && qual[lo] < P.qualClip) lo++;
            while (hi > lo && qual[hi - 1] < P.qualClip) hi--;
        }
        // Trailing/leading N runs carry no information and only widen the band.
        while (lo < hi && raw.charAt(lo) === 'n') lo++;
        while (hi > lo && raw.charAt(hi - 1) === 'n') hi--;

        var seq = raw.substring(lo, hi);
        if (seq.length < P.minReadLen) {
            dropped.push({ name: name, reason: looksProtein
                ? ('only ' + raw.length + ' of ' + letters + ' letters are nucleotide codes — ' +
                   'this looks like a protein sequence, not DNA')
                : ('shorter than minReadLen (' + seq.length + ' < ' + P.minReadLen + ' nt)') });
            continue;
        }
        reads.push({
            idx: reads.length,
            id: name,
            desc: recs[i].desc,
            seq: seq,
            len: seq.length,
            qual: qual ? qual.slice(lo, hi) : null,
            origLen: raw.length,
            clipLeft: lo,
            clipRight: raw.length - hi,
            repeatiness: 0        // filled by the seeding stage
        });
    }
    if (dropped.length) {
        notes.push(dropped.length + ' record(s) were not used: ' +
                   dropped.slice(0, 5).map(function (d) { return d.name + ' (' + d.reason + ')'; }).join(', ') +
                   (dropped.length > 5 ? ', …' : '') + '.');
    }
    return { reads: reads, notes: notes, dropped: dropped };
}

/* ───────────────────────── 3. Banded semi-global alignment ───────────────── */

/* The single DP kernel used by both the overlap stage and the consensus stage.
   It is a SEMI-GLOBAL (a.k.a. overlap) alignment with affine gaps: leading and
   trailing gaps are free in BOTH sequences, so it finds a suffix–prefix overlap,
   a containment, or an internal match, whichever scores best.

   Banding: the alignment is restricted to |(j - i) - k0| <= w around the seed
   diagonal k0. Cell (i, j) is stored at band offset o = (j - i) - (k0 - w), so
   the diagonal predecessor is at the same offset in the previous row, the
   vertical predecessor at o+1 in the previous row, and the horizontal one at
   o-1 in the current row. */

// Substitution table over the 26 lower-case letters. Exact IUPAC match scores
// `match`; a disjoint pair scores `mismatch`; a partial overlap (a vs r, or
// anything vs n) scores 0 — neither rewarded nor punished.
var CAP_SUB_CACHE = null;
function capSubTable(P) {
    if (CAP_SUB_CACHE && CAP_SUB_CACHE.m === P.matchScore && CAP_SUB_CACHE.x === P.mismatchScore) {
        return CAP_SUB_CACHE.t;
    }
    var t = new Int16Array(128 * 128);
    var letters = 'acgtrymkwsbdhvn';
    for (var a = 0; a < letters.length; a++) {
        for (var b = 0; b < letters.length; b++) {
            var ca = letters.charAt(a), cb = letters.charAt(b);
            var ba = CAP_BITS[ca], bb = CAP_BITS[cb];
            var v;
            if (ba === bb && (ba === 1 || ba === 2 || ba === 4 || ba === 8)) v = P.matchScore;
            else if (ba & bb) v = 0;
            else v = P.mismatchScore;
            t[ca.charCodeAt(0) * 128 + cb.charCodeAt(0)] = v;
        }
    }
    CAP_SUB_CACHE = { m: P.matchScore, x: P.mismatchScore, t: t };
    return t;
}

// Grow-only scratch buffers, reused across the millions of calls the overlap
// stage makes. Allocating a fresh traceback matrix per pair would dominate the
// runtime and thrash the GC.
var CAP_TB = null, CAP_TB_SIZE = 0;
// Per-column match flags of the traceback, for the local-identity window of
// the junction test (capOverlapPair). Same grow-only discipline as CAP_TB.
var CAP_MF = null, CAP_MF_SIZE = 0;
function capMatchBuf(need) {
    if (CAP_MF_SIZE < need) {
        var buf = new Uint8Array(need);
        CAP_MF = buf;
        CAP_MF_SIZE = need;
    }
    return CAP_MF;
}
// Width, in aligned columns, of the window over which local identity is judged.
// An unrelated stretch this long inside an otherwise good alignment is a
// chimeric junction; the sequencing error of two reads does not produce it
// (see junctionIdentity in CAP_DEFAULTS for the measurements). Narrower than
// 100 so that two chimeras whose junctions lie 30-50 nt apart are still told
// apart; wider than 30 so that a run of errors in a homopolymer-rich stretch
// of two R9 reads is not.
var CAP_JUNCTION_WINDOW = 50;
function capTraceBuf(need) {
    if (CAP_TB_SIZE < need) {
        // Record the new size only AFTER the allocation succeeds. A very wide band
        // on a very long read can ask for hundreds of megabytes in one contiguous
        // block and get a RangeError; if the bookkeeping had already been updated,
        // every later call would pass the size check and be handed the old, too
        // small buffer. Writes past the end of a typed array are silently dropped
        // and reads past it give undefined, so the traceback would decode as a
        // straight diagonal — no exception, no warning, just wrong alignments for
        // the rest of the session.
        var buf = new Uint8Array(need);
        CAP_TB = buf;
        CAP_TB_SIZE = need;
    }
    return CAP_TB;
}

/**
 * Banded semi-global affine alignment of s against t.
 *
 * @param {string} s   first sequence  (rows, index i)
 * @param {string} t   second sequence (cols, index j)
 * @param {number} k0  expected diagonal, j - i, of the true alignment
 * @param {number} w   band half-width
 * @param {object} P   parameters (scoring)
 * @param {boolean} wantOps  build the operation string (needed by the consensus stage)
 * @returns {object|null} alignment, or null when nothing scores above zero
 *   {score, sStart, sEnd, tStart, tEnd, matches, mismatches, gaps, gapRuns,
 *    cols, identity, ops, edgeTouched,
 *    coreSStart, coreSEnd, coreTStart, coreTEnd, coreScore, coreIdentity,
 *    coreMinWindow}
 *   The core* fields describe the maximal-scoring segment of the alignment
 *   (see the traceback), which is what the overhang test is measured on;
 *   coreMinWindow is the lowest identity of any CAP_JUNCTION_WINDOW consecutive
 *   columns inside it (the core's own identity when it is shorter than that).
 *   ops uses the CIGAR convention, read from s's point of view:
 *     'M' aligned pair, 'I' s has a base t lacks, 'D' t has a base s lacks.
 *   It covers [sStart,sEnd) × [tStart,tEnd) only; the free end gaps are not
 *   part of the string (they are implied by the start/end coordinates).
 */
function capAlignBanded(s, t, k0, w, P, wantOps) {
    var n = s.length, m = t.length;
    if (!n || !m) return null;
    if (w > n + m) w = n + m;
    var W = 2 * w + 1;
    var sub = capSubTable(P);
    var go = P.gapOpen + P.gapExtend, ge = P.gapExtend;
    var tb = capTraceBuf((n + 1) * W);

    var prevM = new Int32Array(W), prevX = new Int32Array(W), prevY = new Int32Array(W);
    var curM = new Int32Array(W), curX = new Int32Array(W), curY = new Int32Array(W);
    var swap;

    var bestScore = CAP_NEG, bestI = -1, bestJ = -1, bestSt = 0;

    for (var i = 0; i <= n; i++) {
        var shift = i + k0 - w;                 // j = o + shift
        var jlo = Math.max(0, shift);
        var jhi = Math.min(m, shift + W - 1);
        var base = i * W;
        var o;
        for (o = 0; o < W; o++) { curM[o] = CAP_NEG; curX[o] = CAP_NEG; curY[o] = CAP_NEG; }

        if (jlo <= jhi) {
            var olo = jlo - shift, ohi = jhi - shift;
            var sc = i > 0 ? s.charCodeAt(i - 1) * 128 : 0;
            for (o = olo; o <= ohi; o++) {
                var j = o + shift;
                var dir = 0;
                if (i === 0 && j === 0) {
                    curM[o] = 0;                            // empty alignment
                } else if (i === 0) {
                    curY[o] = 0;                            // free leading gap in s
                } else if (j === 0) {
                    curX[o] = 0;                            // free leading gap in t
                } else {
                    // M — extend the best of the three states on the diagonal
                    var dm = prevM[o], dx = prevX[o], dy = prevY[o];
                    var best = dm, st = 0;
                    if (dx > best) { best = dx; st = 1; }
                    if (dy > best) { best = dy; st = 2; }
                    curM[o] = best < CAP_HALF_NEG ? CAP_NEG : best + sub[sc + t.charCodeAt(j - 1)];
                    dir = st;
                    // X — gap in t, consuming s: comes from (i-1, j) = offset o+1
                    if (o + 1 < W) {
                        var x1 = prevM[o + 1] + go, x2 = prevX[o + 1] + ge;
                        if (x1 >= x2) { curX[o] = x1 < CAP_HALF_NEG ? CAP_NEG : x1; }
                        else { curX[o] = x2 < CAP_HALF_NEG ? CAP_NEG : x2; dir |= 4; }
                    }
                    // Y — gap in s, consuming t: comes from (i, j-1) = offset o-1
                    if (o - 1 >= 0) {
                        var y1 = curM[o - 1] + go, y2 = curY[o - 1] + ge;
                        if (y1 >= y2) { curY[o] = y1 < CAP_HALF_NEG ? CAP_NEG : y1; }
                        else { curY[o] = y2 < CAP_HALF_NEG ? CAP_NEG : y2; dir |= 16; }
                    }
                }
                tb[base + o] = dir;

                // Trailing gaps are free, so any cell on the last row or the last
                // column is a candidate end point.
                if (i === n || j === m) {
                    if (curM[o] > bestScore) { bestScore = curM[o]; bestI = i; bestJ = j; bestSt = 0; }
                    if (curX[o] > bestScore) { bestScore = curX[o]; bestI = i; bestJ = j; bestSt = 1; }
                    if (curY[o] > bestScore) { bestScore = curY[o]; bestI = i; bestJ = j; bestSt = 2; }
                }
            }
        }
        swap = prevM; prevM = curM; curM = swap;
        swap = prevX; prevX = curX; curX = swap;
        swap = prevY; prevY = curY; curY = swap;
    }

    if (bestI < 0 || bestScore <= 0) return null;

    // ── traceback ──
    var ops = wantOps ? [] : null;
    var matches = 0, mismatches = 0, gaps = 0, gapRuns = 0, cols = 0;
    var ci = bestI, cj = bestJ, cst = bestSt, edge = false, lastGap = 0;
    // The maximal-scoring segment of the path is found in the same pass. The
    // alignment is end-space free: it starts on the first row or column and
    // ends on the last, so on each side it runs out to the end of one of the
    // two sequences. Two reads that share only an interspersed repeat are
    // therefore aligned on through the unrelated sequence beyond the repeat,
    // and the overhang that the repeat filter in capOverlapPair() has to see is
    // buried inside the alignment. The segment stops where the reads stop
    // agreeing, so its unaligned ends are those overhangs. Walking back from
    // the end, `run` is the score of the path from the current boundary to the
    // end and `minRun` the lowest such suffix score seen so far, so
    // run − minRun is the best segment that starts here (Kadane's method).
    var run = 0, minRun = 0, minI = bestI, minJ = bestJ, minM = 0, minX = 0, minG = 0, minCols = 0;
    var coreScore = CAP_NEG, coreSStart = bestI, coreSEnd = bestI, coreTStart = bestJ, coreTEnd = bestJ;
    var coreMatches = 0, coreMismatches = 0, coreGaps = 0, coreColLo = 0, coreColHi = 0, step;
    // One flag per column, in traceback order (mf[0] is the alignment's last
    // column), read back below for the local-identity window of the core.
    var mf = capMatchBuf(n + m + 1);
    while (ci > 0 && cj > 0) {
        var co = cj - (ci + k0 - w);
        if (co <= 0 || co >= W - 1) edge = true;
        var d = tb[ci * W + co];
        if (cst === 0) {
            var a = s.charAt(ci - 1), b = t.charAt(cj - 1);
            var ba = CAP_BITS[a], bb = CAP_BITS[b];
            if (ba & bb) { matches++; mf[cols] = 1; } else { mismatches++; mf[cols] = 0; }
            if (ops) ops.push('M');
            cols++; lastGap = 0;
            step = sub[a.charCodeAt(0) * 128 + b.charCodeAt(0)];
            cst = d & 3;
            ci--; cj--;
        } else if (cst === 1) {
            if (ops) ops.push('I');           // s consumed, t gapped: insertion in s
            gaps++; mf[cols] = 0; cols++; if (lastGap !== 1) gapRuns++;
            lastGap = 1;
            cst = (d & 4) ? 1 : 0;
            step = cst === 1 ? ge : go;       // the run's first base pays the opening cost
            ci--;
        } else {
            if (ops) ops.push('D');           // t consumed, s gapped: deletion in s
            gaps++; mf[cols] = 0; cols++; if (lastGap !== 2) gapRuns++;
            lastGap = 2;
            cst = (d & 16) ? 2 : 0;
            step = cst === 2 ? ge : go;
            cj--;
        }
        // (ci, cj) is now the boundary before the column just read; the segment
        // from here to the lowest-scoring suffix boundary is the best one that
        // starts here. Then this boundary becomes a candidate end for the
        // segments that start further back.
        run += step;
        if (run - minRun > coreScore) {
            coreScore = run - minRun;
            coreSStart = ci; coreTStart = cj; coreSEnd = minI; coreTEnd = minJ;
            coreMatches = matches - minM; coreMismatches = mismatches - minX; coreGaps = gaps - minG;
            coreColLo = minCols; coreColHi = cols;
        }
        if (run < minRun) { minRun = run; minI = ci; minJ = cj; minM = matches; minX = mismatches; minG = gaps; minCols = cols; }
    }

    if (ops) ops.reverse();
    var denom = matches + mismatches + gaps;
    var coreDenom = coreMatches + coreMismatches + coreGaps;
    // The worst CAP_JUNCTION_WINDOW consecutive columns of the core. A chimeric
    // junction is a stretch where the two reads are unrelated, so the identity
    // there falls to what random sequence gives (about half), however good the
    // flanks on either side are and however high the overall identity therefore
    // stays. The flanks are what Kadane's segment keeps the stretch inside of:
    // the segment drops a bad stretch only at its ends, never in its middle.
    var minWin = -1;
    if (coreColHi - coreColLo >= CAP_JUNCTION_WINDOW) {
        var wsum = 0;
        for (var wc = coreColLo; wc < coreColHi; wc++) {
            wsum += mf[wc];
            if (wc - coreColLo >= CAP_JUNCTION_WINDOW) wsum -= mf[wc - CAP_JUNCTION_WINDOW];
            if (wc - coreColLo >= CAP_JUNCTION_WINDOW - 1 && (minWin < 0 || wsum < minWin)) minWin = wsum;
        }
    }
    return {
        score: bestScore,
        sStart: ci, sEnd: bestI,
        tStart: cj, tEnd: bestJ,
        matches: matches, mismatches: mismatches, gaps: gaps, gapRuns: gapRuns,
        cols: cols,
        identity: denom ? (100 * matches / denom) : 0,
        ops: ops ? ops.join('') : null,
        edgeTouched: edge,
        // The maximal-scoring segment: where the two sequences actually agree.
        // Its unaligned ends are the overhangs the repeat filter measures.
        coreSStart: coreSStart, coreSEnd: coreSEnd,
        coreTStart: coreTStart, coreTEnd: coreTEnd,
        coreScore: coreScore,
        coreIdentity: coreDenom ? (100 * coreMatches / coreDenom) : 0,
        coreMinWindow: minWin >= 0 ? (100 * minWin / CAP_JUNCTION_WINDOW)
                                   : (coreDenom ? (100 * coreMatches / coreDenom) : 0)
    };
}

/**
 * Align with automatic band widening: if the optimal path pressed against the
 * band edge, the true alignment may lie outside it, so double the band and retry
 * (bounded by maxBandWidth). Without this, a read with a real indel near the seed
 * diagonal silently loses part of its overlap.
 *
 * `gate` decides whether widening is worth it. It matters a great deal: most
 * candidate pairs at the overlap stage are chance k-mer collisions between
 * unrelated reads, and a garbage alignment ALWAYS hugs the band edge — so
 * without a gate every one of them is re-aligned at 2x, 4x, 8x the band, and the
 * pairs that will be thrown away come to dominate the entire run. Widening only
 * pays when what is already there looks like a real overlap that got clipped.
 */
function capAlignAuto(s, t, k0, P, wantOps, gate) {
    var w = P.bandWidth;
    var res = capAlignBanded(s, t, k0, w, P, wantOps);
    var guard = 0;
    while (res && res.edgeTouched && w < P.maxBandWidth && guard++ < 4 &&
           (!gate || gate(res))) {
        w = Math.min(P.maxBandWidth, w * 2);
        res = capAlignBanded(s, t, k0, w, P, wantOps);
    }
    if (!res && !gate && w < P.maxBandWidth) {
        // Nothing at all in the narrow band: one wide retry before giving up.
        // Only for gate-less callers (the consensus stage), where the read is
        // known to belong here and losing it would corrupt the contig.
        res = capAlignBanded(s, t, k0, P.maxBandWidth, P, wantOps);
    }
    return res;
}

/* ───────────────────────── 3b. Fold-back reads ───────────────────────────── */

/* A fold-back read carries a stretch of its own template twice, the second time
   reverse-complemented: the strand turns at a hairpin and is read back, or the
   complementary strand follows the template strand through the pore. The
   pile-based test in capChimeras() cannot see it. The read's body overlaps
   every sibling, and the tail dents the identity of those overlaps too little
   to fail them, so the read is accepted whole and its tail is laid beyond the
   end of the product; two such tails of different length disagree with each
   other there and the disagreement is written out as ambiguity codes. A 2000 nt
   product came out as a 3068 nt contig with 217 IUPAC columns from two such
   reads among 200.

   The signature is inside the read. A word of the copy is the reverse
   complement of a word of the original, and for every such pair of words the
   two positions add up to the same number (the net indels between the copies
   move it by a few). So each read's words are looked up against their own
   reverse complements, the pairs are bucketed by that sum, and the strongest
   bucket is confirmed with the banded aligner between the read and its reverse
   complement — the same test an overlap between two reads has to pass. The
   confirmed alignment gives the two copies as intervals on the read, mirror
   images of each other. Then:
     • the copies lie next to each other, or overlap (a hairpin, or two
       molecules back to back with at most an adapter's length between them):
       the read is cut between them, or at the fold, and the longer side kept;
     • they lie apart (the tail copies the start of the read rather than its
       end): either copy may be the one that continues the sequence between
       them, so only that middle part is kept;
     • what would be kept is shorter than minReadLen: the read is dropped.
   The copy must reach a read end within maxOverhang bases. A template that
   itself holds an inverted repeat shows the same words in every read, but
   inside the read, with sequence going on after the second copy; only a read
   that ends within the repeat looks like a fold-back, and such a read is cut
   back to the start of the second copy, which the other reads still cover.

   Runs before seeding, so the index, the overlaps and the layout only ever
   see the trimmed read. The bases cut off count as clipped in the placement
   table; a dropped read is listed with its name. Returns
   [{idx, id, len, copies:[e0,e1,l0,l1], kept:[from,to], identity, action}],
   coordinates 0-based half-open on the read as given, idx −1 when dropped. */
var CAP_FOLD_BIN = 32;      // bucket width for the position sums, as DIAGBIN in seeding
var CAP_FOLD_LOOP = 100;    // copies this close are one fold (an adapter can sit between them)
var CAP_FOLD_JUNK = 60;     // a stretch scoring this far below its surroundings parts two copies

/* The aligned copy, from the alignment of a read against its reverse complement.
 *
 * When the tail copies the start of the read, the copy and the original both lie
 * on one diagonal of that alignment, with the rest of the read between them —
 * sequence aligned against its own reverse complement, which matches nothing.
 * The maximal-scoring segment takes both copies and the stretch between them
 * together whenever the stretch costs less than a copy scores, so it cannot
 * be read as one copy. The stretch is found as the lowest-scoring part of the
 * segment (Kadane's method for the minimum), and when it scores CAP_FOLD_JUNK
 * below its surroundings the segment is cut there and the part before it is
 * the copy. A hairpin's arms align in one piece, with at most the loop as a
 * dip, and the whole segment is the copy.
 *
 * Returns {s0, s1, t0, t1, identity}: the copy on the read (s) and on the
 * reverse complement (t), with the identity over its columns. */
function capSelfCopy(s, t, aln, P) {
    var sub = capSubTable(P), go = P.gapOpen + P.gapExtend, ge = P.gapExtend;
    var ops = aln.ops || '', i = aln.sStart, j = aln.tStart, last = '';
    var inCore = false, run = 0, maxRun = 0, minSeg = 0, m = 0, c = 0;
    var maxI = aln.coreSStart, maxJ = aln.coreTStart, maxM = 0, maxC = 0;
    var cutI = -1, cutJ = -1, cutM = 0, cutC = 0;
    for (var o = 0; o < ops.length; o++) {
        if (!inCore && i === aln.coreSStart && j === aln.coreTStart) inCore = true;
        if (i === aln.coreSEnd && j === aln.coreTEnd) break;
        var op = ops.charAt(o), step, match = 0;
        if (op === 'M') {
            var a = s.charAt(i), b = t.charAt(j);
            match = (CAP_BITS[a] & CAP_BITS[b]) ? 1 : 0;
            step = sub[a.charCodeAt(0) * 128 + b.charCodeAt(0)];
            i++; j++;
        } else {
            step = last === op ? ge : go;
            if (op === 'I') i++; else j++;
        }
        last = op;
        if (!inCore) continue;
        run += step; m += match; c++;
        if (run > maxRun) { maxRun = run; maxI = i; maxJ = j; maxM = m; maxC = c; }
        if (run - maxRun < minSeg) { minSeg = run - maxRun; cutI = maxI; cutJ = maxJ; cutM = maxM; cutC = maxC; }
    }
    if (minSeg <= -CAP_FOLD_JUNK && cutI >= 0) {
        return { s0: aln.coreSStart, s1: cutI, t0: aln.coreTStart, t1: cutJ,
                 identity: cutC ? 100 * cutM / cutC : 0 };
    }
    return { s0: aln.coreSStart, s1: aln.coreSEnd, t0: aln.coreTStart, t1: aln.coreTEnd,
             identity: aln.coreIdentity };
}

function capFoldbacks(reads, P) {
    var K = P.kmer, W = P.maxOverhang, found = [];
    if (!P.chimeraCheck || !(P.foldbackMin > 0)) return found;
    var mask = K === 15 ? 0x3FFFFFFF : ((1 << (2 * K)) - 1), rcShift = 2 * (K - 1);
    var codes = null, rcodes = null, ok = null, keyT = null, posT = null, nextT = null, cap = 0, tcap = 0;
    var gate = function (r) {
        return r.coreIdentity >= P.minIdentity - 5 &&
               Math.max(r.coreSEnd - r.coreSStart, r.coreTEnd - r.coreTStart) >= P.foldbackMin * 0.5;
    };
    for (var ri = 0; ri < reads.length; ri++) {
        var rd = reads[ri], L = rd.len, seq = rd.seq;
        // Nothing could be cut off and still leave a read.
        if (L < P.foldbackMin + P.minReadLen) continue;
        var nW = L - K + 1, p, c;
        if (cap < nW) { cap = nW; codes = new Int32Array(cap); rcodes = new Int32Array(cap); ok = new Uint8Array(cap); nextT = new Int32Array(cap); }
        // Forward and reverse-complement code of every word; a word spanning an
        // N, or a palindrome (its own reverse complement), is unusable.
        var code = 0, rcode = 0, valid = 0;
        for (p = 0; p < L; p++) {
            c = CAP_CODE[seq.charCodeAt(p)];
            if (c < 0) { valid = 0; code = 0; rcode = 0; } else {
                code = ((code << 2) | c) & mask;
                rcode = ((rcode >>> 2) | ((3 - c) << rcShift)) & mask;
                valid++;
            }
            if (p >= K - 1) {
                var at = p - K + 1;
                codes[at] = code; rcodes[at] = rcode;
                ok[at] = (valid >= K && code !== rcode) ? 1 : 0;
            }
        }
        // Open-addressing hash of the forward codes: keyT holds the code, posT
        // the first position with it, nextT chains the others.
        var size = 1; while (size < 2 * nW) size <<= 1;
        if (tcap < size) { tcap = size; keyT = new Int32Array(tcap); posT = new Int32Array(tcap); }
        for (p = 0; p < size; p++) keyT[p] = -1;
        var hmask = size - 1;
        for (p = 0; p < nW; p++) {
            if (!ok[p]) continue;
            var h = Math.imul(codes[p], 0x9E3779B1) >>> 0 & hmask;
            while (keyT[h] !== -1 && keyT[h] !== codes[p]) h = (h + 1) & hmask;
            if (keyT[h] === -1) { keyT[h] = codes[p]; posT[h] = p; nextT[p] = -1; }
            else { nextT[p] = posT[h]; posT[h] = p; }
        }
        // Every (p, p') with word(p) = rc(word(p')), p < p', bucketed by p + p'.
        var bins = new Map();
        for (var q = 0; q < nW; q++) {
            if (!ok[q]) continue;
            var hq = Math.imul(rcodes[q], 0x9E3779B1) >>> 0 & hmask;
            while (keyT[hq] !== -1 && keyT[hq] !== rcodes[q]) hq = (hq + 1) & hmask;
            if (keyT[hq] === -1) continue;
            for (p = posT[hq]; p !== -1; p = nextT[p]) {
                if (p >= q) continue;
                var bin = Math.floor((p + q) / CAP_FOLD_BIN);
                var e = bins.get(bin);
                if (e) { e.cnt++; e.sum += p + q; } else bins.set(bin, { cnt: 1, sum: p + q });
            }
        }
        if (!bins.size) continue;
        var best = null, bestCnt = 0, bestBin = 0;
        bins.forEach(function (e, b) {
            var cnt = e.cnt, lo = bins.get(b - 1), hi = bins.get(b + 1);
            if (lo) cnt += lo.cnt;
            if (hi) cnt += hi.cnt;
            if (cnt > bestCnt || (cnt === bestCnt && b < bestBin)) { bestCnt = cnt; best = e; bestBin = b; }
        });
        if (bestCnt < P.minSeeds) continue;
        // Diagonal j − i of the read (i) against its reverse complement (j): the
        // word at p' on the read is the word at L − K − p' on the reverse
        // complement, so the pair (p, p') lies on j − i = L − K − (p + p').
        var d = L - K - Math.round(best.sum / best.cnt);
        var rcs = capRC(seq);
        var aln = capAlignAuto(seq, rcs, d, P, true, gate);
        if (!aln) continue;
        var copy = capSelfCopy(seq, rcs, aln, P);
        var cLen = Math.max(copy.s1 - copy.s0, copy.t1 - copy.t0);
        if (cLen < P.foldbackMin || copy.identity < P.minIdentity) continue;
        // The two copies on the read: the aligned stretch of the read itself, and
        // the aligned stretch of the reverse complement mapped back.
        var i0 = copy.s0, i1 = copy.s1, f0 = L - copy.t1, f1 = L - copy.t0;
        var e0, e1, l0, l1;
        if (i0 <= f0) { e0 = i0; e1 = i1; l0 = f0; l1 = f1; } else { e0 = f0; e1 = f1; l0 = i0; l1 = i1; }
        if (l1 < L - W && e0 > W) continue;         // an inverted repeat inside the read
        var from, to, gap = l0 - e1;
        if (gap >= Math.max(CAP_FOLD_LOOP, P.minReadLen)) {
            // Copies apart: keep what lies between. Either copy may be the one the
            // strand re-read, so neither is kept (a read X·Y·rc(X) is molecule X·Y
            // read on into its own complement, or molecule Y·rc(X) preceded by it,
            // and nothing on the read tells which). The alignment places a copy's
            // edge only to within the last bases that failed to match, so the
            // kept stretch is drawn in by one word at each end: measured on a
            // fold-back read the copy began 4 nt earlier than the alignment said,
            // and those 4 nt reached past the end of the product in the consensus.
            from = e1 + K; to = l0 - K;
        } else if (gap >= 0) {
            if (e1 >= L - l0) { from = 0; to = e1; } else { from = l0; to = L; }
        } else {
            // A hairpin with no loop between the copies aligns to its own reverse
            // complement over its whole length, so both copies came out as that one
            // stretch; the fold is its midpoint, and the copies are reported either
            // side of it.
            var cut = Math.round((e0 + e1 + l0 + l1) / 4);
            e1 = cut; l0 = cut;
            if (cut >= L - cut) { from = 0; to = cut; } else { from = cut; to = L; }
        }
        // Reported on the read as given: capParseReads may already have clipped
        // its ends (endTrim, qualClip, terminal N), so the offsets found here are
        // on the clipped sequence and move by that clip.
        var base = rd.clipLeft;
        var rec = { idx: rd.idx, id: rd.id, len: rd.origLen, copies: [e0 + base, e1 + base, l0 + base, l1 + base],
                    kept: [from + base, to + base], identity: copy.identity, action: 'trimmed', read: rd };
        if (to - from < P.minReadLen) {
            rec.action = 'dropped'; rec.kept = null;
        } else {
            rd.seq = seq.substring(from, to);
            rd.len = to - from;
            if (rd.qual) rd.qual = rd.qual.slice(from, to);
            rd.clipLeft += from;
            rd.clipRight += L - to;
        }
        found.push(rec);
    }
    // Dropped reads leave the array; indices are positions, so renumber.
    var drop = new Set();
    found.forEach(function (f) { if (f.action === 'dropped') drop.add(f.read); });
    if (drop.size) {
        var kept = reads.filter(function (r) { return !drop.has(r); });
        reads.length = 0;
        for (var k = 0; k < kept.length; k++) { kept[k].idx = k; reads.push(kept[k]); }
    }
    found.forEach(function (f) { f.idx = f.action === 'dropped' ? -1 : f.read.idx; delete f.read; });
    return found;
}

/* ───────────────────────── 4. k-mer seeding & candidate pairs ────────────── */

/* Reads are indexed by CANONICAL k-mer (the lexicographically smaller of the
   2-bit packing of the word and of its reverse complement) together with a flag
   saying which of the two was smaller. Two occurrences with equal flags are on
   the same strand; with different flags, one read must be reverse-complemented.
   That gives both-strand seeding at the cost of a single forward pass, instead
   of indexing every read twice.

   Palindromic words (code === rcCode) are skipped: their flag is meaningless, so
   they would produce arbitrary orientation calls. With k = 12 they are 4^6 of
   4^12 words, i.e. 0.02% — no loss of sensitivity. */

var CAP_CODE = new Int8Array(128);      // a→0, c→1, g→2, t→3, everything else −1
(function () {
    var i;
    for (i = 0; i < 128; i++) CAP_CODE[i] = -1;
    CAP_CODE[97] = 0; CAP_CODE[99] = 1; CAP_CODE[103] = 2; CAP_CODE[116] = 3;
})();

// LSD radix sort of `order` (indices into `keys`) by the 30-bit key, three
// 10-bit passes. Deterministic and linear — Array#sort with a comparator on a
// million entries is both slower and (across engines) not guaranteed stable.
function capRadixSort(order, keys, n) {
    var tmp = new Uint32Array(n);
    var count = new Uint32Array(1024);
    for (var shift = 0; shift < 30; shift += 10) {
        var c;
        for (c = 0; c < 1024; c++) count[c] = 0;
        for (var i = 0; i < n; i++) count[(keys[order[i]] >>> shift) & 1023]++;
        var sum = 0;
        for (c = 0; c < 1024; c++) { var t = count[c]; count[c] = sum; sum += t; }
        for (var j = 0; j < n; j++) tmp[count[(keys[order[j]] >>> shift) & 1023]++] = order[j];
        order.set(tmp);
    }
    return order;
}

/**
 * Build the k-mer index and enumerate candidate read pairs.
 *
 * @returns {{pairs:Array, stats:object}}  pairs = [{a, b, flip, diag, seeds}]
 *   `flip` is 1 when read b must be reverse-complemented to agree with read a.
 *   `diag` is the estimated j − i diagonal of the alignment of read a against
 *   the (possibly reverse-complemented) read b.
 */
var CAP_DEEP_FACTOR = 4;    // word counts up to this × the median are coverage, not repeat
var CAP_DEEP_SAMPLE = 8;    // one deep-coverage word in this many is used for seeding
// Reads that follow each read as its partners on a deep-coverage word, so each
// read meets about twice this many. Every partner costs a full alignment: for
// 1000 reads of one amplicon, 50 placed 991 reads in 191 s, 20 placed 975 in
// 110 s, 10 placed 963 in 87 s, and all three found every chimera.
var CAP_DEEP_PARTNERS = 10;

async function capCandidatePairs(reads, P, ctx) {
    var K = P.kmer, nReads = reads.length;
    var i, r, seq, L;

    // ── pass 1: count windows so the typed arrays can be sized exactly ──
    var total = 0;
    for (i = 0; i < nReads; i++) { L = reads[i].len; if (L >= K) total += L - K + 1; }
    if (!total) return { pairs: [], stats: { words: 0, seeds: 0, skippedWords: 0, truncated: false } };

    var keys = new Int32Array(total);    // canonical 2-bit code
    // Read index and position are kept in SEPARATE arrays rather than packed into
    // one 32-bit word. Packing as (readIdx << 21 | pos) silently breaks: `<<` is a
    // 32-bit operation, so read 2048 shifts to 0 and every read above that is
    // mis-attributed — and maxReads defaults to 5000, so the corruption happens
    // well inside the supported range, with no error. It also capped read length
    // at 2 Mb. Two arrays cost 4 more bytes per k-mer and have no such limits.
    var locR = new Int32Array(total);    // read index
    var locP = new Int32Array(total);    // position within that read
    var flags = new Uint8Array(total);   // 1 when the reverse-complement was canonical
    var mask = K === 15 ? 0x3FFFFFFF : ((1 << (2 * K)) - 1);
    var rcShift = 2 * (K - 1);
    var w = 0;

    for (i = 0; i < nReads; i++) {
        seq = reads[i].seq; L = reads[i].len;
        if (L < K) continue;
        var code = 0, rcode = 0, valid = 0;
        for (var p = 0; p < L; p++) {
            var c = CAP_CODE[seq.charCodeAt(p)];
            if (c < 0) { valid = 0; code = 0; rcode = 0; continue; }   // an N resets the window
            code = ((code << 2) | c) & mask;
            rcode = (rcode >>> 2) | ((3 - c) << rcShift);
            rcode &= mask;
            valid++;
            if (valid >= K) {
                var start = p - K + 1;
                if (code === rcode) continue;                          // palindrome: unusable flag
                var canon, fl;
                if (code < rcode) { canon = code; fl = 0; } else { canon = rcode; fl = 1; }
                keys[w] = canon;
                locR[w] = i;
                locP[w] = start;
                flags[w] = fl;
                w++;
            }
        }
    }

    // ── sort by canonical code so equal words become adjacent runs ──
    var order = new Uint32Array(w);
    for (i = 0; i < w; i++) order[i] = i;
    capRadixSort(order, keys, w);

    // ── pass 2: walk the runs, emit seeds, bucket them per (pair, diagonal) ──
    // Bucketing by diagonal in 32 nt bins tolerates the indels that shift the
    // true diagonal, while still separating genuinely different placements
    // (which is what keeps two copies of a repeat apart).
    var DIAGBIN = 32;
    var buckets = new Map();       // "a|b|flip|bin" -> {a,b,flip,sum,cnt}
    var seedCount = 0, skippedWords = 0, deepWords = 0, truncated = false;

    // A word seen more than maxSeedOcc times is either a repeat word or simply a
    // word of a deep amplicon: at 200 reads every error-free word of the
    // amplicon is seen about 170 times. Skipping both, as was done, left the
    // reads of a deep amplicon with only the words their errors happen to share
    // — 22 of 200 reads, 92 of 1000 across five amplicons, ended as singletons,
    // stray two-read contigs appeared, and chimera halves found too few partners
    // to be recognised (10 of 30 at 200 reads per amplicon). The two are told
    // apart by how the count compares with the typical count, measured here in
    // a first pass: up to CAP_DEEP_FACTOR times the median is coverage, beyond it
    // a repeat. A deep word pairs each read only with the next `succ` reads that
    // contain it (in read order), which gives every read a fixed number of
    // partners instead of all of them; one deep word in CAP_DEEP_SAMPLE is used,
    // since a pair needs only minSeeds shared words and has hundreds.
    var runSizes = [];
    for (var r0 = 0; r0 < w; ) {
        var r1 = r0 + 1, k0 = keys[order[r0]];
        while (r1 < w && keys[order[r1]] === k0) r1++;
        runSizes.push(r1 - r0);
        r0 = r1;
    }
    var deepCut = P.maxSeedOcc;
    if (runSizes.length) {
        // Median over word OCCURRENCES, i.e. the count a typical base is part of.
        var sizesSorted = runSizes.slice().sort(function (x, y) { return x - y; });
        var half = w / 2, accW = 0, medG = 1;
        for (var si0 = 0; si0 < sizesSorted.length; si0++) {
            accW += sizesSorted[si0];
            if (accW >= half) { medG = sizesSorted[si0]; break; }
        }
        deepCut = Math.max(P.maxSeedOcc, CAP_DEEP_FACTOR * medG);
    }
    var succ = Math.max(5, Math.min(CAP_DEEP_PARTNERS,
                                    Math.floor(P.maxCandidates / (3 * Math.max(1, nReads)))));
    var occ = [];

    // One seed between occurrences ox and oy of the same word. Returns false
    // once the candidate table is full.
    var addSeed = function (ox, oy) {
        var ax = locR[ox], ay = locR[oy];
        if (ax === ay) return true;                         // self-hit
        var flip = flags[ox] === flags[oy] ? 0 : 1;
        if (flip && !P.allowRC) return true;
        // Order the pair as (lo, hi) so a pair is only ever seen once.
        var a, b, pa, pb;
        if (ax < ay) { a = ax; b = ay; pa = locP[ox]; pb = locP[oy]; }
        else { a = ay; b = ax; pa = locP[oy]; pb = locP[ox]; }
        // Position of the word on read b as it will be ALIGNED: if b is
        // flipped, the word at pb on b's forward strand starts at
        // (Lb − K − pb) on the reverse complement.
        var pbAln = flip ? (reads[b].len - K - pb) : pb;
        var diag = pbAln - pa;
        var bin = Math.floor(diag / DIAGBIN);
        var bk = a + '|' + b + '|' + flip + '|' + bin;
        var e = buckets.get(bk);
        if (e) { e.sum += diag; e.cnt++; }
        else {
            if (buckets.size >= P.maxCandidates) { truncated = true; return false; }
            buckets.set(bk, { a: a, b: b, flip: flip, sum: diag, cnt: 1 });
        }
        seedCount++;
        return true;
    };

    for (i = 0; i < nReads; i++) { reads[i].mult = new Uint16Array(reads[i].len); reads[i].pile = [0, 0]; }

    var run = 0, runsSinceTick = 0;
    while (run < w) {
        var end = run + 1;
        var key = keys[order[run]];
        while (end < w && keys[order[end]] === key) end++;
        var g = end - run;
        if (g >= 2) {
            // Record this word's copy number at every position it occurs at. The
            // per-BASE copy number is what tells a repeat from unique sequence:
            // in unique sequence a word is seen about once per read covering it,
            // i.e. about the depth of coverage, whereas inside an n-copy repeat it
            // is seen n times that. capRepeatThreshold() turns the distribution
            // into a cut-off and capOverlapPair() then refuses overlaps that lie
            // wholly inside repeat sequence — the joins that collapse two loci
            // into one wrong contig.
            for (var mq = run; mq < end; mq++) {
                var mo = order[mq];
                var mr = reads[locR[mo]];
                if (mr.mult) mr.mult[locP[mo]] = g > 65535 ? 65535 : g;
            }
            capTerminalPiles(order, run, end, locR, locP, flags, reads, K);
            if (g > deepCut) {
                // Repeat word. Skipped for seeding — pairing every occurrence with
                // every other is quadratic and produces only repeat-induced joins,
                // which the overhang filter would reject anyway. It IS counted, so
                // reads sitting in repeats can be flagged in the report.
                skippedWords++;
                for (var q = run; q < end; q++) reads[locR[order[q]]].repeatiness++;
            } else if (g > P.maxSeedOcc) {
                // Deep-coverage word: each occurrence is paired with the next
                // `succ` ones in read order.
                deepWords++;
                if (((Math.imul(key, 0x9E3779B1) >>> 24) % CAP_DEEP_SAMPLE) === 0) {
                    occ.length = 0;
                    for (var oq = run; oq < end; oq++) occ.push(order[oq]);
                    occ.sort(function (x, y) { return locR[x] - locR[y] || locP[x] - locP[y]; });
                    for (var dx = 0; dx < occ.length && !truncated; dx++) {
                        var dEnd = Math.min(occ.length, dx + 1 + succ);
                        for (var dy = dx + 1; dy < dEnd; dy++) if (!addSeed(occ[dx], occ[dy])) break;
                    }
                }
            } else {
                for (var x = run; x < end && !truncated; x++) {
                    for (var y = x + 1; y < end; y++) if (!addSeed(order[x], order[y])) break;
                }
            }
        }
        run = end;
        // Awaiting on every k-mer run would cost more than the work itself, so
        // only offer the event loop a turn every few thousand runs; ctx.tick()
        // then decides, on elapsed time, whether to actually yield.
        if (ctx && (++runsSinceTick & 4095) === 0) await ctx.tick('seeding', run, w);
    }

    // ── collapse buckets to one candidate per (pair, orientation) ──
    // Neighbouring bins belong to the same alignment; keep the strongest bin per
    // (a, b, flip) and add the counts of the two adjacent bins to it.
    var byPair = new Map();
    buckets.forEach(function (e, bk) {
        var parts = bk.split('|');
        var pk = parts[0] + '|' + parts[1] + '|' + parts[2];
        var bin = parseInt(parts[3], 10);
        var cur = byPair.get(pk);
        if (!cur) byPair.set(pk, { a: e.a, b: e.b, flip: e.flip, best: e, bins: new Map([[bin, e]]) });
        else {
            cur.bins.set(bin, e);
            if (e.cnt > cur.best.cnt || (e.cnt === cur.best.cnt && e.sum / e.cnt < cur.best.sum / cur.best.cnt)) cur.best = e;
        }
    });

    var pairs = [];
    byPair.forEach(function (v, pk) {
        var parts = pk.split('|');
        var bin = Math.floor((v.best.sum / v.best.cnt) / DIAGBIN);
        var seeds = v.best.cnt;
        var lo = v.bins.get(bin - 1), hi = v.bins.get(bin + 1);
        if (lo) seeds += lo.cnt;
        if (hi) seeds += hi.cnt;
        if (seeds < P.minSeeds) return;
        pairs.push({ a: v.best.a, b: v.best.b, flip: v.best.flip,
                     diag: Math.round(v.best.sum / v.best.cnt), seeds: seeds });
    });

    // Deterministic order: strongest seed support first, then by read indices.
    pairs.sort(function (p, q) {
        if (q.seeds !== p.seeds) return q.seeds - p.seeds;
        if (p.a !== q.a) return p.a - q.a;
        if (p.b !== q.b) return p.b - q.b;
        return p.flip - q.flip;
    });

    capSmoothMultiplicity(reads, P);
    return { pairs: pairs, stats: { words: w, seeds: seedCount, skippedWords: skippedWords,
                                    deepWords: deepWords, deepCut: deepCut,
                                    truncated: truncated, repeatThreshold: capRepeatThreshold(reads, P) } };
}

/**
 * Count, for each end of each read, how many OTHER reads end at the same place.
 *
 * PCR amplicons are read end to end, so every read of an amplicon starts and
 * stops where its siblings do; shotgun reads end anywhere. The amplicon
 * exemption in capCoverageStep() relies on this. One word run (all occurrences
 * of one canonical k-mer) is examined at a time: an occurrence near a read end
 * is described by which end of the TEMPLATE that read end is (the read end
 * flipped by the strand the word was read on) and how far the word lies from it.
 * Occurrences with the same template end and distances within CAP_PILE_JITTER
 * belong to reads that end together. Each read end keeps the largest such count
 * over its words, which approximates the number of distinct reads sharing it.
 */
var CAP_PILE_WINDOW = 24;   // words starting this close to a read end are looked at
var CAP_PILE_JITTER = 10;   // read ends this close together count as the same end
var CAP_PILE_MIN = 2;       // other reads needed for a read end to be a shared one

function capTerminalPiles(order, run, end, locR, locP, flags, reads, K) {
    var hits = null;
    for (var q = run; q < end; q++) {
        var o = order[q], rd = reads[locR[o]], p = locP[o];
        var dS = p, dE = rd.len - K - p;
        if (dS >= CAP_PILE_WINDOW && dE >= CAP_PILE_WINDOW) continue;
        if (!hits) hits = [];
        if (dS < CAP_PILE_WINDOW) hits.push({ read: rd, end: 0, side: flags[o], d: dS });
        if (dE < CAP_PILE_WINDOW) hits.push({ read: rd, end: 1, side: flags[o] ^ 1, d: dE });
    }
    if (!hits || hits.length < 2) return;
    hits.sort(function (x, y) { return x.side - y.side || x.d - y.d; });
    // Two pointers over each side's distances: [lo, hi) holds the hits within
    // the jitter of hit i.
    var lo = 0, hi = 0;
    for (var i = 0; i < hits.length; i++) {
        var h = hits[i];
        while (hits[lo].side !== h.side || hits[lo].d < h.d - CAP_PILE_JITTER) lo++;
        if (hi < i) hi = i;
        while (hi < hits.length && hits[hi].side === h.side && hits[hi].d <= h.d + CAP_PILE_JITTER) hi++;
        var others = hi - lo - 1;
        if (others > h.read.pile[h.end]) h.read.pile[h.end] = others;
    }
}

/**
 * Smooth the per-base word copy numbers over a ±k window.
 *
 * A single sequencing error destroys every k-mer that spans it, so inside a
 * repeat about one base in six carries a word with copy number 1 rather than the
 * repeat's true copy number (at 1.5% error and k = 12, P(clean word) = 0.83).
 * Raw copy numbers therefore mark only a fraction of a repeat as repetitive —
 * measured at 3.7% of positions where 6.7% were genuinely repeated — and a
 * repeat-only overlap looks about half unique, which is not enough to act on.
 *
 * Averaging over a window as wide as the damage an error does restores the
 * signal (an isolated blanked word is outvoted by its neighbours) while leaving
 * genuine repeat boundaries where they are. Positions with no word at all — the
 * last k-1 bases of a read, or windows containing an N — are excluded from the
 * average rather than counted as zero, so a read end is not made to look unique.
 */
function capSmoothMultiplicity(reads, P) {
    var w = P.kmer;
    for (var i = 0; i < reads.length; i++) {
        var m = reads[i].mult;
        if (!m || !m.length) continue;
        var n = m.length, out = new Uint16Array(n);
        var sum = 0, cnt = 0, lo = 0, hi = -1, p, a, b;
        for (p = 0; p < n; p++) {
            a = p - w; if (a < 0) a = 0;
            b = p + w; if (b > n - 1) b = n - 1;
            while (hi < b) { hi++; if (m[hi]) { sum += m[hi]; cnt++; } }
            while (lo < a) { if (m[lo]) { sum -= m[lo]; cnt--; } lo++; }
            out[p] = cnt ? Math.round(sum / cnt) : 0;
        }
        reads[i].mult = out;
    }
}

/**
 * Turn the per-base word copy numbers into a "this is repeat sequence" cut-off.
 *
 * The median copy number over all bases is, by construction, the copy number of
 * ordinary unique sequence — it is essentially the depth of coverage, whatever
 * that happens to be, so nothing has to be assumed about the data. Sequence whose
 * words are seen `repeatFactor` times more often than that is present in more
 * copies than once.
 *
 * A floor of 3 keeps the rule quiet on very shallow input, where the median is
 * 1–2 and ordinary sampling noise would otherwise look like a repeat.
 */
function capRepeatThreshold(reads, P) {
    var hist = new Uint32Array(1024), total = 0, i, j, m;
    for (i = 0; i < reads.length; i++) {
        var mult = reads[i].mult;
        if (!mult) continue;
        for (j = 0; j < mult.length; j++) {
            m = mult[j];
            if (!m) continue;                       // position never indexed (an N, or the tail < k)
            hist[m > 1023 ? 1023 : m]++;
            total++;
        }
    }
    if (!total) return Infinity;
    var half = total >> 1, acc = 0, median = 1;
    for (i = 1; i < 1024; i++) { acc += hist[i]; if (acc >= half) { median = i; break; } }
    return Math.max(3, Math.ceil(P.repeatFactor * median));
}

// Fraction of [from, to) on a read whose words exceed the repeat cut-off.
function capRepeatFraction(read, from, to, threshold) {
    if (!read.mult || to <= from || !isFinite(threshold)) return 0;
    var n = 0, tot = 0;
    var hi = Math.min(to, read.mult.length);
    for (var i = Math.max(0, from); i < hi; i++) {
        tot++;
        if (read.mult[i] >= threshold) n++;
    }
    return tot ? n / tot : 0;
}

// Longest run of single-copy positions in [from, to) of a read. Positions with no
// copy number at all (no word covers them) do not count as single-copy.
function capSingleCopyRun(read, from, to, threshold) {
    if (!read.mult) return 0;
    var best = 0, cur = 0, hi = Math.min(to, read.mult.length);
    for (var i = Math.max(0, from); i < hi; i++) {
        var m = read.mult[i];
        if (m && m < threshold) { if (++cur > best) best = cur; } else cur = 0;
    }
    return best;
}

// Does the MATCHING part of a mostly-repeat overlap take in single-copy sequence
// on both reads? Then the two reads share sequence that exists once. The core
// segment is used, not the whole end-free alignment: two reads from different
// copies of a repeat are aligned on past the repeat's end through unrelated
// flank, and that flank is single-copy but does not match. A run of at least
// max(minOverlap, 2k) bases is required, well beyond the few positions at a
// repeat's edge that the copy-number smoothing blurs.
function capSharesSingleCopy(A, B, flip, aln, P, threshold) {
    var need = Math.max(P.minOverlap, 2 * P.kmer);
    if (capSingleCopyRun(A, aln.coreSStart, aln.coreSEnd, threshold) < need) return false;
    var from = flip ? B.len - aln.coreTEnd : aln.coreTStart;
    var to = flip ? B.len - aln.coreTStart : aln.coreTEnd;
    return capSingleCopyRun(B, from, to, threshold) >= need;
}

// Is a mostly-high-copy overlap a join between PCR amplicon reads rather than a
// repeat? Two shapes qualify. A junction between neighbouring amplicons is a
// coverage step (capCoverageStep). Two reads of ONE amplicon, one of them broken
// off inside an amplicon overlap, end together at the amplicon's end and share
// its single-copy body; refusing them while the broken read still joins the
// neighbouring amplicon leaves it in a cluster of its own, and the branch guard
// then reads the whole junction as a repeat branch and cuts the contig there.
// Both shapes need read ends shared with other reads, which shotgun reads almost
// never have, so shotgun data is judged exactly as before.
function capAmpliconJoin(A, B, flip, aln, P, threshold) {
    if (capCoverageStep(A, B, flip, aln, P, threshold)) return true;
    var W = Math.max(P.kmer, P.maxOverhang);
    var bStart = flip ? 1 : 0, bEnd = flip ? 0 : 1;
    var sharedLeft = aln.sStart <= W && aln.tStart <= W &&
                     A.pile && B.pile && A.pile[0] >= CAP_PILE_MIN && B.pile[bStart] >= CAP_PILE_MIN;
    var sharedRight = A.len - aln.sEnd <= W && B.len - aln.tEnd <= W &&
                      A.pile && B.pile && A.pile[1] >= CAP_PILE_MIN && B.pile[bEnd] >= CAP_PILE_MIN;
    return (sharedLeft || sharedRight) && capSharesSingleCopy(A, B, flip, aln, P, threshold);
}

// Is an all-high-copy overlap a coverage step between PCR amplicons rather than
// a repeat? Each of the four read ends is looked at. Where a read runs on past
// the overlap by more than W bases, the W bases next to the overlap must be
// single-copy: in a repeat the high-copy stretch continues there. Where a read
// ends inside the overlap (or within W of its edge), that read end must be
// shared by other reads (capTerminalPiles): amplicon reads end together, while
// a shotgun read that happens to end exactly at a repeat boundary ends alone.
// Both tests are needed. Flanks alone let through 2–9 wrong joins per 30–40 kb
// shotgun set with exact repeats, from reads ending within W of a repeat edge;
// shared ends alone would pass two amplicons that both end inside one repeat.
function capCoverageStep(A, B, flip, aln, P, threshold) {
    var W = Math.max(P.kmer, P.maxOverhang);
    // [from, to) on the read as aligned; mirrored onto the read's own strand when
    // it was aligned reverse-complemented, since mult and pile are stored forward.
    var ends = [
        { read: A, flank: aln.sStart, from: aln.sStart - W, to: aln.sStart, rc: false, end: 0 },
        { read: A, flank: A.len - aln.sEnd, from: aln.sEnd, to: aln.sEnd + W, rc: false, end: 1 },
        { read: B, flank: aln.tStart, from: aln.tStart - W, to: aln.tStart, rc: !!flip, end: flip ? 1 : 0 },
        { read: B, flank: B.len - aln.tEnd, from: aln.tEnd, to: aln.tEnd + W, rc: !!flip, end: flip ? 0 : 1 }
    ];
    for (var i = 0; i < ends.length; i++) {
        var e = ends[i];
        if (e.flank > W) {
            var from = e.rc ? e.read.len - e.to : e.from;
            var to = e.rc ? e.read.len - e.from : e.to;
            if (capRepeatFraction(e.read, from, to, threshold) >= P.repeatFrac) return false;
        } else if (!e.read.pile || e.read.pile[e.end] < CAP_PILE_MIN) {
            return false;
        }
    }
    return true;
}

/* ───────────────────────────── 5. Overlap computation ────────────────────── */

// The local identity below which a stretch inside an overlap is a chimeric
// junction: the value set, or 35 points under minIdentity, never below 40.
function capJunctionIdentity(P) {
    return P.junctionIdentity > 0 ? P.junctionIdentity : Math.max(40, P.minIdentity - 35);
}

/**
 * Align one candidate pair and decide whether it is an acceptable overlap.
 * Returns an overlap record, or {rejected:'…'} describing why it was refused.
 *
 * Coordinate convention for the record:
 *   bSeq  = flip ? reverseComplement(reads[b].seq) : reads[b].seq
 *   the alignment covers reads[a].seq[aStart,aEnd) against bSeq[bStart,bEnd)
 *   offset = aStart − bStart  — the position of bSeq[0] in read a's frame,
 *   which is exactly what the layout stage needs and can be negative.
 */
function capOverlapPair(reads, cand, P, repeatThreshold) {
    var A = reads[cand.a], B = reads[cand.b];
    var bSeq = cand.flip ? capRC(B.seq) : B.seq;
    // Only chase a wider band when this already looks like an overlap that could
    // be accepted. Unrelated reads that share a chance k-mer align at roughly the
    // identity of random sequence and never qualify, so widening for them is pure
    // waste — and there are far more of them than there are real overlaps.
    var widenGate = function (r) {
        return r.identity >= (P.minIdentity - 5) &&
               Math.max(r.sEnd - r.sStart, r.tEnd - r.tStart) >= P.minOverlap * 0.5;
    };
    var aln = capAlignAuto(A.seq, bSeq, cand.diag, P, false, widenGate);
    if (!aln) return { a: cand.a, b: cand.b, flip: cand.flip, rejected: 'no alignment' };

    // ── the repeat filter ──
    // In a genuine suffix–prefix overlap one read runs off each edge of the
    // alignment, so at each edge at least one of the two unaligned overhangs is
    // ~zero. Two reads that share only an interspersed repeat both continue past
    // the repeat on both sides, so min(overhang) is large at both edges.
    //
    // The alignment itself cannot show this. It is end-space free, so it always
    // runs out to the end of one read on each side — for a repeat pair, straight
    // on through the unrelated sequence beyond the repeat — and min(overhang)
    // measured on it is zero by construction. The overhangs are measured on the
    // alignment's maximal-scoring segment instead, which ends where the reads
    // stop agreeing. Only a segment that would itself qualify as an overlap is
    // judged here: a chance seed between unrelated reads also leaves long
    // unaligned ends, but it fails on length or identity below, and that is the
    // right reason to give for it.
    var coreLen = Math.max(aln.coreSEnd - aln.coreSStart, aln.coreTEnd - aln.coreTStart);
    var leftOver = Math.min(aln.coreSStart, aln.coreTStart);
    var rightOver = Math.min(A.len - aln.coreSEnd, bSeq.length - aln.coreTEnd);
    if (coreLen >= P.minOverlap && aln.coreIdentity >= P.minIdentity &&
        (leftOver > P.maxOverhang || rightOver > P.maxOverhang)) {
        return { a: cand.a, b: cand.b, flip: cand.flip, olen: coreLen, identity: aln.coreIdentity,
                 leftOverhang: leftOver, rightOverhang: rightOver,
                 // The shared stretch, for capProvenRepeats (b on its oriented copy).
                 aStart: aln.coreSStart, aEnd: aln.coreSEnd, bStart: aln.coreTStart, bEnd: aln.coreTEnd,
                 rejected: 'unmatched overhang ' + Math.max(leftOver, rightOver) + ' nt > maxOverhang (' +
                           P.maxOverhang + ') — typical of a shared repeat rather than a true join' };
    }

    var olen = Math.max(aln.sEnd - aln.sStart, aln.tEnd - aln.tStart);
    if (olen < P.minOverlap) {
        return { a: cand.a, b: cand.b, flip: cand.flip, olen: olen, identity: aln.identity,
                 rejected: 'overlap too short (' + olen + ' < ' + P.minOverlap + ')' };
    }
    if (aln.identity < P.minIdentity) {
        return { a: cand.a, b: cand.b, flip: cand.flip, olen: olen, identity: aln.identity,
                 rejected: 'identity too low (' + aln.identity.toFixed(1) + '% < ' + P.minIdentity + '%)' };
    }

    // ── the junction filter ──
    // Two chimeric reads made of the same two products in the same order share
    // one product on one side of their junctions and the other product on the
    // other side, with a stretch between the two junctions where one read has
    // the first product and the other the second: unrelated sequence. Judged
    // as a whole the alignment passes — the stretch is a small part of it — and
    // the overhang test sees nothing either, since the matching flanks on both
    // sides keep the stretch inside the maximal segment. Judged locally it
    // cannot pass: no stretch of CAP_JUNCTION_WINDOW columns between two reads
    // of one molecule falls to the identity of random sequence, and the gap
    // run that the shift between the two junctions forces into the alignment
    // counts against the window as well. Read errors stay well clear of the
    // threshold (see junctionIdentity in CAP_DEFAULTS), so the test costs one
    // genuine overlap in tens of thousands at R9 error levels and none at R10.
    // A junction at the very end of an overlap is caught too, when the
    // overhang test let it through as shorter than maxOverhang.
    if (P.chimeraCheck) {
        var jId = capJunctionIdentity(P);
        if (aln.coreMinWindow < jId) {
            return { a: cand.a, b: cand.b, flip: cand.flip, olen: olen, identity: aln.identity,
                     rejected: 'identity falls to ' + aln.coreMinWindow.toFixed(0) + '% over ' + CAP_JUNCTION_WINDOW +
                               ' columns inside the overlap (< ' + jId + '%) — the reads agree on either side ' +
                               'of that stretch but not across it, as at a chimeric junction' };
        }
    }

    // ── the second repeat filter ──
    // The overhang test above catches the easy case, where the two reads diverge
    // on both sides of a shared repeat. It cannot catch the hard one: a read
    // ending inside a repeat and another beginning inside the same repeat look
    // exactly like a clean suffix–prefix join, and accepting it splices two
    // unrelated loci together. What distinguishes them is that the alignment
    // covers nothing but repeat sequence. An overlap that reaches into
    // single-copy sequence at either end is trustworthy; one that does not is a
    // guess, and it is exactly the guess that produces a confident, wrong contig.
    //
    // High copy number alone does not make a repeat, though. Where two tiled PCR
    // amplicons overlap, the overlap is read by both amplicons' reads, so its
    // words are seen twice as often as the rest — exactly the signature of a
    // two-copy repeat — and every join between neighbouring amplicons lies
    // wholly inside it. Five 1.5 kb amplicons with 350 bp overlaps at equal
    // depth came out as three contigs. Two things tell such a junction from a
    // repeat: the high-copy stretch stops where the overlap stops (it was the
    // other amplicon's reads that raised it), and the reads that end inside the
    // overlap end together with their siblings. Such an overlap is accepted
    // (capAmpliconJoin).
    if (isFinite(repeatThreshold)) {
        var fa = capRepeatFraction(A, aln.sStart, aln.sEnd, repeatThreshold);
        // tStart/tEnd are on the ORIENTED copy of B; B.mult is indexed on B's own
        // forward strand, so a flipped read's interval has to be mirrored back.
        var bFrom = cand.flip ? (B.len - aln.tEnd) : aln.tStart;
        var bTo = cand.flip ? (B.len - aln.tStart) : aln.tEnd;
        var fb = capRepeatFraction(B, bFrom, bTo, repeatThreshold);
        if (fa >= P.repeatFrac && fb >= P.repeatFrac &&
            !capAmpliconJoin(A, B, cand.flip, aln, P, repeatThreshold)) {
            return { a: cand.a, b: cand.b, flip: cand.flip, olen: olen, identity: aln.identity,
                     repeatFraction: Math.min(fa, fb),
                     rejected: 'the aligned region is ' + Math.round(100 * Math.min(fa, fb)) +
                               '% high-copy repeat sequence — this join cannot be told from the ' +
                               'same repeat in another location' };
        }
    }

    var offset = aln.sStart - aln.tStart;
    var contained = null;
    if (offset <= P.maxOverhang && offset + bSeq.length >= A.len - P.maxOverhang) contained = 'a_in_b';
    else if (offset >= -P.maxOverhang && offset + bSeq.length <= A.len + P.maxOverhang) contained = 'b_in_a';

    return {
        a: cand.a, b: cand.b, flip: cand.flip,
        aStart: aln.sStart, aEnd: aln.sEnd, bStart: aln.tStart, bEnd: aln.tEnd,
        offset: offset, olen: olen,
        matches: aln.matches, mismatches: aln.mismatches, gaps: aln.gaps,
        identity: aln.identity, score: aln.score,
        leftOverhang: leftOver, rightOverhang: rightOver,
        contained: contained, seeds: cand.seeds, rejected: null
    };
}

/**
 * Align a shard of candidate pairs. This is the one stage worth parallelising:
 * on any realistic input it dominates the wall-clock time (the banded DP is
 * O(pairs × readLen × band)), and every pair is completely independent of every
 * other — no shared state, no ordering requirement. js/cap-worker.js runs this
 * function and nothing else.
 *
 * Purely synchronous so it can be called straight from a Worker's onmessage.
 */
function capComputeOverlaps(reads, pairs, P, repeatThreshold, from, to) {
    var lo = from || 0, hi = to === undefined ? pairs.length : to;
    var out = [], rejected = [];
    for (var i = lo; i < hi; i++) {
        var ov = capOverlapPair(reads, pairs[i], P, repeatThreshold);
        if (ov.rejected) rejected.push(ov); else out.push(ov);
    }
    return { overlaps: out, rejected: rejected };
}

/**
 * The default overlap runner: same thread, chunked so the page can repaint.
 * The page replaces this with a Worker-pool runner (js/cap-pool.js) when
 * Workers are available and the job is big enough to be worth the hand-off.
 */
async function capSerialOverlapRunner(reads, pairs, P, ctx, repeatThreshold) {
    var out = [], rejected = [];
    for (var i = 0; i < pairs.length; i++) {
        var ov = capOverlapPair(reads, pairs[i], P, repeatThreshold);
        if (ov.rejected) rejected.push(ov); else out.push(ov);
        if (ctx && (i & 63) === 0) await ctx.tick('overlaps', i, pairs.length);
    }
    return { overlaps: out, rejected: rejected };
}

/* Chimeric reads.
 *
 * One read can carry two molecules: two amplicons ligated end to end, a
 * fold-back read that runs on into the reverse complement of its own template,
 * or two molecules read back to back through the pore with an adapter between
 * them. Each half overlaps its own reads perfectly well, so the read passes every
 * overlap filter and then joins two unrelated places in the layout. Measured on
 * five tiled 1.5 kb amplicons, two such reads among 202 were enough to turn the
 * 6.1 kb target into one wrong contig.
 *
 * The signature is inside the read. The reads that overlap the first molecule
 * stop where it stops, the reads that overlap the second one start there, and
 * hardly any read runs across the junction, because only another chimera with
 * the same junction could. So a junction is called where
 *   • at least CAP_CHIMERA_MIN overlapping reads end within a few bases of one
 *     another, and at least as many begin within a few bases of one another,
 *     from W before the ends to CAP_CHIMERA_GAP after them (room for an
 *     adapter's worth of sequence between the two molecules);
 *   • those that end are at least half of all the reads overlapping the read
 *     just before the junction, and those that begin at least half of those just
 *     after it: the whole left side stops, the whole right side starts;
 *   • the reads that span the junction number no more than a quarter of the
 *     smaller pile, which lets a few chimeras sharing one junction, and so
 *     spanning it for one another, still be caught.
 * Piles alone are not enough: at a coverage hole, or among 20x short reads, a
 * few reads end and a few begin together by chance, and with only the first and
 * third rules a shotgun sweep set aside genuine reads. Inside an amplicon read
 * the siblings span every position, so none of this can fire there.
 *
 * Returns [{read, at, to, ending, starting, spanning}], at/to on the read's own
 * forward strand. */
var CAP_CHIMERA_GAP = 200;
var CAP_CHIMERA_MIN = 3;

// Number of sorted values in [lo, hi].
function capCountIn(arr, lo, hi) {
    var a = 0, b = arr.length, m;
    while (a < b) { m = (a + b) >> 1; if (arr[m] < lo) a = m + 1; else b = m; }
    var first = a;
    b = arr.length;
    while (a < b) { m = (a + b) >> 1; if (arr[m] <= hi) a = m + 1; else b = m; }
    return a - first;
}

function capChimeras(reads, overlaps, P) {
    var W = Math.max(P.kmer, P.maxOverhang);
    var edge = Math.max(P.minOverlap, W);
    var iv = new Array(reads.length);
    var add = function (r, s, e) { (iv[r] || (iv[r] = [])).push([s, e]); };
    for (var i = 0; i < overlaps.length; i++) {
        var o = overlaps[i], Lb = reads[o.b].len;
        add(o.a, o.aStart, o.aEnd);
        add(o.b, o.flip ? Lb - o.bEnd : o.bStart, o.flip ? Lb - o.bStart : o.bEnd);
    }
    var countIn = capCountIn;
    var found = [];
    for (var r = 0; r < reads.length; r++) {
        var list = iv[r];
        if (!list || list.length < 2 * CAP_CHIMERA_MIN) continue;
        var L = reads[r].len, ends = [], starts = [];
        for (var j = 0; j < list.length; j++) {
            if (list[j][1] >= edge && list[j][1] <= L - edge) ends.push(list[j][1]);
            if (list[j][0] >= edge && list[j][0] <= L - edge) starts.push(list[j][0]);
        }
        if (ends.length < CAP_CHIMERA_MIN || starts.length < CAP_CHIMERA_MIN) continue;
        ends.sort(function (x, y) { return x - y; });
        starts.sort(function (x, y) { return x - y; });
        // A pile is the reads ending (or starting) within 2 × CAP_PILE_JITTER of the
        // first of them. Reads that stop because their molecule stops do so within a
        // few bases of one another; the looser W let the repeat filter fake a pile,
        // since near a repeat it keeps only the overlaps that end before a certain
        // point, and three such ends fell within ±W at 4 % error.
        var J2 = 2 * CAP_PILE_JITTER, hit = null;
        for (j = 0; j < ends.length && !hit; j++) {
            var at = ends[j];
            var ending = countIn(ends, at, at + J2);
            if (ending < CAP_CHIMERA_MIN) continue;
            // The densest pile of starts from just before the ends to an adapter's
            // length after them.
            var bestTo = -1, starting = 0;
            for (var k = 0; k < starts.length; k++) {
                if (starts[k] < at - W || starts[k] > at + CAP_CHIMERA_GAP) continue;
                var c = countIn(starts, starts[k], starts[k] + J2);
                if (c > starting) { starting = c; bestTo = starts[k]; }
            }
            if (starting < CAP_CHIMERA_MIN) continue;
            var lo = Math.min(at, bestTo) - W, hi = Math.max(at, bestTo) + W, spanning = 0;
            var pL = Math.max(0, lo - W), pR = Math.min(L - 1, hi + W), leftDepth = 0, rightDepth = 0;
            for (k = 0; k < list.length; k++) {
                if (list[k][0] < lo && list[k][1] > hi) spanning++;
                if (list[k][0] <= pL && list[k][1] > pL) leftDepth++;
                if (list[k][0] <= pR && list[k][1] > pR) rightDepth++;
            }
            if (4 * spanning <= Math.min(ending, starting) &&
                2 * ending >= leftDepth && 2 * starting >= rightDepth) {
                hit = { read: r, at: at, to: Math.max(at, bestTo), ending: ending, starting: starting, spanning: spanning };
            }
        }
        if (hit) found.push(hit);
    }
    return found;
}

/* Chimeras proven by parting.
 *
 * The pile test sees only accepted overlaps, and a ligation chimera has few of
 * those. The reads of either molecule align with it up to the junction and then
 * go on with their own sequence, which the overhang filter refuses as it would a
 * shared repeat; what it accepts are reads that lie wholly inside one half, and
 * through those the chimera still reaches the layout, now as a container whose
 * other half hangs off the contig. Measured on 400 reads of two 2 kb products
 * with 10 % chimeras: five junction chimeras, each holding the kept middle of a
 * fold-back read, built two contigs extended by 366 and 676 nt.
 *
 * The refusals themselves are the evidence. On the chimera every partner parts
 * at the junction and no accepted overlap runs across it; on a read of either
 * molecule the same partings fall at that molecule's position, but its siblings
 * cross there. A junction is called where
 *   • at least CAP_CHIMERA_MIN partners part within 2 × CAP_PILE_JITTER of one
 *     another, their agreement with the read ending there (they vouch for the
 *     read before the point and contradict it after) or beginning there;
 *   • those partings are at least half of all the overlaps, accepted or refused,
 *     that cover the read just before (or after) the point: the whole side parts;
 *   • the accepted overlaps crossing the point number no more than a quarter of
 *     the partings, as in the pile test, so chimeras sharing a junction and
 *     crossing it for one another are still caught.
 * One side is enough. When the other molecule has no reads in the set nothing
 * can begin at the junction, but the read's continuation is still contradicted
 * by every partner and vouched for by none. A read of one amplicon has none of
 * this: its siblings cross every position.
 *
 * Returns records shaped as capChimeras's, with evidence 'parting'. */
function capPartingChimeras(reads, overlaps, rejected, P) {
    var W = Math.max(P.kmer, P.maxOverhang), J2 = 2 * CAP_PILE_JITTER;
    var acc = new Array(reads.length), cover = new Array(reads.length);
    var pEnd = new Array(reads.length), pStart = new Array(reads.length);
    var push = function (arr, r, v) { (arr[r] || (arr[r] = [])).push(v); };
    var i, o, Lb, bFrom, bTo;
    for (i = 0; i < overlaps.length; i++) {
        o = overlaps[i]; Lb = reads[o.b].len;
        push(acc, o.a, [o.aStart, o.aEnd]);
        push(acc, o.b, [o.flip ? Lb - o.bEnd : o.bStart, o.flip ? Lb - o.bStart : o.bEnd]);
    }
    for (i = 0; i < rejected.length; i++) {
        o = rejected[i];
        if (o.aStart === undefined || o.leftOverhang === undefined) continue;
        var left = o.leftOverhang > P.maxOverhang, right = o.rightOverhang > P.maxOverhang;
        Lb = reads[o.b].len;
        bFrom = o.flip ? Lb - o.bEnd : o.bStart; bTo = o.flip ? Lb - o.bStart : o.bEnd;
        push(cover, o.a, [o.aStart, o.aEnd]);
        push(cover, o.b, [bFrom, bTo]);
        // The partner agrees with a over [aStart, aEnd) and parts at the edge where
        // both reads go on: at aStart the agreement begins, at aEnd it ends. On b's
        // own strand a flipped overlap swaps the two.
        if (left) { push(pStart, o.a, o.aStart); push(o.flip ? pEnd : pStart, o.b, o.flip ? bTo : bFrom); }
        if (right) { push(pEnd, o.a, o.aEnd); push(o.flip ? pStart : pEnd, o.b, o.flip ? bFrom : bTo); }
    }
    var num = function (x, y) { return x - y; };
    var found = [];
    for (var r = 0; r < reads.length; r++) {
        var E = pEnd[r] || [], S = pStart[r] || [];
        if (E.length < CAP_CHIMERA_MIN && S.length < CAP_CHIMERA_MIN) continue;
        E.sort(num); S.sort(num);
        var L = reads[r].len, A = acc[r] || [], C = cover[r] || [], j, k, c, lo, hi, spanning;
        var depthAt = function (p) {
            var n = 0;
            for (var q = 0; q < A.length; q++) if (A[q][0] <= p && A[q][1] > p) n++;
            for (q = 0; q < C.length; q++) if (C[q][0] <= p && C[q][1] > p) n++;
            return n;
        };
        var crossing = function (from, to) {
            var n = 0;
            for (var q = 0; q < A.length; q++) if (A[q][0] < from && A[q][1] > to) n++;
            return n;
        };
        var hit = null;
        // Agreement ends here: the read's continuation to the right is contradicted.
        for (j = 0; j < E.length && !hit; j++) {
            var at = E[j], ending = capCountIn(E, at, at + J2);
            if (ending < CAP_CHIMERA_MIN) continue;
            var bestTo = at, starting = 0;
            for (k = 0; k < S.length; k++) {
                if (S[k] < at - W || S[k] > at + CAP_CHIMERA_GAP) continue;
                c = capCountIn(S, S[k], S[k] + J2);
                if (c > starting) { starting = c; bestTo = S[k]; }
            }
            lo = Math.min(at, bestTo) - W; hi = Math.max(at, bestTo) + W;
            spanning = crossing(lo, hi);
            if (4 * spanning > ending || 2 * ending < depthAt(Math.max(0, lo - W))) continue;
            hit = { read: r, at: at, to: Math.max(at, bestTo), ending: ending, starting: starting,
                    spanning: spanning, evidence: 'parting' };
        }
        // Agreement begins here: the continuation to the left is contradicted.
        for (j = S.length - 1; j >= 0 && !hit; j--) {
            var to = S[j], starting2 = capCountIn(S, to - J2, to);
            if (starting2 < CAP_CHIMERA_MIN) continue;
            var bestAt = to, ending2 = 0;
            for (k = 0; k < E.length; k++) {
                if (E[k] > to + W || E[k] < to - CAP_CHIMERA_GAP) continue;
                c = capCountIn(E, E[k] - J2, E[k]);
                if (c > ending2) { ending2 = c; bestAt = E[k]; }
            }
            lo = Math.min(to, bestAt) - W; hi = Math.max(to, bestAt) + W;
            spanning = crossing(lo, hi);
            if (4 * spanning > starting2 || 2 * starting2 < depthAt(Math.min(L - 1, hi + W))) continue;
            hit = { read: r, at: Math.min(to, bestAt), to: to, ending: ending2, starting: starting2,
                    spanning: spanning, evidence: 'parting' };
        }
        if (hit) found.push(hit);
    }
    return found;
}

/* Repeats proven by divergence.
 *
 * The copy-number test in capOverlapPair can only see a repeat whose words are
 * counted clearly more often than ordinary sequence. A repeat present twice is
 * counted twice as often only on average: at ~9x the two copies together were
 * seen 12-17 times against a cut-off of 17, so the test rarely fired, and a
 * 30 kb shotgun set with a 1 kb repeat at two loci assembled with 63 joins
 * between the copies.
 *
 * The same data hold direct proof. Two reads from different copies align over
 * the repeat and then go separate ways, and that pair is refused by the overhang
 * filter (124 such pairs in that set, every one between the two copies). The
 * stretch they share is therefore present at two places. Such a refusal counts
 * as proof only when
 *   • the reads part over at least CAP_REPEAT_FLANK bases, longer than the
 *     adapter or barcode left on a read end, and the two continuations share
 *     no more than one word (a noisy stretch that ended the alignment early
 *     still shares many), and
 *   • on each read, the point where they part is crossed by an accepted overlap
 *     with some other read: the read's own continuation is real. A chimeric read
 *     parts from its partners at its junction too, but nothing crosses that.
 * The words of every proven stretch then mark that sequence on every read
 * (read.proven), and an accepted overlap lying at least repeatFrac inside proven
 * repeat on both reads, without max(minOverlap, 2k) shared bases outside it, is
 * refused like any other repeat join.
 *
 * Returns the refused overlaps' records. */
var CAP_REPEAT_FLANK = 150;

// Call fn(code, pos) for each word of seq[from, to) — canonical 2-bit code, as in
// seeding, so a word and its reverse complement agree.
function capEachWord(seq, from, to, K, fn) {
    var mask = K === 15 ? 0x3FFFFFFF : ((1 << (2 * K)) - 1), rcShift = 2 * (K - 1);
    var code = 0, rcode = 0, valid = 0;
    for (var p = Math.max(0, from); p < Math.min(to, seq.length); p++) {
        var c = CAP_CODE[seq.charCodeAt(p)];
        if (c < 0) { valid = 0; code = 0; rcode = 0; continue; }
        code = ((code << 2) | c) & mask;
        rcode = ((rcode >>> 2) | ((3 - c) << rcShift)) & mask;
        if (++valid >= K) fn(code < rcode ? code : rcode, p - K + 1);
    }
}

function capProvenRepeats(reads, overlaps, rejected, P) {
    var W = Math.max(P.kmer, P.maxOverhang), K = P.kmer, i;
    var iv = new Array(reads.length);
    var add = function (r, s, e) { (iv[r] || (iv[r] = [])).push(s, e); };
    for (i = 0; i < overlaps.length; i++) {
        var o = overlaps[i], Lb = reads[o.b].len;
        add(o.a, o.aStart, o.aEnd);
        add(o.b, o.flip ? Lb - o.bEnd : o.bStart, o.flip ? Lb - o.bStart : o.bEnd);
    }
    var crossed = function (r, pos) {
        var L = iv[r];
        if (!L) return false;
        for (var j = 0; j < L.length; j += 2) if (L[j] <= pos - W && L[j + 1] >= pos + W) return true;
        return false;
    };

    // Do the two continuations really differ? Unrelated flanks share no words; a
    // pair whose alignment merely ran into a noisy stretch shares many. At 4 %
    // error one such pair made a whole amplicon "proven repeat".
    var flanksDiffer = function (rj, right) {
        var A = reads[rj.a], B = reads[rj.b];
        var bSeq = rj.flip ? capRC(B.seq) : B.seq;
        var fa = right ? A.seq.substring(rj.aEnd, rj.aEnd + CAP_REPEAT_FLANK)
                       : A.seq.substring(Math.max(0, rj.aStart - CAP_REPEAT_FLANK), rj.aStart);
        var fb = right ? bSeq.substring(rj.bEnd, rj.bEnd + CAP_REPEAT_FLANK)
                       : bSeq.substring(Math.max(0, rj.bStart - CAP_REPEAT_FLANK), rj.bStart);
        var seen = new Set(), shared = 0;
        capEachWord(fa, 0, fa.length, K, function (code) { seen.add(code); });
        capEachWord(fb, 0, fb.length, K, function (code) { if (seen.has(code)) shared++; });
        return shared <= 1;
    };

    var marks = [];
    for (i = 0; i < rejected.length; i++) {
        var rj = rejected[i];
        if (rj.aStart === undefined || rj.leftOverhang === undefined) continue;
        var left = rj.leftOverhang >= CAP_REPEAT_FLANK, right = rj.rightOverhang >= CAP_REPEAT_FLANK;
        if (!left && !right) continue;
        var lb = reads[rj.b].len;
        // b's shared stretch on its own strand; flipped, its left end is on the right.
        var bFrom = rj.flip ? lb - rj.bEnd : rj.bStart, bTo = rj.flip ? lb - rj.bStart : rj.bEnd;
        if (left && !(crossed(rj.a, rj.aStart) && crossed(rj.b, rj.flip ? bTo : bFrom) && flanksDiffer(rj, false))) continue;
        if (right && !(crossed(rj.a, rj.aEnd) && crossed(rj.b, rj.flip ? bFrom : bTo) && flanksDiffer(rj, true))) continue;
        marks.push(rj.a, rj.aStart, rj.aEnd, rj.b, bFrom, bTo);
    }
    if (!marks.length) return [];

    var words = new Set();
    var keep = function (code) { words.add(code); };
    for (i = 0; i < marks.length; i += 3) capEachWord(reads[marks[i]].seq, marks[i + 1], marks[i + 2], K, keep);
    reads.forEach(function (rd) {
        var m = null;
        capEachWord(rd.seq, 0, rd.len, K, function (code, p) {
            if (!words.has(code)) return;
            if (!m) m = new Uint8Array(rd.len);
            for (var q = p; q < p + K; q++) m[q] = 1;
        });
        rd.proven = m;
    });

    var need = Math.max(P.minOverlap, 2 * P.kmer);
    // Fraction of [from, to) in proven repeat, and the longest run outside it.
    var inRepeat = function (rd, from, to) {
        var n = 0, run = 0, best = 0;
        for (var p = from; p < to; p++) {
            if (rd.proven[p]) { n++; run = 0; } else if (++run > best) best = run;
        }
        return { frac: to > from ? n / (to - from) : 0, single: best };
    };
    // An amplicon junction inside a repeat: every read end that lies within the
    // overlap is shared by other reads (capTerminalPiles), as PCR products end
    // together and a shotgun read ends alone. Without this, tiled products whose
    // overlap fell inside a repeat copy found elsewhere were split at that
    // junction, which the products themselves resolve. "Shared" here means by at
    // least a quarter as many reads as the read has overlaps: at 20x, two or
    // three shotgun reads end together by chance, and that alone let a join
    // between the copies of a 1 kb repeat through.
    var piledJoin = function (ov, A, B) {
        if (!A.pile || !B.pile) return false;
        var ends = 0;
        var check = function (inner, r, end) {
            if (!inner) return true;
            ends++;
            var p = reads[r].pile[end], deg = iv[r] ? iv[r].length / 2 : 0;
            return p >= CAP_PILE_MIN && 4 * p >= deg;
        };
        return check(ov.aStart <= W, ov.a, 0) && check(A.len - ov.aEnd <= W, ov.a, 1) &&
               check(ov.bStart <= W, ov.b, ov.flip ? 1 : 0) && check(B.len - ov.bEnd <= W, ov.b, ov.flip ? 0 : 1) &&
               ends > 0;
    };
    var refused = [], kept = [];
    for (i = 0; i < overlaps.length; i++) {
        var ov = overlaps[i], A = reads[ov.a], B = reads[ov.b];
        if (A.proven && B.proven && !piledJoin(ov, A, B)) {
            var LB = B.len;
            var ra = inRepeat(A, ov.aStart, ov.aEnd);
            var rb = inRepeat(B, ov.flip ? LB - ov.bEnd : ov.bStart, ov.flip ? LB - ov.bStart : ov.bEnd);
            if (ra.frac >= P.repeatFrac && rb.frac >= P.repeatFrac && !(ra.single >= need && rb.single >= need)) {
                refused.push({ a: ov.a, b: ov.b, flip: ov.flip, olen: ov.olen, identity: ov.identity,
                               repeatFraction: Math.min(ra.frac, rb.frac),
                               rejected: 'the aligned region is ' + Math.round(100 * Math.min(ra.frac, rb.frac)) +
                                         '% sequence that other reads show at a second place (they align over it ' +
                                         'and then part) — this join cannot be told from that copy' });
                continue;
            }
        }
        kept.push(ov);
    }
    overlaps.length = 0;
    Array.prototype.push.apply(overlaps, kept);
    return refused;
}

/* ───────────────────────────── 6. Layout ─────────────────────────────────── */

/* Greedy overlap-layout. Reads are grouped into components; every read in a
   component carries a position and a strand flag in that component's coordinate
   frame. Overlaps are consumed best-first; one that agrees with the current
   layout merges two components (or is redundant), one that disagrees is recorded
   as a conflict and dropped.

   Components hold their members explicitly rather than using a union-find with
   composed offsets. Offsets do NOT compose by addition once a strand flip is
   involved — flipping a component maps a member at (pos, flip) with length L to
   (−(pos + L), flip ^ 1) — and doing that lazily inside find() is the classic
   source of silent orientation bugs. Merging the smaller component into the
   larger keeps the total work O(n log n), which is irrelevant at these sizes
   but costs nothing. */

/* Where does the partner in overlap `o` sit relative to read `r`, and which end of
   `r` does it extend past?
   Returns {end, lo, hi} with end = -1 (neither), 0 (r's 5′ end) or 1 (r's 3′ end),
   and [lo, hi) the partner's interval — both expressed in r's OWN FORWARD frame.

   The frame matters. `o.offset` places read b's ORIENTED copy in read a's forward
   frame, so from b's point of view read a starts at −offset in that same oriented
   frame — which is b's reverse complement whenever o.flip is set. Reporting an end
   without mirroring back would file the same physical end of b under '5′' or '3′'
   depending on which partner happened to flip it, and the branch guard would then
   compare partners that are not on the same end at all. */
function capExtendedEnd(reads, o, r, tol) {
    var La = reads[o.a].len, Lb = reads[o.b].len;
    var lo, hi, L;
    if (r === o.a) {
        lo = o.offset; hi = o.offset + Lb; L = La;
    } else {
        lo = -o.offset; hi = -o.offset + La; L = Lb;
        if (o.flip) {                       // mirror the oriented-b frame back onto b
            var m = Lb - hi;
            hi = Lb - lo;
            lo = m;
        }
    }
    var end = -1;
    if (hi > L + tol) end = 1;
    else if (lo < -tol) end = 0;
    return { end: end, lo: lo, hi: hi };
}

/* Repeat-branch detection.
 *
 * A read whose 3′ end runs into a repeat has several partners that all extend
 * past that end — one per copy of the repeat. Reads from the SAME copy overlap
 * one another; reads from DIFFERENT copies do not, because outside the repeat
 * they are unrelated sequence. So: take every partner that extends a given end
 * of a read, and ask whether each pair of them overlaps. When two partners are
 * predicted to overlap substantially (their placements relative to the read say
 * so) but no such overlap survived the filters, the read cannot tell the copies
 * apart and that end is AMBIGUOUS.
 *
 * The layout then refuses to extend a contig through an ambiguous end. The
 * assembly stops at the repeat boundary and says so, instead of chaining two
 * unrelated loci into one confident-looking, wrong contig — which is precisely
 * what happens without this check.
 *
 * Returns a Set of "readIdx:end" keys.
 */
function capBranchGuard(reads, overlaps, P) {
    var tol = P.maxOverhang;
    var have = new Set();
    var byEnd = new Map();
    var i, o;

    for (i = 0; i < overlaps.length; i++) {
        o = overlaps[i];
        have.add(o.a + ':' + o.b);
        var ea = capExtendedEnd(reads, o, o.a, tol);
        if (ea.end >= 0) push(o.a * 2 + ea.end, { p: o.b, lo: ea.lo, hi: ea.hi, score: o.score });
        var eb = capExtendedEnd(reads, o, o.b, tol);
        if (eb.end >= 0) push(o.b * 2 + eb.end, { p: o.a, lo: eb.lo, hi: eb.hi, score: o.score });
    }
    function push(k, v) {
        var arr = byEnd.get(k);
        if (!arr) { arr = []; byEnd.set(k, arr); }
        arr.push(v);
    }

    var ambiguous = new Set();
    byEnd.forEach(function (arr, k) {
        if (arr.length < 2) return;

        // Weak partners would create phantom clusters, and at very high coverage
        // the pairwise pass would grow quadratically for no benefit. Judge the
        // question on the strong partners only.
        var best = 0, j;
        for (j = 0; j < arr.length; j++) if (arr[j].score > best) best = arr[j].score;
        var strong = arr.filter(function (v) { return v.score >= 0.6 * best; });
        if (strong.length > 20) {
            strong = strong.slice().sort(function (x, y) { return y.score - x.score || x.p - y.p; }).slice(0, 20);
        }
        if (strong.length < 2) return;

        // Cluster the partners by "these two have an overlap of their own".
        // A single missing edge inside an otherwise connected group is just a
        // filtered-out pair; it must not condemn the end. Only a genuine SPLIT —
        // two groups with no edge between them at all — means the read cannot
        // tell one continuation from another.
        var parent = strong.map(function (_, idx) { return idx; });
        var find = function (x) { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
        var x2, y2;
        for (x2 = 0; x2 < strong.length; x2++) {
            for (y2 = x2 + 1; y2 < strong.length; y2++) {
                var lo2 = Math.min(strong[x2].p, strong[y2].p), hi2 = Math.max(strong[x2].p, strong[y2].p);
                if (have.has(lo2 + ':' + hi2)) parent[find(x2)] = find(y2);
            }
        }
        // Is any cross-cluster pair PREDICTED to overlap? If two partners are not
        // expected to share sequence, the absence of an overlap between them says
        // nothing; it is only evidence when they should have had one.
        for (x2 = 0; x2 < strong.length; x2++) {
            for (y2 = x2 + 1; y2 < strong.length; y2++) {
                if (find(x2) === find(y2)) continue;
                var predicted = Math.min(strong[x2].hi, strong[y2].hi) - Math.max(strong[x2].lo, strong[y2].lo);
                if (predicted >= P.minOverlap) {
                    ambiguous.add((k >> 1) + ':' + (k & 1));
                    return;
                }
            }
        }
    });
    return ambiguous;
}

// Where does read b sit, given read a's placement, in a's component frame?
// `offset` and `flip` are read b's placement relative to read a's FORWARD strand.
function capPlace(aPlace, aLen, bLen, offset, flip) {
    if (!aPlace.flip) return { pos: aPlace.pos + offset, flip: flip };
    // a is stored reverse-complemented in this frame, so the whole relationship
    // is mirrored: the interval [offset, offset+bLen) in a's forward frame maps
    // to [aPos + aLen − offset − bLen, aPos + aLen − offset).
    return { pos: aPlace.pos + aLen - offset - bLen, flip: flip ^ 1 };
}

async function capLayout(reads, overlaps, P, ctx) {
    var n = reads.length, i;

    // component id per read; components[cid] = {id, members: Map(read -> {pos,flip})}
    var compOf = new Int32Array(n);
    var comps = new Array(n);
    for (i = 0; i < n; i++) {
        compOf[i] = i;
        comps[i] = { id: i, members: new Map([[i, { pos: 0, flip: 0 }]]), alive: true };
    }

    // Best-first, with a total order so the result never depends on sort stability.
    var sorted = overlaps.slice().sort(function (p, q) {
        if (q.score !== p.score) return q.score - p.score;
        if (q.olen !== p.olen) return q.olen - p.olen;
        if (q.identity !== p.identity) return q.identity - p.identity;
        if (p.a !== q.a) return p.a - q.a;
        if (p.b !== q.b) return p.b - q.b;
        return p.flip - q.flip;
    });

    var accepted = [], conflicts = [];
    var ambiguous = P.repeatGuard ? capBranchGuard(reads, overlaps, P) : new Set();

    for (i = 0; i < sorted.length; i++) {
        var ov = sorted[i];
        var ca = comps[compOf[ov.a]], cb = comps[compOf[ov.b]];
        var pa = ca.members.get(ov.a);
        var want = capPlace(pa, reads[ov.a].len, reads[ov.b].len, ov.offset, ov.flip);

        if (ca === cb) {
            // Both reads are already placed. The overlap is either redundant
            // confirmation or evidence that the layout is wrong here.
            var pb = cb.members.get(ov.b);
            var dPos = Math.abs(pb.pos - want.pos);
            if (pb.flip !== want.flip) {
                ov.role = 'conflict';
                conflicts.push({ ov: ov, reason: 'orientation contradicts the existing layout' });
            } else if (dPos > P.offsetTolerance) {
                ov.role = 'conflict';
                conflicts.push({ ov: ov, reason: 'offset differs from the layout by ' + dPos +
                                                 ' nt (> offsetTolerance ' + P.offsetTolerance + ')' });
            } else {
                ov.redundant = true;
                ov.role = 'confirming';
                accepted.push(ov);
            }
            continue;
        }

        // Joining two contigs through an end that the branch guard flagged would
        // be a guess. Record it and leave the contigs apart — a reported break is
        // recoverable, a silent mis-join is not.
        if (ambiguous.size) {
            var ea = capExtendedEnd(reads, ov, ov.a, P.maxOverhang).end;
            var eb = capExtendedEnd(reads, ov, ov.b, P.maxOverhang).end;
            if ((ea >= 0 && ambiguous.has(ov.a + ':' + ea)) ||
                (eb >= 0 && ambiguous.has(ov.b + ':' + eb))) {
                ov.role = 'conflict';
                conflicts.push({ ov: ov, reason: 'ambiguous extension — this end overlaps reads that do ' +
                                                 'not overlap each other, the signature of a repeat boundary' });
                continue;
            }
        }

        // Merge cb into ca (or the smaller into the larger, transforming whichever
        // one moves). Compute the transform that maps cb's frame onto ca's.
        var pbCur = cb.members.get(ov.b);
        var flipComp = (want.flip !== pbCur.flip) ? 1 : 0;
        var anchorPos = flipComp ? -(pbCur.pos + reads[ov.b].len) : pbCur.pos;
        var delta = want.pos - anchorPos;

        var src, dst, srcFlip, srcDelta;
        if (cb.members.size <= ca.members.size) {
            src = cb; dst = ca; srcFlip = flipComp; srcDelta = delta;
        } else {
            // Move ca instead: invert the transform. cb→ca is x ↦ f(x)+delta, so
            // ca→cb is x ↦ f(x − delta) because f is an involution.
            src = ca; dst = cb; srcFlip = flipComp;
            srcDelta = flipComp ? delta : -delta;
        }
        var members = src.members;
        members.forEach(function (pl, readIdx) {
            var np;
            if (srcFlip) np = { pos: -(pl.pos + reads[readIdx].len) + srcDelta, flip: pl.flip ^ 1 };
            else np = { pos: pl.pos + srcDelta, flip: pl.flip };
            dst.members.set(readIdx, np);
            compOf[readIdx] = dst.id;
        });
        src.alive = false;
        src.members = null;
        ov.role = 'joining';
        accepted.push(ov);
        if (ctx && (i & 255) === 0) await ctx.tick('layout', i, sorted.length);
    }

    // ── materialise contigs ──
    var groups = new Map();
    for (i = 0; i < n; i++) {
        var cid = compOf[i];
        if (!groups.has(cid)) groups.set(cid, []);
        groups.get(cid).push(i);
    }

    var contigs = [], singletons = [];
    var ids = Array.from(groups.keys()).sort(function (x, y) { return x - y; });
    for (var gi = 0; gi < ids.length; gi++) {
        var idxs = groups.get(ids[gi]);
        if (idxs.length === 1) { singletons.push(idxs[0]); continue; }
        var comp = comps[ids[gi]];
        var placed = idxs.map(function (ri) {
            var pl = comp.members.get(ri);
            return { readIdx: ri, pos: pl.pos, flip: pl.flip, len: reads[ri].len };
        });
        // Normalise so the leftmost read starts at 0, then order by position.
        var min = Infinity;
        placed.forEach(function (p) { if (p.pos < min) min = p.pos; });
        placed.forEach(function (p) { p.pos -= min; });
        placed.sort(function (p, q) { return p.pos - q.pos || p.readIdx - q.readIdx; });

        // A component is only a contig if its reads actually form a continuous
        // chain. They always should — every merge came from an overlap — but a
        // conflicting layout can leave a gap, and a contig with a hole in it
        // would produce nonsense consensus. Split at any gap.
        var pieces = [], cur = [placed[0]], reach = placed[0].pos + placed[0].len;
        for (var k = 1; k < placed.length; k++) {
            if (placed[k].pos >= reach) { pieces.push(cur); cur = [placed[k]]; reach = placed[k].pos + placed[k].len; }
            else { cur.push(placed[k]); reach = Math.max(reach, placed[k].pos + placed[k].len); }
        }
        pieces.push(cur);

        for (var pi = 0; pi < pieces.length; pi++) {
            var piece = pieces[pi];
            if (piece.length === 1) { singletons.push(piece[0].readIdx); continue; }
            var base = piece[0].pos, span = 0;
            piece.forEach(function (p) { p.pos -= base; span = Math.max(span, p.pos + p.len); });
            contigs.push({ id: 0, reads: piece, span: span, split: pieces.length > 1 });
        }
    }

    // Largest contig first, ties broken by the smallest read index for determinism.
    contigs.sort(function (p, q) { return q.span - p.span || p.reads[0].readIdx - q.reads[0].readIdx; });
    contigs.forEach(function (c, ci) { c.id = ci + 1; });
    singletons.sort(function (x, y) { return x - y; });

    return { contigs: contigs, singletons: singletons, accepted: accepted, conflicts: conflicts };
}

/* ─────────────────── 7. Column alignment & consensus per contig ──────────── */

/* The layout gives each read an approximate offset. That is not enough to call a
   consensus: indels shift the reads relative to one another, so the reads must be
   brought into a genuine column-wise multiple alignment first.

   Reads are added in position order to a growing column list. Each read is
   aligned to the CURRENT consensus by the same banded semi-global kernel; where
   the read carries an insertion relative to the consensus, a fresh column is
   created and every previously placed read simply has no entry there.

   Columns are sparse (a Map read→cell) because coverage is a handful of reads
   even when the contig is long: a dense reads × columns matrix would be mostly
   empty and, for a few hundred reads, large.

   A cell packs the base together with the position it came from on the read's
   OWN forward strand, so the consensus caller can look the base's PHRED value up
   in the read's quality array and can tell how close to a read end it sat. Both
   matter: an unclipped Sanger tail is exactly where the spurious bases are. */

function capNewColumn() { return new Map(); }

// cell = readPos * 128 + charCode. readPos is on the read's forward strand.
function capCell(ch, alignedIdx, flip, readLen) {
    var pos = flip ? (readLen - 1 - alignedIdx) : alignedIdx;
    if (pos < 0) pos = 0;
    return pos * 128 + ch.charCodeAt(0);
}
function capCellChar(cell) { return String.fromCharCode(cell & 127); }
// Math.floor, not `>> 7`: a cell is an ordinary JS number, so it stays exact well
// past 2^31, but a bitwise shift would truncate it to 32 bits and start returning
// nonsense (or negatives) for positions beyond ~16 Mb into a read.
function capCellPos(cell) { return Math.floor(cell / 128); }

// Majority base of a column, ignoring reads that do not cover it. Used only as
// the alignment reference while building; the reported consensus is computed
// once at the end with the full quality/polymorphism logic.
//
// Called once per column per read added, so on a long contig this runs millions
// of times: it counts into four local variables rather than allocating a tally
// object, which is what it used to do.
function capQuickBase(col) {
    var na = 0, nc = 0, ng = 0, nt = 0, gapN = 0, ambN = 0, tot = 0;
    col.forEach(function (cell) {
        var code = cell & 127;
        tot++;
        if (code === 97) na++;
        else if (code === 99) nc++;
        else if (code === 103) ng++;
        else if (code === 116) nt++;
        else if (code === 45) gapN++;
        else ambN++;                        // an IUPAC code carried by the read
    });
    if (!tot) return '';
    if (gapN * 2 > tot) return '';          // '' = the column is mostly gap
    var best = '', bestN = 0;
    if (na > bestN) { bestN = na; best = 'a'; }
    if (nc > bestN) { bestN = nc; best = 'c'; }
    if (ng > bestN) { bestN = ng; best = 'g'; }
    if (nt > bestN) { bestN = nt; best = 't'; }
    // A column of nothing but ambiguity codes is still a column: give it an 'n' so
    // it stays part of the reference the next read is aligned against, rather than
    // silently disappearing from the coordinate system.
    return bestN ? best : (ambN ? 'n' : '');
}

// The column `n` consensus bases away from column `from` (backwards when n < 0),
// clamped to [0, cols.length]. Columns that are mostly gap carry no consensus base
// and are stepped over without being counted, which is what turns a distance in
// bases — a layout offset, a read length — into one in columns.
function capWalkBases(cols, from, n) {
    var c = Math.max(0, Math.min(cols.length, from)), seen = 0;
    if (n >= 0) {
        while (c < cols.length && seen < n) { if (capQuickBase(cols[c])) seen++; c++; }
    } else {
        while (c > 0 && seen < -n) { c--; if (capQuickBase(cols[c])) seen++; }
    }
    return c;
}

function capConsensusString(cols) {
    var out = [];
    for (var i = 0; i < cols.length; i++) {
        var b = capQuickBase(cols[i]);
        if (b) out.push(b);
    }
    return out.join('');
}

/**
 * Put a read's extra bases (cells, in read order) into columns that other reads
 * opened at the same place (slots, in column order), or say where a new column
 * is needed. A small alignment: a base placed in a slot whose bases agree with
 * it scores 2, in one whose bases differ 0, a new column -1, and a slot may be
 * left unused. Returns one {slot, after} per cell: the slot's column, or -1 and
 * the column the new one follows (-1 = before the first slot).
 */
function capPlaceInSlots(cells, slots, cols) {
    var k = cells.length, m = slots.length, i, j;
    var maj = slots.map(function (c) {
        var n = { a: 0, c: 0, g: 0, t: 0 }, best = '', bn = 0;
        cols[c].forEach(function (cell) { var ch = capCellChar(cell); if (n[ch] !== undefined) n[ch]++; });
        ['a', 'c', 'g', 't'].forEach(function (b) { if (n[b] > bn) { bn = n[b]; best = b; } });
        return best;
    });
    var W = m + 1, S = new Float64Array((k + 1) * W), T = new Uint8Array((k + 1) * W);
    for (i = 0; i <= k; i++) {
        for (j = 0; j <= m; j++) {
            if (!i && !j) continue;
            var best = -Infinity, how = 0;
            if (i && j) { var v = S[(i - 1) * W + j - 1] + (capCellChar(cells[i - 1]) === maj[j - 1] ? 2 : 0); if (v > best) { best = v; how = 1; } }
            if (j) { var u = S[i * W + j - 1]; if (u > best) { best = u; how = 2; } }       // slot unused
            if (i) { var w = S[(i - 1) * W + j] - 1; if (w > best) { best = w; how = 3; } } // new column
            S[i * W + j] = best; T[i * W + j] = how;
        }
    }
    var out = new Array(k);
    i = k; j = m;
    while (i > 0) {
        var t = T[i * W + j];
        if (t === 1) { out[i - 1] = { slot: slots[j - 1], after: -1 }; i--; j--; }
        else if (t === 2) j--;
        else { out[i - 1] = { slot: -1, after: j > 0 ? slots[j - 1] : -1 }; i--; }
    }
    return out;
}

/**
 * Build the column-wise multiple alignment of one contig.
 * @returns {{cols:Array<Map>, rows:Array}}  rows = [{readIdx, flip, first, last}]
 */
async function capContigAlign(reads, contig, P, ctx) {
    var placed = contig.reads;
    var cols = [];
    var rows = [];
    // Layout coordinates and COLUMN coordinates drift apart as gap columns are
    // inserted: every insertion pushes the columns one further along than the
    // layout says. Assuming the two agree (they do only for the first read) makes
    // a long contig fail to stack at all — each read that misses its search band
    // is appended instead of aligned, which pushes the next one further out
    // still, and the contig ends up several times its true length.
    //
    // So each read is anchored on the one placed before it: the previous read's
    // ACTUAL first column plus the layout distance between them. That is a local
    // measurement, re-taken every iteration, so a single bad placement cannot
    // poison the reads after it.
    var lastFirstCol = -1, lastLayoutPos = 0;
    // Indels accumulate along a contig, so give the consensus stage a wider band
    // than the pairwise stage uses.
    var consBand = Math.max(P.bandWidth, 48);
    var consP = capParams(P);
    consP.bandWidth = consBand;
    consP.maxBandWidth = Math.max(P.maxBandWidth, 1024);

    for (var r = 0; r < placed.length; r++) {
        if (ctx) await ctx.tick('aligning contig', r, placed.length);
        var p = placed[r];
        var rlen = reads[p.readIdx].len;
        var seq = p.flip ? capRC(reads[p.readIdx].seq) : reads[p.readIdx].seq;

        if (!cols.length) {
            for (var i = 0; i < seq.length; i++) {
                var c = capNewColumn();
                c.set(p.readIdx, capCell(seq.charAt(i), i, p.flip, rlen));
                cols.push(c);
            }
            rows.push({ readIdx: p.readIdx, flip: p.flip, first: 0, last: cols.length - 1 });
            lastFirstCol = 0; lastLayoutPos = p.pos;
            continue;
        }

        // Build the reference the read will be aligned against — but only over a
        // WINDOW of columns around where the read is expected to land. Rebuilding
        // the whole consensus for every read costs O(reads × columns × coverage),
        // which on a 30 kb contig is tens of millions of column tallies and
        // dominated the run; the read can only align near its layout position, so
        // the rest of the contig is wasted work.
        //
        // The window is anchored on the PREVIOUS read's actual first column (the
        // nearest thing to a local reference) and is generous enough to absorb the
        // drift between layout and column coordinates.
        //
        // Every distance is counted in consensus BASES, not columns. Each read adds
        // columns for its own insertions, so at depth a contig holds far more
        // columns than bases (2.4k for a 1.5 kb amplicon at 200 ONT reads). Sized in
        // columns, the window stopped short of the read's end, the unaligned tail
        // was appended as new columns, the next window fell shorter still, and the
        // contig ran away: 47.9k columns at 400 reads, out of memory at 1000, and
        // 58 consensus errors on a 4 kb amplicon at 200 reads.
        var margin = Math.max(512, 2 * consBand);
        var expCol = lastFirstCol >= 0 ? capWalkBases(cols, lastFirstCol, p.pos - lastLayoutPos)
                                       : capWalkBases(cols, 0, p.pos);
        var winLo = capWalkBases(cols, expCol, -margin);
        var winHi = cols.length;

        // colOfCons maps a consensus index (within the window) to an ABSOLUTE
        // column index, so the alignment result still lands on the right columns.
        // The reference runs from `lo` until `need` consensus bases lie at or past
        // expCol; winHi is left at the first column not read.
        var colOfCons = [];
        var consChars = [];
        var expConsIdx = -1;
        var buildRef = function (lo, need) {
            colOfCons.length = 0; consChars.length = 0; expConsIdx = -1;
            var ci = lo, after = 0;
            for (; ci < cols.length && after < need; ci++) {
                if (expConsIdx < 0 && ci >= expCol) expConsIdx = consChars.length;
                var b = capQuickBase(cols[ci]);
                if (b) {
                    colOfCons.push(ci); consChars.push(b);
                    if (ci >= expCol) after++;
                }
            }
            winHi = ci;
        };
        buildRef(winLo, seq.length + margin);
        // The window can land on a stretch of gap-only columns and yield no
        // reference at all. Widen to the whole contig before concluding anything;
        // only a contig with no consensus-bearing column ANYWHERE is genuinely
        // degenerate. (Resetting on the windowed result would throw away every
        // read already placed — silent data loss for a purely local problem.)
        if (!consChars.length && (winLo > 0 || winHi < cols.length)) {
            winLo = 0;
            buildRef(0, Infinity);
        }
        var cons = consChars.join('');
        if (!cons.length) {
            // Nothing to align against. Start the column list from this read
            // rather than losing it; the reads placed so far had no usable
            // consensus, so there is nothing to preserve.
            cols.length = 0; r--; continue;
        }
        if (expConsIdx < 0) expConsIdx = consChars.length;

        // Expected diagonal j − i: read base 0 should sit at consensus index
        // expConsIdx within the window.
        var k0 = expConsIdx;
        if (k0 > cons.length) k0 = cons.length;
        if (k0 < -seq.length) k0 = -seq.length;
        var aln = capAlignAuto(seq, cons, k0, consP, true);
        // A read that overlaps by less than the layout promised means the band
        // was centred wrongly (a big local indel, or an accumulated drift the
        // running correction has not caught up with). One wide retry is far
        // cheaper than appending the read and derailing every read after it.
        var expected = Math.min(seq.length, cons.length - Math.max(0, k0));
        if (aln && expected >= P.minOverlap &&
            (aln.sEnd - aln.sStart) < Math.min(expected, P.minOverlap)) {
            var wide = capParams(consP);
            wide.bandWidth = Math.min(4096, Math.max(512, Math.ceil(cons.length / 8)));
            wide.maxBandWidth = wide.bandWidth;
            var retry = capAlignBanded(seq, cons, k0, wide.bandWidth, wide, true);
            if (retry && (retry.sEnd - retry.sStart) > (aln.sEnd - aln.sStart)) aln = retry;
        }
        if (!aln) {
            // Should not happen (the layout only joined reads that overlapped),
            // but never drop a read silently: place it by its layout offset.
            aln = { sStart: 0, sEnd: 0, tStart: Math.max(0, Math.min(cons.length, p.pos)),
                    tEnd: Math.max(0, Math.min(cons.length, p.pos)), ops: '' };
        }

        // Walk the alignment, writing this read's bases into existing columns and
        // collecting insertions. Insertions are applied in a single rebuild pass
        // afterwards, so adding a read is O(columns) rather than O(columns × indels).
        var inserts = [];             // {afterCol, chars: []}
        var i2 = aln.sStart;          // index into seq
        var j2 = aln.tStart;          // index into cons
        var ops = aln.ops || '';
        var pendingIns = null;
        var firstCol = -1;

        for (var oi = 0; oi < ops.length; oi++) {
            var op = ops.charAt(oi);
            if (op === 'M') {
                var col = colOfCons[j2];
                cols[col].set(p.readIdx, capCell(seq.charAt(i2), i2, p.flip, rlen));
                if (firstCol < 0) firstCol = col;
                i2++; j2++; pendingIns = null;
            } else if (op === 'I') {
                // The read has a base the consensus does not: an insertion, so a
                // new column has to be created for it.
                var after = j2 > aln.tStart ? colOfCons[j2 - 1] : (colOfCons[j2] !== undefined ? colOfCons[j2] - 1 : -1);
                if (!pendingIns || pendingIns.afterCol !== after) {
                    pendingIns = { afterCol: after, cells: [] };
                    inserts.push(pendingIns);
                }
                pendingIns.cells.push(capCell(seq.charAt(i2), i2, p.flip, rlen));
                i2++;
            } else {
                // gap in the read: the consensus has a base the read does not.
                // The read position recorded is the next base it will contribute,
                // which is what the quality lookup should use for this deletion.
                var dcol = colOfCons[j2];
                cols[dcol].set(p.readIdx, capCell('-', Math.min(i2, rlen - 1), p.flip, rlen));
                if (firstCol < 0) firstCol = dcol;
                j2++; pendingIns = null;
            }
        }

        // Bases of the read outside the aligned region become brand-new columns.
        // With free end gaps the alignment can only start at the reference's
        // beginning or the read's, so a head exists exactly when the read runs off
        // the left of the aligned region — and it belongs immediately before the
        // first aligned column, not at the front of the whole contig (which is
        // where a windowed reference would otherwise put it).
        var headLen = aln.sStart, tailLen = seq.length - aln.sEnd;
        var head = headLen > 0 ? seq.substring(0, headLen) : '';
        var tail = tailLen > 0 ? seq.substring(aln.sEnd) : '';
        var firstAlignedCol = colOfCons[aln.tStart];
        var lastAlignedCol = aln.tEnd > aln.tStart ? colOfCons[aln.tEnd - 1] : firstAlignedCol;
        if (firstAlignedCol === undefined) firstAlignedCol = winLo;
        if (lastAlignedCol === undefined) lastAlignedCol = firstAlignedCol;

        var myFirst = firstCol >= 0 ? firstCol : firstAlignedCol;
        var byAfter = null;

        if (inserts.length || head || tail) {
            // Group every new column by the existing column it follows, keeping
            // path order within a group: head, then the insertions in the order
            // the alignment produced them, then the tail.
            byAfter = new Map();
            var addAfter = function (afterCol, cell) {
                var bucket = byAfter.get(afterCol);
                if (!bucket) { bucket = []; byAfter.set(afterCol, bucket); }
                bucket.push(cell);
            };
            for (var hi2 = 0; hi2 < head.length; hi2++) {
                addAfter(firstAlignedCol - 1, capCell(head.charAt(hi2), hi2, p.flip, rlen));
            }
            for (var ii = 0; ii < inserts.length; ii++) {
                var ins = inserts[ii];
                for (var cci = 0; cci < ins.cells.length; cci++) addAfter(ins.afterCol, ins.cells[cci]);
            }
            for (var ti = 0; ti < tail.length; ti++) {
                addAfter(lastAlignedCol, capCell(tail.charAt(ti), aln.sEnd + ti, p.flip, rlen));
            }

            // Earlier reads may already have opened columns at the same place for
            // the same extra bases. Those columns are out of the reference while
            // most reads lack the base, so this read could not align to them; if
            // it opened columns of its own, a base that most reads carry but the
            // reference happens to lack would be split over one column per read,
            // never win a column, and drop out of the consensus. So the extra
            // bases go into those columns first (capPlaceInSlots), and only the
            // remainder open new ones.
            var groups = Array.from(byAfter.keys());
            for (var gk = 0; gk < groups.length; gk++) {
                var gAfter = groups[gk], gCells = byAfter.get(gAfter), slots = [];
                for (var sc = gAfter + 1; sc < cols.length && !cols[sc].has(p.readIdx) && !capQuickBase(cols[sc]); sc++) {
                    slots.push(sc);
                }
                if (!slots.length) continue;
                byAfter.delete(gAfter);
                var placed2 = capPlaceInSlots(gCells, slots, cols);
                for (var pc = 0; pc < gCells.length; pc++) {
                    if (placed2[pc].slot >= 0) cols[placed2[pc].slot].set(p.readIdx, gCells[pc]);
                    else addAfter(placed2[pc].after < 0 ? gAfter : placed2[pc].after, gCells[pc]);
                }
            }
        }

        if (byAfter && byAfter.size) {
            // A new column lies inside the span of every read already placed on
            // both sides of it, and for those reads it is a deletion. Writing that
            // gap now, not only in the final pass, is what keeps one read's
            // insertion out of the reference the next read is aligned against:
            // capQuickBase ignores reads that have no cell in a column, so a new
            // column holding a single base read as solid consensus. With
            // insertion-rich reads every read added its errors to the reference,
            // the reference outgrew the template, and later reads were aligned
            // against the excess. Measured on 50 R9-level reads of a 1.5 kb
            // amplicon: the reference grew from 1497 to 1748 bases and the
            // consensus kept 68 spurious insertions near one end.
            var rebuilt = [];
            myFirst = -1;
            var spanners = null;
            var pushCell = function (cell) {
                var nc = capNewColumn();
                nc.set(p.readIdx, cell);
                for (var si = 0; si < spanners.length; si += 2) nc.set(spanners[si], spanners[si + 1]);
                if (myFirst < 0) myFirst = rebuilt.length;
                rebuilt.push(nc);
            };
            var flushAt = function (afterCol) {
                var g = byAfter.get(afterCol);
                if (!g) return;
                // Reads with a cell (base or gap) on both sides: [readIdx, gapCell, …].
                spanners = [];
                var left = afterCol >= 0 ? cols[afterCol] : null;
                var right = afterCol + 1 < cols.length ? cols[afterCol + 1] : null;
                if (left && right) {
                    var small = left.size <= right.size ? left : right, other = small === left ? right : left;
                    small.forEach(function (cell, ri) {
                        if (ri === p.readIdx || !other.has(ri)) return;
                        // The gap carries the read position of the base before it,
                        // as the final pass does, for the quality lookup.
                        spanners.push(ri, capCellPos(left.get(ri)) * 128 + 45);
                    });
                }
                for (var gj = 0; gj < g.length; gj++) pushCell(g[gj]);
            };
            flushAt(-1);
            for (var ck = 0; ck < cols.length; ck++) {
                if (myFirst < 0 && cols[ck].has(p.readIdx)) myFirst = rebuilt.length;
                rebuilt.push(cols[ck]);
                flushAt(ck);
            }
            cols = rebuilt;
        }

        // Complete this read's votes now: a column inside its span that it has no
        // cell in (one opened earlier for another read's extra base, which this
        // read was never aligned to) is a deletion for it. With this, and with the
        // gaps written for spanning reads above, every column holds a vote from
        // every read across it, so capQuickBase is a true majority throughout.
        // The read wrote ops.length + head + tail cells; the scan stops once it
        // has met them all.
        var toMeet = (aln.ops || '').length + headLen + tailLen, met = 0, gapPos = 0, pend = [];
        for (var fc = myFirst >= 0 ? myFirst : 0; fc < cols.length && met < toMeet; fc++) {
            var own = cols[fc].get(p.readIdx);
            if (own === undefined) { if (met) pend.push(fc); continue; }
            met++;
            for (var pj = 0; pj < pend.length; pj++) cols[pend[pj]].set(p.readIdx, gapPos * 128 + 45);
            pend.length = 0;
            gapPos = capCellPos(own);
        }
        rows.push({ readIdx: p.readIdx, flip: p.flip, first: 0, last: 0 });
        lastFirstCol = myFirst >= 0 ? myFirst : 0;
        lastLayoutPos = p.pos;
    }

    // Recompute each read's covered span now that all columns exist, and fill
    // internal gaps so "covered but no base" reads as a deletion rather than as
    // absence of coverage.
    var spanFirst = new Map(), spanLast = new Map();
    for (var cc = 0; cc < cols.length; cc++) {
        cols[cc].forEach(function (cell, ri) {
            if (!spanFirst.has(ri)) spanFirst.set(ri, cc);
            spanLast.set(ri, cc);
        });
    }
    for (var ri2 = 0; ri2 < rows.length; ri2++) {
        var rr = rows[ri2];
        rr.first = spanFirst.has(rr.readIdx) ? spanFirst.get(rr.readIdx) : 0;
        rr.last = spanLast.has(rr.readIdx) ? spanLast.get(rr.readIdx) : -1;
        // Columns inside a read's span that it never wrote to are deletions, not
        // absence of coverage. They inherit the read position of the preceding
        // base so the quality lookup still lands somewhere sensible.
        var lastPos = 0;
        for (var q2 = rr.first; q2 <= rr.last; q2++) {
            var have = cols[q2].get(rr.readIdx);
            if (have === undefined) cols[q2].set(rr.readIdx, lastPos * 128 + 45);   // 45 = '-'
            else lastPos = capCellPos(have);
        }
    }
    rows.sort(function (x, y) { return x.first - y.first || x.readIdx - y.readIdx; });
    return { cols: cols, rows: rows };
}

/**
 * Call the consensus of an aligned contig, encoding real polymorphism as IUPAC.
 *
 * A column becomes ambiguous only when a minority base is supported by
 *   • at least polyMinCount reads, AND
 *   • at least polyMinFrac percent of the column's non-gap coverage, AND
 *   • the column has at least polyMinCov coverage.
 * At 1–2% sequencing error, seeing the SAME wrong base in ≥2 of ≥3 reads at the
 * SAME column is already unlikely; requiring a quarter of the reads makes it
 * negligible, while a true heterozygous site (≈50/50) is caught comfortably.
 * Indel polymorphism has no IUPAC letter, so it is reported in the variant table
 * with kind='indel' instead of being folded into the consensus.
 */
/* P(X >= k) for X ~ Binomial(n, p). Coverage is small, so the terms are summed
   directly and iteratively — no logs, no approximation, no overflow. */
function capBinomTailGE(n, k, p) {
    if (k <= 0) return 1;
    if (k > n) return 0;
    if (p <= 0) return 0;
    if (p >= 1) return 1;
    var q = 1 - p, term = Math.pow(q, n), sum = 0, i;
    for (i = 0; i < k; i++) {
        sum += term;
        term = term * ((n - i) / (i + 1)) * (p / q);
    }
    var res = 1 - sum;
    return res < 0 ? 0 : (res > 1 ? 1 : res);
}

// When the bases in a column are equally supported, the one lying further inside
// its read wins: read ends are where base calls are poorest, and where a primer
// sits. Joining PCR products, the primer-derived end of one product overlaps the
// template-derived middle of its neighbour, and a primer mismatch showed up as a
// two-read tie that the alphabet decided — half the primer mismatches went into
// the consensus. Distance counts up to this many bases; beyond it a base is as
// good as any other, and a tie that remains is written as an IUPAC code.
var CAP_TIE_DEPTH = 50;

function capContigConsensus(reads, aligned, P, model) {
    model = model || {};
    // Probability that one read shows a PARTICULAR wrong base at a position. The
    // error rate is measured from the assembly's own overlaps rather than assumed
    // (see capAssemble), and divided by three because an error has three ways to
    // go wrong and only one of them mimics the minority allele in question.
    var pErr = (model.errorRate || 0.01) / 3;
    // Significance level. Left at "auto" it is 1/(total alignment columns), i.e.
    // a Bonferroni correction targeting about ONE false polymorphic call per
    // assembly however long it is. Without it, the plain fraction rule scales
    // badly: on a 30 kb assembly at 8x it produced 23 calls where the sequence
    // had no polymorphism at all, because two coincident errors in eight reads
    // reach 25% of the column.
    var alpha = model.alpha || 1e-6;
    var cols = aligned.cols;
    var consensus = [], columns = [], variants = [];
    var covSum = 0;

    // Weight of one observed base. With a quality file, PHRED is used directly
    // (it is already a log-odds score, so summing PHRED across reads is the
    // right way to combine independent evidence). Without one, every base
    // weighs 1. Either way, bases within `endPenalty` of a read end are ramped
    // down linearly: read ends are where the base calls are least reliable, and
    // in a 2-read column they are what decides a tie.
    var weightOf = function (ri, pos) {
        var w = 1;
        var rd = reads[ri];
        if (P.qualWeight && rd.qual && pos < rd.qual.length) w = Math.max(1, rd.qual[pos]);
        if (P.endPenalty > 0) {
            var dist = Math.min(pos, rd.len - 1 - pos);
            if (dist < P.endPenalty) w *= (dist + 1) / (P.endPenalty + 1);
        }
        return w;
    };

    for (var ci = 0; ci < cols.length; ci++) {
        var col = cols[ci];
        var counts = { a: 0, c: 0, g: 0, t: 0 };
        var wcounts = { a: 0, c: 0, g: 0, t: 0 };
        var inner = { a: 0, c: 0, g: 0, t: 0 };   // tie-break support, see CAP_TIE_DEPTH
        var gap = 0, other = 0, cov = 0, otherBits = 0;

        col.forEach(function (cell, ri) {
            var ch = capCellChar(cell);
            var cpos = capCellPos(cell);
            var w = weightOf(ri, cpos);
            cov++;
            if (ch === '-') { gap++; return; }
            if (counts[ch] !== undefined) {
                counts[ch]++; wcounts[ch] += w;
                inner[ch] += Math.min(cpos, reads[ri].len - 1 - cpos, CAP_TIE_DEPTH);
                return;
            }
            // An IUPAC code in the READ itself (the input may legitimately contain
            // R, Y, N …). It is not evidence for one base, but it is evidence
            // against the bases it excludes, so its weight is split across the
            // bases it does allow. Its bits are also remembered so that a column
            // where EVERY read is ambiguous still produces a base: dropping it
            // would silently shorten the consensus and shift every coordinate
            // after it.
            var bits = CAP_BITS[ch];
            if (!bits) { other++; return; }
            otherBits |= bits;
            other++;
            var n = 0, bb;
            for (bb = 0; bb < 4; bb++) if (bits & (1 << bb)) n++;
            if (!n) return;
            if (bits & 1) wcounts.a += w / n;
            if (bits & 2) wcounts.c += w / n;
            if (bits & 4) wcounts.g += w / n;
            if (bits & 8) wcounts.t += w / n;
        });

        var nonGap = counts.a + counts.c + counts.g + counts.t;
        covSum += cov;

        if (!nonGap) {
            // No read gave an unambiguous base here. If some read gave an
            // ambiguous one, the column still exists and must still produce a
            // base — the union of what the reads allow — otherwise a run of Ns in
            // the input would vanish from the consensus and shorten it.
            if (otherBits && gap * 2 <= cov) {
                var amb = CAP_FROM_BITS[otherBits];
                columns.push({ col: ci, pos: consensus.length, coverage: cov, gap: gap, other: other,
                               counts: shallowCounts(counts), base: amb, isPoly: false, alleles: [] });
                consensus.push(amb);
                continue;
            }
            // Genuinely nothing here (all gap): keep the column for the alignment
            // view, contribute nothing to the sequence.
            columns.push({ col: ci, pos: -1, coverage: cov, gap: gap, other: other, counts: shallowCounts(counts),
                           base: '', isPoly: false, alleles: [] });
            continue;
        }
        if (gap > cov / 2) {
            // Most reads have no base here: the insertion is a minority event.
            // Drop the column from the consensus but record it as an indel variant
            // when enough reads support it.
            if (nonGap >= P.polyMinCount && cov >= P.polyMinCov &&
                (100 * nonGap / cov) >= P.polyMinFrac) {
                variants.push({ col: ci, pos: consensus.length, kind: 'indel',
                                coverage: cov, counts: shallowCounts(counts), gap: gap,
                                code: '', note: 'insertion present in ' + nonGap + '/' + cov + ' reads' });
            }
            columns.push({ col: ci, pos: -1, coverage: cov, gap: gap, other: other, counts: shallowCounts(counts),
                           base: '', isPoly: false, alleles: [] });
            continue;
        }

        // Rank the observed bases by weighted support, then raw counts, then by
        // how far inside their reads they lie, and finally alphabetically so the
        // same input always gives the same order (no dependence on Map iteration
        // order or sort stability).
        var order = ['a', 'c', 'g', 't'].filter(function (b) { return counts[b] > 0; });
        order.sort(function (x, y) {
            if (wcounts[y] !== wcounts[x]) return wcounts[y] - wcounts[x];
            if (counts[y] !== counts[x]) return counts[y] - counts[x];
            if (inner[y] !== inner[x]) return inner[y] - inner[x];
            return x < y ? -1 : 1;
        });
        // Bases that tie with the leader on all three: nothing in the data
        // prefers one, so none is picked (see below).
        var tieBits = 0;
        for (var ti2 = 1; ti2 < order.length; ti2++) {
            var tb = order[ti2];
            if (wcounts[tb] === wcounts[order[0]] && counts[tb] === counts[order[0]] && inner[tb] === inner[order[0]]) {
                tieBits |= CAP_BITS[tb];
            }
        }

        var major = order[0];
        var bits = CAP_BITS[major];
        var alleles = [{ base: major, count: counts[major] }];
        var isPoly = false;
        var minP = 1;
        if (cov >= P.polyMinCov) {
            for (var oi = 1; oi < order.length; oi++) {
                var b = order[oi];
                var frac = 100 * counts[b] / nonGap;
                if (counts[b] < P.polyMinCount || frac < P.polyMinFrac) continue;
                // Could this many reads carry this base by sequencing error alone?
                var pv = capBinomTailGE(nonGap, counts[b], pErr);
                if (pv < minP) minP = pv;
                if (pv > alpha) continue;
                bits |= CAP_BITS[b];
                alleles.push({ base: b, count: counts[b] });
                isPoly = true;
            }
        }
        // A tie that the polymorphism test did not already turn into a code is
        // written as the IUPAC code of the tied bases and listed with kind 'tie':
        // with two reads that disagree there is no telling which is right, and
        // picking the alphabetically first one, as was done before, put an
        // arbitrary base into the consensus without a trace.
        var tied = !isPoly && tieBits;
        var call = isPoly ? CAP_FROM_BITS[bits] : tied ? CAP_FROM_BITS[CAP_BITS[major] | tieBits] : major;
        var pos = consensus.length;
        consensus.push(call);

        if (isPoly) {
            variants.push({ col: ci, pos: pos, kind: 'snp', coverage: cov,
                            counts: shallowCounts(counts), gap: gap, code: call,
                            pValue: minP, alleles: alleles });
        } else if (tied) {
            variants.push({ col: ci, pos: pos, kind: 'tie', coverage: cov,
                            counts: shallowCounts(counts), gap: gap, code: call,
                            note: 'reads split evenly, the bases equally far inside their reads' });
        }
        // Gaps present in a minority of reads are a deletion polymorphism.
        if (gap > 0 && gap >= P.polyMinCount && cov >= P.polyMinCov &&
            (100 * gap / cov) >= P.polyMinFrac) {
            variants.push({ col: ci, pos: pos, kind: 'indel', coverage: cov,
                            counts: shallowCounts(counts), gap: gap, code: call,
                            note: 'deletion in ' + gap + '/' + cov + ' reads' });
        }

        columns.push({ col: ci, pos: pos, coverage: cov, gap: gap, other: other,
                       counts: shallowCounts(counts), base: call,
                       isPoly: isPoly, alleles: alleles });
    }

    var trimmed = capTrimEnds(consensus, columns, variants, P);
    // Mean coverage of what is left: the cut tails were the thin part.
    var covN = cols.length;
    if (trimmed[0] || trimmed[1]) {
        covSum = 0; covN = 0;
        for (var ki = 0; ki < columns.length; ki++) {
            if (columns[ki].trimmed) continue;
            covSum += columns[ki].coverage; covN++;
        }
    }

    return {
        consensus: consensus.join(''),
        columns: columns,
        variants: variants,
        trimmed: trimmed,
        meanCoverage: covN ? covSum / covN : 0
    };
}

function shallowCounts(c) { return { a: c.a, c: c.c, g: c.g, t: c.t }; }

/* Contig ends by coverage.
 *
 * Whatever reaches the layout with an unmatched tail — a fold-back copy shorter
 * than foldbackMin, a junction within maxOverhang of a read's end, adapter left
 * on a read — is laid beyond the last base the other reads support, and the
 * consensus runs on through it on the word of one or two reads. The end of a
 * contig is therefore cut back to the last column carried by at least
 * minEndCov reads, or by CAP_END_COV_FRAC of the contig's median coverage when
 * that is more (at 200 reads a tail carried by five of them is still a tail).
 *
 * Only a thin stretch longer than maxOverhang is cut. Reads of one PCR product
 * end within a few bases of one another, not at one base, so the coverage
 * tapers over the last bases of a genuine end and the final base may be
 * carried by few reads; cutting that taper would shorten every product by a
 * base or two. A tail is hundreds of bases; the taper is not.
 *
 * Applied to the consensus array in place: the trimmed columns lose their
 * consensus position (they stay in the alignment view, as columns with no
 * consensus base), the kept columns and variants are renumbered, and variants
 * in the trimmed part — the ambiguity codes two disagreeing tails produce — go.
 * Returns [bases cut from the start, bases cut from the end]. */
var CAP_END_COV_FRAC = 0.05;

function capTrimEnds(consensus, columns, variants, P) {
    var cut = [0, 0], i;
    if (!(P.minEndCov > 0)) return cut;
    var covs = [];
    for (i = 0; i < columns.length; i++) if (columns[i].pos >= 0) covs.push(columns[i].coverage);
    var n = covs.length;
    if (!n) return cut;
    var sorted = covs.slice().sort(function (x, y) { return x - y; });
    var t = Math.max(P.minEndCov, Math.ceil(CAP_END_COV_FRAC * sorted[n >> 1]));
    // A contig whose body is itself thin is left alone. At one or two reads of
    // coverage every end is carried by one read, and the rule above cut a
    // 4-read, 9.7 kb contig of the bundled Nanopore example down to the 2 kb
    // where three of its reads overlap. The threshold means something only
    // well under the body's coverage: trimming needs a median of at least
    // twice the threshold.
    if (2 * t > sorted[n >> 1]) return cut;
    var lo = 0;
    while (lo < n && covs[lo] < t) lo++;
    if (lo >= n) return cut;                 // no column reaches t: nothing to cut back to
    var hi = n - 1;
    while (hi > lo && covs[hi] < t) hi--;
    if (lo > P.maxOverhang) cut[0] = lo;
    if (n - 1 - hi > P.maxOverhang) cut[1] = n - 1 - hi;
    if (!cut[0] && !cut[1]) return cut;
    var from = cut[0], to = n - cut[1];
    // Columns are cut by position in the alignment, so the gap-only columns
    // (pos −1) inside a cut stretch go with it and those between kept bases stay.
    var ciLo = -1, ciHi = -1;
    for (i = 0; i < columns.length; i++) {
        if (columns[i].pos < 0) continue;
        if (ciLo < 0 && columns[i].pos >= from) ciLo = i;
        if (columns[i].pos < to) ciHi = i;
    }
    for (i = 0; i < columns.length; i++) {
        var col = columns[i];
        if (i < ciLo || i > ciHi) {
            col.pos = -1; col.base = ''; col.isPoly = false; col.alleles = []; col.trimmed = true;
        } else if (col.pos >= 0) {
            col.pos -= from;
        }
    }
    var keptVariants = [];
    for (i = 0; i < variants.length; i++) {
        var v = variants[i];
        if (v.pos < from || v.pos >= to) continue;
        v.pos -= from;
        keptVariants.push(v);
    }
    variants.length = 0;
    Array.prototype.push.apply(variants, keptVariants);
    consensus.splice(to, consensus.length - to);
    consensus.splice(0, from);
    return cut;
}

/* ───────────────────────────── 8. Orchestration ──────────────────────────── */

/* A run context: reports progress, hands the event loop a turn when a stage has
   been hogging it, and lets a run be cancelled. tick() yields only when enough
   time has passed, so it can be called from inside a hot loop without the
   per-call setTimeout clamp (≈4 ms once nested) dominating the runtime. */
function capContext(opts) {
    opts = opts || {};
    var onProgress = opts.onProgress || null;
    var slice = opts.sliceMs || 30;
    var last = capNow();
    return {
        cancelled: function () { return !!(opts.isCancelled && opts.isCancelled()); },
        tick: function (stage, i, total) {
            var now = capNow();
            if (now - last < slice) return Promise.resolve();
            last = now;
            if (onProgress) onProgress(stage, i, total);
            if (opts.isCancelled && opts.isCancelled()) return Promise.reject(new Error('cancelled'));
            return capYield();
        },
        stage: function (name) { if (onProgress) onProgress(name, 0, 0); }
    };
}

/* Hand the event loop one turn.
 *
 * NOT setTimeout(0): browsers clamp nested timers to ~4 ms, and a tab in the
 * background to ~1000 ms. Since a long assembly yields every ~30 ms of work,
 * that turns a 10-second run in a background tab into a run of many minutes —
 * it looks exactly like a hang, and switching away from the tab is precisely
 * what a user does while waiting. A MessageChannel message is not clamped, so
 * the yield stays cheap whether or not the tab is visible. setTimeout remains
 * as the fallback for environments without MessageChannel.
 */
function capYield() {
    if (typeof MessageChannel === 'function') {
        return new Promise(function (resolve) {
            var ch = new MessageChannel();
            ch.port1.onmessage = function () { ch.port1.close(); resolve(); };
            ch.port2.postMessage(0);
        });
    }
    return new Promise(function (resolve) {
        if (typeof setTimeout === 'function') setTimeout(resolve, 0);
        else resolve();
    });
}

/**
 * Run the whole pipeline.
 *
 * Everything except the overlap stage runs here, chunked so a long run never
 * freezes the page. The overlap stage — the expensive one — is delegated to
 * `opts.overlapRunner`, which is either the same-thread chunked runner (default,
 * and what Node and a file:// page use) or the Worker-pool runner installed by
 * js/cap-pool.js. Both produce identical results; only the wall-clock differs.
 *
 * @param {object} params  {input, quality, …CAP_DEFAULTS}
 * @param {object} [opts]  {onProgress(stage,i,total), overlapRunner, isCancelled, threads}
 * @returns {Promise<object>} the assembly result
 */
async function capAssemble(params, opts) {
    opts = opts || {};
    var ctx = capContext(opts);
    var overlapRunner = opts.overlapRunner || capSerialOverlapRunner;
    var t0 = capNow();
    var P = capParams(params);
    var parsed = capParseReads(params.input || '', params.quality || '', P);
    var reads = parsed.reads;
    var notes = parsed.notes.slice();

    // Fold-back reads are cut back to one copy of their sequence before anything
    // is seeded, so the index and the overlaps never see the second copy.
    var foldbacks = capFoldbacks(reads, P);
    if (foldbacks.length) {
        var fDropped = foldbacks.filter(function (f) { return f.action === 'dropped'; }).length;
        notes.push(foldbacks.length + ' read(s) end in the reverse complement of their own sequence (fold-back reads: ' +
                   'the strand was read on into its complement). ' +
                   (foldbacks.length - fDropped) + ' were cut back to one copy' +
                   (fDropped ? ' and ' + fDropped + ' dropped because what remained was shorter than minReadLen' : '') +
                   '. They are listed under Fold-back reads.');
    }

    if (reads.length < 1) {
        return capEmptyResult(P, notes, 'No usable sequences were found in the input.');
    }
    if (reads.length < 2) {
        // Not an error — the single fragment is reported as a singleton — but say
        // so plainly, because "0 contigs" on its own reads like a failure.
        notes.push('Only one usable fragment was found, so there is nothing to assemble it with. ' +
                   'It is reported below as a singleton. Assembly needs at least two overlapping fragments.');
    }
    // Mean quality per read, used to weight consensus votes.
    for (var i = 0; i < reads.length; i++) {
        if (reads[i].qual && reads[i].qual.length) {
            var s = 0;
            for (var q = 0; q < reads[i].qual.length; q++) s += reads[i].qual[q];
            reads[i].meanQual = s / reads[i].qual.length;
        } else {
            reads[i].meanQual = 1;
        }
    }

    ctx.stage('seeding');
    var cand = await capCandidatePairs(reads, P, ctx);
    if (cand.stats.truncated) {
        notes.push('Candidate-pair list hit the maxCandidates limit (' + P.maxCandidates +
                   '); some overlaps were not evaluated. Raise maxCandidates or the k-mer length.');
    }
    if (cand.stats.skippedWords) {
        notes.push(cand.stats.skippedWords + ' k-mer(s) occurred in more than ' + cand.stats.deepCut +
                   ' places and were treated as repeat words (not used for seeding).');
    }

    ctx.stage('overlaps');
    var ovRes = await overlapRunner(reads, cand.pairs, P, ctx, cand.stats.repeatThreshold);

    // A chimeric read takes no part in the layout: all its overlaps are dropped,
    // and it is reported on its own rather than as a singleton (which would put
    // the joined sequence into the FASTA as if it were a fragment of the target).
    var chimeras = P.chimeraCheck ? capChimeras(reads, ovRes.overlaps, P) : [];
    if (P.chimeraCheck) {
        var flagged = new Set(chimeras.map(function (c) { return c.read; }));
        chimeras.forEach(function (c) { c.evidence = 'piles'; });
        capPartingChimeras(reads, ovRes.overlaps, ovRes.rejected, P).forEach(function (c) {
            if (!flagged.has(c.read)) { flagged.add(c.read); chimeras.push(c); }
        });
        chimeras.sort(function (x, y) { return x.read - y.read; });
    }
    if (chimeras.length) {
        var isChimera = new Set(chimeras.map(function (c) { return c.read; }));
        ovRes.overlaps = ovRes.overlaps.filter(function (o) { return !isChimera.has(o.a) && !isChimera.has(o.b); });
        notes.push(chimeras.length + ' read(s) carry two molecules joined in one read (a ligation chimera, a ' +
                   'fold-back read, or two molecules read back to back) and were set aside: the reads that ' +
                   'overlap one part stop at the junction or part from the read there, those that overlap ' +
                   'the other part start there, and almost none run across it. They are listed under ' +
                   'Chimeric reads.');
    }

    // After the chimeras are out, so a chimera's own continuation cannot vouch
    // for a repeat.
    if (P.repeatGuard) {
        var provenRefused = capProvenRepeats(reads, ovRes.overlaps, ovRes.rejected, P);
        Array.prototype.push.apply(ovRes.rejected, provenRefused);
    }

    ctx.stage('layout');
    var lay = await capLayout(reads, ovRes.overlaps, P, ctx);
    if (chimeras.length) {
        lay.singletons = lay.singletons.filter(function (ri) { return !isChimera.has(ri); });
    }
    lay.chimeras = chimeras;
    lay.foldbacks = foldbacks;

    // Measure the sequencing error rate from the assembly's own overlaps instead
    // of assuming one. A pairwise alignment shows a mismatch when EITHER read is
    // wrong, so the per-read rate is about half the pairwise rate.
    var mmTot = 0, colTot = 0;
    for (var oi2 = 0; oi2 < ovRes.overlaps.length; oi2++) {
        mmTot += ovRes.overlaps[oi2].mismatches;
        colTot += ovRes.overlaps[oi2].matches + ovRes.overlaps[oi2].mismatches;
    }
    var errorRate = colTot ? Math.min(0.25, Math.max(0.0005, (mmTot / colTot) / 2)) : 0.01;

    ctx.stage('consensus');
    // Align every contig first: the significance threshold for a polymorphic
    // column is a correction over the WHOLE assembly, so it cannot be fixed until
    // the total number of columns is known.
    var aligns = [];
    for (var ca2 = 0; ca2 < lay.contigs.length; ca2++) {
        aligns.push(await capContigAlign(reads, lay.contigs[ca2], P, ctx));
    }
    var totalCols = 0;
    aligns.forEach(function (a) { totalCols += a.cols.length; });
    var alpha = P.polyMaxP > 0 ? P.polyMaxP : 1 / Math.max(1, totalCols);
    var model = { errorRate: errorRate, alpha: alpha };

    var contigs = [];
    for (var c = 0; c < lay.contigs.length; c++) {
        var aligned = aligns[c];
        var cons = capContigConsensus(reads, aligned, P, model);
        contigs.push({
            id: lay.contigs[c].id,
            reads: lay.contigs[c].reads,
            rows: aligned.rows,
            cols: aligned.cols,
            consensus: cons.consensus,
            columns: cons.columns,
            variants: cons.variants,
            trimmed: cons.trimmed,
            meanCoverage: cons.meanCoverage,
            split: lay.contigs[c].split
        });
        await ctx.tick('consensus', c + 1, lay.contigs.length);
    }
    contigs.sort(function (p, q) { return q.consensus.length - p.consensus.length || p.id - q.id; });
    contigs.forEach(function (ct, k) { ct.id = k + 1; });

    return capFinish(reads, contigs, lay, ovRes, cand, P, notes, capNow() - t0,
                     opts.threads || 1, model);
}

// Assemble the final result object from already-computed stages. Shared by the
// serial path and the Worker path so both report identically.
function capFinish(reads, contigs, lay, ovRes, cand, P, notes, elapsed, threads, model) {
    var lens = contigs.map(function (c) { return c.consensus.length; }).sort(function (a, b) { return b - a; });
    var totalBp = lens.reduce(function (a, b) { return a + b; }, 0);
    var n50 = 0, acc = 0;
    for (var i = 0; i < lens.length; i++) { acc += lens[i]; if (acc >= totalBp / 2) { n50 = lens[i]; break; } }
    var placed = contigs.reduce(function (a, c) { return a + c.reads.length; }, 0);
    var variants = 0, snps = 0, indels = 0, ties = 0, trimmedBp = 0;
    contigs.forEach(function (c) {
        variants += c.variants.length;
        c.variants.forEach(function (v) {
            if (v.kind === 'snp') snps++;
            else if (v.kind === 'indel') indels++;
            else if (v.kind === 'tie') ties++;
        });
        if (c.trimmed) trimmedBp += c.trimmed[0] + c.trimmed[1];
    });
    var foldbacks = lay.foldbacks || [];
    var foldbacksDropped = foldbacks.filter(function (f) { return f.action === 'dropped'; }).length;

    return {
        ok: true,
        params: P,
        reads: reads,
        contigs: contigs,
        singletons: lay.singletons,
        chimeras: lay.chimeras || [],
        foldbacks: foldbacks,
        overlaps: ovRes.overlaps,
        rejected: ovRes.rejected,
        conflicts: lay.conflicts,
        candidates: cand.pairs.length,
        seedStats: cand.stats,
        notes: notes,
        stats: {
            readsIn: reads.length,
            readsPlaced: placed,
            singletons: lay.singletons.length,
            chimeras: (lay.chimeras || []).length,
            foldbacks: foldbacks.length - foldbacksDropped,
            foldbacksDropped: foldbacksDropped,
            trimmedBp: trimmedBp,
            contigs: contigs.length,
            totalBp: totalBp,
            n50: n50,
            longest: lens.length ? lens[0] : 0,
            meanCoverage: contigs.length
                ? contigs.reduce(function (a, c) { return a + c.meanCoverage; }, 0) / contigs.length : 0,
            overlaps: ovRes.overlaps.length,
            rejected: ovRes.rejected.length,
            conflicts: lay.conflicts.length,
            variants: variants, snps: snps, indels: indels, ties: ties,
            errorRate: model ? model.errorRate : 0,
            alpha: model ? model.alpha : 0,
            junctionIdentity: capJunctionIdentity(P),
            elapsedMs: elapsed,
            threads: threads
        }
    };
}

function capEmptyResult(P, notes, message) {
    return {
        ok: false, message: message, params: P, reads: [], contigs: [], singletons: [], chimeras: [], foldbacks: [],
        overlaps: [], rejected: [], conflicts: [], candidates: 0,
        seedStats: { words: 0, seeds: 0, skippedWords: 0, truncated: false },
        notes: notes,
        stats: { readsIn: 0, readsPlaced: 0, singletons: 0, chimeras: 0, foldbacks: 0, foldbacksDropped: 0,
                 trimmedBp: 0, contigs: 0, totalBp: 0, n50: 0,
                 longest: 0, meanCoverage: 0, overlaps: 0, rejected: 0, conflicts: 0,
                 variants: 0, snps: 0, indels: 0, ties: 0, junctionIdentity: capJunctionIdentity(P),
                 elapsedMs: 0, threads: 1 }
    };
}

function capNow() {
    if (typeof performance !== 'undefined' && performance.now) return performance.now();
    return Date.now();
}

/* ───────────────────────────── 9. Exports ────────────────────────────────── */

// Same guard style as js/repeats.js and js/muscle.js: the file is a plain script
// that works as a browser global, inside a Worker (self), and under Node.
(function (g) {
    var api = {
        assemble: capAssemble,
        params: capParams,
        parseReads: capParseReads,
        parseRecords: capParseFasta,                   // FASTA / FASTQ / tab / plain -> {name, desc, body}
        isFastq: function (text) { return capIsFastq(String(text == null ? '' : text).split(/\r\n?|\n/)); },
        candidatePairs: capCandidatePairs,
        computeOverlaps: capComputeOverlaps,           // synchronous; the Worker entry point
        serialOverlapRunner: capSerialOverlapRunner,
        overlapPair: capOverlapPair,
        foldbacks: capFoldbacks, partingChimeras: capPartingChimeras,
        junctionIdentity: capJunctionIdentity,
        layout: capLayout,
        contigAlign: capContigAlign,
        contigConsensus: capContigConsensus,
        align: capAlignBanded,
        alignAuto: capAlignAuto,
        rc: capRC,
        finish: capFinish,
        context: capContext,
        DEFAULTS: CAP_DEFAULTS,
        IUPAC_FROM_BITS: CAP_FROM_BITS,
        IUPAC_BITS: CAP_BITS
    };
    if (typeof module !== 'undefined' && module.exports) module.exports = api;
    if (g) g.CapJS = api;
})(typeof self !== 'undefined' ? self : (typeof window !== 'undefined' ? window : (typeof globalThis !== 'undefined' ? globalThis : null)));
