'use strict';
/* ────────────────────────────────────────────────────────────────────────────
   cap-report.js — turns an OLC Assembler result into the four artefacts the
   page shows and lets the user download:

     capReportText(result)            the assembly report  (human-readable + TSV tables)
     capReportFasta(result, opts)     consensus sequences  (FASTA)
     capReportPlacements(result)      read placements      (TSV, for downstream use)
     capReportAlignmentText(result)   the read-under-consensus alignment (plain text)
     capReportAlignmentHtml(result)   the same view, coloured, for the browser

   Pure string building — no DOM, no engine internals beyond the documented
   result shape — so it can be unit-tested under Node together with the engine.

   The alignment view follows the convention msa.html established: plain text
   plus inline colour <span>s inside a <pre>, so click-drag selection stays
   inside the box and copies as clean text.
   ──────────────────────────────────────────────────────────────────────────── */
(function (global) {

    function pad(s, n) { s = String(s); return s.length >= n ? s : s + new Array(n - s.length + 1).join(' '); }
    function wrap(s, n) {
        var out = [];
        for (var i = 0; i < s.length; i += n) out.push(s.substring(i, i + n));
        return out.join('\n');
    }
    function esc(s) {
        return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
    }
    function pct(x) { return (isFinite(x) ? x : 0).toFixed(1); }

    /* Where each read sits IN CONSENSUS COORDINATES.
     *
     * contig.reads[].pos is the greedy layout's estimate, made before the reads
     * were stacked into columns. The column alignment then inserts gap columns,
     * so layout positions and consensus positions drift apart along a contig —
     * printing the layout number next to a consensus-based figure (a variant
     * position, the contig length) quietly invites the reader to compare two
     * different coordinate systems. contig.rows carries each read's true first and
     * last COLUMN, and contig.columns maps a column to its consensus position, so
     * the real answer is available; this just looks it up.
     *
     * Returns a Map readIdx -> {start, end} as 1-based inclusive consensus
     * coordinates, falling back to the layout estimate if a read covers no
     * consensus-bearing column at all.
     */
    function consensusSpans(c) {
        var posOfCol = {};
        var i;
        for (i = 0; i < c.columns.length; i++) {
            if (c.columns[i].pos >= 0) posOfCol[c.columns[i].col] = c.columns[i].pos;
        }
        // For a column that carries no consensus base (an insertion seen in a
        // minority of reads), walk outwards to the nearest one that does.
        var nearest = function (col, dir) {
            for (var k = col; k >= 0 && k < c.cols.length; k += dir) {
                if (posOfCol[k] !== undefined) return posOfCol[k];
            }
            return null;
        };
        var out = new Map();
        for (i = 0; i < c.rows.length; i++) {
            var rw = c.rows[i];
            var s = nearest(rw.first, 1);
            var e = nearest(rw.last, -1);
            if (s === null || e === null) out.set(rw.readIdx, null);
            else out.set(rw.readIdx, { start: s + 1, end: e + 1 });
        }
        return out;
    }

    // The same total order capLayout() consumes overlaps in, so every table in the
    // report is reproducible regardless of how the work was distributed.
    function capOverlapOrder(p, q) {
        if (q.score !== p.score) return q.score - p.score;
        if (q.olen !== p.olen) return q.olen - p.olen;
        if (q.identity !== p.identity) return q.identity - p.identity;
        if (p.a !== q.a) return p.a - q.a;
        if (p.b !== q.b) return p.b - q.b;
        return p.flip - q.flip;
    }
    // A rejected overlap may have been refused before it was ever scored, so it is
    // ordered by the read pair itself, which is always present.
    function capRejectOrder(p, q) {
        if (p.a !== q.a) return p.a - q.a;
        if (p.b !== q.b) return p.b - q.b;
        return (p.flip || 0) - (q.flip || 0);
    }

    /* ───────────────────────────── assembly report ───────────────────────── */

    function capReportText(res) {
        if (!res || !res.ok) {
            return 'OLC Assembler\n=============\n\n' +
                   ((res && res.message) || 'The assembly produced no result.') + '\n' +
                   ((res && res.notes && res.notes.length) ? '\nNotes\n-----\n- ' + res.notes.join('\n- ') + '\n' : '');
        }
        var S = res.stats, P = res.params, L = [];
        var reads = res.reads;

        L.push('OLC Assembler report');
        L.push('====================');
        L.push('');
        L.push('Run summary');
        L.push('-----------');
        var row = function (k, v) { L.push(pad(k, 32) + v); };
        row('Reads used', S.readsIn);
        row('Reads placed in contigs', S.readsPlaced);
        row('Singletons (unassembled)', S.singletons);
        if (S.chimeras) row('Chimeric reads (set aside)', S.chimeras);
        if (S.foldbacks || S.foldbacksDropped) {
            row('Fold-back reads trimmed', S.foldbacks + (S.foldbacksDropped ? ' (and ' + S.foldbacksDropped + ' dropped)' : ''));
        }
        row('Contigs', S.contigs);
        row('Total consensus', S.totalBp + ' bp');
        row('Longest contig', S.longest + ' bp');
        row('N50', S.n50 + ' bp');
        row('Mean coverage', pct(S.meanCoverage) + '×');
        row('Candidate pairs aligned', res.candidates);
        row('Overlaps accepted', S.overlaps);
        row('Overlaps rejected', S.rejected);
        row('Layout conflicts', S.conflicts);
        row('Polymorphic sites (SNP)', S.snps);
        row('Indel polymorphisms', S.indels);
        if (S.ties) row('Unresolved ties (IUPAC)', S.ties);
        if (S.trimmedBp) row('Contig ends trimmed', S.trimmedBp + ' bp');
        row('Elapsed', (S.elapsedMs / 1000).toFixed(2) + ' s' +
                       (S.threads > 1 ? ' on ' + S.threads + ' worker threads' : ' on 1 thread'));
        L.push('');
        L.push('Parameters');
        L.push('----------');
        L.push('k-mer ' + P.kmer + ' · min seeds ' + P.minSeeds + ' · min overlap ' + P.minOverlap +
               ' nt · min identity ' + P.minIdentity + '% · max overhang ' + P.maxOverhang + ' nt');
        L.push('scoring: match ' + P.matchScore + ', mismatch ' + P.mismatchScore +
               ', gap open ' + P.gapOpen + ', gap extend ' + P.gapExtend + ' · band ±' + P.bandWidth +
               ' (max ±' + P.maxBandWidth + ')');
        if (P.consGapOpen !== undefined) {
            L.push('consensus: reads fitted with gap open ' + P.consGapOpen + ' · homopolymer runs ' +
                   (P.runVoteMin > 0 ? 'of ≥ ' + P.runVoteMin + ' nt checked against the run lengths in the reads' : 'not checked against the reads'));
        }
        L.push('polymorphism: minor allele ≥ ' + P.polyMinFrac + '% and ≥ ' + P.polyMinCount +
               ' read(s), column coverage ≥ ' + P.polyMinCov +
               (P.qualWeight ? ' · votes weighted by PHRED quality when supplied' : ''));
        L.push('reverse-complement search: ' + (P.allowRC ? 'on' : 'off') +
               ' · chimeric reads: ' + (P.chimeraCheck ? 'set aside' : 'kept') +
               (P.chimeraCheck ? ' (fold-back copies ' + (P.foldbackMin > 0 ? '≥ ' + P.foldbackMin + ' nt trimmed' : 'off') +
                                 (S.junctionIdentity ? ', junctions below ' + S.junctionIdentity + '% local identity refused' : '') + ')' : '') +
               ' · contig ends need ' + (P.minEndCov > 0 ? '≥ ' + P.minEndCov + ' reads (or 5% of the median coverage) in contigs at least twice that deep' : 'no minimum coverage') +
               ' · min read length ' + P.minReadLen + ' nt' +
               (P.endTrim ? ' · end trim ' + P.endTrim + ' nt' : '') +
               (P.qualClip ? ' · quality clip < Q' + P.qualClip : ''));
        L.push('');

        L.push('Contigs');
        L.push('-------');
        L.push('Trimmed5 and Trimmed3 are consensus bases cut from the two ends because too few reads carried them.');
        L.push(['Contig', 'Length(bp)', 'Reads', 'MeanCoverage', 'MinCoverage', 'SNP', 'Indel', 'Forward', 'Reverse',
                'Trimmed5', 'Trimmed3'].join('\t'));
        res.contigs.forEach(function (c) {
            var minCov = Infinity, fwd = 0, rev = 0;
            c.columns.forEach(function (col) { if (col.pos >= 0 && col.coverage < minCov) minCov = col.coverage; });
            c.reads.forEach(function (p) { if (p.flip) rev++; else fwd++; });
            L.push(['Contig' + c.id, c.consensus.length, c.reads.length, pct(c.meanCoverage),
                    isFinite(minCov) ? minCov : 0,
                    c.variants.filter(function (v) { return v.kind === 'snp'; }).length,
                    c.variants.filter(function (v) { return v.kind === 'indel'; }).length,
                    fwd, rev, c.trimmed ? c.trimmed[0] : 0, c.trimmed ? c.trimmed[1] : 0].join('\t'));
        });
        L.push('');

        L.push('Read placement');
        L.push('--------------');
        L.push('Start and End are positions in that contig\'s consensus, 1-based inclusive.');
        L.push(['Contig', 'Read', 'Strand', 'Start', 'End', 'Length(nt)', 'Clipped5', 'Clipped3'].join('\t'));
        res.contigs.forEach(function (c) {
            var spans = consensusSpans(c);
            c.reads.forEach(function (p) {
                var r = reads[p.readIdx];
                var sp = spans.get(p.readIdx);
                L.push(['Contig' + c.id, r.id, p.flip ? '-' : '+',
                        sp ? sp.start : (p.pos + 1), sp ? sp.end : (p.pos + p.len),
                        r.len, r.clipLeft, r.clipRight].join('\t'));
            });
        });
        L.push('');

        L.push('Polymorphic sites');
        L.push('-----------------');
        if (!S.variants) {
            L.push('None. No column had a minority base supported by ≥ ' + P.polyMinCount +
                   ' reads and ≥ ' + P.polyMinFrac + '% of its coverage.');
        } else {
            if (S.ties) {
                L.push('Kind "tie": the reads are split evenly and no base lies further inside its read ' +
                       'than the others, so the consensus carries the IUPAC code of the tied bases. ' +
                       'Not a polymorphism call; more reads, or quality values, resolve it.');
            }
            L.push(['Contig', 'Position', 'Kind', 'Coverage', 'A', 'C', 'G', 'T', 'Gap',
                    'IUPAC', 'MinorFraction(%)', 'Note'].join('\t'));
            res.contigs.forEach(function (c) {
                c.variants.forEach(function (v) {
                    var n = v.counts, tot = n.a + n.c + n.g + n.t;
                    var minor;
                    if (v.kind === 'indel') {
                        // The minority allele of an indel is the gap (for a
                        // deletion) or the inserted base (for an insertion), and
                        // neither is in v.counts — that only tallies A/C/G/T. The
                        // SNP formula therefore reports 0.0 for every indel, which
                        // sorts them all to the bottom of the table exactly when
                        // someone is triaging by how well supported a variant is.
                        minor = v.coverage ? 100 * (v.code ? v.gap : tot) / v.coverage : 0;
                    } else {
                        var srt = [n.a, n.c, n.g, n.t].sort(function (x, y) { return y - x; });
                        minor = tot ? 100 * srt[1] / tot : 0;
                    }
                    L.push(['Contig' + c.id, v.pos + 1, v.kind, v.coverage,
                            n.a, n.c, n.g, n.t, v.gap,
                            (v.code || '-').toUpperCase(), pct(minor), v.note || ''].join('\t'));
                });
            });
        }
        L.push('');

        L.push('Overlaps that passed the filters');
        L.push('--------------------------------');
        L.push('Role: "joining" built a contig, "confirming" agreed with a join already made, ' +
               '"conflict" passed every filter but contradicted the layout and was not used ' +
               '(those are listed with their reasons under Layout conflicts).');
        L.push(['ReadA', 'ReadB', 'Strand', 'Offset', 'Overlap(nt)', 'Identity(%)', 'Score',
                'Mismatches', 'Gaps', 'Contained', 'Role'].join('\t'));
        // Sorted by a TOTAL order, the same one the layout uses. With the Worker
        // pool the overlaps arrive in whatever order the shards finish, so sorting
        // on score alone would leave ties in a run-dependent order and two
        // identical runs would produce different report files.
        res.overlaps.slice().sort(capOverlapOrder).forEach(function (o) {
            L.push([reads[o.a].id, reads[o.b].id, o.flip ? '-' : '+', o.offset, o.olen,
                    pct(o.identity), o.score, o.mismatches, o.gaps,
                    o.contained || '', o.role || 'joining'].join('\t'));
        });
        L.push('');

        if (res.rejected.length) {
            L.push('Rejected overlaps');
            L.push('-----------------');
            L.push('Overlaps that a k-mer seed suggested but that failed a filter. A large ' +
                   '"unmatched overhang" count is the normal signature of an interspersed repeat.');
            L.push(['ReadA', 'ReadB', 'Strand', 'Overlap(nt)', 'Identity(%)', 'Reason'].join('\t'));
            res.rejected.slice().sort(capRejectOrder).forEach(function (o) {
                L.push([reads[o.a].id, reads[o.b].id, o.flip ? '-' : '+',
                        o.olen === undefined ? '' : o.olen,
                        o.identity === undefined ? '' : pct(o.identity), o.rejected].join('\t'));
            });
            L.push('');
        }

        if (res.conflicts.length) {
            L.push('Layout conflicts');
            L.push('----------------');
            L.push('These overlaps passed every filter but contradicted the layout already built ' +
                   'from stronger overlaps. Each one is a place where the data is ambiguous — most ' +
                   'often two copies of a repeat that the reads cannot tell apart.');
            L.push(['ReadA', 'ReadB', 'Strand', 'Overlap(nt)', 'Identity(%)', 'Reason'].join('\t'));
            res.conflicts.slice().sort(function (x, y) { return capRejectOrder(x.ov, y.ov); }).forEach(function (cf) {
                L.push([reads[cf.ov.a].id, reads[cf.ov.b].id, cf.ov.flip ? '-' : '+',
                        cf.ov.olen, pct(cf.ov.identity), cf.reason].join('\t'));
            });
            L.push('');
        }

        if (res.singletons.length) {
            L.push('Singletons');
            L.push('----------');
            L.push('Reads with no acceptable overlap to any other read.');
            L.push(['Read', 'Length(nt)', 'CG(%)', 'RepeatWordHits'].join('\t'));
            res.singletons.forEach(function (ri) {
                var r = reads[ri];
                L.push([r.id, r.len, pct(cg(r.seq)), r.repeatiness].join('\t'));
            });
            L.push('');
        }

        if (res.chimeras && res.chimeras.length) {
            L.push('Chimeric reads');
            L.push('--------------');
            L.push('Reads that carry two molecules: the reads overlapping one part stop at the junction, ' +
                   'those overlapping the other start there, and almost none run across it. They were ' +
                   'left out of the layout and out of the FASTA. Junction = position on the read as given. ' +
                   'Evidence "piles": the overlapping reads end and begin there. Evidence "parting": the ' +
                   'partners agree with the read up to (or from) there and then go on with sequence of ' +
                   'their own, and nothing vouches for the read across it; ReadsEnding and ReadsStarting ' +
                   'then count those partners.');
            L.push(['Read', 'Length(nt)', 'Junction(nt)', 'ReadsEnding', 'ReadsStarting', 'ReadsAcross', 'Evidence'].join('\t'));
            res.chimeras.forEach(function (c) {
                var r = reads[c.read];
                L.push([r.id, r.len, c.to > c.at ? c.at + '-' + c.to : c.at, c.ending, c.starting, c.spanning,
                        c.evidence || 'piles'].join('\t'));
            });
            L.push('');
        }

        if (res.foldbacks && res.foldbacks.length) {
            L.push('Fold-back reads');
            L.push('---------------');
            L.push('Reads that carry a stretch of their own sequence twice, the second time reverse-complemented: ' +
                   'the strand was read on into its complement. Each had the duplicated stretch cut away before ' +
                   'assembly (Kept, on the read as given; the rest counts as clipped in the placement table), or ' +
                   'was dropped when nothing long enough remained. A hairpin keeps its longer arm; copies that lie ' +
                   'apart keep only what lies between them, since either copy may be the one the strand re-read. ' +
                   'Copy1 and Copy2 are the two copies on the read, 1-based.');
            L.push(['Read', 'Length(nt)', 'Copy1', 'Copy2', 'Identity(%)', 'Kept', 'Action'].join('\t'));
            res.foldbacks.forEach(function (f) {
                L.push([f.id, f.len, (f.copies[0] + 1) + '-' + f.copies[1], (f.copies[2] + 1) + '-' + f.copies[3],
                        pct(f.identity), f.kept ? (f.kept[0] + 1) + '-' + f.kept[1] : '-', f.action].join('\t'));
            });
            L.push('');
        }

        var repetitive = reads.filter(function (r) { return r.repeatiness > 0; });
        if (repetitive.length) {
            L.push('Reads containing repeat words');
            L.push('-----------------------------');
            L.push('These reads contain k-mers that occur in more than ' + (res.seedStats.deepCut || P.maxSeedOcc) +
                   ' places. Joins involving them deserve a second look.');
            L.push(['Read', 'RepeatWordHits'].join('\t'));
            repetitive.sort(function (x, y) { return y.repeatiness - x.repeatiness; })
                      .slice(0, 50)
                      .forEach(function (r) { L.push([r.id, r.repeatiness].join('\t')); });
            L.push('');
        }

        if (res.notes.length) {
            L.push('Notes');
            L.push('-----');
            res.notes.forEach(function (n) { L.push('- ' + n); });
            L.push('');
        }
        return L.join('\n');
    }

    function cg(s) {
        var n = 0, tot = 0;
        for (var i = 0; i < s.length; i++) {
            var c = s.charAt(i);
            if (c === 'g' || c === 'c' || c === 's') n++;
            if ('acgt'.indexOf(c) !== -1 || c === 's' || c === 'w') tot++;
        }
        return tot ? 100 * n / tot : 0;
    }

    /* ───────────────────────────── consensus FASTA ───────────────────────── */

    function capReportFasta(res, opts) {
        opts = opts || {};
        var width = opts.width || 60;
        var upper = opts.upper !== false;      // IUPAC codes read better upper-case
        var out = [];
        if (!res || !res.ok) return '';

        res.contigs.forEach(function (c) {
            var snps = c.variants.filter(function (v) { return v.kind === 'snp'; }).length;
            var indels = c.variants.filter(function (v) { return v.kind === 'indel'; }).length;
            var ties = c.variants.filter(function (v) { return v.kind === 'tie'; }).length;
            out.push('>Contig' + c.id + ' length=' + c.consensus.length +
                     ' reads=' + c.reads.length +
                     ' cov=' + pct(c.meanCoverage) + 'x' +
                     ' snp=' + snps + ' indel=' + indels + (ties ? ' tie=' + ties : ''));
            out.push(wrap(upper ? c.consensus.toUpperCase() : c.consensus, width));
        });

        if (opts.includeSingletons !== false && res.singletons.length) {
            res.singletons.forEach(function (ri) {
                var r = res.reads[ri];
                out.push('>' + r.id + ' length=' + r.len + ' reads=1 cov=1.0x unassembled');
                out.push(wrap(upper ? r.seq.toUpperCase() : r.seq, width));
            });
        }
        return out.length ? out.join('\n') + '\n' : '';
    }

    /* ───────────────────────────── placements TSV ────────────────────────── */

    function capReportPlacements(res) {
        if (!res || !res.ok) return '';
        var L = [['contig', 'contig_length', 'read', 'strand', 'contig_start', 'contig_end',
                  'read_length', 'clipped_5prime', 'clipped_3prime'].join('\t')];
        res.contigs.forEach(function (c) {
            var spans = consensusSpans(c);
            c.reads.forEach(function (p) {
                var r = res.reads[p.readIdx];
                var sp = spans.get(p.readIdx);
                L.push(['Contig' + c.id, c.consensus.length, r.id, p.flip ? '-' : '+',
                        sp ? sp.start : (p.pos + 1), sp ? sp.end : (p.pos + p.len),
                        r.len, r.clipLeft, r.clipRight].join('\t'));
            });
        });
        res.singletons.forEach(function (ri) {
            var r = res.reads[ri];
            L.push(['(singleton)', r.len, r.id, '+', 1, r.len, r.len, r.clipLeft, r.clipRight].join('\t'));
        });
        return L.join('\n') + '\n';
    }

    /* ───────────────────────────── alignment view ────────────────────────── */

    /* One block per `width` columns: the consensus on top, then every read that
       covers the block. A base identical to the consensus is shown as a dot, so
       the eye lands straight on the disagreements — which, in an assembly, is the
       whole point of looking at the alignment. */

    // Block width comes straight from a number input, so it can be 0, negative,
    // fractional or absurd. A width of 0 makes the block loop below step by zero
    // and never terminate — the page would hang with no way back — so it is
    // clamped here rather than trusted.
    function blockWidth(w) {
        var n = Math.floor(Number(w));
        if (!isFinite(n) || n < 20) n = 20;
        if (n > 1000) n = 1000;
        return n;
    }

    /* Split the contigs into displayable blocks, within a column budget.
     *
     * A contig that does not fit is CLIPPED, not dropped: dropping it whole meant
     * that a single contig longer than the budget produced an empty view — and
     * the empty view was then reported as "No contigs were assembled", which is
     * the opposite of what happened. Showing the first N columns and saying so is
     * both true and useful. */
    function buildBlocks(res, width, maxCols) {
        var blocks = [];
        var used = 0, truncated = false;

        for (var ci = 0; ci < res.contigs.length; ci++) {
            if (used >= maxCols) { truncated = true; break; }
            var c = res.contigs[ci];
            var nCols = c.cols.length;
            if (used + nCols > maxCols) { nCols = maxCols - used; truncated = true; }
            used += nCols;

            // consensus letter per COLUMN (gap columns show '-')
            var consLine = new Array(nCols), polyAt = new Array(nCols);
            for (var k = 0; k < nCols; k++) {
                var col = c.columns[k];
                consLine[k] = (col && col.base) ? col.base : '-';
                polyAt[k] = !!(col && col.isPoly);
            }

            var rows = c.rows.map(function (rw) {
                var r = res.reads[rw.readIdx];
                return { id: r.id, strand: rw.flip ? '-' : '+', first: rw.first, last: rw.last,
                         cells: c.cols, readIdx: rw.readIdx };
            });

            blocks.push({ contig: c, consLine: consLine, polyAt: polyAt, rows: rows,
                          nCols: nCols, clipped: nCols < c.cols.length });
        }
        return { blocks: blocks, truncated: truncated };
    }

    function truncationNotice(built, res) {
        if (!built.truncated) return '';
        var shown = built.blocks.length;
        return '[ Alignment view truncated: showing ' + shown + ' of ' + res.contigs.length +
               ' contig(s)' + (built.blocks.length && built.blocks[built.blocks.length - 1].clipped
                   ? ', the last one only partly' : '') +
               '. The consensus FASTA and the assembly report cover every contig in full. ]';
    }

    function rowCharAt(row, colIdx, consChar) {
        if (colIdx < row.first || colIdx > row.last) return ' ';
        var cell = row.cells[colIdx].get(row.readIdx);
        if (cell === undefined) return ' ';
        var ch = String.fromCharCode(cell & 127);
        if (ch === '-') return '-';
        return ch === consChar ? '.' : ch;
    }

    function capReportAlignmentText(res, opts) {
        opts = opts || {};
        var width = blockWidth(opts.width || 120);
        var built = buildBlocks(res, width, opts.maxCols || 400000);
        if (!res || !res.ok) return '';
        var L = [];

        built.blocks.forEach(function (B) {
            var c = B.contig;
            var idw = Math.max(12, Math.min(28, Math.max.apply(null, B.rows.map(function (r) { return r.id.length; }).concat([12]))));
            L.push('Contig' + c.id + '  ' + c.consensus.length + ' bp  ·  ' + c.reads.length +
                   ' reads  ·  mean coverage ' + pct(c.meanCoverage) + '×  ·  ' +
                   B.nCols + ' alignment columns');
            L.push(new Array(70).join('─'));
            for (var start = 0; start < B.nCols; start += width) {
                var end = Math.min(start + width, B.nCols);
                L.push('');
                L.push('columns ' + (start + 1) + '–' + end);
                var cline = '';
                for (var k = start; k < end; k++) cline += B.consLine[k];
                L.push(pad('consensus', idw) + '   ' + cline);
                B.rows.forEach(function (row) {
                    if (row.last < start || row.first >= end) return;   // read not in this block
                    var s = '';
                    for (var q = start; q < end; q++) s += rowCharAt(row, q, B.consLine[q]);
                    L.push(pad(row.id.substring(0, idw), idw) + ' ' + row.strand + ' ' + s);
                });
                // marker line: * under polymorphic (IUPAC) columns
                var mline = '';
                for (var m = start; m < end; m++) mline += B.polyAt[m] ? '*' : ' ';
                if (mline.replace(/ /g, '')) L.push(pad('', idw) + '   ' + mline);
            }
            L.push('');
        });
        if (built.truncated) L.push(truncationNotice(built, res));
        return L.join('\n');
    }

    /* Colouring one <span> per base would emit ~150 bytes per column: a 20 kb
       assembly at 5× coverage would be a 15 MB DOM. Runs of identically coloured
       characters — and in an assembly view most of a row is a run of dots — are
       emitted as a single span instead. */
    function runPainter(out) {
        var openStyle = null, buf = '';
        return {
            put: function (ch, style) {
                if (style !== openStyle) { this.flush(); openStyle = style; }
                buf += ch;
            },
            flush: function () {
                if (!buf) { openStyle = null; return; }
                if (openStyle) out.push('<span style="' + openStyle + '">' + buf + '</span>');
                else out.push(buf);
                buf = ''; openStyle = null;
            }
        };
    }

    function capReportAlignmentHtml(res, opts) {
        opts = opts || {};
        var width = blockWidth(opts.width || 120);
        if (!res || !res.ok) return esc((res && res.message) || 'No assembly to show.');
        // A tighter budget than the plain-text view: markup costs bytes, and a
        // browser asked to lay out millions of coloured cells stops being useful
        // long before it runs out of memory.
        var built = buildBlocks(res, width, opts.maxCols || 120000);
        if (!built.blocks.length) {
            // Distinguish 'there is nothing to show' from 'it did not fit'.
            return built.truncated
                ? '<span style="color:#b26a00">' + esc(truncationNotice(built, res)) + '</span>'
                : '<span style="color:#888">No contigs were assembled.</span>';
        }

        var H = [];
        H.push('<span style="color:#888">Key: </span>' +
               '<span style="color:#9e9e9e">.</span><span style="color:#888"> = agrees with the consensus · </span>' +
               '<span style="color:#d32f2f">x</span><span style="color:#888"> = disagrees · </span>' +
               '<span style="color:#c0c0c0">-</span><span style="color:#888"> = gap · </span>' +
               '<span style="color:#7b1fa2;font-weight:600">R</span><span style="color:#888"> = polymorphic consensus (IUPAC)</span>\n\n');

        built.blocks.forEach(function (B) {
            var c = B.contig;
            var idw = Math.max(12, Math.min(28, Math.max.apply(null, B.rows.map(function (r) { return r.id.length; }).concat([12]))));
            H.push('<span style="color:#0d6efd;font-weight:600">Contig' + c.id + '</span>' +
                   '<span style="color:#555">  ' + c.consensus.length + ' bp · ' + c.reads.length +
                   ' reads · mean coverage ' + pct(c.meanCoverage) + '× · ' + B.nCols + ' columns</span>\n');

            for (var start = 0; start < B.nCols; start += width) {
                var end = Math.min(start + width, B.nCols);
                H.push('\n<span style="color:#888">columns ' + (start + 1) + '–' + end + '</span>\n');

                // consensus row
                H.push('<span style="color:#555">' + esc(pad('consensus', idw)) + '   </span>');
                var paint = runPainter(H), k;
                for (k = start; k < end; k++) {
                    var ch = B.consLine[k];
                    if (B.polyAt[k]) paint.put(ch.toUpperCase(), 'color:#7b1fa2;font-weight:700');
                    else if (ch === '-') paint.put('-', 'color:#c0c0c0');
                    else paint.put(ch, 'color:#111');
                }
                paint.flush();
                H.push('\n');

                // read rows
                B.rows.forEach(function (row) {
                    if (row.last < start || row.first >= end) return;
                    H.push('<span style="color:#555">' + esc(pad(row.id.substring(0, idw), idw)) + ' </span>');
                    H.push('<span style="color:' + (row.strand === '-' ? '#b26a00' : '#2e7d32') + '">' + row.strand + '</span> ');
                    var rp = runPainter(H);
                    for (var q = start; q < end; q++) {
                        var cc = rowCharAt(row, q, B.consLine[q]);
                        if (cc === ' ') rp.put(' ', null);
                        else if (cc === '.') rp.put('.', 'color:#b0b0b0');
                        else if (cc === '-') rp.put('-', 'color:#c0c0c0');
                        else rp.put(cc, 'color:#d32f2f;font-weight:600');
                    }
                    rp.flush();
                    H.push('\n');
                });

                var mline = '';
                for (var m = start; m < end; m++) mline += B.polyAt[m] ? '*' : ' ';
                if (mline.replace(/ /g, '')) {
                    H.push('<span style="color:#7b1fa2">' + esc(pad('', idw) + '   ' + mline) + '</span>\n');
                }
            }
            H.push('\n');
        });
        if (built.truncated) {
            H.push('<span style="color:#b26a00">' + esc(truncationNotice(built, res)) + '</span>\n');
        }
        return H.join('');
    }

    /* ───────────────────────────── log / summary ─────────────────────────── */

    function capReportLog(res) {
        if (!res) return '';
        if (!res.ok) return (res.message || 'No result.') + '\n' +
                            (res.notes && res.notes.length ? '\n' + res.notes.join('\n') : '');
        var S = res.stats, L = [];
        L.push('=== OLC Assembler ===');
        L.push('');
        L.push('Reads used:        ' + S.readsIn);
        L.push('Placed / singleton:' + ' ' + S.readsPlaced + ' / ' + S.singletons +
               (S.chimeras ? '  (' + S.chimeras + ' chimeric read(s) set aside)' : ''));
        if (S.foldbacks || S.foldbacksDropped) {
            L.push('Fold-back reads:   ' + S.foldbacks + ' trimmed' + (S.foldbacksDropped ? ', ' + S.foldbacksDropped + ' dropped' : ''));
        }
        L.push('Contigs:           ' + S.contigs + '  (longest ' + S.longest + ' bp, N50 ' + S.n50 + ' bp)');
        L.push('Consensus total:   ' + S.totalBp + ' bp' + (S.trimmedBp ? '  (' + S.trimmedBp + ' bp of thinly covered contig ends trimmed)' : ''));
        L.push('Mean coverage:     ' + pct(S.meanCoverage) + '×');
        L.push('Polymorphic sites: ' + S.snps + ' SNP, ' + S.indels + ' indel' +
               (S.ties ? ', ' + S.ties + ' unresolved tie(s)' : ''));
        L.push('');
        L.push('=== Overlap stage ===');
        L.push('');
        L.push('k-mer words indexed: ' + res.seedStats.words);
        L.push('Seed hits:           ' + res.seedStats.seeds);
        L.push('Repeat words skipped:' + ' ' + res.seedStats.skippedWords +
               (res.seedStats.deepCut ? '  (seen more than ' + res.seedStats.deepCut + ' times)' : ''));
        if (res.seedStats.deepWords) {
            L.push('Deep-coverage words: ' + res.seedStats.deepWords +
                   '  (seen ' + (res.params.maxSeedOcc + 1) + '-' + res.seedStats.deepCut +
                   ' times; each read paired with a fixed number of the others)');
        }
        L.push('Candidate pairs:     ' + res.candidates);
        L.push('Accepted overlaps:   ' + S.overlaps);
        L.push('Rejected overlaps:   ' + S.rejected);
        L.push('Layout conflicts:    ' + S.conflicts);
        L.push('');
        L.push('=== Error model ===');
        L.push('');
        L.push('Sequencing error rate: ' + (100 * (S.errorRate || 0)).toFixed(2) + '%  ' +
               '(measured from the mismatches in this run\'s own overlaps)');
        L.push('Polymorphism cut-off:  p < ' + (S.alpha || 0).toExponential(2) +
               (res.params.polyMaxP > 0 ? '  (set manually)' : '  (auto: about one false call per assembly)'));
        L.push('Repeat cut-off:        word copy number ≥ ' +
               (isFinite(res.seedStats.repeatThreshold) ? res.seedStats.repeatThreshold : 'n/a') +
               '  (' + res.params.repeatFactor + '× the median)');
        L.push('');
        L.push('=== Performance ===');
        L.push('');
        L.push('Elapsed: ' + (S.elapsedMs / 1000).toFixed(2) + ' s');
        L.push('Threads: ' + S.threads + (S.threads > 1 ? ' (Web Workers)' : ' (single thread)'));
        if (res.notes.length) {
            L.push('');
            L.push('=== Notices ===');
            L.push('');
            res.notes.forEach(function (n) { L.push('- ' + n); });
        }
        return L.join('\n');
    }

    var api = {
        text: capReportText,
        fasta: capReportFasta,
        placements: capReportPlacements,
        alignmentText: capReportAlignmentText,
        alignmentHtml: capReportAlignmentHtml,
        log: capReportLog
    };
    if (typeof module !== 'undefined' && module.exports) module.exports = api;
    if (global) global.CapReport = api;

})(typeof window !== 'undefined' ? window : (typeof self !== 'undefined' ? self : globalThis));
