#!/usr/bin/env python
"""Extract plotted data from the vector figures of He & Kruczenski, arXiv:2403.10772v4.

Method (same as references/figure*_*.metadata.json for 2309.12402):
  * PyMuPDF vector paths (page.get_drawings(extended=True)); a marker centre is the
    bounding-box centre of its filled circular path; a line is the vertex list of a stroke.
  * Each axis is calibrated by a least-squares linear fit of >= 2 printed tick labels
    against the PDF coordinate of the tick-mark path; residuals and the position of the
    unlabelled minor ticks are recorded as checks.
      - Mathematica PDFs: tick labels are text spans (font ArialMT), read automatically.
      - Inkscape PDFs (psandeta, D0beforeW, D0afterW): tick labels are glyph outlines.
        The label strings are typed in below from the rendered figure and verified by
        (a) glyph count per label and (b) shape consistency: every glyph assigned the
        same character must have the same normalised outline, different characters
        must differ.
  * Series are identified by exact fill/stroke RGB (and dash pattern); the colour ->
    series map comes from the figure legend when there is one, otherwise from the
    paper text (blue/green/red = the three points of Fig. 3 = ni 44/46/48).
  * Clipping: the effective clip rectangle of every path is reconstructed from the
    clip stack; markers whose centre is outside it (or outside the panel frame) are
    kept with visible=0, line vertices outside it are dropped (count recorded).
  * Node self-check: every marker abscissa is compared with the collocation nodes
    sqrt(s_j), s_j = -20 + 48/(1+cos((j-1/2)pi/50)) [m_pi^2], m_pi = 139.57 MeV.
  * Completeness: every coloured marker / curve in a PDF must be assigned to a series
    (assertion), so nothing data-like is silently skipped.
  * 90-degree crossings (Fig. 4 P1, D0; Fig. 20 D0): linear interpolation in sqrt(s)
    between adjacent markers of one colour.
  * Cross-figure checks (index.json): the same run drawn in two figures (Fig. 4 vs 20b,
    Fig. 8 vs 9) must agree, which tests both calibrations independently.

Inputs: only the figure PDFs of the arXiv v4 source (nothing from any reproduction run).
Writes, next to this script, one CSV + one .metadata.json per panel and index.json.
CSV columns: series,x,y,kind (dot|line), then node_j, dx_node_gev, dx_node_pt, visible,
[y_err_lo, y_err_hi,] segment, pdf_x, pdf_y, marker_diam_pt, path_index.
Usage: python extract_2403_figures.py        (PyMuPDF >= 1.24, numpy; ~10 s)
"""
import csv
import hashlib
import json
import math
import os
import re
import sys

import fitz
import numpy as np

SRC = '/playpen1/shiqiu/gtb2403_prep/paper/src_2403v4/figures'
OUT = os.path.dirname(os.path.abspath(__file__))
MPI_GEV = 0.13957
NODE_M, NODE_NU0 = 50, -20.0
PAPER = 'He & Kruczenski, arXiv:2403.10772v4 (tex: arxiv_submission_3_v10.tex)'

J_NODES = np.arange(1, NODE_M + 1)
S_NODES = NODE_NU0 + 48.0 / (1.0 + np.cos((J_NODES - 0.5) * np.pi / NODE_M))
SQRT_S_NODES_GEV = np.sqrt(S_NODES) * MPI_GEV

BLACK = (0.0, 0.0, 0.0)
JOIN_GAP_PT, JOIN_GAP_REL, JOIN_DY_PT = 0.5, 0.0015, 3.0
DUPLICATES = {}   # series -> markers with an identical centre already read (omitted)


# ----------------------------------------------------------------------------- utilities
def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()


def cdist(a, b):
    return max(abs(x - y) for x, y in zip(a, b))


def is_gray(c, lo=0.3, hi=0.8):
    return c is not None and max(c) - min(c) < 0.02 and lo < c[0] < hi


def rnd(v, n=6):
    return None if v is None else float(round(float(v), n))


class Path:
    __slots__ = ('idx', 'type', 'fill', 'stroke', 'width', 'dashes', 'kinds', 'items',
                 'rect', 'scissor', 'subpaths')

    def __init__(self, idx, d, scissor):
        self.idx = idx
        self.type = d['type']
        self.fill = tuple(d['fill']) if d.get('fill') is not None else None
        self.stroke = tuple(d['color']) if d.get('color') is not None else None
        self.width = d.get('width')
        self.dashes = d.get('dashes') or ''
        self.items = d['items']
        self.kinds = [it[0] for it in d['items']]
        self.rect = fitz.Rect(d['rect'])
        self.scissor = scissor
        self.subpaths = self._subpaths()

    def _subpaths(self):
        subs, cur = [], []
        for it in self.items:
            if it[0] == 'l':
                p1, p2 = it[1], it[2]
            elif it[0] == 'c':
                p1, p2 = it[1], it[4]
            else:
                continue
            if cur and (abs(cur[-1][0] - p1.x) > 1e-6 or abs(cur[-1][1] - p1.y) > 1e-6):
                subs.append(cur)
                cur = []
            if not cur:
                cur.append((p1.x, p1.y))
            cur.append((p2.x, p2.y))
        if cur:
            subs.append(cur)
        return subs

    @property
    def center(self):
        return ((self.rect.x0 + self.rect.x1) / 2, (self.rect.y0 + self.rect.y1) / 2)

    @property
    def dashed(self):
        return bool(re.search(r'\[\s*\d', self.dashes))

    def is_disc(self, dmax=15.0):
        r = self.rect
        return (self.fill is not None and len(self.kinds) >= 4 and all(k == 'c' for k in self.kinds)
                and r.width < dmax and abs(r.width - r.height) < 0.03 * r.width + 0.01)

    def is_segment(self):
        """single straight black-or-coloured segment"""
        return self.type == 's' and self.kinds == ['l']

    def seg(self):
        it = self.items[0]
        return (it[1].x, it[1].y, it[2].x, it[2].y)


def load_paths(pdf):
    page = fitz.open(pdf)[0]
    raw = page.get_drawings(extended=True)
    stack, paths = [], []
    for k, d in enumerate(raw):
        lvl = d.get('level', 0)
        while stack and stack[-1][0] >= lvl:
            stack.pop()
        if d['type'] == 'clip':
            stack.append((lvl, fitz.Rect(d['scissor'])))
            continue
        if d['type'] == 'group':
            continue
        sc = None
        for _, r in stack:
            sc = fitz.Rect(r) if sc is None else (sc & r)
        paths.append(Path(k, d, sc))
    return page, paths


def inside(rect, x, y, tol=0.0):
    return rect is None or (rect.x0 - tol <= x <= rect.x1 + tol and rect.y0 - tol <= y <= rect.y1 + tol)


# ----------------------------------------------------------------------------- calibration
def fit_axis(coords, values):
    coords, values = np.asarray(coords, float), np.asarray(values, float)
    a, b = np.polyfit(coords, values, 1)
    res_pt = coords - (values - b) / a
    return float(a), float(b), [float(r) for r in res_pt]


def minor_check(a, b, majors_pos, minors_pos):
    """distance (pt) of each unlabelled tick from the nearest multiple of the minor step"""
    if not minors_pos or len(majors_pos) < 2:
        return None
    mv = sorted(a * np.asarray(majors_pos) + b)
    step_major = float(np.median(np.diff(mv)))
    # number of minor intervals per major interval from the tick density
    span = abs(a) * (max(majors_pos) - min(majors_pos))
    nin = sum(1 for p in minors_pos if min(majors_pos) < p < max(majors_pos))
    nsub = round(nin / max(len(majors_pos) - 1, 1)) + 1
    step = step_major / nsub
    dev = [abs((a * p + b) / step - round((a * p + b) / step)) * step / abs(a) for p in minors_pos]
    return dict(minor_step=step, n_minor=len(minors_pos), max_dev_pt=max(dev))


def axis_record(name, unit, coords, values, labels, path_idx, minors=None):
    a, b, res = fit_axis(coords, values)
    rec = dict(axis=name, unit=unit, formula=f'{name}_value = slope * pdf_{"x" if name=="x" else "y"} + intercept',
               slope=a, intercept=b, printed_tick_labels=labels, tick_values=[float(v) for v in values],
               tick_pdf_coord=[float(c) for c in coords], tick_path_indices=path_idx,
               tick_fit_residual_pt=res, max_abs_residual_pt=max(abs(r) for r in res),
               pt_per_unit=1.0 / abs(a))
    if minors is not None:
        rec['minor_tick_check'] = minor_check(a, b, coords, minors)
    return rec


# --- Mathematica (text tick labels, Axes style)
def text_spans(page):
    spans = []
    for b in page.get_text('rawdict')['blocks']:
        for line in b.get('lines', []):
            for s in line['spans']:
                t = ''.join(c['c'] for c in s['chars'])
                spans.append(dict(text=t, bbox=fitz.Rect(s['bbox']), font=s['font'], size=s['size']))
    return spans


NUM = re.compile(r'^-?\d+(\.\d+)?$')


def numeric_labels(spans):
    out = []
    for s in spans:
        t = s['text'].strip()
        if s['font'] != 'ArialMT' or not NUM.match(t):
            continue
        bb = fitz.Rect(s['bbox'])
        for m in spans:  # separate minus sign immediately to the left
            if m['text'].strip() in ('-', '−', '–') and 0 <= bb.x0 - m['bbox'].x1 < 4 \
                    and min(bb.y1, m['bbox'].y1) - max(bb.y0, m['bbox'].y0) > 0.5 * bb.height:
                t = '-' + t
                bb = bb | m['bbox']
        out.append(dict(text=t, value=float(t), bbox=bb))
    return out


def mma_axes(paths, spans):
    segs = [p for p in paths if p.is_segment() and p.stroke is not None and cdist(p.stroke, BLACK) < 1e-3]
    hs = [p for p in segs if abs(p.seg()[1] - p.seg()[3]) < 1e-3]
    vs = [p for p in segs if abs(p.seg()[0] - p.seg()[2]) < 1e-3]
    xax = max(hs, key=lambda p: p.rect.width)
    yax = max(vs, key=lambda p: p.rect.height)
    y0 = xax.seg()[1]
    x0 = yax.seg()[0]
    xt = [p for p in vs if p is not yax and p.rect.height < 15 and min(abs(p.rect.y0 - y0), abs(p.rect.y1 - y0)) < 0.05]
    yt = [p for p in hs if p is not xax and p.rect.width < 15 and min(abs(p.rect.x0 - x0), abs(p.rect.x1 - x0)) < 0.05]
    labels = numeric_labels(spans)
    res = {}
    for name, ticks, along in (('x', xt, 0), ('y', yt, 1)):
        L = max((p.rect.height if name == 'x' else p.rect.width) for p in ticks)
        major = [p for p in ticks if (p.rect.height if name == 'x' else p.rect.width) > 0.8 * L]
        minor = [p for p in ticks if p not in major]
        coords, vals, labs, idx = [], [], [], []
        used = set()
        for p in sorted(major, key=lambda p: p.center[along]):
            c = p.center[along]
            best = None
            for k, lab in enumerate(labels):
                bb = lab['bbox']
                if name == 'x' and not (bb.y0 > y0 and bb.y0 - y0 < 40):
                    continue
                if name == 'y' and not (bb.x1 < x0 + 1 and x0 - bb.x1 < 40):
                    continue
                lc = (bb.x0 + bb.x1) / 2 if name == 'x' else (bb.y0 + bb.y1) / 2
                tol = 0.3 * (bb.width if name == 'x' else bb.height)
                if abs(lc - c) < tol and (best is None or abs(lc - c) < best[0]):
                    best = (abs(lc - c), k)
            if best is None:
                continue
            assert best[1] not in used
            used.add(best[1])
            coords.append(c)
            vals.append(labels[best[1]]['value'])
            labs.append(labels[best[1]]['text'])
            idx.append(p.idx)
        assert len(coords) >= 2, f'fewer than two labelled {name} ticks'
        res[name] = (coords, vals, labs, idx, [p.center[along] for p in minor])
    frame = fitz.Rect(min(xax.rect.x0, x0), min(yax.rect.y0, y0), max(xax.rect.x1, x0), max(yax.rect.y1, y0))
    return res, dict(x_axis_path=xax.idx, y_axis_path=yax.idx, x_axis_pdf_y=y0, y_axis_pdf_x=x0), frame


# --- Inkscape (frames, glyph tick labels)
def find_frames(paths):
    blk = [p for p in paths if p.is_segment() and p.stroke is not None and cdist(p.stroke, BLACK) < 1e-3]
    frames = []
    for i in range(len(blk) - 3):
        g = blk[i:i + 4]
        if not all(max(p.rect.width, p.rect.height) > 30 for p in g):
            continue
        xs = [c for p in g for c in (p.seg()[0], p.seg()[2])]
        ys = [c for p in g for c in (p.seg()[1], p.seg()[3])]
        if len({round(x, 2) for x in xs}) == 2 and len({round(y, 2) for y in ys}) == 2 \
                and all(g[k].idx == g[0].idx + k for k in range(4)):
            frames.append((fitz.Rect(min(xs), min(ys), max(xs), max(ys)), [p.idx for p in g]))
    return frames, blk


def glyph_paths(paths):
    return [p for p in paths if p.type == 'f' and p.fill is not None and cdist(p.fill, BLACK) < 1e-3
            and not p.is_disc() and p.rect.height < 20 and p.rect.width < 14]


def cluster_glyphs(glyphs):
    n = len(glyphs)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i in range(n):
        for j in range(i + 1, n):
            a, b = glyphs[i].rect, glyphs[j].rect
            yo = min(a.y1, b.y1) - max(a.y0, b.y0)
            gap = max(a.x0, b.x0) - min(a.x1, b.x1)
            if yo > -0.5 and gap < 0.7 * max(a.height, b.height):
                parent[find(i)] = find(j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(glyphs[i])
    out = []
    for g in groups.values():
        g.sort(key=lambda p: p.rect.x0)
        bb = fitz.Rect(g[0].rect)
        for p in g[1:]:
            bb |= p.rect
        out.append(dict(glyphs=g, bbox=bb))
    return out


def glyph_signature(p):
    """outline points relative to the glyph bbox corner, in units of the glyph height"""
    r = p.rect
    h = max(r.height, 1e-6)
    pts = []
    for it in p.items:
        for q in it[1:]:
            if hasattr(q, 'x'):
                pts += [(q.x - r.x0) / h, (q.y - r.y0) / h]
            elif hasattr(q, 'x0'):
                pts += [(q.x0 - r.x0) / h, (q.y0 - r.y0) / h, (q.x1 - r.x0) / h, (q.y1 - r.y0) / h]
    return (tuple(it[0] for it in p.items), r.width / h, np.array(pts))


def check_glyph_chars(assigned):
    """assigned: list of (char, Path). Same char => same outline; different chars => different.
    Tolerance 0.01 glyph heights for outline glyphs, 0.06 for single-rectangle glyphs ('.', '-'), whose
    normalised form is only the aspect ratio of a < 1 pt high box (rounding 1e-3 pt)."""
    sig = {}
    for ch, p in assigned:
        sig.setdefault(ch, []).append(glyph_signature(p))
    within = {}
    for ch, ss in sig.items():
        d = 0.0
        for s in ss[1:]:
            assert s[0] == ss[0][0] and s[2].shape == ss[0][2].shape, f'glyph outline mismatch for {ch!r}'
            d = max(d, float(np.max(np.abs(s[2] - ss[0][2]))), abs(s[1] - ss[0][1]))
        within[ch] = dict(instances=len(ss), max_normalised_dev=d, n_items=len(ss[0][0]))
        assert d < (0.06 if len(ss[0][0]) == 1 else 0.01), f'same-character glyphs differ for {ch!r}: {d}'
    between = 1e9
    chars = list(sig)
    for i, a in enumerate(chars):
        for b in chars[i + 1:]:
            sa, sb = sig[a][0], sig[b][0]
            if sa[0] == sb[0] and sa[2].shape == sb[2].shape:
                between = min(between, float(np.max(np.abs(sa[2] - sb[2]))), abs(sa[1] - sb[1]) + 0)
    worst = max(v['max_normalised_dev'] for v in within.values())
    assert between > 0.1, f'two different characters have the same outline (min diff {between})'
    return dict(chars=within, min_between_char_dev_same_item_structure=None if between == 1e9 else between,
                max_within_char_dev=worst)


def frame_ticks(blk, frame):
    yt, xt = [], []
    for p in blk:
        x0, y0, x1, y1 = p.seg()
        L = math.hypot(x1 - x0, y1 - y0)
        if L > 5 or L < 0.2:
            continue
        if abs(y0 - y1) < 1e-3 and min(abs(x0 - frame.x0), abs(x1 - frame.x0)) < 0.05 \
                and frame.y0 - 0.1 <= y0 <= frame.y1 + 0.1:
            yt.append((y0, L, p.idx))
        if abs(x0 - x1) < 1e-3 and min(abs(y0 - frame.y1), abs(y1 - frame.y1)) < 0.05 \
                and frame.x0 - 0.1 <= x0 <= frame.x1 + 0.1:
            xt.append((x0, L, p.idx))
    return xt, yt


def inkscape_axis(name, ticks, clusters, frame, label_strings, glyph_log):
    L = max(t[1] for t in ticks)
    major = sorted([t for t in ticks if t[1] > 0.8 * L], key=lambda t: t[0], reverse=(name == 'y'))
    minor = [t[0] for t in ticks if t[1] <= 0.8 * L]
    coords, idx, labs = [], [], []
    for t in major:
        best = None
        for c in clusters:
            bb = c['bbox']
            if name == 'y':
                if not (bb.x1 <= frame.x0 + 0.1 and frame.x0 - bb.x1 < 8):
                    continue
                d = abs((bb.y0 + bb.y1) / 2 - t[0])
            else:
                if not (bb.y0 >= frame.y1 - 0.1 and bb.y0 - frame.y1 < 10):
                    continue
                d = abs((bb.x0 + bb.x1) / 2 - t[0])
            if d < 3.0 and (best is None or d < best[0]):
                best = (d, c)
        if best is None:
            continue
        coords.append(t[0])
        idx.append(t[2])
        labs.append(best[1])
    # majors ordered bottom->top (y) or left->right (x); labels typed in the same order
    assert len(labs) == len(label_strings), f'{name}: {len(labs)} labelled majors vs typed {label_strings}'
    for c, s in zip(labs, label_strings):
        assert len(c['glyphs']) == len(s), f'label {s!r}: {len(c["glyphs"])} glyphs'
        for ch, g in zip(s, c['glyphs']):
            glyph_log.append((ch, g))
    vals = [float(s) for s in label_strings]
    return coords, vals, list(label_strings), idx, minor


# ----------------------------------------------------------------------------- data elements
def node_check(xgev, a_x):
    k = int(np.argmin(np.abs(SQRT_S_NODES_GEV - xgev)))
    dx = xgev - SQRT_S_NODES_GEV[k]
    return int(J_NODES[k]), float(dx), float(dx / abs(a_x))


def collect_errorbars(paths, dots_xy, colour):
    """vertical gray single segments through a gray marker centre -> (ylo, yhi) pdf"""
    vs = [p for p in paths if p.is_segment() and p.stroke is not None and cdist(p.stroke, colour) < 1e-3
          and abs(p.seg()[0] - p.seg()[2]) < 1e-3]
    used = set()
    out = []
    for (x, y) in dots_xy:
        cand = [p for p in vs if abs(p.seg()[0] - x) < 0.35]
        lo, hi = y, y
        changed = True
        while changed:
            changed = False
            for p in cand:
                a, b = sorted((p.seg()[1], p.seg()[3]))
                if a <= hi + 0.6 and b >= lo - 0.6 and (a < lo or b > hi):
                    lo, hi = min(lo, a), max(hi, b)
                    changed = True
        for p in cand:
            a, b = sorted((p.seg()[1], p.seg()[3]))
            if a >= lo - 1e-6 and b <= hi + 1e-6:
                used.add(p.idx)
        out.append((lo, hi) if hi - lo > 1e-6 else (None, None))
    orphans = [p for p in vs if p.idx not in used]
    return out, orphans


def crossings(xs, ys, level):
    out = []
    for k in range(len(xs) - 1):
        y0, y1 = ys[k] - level, ys[k + 1] - level
        if y0 == 0:
            out.append((xs[k], 'on_node'))
        elif y0 * y1 < 0:
            t = y0 / (y0 - y1)
            out.append((xs[k] + t * (xs[k + 1] - xs[k]), 'up' if y1 > y0 else 'down'))
    return out


# ----------------------------------------------------------------------------- writers
def write_panel(stem, rows, extra_cols, meta):
    cols = ['series', 'x', 'y', 'kind'] + extra_cols
    with open(os.path.join(OUT, stem + '.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([('' if r.get(c) is None else (f'{r[c]:.10g}' if isinstance(r[c], float) else r[c]))
                        for c in cols])
    meta['dataset'] = dict(csv_file=stem + '.csv', columns=cols, row_count=len(rows),
                           CSV_format='10 significant digits preserve the graphical conversion; '
                                      'they are not physical or numerical accuracy.')
    with open(os.path.join(OUT, stem + '.metadata.json'), 'w') as f:
        json.dump(meta, f, indent=1)


def source_block(fig, fname, panel, page):
    pdf = os.path.join(SRC, fname)
    return dict(paper=PAPER, figure=fig, panel=panel, pdf=pdf, pdf_sha256=sha256(pdf),
                pdf_page_size_pt=[page.rect.width, page.rect.height],
                pdf_producer=fitz.open(pdf).metadata.get('creator'))


SCOPE = [
    'Graphical readout of the displayed marker centres and line vertices; not the authors\' numerical arrays.',
    'No fit, smoothing, rescaling, extrapolation or invented samples.',
    'Tick-fit and node residuals describe the PDF coordinate conversion only, not a physics uncertainty.',
    'Extracted independently of any reproduction result (frozen before comparison).',
]


def series_block(name, colour, kind, recs, basis, extra=None):
    xs = [r['x'] for r in recs]
    ys = [r['y'] for r in recs]
    b = dict(series=name, kind=kind, pdf_rgb=[rnd(c, 6) for c in colour] if colour else None,
             identification=basis, count=len(recs),
             x_range=[rnd(min(xs)), rnd(max(xs))] if xs else None,
             y_range=[rnd(min(ys)), rnd(max(ys))] if ys else None,
             source_path_indices=sorted({r['path_index'] for r in recs}))
    if kind == 'dot':
        b['visible_count'] = sum(1 for r in recs if r.get('visible') == 1)
        b['duplicate_centres_omitted'] = DUPLICATES.pop(name, [])
    if extra:
        b.update(extra)
    return b


def node_stats(recs):
    """node self-check of one marker series, plus markers whose outline differs from the series' modal marker
    (diameter off by > 0.005 pt) or which are drawn outside the series' contiguous block of paths"""
    js = [r['node_j'] for r in recs]
    dxg = [abs(r['dx_node_gev']) for r in recs]
    dxp = [abs(r['dx_node_pt']) for r in recs]
    dup = sorted({j for j in js if js.count(j) > 1})
    lo, hi = (min(js), max(js)) if js else (None, None)
    diam = float(np.median([r['marker_diam_pt'] for r in recs])) if recs else None
    idx = sorted(r['path_index'] for r in recs)
    med = float(np.median(idx)) if idx else None
    anomal = []
    for r in recs:
        far = abs(r['path_index'] - med) > 3 * len(recs)
        if abs(r['marker_diam_pt'] - diam) > 0.005 or far:
            anomal.append(dict(node_j=r['node_j'], path_index=r['path_index'], marker_diam_pt=rnd(r['marker_diam_pt'], 4),
                               series_median_diam_pt=rnd(diam, 4), drawn_out_of_series_block=far,
                               dx_node_pt=rnd(r['dx_node_pt'], 4)))
    return dict(nodes_hit=sorted(set(js)), node_range=[lo, hi],
                missing_nodes_in_range=[int(j) for j in range(lo, hi + 1) if j not in js] if js else [],
                nodes_1_to_39_not_drawn=[j for j in range(1, 40) if j not in js],
                nodes_with_two_markers=dup,
                max_abs_dx_gev=rnd(max(dxg), 8) if dxg else None,
                rms_dx_gev=rnd(math.sqrt(sum(d * d for d in dxg) / len(dxg)), 8) if dxg else None,
                max_abs_dx_pt=rnd(max(dxp), 5) if dxp else None,
                marker_anomalies=anomal)


def dot_rows(name, dots, cal_x, cal_y, frame, extra=None):
    ax, bx = cal_x
    ay, by = cal_y
    rows = []
    seen = {}
    for p in dots:
        cx, cy = p.center
        key = (round(cx, 3), round(cy, 3))
        if key in seen:
            DUPLICATES.setdefault(name, []).append(dict(path_index=p.idx, same_center_as=seen[key]))
            continue
        seen[key] = p.idx
        x = ax * cx + bx
        y = ay * cy + by
        j, dg, dp = node_check(x, ax)
        vis = int(inside(p.scissor, cx, cy, 0.01) and inside(frame, cx, cy, 0.05))
        r = dict(series=name, x=x, y=y, kind='dot', node_j=j, dx_node_gev=dg, dx_node_pt=dp,
                 pdf_x=cx, pdf_y=cy, marker_diam_pt=p.rect.width, visible=vis, path_index=p.idx)
        if extra:
            r.update(extra(p))
        rows.append(r)
    rows.sort(key=lambda r: r['x'])
    return rows


def line_rows(name, lines, cal_x, cal_y, frame):
    """polylines of one series: identical re-drawn paths are skipped; consecutive pieces are joined when the
    next starts to the right of the previous end within max(JOIN_GAP_PT, JOIN_GAP_REL * plot width) in x and
    JOIN_DY_PT in y (the exporters split curves into pieces separated by one sampling step: 0.08-0.25 pt in
    psandeta, ~0.48 or ~1.0 pt in the Mathematica PDFs); vertices outside the effective clip / plot region are dropped and the
    polyline is split there into numbered segments"""
    ax, bx = cal_x
    ay, by = cal_y
    seen, chains = set(), []
    for p in lines:
        key = tuple((round(x, 3), round(y, 3)) for sp in p.subpaths for (x, y) in sp)
        if key in seen:      # identical polyline drawn twice
            continue
        seen.add(key)
        for sp in p.subpaths:
            pts = [(x, y, p.idx, inside(p.scissor, x, y, 0.3) and inside(frame, x, y, 0.3)) for (x, y) in sp]
            dx = pts[0][0] - chains[-1][-1][0] if chains else -1e9
            dy = pts[0][1] - chains[-1][-1][1] if chains else 1e9
            if -0.01 <= dx <= max(JOIN_GAP_PT, JOIN_GAP_REL * frame.width) and abs(dy) <= JOIN_DY_PT:
                chains[-1] += pts if math.hypot(dx, dy) > 1e-3 else pts[1:]
            else:
                chains.append(pts)
    rows, dropped, seg = [], 0, 0
    for ch in chains:
        cur = []
        for (u, v, idx, ok) in ch + [(None, None, None, False)]:
            if ok:
                cur.append((u, v, idx))
                continue
            if u is not None:
                dropped += 1
            if cur:
                seg += 1
                rows += [dict(series=name, x=ax * u_ + bx, y=ay * v_ + by, kind='line', segment=seg,
                              pdf_x=u_, pdf_y=v_, path_index=i_) for (u_, v_, i_) in cur]
                cur = []
    return rows, dropped, len(seen)


# ----------------------------------------------------------------------------- Inkscape figures
GTB3 = {'GTB_blue_ni44': 'blue', 'GTB_green_ni46': 'green', 'GTB_red_ni48': 'red'}
INK = {
    'psandeta.pdf': dict(
        fig='4', stem='fig04',
        colours={'GTB_blue_ni44': (0.0, 0.0, 1.0), 'GTB_green_ni46': (0.0, 0.658824, 0.0),
                 'GTB_red_ni48': (1.0, 0.0, 0.0)},
        # frames ordered: left column top->bottom, then right column top->bottom
        layout=[('S0', 'delta', '0', '0'), ('S0', 'eta', '0', '0'), ('S2', 'delta', '2', '0'), ('S2', 'eta', '2', '0'),
                ('P1', 'delta', '1', '1'), ('P1', 'eta', '1', '1'), ('D0', 'delta', '0', '2'), ('D0', 'eta', '0', '2'),
                ('D2', 'delta', '2', '2'), ('D2', 'eta', '2', '2'), ('F1', 'delta', '1', '3'), ('F1', 'eta', '1', '3')],
        ylabels={('S0', 'delta'): ['0', '50', '100', '150', '200'], ('D0', 'delta'): ['0', '50', '100', '150', '200'],
                 ('S2', 'delta'): ['-40', '-30', '-20', '-10', '0'], ('D2', 'delta'): ['-6', '-4', '-2', '0'],
                 ('P1', 'delta'): ['0', '50', '100', '150', '200', '250', '300', '350'],
                 ('F1', 'delta'): ['-5', '0', '5']},
        eta_ylabels=['0', '0.5', '1'], xlabels=['0.5', '1.0', '1.5', '2.0'],
        crossings_for=[('P1', 'delta'), ('D0', 'delta')]),
    'D0beforeW.pdf': dict(
        fig='20a', stem='fig20_D0beforeW',
        colours={'GTB_blue_ni44': (0.0, 0.0, 0.992157), 'GTB_green_ni46': (0.0, 0.658824, 0.0),
                 'GTB_red_ni48': (0.992157, 0.0, 0.0)},
        layout=[('D0', 'delta', '0', '2'), ('D0', 'eta', '0', '2')],
        ylabels={('D0', 'delta'): ['0', '50', '100', '150', '200']},
        eta_ylabels=['0', '0.5', '1'], xlabels=['0.5', '1.0', '1.5', '2.0'],
        crossings_for=[('D0', 'delta')]),
    'D0afterW.pdf': dict(
        fig='20b', stem='fig20_D0afterW',
        colours={'GTB_blue_ni44': (0.0, 0.0, 0.992157), 'GTB_green_ni46': (0.0, 0.658824, 0.0),
                 'GTB_red_ni48': (0.992157, 0.0, 0.0)},
        layout=[('D0', 'delta', '0', '2'), ('D0', 'eta', '0', '2')],
        ylabels={('D0', 'delta'): ['0', '50', '100', '150', '200']},
        eta_ylabels=['0', '0.5', '1'], xlabels=['0.5', '1.0', '1.5', '2.0'],
        crossings_for=[('D0', 'delta')]),
}


def do_inkscape(fname, cfg, index):
    page, paths = load_paths(os.path.join(SRC, fname))
    frames, blk = find_frames(paths)
    assert len(frames) == len(cfg['layout']), f'{fname}: {len(frames)} frames'
    frames.sort(key=lambda f: (round(f[0].x0 / 50), f[0].y0))
    glyphs = glyph_paths(paths)
    clusters = cluster_glyphs(glyphs)
    glyph_log = []
    panels = {}
    for (fr, fidx), (wave, q, iso, ell) in zip(frames, cfg['layout']):
        xt, yt = frame_ticks(blk, fr)
        ylab = cfg['ylabels'][(wave, q)] if q == 'delta' else cfg['eta_ylabels']
        yc = inkscape_axis('y', yt, clusters, fr, ylab, glyph_log)
        panels[(wave, q)] = dict(frame=fr, frame_paths=fidx, y=yc, x=None, xt=xt, iso=iso, ell=ell)
        if xt:
            panels[(wave, q)]['x'] = inkscape_axis('x', xt, clusters, fr, cfg['xlabels'], glyph_log)
    glyph_check = check_glyph_chars(glyph_log)
    # delta panels have no x ticks: use the x calibration of the eta panel drawn directly below
    for (wave, q), P in panels.items():
        if P['x'] is None:
            E = panels[(wave, 'eta')]
            assert abs(E['frame'].x0 - P['frame'].x0) < 0.01 and abs(E['frame'].x1 - P['frame'].x1) < 0.01
            P['x_from'] = f'{wave} eta panel (identical frame x-extent: dx0={E["frame"].x0-P["frame"].x0:+.4f} pt, ' \
                          f'dx1={E["frame"].x1-P["frame"].x1:+.4f} pt)'
            P['x'] = E['x']
            P['x_shared'] = True
    # element -> panel by containment of the centre in a frame grown by 4 pt
    def panel_of(x, y):
        best = None
        for key, P in panels.items():
            f = P['frame']
            d = max(f.x0 - x, x - f.x1, f.y0 - y, y - f.y1, 0)
            if d < 4 and (best is None or d < best[0]):
                best = (d, key)
        return None if best is None else best[1]

    dots = [p for p in paths if p.is_disc(8.0)]
    gstroke = [p for p in paths if p.type == 's' and is_gray(p.stroke)]
    # every marker must be a series colour or gray; no coloured curves exist in these PDFs
    odd = [p.idx for p in dots if not is_gray(p.fill) and min(cdist(p.fill, c) for c in cfg['colours'].values()) > 1e-3]
    odd += [p.idx for p in paths if p.type == 's' and p.stroke is not None and not is_gray(p.stroke)
            and cdist(p.stroke, BLACK) > 1e-3]
    assert not odd, f'{fname}: unclassified coloured paths {odd}'
    unassigned = []
    results = {}
    for key, P in panels.items():
        wave, q = key
        xc = axis_record('x', 'GeV (sqrt s)', *P['x'][:4], minors=P['x'][4])
        yc = axis_record('y', 'degrees' if q == 'delta' else 'dimensionless', *P['y'][:4], minors=P['y'][4])
        cal_x, cal_y = (xc['slope'], xc['intercept']), (yc['slope'], yc['intercept'])
        fr = P['frame']
        rows, groups = [], []
        mine = [p for p in dots if panel_of(*p.center) == key]
        for name, col in cfg['colours'].items():
            sel = [p for p in mine if cdist(p.fill, col) < 1e-3]
            r = dot_rows(name, sel, cal_x, cal_y, fr)
            rows += r
            extra = dict(node_check=node_stats(r))
            if q == 'eta':
                extra['max_abs_eta_minus_1'] = rnd(max(abs(x['y'] - 1) for x in r), 6) if r else None
            groups.append(series_block(name, col, 'dot', r,
                                       f'colour {GTB3[name]} = Fig. 3 point of the same colour (tex l.854: "the blue, '
                                       f'green and red points correspond to the points in fig. shape of the same color"); '
                                       f'ni={name[-2:]} from its Fig. 3 position F0 = ni/50 x 0.1632', extra))
        gd = [p for p in mine if is_gray(p.fill)]
        gcol = gd[0].fill if gd else None
        if gd:
            bars, orph = collect_errorbars([p for p in gstroke if panel_of(*p.center) == key],
                                           [p.center for p in gd], gcol)
            barmap = {p.idx: b for p, b in zip(gd, bars)}

            def ext(p, barmap=barmap):
                lo, hi = barmap[p.idx]
                return dict(y_err_lo=None if hi is None else cal_y[0] * hi + cal_y[1],
                            y_err_hi=None if lo is None else cal_y[0] * lo + cal_y[1])
            r = dot_rows('exp_data_gray', gd, cal_x, cal_y, fr, ext)
            for x in r:
                x['node_j'] = x['dx_node_gev'] = x['dx_node_pt'] = None
            rows += r
            groups.append(series_block('exp_data_gray', gcol, 'dot', r,
                                       'gray markers = experimental data [Protopopescu 73, Losty 74, Hyams 75] '
                                       '(caption, tex l.854 "gray dots"); marker diameter (CSV) may separate data sets, '
                                       'not attributed here',
                                       dict(marker_diameters_pt=sorted({round(p.rect.width, 2) for p in gd}),
                                            with_vertical_errorbar=sum(1 for x in r if x.get('y_err_lo') is not None),
                                            orphan_vertical_gray_segments=len(orph))))
        # lines: gray polylines with >= 3 vertices that are not straight reference lines
        lines = [p for p in gstroke if panel_of(*p.center) == key and sum(len(s) for s in p.subpaths) >= 3]
        straight = [p for p in paths if p.is_segment() and p.rect.width > 0.8 * fr.width and p.rect.height < 1e-3
                    and fr.y0 - 1 <= p.rect.y0 <= fr.y1 + 1 and p.idx not in P['frame_paths']
                    and fr.x0 - 0.2 * fr.width <= p.rect.x0 and p.rect.x1 <= fr.x1 + 1]
        if lines:
            lr, dropped, npaths = line_rows('PY_fit_gray', lines, cal_x, cal_y, fr)
            rows += lr
            groups.append(series_block('PY_fit_gray', lines[0].stroke, 'line', lr,
                                       'gray curve = phenomenological fit [Pelaez:2004vs] (tex l.854 "gray lines"; S0 '
                                       'with K Kbar removed, l.855)', dict(n_paths=npaths, vertices_dropped_outside_clip=dropped)))
        ref_lines = [dict(path_index=p.idx, pdf_y=p.rect.y0, value=cal_y[0] * p.rect.y0 + cal_y[1],
                          rgb=[rnd(c, 4) for c in p.stroke]) for p in straight]
        rows.sort(key=lambda r: (r['series'], r['kind'], r.get('segment') or 0, r['x']))
        stem = f'{cfg["stem"]}_{wave}_{q}' if len(panels) > 2 else f'{cfg["stem"]}_{q}'
        sym = ('delta' if q == 'delta' else 'eta') + f'^{P["iso"]}_{P["ell"]}'
        meta = dict(source=source_block(cfg['fig'], fname, f'{wave} {q} ({sym})', page),
                    extraction=dict(
                        method='PyMuPDF get_drawings(extended=True). Marker centre = bbox centre of filled circular '
                               'path; line = stroke vertices; panel = enclosing frame (4 black segments). y ticks '
                               'are black segments on the frame left edge; tick labels are glyph outlines whose '
                               'strings were typed from the rendered figure and verified by glyph count and '
                               'outline consistency (glyph_label_check).',
                        frame_pdf=[rnd(v, 4) for v in fr], frame_path_indices=P['frame_paths'],
                        frame_range_physical=dict(x=[rnd(cal_x[0] * fr.x0 + cal_x[1]), rnd(cal_x[0] * fr.x1 + cal_x[1])],
                                                  y=[rnd(cal_y[0] * fr.y1 + cal_y[1]), rnd(cal_y[0] * fr.y0 + cal_y[1])]),
                        x_calibration=xc, x_calibration_source=P.get('x_from', 'own x ticks'),
                        y_calibration=yc, glyph_label_check=glyph_check,
                        reference_lines_excluded=ref_lines,
                        node_formula='sqrt(s_j)*0.13957 GeV, s_j=-20+48/(1+cos((j-1/2)pi/50)), j=1..50'),
                    groups=groups, scope=SCOPE + [
                        'visible=0 marks a marker whose centre lies outside the panel frame/clip (not seen in the '
                        'paper); it is kept because its position is still defined by the same calibration.',
                        'y_err_lo/y_err_hi: extent of the vertical gray error bar through the marker (if any).',
                        'node_check.nodes_1_to_39_not_drawn: markers absent from the PDF. The figure software omits '
                        'points outside the plot range, so these values presumably lie outside frame_range_physical.y '
                        '(e.g. D0 at low energy, alternate nodes, delta presumably slightly below 0); not recoverable.'])
        if key in cfg['crossings_for']:
            cr = {}
            for name in cfg['colours']:
                rr = sorted([r for r in rows if r['series'] == name and r['kind'] == 'dot'], key=lambda r: r['x'])
                xs, ys = [r['x'] for r in rr], [r['y'] for r in rr]
                c90 = crossings(xs, ys, 90.0)
                c270 = crossings(xs, ys, 270.0)
                first = next((c for c in c90 if c[1] == 'up'), None)
                k = None if first is None else max(i for i, x in enumerate(xs) if x <= first[0])
                kmax = int(np.argmax(ys))
                cr[name] = dict(first_up_crossing_90deg_gev=rnd(first[0], 6) if first else None,
                                max_marker=[rr[kmax]['node_j'], rnd(xs[kmax], 6), rnd(ys[kmax], 4)],
                                bracketing_nodes=None if k is None else [rr[k]['node_j'], rr[k + 1]['node_j']],
                                bracketing_points=None if k is None else [[rnd(xs[k], 6), rnd(ys[k], 4)],
                                                                          [rnd(xs[k + 1], 6), rnd(ys[k + 1], 4)]],
                                all_90deg=[[rnd(c[0], 6), c[1]] for c in c90],
                                all_270deg=[[rnd(c[0], 6), c[1]] for c in c270])
            for name in ('PY_fit_gray',):
                lr = [r for r in rows if r['series'] == name]
                if lr:
                    out = []
                    for sgm in sorted({r['segment'] for r in lr}):
                        rr = sorted([r for r in lr if r['segment'] == sgm], key=lambda r: r['x'])
                        out += [[rnd(c[0], 6), c[1]] for c in crossings([r['x'] for r in rr], [r['y'] for r in rr], 90.0)]
                    cr[name] = dict(all_90deg=out)
            meta['crossings'] = dict(method='linear interpolation in sqrt(s) between adjacent markers of the same '
                                            'colour (sorted by sqrt s), level 90 deg (and 270 deg)', by_series=cr)
        meta['groups_summary'] = {g['series']: g['count'] for g in groups}
        cols = ['node_j', 'dx_node_gev', 'dx_node_pt', 'visible', 'y_err_lo', 'y_err_hi', 'segment',
                'pdf_x', 'pdf_y', 'marker_diam_pt', 'path_index']
        write_panel(stem, rows, cols, meta)
        results[stem] = meta
        index.append(dict(stem=stem, figure=cfg['fig'], file=fname, panel=f'{wave} {q}',
                          series={g['series']: g['count'] for g in groups},
                          x_tick_resid_pt=xc['max_abs_residual_pt'], y_tick_resid_pt=yc['max_abs_residual_pt'],
                          node_check={g['series']: g.get('node_check', {}).get('max_abs_dx_gev')
                                      for g in groups if g['series'].startswith('GTB')},
                          crossings=meta.get('crossings')))
    # anything colour-coded that was not assigned to a panel?
    for p in dots:
        if panel_of(*p.center) is None:
            unassigned.append(p.idx)
    assert not unassigned, f'{fname}: markers outside every frame: {unassigned}'
    # cross-check the shared x calibration: same node, same colour -> same pdf x in delta and eta panels
    xchk = {}
    for wave in {k[0] for k in panels}:
        st_d = f'{cfg["stem"]}_{wave}_delta' if len(panels) > 2 else f'{cfg["stem"]}_delta'
        st_e = f'{cfg["stem"]}_{wave}_eta' if len(panels) > 2 else f'{cfg["stem"]}_eta'
        dev = 0.0
        rd = list(csv.DictReader(open(os.path.join(OUT, st_d + '.csv'))))
        re_ = list(csv.DictReader(open(os.path.join(OUT, st_e + '.csv'))))
        for name in cfg['colours']:
            ed = {int(r['node_j']): float(r['x']) for r in re_ if r['series'] == name}
            for r in rd:
                if r['series'] == name and int(r['node_j']) in ed:
                    dev = max(dev, abs(float(r['x']) - ed[int(r['node_j'])]))
        xchk[wave] = dev
    for stem, meta in results.items():
        wave = meta['source']['panel'].split()[0]
        meta['extraction']['delta_vs_eta_same_node_max_dx_gev'] = rnd(xchk[wave], 9)
        with open(os.path.join(OUT, stem + '.metadata.json'), 'w') as f:
            json.dump(meta, f, indent=1)
    return results


# ----------------------------------------------------------------------------- Mathematica figures
MMA_BLUE, MMA_GREEN, MMA_RED = (0.0, 0.0, 0.999985), (0.0, 0.666656, 0.0), (0.999985, 0.0, 0.0)
C0, C1, C2 = (0.368408, 0.506775, 0.709793), (0.880707, 0.611038, 0.142044), (0.560165, 0.691559, 0.194885)
PT_TEXT = 'colour = point of Fig. 3 (blue ni=44, green ni=46, red ni=48; paper text / caption)'


def legend_rows(spans, paths):
    """legend swatches (coloured squares or short horizontal line samples) and the text to the right of each.
    Main-size spans are attached to the nearest swatch row; smaller (sub/superscript) spans to the row whose
    main text they overlap most in height."""
    sw = [p for p in paths if (p.type in ('f', 'fs') and p.kinds == ['re'] and p.fill is not None
                               and cdist(p.fill, BLACK) > 1e-3)
          or (p.is_segment() and p.stroke is not None and p.rect.height < 1e-3 and 10 < p.rect.width < 60)]
    if not sw:
        return [], set()
    ys = sorted(p.center[1] for p in sw)
    gaps = [b - a for a, b in zip(ys, ys[1:]) if b - a > 1]
    spacing = min(gaps) if gaps else 30.0
    xmin = min(p.rect.x1 for p in sw) - 0.5
    cand = [s for s in spans if not (s['font'] == 'ArialMT' and NUM.match(s['text'].strip()))
            and s['bbox'].x0 >= xmin and s['bbox'].x0 - xmin < 300
            and ys[0] - spacing < (s['bbox'].y0 + s['bbox'].y1) / 2 < ys[-1] + spacing]
    if not cand:
        return [], set()
    big = max(s['size'] for s in cand)
    rows = {p.idx: [] for p in sw}
    for s in cand:
        if s['size'] < 0.9 * big:
            continue
        cy = (s['bbox'].y0 + s['bbox'].y1) / 2
        p = min(sw, key=lambda p: abs(cy - p.center[1]))
        if abs(cy - p.center[1]) < 0.75 * spacing:
            rows[p.idx].append(s)
    band = {k: (min(s['bbox'].y0 for s in v), max(s['bbox'].y1 for s in v)) for k, v in rows.items() if v}
    for s in cand:
        if s['size'] >= 0.9 * big:
            continue
        ov = {k: min(b[1], s['bbox'].y1) - max(b[0], s['bbox'].y0) for k, b in band.items()}
        k = max(ov, key=ov.get)
        if ov[k] > 0:
            rows[k].append(s)
    out = []
    for p in sw:
        ss = sorted(rows[p.idx], key=lambda s: s['bbox'].x0)
        if not ss:
            continue
        out.append(dict(text=''.join(s['text'] for s in ss).replace(' ', ''), swatch_path=p.idx,
                        rgb=p.fill if p.type in ('f', 'fs') else p.stroke, dashed=p.dashed if p.type == 's' else None,
                        swatch='square' if p.type in ('f', 'fs') else 'line', swatch_pdf_y=p.center[1]))
    return out, {p.idx for p in sw}


MMA = [
    dict(stem='fig05_forwardRe', fig='5a', file='forwardRe.pdf', ylab='Re T(s,t=0) (dimensionless, as plotted)',
         legend={'ReT0+(s,t=0)': 'T0+', 'ReT00(s,t=0)': 'T00', 'ReTIt=1(s,t=0)': 'TIt1'},
         roles={'dot': '_GTB_red_ni48', 'solid': '_GTB_red_ni48_interp', 'dashed': '_PY'},
         basis='legend swatch colour -> amplitude; caption: points + solid lines = red point of Fig. 3, '
               'dashed = Pelaez-Yndurain with K Kbar removed'),
    dict(stem='fig05_forwardIm', fig='5b', file='forwardIm.pdf', ylab='Im T(s,t=0) (dimensionless, as plotted)',
         legend={'ImT0+(s,t=0)': 'T0+', 'ImT00(s,t=0)': 'T00', 'ImTIt=1(s,t=0)': 'TIt1'},
         roles={'dot': '_GTB_red_ni48', 'solid': '_GTB_red_ni48_interp', 'dashed': '_PY'},
         basis='legend swatch colour -> amplitude; caption: points + solid lines = red point of Fig. 3, '
               'dashed = Pelaez-Yndurain with K Kbar removed'),
    dict(stem='fig06_S0ff', fig='6', file='S0ff.pdf', ylab='|F^0_0| (dimensionless)',
         colours={MMA_BLUE: 'GTB_blue_ni44', MMA_GREEN: 'GTB_green_ni46', MMA_RED: 'GTB_red_ni48'},
         roles={'dot': '', 'solid': '_interp'}, basis=PT_TEXT),
    dict(stem='fig08_P1FF', fig='8', file='P1FF.pdf', ylab='|F^1_1| (dimensionless)',
         colours={MMA_BLUE: 'GTB_blue_ni44', MMA_RED: 'GTB_red_ni48',
                  (0.333328, 0.333328, 0.999985): 'GTB_blue_ni44', (0.999985, 0.333328, 0.333328): 'GTB_red_ni48'},
         roles={'dot': '', 'solid': '_GS_fit', 'dashed': '_BW_fit'},
         legend_check={'Breit-Wigner': True, 'Gounaris-Sakurai': False},
         basis=PT_TEXT + '; legend: dashed = Breit-Wigner (5.2), solid = Gounaris-Sakurai (5.3)'),
    dict(stem='fig09_eehadron', fig='9', file='eehadron.pdf', ylab='|F^1_1| (dimensionless; green = rescaled R)',
         colours={MMA_BLUE: 'GTB_blue_ni44', BLACK: 'GS_PDG_m770_G149', MMA_GREEN: 'PDG_R_ratio_rescaled'},
         roles={'dot': '', 'solid': '_fit_line'}, black_role='',
         basis='legend: blue line "GTB", black line "m_rho=770 MeV, Gamma_rho=149 MeV"; caption: blue = "one of the '
               'fits (blue) in fig. 8", black = Gounaris-Sakurai with PDG values, green = PDG e+e- -> hadrons '
               'rescaled; blue markers = blue point (as Fig. 8)'),
    dict(stem='fig10_rho1450twopole', fig='10', file='rho1450twopole.pdf', ylab='|calF^1_1|^2 (dimensionless)',
         colours={MMA_BLUE: 'GTB_blue_ni44', MMA_RED: 'GTB_red_ni48'},
         roles={'dot': '', 'solid': '_twoBW_fit'}, basis=PT_TEXT + '; lines = two-Breit-Wigner fit (5.9)'),
    dict(stem='fig11_D0FF', fig='11', file='D0FF.pdf', ylab='|F^0_2| (dimensionless)',
         colours={MMA_BLUE: 'GTB_blue_ni44', MMA_RED: 'GTB_red_ni48', BLACK: 'BW_PDG_f2'},
         roles={'dot': '', 'solid': '_BW_fit'}, black_role='',
         basis=PT_TEXT + '; caption: blue/red lines = BW fits, black = BW with PDG f2 parameters'),
]
for _f, _w, _lab in (('00', 'S0', 'r^0_0'), ('11', 'P1', 'r^1_1'), ('02', 'D0', 'r^0_2')):
    for _ab, _q in (('A', f'|F| {_w}'), ('B', f'alpha {_w} (degrees)')):
        MMA.append(dict(stem=f'fig18_Fasymp{_f}{_ab}', fig='18', file=f'Fasymp{_f}{_ab}.pdf', ylab=_q,
                        legend='auto', roles={'dot': ''},
                        basis='legend swatch colour -> r value (caption: blue = minimal feasible, orange = used, '
                              'green = relaxed); GTB point not stated in the paper'))
for _f, _w in (('00', 'S0'), ('11', 'P1'), ('02', 'D0')):
    for _k in ('SRerr', 'nSR'):
        MMA.append(dict(stem=f'fig19_{_k}{_f}', fig='19', file=f'{_k}{_f}.pdf', ylab=f'|F| {_w}',
                        legend='auto', roles={'dot': ''},
                        basis='legend swatch colour -> epsilon or n_SR value; GTB point not stated in the paper'))


def legend_series_name(text):
    t = text.replace('⨯', 'x').replace('×', 'x').replace('ϵ', 'eps').replace('ε', 'eps')
    t = re.sub(r'10-(\d)', r'10^-\1', t)
    return re.sub(r'[^A-Za-z0-9=.^+x_-]', '', t)


def do_mma(cfg, index):
    fname = cfg['file']
    page, paths = load_paths(os.path.join(SRC, fname))
    spans = text_spans(page)
    ax, axinfo, frame = mma_axes(paths, spans)
    xc = axis_record('x', 'GeV (sqrt s)', *ax['x'][:4], minors=ax['x'][4])
    yc = axis_record('y', cfg['ylab'], *ax['y'][:4], minors=ax['y'][4])
    cal_x, cal_y = (xc['slope'], xc['intercept']), (yc['slope'], yc['intercept'])
    legend, swatch_idx = legend_rows(spans, paths)
    colours = dict(cfg.get('colours', {}))
    if cfg.get('legend') == 'auto':
        for L in legend:
            colours[tuple(L['rgb'])] = legend_series_name(L['text'])
    elif isinstance(cfg.get('legend'), dict):
        for txt, nm in cfg['legend'].items():
            hit = [L for L in legend if L['text'] == txt.replace(' ', '')]
            assert len(hit) == 1, f'{fname}: legend row {txt!r} not found in {[L["text"] for L in legend]}'
            colours[tuple(hit[0]['rgb'])] = nm
    if cfg.get('legend_check'):          # legend line style (dash) per text, P1FF
        for txt, dashed in cfg['legend_check'].items():
            hit = [L for L in legend if L['text'] == txt]
            assert hit and all(L['dashed'] == dashed for L in hit), f'{fname}: legend style for {txt}'
    if fname == 'eehadron.pdf':
        hit = {L['text']: L for L in legend}
        assert cdist(hit['GTB']['rgb'], MMA_BLUE) < 1e-3, hit
        assert cdist(hit['m\u03c1=770MeV,\u0393\u03c1=149MeV']['rgb'], BLACK) < 1e-3, hit

    def cname(c):
        for k, v in colours.items():
            if c is not None and cdist(c, k) < 2e-3:
                return v
        return None

    axis_paths = {axinfo['x_axis_path'], axinfo['y_axis_path']}
    rows, groups, unclassified = [], [], []
    dots = [p for p in paths if p.is_disc() and p.idx not in swatch_idx]
    by = {}
    for p in dots:
        n = cname(p.fill)
        if n is None:
            unclassified.append(p.idx)
            continue
        by.setdefault((n + cfg['roles']['dot'], p.fill), []).append(p)
    for (n, col), sel in by.items():
        r = dot_rows(n, sel, cal_x, cal_y, frame)
        rows += r
        extra = dict(node_check=node_stats(r)) if not n.startswith('PDG') else {}
        if n.startswith('PDG'):
            for x in r:
                x['node_j'] = x['dx_node_gev'] = x['dx_node_pt'] = None
        extra['marker_diameters_pt'] = sorted({round(p.rect.width, 2) for p in sel})
        groups.append(series_block(n, col, 'dot', r, cfg['basis'], extra))
    lines = [p for p in paths if p.type == 's' and p.idx not in swatch_idx and p.idx not in axis_paths
             and p.stroke is not None and (cdist(p.stroke, BLACK) > 1e-3 or sum(len(q) for q in p.subpaths) >= 6)]
    byl = {}
    for p in lines:
        n = cname(p.stroke)
        if n is None:
            unclassified.append(p.idx)
            continue
        role = cfg['roles']['dashed' if p.dashed else 'solid']
        if p.stroke is not None and cdist(p.stroke, BLACK) < 1e-3:
            role = cfg.get('black_role', role)
        byl.setdefault((n + role, p.stroke, p.dashed), []).append(p)
    for (n, col, dashed), sel in byl.items():
        lr, dropped, npaths = line_rows(n, sel, cal_x, cal_y, frame)
        rows += lr
        groups.append(series_block(n, col, 'line', lr, cfg['basis'],
                                   dict(dashed=dashed, dash_pattern=sel[0].dashes, n_paths=npaths,
                                        n_drawn_paths=len(sel), vertices_dropped_outside_clip=dropped)))
    # completeness: every coloured fill is a classified marker or a legend swatch, every coloured stroke a
    # classified curve or a legend sample; remaining black strokes are axes/ticks (<= 2 vertices)
    done = {p.idx for v in by.values() for p in v} | {p.idx for v in byl.values() for p in v} | swatch_idx
    for p in paths:
        if p.idx in done or p.idx in axis_paths:
            continue
        if p.fill is not None and cdist(p.fill, BLACK) > 1e-3:
            unclassified.append(p.idx)
        if p.type == 's' and p.stroke is not None and (cdist(p.stroke, BLACK) > 1e-3 or
                                                       sum(len(q) for q in p.subpaths) > 2):
            unclassified.append(p.idx)
    assert not unclassified, f'{fname}: unclassified coloured paths {sorted(set(unclassified))}'
    rows.sort(key=lambda r: (r['series'], r['kind'], r.get('segment') or 0, r['x']))
    meta = dict(source=source_block(cfg['fig'], fname, 'single', page),
                extraction=dict(
                    method='PyMuPDF get_drawings(extended=True). Marker centre = bbox centre of filled circular path; '
                           'line = stroke vertices inside the effective clip rectangle; axes = longest black '
                           'horizontal/vertical segments; ticks = black segments touching an axis; tick labels = '
                           'ArialMT text spans centred on the major ticks.',
                    axes=axinfo, plot_region_pdf=[rnd(v, 4) for v in frame],
                    x_calibration=xc, y_calibration=yc,
                    legend=[dict(text=L['text'], swatch=L['swatch'], swatch_path_index=L['swatch_path'],
                                 rgb=[rnd(c, 6) for c in L['rgb']], dashed=L['dashed']) for L in legend],
                    node_formula='sqrt(s_j)*0.13957 GeV, s_j=-20+48/(1+cos((j-1/2)pi/50)), j=1..50'),
                groups=groups, scope=SCOPE + ['visible=0: marker centre outside the effective clip region.'])
    meta['groups_summary'] = {g['series']: g['count'] for g in groups}
    cols = ['node_j', 'dx_node_gev', 'dx_node_pt', 'visible', 'segment', 'pdf_x', 'pdf_y', 'marker_diam_pt',
            'path_index']
    write_panel(cfg['stem'], rows, cols, meta)
    index.append(dict(stem=cfg['stem'], figure=cfg['fig'], file=fname,
                      series={g['series']: g['count'] for g in groups},
                      x_tick_resid_pt=xc['max_abs_residual_pt'], y_tick_resid_pt=yc['max_abs_residual_pt'],
                      node_check={g['series']: g['node_check']['max_abs_dx_gev'] for g in groups
                                  if 'node_check' in g}))
    return meta


# ----------------------------------------------------------------------------- cross-figure checks
def _dots(stem, series):
    with open(os.path.join(OUT, stem + '.csv')) as f:
        return {int(r['node_j']): float(r['y']) for r in csv.DictReader(f)
                if r['series'] == series and r['kind'] == 'dot' and r['node_j']}


def _line(stem, series):
    with open(os.path.join(OUT, stem + '.csv')) as f:
        rr = [(float(r['x']), float(r['y'])) for r in csv.DictReader(f) if r['series'] == series]
    rr.sort()
    return np.array([r[0] for r in rr]), np.array([r[1] for r in rr])


def _cmp(a, b):
    k = sorted(set(a) & set(b))
    return dict(common_nodes=len(k), max_abs_diff=rnd(max(abs(a[j] - b[j]) for j in k), 7))


def cross_checks():
    """the same run drawn in two figures must give the same numbers (independent check of both calibrations);
    informational comparisons of the Appendix-B baseline runs are listed without interpretation"""
    out = []
    for c in ('GTB_blue_ni44', 'GTB_green_ni46', 'GTB_red_ni48'):
        for q in ('delta', 'eta'):
            out.append(dict(check=f'Fig.4 D0 {q} vs Fig.20b (after Watson) D0 {q}, {c}',
                            **_cmp(_dots(f'fig04_D0_{q}', c), _dots(f'fig20_D0afterW_{q}', c))))
    out.append(dict(check='Fig.8 vs Fig.9 blue markers |F^1_1|',
                    **_cmp(_dots('fig08_P1FF', 'GTB_blue_ni44'), _dots('fig09_eehadron', 'GTB_blue_ni44'))))
    x9, y9 = _line('fig09_eehadron', 'GTB_blue_ni44_fit_line')
    fit = {}
    for s8 in ('GTB_blue_ni44_BW_fit', 'GTB_blue_ni44_GS_fit', 'GTB_red_ni48_BW_fit', 'GTB_red_ni48_GS_fit'):
        x8, y8 = _line('fig08_P1FF', s8)
        m = (x8 > x9.min()) & (x8 < x9.max())
        d = np.interp(x8[m], x9, y9) - y8[m]
        fit[s8] = dict(rms=rnd(np.sqrt(np.mean(d * d)), 5), max_abs=rnd(np.max(np.abs(d)), 5))
    best = min(fit, key=lambda k: fit[k]['rms'])
    out.append(dict(check='Fig.9 blue "GTB" line vs the four Fig.8 fit curves (interpolated, common range)',
                    result=fit, closest=best))
    mpath = os.path.join(OUT, 'fig09_eehadron.metadata.json')
    with open(mpath) as f:
        m9 = json.load(f)
    for g in m9['groups']:
        if g['series'] == 'GTB_blue_ni44_fit_line':
            g['identified_as'] = f'Fig. 8 {best} (curve comparison: {fit})'
    with open(mpath, 'w') as f:
        json.dump(m9, f, indent=1)
    base = {'00': ('r00=5x10^-2', ['eps00=5x10^-8', 'eps00=1x10^-7', 'eps00=5x10^-7'], 'fig06_S0ff'),
            '11': ('r11=2', ['eps11=1x10^-6', 'eps11=6x10^-6', 'eps11=1x10^-5'], 'fig08_P1FF'),
            '02': ('r02=6', ['eps02=3x10^-6', 'eps02=5x10^-6', 'eps02=7x10^-6'], 'fig11_D0FF')}
    for w, (r_used, eps, ff) in base.items():
        n3 = _dots(f'fig19_nSR{w}', 'nSR=3')
        out.append(dict(check=f'informational: Fig.18 {r_used} (orange, "used") vs Fig.19 nSR=3, wave {w}',
                        **_cmp(_dots(f'fig18_Fasymp{w}A', r_used), n3)))
        for e in eps:
            out.append(dict(check=f'informational: Fig.19 {e} vs Fig.19 nSR=3, wave {w}',
                            **_cmp(_dots(f'fig19_SRerr{w}', e), n3)))
        for c in ('GTB_blue_ni44', 'GTB_green_ni46', 'GTB_red_ni48'):
            g = _dots(ff, c)
            if g:
                out.append(dict(check=f'informational: Fig.19 nSR=3 vs {ff} {c}, wave {w}', **_cmp(n3, g)))
    return out


def raster_note(fname):
    pdf = os.path.join(SRC, fname)
    if not os.path.exists(pdf):
        return dict(file_present=False)
    page = fitz.open(pdf)[0]
    return dict(file_present=True, embedded_images=len(page.get_images()), vector_paths=len(page.get_drawings()),
                text_chars=len(page.get_text().strip()))


# ----------------------------------------------------------------------------- main
def main():
    index = []
    for fname, cfg in INK.items():
        do_inkscape(fname, cfg, index)
    for cfg in MMA:
        do_mma(cfg, index)
    skipped = [
        dict(figure='7', file='rhopole.pdf', reason='raster image embedded in the PDF (3-D |S^1_1| on sheet II); '
                                                    'stale since v1 per inventory; not extracted',
             pdf_content=raster_note('rhopole.pdf')),
        dict(figure='17', file='psandetaconverge.png', reason='raster PNG only (4997x6631); not extracted'),
        dict(figure='3', file='shape.pdf', reason='already extracted: ../../gtb2403_prep/paper/fig3_shape_dots_extracted.txt'),
        dict(figure='1,2,12-16', file='schematics/diagrams', reason='no data'),
    ]
    with open(os.path.join(OUT, 'index.json'), 'w') as f:
        json.dump(dict(paper=PAPER, script=os.path.basename(__file__), script_sha256=sha256(__file__),
                       source_dir=SRC, mpi_gev=MPI_GEV, node_rule='M=50, nu0=-20',
                       nodes=[dict(j=int(j), s_mpi2=rnd(sv, 9), sqrt_s_gev=rnd(e, 9))
                              for j, sv, e in zip(J_NODES[:40], S_NODES[:40], SQRT_S_NODES_GEV[:40])],
                       panels=index, cross_checks=cross_checks(), skipped=skipped), f, indent=1)
    print(f'{len(index)} panels written to {OUT}')


if __name__ == '__main__':
    sys.exit(main())
