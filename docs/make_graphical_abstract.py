#!/usr/bin/env python3
"""
Graphical abstract — the two mouse TCR construction methods in TCRAFT.

Both pipelines consume the same AIRR input and share references, codon optimisation and
FR4 handling; they diverge only in how the finished cassette is delivered.  The figure is
laid out as a decision guide: shared inputs on top, the fork in the middle, and the
criteria for choosing each route at the bottom.

Usage:  uv run python docs/make_graphical_abstract.py
Output: docs/graphical_abstract.pdf / .png
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUT = Path(__file__).parent

# Shared with make_library_figure.py / make_assembly_figure.py — keep in visual sync.
C = dict(
    trbv='#3B82F6', cdr3b='#F97316', trbc='#0F766E', t2a='#A855F7',
    trav='#60A5FA', cdr3a='#F87171', trac='#059669', grey='#6B7280',
    txt='#1F2937', muted='#64748B',
)
HDR = '#B45309'      # knock-in route accent
LIB = '#1D4ED8'      # plasmid-library route accent
SHARED = '#0F766E'   # shared-input accent


def box(ax, x0, y0, w, h, fc, ec, lw=1.3, r=0.012, z=2, alpha=1.0):
    ax.add_patch(mpatches.FancyBboxPatch(
        (x0, y0), w, h, boxstyle=f'round,pad=0,rounding_size={r}',
        facecolor=fc, edgecolor=ec, linewidth=lw, alpha=alpha,
        transform=ax.transAxes, zorder=z, clip_on=False))


def text(ax, x, y, s, size=9, color=None, weight='normal', ha='center', va='center',
         style='normal', family=None, z=5):
    ax.text(x, y, s, fontsize=size, color=color or C['txt'], fontweight=weight,
            ha=ha, va=va, style=style, fontfamily=family, transform=ax.transAxes, zorder=z)


def arrow(ax, x0, y0, x1, y1, color, lw=2.0, ls='-'):
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                xycoords='axes fraction', textcoords='axes fraction', zorder=4,
                arrowprops=dict(arrowstyle='-|>', color=color, lw=lw, linestyle=ls,
                                mutation_scale=16, shrinkA=2, shrinkB=2))


def cassette(ax, x0, x1, y, h, segs, fontsize=5.6):
    """Small coloured cassette bar; segs = [(palette_key, weight, label)]."""
    tot, cur = sum(s[1] for s in segs), x0
    for ck, n, lab in segs:
        w = (n / tot) * (x1 - x0)
        ax.add_patch(mpatches.Rectangle((cur, y), w, h, facecolor=C.get(ck, '#CBD5E1'),
                                        edgecolor='white', linewidth=0.6,
                                        transform=ax.transAxes, zorder=3, clip_on=False))
        if lab and w > 0.030:
            ax.text(cur + w / 2, y + h / 2, lab, ha='center', va='center',
                    fontsize=fontsize, color='white', fontweight='bold',
                    transform=ax.transAxes, zorder=4)
        cur += w


CASSETTE = [('trbv', 16, 'TRBV'), ('cdr3b', 10, 'CDR3β'), ('trbc', 17, 'TRBC'),
            ('t2a', 7, 'T2A'), ('trav', 21, 'TRAV'), ('cdr3a', 10, 'CDR3α'),
            ('trac', 19, 'TRAC')]


def bullets(ax, x, y_top, items, color, size=8.0, dy=0.030):
    for i, s in enumerate(items):
        text(ax, x, y_top - i * dy, '•', size=size, color=color, ha='left')
        text(ax, x + 0.016, y_top - i * dy, s, size=size, ha='left', color=C['txt'])


def build(ax):
    # ── Shared input ────────────────────────────────────────────────────────
    box(ax, 0.255, 0.890, 0.49, 0.062, '#ECFDF5', SHARED, lw=1.6)
    text(ax, 0.50, 0.932, 'AIRR rearrangement file  (e.g. Cell Ranger VDJ)',
         size=10.5, weight='bold', color=SHARED)
    text(ax, 0.50, 0.907, 'paired α/β clonotypes  ·  dual-α and dual-β resolved by read support',
         size=8.0, color=C['muted'])

    text(ax, 0.50, 0.862,
         'shared:  mouse IMGT references  ·  codon optimisation  ·  J-segment FR4  ·  '
         'cysteine-modified constant regions',
         size=8.2, color=C['muted'], style='italic')

    # one TCR cassette, built the same way for both routes
    cassette(ax, 0.295, 0.705, 0.805, 0.034, CASSETTE)
    text(ax, 0.50, 0.787, 'one TCR cassette, identical in both routes',
         size=8.0, color=C['muted'])

    # ── Fork ────────────────────────────────────────────────────────────────
    arrow(ax, 0.47, 0.775, 0.25, 0.735, HDR)
    arrow(ax, 0.53, 0.775, 0.75, 0.735, LIB)

    for x0, accent, name, sub in (
            (0.035, HDR, 'TCRAFT-hdr', 'knock-in at the endogenous locus'),
            (0.525, LIB, 'TCRAFT-library', 'pooled retroviral plasmid library')):
        box(ax, x0, 0.672, 0.44, 0.058, '#FFFFFF', accent, lw=1.8)
        text(ax, x0 + 0.22, 0.712, name, size=11.5, weight='bold', color=accent,
             family='monospace')
        text(ax, x0 + 0.22, 0.689, sub, size=8.4, color=C['muted'])

    # ── Parts ordered ───────────────────────────────────────────────────────
    text(ax, 0.255, 0.645, 'DNA to order', size=8.6, weight='bold', color=C['muted'])
    text(ax, 0.745, 0.645, 'DNA to order', size=8.6, weight='bold', color=C['muted'])

    box(ax, 0.035, 0.505, 0.44, 0.128, '#FFFBEB', '#FCD34D', lw=1.0)
    bullets(ax, 0.055, 0.606, [
        'N  LH–TRBV gene blocks (one per TRBV)',
        '1  TRAC–RH block  ·  1  TRBC–T2A block',
        'N  variable fragments, one per clonotype',
        'all linear, synthesised fresh',
    ], HDR)

    box(ax, 0.525, 0.505, 0.44, 0.128, '#EFF6FF', '#93C5FD', lw=1.0)
    bullets(ax, 0.545, 0.606, [
        '24  TRBV donor plasmids  —  built once, reused',
        '1  TRAC backbone  ·  1  TRBC–T2A donor',
        'N  variable fragments, one per clonotype',
        'vector domesticated once, no PCR',
    ], LIB)

    # ── Assembly ────────────────────────────────────────────────────────────
    arrow(ax, 0.255, 0.497, 0.255, 0.455, HDR)
    arrow(ax, 0.745, 0.497, 0.745, 0.455, LIB)

    box(ax, 0.035, 0.378, 0.44, 0.072, '#FFFFFF', '#E5E7EB', lw=1.0)
    text(ax, 0.255, 0.428, 'two-step Golden Gate', size=9.0, weight='bold')
    text(ax, 0.255, 0.400, 'BbsI  →  BsmBI   ·   all TRBV in one pool',
         size=8.2, color=C['muted'], family='monospace')

    box(ax, 0.525, 0.378, 0.44, 0.072, '#FFFFFF', '#E5E7EB', lw=1.0)
    text(ax, 0.745, 0.428, 'two-step Golden Gate', size=9.0, weight='bold')
    text(ax, 0.745, 0.400, 'BbsI  →  BsmBI   ·   one reaction per TRBV pool',
         size=8.2, color=C['muted'], family='monospace')

    # ── Product ─────────────────────────────────────────────────────────────
    arrow(ax, 0.255, 0.370, 0.255, 0.330, HDR)
    arrow(ax, 0.745, 0.370, 0.745, 0.330, LIB)

    text(ax, 0.255, 0.312, 'linear HDR template', size=9.0, weight='bold', color=HDR)
    cassette(ax, 0.085, 0.425, 0.258, 0.030,
             [('grey', 13, 'LH'), *CASSETTE, ('grey', 13, 'RH')], fontsize=5.0)

    text(ax, 0.745, 0.312, 'circular retroviral plasmid pool', size=9.0, weight='bold',
         color=LIB)
    ax.add_patch(mpatches.Circle((0.745, 0.247), 0.040, facecolor='none',
                                 edgecolor=LIB, linewidth=5.0,
                                 transform=ax.transAxes, zorder=3))
    text(ax, 0.745, 0.247, 'TCR\nIRES\nGFP', size=5.8, color=LIB, weight='bold')

    # ── Delivery ────────────────────────────────────────────────────────────
    arrow(ax, 0.255, 0.240, 0.255, 0.206, HDR)
    arrow(ax, 0.745, 0.203, 0.745, 0.206, LIB)

    text(ax, 0.255, 0.188, 'Cas9-RNP electroporation  →  TRAC or TRBC1 locus',
         size=8.6, weight='bold')
    text(ax, 0.745, 0.188, 'retroviral packaging  →  transduction',
         size=8.6, weight='bold')

    # ── Decision criteria ───────────────────────────────────────────────────
    box(ax, 0.035, 0.022, 0.44, 0.146, '#FFF7ED', HDR, lw=1.5)
    text(ax, 0.255, 0.151, 'choose this when you need', size=8.8, weight='bold', color=HDR)
    bullets(ax, 0.055, 0.119, [
        'expression under the endogenous promoter',
        'single-copy, physiological TCR levels',
        'the native locus intact — in vivo work',
        'mature T cells (TRAC) or DN3 thymocytes (TRBC1)',
    ], HDR, size=7.8, dy=0.026)

    box(ax, 0.525, 0.022, 0.44, 0.146, '#EFF6FF', LIB, lw=1.5)
    text(ax, 0.745, 0.151, 'choose this when you need', size=8.8, weight='bold', color=LIB)
    bullets(ax, 0.545, 0.119, [
        'library scale — thousands of TCRs in a pool',
        'pooled screening with a selectable readout',
        'plasmids reused across experiments',
        'cell lines or primary cells, no Cas9 needed',
    ], LIB, size=7.8, dy=0.026)


def main():
    fig = plt.figure(figsize=(12.0, 9.0))
    ax = fig.add_axes([0.0, 0.0, 1.0, 0.955])
    ax.set_xlim(0, 1), ax.set_ylim(0, 1), ax.axis('off')
    build(ax)
    fig.suptitle('TCRAFT — two routes from a sequenced repertoire to expressed mouse TCRs',
                 fontsize=14, fontweight='bold', y=0.975)
    for ext in ('pdf', 'png'):
        fig.savefig(OUT / f'graphical_abstract.{ext}', dpi=200,
                    bbox_inches='tight', facecolor='white')
    print(f'wrote {OUT}/graphical_abstract.pdf and .png')


if __name__ == '__main__':
    main()
