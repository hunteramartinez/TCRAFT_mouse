#!/usr/bin/env python3
"""
Assembly diagram for the pooled TCR plasmid library (TCRAFT-library).

Panel A — the four components and what is a plasmid vs. synthesised per clonotype
Panel B — why the assembly must be two steps (the circular topology constraint)
Panel C — step 1 (BbsI, 3-piece) and step 2 (BsmBI, insert TRBC-T2A)

Reuses the palette and drawing helpers from make_assembly_figure.py so the two
figures read as one set.

Usage:  uv run python docs/make_library_figure.py
Output: docs/library_figure.pdf / .png
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUT = Path(__file__).parent

# Palette and helpers are kept local rather than imported from make_assembly_figure.py:
# that script is flat top-level code with no __main__ guard, so importing it would rebuild
# and overwrite assembly_figure.pdf/.png as a side effect.  Keep the two in visual sync.
C = dict(
    lh='#6B7280', p2a='#84CC16', trbv='#3B82F6', cdr3b='#F97316',
    trbc='#0F766E', t2a='#A855F7', trav='#60A5FA', cdr3a='#F87171',
    trac='#059669', rh='#6B7280', bsmbi_site='#92400E', bbsi='#1E40AF',
    txt='#1F2937',
)


def draw_bar(ax, segs, x0, x1, y, h, frag_label='', fontsize=7.0):
    tw, tot, cur = x1 - x0, sum(s[1] for s in segs), x0
    for ck, n, txt in segs:
        w = (n / tot) * tw
        ax.add_patch(mpatches.FancyBboxPatch(
            (cur, y), w, h, boxstyle='square,pad=0', facecolor=C.get(ck, '#CBD5E1'),
            edgecolor='white', linewidth=0.7, transform=ax.transAxes, zorder=2, clip_on=True))
        if txt and w > 0.025:
            fs = max(5.2, fontsize - max(0, len(txt) - 8) * 0.35)
            ax.text(cur + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=fs,
                    color='white', fontweight='bold', transform=ax.transAxes,
                    zorder=3, clip_on=True)
        cur += w
    if frag_label:
        ax.text(x0 - 0.008, y + h / 2, frag_label, ha='right', va='center', fontsize=8,
                fontweight='bold', color=C['txt'], transform=ax.transAxes)


def oh_chip(ax, x0, y, h, seq, bg, ec_color, fontsize=8.0):
    CW, th, ty = 0.048, h * 1.08, y - h * 0.04
    ax.add_patch(mpatches.FancyBboxPatch(
        (x0, ty), CW, th, boxstyle='round,pad=0.005', facecolor=bg,
        edgecolor=ec_color, linewidth=1.0, transform=ax.transAxes, zorder=5))
    ax.text(x0 + CW / 2, ty + th / 2, seq, ha='center', va='center', fontsize=fontsize,
            fontweight='bold', fontfamily='monospace', color='#1E1B4B',
            transform=ax.transAxes, zorder=6)
    return x0 + CW / 2


def junction_label(ax, x, y, text, color, fontsize=6.8):
    ax.text(x, y, text, ha='center', va='top', fontsize=fontsize, color=color,
            fontweight='bold', transform=ax.transAxes)


def ligation_arrow(ax, x, y_from, y_to, label=''):
    ax.annotate('', xy=(x, y_to), xytext=(x, y_from), xycoords='axes fraction',
                textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='#374151', lw=1.7, mutation_scale=12))
    if label:
        ax.text(x + 0.012, (y_from + y_to) / 2, label, ha='left', va='center',
                fontsize=8, color='#374151', transform=ax.transAxes)


def panel_label(ax, letter, title):
    ax.text(0.0, 0.99, letter, transform=ax.transAxes, fontsize=13, fontweight='bold',
            va='bottom', color=C['txt'])
    ax.text(0.028, 0.99, title, transform=ax.transAxes, fontsize=9.6, va='bottom',
            color=C['txt'])

OH = dict(v1='CACC', b='TCCA', a2='AAGA', b2='AGAA')
PLASMID = '#DBEAFE'      # fill for "comes from a plasmid"
SYNTH = '#FEF3C7'        # fill for "synthesised per clonotype"


def _legend_chip(ax, x, y, color, label, w=0.022, h=0.020):
    ax.add_patch(mpatches.FancyBboxPatch(
        (x, y), w, h, boxstyle='round,pad=0.003', facecolor=color,
        edgecolor='#64748B', linewidth=0.8, transform=ax.transAxes, zorder=4))
    ax.text(x + w + 0.008, y + h / 2, label, ha='left', va='center',
            fontsize=7.6, color=C['txt'], transform=ax.transAxes)


def panel_a(ax):
    panel_label(ax, 'A', 'Components — 26 plasmids + one synthesised fragment per clonotype')
    rows = [
        (0.80, 'TRAC backbone  ×1',
         [('lh', 14, 'Kozak'), ('rh', 26, 'ccdB stuffer'), ('trac', 24, 'TRAC'),
          ('lh', 10, 'stop'), ('t2a', 26, 'IRES–GFP'), ('rh', 22, 'LTR / ori / AmpR')]),
        (0.60, 'TRBV donors  ×24',
         [('trbv', 100, 'TRBV (codon-optimised, ATG → Cys104)')]),
        (0.40, 'TRBC-T2A donor  ×1',
         [('trbc', 72, 'TRBC (Cys)'), ('t2a', 28, 'furin-GSG-T2A')]),
        (0.18, 'synthesised  ×N',
         [('cdr3b', 20, 'CDR3β+FR4'), ('trbc', 7, ''), ('t2a', 7, ''),
          ('trav', 46, 'TRAV (full, with leader)'), ('cdr3a', 20, 'CDR3α+FR4')]),
    ]
    for y, label, segs in rows:
        draw_bar(ax, segs, 0.31, 0.97, y, 0.085, frag_label=label)

    # Provenance bracket rather than a colour legend: the bar colours already encode which
    # TCR part each segment is, so a second colour key would contradict them.
    ax.plot([0.105, 0.095, 0.095, 0.105], [0.843, 0.843, 0.443, 0.443],
            color='#64748B', lw=1.1, transform=ax.transAxes, clip_on=False)
    ax.text(0.087, 0.643, 'plasmids\n(built once)', ha='right', va='center',
            fontsize=7.6, color='#475569', transform=ax.transAxes)
    ax.plot([0.105, 0.095, 0.095, 0.105], [0.308, 0.308, 0.218, 0.218],
            color='#92400E', lw=1.1, transform=ax.transAxes, clip_on=False)
    ax.text(0.087, 0.263, 'synthesised\nper clonotype', ha='right', va='center',
            fontsize=7.6, color='#92400E', fontweight='bold', transform=ax.transAxes)
    ax.text(0.64, 0.105, 'mouse has 119 TRAV genes vs 44 in human — too many to hold as '
            'plasmids, so TRAV moves into synthesis',
            ha='center', va='top', fontsize=7.6, style='italic', color='#475569',
            transform=ax.transAxes)


def panel_b(ax):
    panel_label(ax, 'B', 'Why two steps: the variable regions are never adjacent')
    cx, cy, r = 0.30, 0.47, 0.30
    arcs = [('TRBV', 90, 140, 'trbv'), ('CDR3β+FR4', 140, 175, 'cdr3b'),
            ('TRBC·T2A', 175, 240, 'trbc'), ('TRAV·CDR3α', 240, 310, 'trav'),
            ('TRAC', 310, 360, 'trac'), ('vector', 0, 90, 'rh')]
    for label, a0, a1, ck in arcs:
        ax.add_patch(mpatches.Wedge((cx, cy), r, a0, a1, width=0.085,
                                    facecolor=C.get(ck, '#CBD5E1'), edgecolor='white',
                                    linewidth=1.2, transform=ax.transAxes, zorder=2))
        mid = np.radians((a0 + a1) / 2)
        # Wedge and labels share the same axes-fraction space, so one radius serves both;
        # scaling y separately here would pull every label off its arc.
        rl = r + 0.055
        ax.text(cx + rl * np.cos(mid), cy + rl * np.sin(mid), label,
                ha='center', va='center', fontsize=7.2, fontweight='bold',
                color=C['txt'], transform=ax.transAxes)

    ax.text(0.68, 0.72,
            'A TCR ORF alternates V–CDR3–C, V–CDR3–C, so the two\n'
            'variable regions are always separated by a constant one.',
            ha='left', va='top', fontsize=8.2, color=C['txt'], transform=ax.transAxes)
    ax.text(0.68, 0.50,
            'The synthesised fragment must stay ONE molecule —\n'
            'cut it open in a pool and α/β pairing scrambles.\n'
            'But one molecule cannot occupy two non-adjacent arcs.',
            ha='left', va='top', fontsize=8.2, color='#92400E', transform=ax.transAxes)
    ax.text(0.68, 0.26,
            '→ Step 1 closes the circle with the fragment intact,\n'
            '   locking pairing.  Step 2 then cuts the now-tethered\n'
            '   middle and inserts TRBC-T2A.',
            ha='left', va='top', fontsize=8.2, fontweight='bold',
            color='#166534', transform=ax.transAxes)


def panel_c(ax):
    panel_label(ax, 'C', 'Step 1 — BbsI, 3-piece   ·   Step 2 — BsmBI, insert TRBC-T2A')
    y1, h = 0.74, 0.075
    draw_bar(ax, [('rh', 100, '① TRAC backbone (ccdB dropped out)')], 0.08, 0.34, y1, h)
    draw_bar(ax, [('trbv', 100, '② TRBV donor')], 0.392, 0.57, y1, h)
    draw_bar(ax, [('cdr3b', 26, 'CDR3β'), ('trbc', 8, ''), ('t2a', 8, ''),
                  ('trav', 38, 'TRAV'), ('cdr3a', 20, 'CDR3α')],
             0.622, 0.90, y1, h)
    ax.text(0.761, y1 + h + 0.03, '④ synthesised fragment', ha='center', va='bottom',
            fontsize=7.8, fontweight='bold', color='#92400E', transform=ax.transAxes)
    for x, oh in ((0.34, OH['v1']), (0.57, 'oh_1a'), (0.90, OH['b'])):
        oh_chip(ax, x, y1, h, oh, '#E0E7FF', C['bbsi'])
    junction_label(ax, 0.364, y1 - 0.022, 'vector→TRBV', C['bbsi'])
    junction_label(ax, 0.594, y1 - 0.022, 'TRBV→CDR3β', C['bbsi'])
    junction_label(ax, 0.924, y1 - 0.022, 'CDR3α→TRAC', C['bbsi'])

    ligation_arrow(ax, 0.52, y1 - 0.07, 0.50, 'BbsI  →  circle; α/β pairing locked')

    y2 = 0.355
    draw_bar(ax, [('trbv', 17, 'TRBV'), ('cdr3b', 13, 'CDR3β')], 0.10, 0.352, y2, h,
             frag_label='step 1')
    draw_bar(ax, [('trav', 26, 'TRAV'), ('cdr3a', 13, 'CDR3α'), ('trac', 15, 'TRAC'),
                  ('rh', 13, 'vector')], 0.502, 0.95, y2, h)
    for x, oh in ((0.352, OH['a2']), (0.454, OH['b2'])):
        oh_chip(ax, x, y2, h, oh, '#FEF3C7', C['bsmbi_site'])
    junction_label(ax, 0.403, y2 - 0.022, 'BsmBI opens the tethered middle',
                   C['bsmbi_site'])

    ligation_arrow(ax, 0.52, y2 - 0.07, 0.165, 'BsmBI  +  ③ TRBC-T2A')

    y3 = 0.045
    draw_bar(ax, [('trbv', 15, 'TRBV'), ('cdr3b', 10, 'CDR3β'), ('trbc', 17, 'TRBC'),
                  ('t2a', 7, 'T2A'), ('trav', 21, 'TRAV'), ('cdr3a', 10, 'CDR3α'),
                  ('trac', 12, 'TRAC'), ('rh', 8, 'vector')],
             0.10, 0.95, y3, h, frag_label='final')


def main():
    fig = plt.figure(figsize=(11.5, 13.0))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.05, 1.35], hspace=0.18,
                          left=0.045, right=0.985, top=0.955, bottom=0.025)
    for i, fn in enumerate((panel_a, panel_b, panel_c)):
        ax = fig.add_subplot(gs[i])
        ax.set_xlim(0, 1), ax.set_ylim(0, 1), ax.axis('off')
        fn(ax)
    fig.suptitle('TCRAFT-library — pooled mouse TCR plasmid library assembly',
                 fontsize=13, fontweight='bold', y=0.985)
    for ext in ('pdf', 'png'):
        fig.savefig(OUT / f'library_figure.{ext}', dpi=200,
                    bbox_inches='tight', facecolor='white')
    print(f'wrote {OUT}/library_figure.pdf and .png')


if __name__ == '__main__':
    main()
