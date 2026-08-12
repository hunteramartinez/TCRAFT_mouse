#!/usr/bin/env python3
"""
Golden Gate assembly diagram for mouse TRAC TCR HDR.
Panel A — gene block anatomy (detailed segments)
Panel B — Step 1: BbsI 3-fragment assembly (schematic staircase)
Panel C — Step 2: BsmBI TRBC-T2A insertion

Usage:  .venv/bin/python3 docs/make_assembly_figure.py
Output: docs/assembly_figure.pdf / .png
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

OUT = Path(__file__).parent

# ── Palette ───────────────────────────────────────────────────────────────────
C = dict(
    lh='#6B7280', p2a='#84CC16', trbv='#3B82F6', cdr3b='#F97316',
    trbc='#0F766E', t2a='#A855F7', trav='#60A5FA', cdr3a='#F87171',
    trac='#059669', rh='#6B7280', bsmbi_site='#92400E',
    bbsi='#1E40AF',   # BbsI enzyme colour
    bsmbi='#B45309',  # BsmBI enzyme colour
    chip1='#BFDBFE',  # pale-blue chip  (BbsI OHs)
    chip2='#FDE68A',  # pale-amber chip (BsmBI OHs)
    txt='#111827', muted='#6B7280',
)
plt.rcParams.update({'font.family': 'DejaVu Sans',
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.spines.left': False, 'axes.spines.bottom': False})

OH = {'1a': 'CCAG', '1b': 'GGAA', '2a': 'GGAA', '2b': 'AGGA'}

# ── Segment definitions  (color_key, relative_nt, label) ─────────────────────
# Anatomy Panel – detailed sub-segments
ANA_LH   = [('lh',300,'LH arm (~300 bp)'),('p2a',66,'P2A'),('trbv',360,'TRBV16')]
ANA_VAR  = [('trbv',10,''),('cdr3b',30,'CDR3β'),('trbc',10,''),
            ('bsmbi_site',15,'2a'),('bsmbi_site',15,'2b'),
            ('t2a',10,''),('trav',330,'TRAV9N-3'),('cdr3a',24,'CDR3α')]
ANA_TRAC = [('trac',320,'TRAC-LVL-Cys'),('rh',634,'RH arm + BGH polyA')]
ANA_TRBC = [('trbc',450,'TRBC1-Cys'),('t2a',75,'T2A')]

# Assembly panels – simplified schematic bars (equal visual weight)
SCH_LH   = [('lh',1.2,'LH arm'),('p2a',0.3,'P2A'),('trbv',1.5,'TRBV16')]
SCH_VAR  = [('cdr3b',0.7,'CDR3β'),('bsmbi_site',0.5,'2a·2b'),
            ('trav',1.4,'TRAV9N-3'),('cdr3a',0.4,'CDR3α')]
SCH_TRAC = [('trac',1.2,'TRAC-LVL-Cys'),('rh',1.8,'RH arm + BGH')]
SCH_TRBC = [('trbc',1.5,'TRBC1-Cys'),('t2a',0.5,'T2A')]

FULL_PRODUCT = [
    ('lh',300,'LH'),('p2a',66,'P2A'),('trbv',360,'TRBV16'),
    ('cdr3b',30,'CDR3β'),('trbc',450,'TRBC1-Cys'),('t2a',75,'T2A'),
    ('trav',330,'TRAV9N-3'),('cdr3a',24,'CDR3α'),
    ('trac',320,'TRAC-LVL-Cys'),('rh',634,'BGH + RH'),
]
STEP1_PRODUCT = [
    ('lh',300,'LH'),('p2a',66,'P2A'),('trbv',360,'TRBV16'),
    ('cdr3b',30,'CDR3β'),('bsmbi_site',30,'2a·2b'),
    ('trav',330,'TRAV9N-3'),('cdr3a',24,'CDR3α'),
    ('trac',320,'TRAC'),('rh',634,'BGH + RH'),
]
S2_LEFT  = [('lh',1.0,'LH'),('p2a',0.2,'P2A'),('trbv',1.3,'TRBV16'),('cdr3b',0.5,'CDR3β')]
S2_RIGHT = [('trav',1.3,'TRAV9N-3'),('cdr3a',0.4,'CDR3α'),
            ('trac',1.0,'TRAC-LVL-Cys'),('rh',1.5,'BGH+RH')]


# ── Primitives ────────────────────────────────────────────────────────────────

def draw_bar(ax, segs, x0, x1, y, h, frag_label='', fontsize=7.0):
    tw  = x1 - x0
    tot = sum(s[1] for s in segs)
    cur = x0
    for ck, n, txt in segs:
        w  = (n / tot) * tw
        fc = C.get(ck, '#CBD5E1')
        ax.add_patch(mpatches.FancyBboxPatch(
            (cur, y), w, h, boxstyle='square,pad=0',
            facecolor=fc, edgecolor='white', linewidth=0.7,
            transform=ax.transAxes, zorder=2, clip_on=True))
        if txt and w > 0.025:
            fs = max(5.2, fontsize - max(0, len(txt)-8)*0.35)
            ax.text(cur+w/2, y+h/2, txt, ha='center', va='center',
                    fontsize=fs, color='white', fontweight='bold',
                    transform=ax.transAxes, zorder=3, clip_on=True)
        cur += w
    if frag_label:
        ax.text(x0-0.008, y+h/2, frag_label, ha='right', va='center',
                fontsize=8, fontweight='bold', color=C['txt'],
                transform=ax.transAxes)


def oh_chip(ax, x0, y, h, seq, bg, ec_color, side='right', fontsize=8.0):
    """Draw a sticky-end chip.  side='right' → chip extends rightward from x0."""
    CW = 0.048
    if side == 'left':
        x0 = x0 - CW
    th = h * 1.08
    ty = y - h * 0.04
    ax.add_patch(mpatches.FancyBboxPatch(
        (x0, ty), CW, th, boxstyle='round,pad=0.005',
        facecolor=bg, edgecolor=ec_color, linewidth=1.0,
        transform=ax.transAxes, zorder=5))
    ax.text(x0+CW/2, ty+th/2, seq, ha='center', va='center',
            fontsize=fontsize, fontweight='bold', fontfamily='monospace',
            color='#1E1B4B', transform=ax.transAxes, zorder=6)
    return x0 + CW/2   # return centre


def junction_label(ax, x, y, text, color, fontsize=6.8):
    ax.text(x, y, text, ha='center', va='top', fontsize=fontsize,
            color=color, fontweight='bold', transform=ax.transAxes)


def ligation_arrow(ax, x, y_from, y_to, label=''):
    ax.annotate('', xy=(x, y_to), xytext=(x, y_from),
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='#374151',
                                lw=1.7, mutation_scale=12))
    if label:
        ax.text(x+0.012, (y_from+y_to)/2, label, ha='left', va='center',
                fontsize=8, color='#374151', transform=ax.transAxes)


def dashed_connect(ax, x1, y1, x2, y2, color):
    """Step-line connecting two OH chips."""
    ym = (y1+y2)/2
    ax.plot([x1,x1,x2,x2],[y1,ym,ym,y2], '--', color=color,
            lw=1.0, alpha=0.7, transform=ax.transAxes, zorder=1)


def panel_label(ax, letter, title):
    ax.text(0.0, 0.99, letter, transform=ax.transAxes,
            fontsize=11, fontweight='bold', color=C['txt'], va='top')
    ax.text(0.028, 0.99, title, transform=ax.transAxes,
            fontsize=9.5, fontweight='bold', color=C['txt'], va='top')


# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 10.5), facecolor='white')
gs  = fig.add_gridspec(3, 1, hspace=0.55,
                        top=0.94, bottom=0.05, left=0.005, right=0.995)
axes = [fig.add_subplot(gs[i]) for i in range(3)]
for ax in axes:
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_aspect('auto'); ax.axis('off')

fig.text(0.5, 0.977,
         'Mouse TRAC TCR HDR — Two-Step Golden Gate Assembly'
         '  (2D2: TRBV16 / TRAV9N-3)',
         ha='center', va='top', fontsize=11.5, fontweight='bold', color=C['txt'])


# ══════════════════════════════════════════════════════════════════════════════
# Panel A — Gene Block Anatomy
# ══════════════════════════════════════════════════════════════════════════════
ax = axes[0]
panel_label(ax, 'A', '  Gene Block Anatomy')

H   = 0.175
TAB = 0.050   # chip width
X0, X1 = 0.18, 0.90   # bar span (X0 leaves room for fragment labels)

# Each row: (name, segs, y, left_oh, left_chip, right_oh, right_chip, enz_col)
rows = [
    ('LH-TRBV16',      ANA_LH,   0.77, None,       None,       'CCAG', C['chip1'], C['bbsi']),
    ('Variable (2D2)', ANA_VAR,  0.53, 'CCAG',     C['chip1'], 'GGAA', C['chip1'], C['bbsi']),
    ('TRAC-RH',        ANA_TRAC, 0.29, 'GGAA',     C['chip1'], None,   None,        C['bbsi']),
    ('TRBC-T2A',       ANA_TRBC, 0.05, 'GGAA',     C['chip2'], 'AGGA', C['chip2'], C['bsmbi']),
]

for (name, segs, y, ltab, lchip, rtab, rchip, ecol) in rows:
    # Adjust bar X0 if there's a left chip
    bx0 = X0 + TAB + 0.004 if ltab else X0
    draw_bar(ax, segs, bx0, X1, y, H, frag_label=name)
    if ltab:
        cx = oh_chip(ax, bx0-TAB-0.004, y, H, ltab, lchip, ecol, side='right')
        ax.text(cx, y-0.035, ltab, ha='center', va='top', fontsize=5.5,
                color=ecol, fontfamily='monospace', fontweight='bold',
                transform=ax.transAxes)
    if rtab:
        cx = oh_chip(ax, X1+0.004, y, H, rtab, rchip, ecol, side='right')
        ax.text(cx, y-0.035, rtab, ha='center', va='top', fontsize=5.5,
                color=ecol, fontfamily='monospace', fontweight='bold',
                transform=ax.transAxes)

# BsmBI site tick marks on Variable bar
var_y   = 0.53
bx0_var = X0 + TAB + 0.004
var_w   = X1 - bx0_var
frac_2a = (10+30+10) / (10+30+10+15+15+10+330+24)
frac_2b = (10+30+10+15) / (10+30+10+15+15+10+330+24)
for frac, lbl in [(frac_2a, f'2a:{OH["2a"]}'), (frac_2b, f'2b:{OH["2b"]}')]:
    xb = bx0_var + frac * var_w
    ax.plot([xb,xb], [var_y+H, var_y+H+0.04], '--', lw=0.8,
            color=C['bsmbi'], alpha=0.7, transform=ax.transAxes)
    ax.text(xb, var_y+H+0.045, lbl, ha='center', va='bottom', fontsize=5.8,
            color=C['bsmbi'], fontweight='bold', fontfamily='monospace',
            transform=ax.transAxes)

# Enzyme callout labels
ax.text(X1+TAB+0.012, 0.77+H+0.01, '← BbsI', ha='left', va='bottom',
        fontsize=6.5, color=C['bbsi'], transform=ax.transAxes)
ax.text(X1+TAB+0.012, 0.05+H+0.01, '← BsmBI', ha='left', va='bottom',
        fontsize=6.5, color=C['bsmbi'], transform=ax.transAxes)

# Legend
legend = [('LH/RH arm',C['lh']),('P2A leader',C['p2a']),('TRBV gene',C['trbv']),
          ('CDR3β',C['cdr3b']),('TRBC1-Cys',C['trbc']),('T2A linker',C['t2a']),
          ('TRAV gene',C['trav']),('CDR3α',C['cdr3a']),
          ('TRAC-LVL-Cys',C['trac']),('Restriction site',C['bsmbi_site'])]
lx, ly = 0.18, -0.03
for i,(lbl,fc) in enumerate(legend):
    xi = lx + i*0.082
    ax.add_patch(mpatches.FancyBboxPatch((xi,ly),0.014,0.055,
        boxstyle='square,pad=0', facecolor=fc, edgecolor='white',
        linewidth=0.3, transform=ax.transAxes))
    ax.text(xi+0.017, ly+0.028, lbl, ha='left', va='center',
            fontsize=5.6, color=C['txt'], transform=ax.transAxes)


# ══════════════════════════════════════════════════════════════════════════════
# Panel B — Step 1: BbsI Assembly  (staircase, schematic widths)
# ══════════════════════════════════════════════════════════════════════════════
ax = axes[1]
panel_label(ax, 'B',
    '  Step 1 — BbsI Digest + Ligation  '
    '( LH-TRBV16 + Variable + TRAC-RH → ~2.1 kb )')

H    = 0.165
TAB  = 0.052
GAP  = 0.005

# Layout (schematic, non-proportional):
#  LH-TRBV16  [─────────────────] [CCAG|
#                           |CCAG] [──Variable──] [GGAA|
#                                          |GGAA] [──TRAC-RH──────]
#
# x positions chosen so matching chips are ADJACENT:
J1x = 0.46    # junction 1 (oh_1a) x-centre  → LH right chip 0.46..0.512, Var left chip 0.512..0.564
J2x = 0.72    # junction 2 (oh_1b) x-centre  → Var right chip 0.72..0.772, TRAC left 0.772..0.824

LH_x0, LH_x1   = 0.06,          J1x
VAR_x0, VAR_x1 = J1x + TAB + GAP, J2x
TRAC_x0, TRAC_x1 = J2x + TAB + GAP, 0.97

y_lh   = 0.80
y_var  = 0.58
y_trac = 0.36

draw_bar(ax, SCH_LH,   LH_x0,   LH_x1,   y_lh,   H, frag_label='LH-TRBV16')
draw_bar(ax, SCH_VAR,  VAR_x0,  VAR_x1,  y_var,  H)
draw_bar(ax, SCH_TRAC, TRAC_x0, TRAC_x1, y_trac, H, frag_label='TRAC-RH')
ax.text((VAR_x0+VAR_x1)/2, y_var-0.03, 'Variable (2D2)',
        ha='center', va='top', fontsize=8, fontweight='bold',
        color=C['txt'], transform=ax.transAxes)

# OH chips at junction 1 (oh_1a = CCAG, BbsI)
c1r = oh_chip(ax, LH_x1  + GAP, y_lh,  H, OH['1a'], C['chip1'], C['bbsi'])  # LH right
c1l = oh_chip(ax, VAR_x0 - TAB, y_var, H, OH['1a'], C['chip1'], C['bbsi'])  # Var left
dashed_connect(ax, c1r, y_lh, c1l, y_var+H, C['bbsi'])
ax.text((c1r+c1l)/2+0.02, (y_lh+y_var+H)/2,
        f'oh_1a = {OH["1a"]}', ha='left', va='center',
        fontsize=7.5, color=C['bbsi'], fontweight='bold', transform=ax.transAxes)

# OH chips at junction 2 (oh_1b = GGAA, BbsI)
c2r = oh_chip(ax, VAR_x1  + GAP, y_var,  H, OH['1b'], C['chip1'], C['bbsi'])  # Var right
c2l = oh_chip(ax, TRAC_x0 - TAB, y_trac, H, OH['1b'], C['chip1'], C['bbsi'])  # TRAC left
dashed_connect(ax, c2r, y_var, c2l, y_trac+H, C['bbsi'])
ax.text((c2r+c2l)/2+0.02, (y_var+y_trac+H)/2,
        f'oh_1b = {OH["1b"]}', ha='left', va='center',
        fontsize=7.5, color=C['bbsi'], fontweight='bold', transform=ax.transAxes)

# Internal BsmBI site marks on Variable bar (2a, 2b)
vw = VAR_x1 - VAR_x0
tot_var = sum(s[1] for s in SCH_VAR)
x2a = VAR_x0 + (0.7/tot_var)*vw
x2b = VAR_x0 + (0.7+0.5)/tot_var*vw
for xb, lbl in [(x2a, f'2a:{OH["2a"]}'), (x2b, f'2b:{OH["2b"]}')]:
    ax.plot([xb,xb],[y_var+H, y_var+H+0.04], '--', lw=0.8,
            color=C['bsmbi'], alpha=0.7, transform=ax.transAxes)
    ax.text(xb, y_var+H+0.045, lbl, ha='center', va='bottom',
            fontsize=6, color=C['bsmbi'], fontweight='bold',
            fontfamily='monospace', transform=ax.transAxes)

# Ligation arrow + Step 1 product
ligation_arrow(ax, 0.50, y_trac-0.02, y_trac-0.14, 'BbsI + T4 Ligase')
draw_bar(ax, STEP1_PRODUCT, 0.06, 0.97, 0.04, H,
         frag_label='Step 1\nproduct\n~2.1 kb')

# Mark BsmBI sites on product
pw   = 0.91
ptot = sum(s[1] for s in STEP1_PRODUCT)
poff = 0.06
cum  = 0
for ck, n, _ in STEP1_PRODUCT:
    if ck == 'bsmbi_site':
        xm = poff + (cum + n/2)/ptot * pw
        ax.plot([xm,xm],[0.04,0.04+H], '-', lw=1.3,
                color=C['bsmbi'], alpha=0.8, transform=ax.transAxes, zorder=4)
        ax.text(xm, 0.04+H+0.01, f'2a·2b\n{OH["2a"]}·{OH["2b"]}',
                ha='center', va='bottom', fontsize=5.8, color=C['bsmbi'],
                fontweight='bold', transform=ax.transAxes)
    cum += n


# ══════════════════════════════════════════════════════════════════════════════
# Panel C — Step 2: BsmBI Assembly  (insertion diagram)
# ══════════════════════════════════════════════════════════════════════════════
ax = axes[2]
panel_label(ax, 'C',
    '  Step 2 — BsmBI Digest + Ligation  '
    '( Step 1 product + TRBC-T2A → ~2.5 kb )')

H   = 0.165
TAB = 0.052
GAP = 0.005

# Step 1 product shown pre-cut at 2a and 2b:
#   LEFT  piece [LH─P2A─TRBV─CDR3β] [GGAA chip] | gap | [AGGA chip] [TRAV─CDR3α─TRAC─RH]
#   TRBC-T2A fragment              [GGAA chip] [TRBC-T2A bar] [AGGA chip]

J2a = 0.40   # 2a junction x
J2b = 0.64   # 2b junction x

LEFT_x0,  LEFT_x1  = 0.06,          J2a
RIGHT_x0, RIGHT_x1 = J2b + TAB + GAP, 0.97
TRBC_x0,  TRBC_x1  = J2a + TAB + GAP, J2b

y_prod  = 0.80
y_trbc  = 0.52

draw_bar(ax, S2_LEFT,  LEFT_x0,  LEFT_x1,  y_prod, H,
         frag_label='Step 1\nproduct\n(BsmBI cut)')
draw_bar(ax, S2_RIGHT, RIGHT_x0, RIGHT_x1, y_prod, H)
draw_bar(ax, SCH_TRBC, TRBC_x0,  TRBC_x1,  y_trbc, H,
         frag_label='TRBC-T2A')

# Bracket showing BsmBI-cleaved gap on the Step 1 product row
ax.annotate('', xy=(J2b+TAB, y_prod+H/2), xytext=(J2a, y_prod+H/2),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='<->', color=C['bsmbi'], lw=1.0))
ax.text((J2a + J2b+TAB)/2, y_prod+H/2+0.025,
        'BsmBI cuts\n(2a + 2b)', ha='center', va='bottom',
        fontsize=6.5, color=C['bsmbi'], transform=ax.transAxes)

# OH chips — 2a junction (GGAA)
c2aL = oh_chip(ax, LEFT_x1+GAP,   y_prod, H, OH['2a'], C['chip2'], C['bsmbi'])  # product left piece right
c2aR = oh_chip(ax, TRBC_x0-TAB,   y_trbc, H, OH['2a'], C['chip2'], C['bsmbi'])  # TRBC-T2A left
dashed_connect(ax, c2aL, y_prod, c2aR, y_trbc+H, C['bsmbi'])
ax.text(c2aL+0.015, (y_prod+y_trbc+H)/2+0.02,
        f'oh_2a = {OH["2a"]}', ha='left', va='center',
        fontsize=7.5, color=C['bsmbi'], fontweight='bold', transform=ax.transAxes)

# OH chips — 2b junction (AGGA)
c2bL = oh_chip(ax, TRBC_x1+GAP,   y_trbc, H, OH['2b'], C['chip2'], C['bsmbi'])  # TRBC-T2A right
c2bR = oh_chip(ax, RIGHT_x0-TAB,  y_prod, H, OH['2b'], C['chip2'], C['bsmbi'])  # product right piece left
dashed_connect(ax, c2bL, y_trbc+H, c2bR, y_prod, C['bsmbi'])
ax.text(c2bL+0.015, (y_prod+y_trbc+H)/2+0.02,
        f'oh_2b = {OH["2b"]}', ha='left', va='center',
        fontsize=7.5, color=C['bsmbi'], fontweight='bold', transform=ax.transAxes)

# Ligation arrow + final product
ligation_arrow(ax, 0.50, y_trbc-0.02, y_trbc-0.14, 'BsmBI + T4 Ligase')
draw_bar(ax, FULL_PRODUCT, 0.06, 0.97, 0.05, H,
         frag_label='Assembled\ntemplate\n~2.5 kb')
ax.text(0.975, 0.05+H/2, '→ PCR with\n  HDR_F / HDR_R',
        ha='left', va='center', fontsize=7, color='#374151',
        transform=ax.transAxes)

# Footer
fig.text(0.5, 0.015,
         f'BbsI OHs (Step 1):   oh_1a = {OH["1a"]}  ·  oh_1b = {OH["1b"]}'
         f'               '
         f'BsmBI OHs (Step 2):  oh_2a = {OH["2a"]}  ·  oh_2b = {OH["2b"]}',
         ha='center', va='bottom', fontsize=8, color='#374151',
         fontfamily='monospace')

# ── Save ──────────────────────────────────────────────────────────────────────
for ext in ('pdf', 'png'):
    fig.savefig(OUT / f'assembly_figure.{ext}',
                dpi=300 if ext == 'png' else 150,
                bbox_inches='tight', facecolor='white')
    print(f'Saved → {OUT}/assembly_figure.{ext}')
plt.close(fig)
