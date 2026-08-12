#!/usr/bin/env python3
"""
Build a great_tables overhang-alignment table showing the exact junction
sequences at each golden gate cut site for the 2D2 clone.

Source: assembled.gb (no restriction-site artefacts) + individual gene blocks.
Output: docs/assembly_table.html

Usage: .venv/bin/python3 docs/make_assembly_table.py
"""
from pathlib import Path
from Bio import SeqIO
from Bio.Restriction import BbsI, BsmBI, RestrictionBatch, Analysis
from Bio.Seq import Seq
import great_tables as gt
import polars as pl

OUT  = Path(__file__).parent
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
rc   = lambda s: s.translate(COMP)[::-1]
CTX  = 14   # flanking nucleotides shown on each side of the overhang

# ── Load sequences ─────────────────────────────────────────────────────────────
gb_dir   = Path('custom_tcrs/gb')
asm_gb   = gb_dir / 'clonotypes/2D2/assembled.gb'   # final product, no GG sites
var_gb   = gb_dir / 'clonotypes/2D2/variable.gb'
trbc_gb  = gb_dir / 'constants/TRBC-T2A.gb'

asm_rec  = SeqIO.read(asm_gb,  'genbank')
var_rec  = SeqIO.read(var_gb,  'genbank')
trbc_rec = SeqIO.read(trbc_gb, 'genbank')

asm_seq  = str(asm_rec.seq).upper()
var_seq  = str(var_rec.seq).upper()
trbc_seq = str(trbc_rec.seq).upper()

# ── Find feature boundaries in assembled.gb ────────────────────────────────────
def feat_bounds(record, label_fragment):
    for f in record.features:
        lbl = f.qualifiers.get('label', [''])[0]
        if label_fragment.lower() in lbl.lower():
            return int(f.location.start), int(f.location.end)
    return None, None

trbv_s,  trbv_e  = feat_bounds(asm_rec, 'TRBV')
cdr3b_s, cdr3b_e = feat_bounds(asm_rec, 'CDR3β')
trbc_s,  trbc_e  = feat_bounds(asm_rec, 'TRBC')
t2a_s,   t2a_e   = feat_bounds(asm_rec, 'T2A')
trav_s,  trav_e  = feat_bounds(asm_rec, 'TRAV')
cdr3a_s, cdr3a_e = feat_bounds(asm_rec, 'CDR3α')
trac_s,  trac_e  = feat_bounds(asm_rec, 'TRAC')

print("Assembled.gb feature positions:")
for name, s, e in [('TRBV',trbv_s,trbv_e),('CDR3b',cdr3b_s,cdr3b_e),
                   ('TRBC',trbc_s,trbc_e),('T2A',t2a_s,t2a_e),
                   ('TRAV',trav_s,trav_e),('CDR3a',cdr3a_s,cdr3a_e),
                   ('TRAC',trac_s,trac_e)]:
    print(f"  {name}: {s}..{e}")


# ── Cut-site positions from BioPython on the gene blocks ──────────────────────
def cuts(seq, enz):
    rb = RestrictionBatch([enz])
    return Analysis(rb, Seq(seq)).full().get(enz, [])

var_bbsi  = cuts(var_seq, BbsI)    # [5, 497] expected
var_bsmbi = cuts(var_seq, BsmBI)   # [85, 104] expected
trbc_bsmbi= cuts(trbc_seq, BsmBI)  # two BsmBI sites

print(f"\nBbsI cuts in variable fragment:  {var_bbsi}")
print(f"BsmBI cuts in variable fragment: {var_bsmbi}")
print(f"BsmBI cuts in TRBC-T2A:          {trbc_bsmbi}")


# ── Extract junction context from the ASSEMBLED sequence ──────────────────────
# IMPORTANT: The OH sits WITHIN the upstream gene (not at the domain boundary).
# For the table we show context centred on the DOMAIN BOUNDARY (where genes/domains
# change), so the reader sees the amino-acid junction directly.
# We then separately report where the OH lies relative to that boundary.

CODON = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

# Reading frame: LH arm ends at complete TCT codon; P2A starts at position 302 on a clean codon boundary.
FRAME_ORIGIN = 302   # position 302 is always at a codon boundary


def nt_to_aa(nt):
    return ''.join(CODON.get(nt[i:i+3], '?') for i in range(0, len(nt)-len(nt)%3, 3))


def domain_junction(seq, boundary, ctx=CTX, frame_origin=FRAME_ORIGIN):
    """
    Return a view centred on the domain boundary position (0-based).
    The context window is extended left/right to start on a codon boundary
    so the protein translation aligns cleanly.

    Returns: (left_nt, right_nt, left_aa, right_aa, oh_offset_from_boundary)
      left_nt / right_nt: CTX nt each side, start aligned to frame
      left_aa / right_aa: translated left/right context
    """
    # Align left edge to nearest codon boundary at or before (boundary - ctx)
    raw_left_start = boundary - ctx
    codon_rem = (raw_left_start - frame_origin) % 3
    left_start = raw_left_start - codon_rem   # step back to codon boundary

    left_nt  = seq[left_start : boundary]
    right_nt = seq[boundary   : boundary + ctx]

    # Trim left_nt to multiple of 3 from the right
    trim = len(left_nt) % 3
    left_nt_aligned = left_nt[trim:]   # drop leading fractional codon
    left_aa  = nt_to_aa(left_nt_aligned)
    right_aa = nt_to_aa(right_nt[:len(right_nt) - len(right_nt) % 3])

    return left_nt[-ctx:], right_nt[:ctx], left_aa, right_aa


def find_oh_near(seq, boundary, oh_seq, search=80):
    """
    Find the occurrence of oh_seq closest to boundary in a ±search window.
    Some OHs sit downstream of their boundary (e.g. oh_1b is the first few nt of
    TRAC, positive offset from the CDR3α/TRAC boundary) so we search both sides.
    Returns the position of the OH closest to boundary.
    """
    window_start = max(0, boundary - search)
    window = seq[window_start : boundary + search]
    candidates = []
    idx = 0
    while True:
        idx = window.find(oh_seq, idx)
        if idx == -1:
            break
        abs_pos = window_start + idx
        candidates.append((abs(abs_pos - boundary), abs_pos))
        idx += 1
    if not candidates:
        raise ValueError(f"OH '{oh_seq}' not found within {search} nt of position {boundary}")
    # Return position closest to boundary
    return min(candidates, key=lambda t: t[0])[1]


EXP = {'1a': 'CCAG', '1b': 'GGAA', '2a': 'GGAA', '2b': 'AGGA'}

# Build junction views centred on DOMAIN BOUNDARIES
junctions = {
    '1a': (trbv_e,  EXP['1a'], 'TRBV16 ends (Cys104)',       'CDR3β body begins'),
    '1b': (trac_s,  EXP['1b'], 'CDR3α body ends',             'TRAC-LVL-Cys begins'),
    '2a': (trbc_s,  EXP['2a'], 'CDR3β body ends',             'TRBC-LVL-Cys begins'),
    '2b': (trav_s,  EXP['2b'], 'T2A ends (oh_2b from T2A)',   'TRAV begins'),
}

jxn_data = {}
print("\nDomain-boundary junction context (assembled product):")
for jid, (boundary, oh, left_desc, right_desc) in junctions.items():
    ln, rn, la, ra = domain_junction(asm_seq, boundary)
    oh_abs  = find_oh_near(asm_seq, boundary, oh)
    oh_rel  = oh_abs - boundary   # negative = OH is upstream of boundary
    print(f"  {jid}  boundary={boundary}"
          f"  OH@{oh_abs}({oh_rel:+d} from boundary)"
          f"  left=...{ln}|  |{rn}...  right")
    print(f"        AA: ...{la}|  |{ra}...")
    jxn_data[jid] = dict(
        left_nt=ln, right_nt=rn, left_aa=la, right_aa=ra,
        oh=oh, oh_rel=oh_rel,
        left_desc=left_desc, right_desc=right_desc,
    )

# Verify OHs match expected
for jid, d in jxn_data.items():
    oh_found = asm_seq[find_oh_near(asm_seq,
                    junctions[jid][0], d['oh']) : find_oh_near(asm_seq,
                    junctions[jid][0], d['oh']) + 4]
    assert oh_found == d['oh'], f"OH mismatch at {jid}: {oh_found} vs {d['oh']}"
print("\nAll OH sequences verified ✓")


# ── Build the DataFrame ────────────────────────────────────────────────────────
# Context is centred on DOMAIN BOUNDARIES (where biological features change),
# not on the OH position.  The OH sits upstream of the boundary by oh_rel nt.

rows = [
    dict(step='1', enzyme='BbsI',  jxn='oh_1a\n(TRBV-specific)',
         left_block='LH-TRBV16', right_block='Variable (2D2)',
         **jxn_data['1a'],
         note=f'oh_1a is {EXP["1a"]} nt before end of TRBV16\n'
              f'TRBV encodes up to Cys104 (=C of CDR3β)\n'
              f'CDR3β body stripped of leading C in design'),
    dict(step='1', enzyme='BbsI',  jxn='oh_1b\n(constant)',
         left_block='Variable (2D2)', right_block='TRAC-RH',
         **jxn_data['1b'],
         note=f'oh_1b is {EXP["1b"]} — start of TRAC-LVL-Cys\n'
              f'J-C split codon (D) reconstituted at this junction\n'
              f'Same OH for every TCR clonotype'),
    dict(step='2', enzyme='BsmBI', jxn='oh_2a\n(constant)',
         left_block='Step-1 product\n(CDR3β end)', right_block='TRBC-T2A',
         **jxn_data['2a'],
         note=f'oh_2a is {EXP["2a"]} — start of TRBC1-Cys\n'
              f'J-C split codon (E) reconstituted at this junction\n'
              f'Same OH for every TCR clonotype'),
    dict(step='2', enzyme='BsmBI', jxn='oh_2b\n(constant)',
         left_block='TRBC-T2A', right_block='Step-1 product\n(TRAV start)',
         **jxn_data['2b'],
         note=f'oh_2b derived from WITHIN T2A sequence\n'
              f'Ensures T2A is fully reconstituted (no extra nt)\n'
              f'Placing an artificial OH here would +1 frameshift TRAV'),
]

df = pl.DataFrame({
    'Step'         : [r['step']       for r in rows],
    'Enzyme'       : [r['enzyme']     for r in rows],
    'Junction'     : [r['jxn']        for r in rows],
    '← Upstream'  : [r['left_block'] for r in rows],
    'Upstream AA'  : [r['left_aa']    for r in rows],
    'Upstream nt'  : [r['left_nt']    for r in rows],
    '|'            : ['←↓→']*4,
    'Downstream nt': [r['right_nt']   for r in rows],
    'Downstream AA': [r['right_aa']   for r in rows],
    'Downstream →' : [r['right_block']for r in rows],
    'OH & position': [f'{r["oh"]}\n({r["oh_rel"]:+d} nt from boundary)'
                      for r in rows],
    'Design note'  : [r['note']       for r in rows],
})


# ── Build great_tables ─────────────────────────────────────────────────────────
BBSI_DK  = '#1E40AF'; BBSI_LT  = '#DBEAFE'
BSMBI_DK = '#92400E'; BSMBI_LT = '#FEF3C7'
MONO     = gt.google_font('JetBrains Mono')

SEQ_COLS = ['Upstream nt', 'Downstream nt']
AA_COLS  = ['Upstream AA', 'Downstream AA']

tbl = (
    gt.GT(df, groupname_col='Step')
    .tab_header(
        title='Golden Gate Junction Sequences — 2D2 (TRBV16 / TRAV9N-3)',
        subtitle=(
            'Context is centred on each DOMAIN BOUNDARY (where biological features change), '
            'not on the OH position. '
            'Reading left-to-right across the ↓ column gives the continuous assembled sequence. '
            'Upstream AA / Downstream AA show the protein translation on each side.'
        ),
    )

    # Column spanners
    .tab_spanner(label='Upstream of boundary',
                 columns=['← Upstream', 'Upstream AA', 'Upstream nt'])
    .tab_spanner(label='Downstream of boundary',
                 columns=['Downstream nt', 'Downstream AA', 'Downstream →'])

    # Hide Step/Enzyme (already in row-group labels)
    .cols_hide(columns=['Step', 'Enzyme'])

    # Alignment
    .cols_align(align='right',  columns=['Upstream AA', 'Upstream nt'])
    .cols_align(align='center', columns=['|', 'OH & position', 'Junction'])
    .cols_align(align='left',   columns=['Downstream nt', 'Downstream AA'])

    # Column widths
    .cols_width({
        'Junction': '105px', '← Upstream': '115px',
        'Upstream AA': '110px', 'Upstream nt': '155px',
        '|': '40px',
        'Downstream nt': '155px', 'Downstream AA': '110px',
        'Downstream →': '115px',
        'OH & position': '90px', 'Design note': '210px',
    })

    # Monospace + bold for nucleotide sequence columns
    .tab_style(
        style=gt.style.text(font=MONO, size='12px', weight='bold', color='#1E293B'),
        locations=gt.loc.body(columns=SEQ_COLS),
    )
    # Monospace (normal weight) for AA columns
    .tab_style(
        style=gt.style.text(font=MONO, size='11px', color='#374151'),
        locations=gt.loc.body(columns=AA_COLS),
    )

    # Light fill on sequence/AA columns
    .tab_style(
        style=gt.style.fill(color='#F0F9FF'),
        locations=gt.loc.body(columns=SEQ_COLS + ['Upstream AA']),
    )
    .tab_style(
        style=gt.style.fill(color='#F0FFF4'),
        locations=gt.loc.body(columns=['Downstream nt', 'Downstream AA']),
    )

    # OH & position — BbsI rows (0,1) blue; BsmBI rows (2,3) amber
    .tab_style(
        style=[gt.style.fill(color=BBSI_LT),
               gt.style.text(color=BBSI_DK, weight='bold',
                              font=MONO, size='13px')],
        locations=gt.loc.body(columns='OH & position', rows=[0, 1]),
    )
    .tab_style(
        style=[gt.style.fill(color=BSMBI_LT),
               gt.style.text(color=BSMBI_DK, weight='bold',
                              font=MONO, size='13px')],
        locations=gt.loc.body(columns='OH & position', rows=[2, 3]),
    )

    # Centre divider column
    .tab_style(
        style=[gt.style.fill(color='#F1F5F9'),
               gt.style.text(color='#94A3B8', size='16px')],
        locations=gt.loc.body(columns='|'),
    )

    # Row-group header colours
    .tab_style(
        style=[gt.style.fill(color='#EFF6FF'),
               gt.style.text(color=BBSI_DK, weight='bold')],
        locations=gt.loc.row_groups(rows=[0, 1]),
    )
    .tab_style(
        style=[gt.style.fill(color='#FFFBEB'),
               gt.style.text(color=BSMBI_DK, weight='bold')],
        locations=gt.loc.row_groups(rows=[2, 3]),
    )

    .tab_source_note(
        'Nucleotide context = 14 nt each side of the domain boundary in the final assembled product '
        '(no restriction-site artefacts).  '
        'AA translation is in reading frame from the clean codon boundary at P2A start (position 302).  '
        'OH position: negative = OH lies that many nt upstream of the domain boundary within the upstream gene.  '
        'oh_1b = oh_2a = GGAA in this design — not a conflict; they are used in separate enzyme reactions.'
    )
    .opt_table_font(font=gt.google_font('Inter'))
    .opt_stylize(style=3)
)

html_out = OUT / 'assembly_table.html'
html_out.write_text(tbl.as_raw_html())
print(f'\nSaved → {html_out}')
print('Open in a browser to view; use browser Print → Save as PDF to export.')
