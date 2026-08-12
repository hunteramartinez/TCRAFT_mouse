#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Design PCR primers for the mouse TRAC left and right homology arms.
Intended for HiFi polymerase (Q5, KAPA HiFi, etc.).

Criteria:
  - Length 18-28 nt
  - At least 1 GC clamp at 3' end
  - No homopolymer run >=5
  - Wallace Tm target: 60-72 C (corresponds to ~55-65 C actual Tm for HiFi)
  - Prefer primer closest to midpoint of target Tm range
"""
from pathlib import Path
import math

_REF = Path(__file__).parent / 'TCRAFT' / 'references' / 'mouse' / 'mm_trac_HDR_sequences.fa'
_OUT = Path(__file__).parent / 'custom_tcrs' / 'ha_primers.tsv'

_COMP = str.maketrans('ACGTacgt', 'TGCAtgca')

def rev_comp(s):
    return s.translate(_COMP)[::-1]

def gc_count(s):
    return sum(1 for b in s.upper() if b in 'GC')

def tm_wallace(s):
    """Simple Wallace-rule Tm (°C). Overestimates for long oligos but useful for ranking."""
    s = s.upper()
    gc = gc_count(s)
    return 4 * gc + 2 * (len(s) - gc)

def longest_hp(s):
    if not s:
        return 0
    best = cur = 1
    for i in range(1, len(s)):
        cur = cur + 1 if s[i].upper() == s[i-1].upper() else 1
        best = max(best, cur)
    return best

def parse_fasta(path):
    seqs, hdr, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if hdr:
                    seqs[hdr] = ''.join(parts)
                hdr, parts = line[1:].strip(), []
            else:
                parts.append(line)
    if hdr:
        seqs[hdr] = ''.join(parts)
    return seqs

def design_primer(seq, reverse=False, min_len=18, max_len=28, tm_lo=60, tm_hi=72):
    """
    Slide over primer lengths and pick the best by:
      1. GC clamp at 3' end (required)
      2. No homopolymer ≥5 (required)
      3. Wallace Tm in [tm_lo, tm_hi]  (preferred; falls back to nearest)
    If reverse=True the primer anneals to the 3' end of seq (rev-comp is returned).
    """
    seq = seq.upper()
    candidates = []
    for n in range(min_len, max_len + 1):
        region = seq[-n:] if reverse else seq[:n]
        primer = rev_comp(region) if reverse else region
        tm = tm_wallace(primer)
        gc = gc_count(primer)
        gc_pct = gc / n * 100
        clamp = primer[-1] in 'GC'
        hp = longest_hp(primer)
        candidates.append(dict(
            primer=primer, length=n, tm=round(tm, 1),
            gc_pct=round(gc_pct, 1), clamp=clamp, hp=hp
        ))

    # Hard filters
    valid = [c for c in candidates if c['clamp'] and c['hp'] < 5]
    if not valid:
        valid = [c for c in candidates if c['clamp']]  # relax HP
    if not valid:
        valid = candidates                               # relax clamp too (should never happen)

    # Prefer in-range Tm
    in_range = [c for c in valid if tm_lo <= c['tm'] <= tm_hi]
    pool = in_range if in_range else valid
    # Among candidates in pool, minimise |Tm - midpoint|; break ties by shorter primer
    mid = (tm_lo + tm_hi) / 2
    return min(pool, key=lambda c: (abs(c['tm'] - mid), c['length']))

def main():
    seqs = parse_fasta(_REF)
    lh = seqs['mm_trac_left_homology_arm'].replace('-', '').upper()
    rh = seqs['mm_trac_right_homology_arm'].replace('-', '').upper()

    print(f"LH arm: {len(lh)} bp")
    print(f"  5' end: {lh[:30]}...")
    print(f"  3' end: ...{lh[-30:]}")
    print(f"RH arm: {len(rh)} bp")
    print(f"  5' end: {rh[:30]}...")
    print(f"  3' end: ...{rh[-30:]}")
    print()

    # Two outer primers that amplify the FULLY ASSEMBLED HDR template (post-golden-gate).
    # Template for PCR = assembled golden gate product.
    # Product = complete linear HDR template [LH--insert--RH] for RNP electroporation.
    rows = [
        ('HDR_F', 'Left_HA_5prime',  'Forward',
         'Anchors at very 5-prime end of LH arm',
         design_primer(lh, reverse=False)),
        ('HDR_R', 'Right_HA_3prime', 'Reverse',
         'Anchors at very 3-prime end of RH arm',
         design_primer(rh, reverse=True)),
    ]

    hdr_line = f"{'Name':<7}  {'Primer (5->3)':<30}  {'Len':>3}  {'Tm_W':>5}  {'GC%':>5}  {'Clamp'}"
    print(hdr_line)
    print('-' * len(hdr_line))
    for name, region, direction, note, p in rows:
        print(f"{name:<7}  {p['primer']:<30}  {p['length']:>3}  {p['tm']:>5}  {p['gc_pct']:>5}  {p['clamp']}")
        print(f"         ({note})")

    # Save TSV
    _OUT.parent.mkdir(exist_ok=True)
    with open(_OUT, 'w') as fh:
        fh.write("name\tanchor_region\tdirection\tnote\tsequence_5to3\tlength_nt\ttm_wallace_C\tgc_pct\tgc_clamp\n")
        for name, region, direction, note, p in rows:
            fh.write(f"{name}\t{region}\t{direction}\t{note}\t{p['primer']}\t{p['length']}\t{p['tm']}\t{p['gc_pct']}\t{p['clamp']}\n")
    print(f"\nSaved -> {_OUT}")

if __name__ == '__main__':
    main()
