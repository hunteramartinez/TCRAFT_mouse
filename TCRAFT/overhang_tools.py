"""
Overhang-set scoring and search for TCRAFT-style pooled Golden Gate assembly.

Scoring follows Gaglione et al. (Immunity 2026), Methods "Overhang design":
  p(O | S)  = 1/2 * [ r(O, RC O) / sum_{s in S} r(O, s)  +  r(O, RC O) / sum_{s in S} r(RC O, s) ]
  f1(W)     = min_i  p(OL_i | T) * p(OR_i | T)        (worst-case correct circular assembly)
  f2(W)     = prod_i p(OL_i | L) * p(OR_i | R)        (orthogonality regulariser)
  f(W)      = w*f1 + (1-w)*f2,  w = 0.5
where L / R are the sets of left / right overhangs and their reverse complements and T = L u R.
r() comes from the NEB 4-nt ligation-fidelity tables (references/neb_tables), symmetrised.

Each "vector" i carries a (left, right) overhang pair; pairs must be unique within a pool.
"""
import random
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

NEB_DIR = Path(__file__).resolve().parent / "references" / "neb_tables"
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def rc(s):
    return s.translate(_COMP)[::-1]


class Fidelity:
    """NEB ligation matrix for one enzyme with fast p()/f1/f2 evaluation."""

    def __init__(self, enzyme):
        m = pd.read_csv(NEB_DIR / f"{enzyme}_ligation_efficiency.csv", index_col=0)
        m.columns = [c.strip() for c in m.columns]
        m.index = [i.strip() for i in m.index]
        m = m.loc[m.index, m.index]  # same ordering on both axes
        self.enzyme = enzyme
        self.labels = list(m.index)
        self.idx = {o: i for i, o in enumerate(self.labels)}
        r = m.values.astype(float)
        self.R = (r + r.T) / 2.0  # NEB matrix is symmetric only up to noise
        self.rc_idx = np.array([self.idx[rc(o)] for o in self.labels])
        self.wc = self.R[np.arange(len(self.labels)), self.rc_idx]  # r(O, RC O)

    def _bg(self, ohs):
        """Row-sum of r over background S = unique(ohs u RC(ohs)) for every overhang."""
        s = set()
        for o in ohs:
            s.add(self.idx[o])
            s.add(self.idx[rc(o)])
        return self.R[:, sorted(s)].sum(axis=1)

    def p_all(self, bg_sum):
        """p(O|S) for every one of the 256 overhangs given background row-sums."""
        a, b = bg_sum, bg_sum[self.rc_idx]
        with np.errstate(divide="ignore", invalid="ignore"):
            return 0.5 * self.wc * (1.0 / a + 1.0 / b)

    def score(self, lefts, rights, w=0.5):
        """-> dict(f, f1, f2, p_vec) for aligned lists of left / right overhangs."""
        li = np.array([self.idx[o] for o in lefts])
        ri = np.array([self.idx[o] for o in rights])
        pT = self.p_all(self._bg(list(lefts) + list(rights)))
        pL = self.p_all(self._bg(lefts))
        pR = self.p_all(self._bg(rights))
        p_vec = pT[li] * pT[ri]
        f1 = float(p_vec.min())
        f2 = float(np.prod(pL[li] * pR[ri]))
        return dict(f=w * f1 + (1 - w) * f2, f1=f1, f2=f2, p_vec=p_vec)


def is_palindromic(o):
    return o == rc(o)


def optimize(fid, left_cands, right_cands, seed=0, w=0.5,
             stage1=(8, 10_000, 1e-4), stage2=(8, 5_000, 1e-3), patience=1_000, init=None):
    """
    Two-stage Metropolis search (as in the paper): 8 random starts x 10k iters @ T=1e-4, then
    8 copies of the best x 5k iters @ T=1e-3.

    left_cands / right_cands: list (one entry per vector) of lists of allowed 4-mers.
    Pairs (left, right) must stay unique across vectors.
    Returns (lefts, rights, score_dict).
    """
    rng = random.Random(seed)
    n = len(left_cands)

    def rand_init():
        for _ in range(10_000):
            L = [rng.choice(c) for c in left_cands]
            R = [rng.choice(c) for c in right_cands]
            if len(set(zip(L, R))) == n:
                return L, R
        raise RuntimeError("could not draw a unique-pair starting set")

    def run(L, R, iters, T):
        cur = fid.score(L, R, w)["f"]
        best, best_state, stall = cur, (list(L), list(R)), 0
        for _ in range(iters):
            i = rng.randrange(n)
            side = rng.random() < 0.5
            L2, R2 = list(L), list(R)
            if side:
                L2[i] = rng.choice(left_cands[i])
            else:
                R2[i] = rng.choice(right_cands[i])
            if len(set(zip(L2, R2))) < n:
                continue
            new = fid.score(L2, R2, w)["f"]
            if new >= cur or rng.random() < np.exp((new - cur) / T):
                L, R, cur = L2, R2, new
            if cur > best + 1e-12:
                best, best_state, stall = cur, (list(L), list(R)), 0
            else:
                stall += 1
                if stall >= patience:
                    break
        return best, best_state

    starts = [init] if init else [rand_init() for _ in range(stage1[0])]
    results = [run(list(L), list(R), stage1[1], stage1[2]) for L, R in starts]
    best, (L, R) = max(results, key=lambda t: t[0])
    results2 = [run(list(L), list(R), stage2[1], stage2[2]) for _ in range(stage2[0])]
    best, (L, R) = max(results2 + [(best, (L, R))], key=lambda t: t[0])
    return L, R, fid.score(L, R, w)


# ---------------------------------------------------------------------------
# Candidate overhangs: 4-mers a sequence can carry via synonymous recoding only
# ---------------------------------------------------------------------------
from Bio.Data.CodonTable import standard_dna_table  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

SYN = {}
for _cod, _aa in standard_dna_table.forward_table.items():
    SYN.setdefault(_aa, []).append(_cod.lower())

# BsaI, BsmBI, BbsI, SapI, XbaI, EcoRI, XhoI, EcoRV, NotI (same set Generate.py screens for)
SITES = ["GGTCTC", "CGTCTC", "GAAGAC", "GCTCTTC", "TCTAGA", "GAATTC", "CTCGAG", "GATATC", "GCGGCCGC"]


def has_site(s):
    s = s.upper()
    return any(x in s or rc(x) in s for x in SITES)


_HOMOPOLY = ("AAAAA", "CCCCC", "GGGGG", "TTTTT")


def count_problems(s):
    """Number of screened restriction-site hits (either strand) + homopolymer runs >= 5."""
    s = s.upper()
    n = 0
    for x in SITES:
        for y in {x, rc(x)}:
            i = s.find(y)
            while i != -1:
                n += 1
                i = s.find(y, i + 1)
    return n + sum(s.count(h) for h in _HOMOPOLY)


def clean_nt(nt, usage, protect=(), sites_only=False):
    """
    Remove screened restriction sites (and homopolymers >= 5 unless sites_only) from an
    in-frame coding sequence using synonymous codon swaps only.
    usage: {aa: [[codon, freq], ...]} (e.g. mouse_codon_usage.json), used to prefer common codons.
    protect: codon indices that must not change (e.g. an already-chosen overhang window).
    Greedy: repeatedly apply the swap that most reduces the problem count, ties -> higher usage.
    """
    freq = {c.lower(): f for aa, lst in usage.items() for c, f in lst}
    cost = (lambda s: len([1 for x in SITES for y in {x, rc(x)} if y in s.upper()])) if sites_only else count_problems
    codons = [nt[i:i + 3].lower() for i in range(0, len(nt) // 3 * 3, 3)]
    prot = str(Seq("".join(codons)).translate())
    protect = set(protect)
    cur = cost("".join(codons))
    while cur > 0:
        best = None
        seq = "".join(codons)
        for i, aa in enumerate(prot):
            if i in protect or len(SYN[aa]) < 2:
                continue
            lo, hi = max(0, (i - 5) * 3), (i + 6) * 3
            before = cost(seq[lo:hi])
            for alt in SYN[aa]:
                if alt == codons[i]:
                    continue
                trial = seq[:i * 3] + alt + seq[i * 3 + 3:]
                gain = before - cost(trial[lo:hi])
                if gain > 0 and (best is None or (gain, freq.get(alt, 0)) > best[0]):
                    best = ((gain, freq.get(alt, 0)), i, alt)
        if best is None:
            raise ValueError("clean_nt: could not remove all problems with synonymous swaps")
        codons[best[1]] = best[2]
        cur = cost("".join(codons))
    out = "".join(codons)
    assert str(Seq(out).translate()) == prot
    return out


def window_candidates(nt, lo, hi):
    """
    nt: in-frame coding sequence (any valid back-translation; only its protein and local
        context are used). Candidate 4-mers must lie fully inside nt[lo:hi].
    Returns {4mer: [(start_nt, first_codon_idx, recoded_span_nt), ...]}.
    Excludes palindromic 4-mers and any recoding that creates a screened restriction site.
    """
    prot = str(Seq(nt[: len(nt) // 3 * 3]).translate())
    out = {}
    for s in range(lo, hi - 3):
        c0, c1 = s // 3, (s + 3) // 3
        left_ctx = nt[max(0, (c0 - 2) * 3): c0 * 3]
        right_ctx = nt[(c1 + 1) * 3: (c1 + 3) * 3]
        for combo in product(*[SYN[prot[c]] for c in range(c0, c1 + 1)]):
            span = "".join(combo)
            if count_problems(left_ctx + span + right_ctx):
                continue
            oh = span[s - c0 * 3: s - c0 * 3 + 4].upper()
            if is_palindromic(oh):
                continue
            out.setdefault(oh, []).append((s, c0, span))
    return out
