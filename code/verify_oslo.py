"""Exact check of Theorem A, Proposition 6.1 and Corollary B for the Oslo rice pile.

Usage: python3 verify_oslo.py [Lmax]   (default Lmax = 5; Lmax = 6 takes a few minutes)

Configurations are slope vectors in {0,1,2}^L. The operator a_x adds a particle
at site x and relaxes with Dhar's rule: slope 0 -> 1; slope 1 -> 2 with
probability p, otherwise topple to 0; slope 2 -> topple to 1. A toppling at x
sends one particle to x-1 and one to x+1; site 0 absorbs, site L+1 means L.
All arithmetic is exact (fractions), for two values of p.
"""
import itertools
import sys
from fractions import Fraction
from functools import lru_cache


def build(L, p):
    """Return (confs, apply, mono) for the pile on L sites with probability p."""
    q = 1 - p

    @lru_cache(maxsize=None)
    def a(i, c):
        """Law of a_i applied to the configuration c, as a tuple of (config, prob)."""
        if i == 0:
            return ((c, Fraction(1)),)
        if i == L + 1:
            i = L
        s = c[i - 1]
        setv = lambda d, v: d[:i - 1] + (v,) + d[i:]

        def topple(c0):
            res = {}
            for k, v in a(i - 1, c0):
                for k2, v2 in a(i + 1, k):
                    res[k2] = res.get(k2, 0) + v * v2
            return res

        if s == 0:
            return ((setv(c, 1), Fraction(1)),)
        out = {}
        if s == 1:
            out[setv(c, 2)] = p
            for k, v in topple(setv(c, 0)).items():
                out[k] = out.get(k, 0) + q * v
        else:
            out = topple(setv(c, 1))
        return tuple(out.items())

    confs = list(itertools.product(range(3), repeat=L))

    def apply(i, vec):
        new = {}
        for c, w in vec.items():
            for k, v in a(i, c):
                new[k] = new.get(k, 0) + w * v
        return {k: v for k, v in new.items() if v != 0}

    def mono(m):
        vec = {(0,) * L: Fraction(1)}
        for i in range(L):
            for _ in range(m[i]):
                vec = apply(i + 1, vec)
        return vec

    for c in confs:
        assert sum(apply(1, {c: Fraction(1)}).values()) == 1
    return confs, apply, mono


def run(L, p):
    confs, apply, mono = build(L, p)
    two, one = (2,) * L, (1,) * L
    Psi = mono(two)
    assert apply(1, Psi) == Psi and apply(1, {two: Fraction(1)}) == Psi
    vu = mono((1,) + (2,) * (L - 1))
    N = L * (L + 1)
    lin = lambda x, X, y, Y: {k: v for k in set(X) | set(Y)
                              if (v := x * X.get(k, 0) + y * Y.get(k, 0)) != 0}
    # Theorem A
    K = {}
    for s in confs:
        g = sum((j + 1) * s[j] for j in range(L))
        K[s] = N - g + (s[-1] == 2)
        r = 0
        while r < L and s[L - 1 - r] == 2:
            r += 1
        x = {s: Fraction(1)}
        for k in range(K[s] - 1):
            assert x != Psi
            x = apply(1, x)
        target = vu if r == 0 else lin(p ** (-r), Psi, 1 - p ** (-r), vu)
        assert x == target and x != Psi, s
        assert apply(1, x) == Psi, s
    # Proposition 6.1 and Lemma 6.2
    def inR(c):
        for i in range(L):
            if c[i] == 0:
                j = i + 1
                while j < L and c[j] == 1:
                    j += 1
                if j >= L or c[j] != 2:
                    return False
        return True
    rec = set(Psi)
    assert all(inR(c) for c in rec) and one in rec
    assert all(K[c] <= N // 2 for c in confs if inR(c))
    # Corollary B
    assert max(K.values()) == N and [c for c in confs if K[c] == N] == [(0,) * L]
    assert max(K[c] for c in rec) == N // 2 and K[one] == N // 2
    dhar_fails = sorted(c for c in rec if c[-1] == 2)
    assert dhar_fails, "a recurrent configuration with last slope 2 exists"
    return len(rec), len(dhar_fails), set(rec) == {c for c in confs if inR(c)}


if __name__ == "__main__":
    Lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for L in range(1, Lmax + 1):
        for p in (Fraction(3, 7), Fraction(4, 5)):
            nrec, nfail, equal = run(L, p)
            print(f"L={L} p={p}: Theorem A holds for all {3 ** L} configurations; "
                  f"{nrec} recurrent, all in R (supp Psi = R: {equal}); "
                  f"identity (1) fails for {nfail} recurrent configurations")
    print("all checks passed")
