"""Draw fig_relaxation.png: K(sigma) against g(sigma) for all 3^5 configurations at L = 5.

Recurrent configurations are those charged by the exact steady state Psi = a_1^2...a_L^2|0>.
Usage: python3 make_figure.py [output.png]   (needs matplotlib)
"""
import sys
from collections import Counter
from fractions import Fraction

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from verify_oslo import build

L = 5
N = L * (L + 1)
confs, apply, mono = build(L, Fraction(1, 2))
rec = set(mono((2,) * L))
pts = Counter()
for c in confs:
    g = sum((j + 1) * c[j] for j in range(L))
    K = N - g + (c[-1] == 2)
    pts[(g, K, c[-1] == 2, c in rec)] += 1

fig, ax = plt.subplots(figsize=(7.2, 4.3))
gs = list(range(0, N + 1))
ax.plot(gs, [N - g for g in gs], color="0.6", lw=1, zorder=1, label=r"$K=L(L+1)-g$")
ax.plot(gs, [N + 1 - g for g in gs], color="0.6", lw=1, ls=":", zorder=1, label=r"$K=L(L+1)-g+1$")
ax.axhline(N // 2, color="#b2182b", ls="--", lw=1.2, zorder=1, label=r"$L(L+1)/2=15$")
for (g, K, last2, isrec), n in pts.items():
    col = "#e08214" if last2 else "#2166ac"
    ax.scatter([g], [K], s=18 + 7 * n, marker="s" if last2 else "o",
               facecolors=col if isrec else "white", edgecolors=col, linewidths=1.1, zorder=3)
ax.scatter([], [], marker="o", facecolors="#2166ac", edgecolors="#2166ac", label=r"$\sigma_5\neq2$, recurrent")
ax.scatter([], [], marker="o", facecolors="white", edgecolors="#2166ac", label=r"$\sigma_5\neq2$, transient")
ax.scatter([], [], marker="s", facecolors="#e08214", edgecolors="#e08214", label=r"$\sigma_5=2$, recurrent")
ax.scatter([], [], marker="s", facecolors="white", edgecolors="#e08214", label=r"$\sigma_5=2$, transient")
ax.annotate("empty pile: $K=30$", xy=(0, 30), xytext=(4, 28.5), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.set_xlabel(r"number of grains $g(\sigma)$")
ax.set_ylabel(r"relaxation time $K(\sigma)$")
ax.set_xlim(-1, 31)
ax.set_ylim(-0.5, 31.5)
ax.legend(fontsize=7.5, loc="upper right", frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
out = sys.argv[1] if len(sys.argv) > 1 else "fig_relaxation.png"
plt.savefig(out, dpi=200)
print(f"{len(rec)} recurrent; largest K on recurrent configurations = "
      f"{max(K for (g, K, l, r) in pts if r)}; written {out}")
