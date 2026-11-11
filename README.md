# The Oslo Rice Pile Relaxes Exactly in L(L+1) Steps

Paper (`paper/main.tex`, `paper/main.pdf`) and the program that rechecks its statements in exact rational arithmetic for small piles.

```
cd code && python3 verify_oslo.py 5   # Theorem A, Proposition 6.1, Corollary B for L <= 5; no dependencies
python3 verify_oslo.py 6               # the same for L <= 6, as stated in the paper; about 20 minutes
python3 make_figure.py                # fig_relaxation.png; needs matplotlib
```

DOI (all versions): [10.5281/zenodo.24123896](https://doi.org/10.5281/zenodo.24123896)

Author: Deep Bhattacharjee (GitHub-related contact: egspl@proton.me). License: MIT.
