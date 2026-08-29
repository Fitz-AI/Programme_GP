# Reproducing the counterexample of `note179.tex`

This directory contains the PARI/GP scripts referenced in §8
("Computational scripts"), the arXiv note refuting
Conjecture 10.2.8 of the author's 2016 thesis for $n=9$ via an explicit
pencil (draw no. 179 of campaign E03).

## Requirements

- PARI/GP **2.17.4** (the version these scripts were tested with; earlier
  2.1x releases should work but have not been checked).
- Two files already present elsewhere in this repository:
  `../libre/qfsolve.gp` (Denis Simon, for `QfWittinvariant`) and
  `../programme/Changements_de_bases.gp` (for `Hyperbolique`).
  `audit_complet.gp` and `certificat_decomposition.gp` `read()` them via
  these relative paths, so this directory must stay a sibling of `libre/`
  and `programme/` at the repository root.


## The six scripts and what each one establishes

| Script | Establishes | Cited as |
|---|---|---|
| `check_singularite.gp` | The pair $(A,B)$ satisfies Condition 1: $\det A,\det B\ne0$, $\gcd(\Delta,\Delta')=1$. | Proposition 3.1 |
| `preuve_v2.gp` | $\det(\lambda A+\mu B)$ is odd at the three points of $\mathbb{P}^1(\mathbb{F}_2)$, hence at every primitive integer pair. | Lemma 4.1 |
| `verif_12classes.gp` | The 48 primitive pairs mod 8 fall into exactly 12 orbits under the unit action, and the 12 representatives used cover them without duplication. | Corollary 4.3 |
| `certificat_decomposition.gp` | For each of the 12 representatives, an exact rational $P$ with $\det P\ne0$ such that $PQP^{\mathsf t}=\mathbb{H}^3\perp R$ genuinely block-diagonal, and $R$ integral with odd determinant. | Equation (2), Table 1 |
| `audit_complet.gp` | Stage 2's exhaustive mod-64 enumeration and the Hasse-invariant cross-check at the 12 classes, plus the two positive controls (36 verdicts, 36 agreements). | Table 1, Remark 5.3 |
| `arbitre_mod64.gp` | Standalone validation of the mod-64 exhaustive-enumeration routine on three cases with known answers (one anisotropic, two isotropic, one non-diagonal). | Lemma 5.1 |



