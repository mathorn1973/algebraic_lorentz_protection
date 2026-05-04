# algebraic-lorentz-protection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20029689.svg)](https://doi.org/10.5281/zenodo.20029689)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC--BY%204.0-lightgrey.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org)

Verification scripts for the paper

> A. M. Thorn,
> *Dense cyclotomic subgroups of the Lorentz group and the
> Collins--Perez--Sudarsky obstruction*,
> 2026.

## Paper

The full preprint and its LaTeX source are included in this repository:

- **PDF**:   [`Thorn2026_lorentz_density.pdf`](Thorn2026_lorentz_density.pdf)
- **Source**: [`Thorn2026_lorentz_density.tex`](Thorn2026_lorentz_density.tex)

The paper bibliography references this verification suite by its
Zenodo DOI [10.5281/zenodo.20029689](https://doi.org/10.5281/zenodo.20029689).
Both the paper and the code are archived together on Zenodo via this
concept DOI; each tagged GitHub release produces a fresh version DOI
that resolves to that exact snapshot.

To recompile the paper from source, any modern LaTeX engine works.
We used [Tectonic](https://tectonic-typesetting.github.io/) 0.15.0:

```
tectonic Thorn2026_lorentz_density.tex
```

## Run

```
python3 algebraic_lorentz_protection.py
```

## Expected output

```
11/11 PASS
```

A full sample run is provided in `verification_output.txt`.

## What is verified

The script verifies eleven independent claims used in the paper:

1. Elliptic commutator: the commutator
   `C = B_J B_g B_J^{-1} B_g^{-1}` of two conjugate boosts of rapidity
   `log phi` whose axes meet at angle `arccos(1/sqrt 5)` has trace
   `2 + 2 c_omega` with `c_omega = (31 + 8 sqrt 5)/50`. Verified at
   200-digit precision.

2. Minimal polynomial: `c_omega` satisfies `2500 x^2 - 3100 x + 641 = 0`
   over Q, with discriminant `800^2 * 5` non-square, so `c_omega` has
   degree two over Q and lies in `Q(sqrt 5)`.

3. Cyclotomic cosine exclusion: `omega/pi` is irrational. The eight
   candidate cosines from cyclotomic orders `m in {5, 8, 10, 12}` are
   each numerically distinct from `c_omega`.

4. Molien series: Klein's formula for the icosahedral group `A_5`
   acting on `R^3`,

       M(t) = (1 + t^15) / ((1 - t^2)(1 - t^6)(1 - t^10)),

   is verified by comparing coefficients to the dimensions
   `dim((S^d R^3)^{A_5})` for `d = 0..16`.

5. Anisotropy degree: `dim((S^2 R^3)^{A_5}) = dim((S^4 R^3)^{A_5}) = 1`,
   `dim((S^6 R^3)^{A_5}) = 2`. The first non-isotropic invariant is at
   degree six.

6. Cubic dispersion exponent: spread of free dispersion over the
   momentum sphere fits `|k|^4` for cubic `O_h`-symmetric couplings.
   Fitted exponent within 0.05 of 4.

7. Icosahedral dispersion exponent: same fit gives `|k|^6` for
   `A_5`-symmetric couplings. Fitted exponent within 0.05 of 6.

8. Galois equivariance: `sigma_k * M_J * sigma_k^{-1} = M_{sigma_k(J)}`
   holds in the regular representation of `Z[zeta_5]` for all
   `k = 2, 3, 4`.

9. Elliptic classification of C: `C^T eta C = eta`, `det C = 1`,
   `rank(C - I) = 2` (no Jordan block on the eigenvalue 1), and
   the restriction of the Lorentz form to `ker(C - I)` has Gram
   matrix `diag(-1, (13 + sqrt 5)/2)` of signature (1, 1). Hence
   `C` is conjugate in `SO+(3,1)^0` to a spatial rotation by angle
   `omega`.

10. Explicit `A_5` element: the rotation matrix

         g = [[-1/2,         (1-sqrt 5)/4,  (1+sqrt 5)/4],
              [(-1+sqrt 5)/4, (1+sqrt 5)/4, 1/2         ],
              [-(1+sqrt 5)/4, 1/2,          (1-sqrt 5)/4]]

    satisfies `g^T g = I`, `det g = 1`, `g^3 = I`, `g != I`,
    `g != I` after one or two iterations, `g(v_1) = v_2` for the
    icosahedron vertices `v_1 = (0, 1, phi)`, `v_2 = (1, phi, 0)`,
    and the angle between them satisfies `cos alpha = 1/sqrt 5`.

11. Group membership: `g` permutes the set of twelve standard
    icosahedron vertices `{(0, +-1, +-phi), (+-1, +-phi, 0),
    (+-phi, 0, +-1)}` with orbit structure `[3, 3, 3, 3]`. Combined
    with `det g = 1`, this forces `g` into the rotational icosahedral
    group `I = A_5`.

## Tested environment

```
Python  3.10.12
sympy   1.14.0
mpmath  1.3.0
numpy   2.2.6
```

The verifications use exact rational and sympy-symbolic arithmetic
where possible, and 200-digit `mpmath` floating point for the trace
calculation. No claim depends on a fragile numerical fit. A
runtime-fitted dispersion exponent (claims 6 and 7) is the only
non-symbolic part and is given a generous tolerance of 0.05.

## Reproducibility

The full output of one run on the tested environment is saved in
`verification_output.txt` for reference.

Total runtime under one second on a modern CPU.

## Citation

If you use this verification suite, please cite both the article and
the code. Machine-readable metadata is in `CITATION.cff` (used
automatically by the GitHub "Cite this repository" button).

Plain-text suggested citation for the code:

    A. M. Thorn (2026), algebraic-lorentz-protection: verification
    scripts. Zenodo. https://doi.org/10.5281/zenodo.20029689

Plain-text suggested citation for the article:

    A. M. Thorn (2026), Dense cyclotomic subgroups of the Lorentz
    group and the Collins-Perez-Sudarsky obstruction. Preprint,
    archived with the verification code at
    https://doi.org/10.5281/zenodo.20029689

## License

This work is licensed under Creative Commons Attribution 4.0
International (CC BY 4.0). See `LICENSE` for the short notice and
links to the full legal text. You are free to share and adapt this
material for any purpose, including commercially, with attribution.
