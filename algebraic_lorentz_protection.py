"""
algebraic_lorentz_protection.py

Verification scripts for the paper

  "Dense cyclotomic subgroups of the Lorentz group and the
   Collins-Perez-Sudarsky obstruction"
  A. M. Thorn, 2026

Run end to end:

  python3 algebraic_lorentz_protection.py

Dependencies: sympy, numpy, mpmath.

Output: PASS/FAIL for eleven claims. All eleven pass in under one minute.

  [1]  Elliptic commutator: C has eigenvalues {1, 1, e^iw, e^-iw}
       with cos(w) = (31 + 8 sqrt 5) / 50.
  [2]  Minimal polynomial: 2500 x^2 - 3100 x + 641 vanishes at cos(w).
  [3]  Cyclotomic cosine exclusion: w/pi is irrational.
  [4]  Molien series for icosahedral group A_5 acting on R^3.
  [5]  First anisotropic A_5-invariant on R^3 is at degree 6.
  [6]  Cubic spatial coupling: dispersion anisotropy ~ |k|^4.
  [7]  Icosahedral spatial coupling: dispersion anisotropy ~ |k|^6.
  [8]  Galois equivariance: sigma_k . M_J . sigma_k^{-1} == M_{sigma_k(J)}.
  [9]  Elliptic classification: rank(C-I) = 2 and the fixed plane carries
       a non-degenerate Lorentz form of signature (1,1).
  [10] Explicit A_5 element: order-3 rotation g with g(v_1) = v_2 for
       icosahedron vertices v_1 = (0,1,phi), v_2 = (1,phi,0).
  [11] Group membership: g permutes the set of 12 standard icosahedron
       vertices, hence g lies in the rotational icosahedral group I.
"""

import sympy as sp
import numpy as np
import mpmath as mp


# ============================================================================
# Section 1.  Elliptic commutator and rotation cosine
# ============================================================================
#
# Build the Lorentz boost B(n, rho) of rapidity rho along unit vector n,
# then form the commutator C = B_J B_g B_J^-1 B_g^-1, where B_J is along +z
# and B_g is along the direction making angle alpha with +z, with cos(alpha)
# = 1/sqrt(5) (the angle between adjacent icosahedral vertex axes), and
# rho = log(phi).
#
# C is a pure elliptic Lorentz element with eigenvalues {1, 1, e^iw, e^-iw}
# and cos(w) = (Tr(C) - 2)/2.  Verified at 200 digits to match the closed
# form (31 + 8 sqrt 5) / 50.

mp.mp.dps = 200

phi_n = (1 + mp.sqrt(5)) / 2
cosh_rho = mp.sqrt(5) / 2
sinh_rho = mp.mpf(1) / 2
cos_alpha = 1 / mp.sqrt(5)
sin_alpha = mp.sqrt(1 - cos_alpha**2)

def boost_n(nx, ny, nz, ch, sh):
    """4x4 Lorentz boost in mpmath, axis (nx, ny, nz) unit, cosh = ch, sinh = sh."""
    M = mp.matrix(4, 4)
    M[0, 0] = ch
    for i, ni in enumerate([nx, ny, nz]):
        M[0, i+1] = sh * ni
        M[i+1, 0] = sh * ni
    n = [nx, ny, nz]
    for i in range(3):
        for j in range(3):
            d = mp.mpf(1) if i == j else mp.mpf(0)
            M[i+1, j+1] = d + (ch - 1) * n[i] * n[j]
    return M

B_J = boost_n(0, 0, 1, cosh_rho, sinh_rho)
B_g = boost_n(sin_alpha, 0, cos_alpha, cosh_rho, sinh_rho)
B_J_inv = boost_n(0, 0, 1, cosh_rho, -sinh_rho)
B_g_inv = boost_n(sin_alpha, 0, cos_alpha, cosh_rho, -sinh_rho)

C = B_J * B_g * B_J_inv * B_g_inv

# Trace and rotation cosine
trace_C = C[0, 0] + C[1, 1] + C[2, 2] + C[3, 3]
cos_omega_num = (trace_C - 2) / 2

# Verify Lorentz invariance: C^T eta C = eta
eta = mp.matrix([[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]])
gram = C.T * eta * C - eta
gram_err = max(abs(gram[i,j]) for i in range(4) for j in range(4))

# Closed form: cos(omega) = (31 + 8 sqrt 5) / 50
closed_form = (mp.mpf(31) + mp.mpf(8) * mp.sqrt(5)) / mp.mpf(50)
diff_1 = abs(cos_omega_num - closed_form)

print(f"[1]  Elliptic commutator rotation cosine")
print(f"     Tr(C)            = {mp.nstr(trace_C, 30)}")
print(f"     cos(omega)        = {mp.nstr(cos_omega_num, 30)}")
print(f"     closed form       = {mp.nstr(closed_form, 30)}")
print(f"     |measured - closed| = {mp.nstr(diff_1, 5)}")
print(f"     |C^T eta C - eta|   = {mp.nstr(gram_err, 5)}  (Lorentz check)")
pass_1 = (diff_1 < mp.mpf(10)**(-150)) and (gram_err < mp.mpf(10)**(-150))
print(f"     ==> {'PASS' if pass_1 else 'FAIL'}")
print()


# ============================================================================
# Section 2.  Minimal polynomial of cos(omega) over Q
# ============================================================================
#
# Verify symbolically that the closed form
#
#   cos(omega) = (31 + 8 sqrt 5) / 50
#
# satisfies the irreducible quadratic
#
#   2500 x^2 - 3100 x + 641 = 0
#
# over Q.  Its discriminant 800^2 * 5 is non-square in Q, so cos(omega)
# has degree 2 over Q and lies in Q(sqrt 5).

sqrt5 = sp.sqrt(5)
cw_sym = sp.Rational(31, 50) + sp.Rational(8, 50) * sqrt5

a, b, c = 2500, -3100, 641
poly_value = sp.simplify(a * cw_sym**2 + b * cw_sym + c)
disc = b**2 - 4 * a * c

print(f"[2]  Minimal polynomial of cos(omega)")
print(f"     polynomial:         {a} x^2 + ({b}) x + {c}")
print(f"     discriminant:       {disc} = 800^2 * 5  (non-square in Q)")
print(f"     poly evaluated at closed form: {poly_value}")
pass_2 = (poly_value == 0) and (disc == 800**2 * 5)
print(f"     ==> {'PASS' if pass_2 else 'FAIL'}")
print()


# ============================================================================
# Section 3.  Cyclotomic cosine exclusion: omega/pi is irrational
# ============================================================================
#
# If omega = p pi / q in lowest terms with q >= 3, then 2 cos omega lives
# in the maximal real subfield of Q(zeta_{2q}), of degree phi(2q)/2 over
# Q.  For c_omega to have degree 2, phi(2q) = 4 forces q in {4, 5, 6}.
# The 8 candidate cosines are:
#   +/- sqrt(2)/2          (q = 4),
#   +/- sqrt(3)/2          (q = 6),
#   +/- phi/2 = +/- (1+sqrt 5)/4    (q = 5),
#   +/- (phi-1)/2 = +/- (sqrt 5 - 1)/4    (q = 5).
# Cases q = 1, 2, 3 give degree-1 cosines (+/-1, 0, +/-1/2) and are
# excluded since c_omega has degree 2.  We check numerical distinctness
# of c_omega from each candidate.

candidates = {
    "+sqrt(2)/2": mp.sqrt(2) / 2,
    "-sqrt(2)/2": -mp.sqrt(2) / 2,
    "+sqrt(3)/2": mp.sqrt(3) / 2,
    "-sqrt(3)/2": -mp.sqrt(3) / 2,
    "+phi/2":     phi_n / 2,
    "-phi/2":    -phi_n / 2,
    "+(phi-1)/2": (phi_n - 1) / 2,
    "-(phi-1)/2": -(phi_n - 1) / 2,
}

print(f"[3]  Cyclotomic cosine exclusion of rational omega/pi")
all_excluded = True
for name, val in candidates.items():
    diff = abs(closed_form - val)
    excluded = diff > mp.mpf(10)**(-50)
    all_excluded = all_excluded and excluded
    print(f"     {name:<14} diff = {mp.nstr(diff, 5)}  {'excluded' if excluded else '*** MATCH ***'}")
pass_3 = all_excluded
print(f"     ==> {'PASS' if pass_3 else 'FAIL'}")
print()


# ============================================================================
# Section 4.  Molien series for I = A_5 acting on R^3
# ============================================================================
#
# Klein 1884:  M_I(t) = (1 + t^15) / ((1 - t^2)(1 - t^6)(1 - t^{10})).
# Expand as power series and read off coefficients.

t = sp.symbols('t')
M_series = (1 + t**15) / ((1 - t**2) * (1 - t**6) * (1 - t**10))
expanded = sp.series(M_series, t, 0, 17).removeO()
expanded_poly = sp.Poly(expanded, t)

molien_dims = {d: int(expanded_poly.nth(d)) for d in range(17)}

# Klein 1884 prediction for d = 0, 2, 4, 6, 8, 10, 12, 14, 15, 16:
# dim((S^d R^3)^I) = 1, 1, 1, 2, 2, 3, 4, 4, 1, 5
expected_dims = {
    0: 1, 2: 1, 4: 1, 6: 2, 8: 2, 10: 3, 12: 4, 14: 4, 15: 1, 16: 5,
    1: 0, 3: 0, 5: 0, 7: 0, 9: 0, 11: 0, 13: 0,
}

print(f"[4]  Molien series for I on R^3")
print(f"     M_I(t) = (1 + t^15) / ((1 - t^2)(1 - t^6)(1 - t^10))")
print(f"     coefficients d = 0..16:")
for d in range(17):
    star = "" if molien_dims[d] == expected_dims[d] else "  *** MISMATCH ***"
    print(f"       a_{d:<2} = {molien_dims[d]}   (expected {expected_dims[d]}){star}")
pass_4 = all(molien_dims[d] == expected_dims[d] for d in range(17))
print(f"     ==> {'PASS' if pass_4 else 'FAIL'}")
print()

# Standalone check of Klein's invariant degrees (the dimensions for d = 2, 4, 6
# determine that anisotropy starts at degree 6).
key_dims_correct = (
    molien_dims[2] == 1 and       # only |k|^2 invariant
    molien_dims[4] == 1 and       # only (|k|^2)^2 invariant
    molien_dims[6] == 2           # |k|^6 plus Klein form K_6
)
print(f"[5]  Anisotropy degree of I-invariants on R^3")
print(f"     dim((S^2 R^3)^I) = {molien_dims[2]}  (only |k|^2)")
print(f"     dim((S^4 R^3)^I) = {molien_dims[4]}  (only (|k|^2)^2: leading isotropy preserved)")
print(f"     dim((S^6 R^3)^I) = {molien_dims[6]}  (|k|^6 plus Klein form K_6: first anisotropy)")
pass_5 = key_dims_correct
print(f"     ==> {'PASS' if pass_5 else 'FAIL'}")
print()


# ============================================================================
# Section 6.  Free dispersion exponents: cubic |k|^4 vs icosahedral |k|^6
# ============================================================================
#
# E(k) = sum_e (1 - cos(k . e)).  Sample 5000 directions on the unit sphere
# (Fibonacci spiral); fit the angular spread to a power law in |k|.
# This is a separate spatial-symmetry check; it does NOT directly verify
# the Lorentz protection corollary of the paper.  It verifies that
# I-invariant spatial couplings have first anisotropy at degree 6,
# whereas O_h-invariant cubic couplings have anisotropy already at degree 4.

phi_np = (1 + 5**0.5) / 2
cubic_couplings = np.array([
    [1, 0, 0], [-1, 0, 0],
    [0, 1, 0], [0, -1, 0],
    [0, 0, 1], [0, 0, -1],
], dtype=float)
cubic_couplings /= np.linalg.norm(cubic_couplings, axis=1, keepdims=True)

icos_couplings = np.array([
    [0,  1,  phi_np], [0,  1, -phi_np], [0, -1,  phi_np], [0, -1, -phi_np],
    [1,  phi_np, 0],  [1, -phi_np, 0],  [-1, phi_np, 0],  [-1, -phi_np, 0],
    [phi_np, 0,  1],  [phi_np, 0, -1],  [-phi_np, 0, 1],  [-phi_np, 0, -1],
], dtype=float)
icos_couplings /= np.linalg.norm(icos_couplings, axis=1, keepdims=True)

N_dirs = 5000
indices = np.arange(0, N_dirs, dtype=float) + 0.5
phi_angle = np.arccos(1 - 2 * indices / N_dirs)
theta_angle = np.pi * (1 + 5**0.5) * indices
dirs = np.stack([
    np.cos(theta_angle) * np.sin(phi_angle),
    np.sin(theta_angle) * np.sin(phi_angle),
    np.cos(phi_angle),
], axis=-1)

k_values = np.array([0.05, 0.1, 0.15, 0.2, 0.3])
cubic_spreads = []
icos_spreads = []
for kmag in k_values:
    kv = kmag * dirs
    Ec = np.sum(1.0 - np.cos(kv @ cubic_couplings.T), axis=-1)
    Ei = np.sum(1.0 - np.cos(kv @ icos_couplings.T), axis=-1)
    cubic_spreads.append(Ec.max() - Ec.min())
    icos_spreads.append(Ei.max() - Ei.min())
cubic_spreads = np.array(cubic_spreads)
icos_spreads = np.array(icos_spreads)

cubic_p, cubic_a = np.polyfit(np.log(k_values), np.log(cubic_spreads), 1)
icos_p, icos_a = np.polyfit(np.log(k_values), np.log(icos_spreads), 1)

print(f"[6]  Cubic coupling: dispersion anisotropy spread ~ |k|^4")
print(f"     fit: spread = {np.exp(cubic_a):.4f} . |k|^{cubic_p:.4f}")
print(f"     predicted exponent: 4  (Molien dim((S^4 R^3)^O_h) = 2)")
pass_6 = abs(cubic_p - 4.0) < 0.05
print(f"     ==> {'PASS' if pass_6 else 'FAIL'}")
print()

print(f"[7]  Icosahedral coupling: dispersion anisotropy spread ~ |k|^6")
print(f"     fit: spread = {np.exp(icos_a):.4e} . |k|^{icos_p:.4f}")
print(f"     predicted exponent: 6  (Molien dim((S^6 R^3)^I) = 2)")
pass_7 = abs(icos_p - 6.0) < 0.05
print(f"     ==> {'PASS' if pass_7 else 'FAIL'}")
print()


# ============================================================================
# Section 8.  Galois equivariance of M_J
# ============================================================================
#
# Build M_zeta as the companion matrix of zeta_5 acting on Z[zeta_5] in basis
# (1, zeta, zeta^2, zeta^3) using zeta^4 = -1 - zeta - zeta^2 - zeta^3.
# Build sigma_k for k = 2, 3, 4 and verify the substrate identities:
#   M_zeta^5 == I
#   sigma_k . M_J . sigma_k^{-1} == M_{sigma_k(J)} = I + M_zeta^{2k}.
#
# This is the algebraic foundation for the Galois-equivariant choice of J
# and the conjugate-magnitude pattern (1/phi, phi, phi, 1/phi).

M_zeta = sp.Matrix([
    [0, 0, 0, -1],
    [1, 0, 0, -1],
    [0, 1, 0, -1],
    [0, 0, 1, -1],
])
I4 = sp.eye(4)
M_J_sym = I4 + M_zeta**2

def reduce_zeta(j):
    """Express zeta^j in basis (1, zeta, zeta^2, zeta^3)."""
    j = j % 5
    if j < 4:
        v = [0, 0, 0, 0]
        v[j] = 1
        return v
    return [-1, -1, -1, -1]

def galois_sigma(k):
    """Matrix of sigma_k: zeta -> zeta^k in the basis (1, zeta, zeta^2, zeta^3)."""
    M = sp.zeros(4, 4)
    for col in range(4):
        v = reduce_zeta(col * k)
        for row in range(4):
            M[row, col] = v[row]
    return M

# Check 1: M_zeta^5 = I
M_zeta_5 = M_zeta**5
zeta5_eq_I = (M_zeta_5 == I4)

# Check 2: sigma_k . M_J . sigma_k^{-1} = I + M_zeta^{2k} for k = 2, 3, 4
galois_checks = []
for k in [2, 3, 4]:
    sk = galois_sigma(k)
    sk_inv = sk**-1
    LHS = sk * M_J_sym * sk_inv
    RHS = I4 + M_zeta**(2 * k)
    galois_checks.append(((LHS - RHS) == sp.zeros(4, 4), k))

print(f"[8]  Galois equivariance:  sigma_k . M_J . sigma_k^{{-1}} == M_{{sigma_k(J)}}")
print(f"     M_zeta^5 == I:                              {zeta5_eq_I}")
for ok, k in galois_checks:
    print(f"     sigma_{k} . M_J . sigma_{k}^-1 == I + M_zeta^{2*k}: {ok}")
pass_8 = zeta5_eq_I and all(ok for ok, _ in galois_checks)
print(f"     ==> {'PASS' if pass_8 else 'FAIL'}")
print()


# ============================================================================
# Section 9.  Elliptic classification of C: rank, fixed plane signature
# ============================================================================
#
# Eigenvalues alone do not classify a Lorentz transformation; one must
# also verify that the eigenvalue 1 has no Jordan block (otherwise C
# would be parabolic, not purely elliptic).  We check three conditions
# in exact arithmetic over Q(sqrt 5):
#
#   (a) C is Lorentz: C^T eta C = eta
#   (b) eigenvalue 1 has 2-dimensional fixed space: rank(C - I) = 2
#   (c) restriction of eta to ker(C - I) is non-degenerate Lorentzian:
#       the 2x2 Gram matrix has signature (1, 1), equivalently negative
#       determinant.
#
# Together these prove C is conjugate in SO+(3,1)^0 to a spatial rotation
# by angle omega about the spacelike axis perpendicular to the fixed
# (timelike + spacelike) plane.  This is the foundation for the closure
# of <C> being a one-parameter elliptic subgroup.

import sympy as sp

s5_sym = sp.sqrt(5)
ch_sym = s5_sym / 2
sh_sym = sp.Rational(1, 2)
ca_sym = 1 / s5_sym
sa_sym = sp.simplify(sp.sqrt(1 - ca_sym**2))

def boost_sym(nx, ny, nz, ch, sh):
    """4x4 sympy boost matrix in (t, x, y, z), axis (nx, ny, nz) unit."""
    M = sp.zeros(4, 4)
    M[0, 0] = ch
    n = [nx, ny, nz]
    for i in range(3):
        M[0, i+1] = sh * n[i]
        M[i+1, 0] = sh * n[i]
    for i in range(3):
        for j in range(3):
            d = 1 if i == j else 0
            M[i+1, j+1] = d + (ch - 1) * n[i] * n[j]
    return M

B1_sym = boost_sym(0, 0, 1, ch_sym, sh_sym)
B2_sym = boost_sym(sa_sym, 0, ca_sym, ch_sym, sh_sym)
B1i_sym = boost_sym(0, 0, 1, ch_sym, -sh_sym)
B2i_sym = boost_sym(sa_sym, 0, ca_sym, ch_sym, -sh_sym)
C_sym = sp.simplify(B1_sym * B2_sym * B1i_sym * B2i_sym)

eta_sym = sp.diag(1, -1, -1, -1)

# Test (a): Lorentz
gram_C = sp.simplify(C_sym.T * eta_sym * C_sym - eta_sym)
is_lorentz_sym = (gram_C == sp.zeros(4, 4))

# Test (a'): det = 1
det_C_sym = sp.simplify(C_sym.det())
is_det_one = (det_C_sym == 1)

# Test (b): rank(C - I) = 2
CmI = sp.simplify(C_sym - sp.eye(4))
rank_CmI = CmI.rank()

# Test (c): fixed plane Gram has signature (1, 1)
fixed = CmI.nullspace()
G_fixed = sp.Matrix([
    [sp.simplify((fixed[i].T * eta_sym * fixed[j])[0, 0]) for j in range(len(fixed))]
    for i in range(len(fixed))
])
det_G_sym = sp.simplify(G_fixed.det())
is_signature_11 = (sp.N(det_G_sym) < 0)

print(f"[9]  Elliptic classification of C")
print(f"     C^T eta C - eta == 0:    {is_lorentz_sym}")
print(f"     det(C) == 1:             {is_det_one}")
print(f"     rank(C - I) == 2:        {rank_CmI == 2}  (rank = {rank_CmI})")
print(f"     dim ker(C - I) == 2:     {len(fixed) == 2}")
print(f"     fixed-plane Gram det:    {det_G_sym}")
print(f"     signature (1, 1):        {is_signature_11}  (det < 0)")
pass_9 = is_lorentz_sym and is_det_one and (rank_CmI == 2) and (len(fixed) == 2) and is_signature_11
print(f"     ==> {'PASS' if pass_9 else 'FAIL'}")
print()


# ============================================================================
# Section 10.  Explicit A_5 element sending v_1 -> v_2
# ============================================================================
#
# Use the standard icosahedron vertices on R^3:
#
#   v_1 = (0, 1, phi),    v_2 = (1, phi, 0).
#
# The explicit order-three rotation
#
#       (-1/2          (1-sqrt5)/4   (1+sqrt5)/4 )
#   g = ( (-1+sqrt5)/4  (1+sqrt5)/4  1/2          )
#       ( -(1+sqrt5)/4  1/2          (1-sqrt5)/4 )
#
# satisfies g in SO(3) (rotation, det = 1, orthogonal), g^3 = I,
# g v_1 = v_2.  The angle alpha between v_1 and v_2 satisfies
# cos alpha = (v_1 . v_2) / (||v_1|| ||v_2||) = 1/sqrt(5).  This
# realizes the abstract construction in section 1 with a concrete
# A_5 element, removing any ambiguity about whether the chosen
# axis configuration lives in A_5 at all.

g_A5 = sp.Matrix([
    [sp.Rational(-1, 2),     (1 - s5_sym) / 4,         (1 + s5_sym) / 4],
    [(-1 + s5_sym) / 4,      (1 + s5_sym) / 4,         sp.Rational(1, 2)],
    [-(1 + s5_sym) / 4,      sp.Rational(1, 2),        (1 - s5_sym) / 4],
])

phi_sym = (1 + s5_sym) / 2
v1 = sp.Matrix([0, 1, phi_sym])
v2 = sp.Matrix([1, phi_sym, 0])

# Check orthogonal: g^T g = I
gTg = sp.simplify(g_A5.T * g_A5)
is_orthogonal = (gTg == sp.eye(3))

# Check det g = 1 (rotation)
det_g = sp.simplify(g_A5.det())
is_rotation = (det_g == 1)

# Check g^3 = I (order 3)
g3 = sp.simplify(g_A5 ** 3)
is_order_3 = (g3 == sp.eye(3))

# Check g^2 != I (so order is exactly 3)
g2 = sp.simplify(g_A5 ** 2)
is_not_order_2 = (g2 != sp.eye(3))

# Check g v_1 = v_2
gv1 = sp.simplify(g_A5 * v1)
sends_v1_to_v2 = (sp.simplify(gv1 - v2) == sp.zeros(3, 1))

# Verify cos alpha = 1/sqrt(5)
norm_sq = sp.simplify(v1.dot(v1))
inner = sp.simplify(v1.dot(v2))
# norm_sq = phi + 2; inner = phi
norm_sq_expected = phi_sym + 2
inner_expected = phi_sym
norm_sq_correct = (sp.simplify(norm_sq - norm_sq_expected) == 0)
inner_correct = (sp.simplify(inner - inner_expected) == 0)

cos_alpha_real = sp.simplify(inner / sp.sqrt(norm_sq * sp.simplify(v2.dot(v2))))
expected_cos = 1 / s5_sym
cos_alpha_correct = (sp.simplify(cos_alpha_real - expected_cos) == 0)

print(f"[10] Explicit A_5 element g sending v_1 -> v_2")
print(f"     g^T g == I:                          {is_orthogonal}")
print(f"     det(g) == 1:                         {is_rotation}")
print(f"     g^3 == I (order 3):                  {is_order_3}")
print(f"     g^2 != I (so order is exactly 3):    {is_not_order_2}")
print(f"     g v_1 == v_2:                        {sends_v1_to_v2}")
print(f"     ||v_1||^2 = ||v_2||^2 = phi + 2:     {norm_sq_correct}")
print(f"     v_1 . v_2 == phi:                    {inner_correct}")
print(f"     cos alpha == 1/sqrt(5):              {cos_alpha_correct}")
pass_10 = (is_orthogonal and is_rotation and is_order_3 and is_not_order_2
           and sends_v1_to_v2 and norm_sq_correct and inner_correct
           and cos_alpha_correct)
print(f"     ==> {'PASS' if pass_10 else 'FAIL'}")
print()


# ============================================================================
# Section 11.  Group membership: g permutes 12 icosahedron vertices
# ============================================================================
#
# The conditions g^T g = I, det g = 1, g^3 = I, g v_1 = v_2 from test [10]
# are necessary for g in A_5 but not sufficient.  A reviewer can ask:
# "why is this matrix in I = A_5, not merely in SO(3)?"  The clean answer
# is to verify that g preserves the set of 12 standard icosahedron
# vertices.  The rotational symmetry group of the standard icosahedron is
# precisely I = A_5, so g . V = V (as a set) implies g in A_5.
#
# Standard icosahedron vertex set:
#   (0, +/- 1, +/- phi),
#   (+/- 1, +/- phi, 0),
#   (+/- phi, 0, +/- 1).
#
# Twelve vertices in total.  An order-3 element of A_5 acts on V with
# orbit structure 4 x 3-cycles (no fixed vertices), since A_5 acts
# transitively on V with stabilizer C_5.

vertex_set = []
for s1 in [1, -1]:
    for s2 in [1, -1]:
        vertex_set.append(sp.Matrix([0, s1, s2 * phi_sym]))
        vertex_set.append(sp.Matrix([s1, s2 * phi_sym, 0]))
        vertex_set.append(sp.Matrix([s1 * phi_sym, 0, s2]))

def canonical_vertex(v):
    """Canonicalize a Matrix vertex to a hashable tuple of simplified entries."""
    return tuple(sp.simplify(v[i, 0]) for i in range(3))

V_set = set(canonical_vertex(v) for v in vertex_set)
gV_set = set(canonical_vertex(sp.simplify(g_A5 * v)) for v in vertex_set)

set_equal = (V_set == gV_set)
n_vertices = len(V_set)
n_images = len(gV_set)

# Verify orbit structure: |V| = 12, all orbits should have length dividing 3
orbits = []
remaining = set(V_set)
while remaining:
    v_canon = next(iter(remaining))
    orbit = [v_canon]
    v_mat = sp.Matrix([sp.sympify(x) for x in v_canon])
    for _ in range(2):
        v_mat = sp.simplify(g_A5 * v_mat)
        c = canonical_vertex(v_mat)
        if c == v_canon:
            break
        orbit.append(c)
    for x in orbit:
        remaining.discard(x)
    orbits.append(orbit)

orbit_sizes = sorted([len(o) for o in orbits])
expected_sizes = [3, 3, 3, 3]
orbit_correct = (orbit_sizes == expected_sizes)

print(f"[11] Group membership of g via vertex permutation")
print(f"     |V| (12 standard icosahedron vertices):    {n_vertices}")
print(f"     |g . V|:                                    {n_images}")
print(f"     g permutes V (g . V == V):                  {set_equal}")
print(f"     orbit sizes (expect [3, 3, 3, 3]):          {orbit_sizes}")
print(f"     ==> {'PASS' if (set_equal and n_vertices == 12 and orbit_correct) else 'FAIL'}")
pass_11 = set_equal and n_vertices == 12 and orbit_correct
print()


# ============================================================================
# Final tally
# ============================================================================

results = [
    ("[1]  elliptic commutator cos(omega)",     pass_1),
    ("[2]  minimal polynomial of cos(omega)",   pass_2),
    ("[3]  cyclotomic cosine exclusion",        pass_3),
    ("[4]  Molien series coefficients",         pass_4),
    ("[5]  anisotropy degree of I-invariants",  pass_5),
    ("[6]  cubic dispersion exponent ~ 4",      pass_6),
    ("[7]  icos dispersion exponent ~ 6",       pass_7),
    ("[8]  Galois equivariance of M_J",         pass_8),
    ("[9]  elliptic classification of C",       pass_9),
    ("[10] explicit A_5 element",               pass_10),
    ("[11] g permutes 12 icosahedron vertices", pass_11),
]
n_pass = sum(1 for _, ok in results if ok)
n_total = len(results)

print("=" * 72)
for label, ok in results:
    print(f"  {label:<46}  {'PASS' if ok else 'FAIL'}")
print("=" * 72)
print(f"  {n_pass}/{n_total} PASS")
print()
