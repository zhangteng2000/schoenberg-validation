"""Exact symbolic checks for the elementary negative-middle residue proof.

The analytic inequalities and domains are stated and proved in
turn2_negative_center_monotonicity_proof.tex.  This short checker verifies
only the polynomial/rational identities and endpoint constants used there;
it does not replace any mathematical step with numerical sampling.
"""

import sympy as sp

t, A, B, z, zdot, K = sp.symbols("t A B z zdot K")
lhs = -zdot * z * (t + A) * (t + B) / (K * t * (t + z))
rhs = (
    -z * zdot / K
    - A * B * zdot / (K * t)
    + zdot * (z - A) * (z - B) / (K * (t + z))
)
assert sp.cancel(lhs - rhs) == 0
print("Exact critical-residue partial fraction identity: PASS")

h, k, q, p = sp.symbols("h k q p")
n = h + k + 1
d = h - k
D0 = h * (k + q**2) + (k - q)**2
F = (p + 1) * D0 - 2 * h * q * (q - 1)
G = (k - q) * (n * p + h + k + k * q) - h**2 * q
bracket1 = (h + 1) * q**2 + (d + 1) * q - k
bracket2 = (h + k) * (d + 1) + 2 * d - (d - 1) * q
assert sp.expand(F - G - p * bracket1 - q * bracket2) == 0
assert sp.expand(bracket1.subs(q, 1) - 2 * (d + 1)) == 0
assert sp.expand(bracket2.subs(q, k) - h * (d + 3)) == 0
print("Exact F-G decomposition and its two positive endpoint bounds: PASS")

P, Q, w = sp.symbols("P Q w")
bound = w * (2 * w + A + B) * (A - w) * (B - w)
bound /= (A + w)**2 * (B + w)**2
sub = {A: w * (1 + P) / (1 - P), B: w * (1 + Q) / (1 - Q)}
assert sp.cancel(bound.subs(sub) - P * Q * (2 - P - Q) / 2) == 0
s = sp.symbols("s")
maxpoly = s**2 * (2 - s) / 8
assert maxpoly.subs(s, sp.Rational(4, 3)) == sp.Rational(4, 27)
print("Exact positive-critical 4/27 substitution and maximum: PASS")

assert sp.Rational(5, 7) - sp.Rational(1, 2) - sp.Rational(4, 27) == sp.Rational(25, 378)
assert sp.Rational(3, 5) - sp.Rational(4, 27) == sp.Rational(61, 135)
print("Exact strict monotonicity margins: 25/378 (n>=7), 61/135 (n=5)")
print("All elementary algebraic checks passed.")

