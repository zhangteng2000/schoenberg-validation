"""Exact symbolic reconstruction of cubic-tail polynomial certificates.

Usage from the manuscript directory:
    python certificates/verify_odd_ge_three.py
    python certificates/verify_odd_ge_three.py --smoke

Requires SymPy. Reads only odd_ge_3_allocations.json, whose values are the
43 allocation rows printed in the appendix. All inequalities are checked
as exact coefficient identities, without floating-point sampling.
"""
import json
import sys
from pathlib import Path
from math import comb
import sympy as S
from sympy.parsing.sympy_parser import (
    parse_expr, standard_transformations, implicit_multiplication_application
)

K, H, v = S.symbols('K H v', nonnegative=True)
transforms = standard_transformations + (implicit_multiplication_application,)


def parse_allocation(expr):
    return parse_expr(expr, local_dict={'K': K, 'H': H},
                      transformations=transforms, evaluate=True)


def is_bernstein_nonnegative(expr, label, degree=None):
    """Check all K,H monomial coefficients of every v Bernstein coefficient."""
    p = S.Poly(S.expand(expr), v)
    if p.is_zero:
        return 0
    degree = p.degree() if degree is None else degree
    assert degree >= p.degree(), label
    # The normalized Bernstein coefficients are the following sums.
    count = 0
    for j in range(degree + 1):
        beta = sum(S.Rational(comb(j, i), comb(degree, i)) * p.nth(i)
                   for i in range(j + 1))
        q = S.Poly(S.expand(beta), K, H)
        for mon, coeff in q.terms():
            assert coeff >= 0, (label, j, mon, coeff)
            if coeff > 0:
                count += 1
    return count


def unbalanced(sign, allocation, chosen=None):
    k, h, m = K+1, H+1, K+H+2
    n, ell = 2*m+1, K+2*H+3
    tau = k*v if sign == '+' else -ell*v
    A, B, C = n*ell+tau, (n-1)*tau, n*k-tau
    SS = n*(k-ell)+(n-3)*tau
    DD = ((n-1)*tau-n*(k-ell))**2 + 4*n*k*ell
    counts = 0
    processed = 0
    for jstr, row in allocation.items():
        j = int(jstr)
        if chosen is not None and j not in chosen:
            continue
        alpha = [parse_allocation(s) for s in row]
        e = [(j >> a) & 1 for a in range(6)]
        gamma_c = e[0]*h*(m+1)**3 + e[1]*(h-1)*m**3
        gamma_u = e[2]*m*(m+1)**3 + e[3]*(m+1)*m**3
        gamma_v = e[4]*m*(m+1)**3 + e[5]*(m+1)*m**3
        q0 = (A**3 * (alpha[0]*(m+1)**3 + alpha[1]*m**3 + alpha[2])
              + (1 if sign == '+' else -1)*B**3*
              (alpha[3]*(m+1)**3 + alpha[4]*m**3 + alpha[5])
              + C**3*(alpha[6] - gamma_c))
        if gamma_u == 0 and gamma_v == 0:
            counts += is_bernstein_nonnegative(q0, (sign, j, 'Q0'))
            processed += 1
            continue
        rp = 8*q0 - (gamma_v-gamma_u)*SS*(SS**2+3*DD)
        tp = (gamma_u+gamma_v)*(DD+3*SS**2)
        ep = rp**2 - DD*tp**2
        counts += is_bernstein_nonnegative(rp, (sign, j, 'R+'))
        counts += is_bernstein_nonnegative(ep, (sign, j, 'E+'))
        processed += 2
        if sign == '-' and gamma_v != 0:
            rm = 8*q0 + (gamma_u+gamma_v)*SS*(SS**2+3*DD)
            counts += is_bernstein_nonnegative(
                rm, (sign, j, 'R-'),
                degree=4 if j == 20 else None)
            processed += 1
            if S.expand(gamma_u-gamma_v) != 0:
                tm = (gamma_u-gamma_v)*(DD+3*SS**2)
                em = rm**2 - DD*tm**2
                counts += is_bernstein_nonnegative(em, (sign, j, 'E-'))
                processed += 1
    return processed, counts


def balanced():
    m = K + 2
    A, C, B = 2*m+1+v, 2*m+1-v, 2*m*v
    ss = (m-1)*v
    dd = m**2*v**2+2*m+1
    a0, b0 = m*(m+1)**3, (m+1)*m**3
    p3 = m**3*(A**3+C**3) + ((m-1)*(m+1)**3+m**4)*B**3
    rows = [
        (m**3*A**3+(m-1)*(m+1)**3*B**3, a0, 0),
        (m**3*(A**3+C**3)+(m-1)*((m+1)**3+m**3)*B**3, a0, a0),
        (p3, a0+b0, 0),
        (p3+m*A**3, a0+b0, a0),
        ((m**3+m)*(A**3+C**3)+((m-1)*(m+1)**3+m**4+1)*B**3,
         a0+b0, a0+b0)
    ]
    total = 0
    for j, (P, c_v, c_u) in enumerate(rows):
        RR = P - (c_v-c_u)*ss*(3*dd+ss**2)
        TT = (c_v+c_u)*(dd+3*ss**2)
        EE = RR**2 - dd*TT**2
        total += is_bernstein_nonnegative(RR, ('balanced', j, 'R'))
        total += is_bernstein_nonnegative(EE, ('balanced', j, 'E'))
    return 10, total


def main():
    data = json.loads((Path(__file__).parent/'odd_ge_3_allocations.json').read_text())
    assert len(data['+']) == 26 and len(data['-']) == 17
    if '--smoke' in sys.argv:
        for sign in ['+', '-']:
            print(sign, unbalanced(sign, data[sign], {1, 3}), flush=True)
        print('PASS smoke', flush=True)
    else:
        pos = unbalanced('+', data['+'])
        print('Nonnegative-middle:', pos, flush=True)
        assert pos[0] == 50
        neg = unbalanced('-', data['-'])
        print('Nonpositive-middle:', neg, flush=True)
        assert neg[0] == 44
        bal = balanced()
        print('Balanced:', bal, flush=True)
        assert bal[0] == 10
        print('PASS 104 exact identities', flush=True)


if __name__ == '__main__':
    main()
