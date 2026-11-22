"""Exact rational certificate for H_p'(1/17)>0, 3/2<=p<2.

Only Fraction arithmetic is used.  The logarithms are enclosed by
the positive atanh series with an explicit geometric remainder.
"""
from fractions import Fraction as F
from math import factorial


def log_from_atanh(z, terms=24):
    lo = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    hi = lo + 2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
    return lo, hi


alpha = log_from_atanh(F(1,33))  # log(17/16)
log2 = log_from_atanh(F(1,3))
beta = (4*log2[0]+alpha[0], 4*log2[1]+alpha[1])  # log 17
terms = [
    (F(4928,4913), -F(288,4913), alpha),
    (F(9536,83521), -F(288,4913),
     (alpha[0]+beta[0],alpha[1]+beta[1])),
    (-F(576,4913), F(288,4913), beta),
    (F(1,83521), F(0), (2*beta[0],2*beta[1])),
]


def coefficient_interval(j):
    lo = hi = F(0)
    for a,b,(l,u) in terms:
        for c,k in ((a,j),(b,j-1)):
            x,y=l**k/F(factorial(k)),u**k/F(factorial(k))
            lo += c*(x if c>=0 else y)
            hi += c*(y if c>=0 else x)
    return lo,hi


assert sum(a for a,b,interval in terms)==1
assert coefficient_interval(1)[0]>F(49,100000)
for j in range(2,6):
    assert coefficient_interval(j)[0]>0
for j in range(6,12):
    assert coefficient_interval(j)[0]>-F(1,1000)

# For 0<t<=1/2, the omitted degree>=12 part divided by t.
simple_bounds=[F(1,16),F(3),F(3),F(6)]
tail=F(0)
for (a,b,interval),L in zip(terms,simple_bounds):
    assert interval[1]<L
    tail += abs(a)*L**12/F(factorial(12))*F(1,2)**11/(1-L/F(26))
    tail += abs(b)*L**11/F(factorial(11))*F(1,2)**11/(1-L/F(24))
assert tail<F(1,1000000)
margin=F(49,100000)-F(1,16000)-F(1,1000000)
assert margin>0
print('PASS: all exact logarithmic, Taylor-coefficient, and tail bounds.')
print('Certified lower bound J(t)/t >',margin)

