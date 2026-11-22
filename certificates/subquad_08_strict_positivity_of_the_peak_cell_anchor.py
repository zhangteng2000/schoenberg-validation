"""Exact rational checks for strict peak-cell left-anchor positivity.

The analytic proof is in turn3_norm_peak_anchor_positive.tex.
This script checks every nontrivial numerical constant there and the
finite exponential certificate for r_p<1/25 when 2/5<=q<=3/5.
No floating-point assertion or numerical sampling is used.
"""
from fractions import Fraction as F
from math import factorial


def log_atanh(z, terms=30):
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
    hi=lo+2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
    return lo,hi


# Elementary exponential/logarithmic comparisons used in the proof.
assert F(25,2)**3>F(9,2)**5
assert F(9,200)**3<F(2,5)**10
assert F(22,5)**5<F(12)**3
assert F(47,20)**5<F(17,2)**2
e_upper=F(8,3)+F(5,96)
assert e_upper<F(11,4)
assert F(11,4)**2<F(17,2)
assert sum((F(5,2)**j/F(factorial(j)) for j in range(7)),F(0))>12

assert F(1973,2000)*3-F(3,7)-1-F(23,20)-F(7,3)*F(9,200)-2*F(9,200)**2==F(38063,140000)>0
assert F(121,125)*3-F(2,3)-F(2,5)-F(6,5)-F(7,3)*F(2,25)-2*F(2,25)**2==F(821,1875)>0

for R,q,Y,target in [(F(1,12),F(2,5),F(22,5),F(571,1350)),(F(2,17),F(3,5),F(47,20),F(3421,57800))]:
    aa=1-q
    assert (q-2*R*R)/aa-R>0
    val=(q-2*R*R)*(Y-1)/aa-R*Y-(q+2)/2-(6+2*q-q*q)*R/3-2*R*R
    assert val==target>0
assert F(1,2)-F(10,3)*F(1,12)-6*F(1,12)**2==F(13,72)>0
assert F(1,2)-F(10,3)*F(2,17)-6*F(2,17)**2==F(43,1734)>0

# Middle exponent location certificate.
r=F(1,25)
alpha=log_atanh(F(1,49))
log2=log_atanh(F(1,3))
log54=log_atanh(F(1,9))
beta=(4*log2[0]+2*log54[0],4*log2[1]+2*log54[1])
terms=[
 ((1-r)*(1+r+2*r*r),-(1-r)*(r+r*r),alpha),
 ((1-r)*r*(2+r+r*r),-(1-r)*r*(1+r),(alpha[0]+beta[0],alpha[1]+beta[1])),
 (-2*r*(1-r*r),r*(1-r*r),beta),
 (r**4,F(0),(2*beta[0],2*beta[1])),
]


def coefficient_interval(j):
    lo=hi=F(0)
    for a,b,(l,u) in terms:
        for cc,k in [(a,j),(b,j-1)]:
            x,y=l**k/F(factorial(k)),u**k/F(factorial(k))
            lo+=cc*(x if cc>=0 else y)
            hi+=cc*(y if cc>=0 else x)
    return lo,hi


assert sum(a for a,b,iv in terms)==1
assert coefficient_interval(1)[0]>-F(3,5000)
assert coefficient_interval(2)[0]>0
assert coefficient_interval(3)[0]>F(17,5000)
assert coefficient_interval(4)[0]>F(3,1000)
for j in [5,6]:
    assert coefficient_interval(j)[0]>0
for j in range(7,12):
    assert coefficient_interval(j)[0]>-F(3,5000)
T=F(3,5)
tail=F(0)
for (a,b,iv),L in zip(terms,[F(1,24),F(10,3),F(13,4),F(13,2)]):
    assert iv[1]<L
    tail+=abs(a)*L**12/F(factorial(12))*T**11/(1-L*T/F(13))
    tail+=abs(b)*L**11/F(factorial(11))*T**11/(1-L*T/F(12))
assert tail<F(1,100000)
margin=-F(3,5000)+F(17,5000)*F(2,5)**2+F(3,1000)*F(2,5)**3-F(3,5000)*T**6/(1-T)-F(1,100000)
assert margin==F(3501,62500000)>0
print('PASS: all exact constants and the middle-location Taylor certificate.')
print('Certified J/(1-q) lower bound:',margin)

