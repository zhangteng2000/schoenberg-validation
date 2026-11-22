"""Exact rational certificate: H_p'(1/18)>0 for 37/25<=p<2.
Together with r_p<q^2/(4+q^2) for q<=12/25, this proves
r_p<1/18 and theta_p<1/19 for all 1<p<2.
All arithmetic below is rational; no floating-point assertion is used.
"""
from fractions import Fraction as F
from math import factorial


def log_from_atanh(z, terms=24):
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
    hi=lo+2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
    return lo,hi

r=F(1,18)
alpha=log_from_atanh(F(1,35)) # log(18/17)
log2=log_from_atanh(F(1,3))
log98=log_from_atanh(F(1,17))
beta=(4*log2[0]+log98[0],4*log2[1]+log98[1]) # log18
terms=[
 ((1-r)*(1+r+2*r*r),-(1-r)*(r+r*r),alpha),
 ((1-r)*r*(2+r+r*r),-(1-r)*r*(1+r),(alpha[0]+beta[0],alpha[1]+beta[1])),
 (-2*r*(1-r*r),r*(1-r*r),beta),
 (r**4,F(0),(2*beta[0],2*beta[1])),
]

def coefficient_interval(j):
    lo=hi=F(0)
    for a,b,(l,u) in terms:
        for cc,k in ((a,j),(b,j-1)):
            x,y=l**k/F(factorial(k)),u**k/F(factorial(k))
            lo+=cc*(x if cc>=0 else y)
            hi+=cc*(y if cc>=0 else x)
    return lo,hi

assert sum(a for a,b,iv in terms)==1
assert coefficient_interval(1)[0]>F(1,5000)
for j in range(2,6):
    assert coefficient_interval(j)[0]>0
for j in range(6,12):
    assert coefficient_interval(j)[0]>-F(1,1000)

T=F(13,25)
tail=F(0)
for (a,b,iv),L in zip(terms,[F(1,17),F(3),F(3),F(6)]):
    assert iv[1]<L
    tail+=abs(a)*L**12/F(factorial(12))*T**11/(1-L*T/F(13))
    tail+=abs(b)*L**11/F(factorial(11))*T**11/(1-L*T/F(12))
assert tail<F(1,500000)
margin=F(1,5000)-F(1,1000)*T**5/(1-T)-F(1,500000)
assert margin>F(1,10000)
assert F(12,25)**2/(4+F(12,25)**2)<r
print('PASS: exact Taylor bounds imply H_p\'(1/18)>0 for 37/25<=p<2.')
print('Certified lower bound J(1-q)/(1-q) >',margin)
print('Together with the existing small-q bound: r_p<1/18, theta_p<1/19 for every 1<p<2.')

