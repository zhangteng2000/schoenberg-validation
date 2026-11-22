"""Exact rational proof r_p > (p-1)^2/32, 1<p<2.
The only transcendental estimates use positive atanh log series.
"""
from fractions import Fraction as F
from math import factorial

assert (F(5,4)**8/F(32)) < F(3,4)**5
assert F(9,8)*F(3,4)+F(56,27)/32+F(1,32)==F(203,216)<1

def log_interval(z, terms=24):
    a=2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
    return a,a+2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
a=log_interval(F(1,63))
b0=log_interval(F(1,3)); b=(5*b0[0],5*b0[1])
r=F(1,32)
terms=[((1-r)*(1+r+2*r*r),-(1-r)*(r+r*r),a),
       ((1-r)*r*(2+r+r*r),-(1-r)*r*(1+r),(a[0]+b[0],a[1]+b[1])),
       (-2*r*(1-r*r),r*(1-r*r),b),
       (r**4,F(0),(2*b[0],2*b[1]))]
assert sum(x for x,y,I in terms)==1

def coeff(j):
    lo=hi=F(0)
    for x,y,(l,u) in terms:
        for c,k in [(x,j),(y,j-1)]:
            L=l**k/F(factorial(k)); U=u**k/F(factorial(k))
            lo+=c*(L if c>=0 else U)
            hi+=c*(U if c>=0 else L)
    return lo,hi
assert coeff(1)[1]<-F(1,1600)
assert coeff(2)[1]<0
for j in range(3,12): assert coeff(j)[1]<F(1,400)
# Bound the entire omitted degree >=12 tail, divided by e, 0<e<=1/5.
Lbounds=[F(1,31),F(4),F(4),F(7)]
tail=F(0); E=F(1,5)
for (x,y,I),L in zip(terms,Lbounds):
    assert I[1]<L
    tail+=abs(x)*L**12/F(factorial(12))*E**11/(1-L*E/F(13))
    tail+=abs(y)*L**11/F(factorial(11))*E**11/(1-L*E/F(12))
assert tail<F(1,1000000)
margin=F(1,1600)-F(1,400)*E**2/(1-E)-F(1,1000000)
assert margin>0
print('PASS: exact lower location r_p > q^2/32.')
print('For q<=4/5, K(q^2/32)/q^2 < -13/216.')
print('For q>=4/5, J(1/32)/(1-q) < -', margin)

