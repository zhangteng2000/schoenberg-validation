"""Exact rational polynomial certificate for the three peak derivative signs.

Reconstructs all expressions from scratch; no cached symbolic object,
floating-point computation, or numerical sampling is used.
Dependencies: SymPy, Python's fractions and math standard modules.
"""
import sympy as s
from fractions import Fraction as F
from math import comb, prod

r,q,X,Z=s.symbols('r q X Z')
u,v,w=s.symbols('u v w')


def derivative(f):
    return (1+r)**2*(s.diff(f,r)+q*X/r*s.diff(f,X)-q*Z/(1-r)*s.diff(f,Z))


def bernstein(poly, variables):
    pp=s.Poly(poly,*variables)
    degrees=pp.degree_list()
    coeff={mo:F(int(co.p),int(co.q)) for mo,co in pp.terms()}
    for axis,degree in enumerate(degrees):
        transformed={}
        for mo,co in coeff.items():
            i=mo[axis]
            for h in range(i,degree+1):
                target=list(mo)
                target[axis]=h
                target=tuple(target)
                transformed[target]=transformed.get(target,F(0))+co*F(comb(h,i),comb(degree,i))
        coeff=transformed
    assert all(co>=0 for co in coeff.values())
    total=prod(d+1 for d in degrees)
    positive=sum(co>0 for co in coeff.values())
    assert positive>0
    return degrees,total,positive,total-positive,coeff


H=(1+r)*(1+r*X-(1-r)*Z)/(r*(1+X))
Np=X*r**3+Z*q*r*r+Z*q*r+Z*r**3+Z*r+2*Z+q*r*r+q*r-r-2
T=(1+r)**2*Np/(r**3*(1+X))
H1=s.factor(derivative(H))
H2=s.factor(derivative(H1))
H3=s.factor(derivative(H2))
T1=s.factor(derivative(T))
nums=[
 s.cancel(r**3*(1+X)**3*(r-1)*H2),
 s.cancel(r**4*(1+X)**4*(r-1)**2*H3),
 s.cancel(r**4*(1+X)**2*(r-1)*T1),
]
names=['J_second','minus_J_third','T_first']

alpha=1-q
Xsub=r*(1+2*alpha+2*alpha**2)+alpha*(1-3*r-2*alpha*r)*w
Z6=s.Integer(1)
aj=q
for j in range(1,7):
    if j>1:
        aj=aj*(j-1-q)/j
    Z6-=aj*r**j

for name,N,nu,e in zip(names,nums,[3,4,4],[1,1,3]):
    assert s.degree(N,Z)==1
    # Strictly positive coefficient of Z on the larger box.
    K=s.cancel(s.diff(N,Z)/(1+r)**nu)
    stats=bernstein(K.subs({q:u,r:s.Rational(2,17)*v,X:w}),(u,v,w))
    assert stats[2]==stats[1]  # Every Bernstein coefficient is positive.
    print(name,'Z coefficient:',stats[:4],flush=True)

    # The three exact numerator upper bounds.
    P=s.cancel((-N.subs(Z,Z6)/r**e).subs(X,Xsub)/(1-q))
    P=s.Poly(P,q,r,w).as_expr()  # Exact polynomial division is required.
    for label,ql,qu,R in [
        ('middle',s.Rational(2,5),s.Rational(3,5),s.Rational(1,12)),
        ('high',s.Rational(3,5),s.Integer(1),s.Rational(2,17)),
    ]:
        mapped=P.subs({q:ql+(qu-ql)*u,r:R*v})
        stats=bernstein(mapped,(u,v,w))
        print(name,label,':',stats[:4],flush=True)
    low=s.cancel(P.subs(r,q*q*v/2)/q**2)
    low=s.Poly(low,q,v,w).as_expr()
    stats=bernstein(low.subs(q,s.Rational(2,5)*u),(u,v,w))
    print(name,'low:',stats[:4],flush=True)

print('PASS: J_second>0, J_third<0, T_first>0 on all three peak domains.')

