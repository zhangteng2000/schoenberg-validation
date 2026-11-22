"""Exact Bernstein certificate for theta J'''+J''/2<0 on peak domains."""
import sympy as s
from fractions import Fraction as F
from math import comb, prod

r,q,X,Z=s.symbols('r q X Z')
u,v,w=s.symbols('u v w')


def derivative(f):
    return (1+r)**2*(s.diff(f,r)+q*X/r*s.diff(f,X)-q*Z/(1-r)*s.diff(f,Z))


def bernstein(poly):
    pp=s.Poly(poly,u,v,w)
    degrees=pp.degree_list()
    coeff={mo:F(int(co.p),int(co.q)) for mo,co in pp.terms()}
    for axis,degree in enumerate(degrees):
        out={}
        for mo,co in coeff.items():
            i=mo[axis]
            for h in range(i,degree+1):
                target=list(mo)
                target[axis]=h
                target=tuple(target)
                out[target]=out.get(target,F(0))+co*F(comb(h,i),comb(degree,i))
        coeff=out
    assert all(z>=0 for z in coeff.values())
    total=prod(d+1 for d in degrees)
    positive=sum(z>0 for z in coeff.values())
    assert positive>0
    return degrees,total,positive,total-positive


H=(1+r)*(1+r*X-(1-r)*Z)/(r*(1+X))
H1=s.factor(derivative(H))
H2=s.factor(derivative(H1))
H3=s.factor(derivative(H2))
weighted=s.factor(r/(1+r)*H3+s.Rational(1,2)*H2)
N=s.cancel(2*r**3*(1+X)**4*(r-1)**2*weighted)
assert s.degree(N,Z)==1
K=s.cancel(s.diff(N,Z)/(1+r)**3)
stats=bernstein(K.subs({q:u,r:s.Rational(2,17)*v,X:w}))
assert stats[1]==stats[2]
print('Z coefficient:',stats,flush=True)

alpha=1-q
Z6=s.Integer(1)
aj=q
for j in range(1,7):
    if j>1:
        aj=aj*(j-1-q)/j
    Z6-=aj*r**j
Xold=r*(1+2*alpha+2*alpha**2)+alpha*(1-3*r-2*alpha*r)*w
P=s.cancel((-N.subs(Z,Z6)/r).subs(X,Xold)/(1-q))
P=s.Poly(P,q,r,w).as_expr()
for label,ql,qu,R in [
    ('middle',s.Rational(2,5),s.Rational(3,5),s.Rational(1,12)),
    ('upper_middle',s.Rational(3,5),s.Rational(9,10),s.Rational(2,17)),
]:
    stats=bernstein(P.subs({q:ql+(qu-ql)*u,r:R*v}))
    print(label,stats,flush=True)
low=s.Poly(s.cancel(P.subs(r,q*q*v/2)/q**2),q,v,w).as_expr()
print('low',bernstein(low.subs(q,s.Rational(2,5)*u)),flush=True)

R=F(2,17)
assert R<F(9,10)**10
assert R**9<F(251,1700)**10
assert F(3,10)-2*R-F(1,5)*R>0
Xnew=r*(1+2*alpha+2*alpha**2)+alpha*(s.Rational(3,10)-2*r-2*alpha*r)*w
Pnew=s.cancel((-N.subs(Z,Z6)/r).subs(X,Xnew)/(1-q))
Pnew=s.Poly(Pnew,q,r,w).as_expr()
print('high',bernstein(Pnew.subs({q:s.Rational(9,10)+u/10,r:s.Rational(2,17)*v})),flush=True)
assert 12**2*2<17**2
print('PASS: theta J_third+J_second/2<0 throughout the full peak domain.')

