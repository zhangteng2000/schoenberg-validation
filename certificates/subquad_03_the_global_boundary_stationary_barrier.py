"""Exact standalone certificate for the boundary stationary barrier.
No numeric optimization or floating-point arithmetic is used.
"""
import sympy as S
from fractions import Fraction
from math import comb
from itertools import product
r,s,q,Y,Z=S.symbols('r s q Y Z',positive=True)
U,V,Q=S.symbols('U V Q')
D=1+s-r
I=(s-r)/D*(1+Y+q*(1+s)*Y/s)/(1+Y)+(1+s)*(r*Y+D-(1+s-2*r)*Z)/(r*D**2*(1+Y))
Ir=S.factor(S.diff(I,r)-q*Z/(1-r)*S.diff(I,Z))
Is=S.factor(S.diff(I,s)+q*Y/s*S.diff(I,Y))
Z0=S.factor(-Ir.subs(Z,0)/S.diff(Ir,Z))
A=(1-r+q*r)*(1+s-r)*(1+s-2*r)+2*r*r*(1-r)
M=S.factor(-s*(1+Y)*A/r*(I-q-1).subs(Z,Z0))
N=S.factor(s*s*(1+Y)**2*(1+s)*A/r*Is.subs(Z,Z0))
assert S.Poly(M,Y).degree()==1
assert S.Poly(N,Y).degree()==2
assert S.denom(M)==1 and S.denom(N)==1
M1=S.factor(S.diff(M,Y));M0=S.factor(M.subs(Y,0));N2=S.factor(S.diff(N,Y,2)/2)
Y0=-M0/M1
H=S.factor(M1*S.diff(N,Y).subs(Y,Y0)/(q*(1+s)))
assert S.denom(H)==1
B=2*q*r-q*s-q+s-1
C=2*q*r*r-3*q*r*s-q*r+q*s*s+q*s-4*r*r+3*r*s+5*r-s*s-s-2
N0_expected=q*q*s*(q-2)*(q-1)*(1+s)**2*B*C*A/M1**2
assert S.cancel(N.subs(Y,Y0)-N0_expected)==0
assert S.factor(C.subs(r,s))==-2*(s-1)**2
assert S.factor(S.diff(C,r).subs(r,s))==(q-5)*(s-1)
# Exact identities used to obtain q > 2s/(1+s) from M>=0 and Y<=1.
F=S.Poly(M1+M0,q)
assert S.factor(F.coeff_monomial(q*q))==(2*r-s-1)*(r*s-r-s*s-s)
F0_expected=-4*s*(1-r)
assert S.factor((M1+M0).subs(q,0)-F0_expected)==0
B0=4*r*r+3*r*s*s-3*r-s**3-2*s*s-s
assert S.cancel((M1+M0).subs(q,2*s/(1+s))+2*s*(s-1)*B0/(s+1)**2)==0
assert S.factor(B0.subs(r,0))==-s*(s+1)**2
assert S.factor(B0.subs(r,s))==2*s*(s-1)*(s+2)
# M0 increases in r on 0<=r<=s<1; its endpoint quotient factors negatively.
assert S.factor(S.diff(M0,r).subs(r,s)/s)==q*q*s-5*q*s-3*q*q+5*q+2
T=S.cancel(M0.subs(r,s)/s)
assert S.factor(T-(s-1)*(2+q-q*q-2*q*s))==0


def bern(poly,vars):
 P=S.Poly(poly,*vars);degs=P.degree_list();co={k:Fraction(v) for k,v in P.terms()}
 for ax,d in enumerate(degs):
  new={};other=[i for i in range(len(vars)) if i!=ax]
  for inds in product(*[range(degs[i]+1) for i in other]):
   base=dict(zip(other,inds))
   for i in range(d+1):
    val=Fraction(0)
    for j in range(i+1):
     key=tuple(j if z==ax else base[z] for z in range(len(vars)))
     val+=co.get(key,0)*Fraction(comb(i,j),comb(d,j))
    if val:
     key=tuple(i if z==ax else base[z] for z in range(len(vars)));new[key]=val
  co=new
 return degs,co

polys={'M1':M1,'N2':N2,'H':H}
expected={(0,'M1'):23,(0,'N2'):47,(0,'H'):225}
for branch in [0]:
 qq=Q
 smax=q/(2-q) if branch==0 else (2-q)/(4+q)
 for name,P in polys.items():
  expr=S.factor(P.subs(r,s*U).subs(s,smax*V).subs(q,qq));num,den=S.fraction(expr)
  assert not den.has(U,V)
  if den.subs(Q,S.Rational(1,2))<0:num,den=-num,-den
  dd,db=bern(den,(Q,))
  assert all(db.get((j,),0)>0 for j in range(dd[0]+1))
  pp=S.Poly(num,Q,U,V)
  mins=tuple(min(m[j] for m in pp.monoms()) for j in range(3))
  num=S.cancel(num/(Q**mins[0]*U**mins[1]*V**mins[2]))
  assert S.denom(num)==1
  ds,bb=bern(num,(Q,U,V))
  assert len(bb)==expected[(branch,name)]
  assert all(v>=0 for v in bb.values())
  # The physical Q,U,V all belong to (0,1). Every Bernstein basis
  # function is then positive, so one positive coefficient ensures strictness.
  assert any(v>0 for v in bb.values())
  print('PASS',branch,name,'denominator',S.factor(den),'stripped',mins,'degrees',ds,'nonzero Bernstein coefficients',len(bb),flush=True)
print('PASS: every symbolic identity, denominator, coefficient sign, and strictness check.',flush=True)

