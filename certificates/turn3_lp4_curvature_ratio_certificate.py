"""Exact curvature-ratio certificate on 1/20<=q<=1.

Every sign decision uses integer directed intervals; all symbolic expressions
are reconstructed from the defining H and right-anchor formulas.
No cached objects, sampled assertions, or floating-point signs are used.
"""
import sympy as s
from fractions import Fraction as F
from math import factorial,comb,prod
import time
Q=F
S=10**60
class I:
 def __init__(self,x=0,hi=None,raw=False):
  if raw:self.lo,self.hi=x,hi;return
  if isinstance(x,I):self.lo,self.hi=x.lo,x.hi;return
  x=Q(x);self.lo=(x.numerator*S)//x.denominator;self.hi=-((-x.numerator*S)//x.denominator)
 @classmethod
 def bounds(cls,a,b):return cls(I(a).lo,I(b).hi,True)
 def __add__(self,o):
  o=I(o);return I(self.lo+o.lo,self.hi+o.hi,True)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo,True)
 def __sub__(self,o):return self+-I(o)
 def __rsub__(self,o):return I(o)+-self
 def __mul__(self,o):
  o=I(o);v=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi];return I(min(v)//S,-((-max(v))//S),True)
 __rmul__=__mul__
 def __truediv__(self,o):
  o=I(o);assert not(o.lo<=0<=o.hi)
  a=Q(S*S,o.lo);b=Q(S*S,o.hi);lo=min(a,b);hi=max(a,b)
  rec=I(lo.numerator//lo.denominator,-((-hi.numerator)//hi.denominator),True)
  return self*rec
 def __rtruediv__(self,o):return I(o)/self
 def __pow__(self,k):
  assert isinstance(k,int) and k>=0
  out=I(1);base=self
  while k:
   if k&1:out=out*base
   k//=2
   if k:base=base*base
  return out
 def __abs__(self):
  if self.lo>=0:return self
  if self.hi<=0:return -self
  return I(0,max(-self.lo,self.hi),True)


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

r,q,X,Z=s.symbols('r q X Z');u,v,w=s.symbols('u v w')
def derivative(f):return (1+r)**2*(s.diff(f,r)+q*X/r*s.diff(f,X)-q*Z/(1-r)*s.diff(f,Z))
H=(1+r)*(1+r*X-(1-r)*Z)/(r*(1+X))
Np=X*r**3+Z*q*r*r+Z*q*r+Z*r**3+Z*r+2*Z+q*r*r+q*r-r-2
T=(1+r)**2*Np/(r**3*(1+X))
H2=s.factor(derivative(s.factor(derivative(H))))
N=s.cancel((s.Rational(9,4)*T-r/(1+r)*H2)*r**3*(1+X)**3*(1-r))
K=s.cancel(s.diff(N,Z)/(1+r)**2)
bs=bernstein(K.subs({q:u,r:s.Rational(2,17)*v,X:w}),(u,v,w))
assert bs[2]==bs[1]
print('PASS: positive Z coefficient',bs[:4],flush=True)
# Rigorous binomial lower bound, since its positive coefficients decrease.
Zlo=s.Integer(1);aj=q
for j in range(1,10):
 if j>1:aj=aj*(j-1-q)/j
 Zlo-=aj*r**j
aj=aj*(10-1-q)/10
Zlo-=aj*r**10/(1-r)
V=s.symbols('V')
R=s.cancel(14515200*(1-r)*N.subs({Z:Zlo,X:r*(1+(1-q)*V)})/((1-q)*r**3*(1+r)**2))
R=s.Poly(R,r,q,V)
assert R.domain==s.ZZ
assert R.degree_list()==(14,13,3)
assert len(R.terms())==412
co={mo:int(cc) for mo,cc in R.terms()}
rows=[]
for ir in range(R.degree(r)+1):
 qr=[]
 for iq in range(R.degree(q)+1):
  vs=[co.get((ir,iq,iv),0) for iv in range(R.degree(V)+1)]
  while vs and vs[-1]==0:vs.pop()
  qr.append(vs)
 while qr and not qr[-1]:qr.pop()
 rows.append(qr)
while rows and not rows[-1]:rows.pop()
def hh(arr,x):
 out=I(0)
 for cc in reversed(arr):out=out*x+cc
 return out
def eval_R(rb,qb,vb):return hh([hh([hh(vs,vb) for vs in qr],qb) for qr in rows],rb)

def logunit_fast(xx):
 xx=F(xx);assert 1<=xx<=2
 z=I((xx-1)/(xx+1));zz=z*z
 ans=I(0)
 for j in range(79,-1,-1):ans=ans*zz+F(2,2*j+1)
 ans=ans*z
 tail=2*z**161/(161*(1-zz))
 return I(ans.lo,ans.hi+tail.hi,True)
LOG2=logunit_fast(F(2))
def logrec_point(rr):
 xx=1/F(rr);k=0
 while xx>2:xx/=2;k+=1
 while xx<1:xx*=2;k-=1
 return k*LOG2+logunit_fast(xx)
def logrec_ball(rb):
 assert rb.lo>0
 lo=logrec_point(F(rb.hi,S));hi=logrec_point(F(rb.lo,S))
 return I(lo.lo,hi.hi,True)
def Vball(qb,rb):
 LL=logrec_ball(rb);xx=(1-qb)*LL
 assert xx.lo>=0
 ans=I(0)
 for j in range(49,-1,-1):ans=ans*xx+F(1,factorial(j+1))
 ee=(xx.hi+S-1)//S
 tail=xx**50/factorial(51)*3**ee
 return LL*(ans+I(0,tail.hi,True))
def box(qlo,qhi,vlo,vhi,kind):
 qb=I.bounds(qlo,qhi);vb=I.bounds(vlo,vhi)
 lb=qb*qb/(8+9*qb+24*qb*qb)
 ub=qb*qb/(4+qb*qb) if kind=='low' else (I(F(1,25)) if kind=='mid' else I(F(1,18)))
 rb=lb+(ub-lb)*vb
 if rb.lo<=0 or rb.hi>=S:return None
 return eval_R(rb,qb,Vball(qb,rb))
stack=[]
for qlo,qhi,kind in [(F(1,20),F(2,5),'low'),(F(2,5),F(3,5),'mid'),(F(3,5),F(1),'high')]:
 stack.append((qlo,qhi,F(0),F(1),kind,0,0))
count=good=0;minimum=None;start=time.time()
while stack:
 qlo,qhi,vlo,vhi,kind,dq,dv=stack.pop();count+=1
 vv=box(qlo,qhi,vlo,vhi,kind)
 if vv is not None and vv.lo>0:
  good+=1
  minimum=vv.lo if minimum is None else min(minimum,vv.lo)
 else:
  assert dq+dv<32,(qlo,qhi,vlo,vhi,kind)
  if dq<=dv+2:
   mid=(qlo+qhi)/2
   stack.extend([(qlo,mid,vlo,vhi,kind,dq+1,dv),(mid,qhi,vlo,vhi,kind,dq+1,dv)])
  else:
   mid=(vlo+vhi)/2
   stack.extend([(qlo,qhi,vlo,mid,kind,dq,dv+1),(qlo,qhi,mid,vhi,kind,dq,dv+1)])
print('PASS: full compact curvature ratio; boxes',good,'nodes',count,flush=True)
assert minimum>200*S
print('Strict normalized lower bound:',F(minimum,S),flush=True)
print('Elapsed seconds (diagnostic only):',time.time()-start,flush=True)

