"""Pure integer/rational certificate. No floating-point operation enters a sign test."""
from fractions import Fraction as Q
from math import factorial,comb
from functools import lru_cache
S=10**70
ORDER=40
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

def logunit(x):
 x=Q(x);assert 1<=x<=2;z=(x-1)/(x+1);M=55
 p=sum((2*z**(2*j+1)/(2*j+1) for j in range(M)),Q(0))
 tail=2*z**(2*M+1)/((2*M+1)*(1-z*z))
 return I.bounds(p,p+tail)
LOG2=logunit(Q(2))
@lru_cache(None)
def logrec(r):
 x=1/Q(r);m=0
 while x>2:x/=2;m+=1
 while x<1:x*=2;m-=1
 return m*LOG2+logunit(x)

def Nterms(r):
 r=Q(r);return [(2-r,r,I(0)),(r**3*(1-2*r),-r**3,logrec(r)),(-(1-r)*(2+r+r*r+2*r**3),(1-r)*(r+r*r),logrec(1-r))]

def Hterms(r,s):
 def first(x):return [(1+x,I(0)),((1+x)*x*x,logrec(x)),(-(1+x)*(1-x)**2,logrec(1-x))]
 def second(x):return [(x,I(0)),(x*x,logrec(x))]
 out=[]
 for u,l in first(r):
  for v,m in second(s):out.append((u*v,Q(0),l+m))
 for u,l in first(s):
  for v,m in second(r):out.append((-u*v,Q(0),l+m))
 return out

def coeffs(terms,b):
 assert sum(a for a,_,_ in terms)==0
 cs=[]
 for j in range(1,ORDER+1):
  cs.append(sum((I(a)*l**j/factorial(j)+I(bb)*l**(j-1)/factorial(j-1) for a,bb,l in terms),I(0)))
 B=I(b);tail=I(0)
 for a,bb,l in terms:
  x=l*B;assert x.lo>=0 and x.hi<8*S
  # exp(x)<exp(8)<3**8, since e<3.
  tail+=(abs(I(a))*x**(ORDER+1)/factorial(ORDER+1)/B+abs(I(bb))*x**ORDER/factorial(ORDER))*3**8
 return cs,tail

def bernstein(cs,a,b):
 deg=len(cs)-1;a=I(a);width=I(b)-a
 ps=[sum((cs[i]*comb(i,j)*a**(i-j)*width**j for i in range(j,deg+1)),I(0)) for j in range(deg+1)]
 return [sum((ps[j]*Q(comb(k,j),comb(deg,j)) for j in range(k+1)),I(0)) for k in range(deg+1)]

def certify(name,terms,alo,ahi):
 cs,tail=coeffs(terms,ahi);stack=[(alo,ahi,0)];cells=0;lower=None
 while stack:
  a,b,depth=stack.pop();bs=bernstein(cs,a,b);vals=[x-tail for x in bs]
  if all(x.lo>0 for x in vals):
   cells+=1;v=min(x.lo for x in vals);lower=v if lower is None else min(lower,v)
  else:
   assert depth<12,(name,a,b,min(x.lo for x in vals))
   mid=(a+b)/2;stack.extend([(a,mid,depth+1),(mid,b,depth+1)])
 # Decimal formatting is diagnostic only; every decision above used integers.
 print('PASS',name,'cells',cells,'strict lower >',Q(lower,S),flush=True)

bands=[(Q(3,10),Q(2,5),Q(7,200)),(Q(2,5),Q(1,2),Q(9,200)),(Q(1,2),Q(3,5),Q(11,200)),(Q(3,5),Q(7,10),Q(3,50)),(Q(7,10),Q(4,5),Q(67,1000))]
for lo,hi,r in bands:
 certify('N '+str((lo,hi,r)),Nterms(r),1-hi,1-lo)
 certify('H-half '+str((lo,hi,r)),Hterms(r,r/(2+r)),1-hi,1-lo)
r=Q(5,69)
certify('N high',Nterms(r),Q(0),Q(1,5))
certify('H integer29',Hterms(r,Q(1,28)),Q(0),Q(1,5))
certify('H k>=2',Hterms(r,2*r/(3+r)),Q(0),Q(1,5))
assert Q(200,9)**7>Q(26,3)**10
assert -2+Q(5,6)+Q(128,225)==-Q(269,450)
print('ALL PURE RATIONAL BRANCH CERTIFICATES PASS',flush=True)

