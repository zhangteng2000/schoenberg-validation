from flint import fmpq_mpoly_ctx
import json,time
ctx=fmpq_mpoly_ctx.get(('C','B','T','Q'))
class Rat:
 def __init__(self,n,d=1):
  if isinstance(n,Rat):self.n,self.d=n.n,n.d;return
  n=ctx.constant(n) if isinstance(n,(int,str)) else n;d=ctx.constant(d) if isinstance(d,(int,str)) else d
  g=n.gcd(d);self.n=n//g;self.d=d//g
  lc=self.d.leading_coefficient();self.n=self.n/lc;self.d=self.d/lc
 def __add__(self,o):
  o=Rat(o);g=self.d.gcd(o.d);ad=o.d//g;bd=self.d//g;return Rat(self.n*ad+o.n*bd,self.d*ad)
 __radd__=__add__
 def __neg__(self):return Rat(-self.n,self.d)
 def __sub__(self,o):return self+-Rat(o)
 def __rsub__(self,o):return Rat(o)+-self
 def __mul__(self,o):
  o=Rat(o);g=self.n.gcd(o.d);h=o.n.gcd(self.d);return Rat((self.n//g)*(o.n//h),(self.d//h)*(o.d//g))
 __rmul__=__mul__
 def __truediv__(self,o):o=Rat(o);return self*Rat(o.d,o.n)
 def __rtruediv__(self,o):return Rat(o)/self
 def __pow__(self,i):return Rat(self.n**i,self.d**i) if i>=0 else Rat(self.d**(-i),self.n**(-i))
C,Bb,Tt,q=[Rat(g) for g in ctx.gens()]
c=(1+C)/2;t=(1+q+(2-q)*Tt)/3;b=t+(1-t)*Bb;x=1-(1-c)*(b-t)/(b+c)
u=x-t;v=b+t;A=1-t;B=1+t;E=1+x;L=1-b;d=x+b;W=1-x+b-t;e=1-x
K=2*x*x+x-b;J=q*x*E*d/u-K;J0=K+q*d*(b-x*x)/v;h1=-L*L*J/d**3;h0=1+q*L/v+L*L*J0/d**3
wa=(1+b)*(b-t)/(W*d*E);wb=e*(x+t)/(W*d*L);w1=(c-x*b)/(E*L);P0=b*wb;P1=x*wa
T0=q*b*wb/v-e*(1+b)/W*(x-b+q*d*b/v)/d**3
T1=-q*x*wa/u-e*(1+b)/W*(x-b-q*d*x/u)/d**3
M=A-q*e
# Verify the spectral moment identities and the two Hermite coefficients.
def zero(r): return not r.n
assert zero(wa+wb+w1-1)
assert zero(-x*wa+b*wb+w1)
assert zero(x*x*wa+b*b*wb+w1-c)
assert zero(-x**3*wa+b**3*wb+w1+(1-c)*t)
assert zero(c-(x*b-t)/W)
da=-1-q*x/u;sec=-x/d;dc2=(sec-da)/d;dc3=(da-2*sec)/d**2
assert zero(h1-(x+da*E+dc2*E**2+dc3*E**2*L))
da=0;db=1+q*b/v;sec=b/d;dc2=sec/d;dc3=(db-2*sec)/d**2
assert zero(h0-(dc2*E**2+dc3*E**2*L))
assert zero(T0-q*P0/B-e/(W*d**3*v*B)*((1+b)*B*v*(b-x)-q*b*d*e*(d+B)))
print('Derived rational right-Radau components.',flush=True)
UB=1+q*L/v;LB=B/(B-q*L)
def transform(poly):
 for j in range(4):
  deg=poly.degrees()[j];groups=[{} for _ in range(deg+1)]
  for ex,co in poly.terms():
   ee=list(ex);k=ee[j];ee[j]=0;groups[k][tuple(ee)]=co
  groups=[ctx.from_dict(g) for g in groups];res=groups[deg]
  for k in range(deg-1,-1,-1):res=ctx.gen(j)*res+groups[k]*(1+ctx.gen(j))**(deg-k)
  poly=res
 return poly
for label,H in [('upper',UB),('lower',LB)]:
 X=(H-h0)/h1
 core=A*X*(T0+T1*X)+q*(P0+P1*X)*X+H*(M*(T0+T1*X)-q*M*(P0+P1*X)/B+q*w1*X*(1+A/B))
 signs=[]
 for side,p in [('numerator',core.n),('denominator',core.d)]:
  cf,factors=p.factor();sgn=1 if cf>0 else -1
  for fac,ex in factors:
   if ex%2==0:
    if side=='numerator':
     assert fac==ctx.gen(3) or fac==-ctx.gen(3)
     # The only even numerator factor is a power of Q=q>0.
    continue
   tp=transform(fac);co=tp.coeffs()
   sign=1 if all(z>=0 for z in co) else -1 if all(z<=0 for z in co) else 0
   assert sign!=0,(label,side,'mixed Bernstein signs',fac.degrees())
   assert any(z!=0 for z in co)
   # Also retain strictness at c=1/2, corresponding to C=0.
   assert any(z!=0 and ee[0]==0 for ee,z in tp.terms())
   sgn*=sign
   if len(fac)>1000:print(label,side,'large factor degree',fac.degrees(),'strict-sign coefficients',len(tp),'sign',sign,flush=True)
  signs.append(sgn)
 assert signs[0]==signs[1],(label,signs)
 print('PASS',label,'endpoint core > 0; numerator and denominator signs',signs,flush=True)
print('PASS: exact right-Radau joint-stationary barrier certificate.',flush=True)

