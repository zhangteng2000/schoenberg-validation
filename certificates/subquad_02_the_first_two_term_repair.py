"""Fresh exact certificate for the first-two-term repair of the zero-family series."""
import sympy as S, math
from pathlib import Path
r,s,q,X=S.symbols('r s q X')
rp=(1+r)*(1+r-2*r*s)/s**2
D=lambda f:S.diff(f,s)+rp*S.diff(f,r)+q*X*rp/r*S.diff(f,X)
F1=s**2*(X+r)/(1+X)
F2=s**3*(X+r**2)/(1+X)
expr=S.factor(D(D(F1+(1-q)*F2/3)))
num,den=S.fraction(expr)
assert S.expand(den-3*r**2*s**2*(X+1)**3)==0
R,e,u,t=S.symbols('R e u t')
ld=1+e*R;xn=1+R+e*(1-e)*R**2*t
# exact polynomial substitutions in two separate steps
p0=S.Poly(S.expand(num),r,s,q,X)
N=0
for (ir,iss,iq,ix),co in p0.terms():
 N+=co*(1+R)**(ir+3-iss)*u**iss*(1-e)**iq*xn**ix*ld**(3-ix)
N=S.Poly(S.expand(N),R,e,u,t)
print('degrees',N.degree_list(),'terms',len(N.terms()),flush=True)
Nexpr=S.factor(N.as_expr())
print('factors',str(Nexpr)[:160],flush=True)
# q=1 is identically zero, remove e.
assert S.expand(N.as_expr().subs(e,0))==0
P=S.Poly(S.cancel(N.as_expr()/e),R,e,u,t)
neg=[];total=0;mint=None;constant_block=[]
for kr in range(P.degree(R)+1):
 pp=S.Poly(P.as_expr().coeff(R,kr),e,u,t)
 if pp.is_zero:continue
 ds=pp.degree_list(); terms=pp.terms()
 for ie in range(ds[0]+1):
  for iu in range(ds[1]+1):
   for it in range(ds[2]+1):
    b=0
    for (ke,ku,kt),co in terms:
     if ke<=ie and ku<=iu and kt<=it:
      b+=co*S.Rational(math.comb(ie,ke),math.comb(ds[0],ke))*S.Rational(math.comb(iu,ku),math.comb(ds[1],ku))*S.Rational(math.comb(it,kt),math.comb(ds[2],kt))
    total+=1
    if kr==0:constant_block.append(b)
    if mint is None or b<mint:mint=b
    if b<0:neg.append((kr,ie,iu,it,b))
 print('Rdegree',kr,'deg',ds,'negative_so_far',len(neg),flush=True)
print('TOTAL',total,'NEGATIVE',len(neg),'MIN',mint,'samples',neg[:30],flush=True)
assert total==1072
assert not neg
assert min(constant_block)==16
print('STRICT_CONSTANT_BLOCK',constant_block,flush=True)
print('ALL EXACT ASSERTIONS PASSED',flush=True)
Path('finite_zero_first_two_certificate_output.txt').write_text(str({'total':total,'negative':len(neg),'min':str(mint),'samples':[str(v) for v in neg[:100]]}))

