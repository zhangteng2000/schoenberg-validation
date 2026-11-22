"""Exact fresh verifier for logarithmic Euler curvature compensation.

All expressions are reconstructed from the stated rational formulas.
There are no cache or pickle inputs. SymPy rational arithmetic only.
This certificate proves a support restriction, not the full conjecture.
"""
import sympy as S
x,m,z,M,v=S.symbols('x m z M v', nonnegative=True)
a=m+1;b=m;n=2*m+1;al=1/a;be=1/b
h=b/(x+al)**2+a/(x+be)**2
G=lambda t:1/t-1/(x+t)
e=a*G(al)+b*G(be)-n*G(S.Integer(1))-x*h
L=S.cancel(S.diff(e,x,2)-S.diff(h,x,2)/S.diff(h,x)*S.diff(e,x))
K=6*a*b*(be-al)**2/((x+al)*(x+be)*(b*(x+be)**3+a*(x+al)**3))

def region_certificate(q, region, m_base=None, m_value=None):
    if region=='small':
        MM=S.Integer(1);nu=z;de=a*(1+z)
    elif region=='middle':
        MM=x/al;nu=b+a*z;de=a*b*(1+z)
    else:
        MM=x*x/(al*be);nu=1+b*z;de=b
    expr=S.factor(S.cancel(L-q*(n-1)*K*MM))
    numer,denom=S.fraction(expr)
    # The cancelled denominator has positive coefficients in m,x.
    dd=S.Poly(S.expand(denom),x,m)
    assert all(c>0 for _,c in dd.terms())
    poly=S.Poly(numer,x);d=poly.degree()
    transformed=S.expand(sum(poly.nth(j)*nu**j*de**(d-j) for j in range(d+1)))
    transformed=S.expand(transformed.subs(m,m_base+M if m_value is None else m_value))
    final=S.Poly(transformed,z,M)
    assert all(c>=0 for _,c in final.terms())
    assert final.coeff_monomial(1)>0
    print('scalar',region,'m_base',m_base,'m_value',m_value,
          'terms',len(final.terms()),'constant',final.coeff_monomial(1),flush=True)

q=1+n*n/(8*a*a*b*b)
for reg in ['small','middle','large']:
    region_certificate(q,reg,m_base=3)
# m=2: the original reflected estimate suffices through x=beta.
for reg in ['small','middle']:
    region_certificate(q,reg,m_value=S.Integer(2))
# In the large region use the stronger Cauchy-weighted kernel estimate.
region_certificate(S.Rational(25,24),'large',m_value=S.Integer(2))

# Verify the stronger negative-kernel estimate at m=2, after alpha=1.
k=S.Rational(3,2);qq=S.Rational(25,24)
A=(k-1)**2*(k+1)
P1=4*k*x**4+(8*k*k+8*k)*x**3+(12*k**3+12*k)*x*x+(4*k**4+4*k**3+4*k*k+4*k)*x+k**5+2*k**3+k
P2=(2*k+2)*x**4+16*k*x**3+(12*k*k+12*k)*x*x+(8*k**3+8*k)*x+k**4+k**3+k*k+k
P=-A*x**4+v*P1+v*v*P2
N=S.expand(qq*(k-1)**2*x*x*(1+v)*(k+v)*(x+v)**2-k*(k-1)**2*(x+v)**4+2*v*P)
poly=S.Poly(S.expand(N.subs(x,k+z)),z)
for j in range(poly.degree()+1):
    cf=S.Poly(poly.nth(j),v)
    c0,c1,c2=[cf.nth(i)for i in range(3)]
    assert c0>0 and c2>0 and 4*c0*c2-c1*c1>0
    assert all(cf.nth(i)>=0 for i in range(3,cf.degree()+1))
    print('n5 kernel z-degree',j,'quadratic discriminant margin',
          S.factor(4*c0*c2-c1*c1),flush=True)
print('ALL LOG CURVATURE COMPENSATION CERTIFICATES PASSED',flush=True)
"""Fresh exact verification of both log reflected-kernel identities."""
import sympy as S
x,v,k,t=S.symbols('x v k t',positive=True)
G=lambda xx,tt:1/tt-1/(xx+tt)
Dplus=-1/((x+t)**2*(v+t))
Dminus=(-G(v,t)-G(x,t)+(v+x)/(x+t)**2)/(v+x)**2
single_difference=-2*v*(x*x-t*v)/(t*(t+v)*(x+t)**2*(x+v)**2)
assert S.cancel(Dminus-Dplus-single_difference)==0
B=(x+k)**3+k*(x+1)**3
Gamma_num=3*((x+k)**4+k*(x+1)**4)
Gamma_den=(x+1)*(x+k)*B
Knum=6*k*(k-1)**2
Kden=(x+1)*(x+k)*B
P1=4*k*x**4+(8*k*k+8*k)*x**3+(12*k**3+12*k)*x*x+(4*k**4+4*k**3+4*k*k+4*k)*x+k**5+2*k**3+k
P2=(2*k+2)*x**4+16*k*x**3+(12*k*k+12*k)*x*x+(8*k**3+8*k)*x+k**4+k**3+k*k+k
P=-(k-1)**2*(k+1)*x**4+v*P1+v*v*P2

def verify_operator(N,D,RN,RD):
    pp=lambda z:S.Poly(z,x,v,k)
    N,D,gn,gd,RN,RD=map(pp,[N,D,Gamma_num,Gamma_den,RN,RD])
    A=N.diff(x)*D-N*D.diff(x)
    LN=(A.diff(x)*D-2*A*D.diff(x))*gd+A*D*gn
    LD=D**3*gd
    assert (LN*RD-RN*LD).is_zero

Dp_num=-((k+v)*(x+k)**2+k*(1+v)*(x+1)**2)
Dp_den=(1+v)*(k+v)*(x+1)**2*(x+k)**2
verify_operator(Dp_num,Dp_den,-Knum,Kden*(1+v)*(k+v))
Dm_num=-2*v*((x*x-v)*(k+v)*(x+k)**2+(x*x-k*v)*(1+v)*(x+1)**2)
Dm_den=(x+v)**2*(1+v)*(k+v)*(x+1)**2*(x+k)**2
verify_operator(Dm_num,Dm_den,12*v*P,(k+v)*(k+x)*(v+1)*(x+1)*B*(x+v)**4)
print('BOTH REFLECTED LOG KERNEL IDENTITIES PASSED',flush=True)

