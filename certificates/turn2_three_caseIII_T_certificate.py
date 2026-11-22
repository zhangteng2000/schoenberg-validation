"""Exact certificate for the unbalanced three-value rational trace inequality.

No floating point arithmetic or saved intermediate data is used.
Run with Python 3 and SymPy.  k>=1, d>=1, r in [0,1].
A nonnegative power expansion in k-1,d-1 and Bernstein expansion in r
certifies each spectral coefficient of the two rational square-root bounds.
"""
import math
import sympy as S

k,d,r,s=S.symbols('k d r s')
m=k+d; n=2*m+1; ell=k+2*d; a=m+1; b=m
A=n*ell+k*r; B=2*m*k*r; C=k*(n-r)
assert S.expand(-k*A+B+ell*C)==0
M=S.cancel((A*C+ell*A*B-k*B*C)/n).expand()
V2=S.expand((B+C-A)**2+2*M)
Disc=S.expand((B+C-A)**2+4*M)
L0=S.expand(B+A-C); H0=S.expand(2*m*n-B); E=ell+k*r
delta=S.expand(Disc-L0**2)
assert S.expand(delta-2*(H0-L0)*E)==0
assert S.expand(delta-4*k*(1-r)*n*E)==0

# The nine persistent-root contributions after canceling multiplicities.
rows=[(-d,a,A),(-d-1,b,A),(-k,1,A),
      (d,a,C),(d-1,b,C),(-ell,1,C),
      (-(m-1),a,B),(-m,b,B),(-1,1,B)]
domain=S.ZZ.poly_ring(k,d,r)
Dr=S.Poly(1,s,domain=domain); Nr=S.Poly(0,s,domain=domain)
for weight,scale,value in rows:
    atom=S.expand(scale*value)
    den=S.Poly(s+atom,s,domain=domain)
    Nr=Nr*den+Dr*S.Poly(S.expand(weight*atom**2),s,domain=domain)
    Dr=Dr*den

def endpoint_coefficients(L):
    Dc=(s**4+(a+b)*L*s**3+((a*a+b*b)*M+a*b*Disc)*s*s
        +a*b*(a+b)*M*L*s+a*a*b*b*M*M)
    Nc=((m*a*a+(m+1)*b*b)*V2*s**3
        +L*(m*a*a*(b*V2+a*M)+(m+1)*b*b*(a*V2+b*M))*s*s
        +(m*a*a*(b*b*M*V2+a*b*M*Disc)
          +(m+1)*b*b*(a*a*M*V2+a*b*M*Disc))*s
        +a*a*b*b*M*M*(m*a+(m+1)*b)*L)
    N=Nr*S.Poly(Dc,s,domain=domain)+Dr*S.Poly(Nc,s,domain=domain)
    assert N.degree()==11
    return [S.Poly(N.nth(j),k,d,r) for j in range(12)]

Nlo=endpoint_coefficients(L0); Nhi=endpoint_coefficients(H0)

def exact_check(poly):
    """Return count of positive coefficients; assert all others are zero.

    First shift k=K+1,d=D+1 by exact binomial transforms. For degree R
    in r, the coefficient of r^t(1-r)^(R-t) is
       sum_{i<=t} c_i * binom(R-i,t-i).
    This is the usual Bernstein coefficient multiplied by binom(R,t).
    The positive factor does not change its sign.
    """
    K,D,R=poly.degree_list()
    arrays=[[[0]*(D+1) for _ in range(K+1)] for _ in range(R+1)]
    for (i,j,t),coef in poly.as_dict().items(): arrays[t][i][j]=int(coef)
    for t in range(R+1):
        tmp=[[0]*(D+1) for _ in range(K+1)]
        out=[[0]*(D+1) for _ in range(K+1)]
        for i in range(K+1):
            for j in range(D+1):
                c=arrays[t][i][j]
                if c:
                    for jj in range(j+1): tmp[i][jj]+=c*math.comb(j,jj)
        for i in range(K+1):
            for j in range(D+1):
                c=tmp[i][j]
                if c:
                    for ii in range(i+1): out[ii][j]+=c*math.comb(i,ii)
        arrays[t]=out
    positive=0
    for t in range(R+1):
        for i in range(K+1):
            for j in range(D+1):
                c=sum(arrays[h][i][j]*math.comb(R-h,t-h) for h in range(t+1))
                assert c>=0, (i,j,t,c)
                positive+=(c>0)
    return (K,D,R),positive

families=[('upper',L0-E,E),
          ('lower',4*L0**2+delta-4*L0*E,4*L0*E)]
total=0
for name,f0,f1 in families:
    P0=S.Poly(f0,k,d,r); P1=S.Poly(f1,k,d,r)
    counts=[]
    for j in range(12):
        poly=P0*Nlo[j]+P1*Nhi[j]
        degrees,count=exact_check(poly)
        counts.append(count)
        if j==11:
            base=S.Poly(poly.as_expr().subs({k:1,d:1}),r)
            R=base.degree()
            assert all(sum(base.nth(i)*math.comb(R-i,t-i) for i in range(t+1))>0 for t in range(R+1))
        print(name,j,'degrees',degrees,'positive',count,flush=True)
    total+=sum(counts)
    print(name,'counts',counts,'sum',sum(counts),flush=True)
assert total==140107
print('VERIFIED: all 140107 nonzero exact coefficients are positive.',flush=True)

