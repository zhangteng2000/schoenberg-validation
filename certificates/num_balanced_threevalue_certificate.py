"""Exact certificate for the BALANCED THREE-VALUE SUBFAMILY only.

For n=2m+1, roots are (-m-t)^m, (2m*t)^1, (m-t)^m,
with |t|<=m/n. Set v=1-(n*t/m)^2 in [0,1].

After cancelling common positive factors and scaling the determinant variable,
the cross-determinant difference is
    A^(m-1) D - B^(m-1) E.
All polynomials below have nonnegative coefficients in Y for m>=2, v in[0,1],
except that H=A-B and J need verification.

For m>=3, H>=0 and
    J=B^2(D-E)+(m-1)BHD+binom(m-1,2)H^2D >=0
coefficientwise in Y imply the target by the ordinary binomial expansion:
    A^(m-1)D-B^(m-1)E >= B^(m-3) J.
The certificates for H/v and J/v^2 are obtained by writing m=r+3,
expanding in nonnegative powers of r,Y, and representing every remaining
polynomial in v in its exact Bernstein basis on [0,1]. All Bernstein
coefficients are nonnegative rational numbers. m=2 is checked directly.

This is NOT a certificate for arbitrary roots or arbitrary three-value
multiplicities. SymPy is used only for exact integer/rational algebra.
"""
import sympy as sp
m,Y,v,r=sp.symbols('m Y v r')
n=2*m+1;a=m+1;u2=1-v
S=(2*m*m-2*m+1)*u2+n
T=(2*m-1)*u2+n
Q=n*n+u2
R=n*n-u2
P=lambda e:sp.Poly(e,m,Y,v,domain=sp.QQ)
LA=P(m**4*Y**2+2*a*a*m*m*S*Y+a**4*T*T)
LB=P(Y*Y+2*S*Y+T*T)
R1=P(m**4*Y**2+2*m*m*Q*Y+R*R)
RB=P(Y*Y+2*Q*Y+R*R)
Fa=P(Y+4*a*a*u2)
Fb=P(Y+4*m*m*u2)
F1=P(Y+4*u2)
A=LA*LB;B=R1*Fa*Fb
D=LA*LB*LB;E=RB*R1*Fb*F1
H=A-B
J=B*B*(D-E)+P(m-1)*B*H*D+P(sp.Rational(1,2)*(m-1)*(m-2))*H*H*D

def verify_bernstein(Z,vanishing_order,shift=None):
    if shift is None:
        ZZ=sp.Poly(Z.as_expr().subs(m,2),Y,v,domain=sp.QQ)
        groups={}
        for (iy,iv),c in ZZ.terms():
            assert iv>=vanishing_order
            groups.setdefault((iy,),{})[iv-vanishing_order]=c
    else:
        ZZ=sp.Poly(Z.as_expr().subs(m,r+shift),Y,v,r,domain=sp.QQ)
        groups={}
        for (iy,iv,ir),c in ZZ.terms():
            assert iv>=vanishing_order
            groups.setdefault((iy,ir),{})[iv-vanishing_order]=c
    count=0
    for key,co in groups.items():
        deg=max(co)
        for k in range(deg+1):
            b=sum(co.get(j,0)*sp.binomial(k,j)/sp.binomial(deg,j)
                  for j in range(k+1))
            assert b>=0,(key,k,b)
            count+=1
    return len(groups),count

assert H.eval(v,0).is_zero
assert J.eval(v,0).is_zero
print('H/v, m>=3:',verify_bernstein(H,1,shift=3))
print('J/v^2, m>=3:',verify_bernstein(J,2,shift=3))
Z2=A*D-B*E
print('Target/v^2, m=2:',verify_bernstein(Z2,2,shift=None))
# The quadratic binomial truncation is exact at m=3.
assert (A*A*D-B*B*E-J).eval(m,3).is_zero
print('All exact nonnegativity assertions passed.')

