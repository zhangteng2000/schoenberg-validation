"""Fresh exact scalar certificate for minority averaging, central critical>0.

This verifies both one-scale bounds used under the logarithmic minimum's
scale constraint. It uses only polynomial division and positive coefficients;
no cached inputs, floating-point computations, or optimization are involved.
"""
import sympy as S
k,h,p,z,K,H=S.symbols('k h p z K H',nonnegative=True)
n=k+h+1;c=(n-2)/n;d=(k-1)/k;d3=(n-3)/n
F=h*z*z+((k+1)*p+k-h)*z-p-k/n
Mnum=3*c*(p+z)-d*(1-z)*(5*p+4*z+1)+d3*(2*p+1)

def positive_coefficients(expr,kbase,variables,label):
    shifted=S.cancel(expr.subs(h,k+2+H).subs(k,kbase+K))
    numerator,denominator=S.fraction(shifted)
    pp=S.Poly(S.expand(numerator),*variables)
    dd=S.Poly(S.expand(denominator),*variables)
    assert all(coef>=0 for _,coef in pp.terms()), label
    assert all(coef>0 for _,coef in dd.terms()), (label,'denominator')
    print(label,'positive numerator terms',len(pp.terms()),flush=True)
    return pp

qm,rm=S.div(Mnum,F,z)
assert S.cancel(qm-4*(k-1)/(h*k))==0
positive_coefficients(qm,1,(p,z,K,H),'M quotient')
P=positive_coefficients(rm,1,(p,z,K,H),'M remainder')
assert P.coeff_monomial(1)>0

for name,kappa in [('first',S.Integer(1)),('second',(n-1)/(n+1))]:
    T=p/kappa;M=Mnum/(4*(p+z));q0=d*(1-z)
    phi=q0*T*(T+p)+M*(T+p+1)*(T+p)-kappa*T*(T+p+z)
    endpoint=S.cancel(4*(p+z)*kappa*kappa*phi/p)
    quotient,remainder=S.div(endpoint,F,z)
    assert S.Poly(remainder,z).degree()<=1
    positive_coefficients(quotient,2,(p,K,H),name+' endpoint quotient')
    positive_coefficients(S.diff(remainder,z),2,(p,K,H),name+' remainder slope')
    P=positive_coefficients(remainder.subs(z,1/(k+1)),2,(p,K,H),name+' remainder at secular lower bound')
    assert P.coeff_monomial(1)>0
    # For k=1 the unscaled Phi decreases as z increases; z<=1 suffices.
    at_one=S.factor(endpoint.subs({k:1,z:1}))
    shifted=S.cancel(at_one.subs(h,3+H))
    nn,dd=S.fraction(shifted)
    pp=S.Poly(S.expand(nn),p,H);dp=S.Poly(S.expand(dd),p,H)
    assert all(coef>0 for _,coef in pp.terms())
    assert all(coef>0 for _,coef in dp.terms())
    print(name,'k=1 endpoint',at_one,flush=True)
print('BOTH CENTRAL-POSITIVE ONE-SCALE MINORITY AVERAGING BOUNDS VERIFIED',flush=True)

