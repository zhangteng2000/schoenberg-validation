import sympy as S
r,d,s=S.symbols('r d s');m=r+d;h=r+2*d+1;g=2*d+1
terms=[(-d,(m+1)*h),(-d-1,m*h),(-r,h),(d+1,(m+1)*r),(d,m*r),(-h,r),(m,(m+1)*g),(m+1,m*g)]
dom=S.ZZ.poly_ring(r,d);N=S.Poly(0,s,domain=dom);D=S.Poly(1,s,domain=dom)
for w,v in terms:
 pp=S.Poly(s+v,s,domain=dom);N=N*pp+D*S.Poly(w*v*v,s,domain=dom);D*=pp
assert N.degree()==6
for i,cf in enumerate(N.all_coeffs()[::-1]):
 factor=S.factor(cf);pp=S.Poly(cf,r,d).shift_list([1,0]);neg=sum(int(bool(c<0)) for c in pp.coeffs())
 print(i,'negative',neg,'terms',len(pp.terms()),'factor',factor,flush=True)
 assert neg==0
print('ALL COEFFICIENTS NONNEGATIVE EXACTLY',flush=True)

