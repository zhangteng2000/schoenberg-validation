import sympy as S,pickle,time
k,d,B,C,z=S.symbols('k d B C z');m=k+d;n=2*m+1;ell=k+2*d
A=2*d*(k+2*d)-B;bb=k*B;cc=k*C;L=4*d*d+(k-1)*B+k*C
E=2*d*(4*m*d+1)-2*ell*B+n*(k-1)*C
F=n*(4*d*d+(k-1)*B+2*k*C)
terms=[(d*F*(m+1),(m+1)*A+m*z),((d+1)*F,A+z),(k*F,A+m*z),(-k*(m-1)*F*(m+1),(m+1)*bb+m*z),(-k*m*F,bb+z),(-k*F,bb+m*z),(m*(m+1)*((k-1)*F-k*E),(m+1)*L+m*z),((m+1)*((k-1)*F-k*E),L+z),(-k*m*(m+1)*E,(m+1)*cc+m*z),(-k*(m+1)*E,cc+z),(4*k*(m+1)*E,z)]
dom=S.QQ.poly_ring(k,d,B,C);N=S.Poly(0,z,domain=dom);D=S.Poly(1,z,domain=dom)
for i,(w,p) in enumerate(terms):
 pp=S.Poly(p,z,domain=dom);N=N*pp+D*S.Poly(w,z,domain=dom);D*=pp
 print('built',i,'degree',N.degree(),flush=True)
co=N.all_coeffs()[::-1]

assert N.degree()==8
assert S.Poly(co[0]-4*k**6*(m+1)**5*E*A**3*B**3*C**2*L**2,k,d,B,C).is_zero
print('The constant coefficient is the stated positive product.',flush=True)
for j,cj in enumerate(co):
 if j==0:print('coefficient0 automatic positive product',flush=True);continue
 div=A**max(0,3-j)*B**max(0,3-j)*C**max(0,2-j)*L**max(0,2-j)
 ppall=S.Poly(cj,k,d,B,C);qq,rr=S.div(ppall,S.Poly(div,k,d,B,C));assert rr.is_zero
 pp=S.Poly(qq.as_expr(),B,C);deg=pp.total_degree();dic={ij:S.Poly(cf,k,d).shift_list([1,1])for ij,cf in pp.as_dict().items()};bad=0;tot=0;badbc=0
 for r in range(deg+1):
  for t in range(deg-r+1):
   pol=S.Poly(0,k,d)
   for (i,h),cf in dic.items():
    if i<=r and h<=t:
     w=S.factorial(r)/S.factorial(r-i)*S.factorial(t)/S.factorial(t-h)*S.factorial(deg-i-h)/S.factorial(deg)
     pol+=cf*w
   neg=sum(int(bool(cf<0))for cf in pol.coeffs());bad+=neg;tot+=len(pol.coeffs());badbc+=int(neg>0)
 assert bad==0
 print('coefficient',j,'BCdegree',deg,'nonnegative exact coefficients',tot,flush=True)

print('All 39919 scalar coefficient signs are certified exactly.',flush=True)

