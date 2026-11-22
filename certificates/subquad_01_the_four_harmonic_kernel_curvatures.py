import sympy as S,math
s,r,l=S.symbols('s r l');V=(1+r)*(1+r-2*r*s)/s**2
D=lambda z:S.diff(z,s)+V*S.diff(z,r)
Den=2*(r+l)*(1+l*r)
A=s*s*r*(1+l*l+2*l*r)/Den
B=s*s*r*((1+l*l)*r+2*l)/Den
t=s*r;F1=A+B;F2=A*s+B*t
exprs={'F2':D(D(F2)),'at_s_left':D(D(F2))-s*D(D(F1)), 'at_s_right':D(D(B*t))-s*D(D(B)), 'at_t':2*D(B)*D(t)+B*D(D(t))}
R,U=S.symbols('R U')
expected={'F2':206,'at_s_left':206,'at_s_right':185,'at_t':104}
for name,expr in exprs.items():
 e=S.factor(expr);N,Q=S.fraction(e)
 expected_Q=s*(l+r)**(2 if name=='at_t' else 3)*(l*r+1)**(2 if name=='at_t' else 3)
 assert S.expand(Q-expected_Q)==0
 N,Nden=S.fraction(S.together(N.subs(s,U/r)))
 assert Nden==(r if name in ['F2','at_s_left'] else 1), (name,Nden)
 P=S.Poly(S.expand(N.subs(r,1+R)),R,l,U)
 print('BEGIN',name,'den',Q,'degree',P.degree_list(),'terms',len(P.terms()),flush=True)
 neg=[];total=0
 for kr in range(P.degree(R)+1):
  pp=S.Poly(P.as_expr().coeff(R,kr),l,U)
  if pp.is_zero:continue
  dl=pp.degree(l);du=pp.degree(U)
  for il in range(dl+1):
   for iu in range(du+1):
    b=sum(co*S.Rational(math.comb(il,kl),math.comb(dl,kl))*S.Rational(math.comb(iu,ku),math.comb(du,ku)) for (kl,ku),co in pp.terms() if kl<=il and ku<=iu)
    total+=1
    if b<0:neg.append((kr,il,iu,b))
 print('RESULT',name,'total',total,'negative',len(neg),neg[:8],flush=True)
 assert total==expected[name], (name,total)
 assert not neg, (name,neg)
 if not neg:
  open('finite_zero_all_m_'+name+'_expression.txt','w').write(str(e))

