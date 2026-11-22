"""Optional exact verification of the six-coefficient identity in the proof."""
import sympy as S
m,y=S.symbols('m y')
n=2*m+1;alpha=m/(m+1)
H=(y+alpha/n)/(y+alpha)**3+(y-1/n)/(y+1)**3-1/(y+m)**2
b=[4*m**5*(m+2),2*m**3*(6*m**3+19*m**2+8*m-3),2*m**2*(m+1)*(6*m**3+24*m**2+7*m-5),m*(m+1)**2*(4*m**3+30*m**2+8*m-5),(m+1)**2*(8*m**3+12*m**2+4*m-1),(m+1)**3*(2*m+1)]
P=sum(co*y**j for j,co in enumerate(b))
assert S.cancel(H-y*P/((m+y)**2*n*(y+1)**3*(m*y+m+y)**3))==0
print('The elementary six-coefficient rational identity is exact.')

