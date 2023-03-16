import random
var('x y')
bound = 10
f = sum([random.randint(-9,9) * x^i * y^(bound-1-i) for i in range(bound)])
R.<x> = QQ[]


ls = [[(i,f(i,j)) for i in range(bound)] for j in range(bound)]
l = [R.lagrange_polynomial(ls[j]) for j in range(bound)]
m = [[(i,l[i].monomial_coefficient(x^(bound-1-j))) for i in range(bound)] for j in range(bound)]
ms = [R.lagrange_polynomial(m[j]) for j in range(bound)]

var('x y')
print(f(x,y).expand())
sum([ms[i].subs(x=y) * x^(bound-i-1) for i in range(bound)]).expand()

