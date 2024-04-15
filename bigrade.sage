import itertools



# My Macualay2 code
# R = QQ[x_0..x_4,y_0,y_1, Degrees=>{5:{0,1},{1,-2},{1,-3}}]
# f = random({0,2}, R)
# g = random({0,3},R)
# F = y_0 * f + y_1 * g
# J = R/ideal jacobian F
# for p from 0 to 5 list hilbertFunction({p,0}, J)
# basis({0,0}, J)
# basis({1,0}, J)
# basis({2,0}, J)



# Nick's Macaulay2 code
# R = QQ[x_0..x_4,y_0,y_1, Degrees=>{5:{0,1},{1,-2},{1,-3}}]
# f = random({1,0}, R)
# J = R/ideal jacobian f
# for p from 0 to 5 list hilbertFunction({p,0}, J)
# basis({0,0}, J)
# basis({1,0}, J)
# basis({2,0}, J)




# for convenience
# x = gens[0]
# y = gens[1]
# z = gens[2]
# a = gens[3]
# b = gens[4]
# c = gens[5]

# f = a^3+b^3-c^2-x*y*z
# f = x^2 + y^2-z^2
# g = x^3+y^3+z^3

# f_0 = f
# f_1 = g


# poly_list = list(var('f_%d' % i) for i in range(num_poly))


# Z = V(f0,..,fc)
# then canonical bundle of Z w_z \cong O_Z(m) with
# m = sum(di) - n - 1


# J_p,m = H^{n-p,p}_van(X)




# d_i = deg(f_i)

# bigrading is
# deg(x_i) = (0,1)
# deg(y_i) = (1,-deg(f_i))


# number of variables
n = 5

# number of hypersurfaces in complete intersection

num_poly = 2

# x_vars = list('x_%d' % i for i in range(n))
# y_vars = list('y_%d' % i for i in range(num_poly))
p = 13

R.<x_0,x_1,x_2,x_3,x_4,y_1,y_2> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,y_1,y_2> = QQ[]
# R = PolynomialRing(QQ, names= x_vars + y_vars) # M + F call this S
gens = R.gens()

x_vars = gens[:n]
y_vars = gens[n:]

# J_p,m consists of p copies of y_i and sum_d_i copies of x_j
# maybe compute by first all monomials in J_p,m for fixed p, then finding a basis

f = (5/4)*x_0^2+(3/7)*x_0*x_1+(2/3)*x_1^2+4*x_0*x_2+x_1*x_2+(2/3)*x_2^2+(1/5)*x_0*x_3+(1/7)*x_1*x_3+x_2*x_3+5*x_3^2+(3/2)*x_0*x_4+(3/4)*x_1*x_4+x_2*x_4+(1/2)*x_3*x_4+(1/2)*x_4^2
g = (7/8)*x_0^3+(3/2)*x_0^2*x_1+x_0*x_1^2+(9/5)*x_1^3+(5/9)*x_0^2*x_2+(4/9)*x_0*x_1*x_2+(1/10)*x_1^2*x_2+(7/9)*x_0*x_2^2+x_1*x_2^2+9*x_2^3+(10/9)*x_0^2*x_3+(3/8)*x_0*x_1*x_3+(4/3)*x_1^2*x_3+(8/3)*x_0*x_2*x_3+2*x_1*x_2*x_3+(7/8)*x_2^2*x_3+(2/5)*x_0*x_3^2+(10/7)*x_1*x_3^2+(4/9)*x_2*x_3^2+(9/2)*x_3^3+(1/2)*x_0^2*x_4+(3/2)*x_0*x_1*x_4+(2/5)*x_1^2*x_4+(1/3)*x_0*x_2*x_4+x_1*x_2*x_4+(7/9)*x_2^2*x_4+(3/2)*x_0*x_3*x_4+4*x_1*x_3*x_4+(3/10)*x_2*x_3*x_4+(3/4)*x_3^2*x_4+(8/3)*x_0*x_4^2+(3/7)*x_1*x_4^2+(2/5)*x_2*x_4^2+(5/2)*x_3*x_4^2+7*x_4^3

B = [1,x_0*x_1*x_4*y_1, x_0*x_2*x_4*y_1, x_0*x_3^2*y_1, x_0*x_3*x_4*y_1, x_0*x_4^2*y_1, x_1^2*x_4*y_1, x_1*x_2*x_4*y_1, x_1*x_3^2*y_1, x_1*x_3*x_4*y_1, x_1*x_4^2*y_1, x_2^2*x_3*y_1, x_2^2*x_4*y_1, x_2*x_3^2*y_1, x_2*x_3*x_4*y_1, x_2*x_4^2*y_1, x_3^3*y_1, x_3^2*x_4*y_1, x_3*x_4^2*y_1, x_4^3*y_1,x_4^6*y_1^2]

f_0 = f; f_1 = g;
poly_list = [f_0,f_1]
d = [_.degree() for _ in poly_list]
m = sum([f.degree() for f in poly_list]) - n #- 1

F = sum([y_vars[i] * poly_list[i] for i in range(num_poly)])
I  = F.jacobian_ideal()
J = R.quotient_ring(I)



def monomial_degree(m):
#   e looks like  [(1, 1, 1, 1, 1)]
    e = m.exponents()[0]
#   x_i contribute (0, sum(e[:n]))
#   y_i contribute (sum(e[n:],sum([-poly_list[i].degree() for i in range(num_poly)] ))
    return [sum(e[n:]), sum(e[:n])+sum([-d[i] for i in range(num_poly)] )]

# each monomial is of the form m / F^j; returns j
def pole_order(m):
    e = m.exponents()[0]
    return sum(e[n:])

def sigma(g):
    return sum([g.coefficient(monomial)* monomial^p for monomial in g.monomials()])

def frobenius(g,prec=2):
    d = pole_order(g)
    summer = 0
    sigma_g = sigma(g)
    fj = R(1)
    for j in range(prec):
#         cacheing may be appropriate
        numer = binomial(-d,j) * binomial(d+prec-1,d+j) * sigma_g*sigma(fj)
        summer += numer
        fj *= f
    return summer


def frobenius_on_cohom(i,prec = 2):
    g = B[i]*prod(gens)
    g = frobenius(g,prec)
    return R(g *p^(n -2) / prod(gens))