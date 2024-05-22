import itertools


# R = QQ[x_0,x_1,x_2,y_0, Degrees=>{3:{0,1},{1,-3}}]
# f = x_0^3 + x_1^3 + x_2^3 - x_0*x_1*x_2
# f = random({0,2}, R)
# F = y_0 * f
# J = R/ideal jacobian F
# for p from 0 to 1 list hilbertFunction({p,0}, J)
# basis({0,0}, J)
# basis({1,0}, J)


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

# My Macualay2 code
# R = QQ[x_0..x_3,y_0,y_1, Degrees=>{4:{0,1},{1,-2},{1,-2}}]
# f = random({0,2}, R)
# g = random({0,2},R)
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
n = 3

# number of hypersurfaces in complete intersection

num_poly = 1

# x_vars = list('x_%d' % i for i in range(n))
# y_vars = list('y_%d' % i for i in range(num_poly))
p = 17
prec = 2

R.<x_0,x_1,x_2,y_0> = QQ[]

# R.<x_0,x_1,x_2,x_3,y_0,y_1> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,y_1,y_2> = QQ[]
# R = PolynomialRing(QQ, names= x_vars + y_vars) # M + F call this S
gens = R.gens()

x_vars = gens[:n]
y_vars = gens[n:]

# J_p,m consists of p copies of y_i and sum_d_i copies of x_j
# maybe compute by first all monomials in J_p,m for fixed p, then finding a basis
# f = (7/4)*x_0^2+(6/7)*x_0*x_1+x_1^2+(1/5)*x_0*x_2+(3/10)*x_1*x_2+(3/8)*x_2^2+x_0*x_3+(3/2)*x_1*x_3+(3/2)*x_2*x_3+(2/5)*x_3^2
# g = (3/5)*x_0^2+(3/7)*x_0*x_1+(9/7)*x_1^2+(1/7)*x_0*x_2+(7/6)*x_1*x_2+4*x_2^2+2*x_0*x_3+(7/6)*x_1*x_3+(1/4)*x_2*x_3+(3/4)*x_3^2

# f = x_0^3 + x_1^3 + x_2^3 - x_0*x_1*x_2
# f = (5/4)*x_0^2+(3/7)*x_0*x_1+(2/3)*x_1^2+4*x_0*x_2+x_1*x_2+(2/3)*x_2^2+(1/5)*x_0*x_3+(1/7)*x_1*x_3+x_2*x_3+5*x_3^2+(3/2)*x_0*x_4+(3/4)*x_1*x_4+x_2*x_4+(1/2)*x_3*x_4+(1/2)*x_4^2
# g = (7/8)*x_0^3+(3/2)*x_0^2*x_1+x_0*x_1^2+(9/5)*x_1^3+(5/9)*x_0^2*x_2+(4/9)*x_0*x_1*x_2+(1/10)*x_1^2*x_2+(7/9)*x_0*x_2^2+x_1*x_2^2+9*x_2^3+(10/9)*x_0^2*x_3+(3/8)*x_0*x_1*x_3+(4/3)*x_1^2*x_3+(8/3)*x_0*x_2*x_3+2*x_1*x_2*x_3+(7/8)*x_2^2*x_3+(2/5)*x_0*x_3^2+(10/7)*x_1*x_3^2+(4/9)*x_2*x_3^2+(9/2)*x_3^3+(1/2)*x_0^2*x_4+(3/2)*x_0*x_1*x_4+(2/5)*x_1^2*x_4+(1/3)*x_0*x_2*x_4+x_1*x_2*x_4+(7/9)*x_2^2*x_4+(3/2)*x_0*x_3*x_4+4*x_1*x_3*x_4+(3/10)*x_2*x_3*x_4+(3/4)*x_3^2*x_4+(8/3)*x_0*x_4^2+(3/7)*x_1*x_4^2+(2/5)*x_2*x_4^2+(5/2)*x_3*x_4^2+7*x_4^3

f = x_0^3 + x_1^3 + x_2^3 - x_0*x_1*x_2

B = [R(1), x_2^3*y_0]
# B = [1,x_0*x_1*x_4*y_1, x_0*x_2*x_4*y_1, x_0*x_3^2*y_1, x_0*x_3*x_4*y_1, x_0*x_4^2*y_1, x_1^2*x_4*y_1, x_1*x_2*x_4*y_1, x_1*x_3^2*y_1, x_1*x_3*x_4*y_1, x_1*x_4^2*y_1, x_2^2*x_3*y_1, x_2^2*x_4*y_1, x_2*x_3^2*y_1, x_2*x_3*x_4*y_1, x_2*x_4^2*y_1, x_3^3*y_1, x_3^2*x_4*y_1, x_3*x_4^2*y_1, x_4^3*y_1,x_4^6*y_1^2]
# B = [R(1), x_3^2*y_1]

# f_0 = f; f_1 = g;
# poly_list = [f_0,f_1]
poly_list = [f]
d = [_.degree() for _ in poly_list]
m = sum([f.degree() for f in poly_list]) - n #- 1

F = sum([y_vars[i] * poly_list[i] for i in range(num_poly)])
I  = F.jacobian_ideal()
J = R.quotient_ring(I)

xI = R.ideal([_*F.derivative(_) for _ in R.gens()])# + [f])
xJ = R.quotient(xI)

# toString basis({i,0},R)
P_0_mons = [R(1)]
P_1_mons = [x_0^3*y_0, x_0^2*x_1*y_0, x_0^2*x_2*y_0, x_0*x_1^2*y_0, x_0*x_1*x_2*y_0, x_0*x_2^2*y_0, x_1^3*y_0,x_1^2*x_2*y_0, x_1*x_2^2*y_0, x_2^3*y_0]
P_2_mons = [x_0^6*y_0^2, x_0^5*x_1*y_0^2, x_0^5*x_2*y_0^2, x_0^4*x_1^2*y_0^2, x_0^4*x_1*x_2*y_0^2, x_0^4*x_2^2*y_0^2, x_0^3*x_1^3*y_0^2, x_0^3*x_1^2*x_2*y_0^2, x_0^3*x_1*x_2^2*y_0^2, x_0^3*x_2^3*y_0^2, x_0^2*x_1^4*y_0^2,x_0^2*x_1^3*x_2*y_0^2, x_0^2*x_1^2*x_2^2*y_0^2, x_0^2*x_1*x_2^3*y_0^2, x_0^2*x_2^4*y_0^2, x_0*x_1^5*y_0^2,x_0*x_1^4*x_2*y_0^2, x_0*x_1^3*x_2^2*y_0^2, x_0*x_1^2*x_2^3*y_0^2, x_0*x_1*x_2^4*y_0^2, x_0*x_2^5*y_0^2,x_1^6*y_0^2, x_1^5*x_2*y_0^2, x_1^4*x_2^2*y_0^2, x_1^3*x_2^3*y_0^2, x_1^2*x_2^4*y_0^2, x_1*x_2^5*y_0^2,x_2^6*y_0^2]
# P_1_mons = [x_0^3*y_1, x_0^2*x_1*y_1, x_0^2*x_2*y_1, x_0^2*x_3*y_1, x_0^2*x_4*y_1, x_0^2*y_0, x_0*x_1^2*y_1,x_0*x_1*x_2*y_1, x_0*x_1*x_3*y_1, x_0*x_1*x_4*y_1, x_0*x_1*y_0, x_0*x_2^2*y_1, x_0*x_2*x_3*y_1, x_0*x_2*x_4*y_1,x_0*x_2*y_0, x_0*x_3^2*y_1, x_0*x_3*x_4*y_1, x_0*x_3*y_0, x_0*x_4^2*y_1, x_0*x_4*y_0, x_1^3*y_1, x_1^2*x_2*y_1,x_1^2*x_3*y_1, x_1^2*x_4*y_1, x_1^2*y_0, x_1*x_2^2*y_1, x_1*x_2*x_3*y_1, x_1*x_2*x_4*y_1, x_1*x_2*y_0,x_1*x_3^2*y_1, x_1*x_3*x_4*y_1, x_1*x_3*y_0, x_1*x_4^2*y_1, x_1*x_4*y_0, x_2^3*y_1, x_2^2*x_3*y_1, x_2^2*x_4*y_1,x_2^2*y_0, x_2*x_3^2*y_1, x_2*x_3*x_4*y_1, x_2*x_3*y_0, x_2*x_4^2*y_1, x_2*x_4*y_0, x_3^3*y_1, x_3^2*x_4*y_1,x_3^2*y_0, x_3*x_4^2*y_1, x_3*x_4*y_0, x_4^3*y_1, x_4^2*y_0]
# P_2_mons = [x_0^6*y_1^2, x_0^5*x_1*y_1^2, x_0^5*x_2*y_1^2, x_0^5*x_3*y_1^2, x_0^5*x_4*y_1^2, x_0^5*y_0*y_1,x_0^4*x_1^2*y_1^2, x_0^4*x_1*x_2*y_1^2, x_0^4*x_1*x_3*y_1^2, x_0^4*x_1*x_4*y_1^2, x_0^4*x_1*y_0*y_1,x_0^4*x_2^2*y_1^2, x_0^4*x_2*x_3*y_1^2, x_0^4*x_2*x_4*y_1^2, x_0^4*x_2*y_0*y_1, x_0^4*x_3^2*y_1^2,x_0^4*x_3*x_4*y_1^2, x_0^4*x_3*y_0*y_1, x_0^4*x_4^2*y_1^2, x_0^4*x_4*y_0*y_1, x_0^4*y_0^2, x_0^3*x_1^3*y_1^2,x_0^3*x_1^2*x_2*y_1^2, x_0^3*x_1^2*x_3*y_1^2, x_0^3*x_1^2*x_4*y_1^2, x_0^3*x_1^2*y_0*y_1, x_0^3*x_1*x_2^2*y_1^2,x_0^3*x_1*x_2*x_3*y_1^2, x_0^3*x_1*x_2*x_4*y_1^2, x_0^3*x_1*x_2*y_0*y_1, x_0^3*x_1*x_3^2*y_1^2,x_0^3*x_1*x_3*x_4*y_1^2, x_0^3*x_1*x_3*y_0*y_1, x_0^3*x_1*x_4^2*y_1^2, x_0^3*x_1*x_4*y_0*y_1, x_0^3*x_1*y_0^2,x_0^3*x_2^3*y_1^2, x_0^3*x_2^2*x_3*y_1^2, x_0^3*x_2^2*x_4*y_1^2, x_0^3*x_2^2*y_0*y_1, x_0^3*x_2*x_3^2*y_1^2,x_0^3*x_2*x_3*x_4*y_1^2, x_0^3*x_2*x_3*y_0*y_1, x_0^3*x_2*x_4^2*y_1^2, x_0^3*x_2*x_4*y_0*y_1, x_0^3*x_2*y_0^2,x_0^3*x_3^3*y_1^2, x_0^3*x_3^2*x_4*y_1^2, x_0^3*x_3^2*y_0*y_1, x_0^3*x_3*x_4^2*y_1^2, x_0^3*x_3*x_4*y_0*y_1,x_0^3*x_3*y_0^2, x_0^3*x_4^3*y_1^2, x_0^3*x_4^2*y_0*y_1, x_0^3*x_4*y_0^2, x_0^2*x_1^4*y_1^2,x_0^2*x_1^3*x_2*y_1^2, x_0^2*x_1^3*x_3*y_1^2, x_0^2*x_1^3*x_4*y_1^2, x_0^2*x_1^3*y_0*y_1, x_0^2*x_1^2*x_2^2*y_1^2,x_0^2*x_1^2*x_2*x_3*y_1^2, x_0^2*x_1^2*x_2*x_4*y_1^2, x_0^2*x_1^2*x_2*y_0*y_1, x_0^2*x_1^2*x_3^2*y_1^2,x_0^2*x_1^2*x_3*x_4*y_1^2, x_0^2*x_1^2*x_3*y_0*y_1, x_0^2*x_1^2*x_4^2*y_1^2, x_0^2*x_1^2*x_4*y_0*y_1,x_0^2*x_1^2*y_0^2, x_0^2*x_1*x_2^3*y_1^2, x_0^2*x_1*x_2^2*x_3*y_1^2, x_0^2*x_1*x_2^2*x_4*y_1^2,x_0^2*x_1*x_2^2*y_0*y_1, x_0^2*x_1*x_2*x_3^2*y_1^2, x_0^2*x_1*x_2*x_3*x_4*y_1^2, x_0^2*x_1*x_2*x_3*y_0*y_1,x_0^2*x_1*x_2*x_4^2*y_1^2, x_0^2*x_1*x_2*x_4*y_0*y_1, x_0^2*x_1*x_2*y_0^2, x_0^2*x_1*x_3^3*y_1^2,x_0^2*x_1*x_3^2*x_4*y_1^2, x_0^2*x_1*x_3^2*y_0*y_1, x_0^2*x_1*x_3*x_4^2*y_1^2, x_0^2*x_1*x_3*x_4*y_0*y_1,x_0^2*x_1*x_3*y_0^2, x_0^2*x_1*x_4^3*y_1^2, x_0^2*x_1*x_4^2*y_0*y_1, x_0^2*x_1*x_4*y_0^2, x_0^2*x_2^4*y_1^2,x_0^2*x_2^3*x_3*y_1^2, x_0^2*x_2^3*x_4*y_1^2, x_0^2*x_2^3*y_0*y_1, x_0^2*x_2^2*x_3^2*y_1^2,x_0^2*x_2^2*x_3*x_4*y_1^2, x_0^2*x_2^2*x_3*y_0*y_1, x_0^2*x_2^2*x_4^2*y_1^2, x_0^2*x_2^2*x_4*y_0*y_1,x_0^2*x_2^2*y_0^2, x_0^2*x_2*x_3^3*y_1^2, x_0^2*x_2*x_3^2*x_4*y_1^2, x_0^2*x_2*x_3^2*y_0*y_1,x_0^2*x_2*x_3*x_4^2*y_1^2, x_0^2*x_2*x_3*x_4*y_0*y_1, x_0^2*x_2*x_3*y_0^2, x_0^2*x_2*x_4^3*y_1^2,x_0^2*x_2*x_4^2*y_0*y_1, x_0^2*x_2*x_4*y_0^2, x_0^2*x_3^4*y_1^2, x_0^2*x_3^3*x_4*y_1^2, x_0^2*x_3^3*y_0*y_1,x_0^2*x_3^2*x_4^2*y_1^2, x_0^2*x_3^2*x_4*y_0*y_1, x_0^2*x_3^2*y_0^2, x_0^2*x_3*x_4^3*y_1^2,x_0^2*x_3*x_4^2*y_0*y_1, x_0^2*x_3*x_4*y_0^2, x_0^2*x_4^4*y_1^2, x_0^2*x_4^3*y_0*y_1, x_0^2*x_4^2*y_0^2,x_0*x_1^5*y_1^2, x_0*x_1^4*x_2*y_1^2, x_0*x_1^4*x_3*y_1^2, x_0*x_1^4*x_4*y_1^2, x_0*x_1^4*y_0*y_1,x_0*x_1^3*x_2^2*y_1^2, x_0*x_1^3*x_2*x_3*y_1^2, x_0*x_1^3*x_2*x_4*y_1^2, x_0*x_1^3*x_2*y_0*y_1,x_0*x_1^3*x_3^2*y_1^2, x_0*x_1^3*x_3*x_4*y_1^2, x_0*x_1^3*x_3*y_0*y_1, x_0*x_1^3*x_4^2*y_1^2,x_0*x_1^3*x_4*y_0*y_1, x_0*x_1^3*y_0^2, x_0*x_1^2*x_2^3*y_1^2, x_0*x_1^2*x_2^2*x_3*y_1^2,x_0*x_1^2*x_2^2*x_4*y_1^2, x_0*x_1^2*x_2^2*y_0*y_1, x_0*x_1^2*x_2*x_3^2*y_1^2, x_0*x_1^2*x_2*x_3*x_4*y_1^2,x_0*x_1^2*x_2*x_3*y_0*y_1, x_0*x_1^2*x_2*x_4^2*y_1^2, x_0*x_1^2*x_2*x_4*y_0*y_1, x_0*x_1^2*x_2*y_0^2,x_0*x_1^2*x_3^3*y_1^2, x_0*x_1^2*x_3^2*x_4*y_1^2, x_0*x_1^2*x_3^2*y_0*y_1, x_0*x_1^2*x_3*x_4^2*y_1^2,x_0*x_1^2*x_3*x_4*y_0*y_1, x_0*x_1^2*x_3*y_0^2, x_0*x_1^2*x_4^3*y_1^2, x_0*x_1^2*x_4^2*y_0*y_1,x_0*x_1^2*x_4*y_0^2, x_0*x_1*x_2^4*y_1^2, x_0*x_1*x_2^3*x_3*y_1^2, x_0*x_1*x_2^3*x_4*y_1^2, x_0*x_1*x_2^3*y_0*y_1,x_0*x_1*x_2^2*x_3^2*y_1^2, x_0*x_1*x_2^2*x_3*x_4*y_1^2, x_0*x_1*x_2^2*x_3*y_0*y_1, x_0*x_1*x_2^2*x_4^2*y_1^2,x_0*x_1*x_2^2*x_4*y_0*y_1, x_0*x_1*x_2^2*y_0^2, x_0*x_1*x_2*x_3^3*y_1^2, x_0*x_1*x_2*x_3^2*x_4*y_1^2,x_0*x_1*x_2*x_3^2*y_0*y_1, x_0*x_1*x_2*x_3*x_4^2*y_1^2, x_0*x_1*x_2*x_3*x_4*y_0*y_1, x_0*x_1*x_2*x_3*y_0^2,x_0*x_1*x_2*x_4^3*y_1^2, x_0*x_1*x_2*x_4^2*y_0*y_1, x_0*x_1*x_2*x_4*y_0^2, x_0*x_1*x_3^4*y_1^2,x_0*x_1*x_3^3*x_4*y_1^2, x_0*x_1*x_3^3*y_0*y_1, x_0*x_1*x_3^2*x_4^2*y_1^2, x_0*x_1*x_3^2*x_4*y_0*y_1,x_0*x_1*x_3^2*y_0^2, x_0*x_1*x_3*x_4^3*y_1^2, x_0*x_1*x_3*x_4^2*y_0*y_1, x_0*x_1*x_3*x_4*y_0^2,x_0*x_1*x_4^4*y_1^2, x_0*x_1*x_4^3*y_0*y_1, x_0*x_1*x_4^2*y_0^2, x_0*x_2^5*y_1^2, x_0*x_2^4*x_3*y_1^2,x_0*x_2^4*x_4*y_1^2, x_0*x_2^4*y_0*y_1, x_0*x_2^3*x_3^2*y_1^2, x_0*x_2^3*x_3*x_4*y_1^2, x_0*x_2^3*x_3*y_0*y_1,x_0*x_2^3*x_4^2*y_1^2, x_0*x_2^3*x_4*y_0*y_1, x_0*x_2^3*y_0^2, x_0*x_2^2*x_3^3*y_1^2, x_0*x_2^2*x_3^2*x_4*y_1^2,x_0*x_2^2*x_3^2*y_0*y_1, x_0*x_2^2*x_3*x_4^2*y_1^2, x_0*x_2^2*x_3*x_4*y_0*y_1, x_0*x_2^2*x_3*y_0^2,x_0*x_2^2*x_4^3*y_1^2, x_0*x_2^2*x_4^2*y_0*y_1, x_0*x_2^2*x_4*y_0^2, x_0*x_2*x_3^4*y_1^2, x_0*x_2*x_3^3*x_4*y_1^2,x_0*x_2*x_3^3*y_0*y_1, x_0*x_2*x_3^2*x_4^2*y_1^2, x_0*x_2*x_3^2*x_4*y_0*y_1, x_0*x_2*x_3^2*y_0^2,x_0*x_2*x_3*x_4^3*y_1^2, x_0*x_2*x_3*x_4^2*y_0*y_1, x_0*x_2*x_3*x_4*y_0^2, x_0*x_2*x_4^4*y_1^2,x_0*x_2*x_4^3*y_0*y_1, x_0*x_2*x_4^2*y_0^2, x_0*x_3^5*y_1^2, x_0*x_3^4*x_4*y_1^2, x_0*x_3^4*y_0*y_1,x_0*x_3^3*x_4^2*y_1^2, x_0*x_3^3*x_4*y_0*y_1, x_0*x_3^3*y_0^2, x_0*x_3^2*x_4^3*y_1^2, x_0*x_3^2*x_4^2*y_0*y_1,x_0*x_3^2*x_4*y_0^2, x_0*x_3*x_4^4*y_1^2, x_0*x_3*x_4^3*y_0*y_1, x_0*x_3*x_4^2*y_0^2, x_0*x_4^5*y_1^2,x_0*x_4^4*y_0*y_1, x_0*x_4^3*y_0^2, x_1^6*y_1^2, x_1^5*x_2*y_1^2, x_1^5*x_3*y_1^2, x_1^5*x_4*y_1^2, x_1^5*y_0*y_1,x_1^4*x_2^2*y_1^2, x_1^4*x_2*x_3*y_1^2, x_1^4*x_2*x_4*y_1^2, x_1^4*x_2*y_0*y_1, x_1^4*x_3^2*y_1^2,x_1^4*x_3*x_4*y_1^2, x_1^4*x_3*y_0*y_1, x_1^4*x_4^2*y_1^2, x_1^4*x_4*y_0*y_1, x_1^4*y_0^2, x_1^3*x_2^3*y_1^2,x_1^3*x_2^2*x_3*y_1^2, x_1^3*x_2^2*x_4*y_1^2, x_1^3*x_2^2*y_0*y_1, x_1^3*x_2*x_3^2*y_1^2, x_1^3*x_2*x_3*x_4*y_1^2,x_1^3*x_2*x_3*y_0*y_1, x_1^3*x_2*x_4^2*y_1^2, x_1^3*x_2*x_4*y_0*y_1, x_1^3*x_2*y_0^2, x_1^3*x_3^3*y_1^2,x_1^3*x_3^2*x_4*y_1^2, x_1^3*x_3^2*y_0*y_1, x_1^3*x_3*x_4^2*y_1^2, x_1^3*x_3*x_4*y_0*y_1, x_1^3*x_3*y_0^2,x_1^3*x_4^3*y_1^2, x_1^3*x_4^2*y_0*y_1, x_1^3*x_4*y_0^2, x_1^2*x_2^4*y_1^2, x_1^2*x_2^3*x_3*y_1^2,x_1^2*x_2^3*x_4*y_1^2, x_1^2*x_2^3*y_0*y_1, x_1^2*x_2^2*x_3^2*y_1^2, x_1^2*x_2^2*x_3*x_4*y_1^2,x_1^2*x_2^2*x_3*y_0*y_1, x_1^2*x_2^2*x_4^2*y_1^2, x_1^2*x_2^2*x_4*y_0*y_1, x_1^2*x_2^2*y_0^2,x_1^2*x_2*x_3^3*y_1^2, x_1^2*x_2*x_3^2*x_4*y_1^2, x_1^2*x_2*x_3^2*y_0*y_1, x_1^2*x_2*x_3*x_4^2*y_1^2,x_1^2*x_2*x_3*x_4*y_0*y_1, x_1^2*x_2*x_3*y_0^2, x_1^2*x_2*x_4^3*y_1^2, x_1^2*x_2*x_4^2*y_0*y_1,x_1^2*x_2*x_4*y_0^2, x_1^2*x_3^4*y_1^2, x_1^2*x_3^3*x_4*y_1^2, x_1^2*x_3^3*y_0*y_1, x_1^2*x_3^2*x_4^2*y_1^2,x_1^2*x_3^2*x_4*y_0*y_1, x_1^2*x_3^2*y_0^2, x_1^2*x_3*x_4^3*y_1^2, x_1^2*x_3*x_4^2*y_0*y_1, x_1^2*x_3*x_4*y_0^2,x_1^2*x_4^4*y_1^2, x_1^2*x_4^3*y_0*y_1, x_1^2*x_4^2*y_0^2, x_1*x_2^5*y_1^2, x_1*x_2^4*x_3*y_1^2,x_1*x_2^4*x_4*y_1^2, x_1*x_2^4*y_0*y_1, x_1*x_2^3*x_3^2*y_1^2, x_1*x_2^3*x_3*x_4*y_1^2, x_1*x_2^3*x_3*y_0*y_1,x_1*x_2^3*x_4^2*y_1^2, x_1*x_2^3*x_4*y_0*y_1, x_1*x_2^3*y_0^2, x_1*x_2^2*x_3^3*y_1^2, x_1*x_2^2*x_3^2*x_4*y_1^2,x_1*x_2^2*x_3^2*y_0*y_1, x_1*x_2^2*x_3*x_4^2*y_1^2, x_1*x_2^2*x_3*x_4*y_0*y_1, x_1*x_2^2*x_3*y_0^2,x_1*x_2^2*x_4^3*y_1^2, x_1*x_2^2*x_4^2*y_0*y_1, x_1*x_2^2*x_4*y_0^2, x_1*x_2*x_3^4*y_1^2, x_1*x_2*x_3^3*x_4*y_1^2,x_1*x_2*x_3^3*y_0*y_1, x_1*x_2*x_3^2*x_4^2*y_1^2, x_1*x_2*x_3^2*x_4*y_0*y_1, x_1*x_2*x_3^2*y_0^2,x_1*x_2*x_3*x_4^3*y_1^2, x_1*x_2*x_3*x_4^2*y_0*y_1, x_1*x_2*x_3*x_4*y_0^2, x_1*x_2*x_4^4*y_1^2,x_1*x_2*x_4^3*y_0*y_1, x_1*x_2*x_4^2*y_0^2, x_1*x_3^5*y_1^2, x_1*x_3^4*x_4*y_1^2, x_1*x_3^4*y_0*y_1,x_1*x_3^3*x_4^2*y_1^2, x_1*x_3^3*x_4*y_0*y_1, x_1*x_3^3*y_0^2, x_1*x_3^2*x_4^3*y_1^2, x_1*x_3^2*x_4^2*y_0*y_1,x_1*x_3^2*x_4*y_0^2, x_1*x_3*x_4^4*y_1^2, x_1*x_3*x_4^3*y_0*y_1, x_1*x_3*x_4^2*y_0^2, x_1*x_4^5*y_1^2,x_1*x_4^4*y_0*y_1, x_1*x_4^3*y_0^2, x_2^6*y_1^2, x_2^5*x_3*y_1^2, x_2^5*x_4*y_1^2, x_2^5*y_0*y_1,x_2^4*x_3^2*y_1^2, x_2^4*x_3*x_4*y_1^2, x_2^4*x_3*y_0*y_1, x_2^4*x_4^2*y_1^2, x_2^4*x_4*y_0*y_1, x_2^4*y_0^2,x_2^3*x_3^3*y_1^2, x_2^3*x_3^2*x_4*y_1^2, x_2^3*x_3^2*y_0*y_1, x_2^3*x_3*x_4^2*y_1^2, x_2^3*x_3*x_4*y_0*y_1,x_2^3*x_3*y_0^2, x_2^3*x_4^3*y_1^2, x_2^3*x_4^2*y_0*y_1, x_2^3*x_4*y_0^2, x_2^2*x_3^4*y_1^2,x_2^2*x_3^3*x_4*y_1^2, x_2^2*x_3^3*y_0*y_1, x_2^2*x_3^2*x_4^2*y_1^2, x_2^2*x_3^2*x_4*y_0*y_1, x_2^2*x_3^2*y_0^2,x_2^2*x_3*x_4^3*y_1^2, x_2^2*x_3*x_4^2*y_0*y_1, x_2^2*x_3*x_4*y_0^2, x_2^2*x_4^4*y_1^2, x_2^2*x_4^3*y_0*y_1,x_2^2*x_4^2*y_0^2, x_2*x_3^5*y_1^2, x_2*x_3^4*x_4*y_1^2, x_2*x_3^4*y_0*y_1, x_2*x_3^3*x_4^2*y_1^2,x_2*x_3^3*x_4*y_0*y_1, x_2*x_3^3*y_0^2, x_2*x_3^2*x_4^3*y_1^2, x_2*x_3^2*x_4^2*y_0*y_1, x_2*x_3^2*x_4*y_0^2,x_2*x_3*x_4^4*y_1^2, x_2*x_3*x_4^3*y_0*y_1, x_2*x_3*x_4^2*y_0^2, x_2*x_4^5*y_1^2, x_2*x_4^4*y_0*y_1,x_2*x_4^3*y_0^2, x_3^6*y_1^2, x_3^5*x_4*y_1^2, x_3^5*y_0*y_1, x_3^4*x_4^2*y_1^2, x_3^4*x_4*y_0*y_1, x_3^4*y_0^2,x_3^3*x_4^3*y_1^2, x_3^3*x_4^2*y_0*y_1, x_3^3*x_4*y_0^2, x_3^2*x_4^4*y_1^2, x_3^2*x_4^3*y_0*y_1,x_3^2*x_4^2*y_0^2, x_3*x_4^5*y_1^2, x_3*x_4^4*y_0*y_1, x_3*x_4^3*y_0^2, x_4^6*y_1^2, x_4^5*y_0*y_1, x_4^4*y_0^2]


# P_1_mons = [x_0^2*y_0, x_0^2*y_1, x_0*x_1*y_0, x_0*x_1*y_1, x_0*x_2*y_0, x_0*x_2*y_1, x_0*x_3*y_0, x_0*x_3*y_1,x_1^2*y_0, x_1^2*y_1, x_1*x_2*y_0, x_1*x_2*y_1, x_1*x_3*y_0, x_1*x_3*y_1, x_2^2*y_0, x_2^2*y_1, x_2*x_3*y_0,x_2*x_3*y_1, x_3^2*y_0, x_3^2*y_1]
# P_2_mons = [x_0^4*y_0^2, x_0^4*y_0*y_1, x_0^4*y_1^2, x_0^3*x_1*y_0^2, x_0^3*x_1*y_0*y_1, x_0^3*x_1*y_1^2,x_0^3*x_2*y_0^2, x_0^3*x_2*y_0*y_1, x_0^3*x_2*y_1^2, x_0^3*x_3*y_0^2, x_0^3*x_3*y_0*y_1, x_0^3*x_3*y_1^2,x_0^2*x_1^2*y_0^2, x_0^2*x_1^2*y_0*y_1, x_0^2*x_1^2*y_1^2, x_0^2*x_1*x_2*y_0^2, x_0^2*x_1*x_2*y_0*y_1,x_0^2*x_1*x_2*y_1^2, x_0^2*x_1*x_3*y_0^2, x_0^2*x_1*x_3*y_0*y_1, x_0^2*x_1*x_3*y_1^2, x_0^2*x_2^2*y_0^2,x_0^2*x_2^2*y_0*y_1, x_0^2*x_2^2*y_1^2, x_0^2*x_2*x_3*y_0^2, x_0^2*x_2*x_3*y_0*y_1, x_0^2*x_2*x_3*y_1^2,x_0^2*x_3^2*y_0^2, x_0^2*x_3^2*y_0*y_1, x_0^2*x_3^2*y_1^2, x_0*x_1^3*y_0^2, x_0*x_1^3*y_0*y_1, x_0*x_1^3*y_1^2,x_0*x_1^2*x_2*y_0^2, x_0*x_1^2*x_2*y_0*y_1, x_0*x_1^2*x_2*y_1^2, x_0*x_1^2*x_3*y_0^2, x_0*x_1^2*x_3*y_0*y_1,x_0*x_1^2*x_3*y_1^2, x_0*x_1*x_2^2*y_0^2, x_0*x_1*x_2^2*y_0*y_1, x_0*x_1*x_2^2*y_1^2, x_0*x_1*x_2*x_3*y_0^2,x_0*x_1*x_2*x_3*y_0*y_1, x_0*x_1*x_2*x_3*y_1^2, x_0*x_1*x_3^2*y_0^2, x_0*x_1*x_3^2*y_0*y_1, x_0*x_1*x_3^2*y_1^2,x_0*x_2^3*y_0^2, x_0*x_2^3*y_0*y_1, x_0*x_2^3*y_1^2, x_0*x_2^2*x_3*y_0^2, x_0*x_2^2*x_3*y_0*y_1,x_0*x_2^2*x_3*y_1^2, x_0*x_2*x_3^2*y_0^2, x_0*x_2*x_3^2*y_0*y_1, x_0*x_2*x_3^2*y_1^2, x_0*x_3^3*y_0^2,x_0*x_3^3*y_0*y_1, x_0*x_3^3*y_1^2, x_1^4*y_0^2, x_1^4*y_0*y_1, x_1^4*y_1^2, x_1^3*x_2*y_0^2, x_1^3*x_2*y_0*y_1,x_1^3*x_2*y_1^2, x_1^3*x_3*y_0^2, x_1^3*x_3*y_0*y_1, x_1^3*x_3*y_1^2, x_1^2*x_2^2*y_0^2, x_1^2*x_2^2*y_0*y_1,x_1^2*x_2^2*y_1^2, x_1^2*x_2*x_3*y_0^2, x_1^2*x_2*x_3*y_0*y_1, x_1^2*x_2*x_3*y_1^2, x_1^2*x_3^2*y_0^2,x_1^2*x_3^2*y_0*y_1, x_1^2*x_3^2*y_1^2, x_1*x_2^3*y_0^2, x_1*x_2^3*y_0*y_1, x_1*x_2^3*y_1^2, x_1*x_2^2*x_3*y_0^2,x_1*x_2^2*x_3*y_0*y_1, x_1*x_2^2*x_3*y_1^2, x_1*x_2*x_3^2*y_0^2, x_1*x_2*x_3^2*y_0*y_1, x_1*x_2*x_3^2*y_1^2,x_1*x_3^3*y_0^2, x_1*x_3^3*y_0*y_1, x_1*x_3^3*y_1^2, x_2^4*y_0^2, x_2^4*y_0*y_1, x_2^4*y_1^2, x_2^3*x_3*y_0^2,x_2^3*x_3*y_0*y_1, x_2^3*x_3*y_1^2, x_2^2*x_3^2*y_0^2, x_2^2*x_3^2*y_0*y_1, x_2^2*x_3^2*y_1^2, x_2*x_3^3*y_0^2,x_2*x_3^3*y_0*y_1, x_2*x_3^3*y_1^2, x_3^4*y_0^2, x_3^4*y_0*y_1, x_3^4*y_1^2]

def monomial_to_vector(m):
    return list(m.exponents()[0])

def vector_to_monomial(v):
    return prod([gens[i]^v[i] for i in range(n+num_poly)])

P_0_pts = [monomial_to_vector(_) for _ in P_0_mons]
P_1_pts = [monomial_to_vector(_) for _ in P_1_mons]
P_2_pts = [monomial_to_vector(_) for _ in P_2_mons]
# P_3_pts = [monomial_to_vector(_) for _ in P_3_mons]


Pn_minus_1 = P_2_mons
Pn_minus_1_pts = P_2_pts
size_pn_minus_1 = len(Pn_minus_1_pts)

def to_pn_minus_1_basis(g):
    return_vec = size_pn_minus_1 * [0]
    for monomial in g.monomials():
        ind = Pn_minus_1.index(monomial)
        return_vec[ind] = g.monomial_coefficient(monomial)
    return return_vec

def from_pn_minus_1_basis(g_vec):
    return sum([g_vec[i] * Pn_minus_1[i] for i in range(size_pn_minus_1)])


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

# def degree_vector(v):
#     return sum([weights[i] * v[i] for i in range(n+num_poly)]) // fdeg


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
        fj *= F
    return summer


def frobenius_on_cohom(i,prec = 2):
    g = B[i]*prod(gens)
    g = frobenius(g,prec)
    return R(g *p^(n -2) / prod(gens))


lift_dict = {}

def lift_poly(g):
    summer = (n+num_poly) * [0]
    for monomial in g.monomials():

        if not monomial in lift_dict.keys():
            c = xJ(monomial).lift()
            r = (monomial-c).lift(xI)
            lift_dict[monomial] = r
        monomial_lift = lift_dict[monomial]
        c = g.monomial_coefficient(monomial)
        for i in range(n+num_poly):
            summer[i] += c*monomial_lift[i]
    return summer   


def to_uvg(h):
    hdict = h.dict()
    return_list = []


    for etuple in hdict.keys():
        vector = list(etuple)
        c = hdict[etuple]
        # todo: save on divisibility?

        g = 0
        for g_vec in Pn_minus_1_pts:
        # for g_vec in Pn_minus_1_pts:
            if all([vector[i] >= g_vec[i] for i in range(n+num_poly)]):
                g = g_vec
        vector = [vector[i] - g[i] for i in range(n+num_poly)]
        u = vector
        g = c * vector_to_monomial(g)
        return_list.append([u,g])
    return return_list


# set of matrices of size |Pn_minus_1| x |Pn_minus_1|
Ruv_const_dict = {}
Ruv_u_dict = {}

def Ruv_const_helper(v,g):
    gi = lift_poly(vector_to_monomial(v)*g)
    h = sum([gi[i] + gens[i] * (gi[i]).derivative(gens[i]) for i in range(n+num_poly)])
    return to_pn_minus_1_basis(h)

def Ruv_u_helper(v,g):
    gi = lift_poly(vector_to_monomial(v)*g)
    return_ls = [gi[i] for i in range(n+num_poly)]
    return [to_pn_minus_1_basis(a) for a in return_ls]

def compute_Ruv(v):
    print('computing matrix for',v)
    Ruv_u_mat = [list(matrix(size_pn_minus_1)) for i in range(n+num_poly)]
    Ruv_const_mat = list(matrix(size_pn_minus_1))

    for i in range(size_pn_minus_1):
        g = Pn_minus_1[i]
        temp = Ruv_u_helper(v,g)
        for j in range(n):
            print(j)
            Ruv_u_mat[j][i] = temp[j]
        Ruv_const_mat[i] = Ruv_const_helper(v,g)
    Ruv_const_dict[tuple(v)] = matrix(Ruv_const_mat)
    Ruv_u_dict[tuple(v)] = tuple([matrix(Ruv_u_mat[i]) for i in range(n+num_poly)])
    return



def reduce_griffiths_dwork(u,g):
    g_vec = vector(matrix(to_pn_minus_1_basis(g)))
    # todo: speed up!
    while (u not in P_1_pts):
        print(u)
        best_k = -1
        best_v = -1
        for v in P_1_pts:
            k = max(u) # garbage value
            for i in range(n+num_poly):
                if v[i]>0:
                    k = min(k, u[i]//v[i])
            if k > best_k:
                best_k = k
                best_v = v
        v = best_v
        k = best_k
        # todo: check this
        if pole_order(vector_to_monomial([u[i] - k*v[i] for i in range(n+num_poly)])) == 0:
            # print("error: overcounted")
            k -= 1
        # print(u,v,k)

        if not tuple(v) in Ruv_u_dict.keys():
            compute_Ruv(v)
        C = Ruv_const_dict[tuple(v)]
        D = Ruv_u_dict[tuple(v)]

        E = C + sum([u[i]* D[i] for i in range(n+num_poly)])
        F_ = sum([v[i] * D[i] for i in range(n+num_poly)])

        for ind in range(1,k+1):
            # left g -> g A
            g_vec *= (E-ind*F_)
            # g_vec *=  (C + sum([(u[i]-ind*v[i]) * D[i] for i in range(n)]))
        u = [u[i] - k*v[i] for i in range(n+num_poly)]
    g = from_pn_minus_1_basis(vector(g_vec))
    
    return u,g


reduction_dict = {}

frob_matrix = [[0 for i in range(len(B))] for j in range(len(B))]

for i in range(len(B)):
    h = frobenius_on_cohom(i,prec)
    htemp = 0
    for u,g in to_uvg(h):
        print(u,g)
        # todo: can deduce degree u
        # todo: can just keep track of denom as (denom % p^prec) * p^something
        denom = factorial(pole_order(vector_to_monomial(u))-1+n)

        # this is the slow step
        u,g = reduce_griffiths_dwork(u,g)

        htemp += vector_to_monomial(u) * g // denom
    h = htemp
    
    summer = R(0)
#   monomials to reduce (with coefficients)
    monomial_list = [R(h.monomial_coefficient(monomial)) * monomial for monomial in h.monomials()]
    while len(monomial_list) > 0:
        term  = monomial_list.pop()
        if term not in QQ:  
            if len(term.monomials())>1:
                print('error: too many terms')
            monomial = R(term.monomials()[0])
            if not monomial in B:
                if not monomial in reduction_dict.keys():
                    q = J(monomial).lift()
                    r = monomial - q
                    l = r.lift(I)
                    temp = sum([l[i].derivative(gens[i]) for i in range(n+num_poly)])
                    reduction_dict[monomial] = temp + q
                result = reduction_dict[monomial]
                for _ in result.monomials():
                    monomial_list.append(_*term.monomial_coefficient(monomial) * result.monomial_coefficient(_))
            else:
                summer += term.monomial_coefficient(monomial) * monomial * factorial(pole_order(monomial))
        else:
            summer += term
    
    for j in range(len(B)):
        frob_matrix[i][j] = summer.monomial_coefficient(R(B[j])) #% p^prec
        
    print(B[i],summer)

if all([all([a.ord(p)>=0 for a in b]) for b in frob_matrix]):
    frob_matrix = [[a % p^ prec for a in b] for b in frob_matrix]
    print(frob_matrix)
else:
    print("Warning: non-invertible elements encountered")
frob_matrix = matrix(frob_matrix)

poly = R(frob_matrix.characteristic_polynomial())
print(poly)

# todo: mod these based on weil conjectures
poly = sum([(poly.monomial_coefficient(mono) % p^(prec) )* mono if (poly.monomial_coefficient(mono)% p^prec) < (p^prec)//2 else (-p^(prec)+poly.monomial_coefficient(mono)% p^prec)*mono for mono in poly.monomials()])
print(poly)