# number of variables
n = 3

# number of hypersurfaces in complete intersection
num_poly = 2

x_vars = list('x_%d' % i for i in range(n))
y_vars = list('y_%d' % i for i in range(num_poly))

R = PolynomialRing(QQ, names= x_vars + y_vars) # M + F call this S
gens = R.gens()

x_vars = gens[:n]
y_vars = gens[n:]


# for convenience
x = gens[0]
y = gens[1]
z = gens[2]

f = x*y + y^2 - z^2
g = x^3+y^3+z^3

f_0 = f
f_1 = g

poly_list = [f_0,f_1]
# poly_list = list(var('f_%d' % i) for i in range(num_poly))



# d_i = deg(f_i)

# bigrading is
# deg(x_i) = (0,1)
# deg(y_i) = (1,-deg(f_i))

# Z = V(f0,..,fc)
# then canonical bundle of Z w_z \cong O_Z(m) with
# m = sum(di) - n - 1

F = sum([y_vars[i] * poly_list[i] for i in range(num_poly)])
I  = F.jacobian_ideal()
J = R.quotient_ring(I)


