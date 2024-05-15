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

# R = QQ[x_0..x_3,y_0,y_1, Degrees=>{4:{0,1},{1,-2},{1,-2}}]
# f = x_0^2 + 2*x_1^2 - x_2^2 + 7*x_3^2
# g = x_0*x_1 + 2*x_0*x_2 - 4*x_1*x_3+x_2^2
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

p = 2
prec = 4


R.<x_0,x_1,x_2,y_0> = QQ[]
# R.<x_0,x_1,x_2,x_3,y_0,y_1> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,y_0,y_1> = QQ[]

gens = R.gens()

x_vars = gens[:n]
y_vars = gens[n:]

# J_p,m consists of p copies of y_i and sum_d_i copies of x_j
# maybe compute by first all monomials in J_p,m for fixed p, then finding a basis

# f = x_0^3 + x_1^3 + x_2^3 - x_0*x_1*x_2
# f = x_0^2 + x_1^2 + x_2^2 + x_3^2
# g = x_0^2 + 2*x_1^2 + 3*x_2^2 + 4*x_3^2
# smooth mod 5,7,11,13,17,19,23
# not 2,3


f = x_0^2*x_1+x_1^3+x_0*x_1*x_2+x_1*x_2^2+x_2^3
# f = x_0^2+2*x_0*x_1+2*x_1^2-x_0*x_2+x_1*x_2-2*x_0*x_3-x_1*x_3+2*x_2*x_3-2*x_3^2+x_0*x_4+2*x_2*x_4+2*x_3*x_4
# g = -2*x_0^3-2*x_0^2*x_1+x_0*x_1^2-x_1^3+2*x_0^2*x_2-x_0*x_1*x_2-2*x_1^2*x_2-2*x_1*x_2^2+x_0^2*x_3+2*x_0*x_1*x_3+2*x_1^2*x_3+2*x_0*x_2*x_3-x_1*x_2*x_3+x_0*x_3^2+x_2*x_3^2-x_0^2*x_4+2*x_1^2*x_4-2*x_0*x_2*x_4+x_1*x_2*x_4-2*x_2^2*x_4-2*x_2*x_3*x_4+2*x_0*x_4^2-x_2*x_4^2-2*x_4^3


B = [R(1), x_2^3*y_0] # need to be careful to use a basis sage likes
# B = [R(1), x_3^2*y_1]
# B = [R(1),x_0*x_1*x_4*y_1, x_0*x_2*x_4*y_1, x_0*x_3^2*y_1, x_0*x_3*x_4*y_1, x_0*x_4^2*y_1, x_1^2*x_4*y_1,x_1*x_2*x_4*y_1, x_1*x_3^2*y_1, x_1*x_3*x_4*y_1, x_1*x_4^2*y_1, x_2^2*x_3*y_1, x_2^2*x_4*y_1, x_2*x_3^2*y_1,x_2*x_3*x_4*y_1, x_2*x_4^2*y_1, x_3^3*y_1, x_3^2*x_4*y_1, x_3*x_4^2*y_1, x_4^3*y_1,x_4^6*y_1^2]
# B = [R(1),x_0*x_2^2*y_0]


poly_list = [f]
# f_0 = f; f_1 = g;
# poly_list = [f_0,f_1]



d = [_.degree() for _ in poly_list]
m = sum([f.degree() for f in poly_list]) - n #- 1

F = sum([y_vars[i] * poly_list[i] for i in range(num_poly)])

# smoothness check
# assert R.ideal([F.derivative(_) for _ in gens(R)]).ngens() == n + num_poly


I = F.jacobian_ideal()
J = R.quotient_ring(I)

# xI = R.ideal([_*F.derivative(_) for _ in R.gens()])# + [f])
# xJ = R.quotient(xI)

# toString basis({i,0},R)

def monomial_to_vector(m):
    return list(m.exponents()[0])

def vector_to_monomial(v):
    return prod([gens[i]^v[i] for i in range(n+num_poly)])


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
        numer = binomial(-d,j) * binomial(d+prec-1,d+j) * sigma_g*sigma(fj) / factorial(pole_order(sigma_g) + pole_order(sigma(fj))-1)
        summer += numer
        fj *= F
    return summer


def frobenius_on_cohom(i,prec = 2):
    g = R(B[i])*prod(gens)
    g = frobenius(g,prec)
    return R(g *p^(n -2) / prod(gens))


reduction_dict = {}

frob_matrix = [[0 for i in range(len(B))] for j in range(len(B))]

for i in range(len(B)):
    h = frobenius_on_cohom(i,prec)
    summer = R(0)
    monomial_list = [R(h.monomial_coefficient(monomial)) * monomial for monomial in h.monomials()]
    while len(monomial_list) > 0:
        term  = monomial_list[0]
        monomial_list.remove(term)
        print(term)
        if term not in QQ:  
            # print(term)
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
                result = term.monomial_coefficient(monomial)* reduction_dict[monomial]
                # result = sum([_*term.monomial_coefficient(monomial) * result.monomial_coefficient(_) for _ in result.monomials()])
                # for _ in result.monomials():
                    # monomial_list.append(_ * result.monomial_coefficient(_))
                h = sum(monomial_list) + result
                monomial_list = [R(h.monomial_coefficient(monomial)) * monomial for monomial in h.monomials()]
                
            else:
                summer += term
    #             print(term)
        else:
            summer += term
    
    for j in range(len(B)):
        frob_matrix[i][j] = summer.monomial_coefficient(R(B[j])) * factorial(pole_order(B[j]) + num_poly-1) # !!?!?!?!?!/
        
    print(B[i],summer)
    # assert summer == sum([frob_matrix[i][_]*B[_] for _ in range(len(B))])

if all([all([a.ord(p)>=0 for a in b]) for b in frob_matrix]):
    frob_matrix = [[a % p^ prec for a in b] for b in frob_matrix]
    print(frob_matrix)
else:
    print("Warning: non-invertible elements encountered")
frob_matrix = matrix(frob_matrix)

poly = (frob_matrix.characteristic_polynomial())
print(poly)

# todo: mod these based on weil conjectures
poly = sum([(poly.monomial_coefficient(mono) % p^(prec) )* mono if (poly.monomial_coefficient(mono)% p^prec) < (p^prec)//2 else (-p^(prec)+poly.monomial_coefficient(mono)% p^prec)*mono for mono in poly.monomials()])
print(poly)
print(p)