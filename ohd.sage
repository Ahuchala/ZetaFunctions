from sympy.utilities.iterables import multiset_permutations

DEBUG = True

p = 53

prec = 3

if DEBUG:
    assert(is_prime(p))

# todo: modular arithmetic
# R.<x,y,z> = QQ[]
R.<w,x,y,z> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4> = QQ[]


# weights = [3,1,1,1]
weights = [1,1,1,1]
# weights = [1,1,1]
# weights = [1,3,1]

# weights = [11,14,18,20,25] # --> correct from typo in example 7.2

Rgens = R.gens()
n = len(Rgens) #number of variables




# f = x^3 - y^3 - z^3 - x*y*z
# f = x^4 + y^4 + z^4-x*y*z^2+4*x^2*z^2
# f = x^5 + y^5 + z^5 - x*y*z^3
# f = x^4 + y^4 + z^4-x*y*z^2+4*x^2*z^2
# f = y^2 - x^6 - z^6-x^3*z^3
# f = y^2 -(-2*x^6-x^5*z+3*x^4*z^2+x^3*z^3-2*x^2*z^4+x*z^5+3*z^6)
# f = (2)*x^3+3*x^2*y+(4)*x*y^2+(5)*y^3+(5)*x^2*z+x*y*z+(7)*y^2*z+(5)*x*z^2+(7)*y*z^2+(1)*z^3

# f = x_0^8 + x_1^5 * x_2 + x_0^2 * x_1^2 *x_2*x_3 + x_1*x_2^3*x_3 + x_1^2*x_3^3 + x_0*x_1*x_2*x_3*x_4+x_2*x_3*x_4^2

# f =w^2 - x^6 - y^6 - z^6
# f = w^3 + x^3 +y^3 - z^3 - w*x*z+2*y*z^2
f = x^4 + y^4 + z^4 + w^4 - w * x * y *z
# f = -w^4-w^3*x-w^2*x^2-x^4-w^3*y-w^2*x*y-w*x^2*y+x^3*y+w^2*y^2+w*x*y^2+x^2*y^2-w*y^3+y^4+w^3*z-w^2*x*z-x^3*z-w^2*y*z+w*x*y*z-w*y^2*z+x*y^2*z-w^2*z^2-x^2*z^2-w*y*z^2+x*y*z^2-y^2*z^2+y*z^3+z^4
# f= 2*x_0^3+2*x_0*x_1^2+x_1^3+2*x_0^2*x_2-x_0*x_1*x_2+2*x_1^2*x_2+x_0*x_2^2+2*x_1*x_2^2+x_2^3-x_0^2*x_3-x_0*x_1*x_3-x_0*x_2*x_3+x_1*x_2*x_3+2*x_2^2*x_3+x_0*x_3^2-x_1*x_3^2-x_2*x_3^2+2*x_0^2*x_4+2*x_0*x_1*x_4-x_1^2*x_4-2*x_0*x_2*x_4-x_1*x_2*x_4+2*x_2^2*x_4+x_0*x_4^2-x_1*x_4^2-2*x_2*x_4^2
# f = w^4+2*w*x^3-2*x^4-x^3*y-x^2*y^2-y^4+w^3*z-x^3*z-2*w^2*y*z+2*w*x*y*z-x^2*y*z-w*y^2*z+2*x*y^2*z-2*y^3*z-w^2*z^2-2*w*x*z^2+x^2*z^2-2*w*y*z^2+x*y*z^2+y^2*z^2+2*w*z^3+2*x*z^3-2*y*z^3-2*z^4
I = R.ideal([f.derivative(_) for _ in R.gens()])# + [f])
J = R.quotient(I)

xI = R.ideal([_*f.derivative(_) for _ in R.gens()])# + [f])
xJ = R.quotient(xI)



# todo: check smoothness via polytope considerations?
if DEBUG:
    r = str(R.gens()).replace('(','').replace(')','').replace(' ','').split(',')
    X = toric_varieties.WP(weights, base_ring = GF(p), names = r)
    assert(X.subscheme(f).is_smooth())
    print("X is smooth")

# d = f.degree()
# fdegree = d

# i.e. y^3 -> [0,3,0,0] when n=4
def monomial_to_vector(m):
    return list(m.exponents()[0])

def vector_to_monomial(v):
    return prod([Rgens[i]^v[i] for i in range(n)])

d = sum([weights[i] * monomial_to_vector(f.monomials()[0])[i] for i in range(n)])
fdegree = d
if DEBUG:
    assert(all([sum([weights[i]*monomial_to_vector(a)[i] for i in range(n)])==fdegree for a in f.monomials()]))


load("aux_functions.sage")

vertices =list(matrix.identity(n-1))
for i in range(n-1):
    for j in range(n-1):
        vertices[i][j] *= fdegree//weights[i] #todo: double check
vertices = [list(a) for a in vertices]
vertices += [list(0 for _ in range(n-1))]
lattice_polytope = Polyhedron(vertices)
print(vertices)

P1 = lattice_polytope.integral_points()
P1 = [prod([Rgens[i]^mon[i] for i in range(n-1)])*Rgens[n-1]^((fdeg - sum([weights[i]*mon[i] for i in range(n-1)]))//weights[n-1]) for mon in P1]
P1_int = [a for a in P1 if lattice_polytope.interior_contains(a)]
P1_int = [prod([Rgens[i]^mon[i] for i in range(n-1)])*Rgens[n-1]^((fdeg - sum([weights[i]*mon[i] for i in range(n-1)]))//weights[n-1])//(prod(Rgens)) for mon in P1_int]

P1_pts = [monomial_to_vector(a) for a in P1]

lift_dict = {}
B= []

for l in range(1,n):
    Pd = lattice_polytope.dilation(l)
    Pd_int = [a for a in Pd.integral_points() if Pd.interior_contains(a)]
    for mon in Pd_int:
        mon = prod([Rgens[i]^mon[i] for i in range(n-1)])*Rgens[n-1]^((l*fdeg - sum([weights[i]*mon[i] for i in range(n-1)]))//weights[n-1])//(prod(Rgens))
        # check if in I
        c = J(mon).lift()
        r = (mon-c).lift(I)
        if c != 0:
            B += c.monomials()

B = list(Set(B))
print(len(B))
print(B)

Pn_minus_1 = set()
l = n-1+1
Pd = lattice_polytope.dilation(l)
Pd = [a for a in Pd.integral_points()]

for mon in Pd:
    mon = prod([Rgens[i]^mon[i] for i in range(n-1)])*Rgens[n-1]^((l*fdeg - sum([weights[i]*mon[i] for i in range(n-1)]))//weights[n-1])//(prod(Rgens))
    Pn_minus_1.add(mon)
Pn_minus_1.remove(0)
Pn_minus_1_pts = [monomial_to_vector(a) for a in Pn_minus_1]
size_pn_minus_1 = len(Pn_minus_1_pts)

Pn_minus_1_list = list(Pn_minus_1)

def to_pn_minus_1_basis(g):
    return_vec = size_pn_minus_1 * [0]
    for monomial in g.monomials():
        ind = Pn_minus_1_list.index(monomial)
        return_vec[ind] = g.monomial_coefficient(monomial)
    return return_vec

# todo: figure out why I need to cast to R
def from_pn_minus_1_basis(g_vec):
    return sum([g_vec[i] * Pn_minus_1_list[i] for i in range(size_pn_minus_1)])

    # return sum([R(g_vec[i]) * Pn_minus_1_list[i] for i in range(size_pn_minus_1)])

    
for part in Partitions((n-1)*d-sum(weights)+1+d, max_length=n):
#       pad with zeros
    part += [0] * (n - len(part))

    for perm in multiset_permutations(part):

        perm = [perm[i] // weights[i] for i in range(n)]
        m = vector_to_monomial(perm)
        if not m in lift_dict.keys():
#                 use cached lifts
            c = xJ(m).lift()
            r = (m-c).lift(xI)
            lift_dict[m] = r
            
print(len(B), (-1)^(n)*(1/d * ((1-d)^(n)-1)+1))


def sigma(g):
    return sum([g.coefficient(monomial)* monomial^p for monomial in g.monomials()])

def frobenius(g,prec=2):
    d = degree(g)
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
    g = B[i]*prod(Rgens)
    g = frobenius(g,prec)
    return R(g *p^(n -2) / prod(Rgens))

def lift_poly(g):
    summer = n * [0]
    for monomial in g.monomials():

        if not monomial in lift_dict.keys():
            c = xJ(monomial).lift()
            r = (monomial-c).lift(xI)
            lift_dict[monomial] = r
        monomial_lift = lift_dict[monomial]
        c = g.monomial_coefficient(monomial)
        for i in range(n):
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
            if all([vector[i] >= g_vec[i] for i in range(n)]):
                g = g_vec
        vector = [vector[i] - g[i] for i in range(n)]
        u = vector
        g = c * vector_to_monomial(g)
        return_list.append([u,g])
    return return_list


# set of matrices of size |Pn_minus_1| x |Pn_minus_1|
Ruv_const_dict = {}
Ruv_u_dict = {}

def Ruv_const_helper(v,g):
    gi = lift_poly(vector_to_monomial(v)*g)
    h = sum([gi[i] + Rgens[i] * (gi[i]).derivative(Rgens[i]) for i in range(n)])
    return to_pn_minus_1_basis(h)

def Ruv_u_helper(v,g):
    gi = lift_poly(vector_to_monomial(v)*g)
    return_ls = [gi[i] for i in range(n)]
    return [to_pn_minus_1_basis(a) for a in return_ls]

def compute_Ruv(v):
    Ruv_u_mat = [list(matrix(size_pn_minus_1)) for i in range(n)]
    Ruv_const_mat = list(matrix(size_pn_minus_1))

    for i in range(size_pn_minus_1):
        g = Pn_minus_1_list[i]
        temp = Ruv_u_helper(v,g)
        for j in range(n):
            Ruv_u_mat[j][i] = temp[j]
        Ruv_const_mat[i] = Ruv_const_helper(v,g)
    Ruv_const_dict[tuple(v)] = matrix(Ruv_const_mat)
    Ruv_u_dict[tuple(v)] = tuple([matrix(Ruv_u_mat[i]) for i in range(n)])
    return


def reduce_griffiths_dwork(u,g):
    g_vec = vector(matrix(to_pn_minus_1_basis(g)))
    # todo: speed up!
    while (u not in P1_pts):
        print(u)
        best_k = -1
        best_v = -1
        for v in P1_pts:
            k = max(u) # garbage value
            for i in range(n):
                if v[i]>0:
                    k = min(k, u[i]//v[i])
            if k > best_k:
                best_k = k
                best_v = v
        v = best_v
        k = best_k
        if degree_vector([u[i] - k*v[i] for i in range(n)]) == 0:
            # print("error: overcounted")
            k -= 1
        # print(u,v,k)

        if not tuple(v) in Ruv_u_dict.keys():
            compute_Ruv(v)
        C = Ruv_const_dict[tuple(v)]
        D = Ruv_u_dict[tuple(v)]

        E = C + sum([u[i]* D[i] for i in range(n)])
        F = sum([v[i] * D[i] for i in range(n)])

        for ind in range(1,k+1):
            # left g -> g A
            g_vec *= (E-ind*F)
            # g_vec *=  (C + sum([(u[i]-ind*v[i]) * D[i] for i in range(n)]))
        u = [u[i] - k*v[i] for i in range(n)]
    g = from_pn_minus_1_basis(vector(g_vec))
    
    return u,g


reduction_dict = {}

frob_matrix = [[0 for i in range(len(B))] for j in range(len(B))]

for i in range(len(B)):
    h = frobenius_on_cohom(i,prec)
    htemp = 0
    for u,g in to_uvg(h):
        # todo: can deduce degree u
        # todo: can just keep track of denom as (denom % p^prec) * p^something
        denom = factorial(degree(vector_to_monomial(u))-1+n)

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
                    temp = sum([l[i].derivative(Rgens[i]) for i in range(n)])
                    reduction_dict[monomial] = temp + q
                result = reduction_dict[monomial]
                for _ in result.monomials():
                    monomial_list.append(_*term.monomial_coefficient(monomial) * result.monomial_coefficient(_))
            else:
                summer += term.monomial_coefficient(monomial) * monomial * factorial(degree(monomial))
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
