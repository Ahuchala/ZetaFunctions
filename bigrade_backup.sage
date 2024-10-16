# import itertools
# import numpy as np

# USE_CYTHON = True
USE_CYTHON = False

# poly_list = list(var('f_%d' % i) for i in range(num_poly))

# Z = V(f1,..,fc)
# then canonical bundle of Z w_z \cong O_Z(m) with
# m = sum(di) - n - 1

# J_p,m = H^{n-p,p}_van(X)

# d_i = deg(f_i)

# bigrading is
# deg(x_i) = (0,1)
# deg(y_i) = (1,-deg(f_i))


# number of hypersurfaces in complete intersection

# num_poly = 1
num_poly = 2

p = 7
# p = 10193
prec = 2


load("mat_mul.sage")


R.<x_0,x_1,x_2,x_3,y_1,y_2> = QQ[]
# R.<x_0,x_1,x_2,y_1> = QQ[]
# R.<x_0,x_1,x_2,x_3,y_1> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,y_1> = QQ[]

# R.<x_0,x_1,x_2,x_3,y_1,y_2> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,y_1,y_2> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,x_5,y_1,y_2> = QQ[]


gens = R.gens()
n = len(gens) - num_poly-1
# assert len(gens) == num_poly + n+1

x_vars = gens[:n+1]
y_vars = gens[n+1:]

# J_p,m consists of p copies of y_i and sum_d_i copies of x_j
# maybe compute by first all monomials in J_p,m for fixed p, then finding a basis

# f =  -5*x_0^3-x_0*x_1^2+4*x_0^2*x_2+2*x_0*x_1*x_2-3*x_1^2*x_2-2*x_0*x_2^2-6*x_1*x_2^2+5*x_2^3

f = sum([gen^2 for gen in x_vars])
# f = x_0^2*x_1^2 - 4*x_0^3*x_2 - 4*x_1^3*x_2 - 8*x_2^4 + 2*x_0*x_1*x_2*x_3 + x_2^2*x_3^2 - 4*x_0*x_1^2*x_4 - 4*x_0*x_2^2*x_4 - 4*x_3^3*x_4 + 2*x_0*x_1*x_4^3 + 2*x_2*x_3*x_4^2 + x_4^4

# f = x_0^3 + x_1^3 + x_2^3 - x_0*x_1*x_2 + x_3^3
# f = x_0^2 + x_1^2 + x_2^2 + x_3^2
g = x_0^2 + 2*x_1^2 + 3*x_2^2 + 4*x_3^2
# h = x_0*x_1 + x_1*x_2 + x_2*x_3

# f = x_0^3 + x_1^3 + x_2^3 - x_0*x_1*x_2

# f = x_0^2*x_1^2 - 4*x_0^3*x_2 - 4*x_1^3*x_2 - 8*x_2^4 + 2*x_0*x_1*x_2*x_3 + x_2^2*x_3^2 - 4*x_0*x_1^2*x_4 - 4*x_0*x_2^2*x_4 - 4*x_3^3*x_4 + 2*x_0*x_1*x_4^3 + 2*x_2*x_3*x_4^2 + x_4^4
# f = x_0^4+x_0^3*x_1+x_0*x_1^3+x_0^2*x_1*x_2+x_0*x_1^2*x_2+x_1^3*x_2+x_0^2*x_2^2+x_0*x_2^3+x_1*x_2^3+x_0*x_1*x_2*x_3+x_2^3*x_3+x_0^2*x_3^2+x_0*x_1*x_3^2+x_1*x_2*x_3^2+x_0*x_3^3+x_1*x_3^3+x_2*x_3^3
# f = x_0^2+2*x_0*x_1+2*x_1^2-x_0*x_2+x_1*x_2-2*x_0*x_3-x_1*x_3+2*x_2*x_3-2*x_3^2+x_0*x_4+2*x_2*x_4+2*x_3*x_4
# g = -2*x_0^3-2*x_0^2*x_1+x_0*x_1^2-x_1^3+2*x_0^2*x_2-x_0*x_1*x_2-2*x_1^2*x_2-2*x_1*x_2^2+x_0^2*x_3+2*x_0*x_1*x_3+2*x_1^2*x_3+2*x_0*x_2*x_3-x_1*x_2*x_3+x_0*x_3^2+x_2*x_3^2-x_0^2*x_4+2*x_1^2*x_4-2*x_0*x_2*x_4+x_1*x_2*x_4-2*x_2^2*x_4-2*x_2*x_3*x_4+2*x_0*x_4^2-x_2*x_4^2-2*x_4^3
# f = sum([gens[i]^2 for i in range(len(gens[:n+1]))])

poly_list = [f]
if num_poly == 2:
    f_0 = f; f_1 = g;
    poly_list = [f_0,f_1]
# poly_list = [f,g,h]

# todo: check if no rational points?

d = [_.degree() for _ in poly_list]
m = sum([f.degree() for f in poly_list]) - n - 1

F = sum([y_vars[i] * poly_list[i] for i in range(num_poly)])






assert p in Primes(), f'Warning: {p} is not prime'





def monomial_to_vector(m):
    return list(m.exponents()[0])

def vector_to_monomial(v):
    return prod([gens[i]^v[i] for i in range(n+num_poly+1)])


def monomial_degree(m):
#   e looks like  [(1, 1, 1, 1, 1)]
    e = m.exponents()[0]
#   x_i contribute (0, sum(e[:n]))
#   y_i contribute (sum(e[n:],sum([-poly_list[i].degree() for i in range(num_poly)] ))
    return [sum(e[n+1:]),  sum(e[:n+1]) +     sum([e[i+n+1]*-d[i] for i in range(num_poly)] )]

# each monomial is of the form m / F^j; returns j
def pole_order(m):
    e = m.exponents()[0]
    return sum(e[n+1:])


# Macaulay2 code to check smoothness
s = f'''
k = ZZ/{p};
R = k[x_0..x_{n}];
I = ideal {str(tuple(poly_list))[:-2] if len(poly_list) == 1 else str(tuple(poly_list))[:-1] });
J = I + minors({num_poly}, jacobian I);
saturate J == R
'''
assert str(macaulay2(s)) == "true", f'Warning: F not smooth modulo {p}'


# Macaulay2 code to run to compute Griffiths ring
s = f'''
k = ZZ/{p};
R = k[x_0..x_{n},y_1..y_{num_poly}, Degrees=>{{{n+1}:{{0,1}},
{str([(1,-i) for i in d])[1:-1].replace('(','{').replace(')','}')}}}];
F={sum([y_vars[i] * poly_list[i] for i in range(num_poly)])};
J = R/ideal jacobian F;
for i from 0 to {n-num_poly} list toString basis({{i,{m}}}, J)
'''
t = str(macaulay2(s))
t = t.replace("{","").replace("}","").replace("^","**").replace(" ","").split("matrix")[1:]
t = "".join(t).split(",")
B = []
for _ in t:
    if _ != "":
        eval("B.append(R(" + _ +"))")
print(B)
# assert False

max_cohomology_pole_order = max([pole_order(_) for _ in B])

# Macaulay2 code to run to compute P1
# Note that this is actually U_{1,0} rather than U_{1,m}
# since we need U_{1,0}*P_n in P_{n+1}
s = f'''
k = ZZ/{p};
R = k[x_0..x_{n},y_1..y_{num_poly}, Degrees=>{{{n+1}:{{0,1}},
{str([(1,-i) for i in d])[1:-1].replace('(','{').replace(')','}')}}}];
toString basis({{1,{0}}}, R)
'''
t = str(macaulay2(s))
t = t.replace("{","").replace("}","").replace("^","**").replace(" ","").split("matrix")[1:]
t = "".join(t).split(",")
P1 = []
for _ in t:
    if _ != "":
        eval("P1.append(R(" + _ +"))")
# print(P1)

P1_pts = [monomial_to_vector(_) for _ in P1]

# Macaulay2 code to run to compute Pn
s = f'''
k = ZZ/{p};
R = k[x_0..x_{n},y_1..y_{num_poly}, Degrees=>{{{n+1}:{{0,1}},
{str([(1,-i) for i in d])[1:-1].replace('(','{').replace(')','}')}}}];
toString basis({{{n},{-2}}}, R)
'''
# toString basis({{{n},{m}}}, R)

# seems like it should be max_cohomology_pole_order instead of n
t = str(macaulay2(s))
t = t.replace("{","").replace("}","").replace("^","**").replace(" ","").split("matrix")[1:]
t = "".join(t).split(",")
Pn = []
for _ in t:
    if _ != "":
        eval("Pn.append(R(" + _ +"))")
# print(Pn)
size_pn = len(Pn)
Pn_pts = [monomial_to_vector(mon) for mon in Pn]



xI = R.ideal([gen *F.derivative(gen) for gen in gens])


I = F.jacobian_ideal()
J = R.quotient_ring(I)


# returns frobenius of polynomial g
# e.g x + y -> x^p + y^p
def sigma(g):
    g_dict = g.dict()
    g_dict_keys = set(g_dict.keys())
    for key in g_dict_keys:
        g_dict[tuple([k*p for k in key])] = g_dict.pop(key)
    return R(g_dict)
    # return sum([g.coefficient(monomial)* monomial^p for monomial in g.monomials()])

# todo: make this return a dict instead of a polynomial
# returns frobenius of logarithmic form g/f^d * omega
def frobenius(g,prec=2):
    d = pole_order(g)
    summer = 0
    sigma_g = sigma(g)
    fj = R(1)
    for j in range(prec):
#        cacheing may be appropriate
        numer = binomial(-d,j) * binomial(d+prec-1,d+j) * sigma_g*sigma(fj)
        summer += numer
        fj *= F
    return summer


def frobenius_on_cohom(i,prec = 2):
    g = R(B[i])*prod(gens)
    g = frobenius(g,prec)
    return R(g *p^(n -1) / prod(gens))


reduction_dict = {}

frob_matrix = [[0 for i in range(len(B))] for j in range(len(B))]

# in: a polynomial, typically ouput of frobenius_on_cohom
# out: a list of tuples representing polynomials x^u g
#     with g in Pn
def to_ug(frobenius_of_Bi):
    return_list = []
    # todo: make optimal choices of which g

    # (x_0^5*x_1*y_1^2).dict() looks like {(5, 1, 0, 2): 1}
    hdict = frobenius_of_Bi.dict()

    

    while hdict: # is_nonempty
        hdict_keys = hdict.keys()
        hdict_keys_ls = list(hdict_keys)

    # for etuple in hdict.keys():
        etuple = hdict_keys_ls[0]

        vector = list(etuple)
        c = hdict[etuple]
        # todo: save on divisibility?

        # python uses function scope?
        # g = Pn_pts[0]
        for g_vec in Pn_pts:
            if all([vector[i] >= g_vec[i] for i in range(n+num_poly+1)]):
                g = g_vec
                break
        u = [vector[i] - g[i] for i in range(n+num_poly+1)]

        # g = R(0)
        g =  vector_to_monomial(g) * hdict[etuple]
        hdict.pop(etuple)


        # g =  vector_to_monomial(g)
        # assert g != 0
        # assert g in Pn
        # print(g)
        # g *= c
        # for etuple_2 in hdict_keys_ls:
        #     temp_vec = [etuple_2[i]-u[i] for i in range(n+num_poly+1)]
        #     if all([temp_vec[i] >= 0 for i in range(n+num_poly+1)]):
        #         print('popped')
        #         g += hdict[etuple_2] * vector_to_monomial(temp_vec)
        #         hdict.pop(etuple_2)

        return_list.append([u,g])
    return return_list

def to_pn_basis(g):
    return_vec = size_pn * [0]
    for monomial in g.monomials():
        ind = Pn.index(monomial)
        return_vec[ind] = g.monomial_coefficient(monomial)
    return return_vec

def from_pn_basis(g_vec):
    return sum([g_vec[i] * Pn[i] for i in range(size_pn)])

lift_dict = {}

# returns g.lift(I)
def lift_poly(g):
    summer = (n+num_poly+1) * [0]
    for monomial in g.monomials():

        if not monomial in lift_dict.keys():
            c = J(monomial).lift()
            r = (monomial-c).lift(xI)
            term = [0 for _ in range(n+num_poly+1)]
            for i in range(n+num_poly+1):
                term_i = r[i].dict()
                for _ in term_i.keys():
                    term_i[_] = (term_i[_]) % p^(2+prec) # this converts to ZZ/p^prec
                term[i] = R(term_i)
            lift_dict[monomial] = term
        monomial_lift = lift_dict[monomial]
        c = g.monomial_coefficient(monomial)
        for i in range(n+num_poly+1):
            summer[i] += c*monomial_lift[i]
    return summer   

Ruv_const_dict = {}
Ruv_u_dict = {}

def Ruv_const_helper(v,g):
    gi = lift_poly(vector_to_monomial(v)*g)
    h = sum([gi[i] + gens[i] * (gi[i]).derivative(gens[i]) for i in range(n+num_poly+1)])
    # print(h)
    return to_pn_basis(h)

def Ruv_u_helper(v,g):
    print(v,g)
    print(vector_to_monomial(v)*g)
    gi = lift_poly(vector_to_monomial(v)*g)
    # return_ls = [gi[i] for i in range(n+num_poly+1)]
    # return [to_pn_basis(a) for a in return_ls]
    # print(gi)
    return [to_pn_basis(a) for a in gi]


def compute_Ruv(v):
    Ruv_u_mat = [list(matrix(size_pn)) for i in range(n+num_poly+1)]
    Ruv_const_mat = list(matrix(size_pn))

    for i in range(size_pn):
        g = Pn[i]
        temp = Ruv_u_helper(v,g)
        for j in range(n+num_poly+1):
            Ruv_u_mat[j][i] = temp[j]
        Ruv_const_mat[i] = Ruv_const_helper(v,g)
    Ruv_const_dict[tuple(v)] = matrix(Ruv_const_mat)
    Ruv_u_dict[tuple(v)] = tuple([matrix(Ruv_u_mat[i]) for i in range(n+num_poly+1)])
    return

# writes u as u' + kv with k maximal and v in P1
def compute_vk(u):
    best_k = -1
    best_v = -1
    for v in P1_pts:
        k = max(u) # garbage value
        for i in range(n+num_poly+1):
            if v[i]>0:
                k = min(k, u[i]//v[i])
        if k > best_k:
            best_k = k
            best_v = v
    v = best_v
    k = best_k
    # I think I could just count the y_i's
    if pole_order(vector_to_monomial([u[i] - k*v[i] for i in range(n+num_poly+1)])) == 0: #m?
        # print("error: overcounted")
        k -= 1
    return v,k

def reduce_griffiths_dwork(u,g):
    g_vec = vector(matrix(to_pn_basis(g)))

    # todo: work out precise bounds!
    # g_vec = np.matrix(to_pn_basis(g),dtype=np.int32)

    # todo: speed up!
    while(pole_order(vector_to_monomial(u)))>max_cohomology_pole_order:
    # while (u not in P1_pts):
        print(u)

        v,k = compute_vk(u)

        if not tuple(v) in Ruv_u_dict.keys():
            compute_Ruv(v)
        C = Ruv_const_dict[tuple(v)]
        D = Ruv_u_dict[tuple(v)]

        E = C + sum([u[i]* D[i] for i in range(n+num_poly+1)])
        F = sum([v[i] * D[i] for i in range(n+num_poly+1)])

        g_vec = mat_mul(E,F,k,g_vec)
        # for ind in range(1,k+1):
            # left g -> g A
            # g_vec *= np.matrix(E-ind*F,dtype=np.int32)
        u = [u[i] - k*v[i] for i in range(n+num_poly+1)]
    # print(g_vec)
    # g_list = [R(_) for _ in g_list]
    # g = from_pn_basis(g_vec.tolist()[0])
    g = from_pn_basis(g_vec)

    
    return u,g

# todo: remove zero matrices which seem to all correspond to y_i

for i in range(len(B)):
    h = frobenius_on_cohom(i,prec)
    htemp = 0
    for u,g in to_ug(h):
        denom = factorial(pole_order(vector_to_monomial(u))+pole_order(g))

        print(u,g)
        # this is the slow step
        u,g = reduce_griffiths_dwork(u,g)


        htemp += vector_to_monomial(u) * g // denom
    h = htemp

    summer = R(0)
    monomial_list = [R(h.monomial_coefficient(monomial)) * monomial for monomial in h.monomials()]
    while len(monomial_list) > 0:
        term  = monomial_list[0]
        monomial_list.remove(term)
        # print(term)
        if term not in QQ:  
            print(term)
            if len(term.monomials())>1:
                print('error: too many terms')
            monomial = R(term.monomials()[0])
            if not monomial in B:
                if not monomial in reduction_dict.keys():
                    q = J(monomial).lift()
                    r = monomial - q
                    l = r.lift(I)
                    temp = sum([l[i].derivative(gens[i]) for i in range(n+num_poly+1)])
                    reduction_dict[monomial] = temp + q
                # this should be implemented better
                result = term.monomial_coefficient(monomial)* reduction_dict[monomial]
                h = sum(monomial_list) + result
                monomial_list = [R(h.monomial_coefficient(monomial)) * monomial for monomial in h.monomials()]
                
            else:
                summer += term
        else:
            summer += term
    
    for j in range(len(B)):
        frob_matrix[i][j] = summer.monomial_coefficient(R(B[j])) * factorial(pole_order(B[j])+num_poly -1) #% p^prec
        
    print(B[i],summer)

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
print(size_pn)