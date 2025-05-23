# maybe we require f to be nondegenerate too?
# why does this have issues? p =13, f =  -5*x_0^3-x_0*x_1^2+4*x_0^2*x_2+2*x_0*x_1*x_2-3*x_1^2*x_2-2*x_0*x_2^2-6*x_1*x_2^2+5*x_2^3
# R/xI has has hilbert function {1, 3, 6, 7, 6, 3, 2, 2, 2, 2, 2}
# what gives?

# maybe check homological dimension of R / xI and compute accordingly

import subprocess
import re

# import itertools
# import numpy as np

# USE_CYTHON = True
USE_CYTHON = False


# whether to use QQ instead of Z/p^r precision in intermediate arithmetic, mainly for debugging
#USE_RATIONAL_ARITHMETIC = False
#USE_RATIONAL_ARITHMETIC = True

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
# num_poly = 3


# p = Primes().next(2^13)
p = 7
prec = 2# todo: work this out
arithmetic_precision_increase = 3 # todo: work this out too
prec_arithmetic = p^(arithmetic_precision_increase+prec)
arithmetic_ring = Integers(prec_arithmetic)

# load("mat_mul.spyx")
load("mat_mul.sage")
load("hirzebruch.sage")

R.<x_0,x_1,x_2,x_3,y_1,y_2> = QQ[]
# R.<x_0,x_1,x_2,y_1> = QQ[]
# R.<x_0,x_1,x_2,x_3,y_1> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,y_1> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,x_5,y_1> = QQ[]


#R.<x_0,x_1,x_2,x_3,y_1,y_2> = QQ[]
#R.<x_0,x_1,x_2,x_3,x_4,y_1,y_2> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,x_5,y_1,y_2> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,y_1,y_2,y_3> = QQ[]



gens = R.gens()
prod_gens = prod(gens)
n = len(gens) - num_poly-1
# assert len(gens) == num_poly + n+1

x_vars = gens[:n+1]
y_vars = gens[n+1:]

# J_p,m consists of p copies of y_i and sum_d_i copies of x_j
# maybe compute by first all monomials in J_p,m for fixed p, then finding a basis

# f = x_1^2*x_2-x_0^3 - x_0*x_2^2 #degenerate
# f = x_0*x_1 + x_2*x_3 + x_4^2
# g = 2*x_0^3+2*x_0^2*x_1+2*x_0^2*x_2+x_0*x_1*x_2-2*x_1^2*x_2-x_1*x_2^2-x_0*x_1*x_3+x_0*x_2*x_3-x_2^2*x_3+2*x_0*x_3^2+x_2*x_3^2-2*x_0^2*x_4-2*x_0*x_1*x_4+x_0*x_2*x_4-x_2^2*x_4+2*x_0*x_3*x_4-x_1*x_3*x_4+2*x_2*x_3*x_4+2*x_3^2*x_4+2*x_0*x_4^2-2*x_1*x_4^2+2*x_3*x_4^2-2*x_4^3
# f = sum(gen^3 for gen in x_vars)
g = sum(gen^2 for gen in x_vars)
# f = x_0^2*x_1^2 - 4*x_0^3*x_2 - 4*x_1^3*x_2 - 8*x_2^4 + 2*x_0*x_1*x_2*x_3 + x_2^2*x_3^2 - 4*x_0*x_1^2*x_4 - 4*x_0*x_2^2*x_4 - 4*x_3^3*x_4 + 2*x_0*x_1*x_4^3 + 2*x_2*x_3*x_4^2 + x_4^4

# f = x_1^2*x_2 - x_0^3 - 8*x_0*x_2^2 - 17*x_2^3

# f =  -5*x_0^3-x_0*x_1^2+4*x_0^2*x_2+2*x_0*x_1*x_2-3*x_1^2*x_2-2*x_0*x_2^2-6*x_1*x_2^2+5*x_2^3
# g = sum(gen^2 for gen in x_vars)
f = x_0^3 + x_1^3 + x_2^3 - x_0*x_1*x_2 + x_3^3
# g = x_0^2 + 2*x_1^2 + 3*x_2^2 + 4*x_3^2
# g = sum(x_vars[i] * x_vars[(i+1)%(n+1)] for i in range(n+1))+x_0^2
# h = sum(x_vars[i] * x_vars[(i+2)%(n+1)] for i in range(n+1))-x_1^2
# print(f,g,h)
# h = x_0*x_1 + x_1*x_2 + x_2*x_3

# f = x_0^3 + x_1^3 + x_2^3 - x_0*x_1*x_2 + x_3^3

# f = x_0^2*x_1^2 - 4*x_0^3*x_2 - 4*x_1^3*x_2 - 8*x_2^4 + 2*x_0*x_1*x_2*x_3 + x_2^2*x_3^2 - 4*x_0*x_1^2*x_4 - 4*x_0*x_2^2*x_4 - 4*x_3^3*x_4 + 2*x_0*x_1*x_4^3 + 2*x_2*x_3*x_4^2 + x_4^4
# f = x_0^4+x_0^3*x_1+x_0*x_1^3+x_0^2*x_1*x_2+x_0*x_1^2*x_2+x_1^3*x_2+x_0^2*x_2^2+x_0*x_2^3+x_1*x_2^3+x_0*x_1*x_2*x_3+x_2^3*x_3+x_0^2*x_3^2+x_0*x_1*x_3^2+x_1*x_2*x_3^2+x_0*x_3^3+x_1*x_3^3+x_2*x_3^3
# f = x_0^2+2*x_0*x_1+2*x_1^2-x_0*x_2+x_1*x_2-2*x_0*x_3-x_1*x_3+2*x_2*x_3-2*x_3^2+x_0*x_4+2*x_2*x_4+2*x_3*x_4
# g = -2*x_0^3-2*x_0^2*x_1+x_0*x_1^2-x_1^3+2*x_0^2*x_2-x_0*x_1*x_2-2*x_1^2*x_2-2*x_1*x_2^2+x_0^2*x_3+2*x_0*x_1*x_3+2*x_1^2*x_3+2*x_0*x_2*x_3-x_1*x_2*x_3+x_0*x_3^2+x_2*x_3^2-x_0^2*x_4+2*x_1^2*x_4-2*x_0*x_2*x_4+x_1*x_2*x_4-2*x_2^2*x_4-2*x_2*x_3*x_4+2*x_0*x_4^2-x_2*x_4^2-2*x_4^3
# f = sum(gens[i]^2 for i in range(len(gens[:n+1])))

poly_list = [f]
if num_poly == 2:
    f_0 = f; f_1 = g;
    poly_list = [f_0,f_1]
elif num_poly == 3:
    poly_list = [f,g,h]

# todo: check if no rational points?

degree_of_polynomials = [_.degree() for _ in poly_list]
m = sum(f.degree() for f in poly_list) - n - 1

F = sum(y_vars[i] * poly_list[i] for i in range(num_poly))




assert all(poly.is_homogeneous() for poly in poly_list), "Warning: at least one polynomial is not homogeneous"
assert p in Primes(), f'Warning: {p} is not prime'





# todo: maybe use .dict() instead?
def monomial_to_vector(m):
    return list(m.exponents()[0])

def vector_to_monomial(v):
    return prod(gens[i]^v[i] for i in range(n+num_poly+1))

def vector_addition(u,v):
    return tuple(u[i] + v[i] for i in range(n+num_poly+1))

# given constant and a dict of a polynomial g, returns (c * g).dict()
def multiplication_by_scalar(c,g_dict):
    return_dict = {}
    for g_key in g_dict:
        return_dict[g_key] = c * g_dict[g_key]
    return return_dict

# given two polynomial dicts, returns the dict of their product
def polynomial_multiplication(g_dict,h_dict):
    return_dict = {}
    # foil it...
    for g_key in g_dict:
        g_val = g_dict[g_key]
        for h_key in h_dict:
            return_dict[vector_addition(g_key,h_key)] = g_val * h_dict[h_key]
    return return_dict

# divides by prod(gens)
def divide_by_x0xn(g_dict):
    return_dict = {}
    for g_key in g_dict:
        return_dict[tuple(g_key_i-1 for i,g_key_i in enumerate(g_key))] = g_dict[g_key]
    return return_dict


def monomial_degree(monomial):
#   e looks like  [(1, 1, 1, 1, 1)]
    e = monomial.exponents()[0]
#   x_i contribute (0, e[i])
#   y_i contribute (e[i+n+1],-poly_list[i].degree()] ))
    return [sum(e[n+1:]), sum(e[:n+1]) + sum(e[i+n+1]*-degree_of_polynomials[i] for i in range(num_poly))]

# for a monomial in vector form
def pole_order_vector(vec):
    return sum(vec[n+1:])

# for a monomial in iterator form
# warning: this must only contain data of y coordinates
# def pole_order_iterator(_iter):
#     # assert len(list(_iter)) == num_poly, "Warning: monomial entered was not a proper iterator of y variables"
#     return sum(_iter)

# redundant method
def pole_order_monomial(monomial):
    return pole_order(monomial)

# each monomial is of the form m / F^j; returns j
def pole_order(monomial):
    e = monomial.exponents()[0]
    return sum(e[n+1:])

# run Macaulay2 programs


def run_macualay2_program(function_name, args):
    file_name = "compute_griffiths.m2"
    ls = ["M2", "--script", file_name, function_name]
    for arg in args:
        ls.append(str(arg))
    a = subprocess.run(ls,capture_output=True,check=True).stdout
    return a

regex_contains_true = re.compile(r"[true]+")

# does it contain the string true
# since s == "b'true\n'" seems not to work
def macaulay2_smooth_check_to_sage(s):    
    return regex_contains_true.search(s)


s = str(run_macualay2_program("assertSmooth", [n,poly_list,p]))
assert macaulay2_smooth_check_to_sage(s), f'Warning: F not smooth modulo {p}'



# first make ^ into **
regex_remove_hyphens = re.compile(r"(\^)")

# match characters in e.g. x_1**2*y_1
# exclude only 1 word match to rule out word matrix
# reallow decimal since 1 is valid
regex_select_monomials = re.compile(r"([x|y|\d|_|^\*]{2,}|\d+)")

# in: macaulay list of matrices
# out: sage list of matrices
def macaulay2_matrices_to_sage(s):    
    s = re.sub(regex_remove_hyphens,"**",s)
    return re.findall(regex_select_monomials,s)


# Macaulay2 code to run to compute Griffiths ring

a = str(run_macualay2_program("computeGriffithsRing", [n,poly_list,p]))
t = macaulay2_matrices_to_sage(a)
B = [R(_) for _ in t]
print(B)
size_B = len(B)
# assert False

hodge_numbers = h_pq(degree_of_polynomials,n)
hodge_slopes = h_pq_to_hodge_polygon(hodge_numbers)
# assert [pole_order(_) for _ in B] == hodge_slopes, f"Warning: issue with Hodge numbers {[pole_order(_) for _ in B]}, {hodge_slopes}"

max_cohomology_pole_order = max(hodge_slopes)

# Macaulay2 code to run to compute P1
# Note that this is actually U_{1,0} rather than U_{1,m}
# since we need U_{1,0}*P_n in P_{n+1}

a = str(run_macualay2_program("computeP1", [n,poly_list,p]))
t = macaulay2_matrices_to_sage(a)
P1 = [R(_) for _ in t]
# print(P1)

P1_pts = [monomial_to_vector(_) for _ in P1]

a = str(run_macualay2_program("computePn", [n,poly_list,p]))
t = macaulay2_matrices_to_sage(a)
Pn = [R(_) for _ in t]

# print(Pn)
size_pn = len(Pn)
Pn_pts = [monomial_to_vector(mon) for mon in Pn]

xI = R.ideal([gen *F.derivative(gen) for gen in gens])


I = F.jacobian_ideal()
J = R.quotient_ring(I)



def sigma(g,return_as_dict = False):
    return sigma_dict(g.dict(),return_as_dict)


# returns frobenius of polynomial dict g.dict()
# e.g x + y -> x^p + y^p
def sigma_dict(g_dict,return_as_dict = False):
    g_dict_keys = set(g_dict.keys())
    for key in g_dict_keys:
        g_dict[tuple(k*p for k in key)] = g_dict.pop(key) # if Fq is eventually implemented this will change
    if return_as_dict:
        return g_dict
    return R(g_dict)

# todo: make this return a dict instead of a polynomial
# returns frobenius of logarithmic form g/f^d * omega
def frobenius(g,prec=2,return_as_dict = True):
    d = pole_order(g)
    summer_dict = {}
    sigma_g_dict = sigma(g,True)
    fj = R(1)

    for j in range(prec):
#        cacheing may be appropriate
        const = binomial(-d,j) * binomial(d+prec-1,d+j)
        numer =  multiplication_by_scalar(const,polynomial_multiplication(sigma_g_dict,sigma(fj,True)))
        summer_dict |= numer
        # summer += numer
        fj *= F
    if return_as_dict:
        return summer_dict
    return R(summer_dict)


# returns as a dict
def frobenius_on_cohom(i,prec = 2):
    Bi_times_prod_gens = B[i]*prod_gens
    Bi_times_prod_gens = frobenius(Bi_times_prod_gens,prec)
    # return R(g *p^(n -1) / prod(gens))
    return multiplication_by_scalar(p^(n -1),divide_by_x0xn(Bi_times_prod_gens))

reduction_dict = {}

frob_matrix = [[0 for i in range(size_B)] for j in range(size_B)]

# in: a polynomial, typically ouput of frobenius_on_cohom
# out: a list of tuples representing polynomials x^u g
#     with g in Pn
def to_ug(frobenius_of_Bi):
    return_list = []
    # todo: make optimal choices of which g

    # (x_0^5*x_1*y_1^2).dict() looks like {(5, 1, 0, 2): 1}
    # hdict = frobenius_of_Bi.dict()
    hdict = frobenius_of_Bi

    while hdict: # is_nonempty
        hdict_keys = hdict.keys()
        hdict_keys_ls = list(hdict_keys)

    # for etuple in hdict.keys():
        etuple = hdict_keys_ls[0]

        vector = list(etuple)
        # todo: save on divisibility?

        # python uses function scope
        # g = Pn_pts[0]
        for g_vec in Pn_pts:
            if all(vector[i] >= g_vec[i] for i in range(n+num_poly+1)):
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
        #     if all(temp_vec[i] >= 0 for i in range(n+num_poly+1)):
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
    return sum(QQ(g_vec[i]) * Pn[i] for i in range(size_pn))

lift_dict = {}

# returns g.lift(I)
def lift_poly(g):
    summer = (n+num_poly+1) * [0]
    for monomial in g.monomials():

        if monomial not in lift_dict:
            c = J(monomial).lift()
            r = (monomial-c).lift(xI)
            term = [0 for _ in range(n+num_poly+1)]
            for i in range(n+num_poly+1):
                term_i = r[i].dict()
                # if not USE_RATIONAL_ARITHMETIC:
                    # for _ in term_i.keys():
                    # note: this will cause issues for some primes. why?
                        # term_i[_] = (term_i[_]) % p^(prec+arithmetic_precision_increase) # this converts to ZZ/p^prec
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
    h = sum(gi[i] + gens[i] * (gi[i]).derivative(gens[i]) for i in range(n+num_poly+1))
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
    if USE_RATIONAL_ARITHMETIC:
        Ruv_u_mat = [list(matrix(QQ,size_pn,size_pn)) for i in range(n+num_poly+1)]
        Ruv_const_mat = list(matrix(QQ,size_pn,size_pn))

    else:
        Ruv_u_mat = [list(matrix(arithmetic_ring,size_pn,size_pn)) for i in range(n+num_poly+1)]
        Ruv_const_mat = list(matrix(arithmetic_ring,size_pn,size_pn))



    for i in range(size_pn):
        g = Pn[i]
        temp = Ruv_u_helper(v,g)
        for j in range(n+num_poly+1):
            Ruv_u_mat[j][i] = temp[j]
        Ruv_const_mat[i] = Ruv_const_helper(v,g)
    if USE_RATIONAL_ARITHMETIC:
        Ruv_const_dict[tuple(v)] = matrix(Ruv_const_mat)
        Ruv_u_dict[tuple(v)] = tuple(matrix(Ruv_u_mat[i]) for i in range(n+num_poly+1))
    else:
        Ruv_const_dict[tuple(v)] = matrix(arithmetic_ring,Ruv_const_mat)
        Ruv_u_dict[tuple(v)] = tuple(matrix(arithmetic_ring,Ruv_u_mat[i]) for i in range(n+num_poly+1))
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
    if pole_order_vector([u[i] - k*v[i] for i in range(n+num_poly+1)]) == 0: #m?

    # if pole_order_iterator(u[i] - k*v[i] for i in range(n+1,n+num_poly+1)) == 0: #m?
        # print("error: overcounted")
        k -= 1
    return v,k

def reduce_griffiths_dwork(u,g):
    if USE_RATIONAL_ARITHMETIC:
        g_vec = vector(matrix(to_pn_basis(g)))
    else:
        g_vec = vector(matrix(arithmetic_ring,to_pn_basis(g)))


    # todo: work out precise bounds!

    # todo: speed up!
    while(pole_order_vector(u))>max_cohomology_pole_order: #add if poly not in ideal for bad primes
    # while (u not in P1_pts):
        print(u)

        v,k = compute_vk(u)

        if tuple(v) not in Ruv_u_dict:
            compute_Ruv(v)
        C = Ruv_const_dict[tuple(v)]
        D = Ruv_u_dict[tuple(v)]

        E = C + sum(u[i]* D[i] for i in range(n+num_poly+1))
        F = sum(v[i] * D[i] for i in range(n+num_poly+1))

        g_vec = mat_mul(E,F,k,g_vec)
 
        u = [u[i] - k*v[i] for i in range(n+num_poly+1)]

    # if not USE_RATIONAL_ARITHMETIC:
        # g = R(g)
    g = from_pn_basis(g_vec)

    
    return u,g

# todo: remove zero matrices which seem to all correspond to y_i

reduce_to_B_dict = {}

def compute_new_reduce_to_B_key(key):
    monomial_list = [vector_to_monomial(key)]
    summer = R(0)
    while len(monomial_list) > 0:
        term  = monomial_list[0]
        monomial_list.remove(term)
        if term not in QQ:  
            print(term)
            monomial = R(term.monomials()[0])
            if monomial not in B:
                if monomial not in reduction_dict:
                    q = J(monomial).lift()
                    r = monomial - q
                    l = r.lift(I)
                    temp = sum(l[i].derivative(gens[i]) for i in range(n+num_poly+1))
                    reduction_dict[monomial] = temp + q
                # this could be implemented better
                # but ultimately should be done in M2 anyways
                result = term.monomial_coefficient(monomial)* reduction_dict[monomial]
                h = sum(monomial_list) + result
                monomial_list = [R(h.monomial_coefficient(monomial)) * monomial for monomial in h.monomials()]
            else:
                summer += term
        else:
            summer += term
    reduce_to_B_dict[key] = summer

def reduce_to_B(h_dict):
    ans = R(0)
    for key in h_dict:
        if key not in reduce_to_B_dict:
            compute_new_reduce_to_B_key(key)
        ans += reduce_to_B_dict[key] * h_dict[key]
    return ans


for i,B_i in enumerate(B):
    h = frobenius_on_cohom(i,prec)
    htemp = 0
    for u,g in to_ug(h):
        denom =  factorial(pole_order_vector(u)+pole_order(g))

        # denom =  factorial(pole_order_vector(u)+max_cohomology_pole_order+1)

        print(u,g)
        # this is the slow step
        u,g = reduce_griffiths_dwork(u,g)


        htemp += vector_to_monomial(u) * g / denom

    # takes h from coker (x dF/dxi) to coker (dF/dxi)
    summer = reduce_to_B(htemp.dict())

    
    for j,B_j in enumerate(B):
        frob_matrix[i][j] = summer.monomial_coefficient(B_j) * factorial(pole_order(B_j)+num_poly-1) #% p^prec
        
    print(B_i,summer)

if all(all(a.ord(p)>=0 for a in b) for b in frob_matrix):
    frob_matrix = [[a % p^ prec for a in b] for b in frob_matrix]
    print(frob_matrix)
else:
    print("Warning: non-invertible elements encountered")
frob_matrix = matrix(frob_matrix)

poly = frob_matrix.characteristic_polynomial()
print(poly)

# todo: mod these based on weil conjectures
poly = sum((poly.monomial_coefficient(mono) % p^(prec) )* mono if (poly.monomial_coefficient(mono)% p^prec) < (p^prec)//2 else (-p^(prec)+poly.monomial_coefficient(mono)% p^prec)*mono for mono in poly.monomials())
print(poly)
print(p)
print(size_pn)




