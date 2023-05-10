from sympy.utilities.iterables import multiset_permutations
from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL

p = 29
# p = 61
# p = 389

R.<x,y,z> = QQ[]

# weight = [1,1,1,1]
# weight = [1,1,3,1]
# weight = [1,3,1]
weight = [1,1,1]

# R.<w,x,y,z> = QQ[]
# f = w^4 + x^4 + y^4 + z^4
# f = w^4+2*w*x^3-2*x^4-x^3*y-x^2*y^2-y^4+w^3*z-x^3*z-2*w^2*y*z+2*w*x*y*z-x^2*y*z-w*y^2*z+2*x*y^2*z-2*y^3*z-w^2*z^2-2*w*x*z^2+x^2*z^2-2*w*y*z^2+x*y*z^2+y^2*z^2+2*w*z^3+2*x*z^3-2*y*z^3-2*z^4
# f = y^2 - x^6 - z^6 - x^3*z^3
# f = y^2 - x*w*z*(x+w+z)*(x+2*w+z)*(3*x+2*w-z)
# f =x^6 - z^2 + y^6 + w^6
# f = w^6 - x^6 - y^6 - z^2
# f = x^6 - y^2 + z^6
n = len(R.gens())
f = x^3 + y^3 + z^3

# f = x^5 + y^5 + z^5
# f = x^4 + y^4 + z^4

# f = x*y*z + 2*x^2 * y + 3*x^2 * 4*z+ x*y^2 + 5*y^2 * z + 6*x^3 + 7*y^3 + 8*z^3
# f = x^3  + y^3 + z^3-x*y*z
d = f.degree()
fdegree = d

prec = 2
Rgens = R.gens()

# eg vertices = [[1,0],[0,1],[0,0]]
# vertices = f.newton_polytope().vertices()
# vertices = [[2,0],[0,6],[0,0]]
# vertices = [[6,0,0],[0,2,0],[0,0,0],[0,0,6]]

# vertices = [[2,0,0],[0,6,0],[0,0,0],[0,0,6]]

# vertices = scale_by_d(d)
def monomial_to_vector(m):
    return list(R(m).exponents()[0])


def vector_to_monomial(v,l=0):
    if len(v) == n-1:
            return 'error'
    return R(prod([Rgens[i]^v[i] for i in range(n)]))

fweight = sum([weight[i] * monomial_to_vector(f.monomials()[0])[i] for i in range(n)])

vertices = [list(_) for _ in matrix.identity(n-1)] + [(n-1) * [0]]
for i in range(n-1):
    for j in range(n-1):
        vertices[i][j] = vertices[i][j] * fweight//weight[i+1]
print(vertices)

def scale_by_d(scale):
    return [[(scale*_) for _ in __] for __ in vertices]    

# v a point in point in P_dl
def affine_vector_to_monomial(v,l):
    return R(Rgens[0]^(((l*d)-sum([v[i] * weight[i+1] for i in range(n-1)]))//weight[0]) * prod([Rgens[i+1]^v[i] for i in range(n-1)]))


# this only works for weighted projective space
def degree(g):
    if g == R(0):
        return -1
    if g == R(1):
        return 0
#     shouldn't matter which one we pick
    mon = g.monomials()[0]
    return sum([weight[i] * monomial_to_vector(mon)[i] for i in range(n)]) // fweight
#     return (g.degree()//fdegree

I = R.ideal([Rgens[i]*f.derivative(Rgens[i]) for i in range(n)])

J = R.quotient(I)

df = [Rgens[i]*f.derivative(Rgens[i]) for i in range(n)]


lift_dict = {}

# TODO: cache lifts

def compute_primitive_cohomology():
    return_set = []
    for scale in range(1,n):
        P_int = list(LatticePolytope(scale_by_d(scale)).interior_points())
        
        monomial_int = [I.reduce(affine_vector_to_monomial(_,scale)) for _ in P_int]
        
#       use gram schmidt to find a basis of (P_int + If)/(If)
        s = Sequence(monomial_int)
        C, m = s.coefficient_matrix()
        G,M = C.gram_schmidt()
        G = G.echelon_form()
        for _ in G:
#           go from vector basis back to monomials
            a = vector(m).dot_product(_)
#           clean up constants
            return_set.append( a / gcd(a.monomial_coefficient(mono) for mono in a.monomials()))
            
    return (list(return_set))
B = compute_primitive_cohomology()

# def compute_primitive_cohomology():
#     return_set = []
#     for scale in range(1,n):
#         l = list(LatticePolytope(scale_by_d(d*scale)).interior_points())
        
#         temp_return_set = []
#         l = [I.reduce(affine_vector_to_monomial(_,scale)) for _ in l]
#         for _ in l:
#             ls = temp_return_set + [_]
#             s = Sequence(ls)
# #             print(s)
#             C, m = s.coefficient_matrix()
#             if C.kernel().dimension() == 0:
#                 temp_return_set.append(_)
#                 return_set.append(_)
#     return (list(return_set))
# B = compute_primitive_cohomology()


print(len(B), (-1)^(n)*(1/d * ((1-d)^(n)-1)+1))

def sigma(g):
    return sum([g.coefficient(monomial)* monomial^p for monomial in g.monomials()])

def frobenius(g,prec=2):
    d = degree(g)
#     d = g.degree() // fdegree
    summer = 0
    sigma_g = sigma(g)
    fj = R(1)
    for j in range(prec):
        numer = binomial(-d,j) * binomial(d+prec-1,d+j) * sigma_g*sigma(fj)
        summer += numer
        fj *= f
    return summer

def frobenius_on_cohom(i,prec = 2):
    g = B[i]
    g = frobenius(g,prec)
    return R(g * p^(n-2))
#     return R(g * p)


# reduction_dict = {}

# frob_matrix = [[0 for i in range(len(B))] for j in range(len(B))]


# print('computing frobenius action and reduction')


# for i in range(len(B)):
#     g = frobenius_on_cohom(i,prec)
#     print(B[i])
#     ans = R(0)
#     while g != R(0):
# #         print(g)
#         q = J(g).lift()
#         r = (g-q).lift(I)
#         ans += q
#         if all([_ in B for _ in g.monomials()]) or degree(g) == 1:
#             # print(ans) #print(g/fdegree)

#             break
#         right_hand_term = sum([ Rgens[j] * ((r[j]).derivative(Rgens[j])) for j in range(n)])
#         right_hand_term = sum([mono * right_hand_term.monomial_coefficient(mono) / (degree(mono)) for mono in right_hand_term.monomials()])

# #         right_hand_term = sum([mono * right_hand_term.monomial_coefficient(mono) / (mono.degree()//fdegree) for mono in right_hand_term.monomials()])
#         g = right_hand_term
#     for j in range(len(B)):
#         frob_matrix[i][j] = ans.monomial_coefficient(R(B[j])) % p^prec
        
# frob_matrix = matrix(frob_matrix)
# print(frob_matrix)
# print(frob_matrix.characteristic_polynomial())
# print(frob_matrix.characteristic_polynomial() % p^prec)



qr_dict = {}

for scale in range(n+1):
#         P_int = list(LatticePolytope(scale_by_d(d*scale)).points())

    P_int = list(LatticePolytope(scale_by_d(scale)).points())

#     P_int = list(LatticePolytope(scale_by_d(scale)).interior_points())
#     print(P_int)
    monomial_int = [affine_vector_to_monomial(_,scale) for _ in P_int]
    for monomial in monomial_int:
        q = J(monomial).lift()
        r = (monomial - q).lift(I)
        qr_dict[monomial] = [q,r]

R_dict = {}
# dear god make this not recursive

def reduce_monomial(u,h):
    if I.reduce(h) == h:
        return h
    q,r = qr_dict[h]
#     l = degree(h)+degree(vector_to_monomial(u))
    
#     ans = q+sum([(u[i] + 1)*r[i] + Rgens[i] * (r[i]).derivative(Rgens[i]) for i in range(n)])
# WHY THE HELL DID HE WRITE IT AS u[i] + 1?????
#     ans = q+(sum([u[i]*r[i] + Rgens[i] * (r[i]).derivative(Rgens[i]) for i in range(n)]))
    ans = q+sum([Rgens[i] * (r[i]).derivative(Rgens[i]) for i in range(n)])
    ans = sum([ans.monomial_coefficient(monomial) * reduce_monomial(u,monomial) for monomial in ans.monomials()])
    for i in range(n):
        for monomial in r[i].monomials():
            ans += u[i] * r[i].monomial_coefficient(monomial)*reduce_monomial(u,monomial)
    return ans
#     return sum([ans.monomial_coefficient(monomial)*reduce_monomial(u,monomial) for monomial in ans.monomials()])

#     else:
#         dic = h.dict()
#         return sum([prod([u[i] ^ key[i] for i in range(n)]) * reduce_monomial(u,dic[key]) for key in dic.keys()])

    
# for scale in range(1,n):

#     P_int = list(LatticePolytope(scale_by_d(scale)).interior_points())

#     monomial_int = [affine_vector_to_monomial(_,scale) for _ in P_int]
#     for monomial in monomial_int:
#         q,r = qr_dict[monomial]
#         print(q,r)

# var('a b c d')
return_set = []
P1 = [affine_vector_to_monomial(_,1) for _ in LatticePolytope(scale_by_d(1)).points()]
P1 = sorted(P1)
for scale in range(n):
    P_int = list(LatticePolytope(scale_by_d(scale)).points())
    monomial_int = []
    for _ in P_int:
        q,r = qr_dict[affine_vector_to_monomial(_,scale)]
        monomial_int.append(q)
#     monomial_int = [I.reduce(affine_vector_to_monomial(_,scale)) for _ in P_int]

#     should be able to get a monomial basis by adding one monomial at a time

#       use gram schmidt to find a basis of (P_int + If)/(If)
    s = Sequence(monomial_int)
    C, m = s.coefficient_matrix()
    G,M = C.gram_schmidt()
    G = G.echelon_form()
    for _ in G:
#           go from vector basis back to monomials
        a = vector(m).dot_product(_)
#           clean up constants
        return_set.append( a / gcd(a.monomial_coefficient(mono) for mono in a.monomials()))

psi_basis = (list(return_set))
print(psi_basis)
print(len(psi_basis))

sv_dict = {}

S.<u1,u2,u3> = R[]
U = [u1,u2,u3]
if n == 4:
    S.<u1,u2,u3,u4> = R[]
    U = [u1,u2,u3,u4]


for mon in psi_basis:
    for v in P1:
        sv_dict[(mon,v)] = reduce_monomial(U, v*mon)

# assemble info from sv_dict into a dict of matrices


s_dict = {}
len_psi_basis = len(psi_basis)

# this is slow so only compute it for v in P1 when needed
def compute_s_dict_v(v):
    mat = [[S(0) for _ in range(len_psi_basis)] for __ in range(len_psi_basis)]
    for i in range(len_psi_basis):
        row_vec = sv_dict[(psi_basis[i],v)]
        if row_vec in R:
            for monomial in row_vec.monomials():
                mat[i][psi_basis.index(monomial)] = mat[i][psi_basis.index(monomial)] + row_vec.monomial_coefficient(monomial)
        else:
            dic = row_vec.dict()
            for term in dic.keys():
                for monomial in R(dic[term]).monomials():
                    mat[i][psi_basis.index(monomial)] = mat[i][psi_basis.index(monomial)] + R(dic[term]).monomial_coefficient(monomial) * prod([U[_]^term[_] for _ in range(n)])
    s_dict[v] = matrix(mat)

# does this do anything

# R.<x,y,z> = QQ[]
# if n == 4:
#     R.<w,x,y,z> = QQ[]

def reduce_uv(monomial,coeff):
    u = monomial_to_vector(monomial)
    g = [R(0)] * len_psi_basis
    g[0] = coeff
    g = matrix(g)

    for v in P1:
        v_vec = monomial_to_vector(v)
        if not v in s_dict.keys():
            compute_s_dict_v(v)
        mat = (s_dict[v])

        while all([u[i]>=v_vec[i] for i in range(n)]):

            for ii in range(n):
                u[ii] -= v_vec[ii]
            if n == 3:
                # print(u)
                m = mat.subs(u1=u[0],u2=u[1],u3=u[2])
                # print(m)
            elif n == 4:
                m = mat.subs(u1=u[0],u2=u[1],u3=u[2],u4=u[3])
            g = g*m# not transpose()

    h = (sum([g[0][ii] * psi_basis[ii] for ii in range(len_psi_basis)]))
    return h

frob_matrix = [[0 for i in range(len(B))] for j in range(len(B))]


for i in range(len(B)):
    ans =R(0)
    fro = frobenius_on_cohom(i,prec)
    print(B[i])
    for monomial in fro.monomials():
        h = reduce_uv(monomial,fro.monomial_coefficient(monomial)) / factorial(degree(monomial)-1)
        print(h)
        for b in B:
            ans += b * (h.monomial_coefficient(b)%p^prec )
        for j in range(len(B)):
            frob_matrix[i][j] = ans.monomial_coefficient(R(B[j])) % p^prec
        
frob_matrix = matrix(frob_matrix)
print(frob_matrix)
print(frob_matrix.characteristic_polynomial())
print(frob_matrix.characteristic_polynomial() % p^prec)