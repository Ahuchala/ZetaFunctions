from sympy.utilities.iterables import multiset_permutations
from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL

p = 7
# p = 61

# R.<x,y,z> = QQ[]

weight = [1,1,1,1]
# weight = [1,1,3,1]
# weight = [1,3,1]
# weight = [1,1,1]

R.<w,x,y,z> = QQ[]
# f = w^4 + x^4 + y^4 + z^4
f = w^4+2*w*x^3-2*x^4-x^3*y-x^2*y^2-y^4+w^3*z-x^3*z-2*w^2*y*z+2*w*x*y*z-x^2*y*z-w*y^2*z+2*x*y^2*z-2*y^3*z-w^2*z^2-2*w*x*z^2+x^2*z^2-2*w*y*z^2+x*y*z^2+y^2*z^2+2*w*z^3+2*x*z^3-2*y*z^3-2*z^4
# f = y^2 - x^6 - z^6 - x^3*z^3
# f = y^2 - x*w*z*(x+w+z)*(x+2*w+z)*(3*x+2*w-z)
# f =x^6 - z^2 + y^6 + w^6
# f = w^6 - x^6 - y^6 - z^2
n = len(R.gens())

# f = x^5 + y^5 + z^5
# f = x^4 + y^4 + z^4

# f = x*y*z + 2*x^2 * y + 3*x^2 * 4*z+ x*y^2 + 5*y^2 * z + 6*x^3 + 7*y^3 + 8*z^3
# f = x^3  + y^3 + z^3-x*y*z
d = f.degree()
fdegree = d

prec = 2
Rgens = R.gens()

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
print('Computing primitive cohomology')

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
    # return R(g * p^(n-2))
    return R(g * p)


# reduction_dict = {}

frob_matrix = [[0 for i in range(len(B))] for j in range(len(B))]


print('computing frobenius action and reduction')


for i in range(len(B)):
    g = frobenius_on_cohom(i,prec)
    print(B[i])
    ans = R(0)
    while g != R(0):
#         print(g)
        q = J(g).lift()
        r = (g-q).lift(I)
        ans += q
        if all([_ in B for _ in g.monomials()]) or degree(g) == 1:
#             print(ans) #print(g/fdegree)

            break
        right_hand_term = sum([ Rgens[j] * ((r[j]).derivative(Rgens[j])) for j in range(n)])
        right_hand_term = sum([mono * right_hand_term.monomial_coefficient(mono) / (degree(mono)) for mono in right_hand_term.monomials()])

#         right_hand_term = sum([mono * right_hand_term.monomial_coefficient(mono) / (mono.degree()//fdegree) for mono in right_hand_term.monomials()])
        g = right_hand_term
    for j in range(len(B)):
        frob_matrix[i][j] = ans.monomial_coefficient(R(B[j])) % p^prec
        
frob_matrix = matrix(frob_matrix)
print(frob_matrix)
print(frob_matrix.characteristic_polynomial() %p^(prec))

poly = R(frob_matrix.characteristic_polynomial())
print(poly)
poly = sum([(poly.monomial_coefficient(mono) % p^(prec) )* mono if (poly.monomial_coefficient(mono)% p^prec) < (p^prec)//2 else (-p^(prec)+poly.monomial_coefficient(mono)% p^prec)*mono for mono in poly.monomials()])
print(poly)

for r in range(1,len(B)+1):
    print(sum([p^(_*r) for _ in range(3)]) + (-1)^n * (frob_matrix^r).trace())