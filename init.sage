load aux_functions.sage

p = 7


R.<x,y,z> = QQ[]
weight = [1,3,1]

f = y^2 - x^6 - z^6 - x^3*z^3
n = len(R.gens())

d = f.degree()
fdegree = d
prec = 1
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

print(df)