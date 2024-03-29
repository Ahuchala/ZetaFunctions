

# This file was *autogenerated* from the file init.sage
from sage.all_cmdline import *   # import sage library

_sage_const_7 = Integer(7); _sage_const_1 = Integer(1); _sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_6 = Integer(6); _sage_const_0 = Integer(0)# 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
sage.repl.load.load(sage.repl.load.base64.b64decode("YXV4X2Z1bmN0aW9ucy5zYWdl"),globals(),False)


p = _sage_const_7 


R = QQ['x, y, z']; (x, y, z,) = R._first_ngens(3)
weight = [_sage_const_1 ,_sage_const_3 ,_sage_const_1 ]

f = y**_sage_const_2  - x**_sage_const_6  - z**_sage_const_6  - x**_sage_const_3 *z**_sage_const_3 
n = len(R.gens())

d = f.degree()
fdegree = d
prec = _sage_const_1 
Rgens = R.gens()

def monomial_to_vector(m):
    return list(R(m).exponents()[_sage_const_0 ])


def vector_to_monomial(v,l=_sage_const_0 ):
    if len(v) == n-_sage_const_1 :
        return 'error'
    return R(prod([Rgens[i]**v[i] for i in range(n)]))

fweight = sum([weight[i] * monomial_to_vector(f.monomials()[_sage_const_0 ])[i] for i in range(n)])

vertices = [list(_) for _ in matrix.identity(n-_sage_const_1 )] + [(n-_sage_const_1 ) * [_sage_const_0 ]]
for i in range(n-_sage_const_1 ):
    for j in range(n-_sage_const_1 ):
        vertices[i][j] = vertices[i][j] * fweight//weight[i+_sage_const_1 ]
print(vertices)

def scale_by_d(scale):
    return [[(scale*_) for _ in __] for __ in vertices]    

# v a point in point in P_dl
def affine_vector_to_monomial(v,l):
    return R(Rgens[_sage_const_0 ]**(((l*d)-sum([v[i] * weight[i+_sage_const_1 ] for i in range(n-_sage_const_1 )]))//weight[_sage_const_0 ]) * prod([Rgens[i+_sage_const_1 ]**v[i] for i in range(n-_sage_const_1 )]))


# this only works for weighted projective space
def degree(g):
    if g == R(_sage_const_0 ):
        return -_sage_const_1 
    if g == R(_sage_const_1 ):
        return _sage_const_0 
#     shouldn't matter which one we pick
    mon = g.monomials()[_sage_const_0 ]
    return sum([weight[i] * monomial_to_vector(mon)[i] for i in range(n)]) // fweight
#     return (g.degree()//fdegree

I = R.ideal([Rgens[i]*f.derivative(Rgens[i]) for i in range(n)])

J = R.quotient(I)

df = [Rgens[i]*f.derivative(Rgens[i]) for i in range(n)]

print(df)

