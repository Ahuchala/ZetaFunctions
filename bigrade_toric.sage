import itertools

# poly_list = list(var('f_%d' % i) for i in range(num_poly))

# Z = V(f0,..,fc)
# then canonical bundle of Z w_z \cong O_Z(m) with
# m = sum(di) - n - 1

# J_p,m = H^{n-p,p}_van(X)

# d_i = deg(f_i)

# bigrading is
# deg(x_i) = (0,1)
# deg(y_i) = (1,-deg(f_i))


# number of hypersurfaces in complete intersection

num_poly = 2

p = 5
prec = 2


R.<x_0,x_1,x_2,x_3,y_1,y_2> = QQ[]
# R.<x_0,x_1,x_2,x_3,y_1> = QQ[]
# R.<x_0,x_1,x_2,x_3,y_1,y_2> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,y_1,y_2> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4,x_5,y_1,y_2> = QQ[]


gens = R.gens()
prod_gens = prod(gens)
n = len(gens) - num_poly-1
# assert len(gens) == num_poly + n+1

x_vars = gens[:n+1]
y_vars = gens[n+1:]

# J_p,m consists of p copies of y_i and sum_d_i copies of x_j
# maybe compute by first all monomials in J_p,m for fixed p, then finding a basis

# f = x_0^3 + x_1^3 + x_2^3 - x_0*x_1*x_2 #+ x_3^3
f = x_0^2 + x_1^2 + x_2^2 + x_3^2
g = x_0^2 + 2*x_1^2 + 3*x_2^2 + 4*x_3^2
# h = x_0*x_1 + x_1*x_2 + x_2*x_3


# f = x_0^4+x_0^3*x_1+x_0*x_1^3+x_0^2*x_1*x_2+x_0*x_1^2*x_2+x_1^3*x_2+x_0^2*x_2^2+x_0*x_2^3+x_1*x_2^3+x_0*x_1*x_2*x_3+x_2^3*x_3+x_0^2*x_3^2+x_0*x_1*x_3^2+x_1*x_2*x_3^2+x_0*x_3^3+x_1*x_3^3+x_2*x_3^3
# f = x_0^2+2*x_0*x_1+2*x_1^2-x_0*x_2+x_1*x_2-2*x_0*x_3-x_1*x_3+2*x_2*x_3-2*x_3^2+x_0*x_4+2*x_2*x_4+2*x_3*x_4
# g = -2*x_0^3-2*x_0^2*x_1+x_0*x_1^2-x_1^3+2*x_0^2*x_2-x_0*x_1*x_2-2*x_1^2*x_2-2*x_1*x_2^2+x_0^2*x_3+2*x_0*x_1*x_3+2*x_1^2*x_3+2*x_0*x_2*x_3-x_1*x_2*x_3+x_0*x_3^2+x_2*x_3^2-x_0^2*x_4+2*x_1^2*x_4-2*x_0*x_2*x_4+x_1*x_2*x_4-2*x_2^2*x_4-2*x_2*x_3*x_4+2*x_0*x_4^2-x_2*x_4^2-2*x_4^3
# f = sum([gens[i]^2 for i in range(len(gens[:n+1]))])

# poly_list = [f]
f_0 = f; f_1 = g;
poly_list = [f_0,f_1]
# poly_list = [f,g,h]

# todo: check if no rational points?

d = [_.degree() for _ in poly_list]
m = sum([f.degree() for f in poly_list]) - n - 1

F = sum([y_vars[i] * poly_list[i] for i in range(num_poly)])


assert p in Primes(), f'Warning: {p} is not prime'


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
R = QQ[x_0..x_{n},y_1..y_{num_poly}, Degrees=>{{{n+1}:{{0,1}},
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
        eval("B.append(R(" + _ +")*prod(gens))")
print(s)
assert False

# Macaulay2 code to run to compute P1
s = f'''
R = QQ[x_0..x_{n},y_1..y_{num_poly}, Degrees=>{{{n+1}:{{0,1}},
{str([(1,-i) for i in d])[1:-1].replace('(','{').replace(')','}')}}}];
toString basis({{1,{m}}}, R)
'''
t = str(macaulay2(s))
t = t.replace("{","").replace("}","").replace("^","**").replace(" ","").split("matrix")[1:]
t = "".join(t).split(",")
P1 = []
for _ in t:
    if _ != "":
        eval("P1.append(R(" + _ +")*prod(gens))")
print(P1)


# Macaulay2 code to run to compute Pn
s = f'''
R = QQ[x_0..x_{n},y_1..y_{num_poly}, Degrees=>{{{n+1}:{{0,1}},
{str([(1,-i) for i in d])[1:-1].replace('(','{').replace(')','}')}}}];
toString basis({{{n},{m}}}, R)
'''
t = str(macaulay2(s))
t = t.replace("{","").replace("}","").replace("^","**").replace(" ","").split("matrix")[1:]
t = "".join(t).split(",")
Pn = []
for _ in t:
    if _ != "":
        eval("Pn.append(R(" + _ +")*prod(gens))")
print(Pn)
size_pn = len(Pn)

# s += "for p from 0 to " + str(n+2) + " list hilbertFunction({p,"  + str(m) + "}, J);"

# todo: go ahead and convert this into a basis sage likes


I = R.ideal([gen * F.derivative(gen) for gen in gens])
# I = F.jacobian_ideal()
J = R.quotient_ring(I)

B = sum([J(_).lift().monomials() for _ in B],[])
print(B)


def monomial_to_vector(m):
    return list(m.exponents()[0])

def vector_to_monomial(v):
    return prod([gens[i]^v[i] for i in range(n+num_poly+1)])

P1_pts = [monomial_to_vector(_) for _ in P1]
Pn_pts = [monomial_to_vector(mon) for mon in Pn]

def monomial_degree(m):
#   e looks like  [(1, 1, 1, 1, 1)]
    e = m.exponents()[0]
#   x_i contribute (0, sum(e[:n]))
#   y_i contribute (sum(e[n:],sum([-poly_list[i].degree() for i in range(num_poly)] ))
    return [sum(e[n+1:]), sum(e[:n+1])+sum([-d[i] for i in range(num_poly)] )]

# each monomial is of the form m / F^j; returns j
def pole_order(m):
    e = m.exponents()[0]
    return sum(e[n+1:])#  - num_poly

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
        numer = binomial(-d,j) * binomial(d+prec-1,d+j) * sigma_g*sigma(fj) / factorial(pole_order(sigma_g) + pole_order(sigma(fj))+num_poly-2)
        summer += numer
        fj *= F
    return summer


def frobenius_on_cohom(i,prec = 2):
    g = R(B[i])#*prod(gens)
    g = frobenius(g,prec)
    return R(g *p^(n -1) )#/ prod(gens))


reduction_dict = {}

frob_matrix = [[0 for i in range(len(B))] for j in range(len(B))]

for i in range(len(B)):
    h = frobenius_on_cohom(i,prec)
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
                    temp = sum([gens[i] * (l[i].derivative(gens[i])) for i in range(n+num_poly+1)])
                    reduction_dict[monomial] = temp + q
                result = term.monomial_coefficient(monomial)* reduction_dict[monomial]
                # result = sum([_*term.monomial_coefficient(monomial) * result.monomial_coefficient(_) for _ in result.monomials()])
                # for _ in result.monomials():
                    # monomial_list.append(_ * result.monomial_coefficient(_))
                h = sum(monomial_list) + result
                monomial_list = [R(h.monomial_coefficient(monomial)) * monomial for monomial in h.monomials()]
                
            else:
                summer += term
        else:
            summer += term
    
    for j in range(len(B)):
        frob_matrix[i][j] = summer.monomial_coefficient(R(B[j])) * factorial(pole_order(B[j])+num_poly -2) #% p^prec
        
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