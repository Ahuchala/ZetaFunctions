from sympy.utilities.iterables import multiset_permutations


use_macaulay = False

n = 4 #number of variables

p = 17

prec =3

# K = Qp(p,8)

R.<x,y,z> = QQ[]
# R.<w,x,y,z> = QQ[]
# R.<x_0,x_1,x_2,x_3,x_4> = QQ[]

Rgens = R.gens()



# f = x^3 + y^3 + z^3 - x*y*z
# f = x^4 + y^4 + z^4

f = x^4 + y^4 + z^4-x*y*z^2+4*x^2*z^2

# f = (2)*x^3+3*x^2*y+(4)*x*y^2+(5)*y^3+(5)*x^2*z+x*y*z+(7)*y^2*z+(5)*x*z^2+(7)*y*z^2+(1)*z^3
# f = x^6 + y^6 + z^6 - x^4*y*z
# f = x^9 + y^9 + z^9
# f = w^4 + x^4 + y^4  + z^4
# f = x^2 + x*y + y^2 + z^2
# f= 2*x_0^3+2*x_0*x_1^2+x_1^3+2*x_0^2*x_2-x_0*x_1*x_2+2*x_1^2*x_2+x_0*x_2^2+2*x_1*x_2^2+x_2^3-x_0^2*x_3-x_0*x_1*x_3-x_0*x_2*x_3+x_1*x_2*x_3+2*x_2^2*x_3+x_0*x_3^2-x_1*x_3^2-x_2*x_3^2+2*x_0^2*x_4+2*x_0*x_1*x_4-x_1^2*x_4-2*x_0*x_2*x_4-x_1*x_2*x_4+2*x_2^2*x_4+x_0*x_4^2-x_1*x_4^2-2*x_2*x_4^2
# f = w^4+2*w*x^3-2*x^4-x^3*y-x^2*y^2-y^4+w^3*z-x^3*z-2*w^2*y*z+2*w*x*y*z-x^2*y*z-w*y^2*z+2*x*y^2*z-2*y^3*z-w^2*z^2-2*w*x*z^2+x^2*z^2-2*w*y*z^2+x*y*z^2+y^2*z^2+2*w*z^3+2*x*z^3-2*y*z^3-2*z^4
I = R.ideal([f.derivative(_) for _ in R.gens()])# + [f])
J = R.quotient(I)

xI = R.ideal([_*f.derivative(_) for _ in R.gens()])# + [f])
xJ = R.quotient(xI)

d = f.degree()
fdegree = d

# i.e. y^3 -> [0,3,0,0] when n=4
def monomial_to_vector(m):
    return list(m.exponents()[0])

def vector_to_monomial(v):
    return prod([Rgens[i]^v[i] for i in range(n)])

lift_dict = {}
for l in range(n):
    if (l*d-n)>=0:
        print('computing lifts of degree ' + str(l*d-n))
        for part in Partitions(l*d-n, max_length=n):
    #       pad with zeros
            part += [0] * (n - len(part))

            for perm in multiset_permutations(part):
                m = vector_to_monomial(perm)
    #                 use cached lifts
                c = J(m).lift()
                r = (m-c).lift(I)
                lift_dict[m] = (c,r)

if use_macaulay:
    print('computing lifts of degree ' + str((n-1)*d-(n-1)+d) + ' using macaulay')
    s = ''
    # toString (basis(''' + str((n-1)*d-(n-1)+d) + ''',R) % M), toString (basis(''' + str((n-1)*d-(n-1)+d) + ''',R) // M)
    if n == 3:
        s = macaulay2('''
        R = QQ[x..z]; 
        f = ''' + str(f) + ''';
        M = matrix{ for v in gens R list v*diff(v,f) };
        toString basis(''' + str((n-1)*d-(n-1)+d) + ''',R), toString (basis(''' + str((n-1)*d-(n-1)+d) + ''',R) // M)
        ''')
    elif n == 4:
        s = macaulay2('''
        R = QQ[w..z]; 
        f = ''' + str(f) + ''';
        M = matrix{for v in gens R list v*diff(v,f) };
        toString basis(''' + str((n-1)*d-(n-1)+d) + ''',R), toString (basis(''' + str((n-1)*d-(n-1)+d) + ''',R) // M)
        ''')
    # s = [str(_).replace("matrix {{","").replace("}","") for _ in s]
    # s = ",".join(s)
    s = str(s)
    s = s[1:-1]
    s = s.replace('matrix {{','').replace(' ','')
    b,s = s.split("}}")[:-1]
    b = b.split(',')
    # s.split('},{')
    # s = [_.split(',') for _ in s]
    s = s.split('},{')
    s = [_.split(',') for _ in s]
    s[0] = s[0][1:]
    for i in range(len(b)):
        lift_dict[R(b[i])] = [R(_[i]) for _ in s]
        # print(b[i],[_[i] for _ in s])
    # print(s)

else:
    print('computing lifts of degree ' + str((n-1)*d-(n-1)+d) + ' using sage')
    for part in Partitions((n-1)*d-(n-1)+d, max_length=n):
    #       pad with zeros
        part += [0] * (n - len(part))

        for perm in multiset_permutations(part):
            m = vector_to_monomial(perm)
    #                 use cached lifts
            c = xJ(m).lift()
            r = (m-c).lift(xI)
            lift_dict[m] = r
            


def compute_primitive_cohomology():

#     R = QQ[w..z]
#     f = ...
#     J = R/ideal jacobian f
#     toString basis(4,J)
    B = set()
    for l in range(n):
#         print('computing Jacobian basis of rank ' + str(l*d))
        if (l*d-n)>=0:
            for part in Partitions(l*d-n, max_length=n):
        #       pad with zeros
                part += [0] * (n - len(part))

                for perm in multiset_permutations(part):
                    m = vector_to_monomial(perm)
                    for monomial in J(m).lift():
                        B.add(monomial[1])
    #                 use cached lifts
#                     c,r = lift_dict[m]
#                     for monomial in c.monomials():
#                         B.add(monomial)
    return list(B)
B = compute_primitive_cohomology()

print(len(B), (-1)^(n)*(1/d * ((1-d)^(n)-1)+1))


def sigma(g):
    return sum([g.coefficient(monomial)* monomial^p for monomial in g.monomials()])

def frobenius(g,prec=2):
#     m = g.degree()//f.degree()
    d = g.degree() // f.degree()
    summer = 0
    sigma_g = sigma(g)
    fj = R(1)
    for j in range(prec):
#         summer += binomial(-d,j) * binomial(d+N-1,d+j) * sigma(g)*sigma(qf^j)
#         cacheing may be appropriate
        numer = binomial(-d,j) * binomial(d+prec-1,d+j) * sigma_g*sigma(fj)
        summer += numer #/(f^(p*(d+j)))
        fj *= f
    return summer

def frobenius_on_cohom(i,prec = 2):
    g = B[i]*prod(Rgens)
    g = frobenius(g,prec)
    return R(g *p^(n -2) / prod(Rgens))

def lift_poly(g):
    summer = n * [0]
    for monomial in g.monomials():
        monomial_lift = lift_dict[monomial]
        c = g.monomial_coefficient(monomial)
        for i in range(n):
            summer[i] += c*monomial_lift[i]
    return summer

def Ruv(u,v,g):
    gi = lift_poly(vector_to_monomial(v)*g)
    deg_g = (n-1)*fdegree-(n-1)
    m = (sum(u) +deg_g+ n) // fdegree

#     m = (sum(u) + sum(v) - g.degree() + n-1) // fdegree - 1
#     print(gi)
    h = sum([(u[i] +1)*gi[i] + Rgens[i] * (gi[i]).derivative(Rgens[i]) for i in range(n)])
    v = n * [0]
    for i in range(n):
        while u[i] > 0:
            u[i] = u[i] - 1
            v[i] = v[i] + 1
            if sum(v) == fdegree:
                return (u,v,h/m)

# this is hacky but only done once per basis element
def to_uvg(h):
    hdict = h.dict()
    return_list = []
    for etuple in hdict.keys():
        vector = list(etuple)
        c = hdict[etuple]
        
#         first check if it's divisible by something in the list
        already_divisible = False
        for j in range(len(return_list)):
            if not already_divisible:
                known_mono_tuple = return_list[j]
                u,v,g = known_mono_tuple
                if sum(u) + fdegree + g.degree() == sum(etuple):
                    if sum([etuple[i] >= u[i] + v[i] for i in range(n)]) == n:
                        g += c * vector_to_monomial([etuple[i] - u[i] - v[i] for i in range(n)])
                        return_list[j] = [u,v,g]
                        already_divisible = True
        if not already_divisible:
    
            
            d = (sum(vector) -n)// f.degree()

            u = n * [0]
            v = n * [0]

            for i in range(n):
                while vector[i] > 0:
                    vector[i] = vector[i] -1
                    v[i] = v[i] + 1
                    if sum(v) == fdegree:
                        u = vector
                        vector = n * [0]

            g_vec = n * [0]
            vector = u
            for i in range(n):
                while vector[i] > 0:
                    vector[i] = vector[i] - 1
                    g_vec[i] = g_vec[i] + 1
                    if sum(g_vec) == (n-1) * fdegree - (n -1):
                        u = vector
                        vector = n * [0]
            g = c * vector_to_monomial(g_vec)
            return_list.append([u,v,g])
    return return_list


reduction_dict = {}

frob_matrix = [[0 for i in range(len(B))] for j in range(len(B))]

# sum([Qx(h).lift(J)[i] * [m  * qf.derivative(xi) for xi in qxgens][i] for i in range(3)]) == Qx(h)
for i in range(len(B)):
    h = frobenius_on_cohom(i,prec)
    print(h)
    htemp = 0
    for u,v,g in to_uvg(h):
        while sum(u) > fdegree:
            u,v,g = Ruv(u,v,g)
        htemp += vector_to_monomial(u) * vector_to_monomial(v) * g
    h = htemp
    
    summer = R(0)
#   monomials to reduce (with coefficients)
    monomial_list = [R(h.monomial_coefficient(monomial)) * monomial for monomial in h.monomials()]
    while len(monomial_list) > 0:
        term  = monomial_list.pop()
#         print(term,term.monomials()[0])
        if term not in QQ:
            if len(term.monomials())>1:
                print('error: too many terms')
            monomial = R(term.monomials()[0])
            if not monomial in B:
                if not monomial in reduction_dict.keys():
                    q = J(monomial).lift()
                    r = monomial - q
                    l = r.lift(I)
                    m = (monomial.degree() +n)// f.degree()
#                     print(monomial,m)
                    temp = sum([l[i].derivative(Rgens[i]) for i in range(n)])/(m-1) #m-1? m+1? coherent for m+1?
                    reduction_dict[monomial] = temp + q
                result = reduction_dict[monomial]
                for _ in result.monomials():
                    monomial_list.append(_*term.monomial_coefficient(monomial) * result.monomial_coefficient(_))
            else:
                summer += term.monomial_coefficient(monomial) * monomial
        else:
            summer += term
    
    for j in range(len(B)):
        frob_matrix[i][j] = summer.monomial_coefficient(R(B[j])) % p^prec
        
    print(B[i],summer)

frob_matrix = matrix(frob_matrix)
print(frob_matrix)
print(frob_matrix.characteristic_polynomial() %p^(prec))

poly = R(frob_matrix.characteristic_polynomial())
print(poly)
poly = sum([(poly.monomial_coefficient(mono) % p^(prec) )* mono if (poly.monomial_coefficient(mono)% p^prec) < (p^prec)//2 else (-p^(prec)+poly.monomial_coefficient(mono)% p^prec)*mono for mono in poly.monomials()])
print(poly)

for r in range(1,len(B)+1):
    print(sum([p^(_*r) for _ in range(3)]) + (-1)^n * (frob_matrix^r).trace())