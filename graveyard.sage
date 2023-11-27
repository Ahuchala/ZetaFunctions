use_macaulay = False

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



# def Ruv_helper(u,v,g):
#     gi = lift_poly(vector_to_monomial(v)*g)
#     h = sum([(u[i] +1)*gi[i] + Rgens[i] * (gi[i]).derivative(Rgens[i]) for i in range(n)])
#     return h

def Ruv(u,v,g):
    if not tuple(v) in Ruv_u_dict.keys():
        compute_Ruv(v)
    return g * (Ruv_const_dict[tuple(v)] + sum([u[i] * Ruv_u_dict[tuple(v)][i] for i in range(n)]))

assert(ProjectiveSpace(R.change_ring(GF(p))).subscheme(f).is_smooth())


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