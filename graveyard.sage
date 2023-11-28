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




# break up h into chunks of monomials in connected components
# for now, just rule out the ones with no neighbors
def to_uvg_helper(h):
    hdict = h.dict()
    hdict_keys = hdict.keys()
    keys_to_pop = set()
    return_list = []

    degrees = set([degree_vector(a) for a in hdict_keys])
    for deg in degrees:
        mon_keys = set([a for a in hdict_keys if degree_vector(a) == deg])


        for h in mon_keys:
            KEEP_GOING = True # necessary?
            for h2 in mon_keys:
                
                # no other monomials could share a factor with h
                if KEEP_GOING and sum([degree_vector([min(h[i],h2[i]) for i in range(n)]) >= deg - (n- 1)]) == 1:
                    keys_to_pop.add(h)
                    
                    for v in Pn_minus_1_pts:
                        if all([h[i]>=v[i] for i in range(n)]):
                            if KEEP_GOING:
                                return_list.append([[h[i]-v[i] for i in range(n)], hdict[h]*vector_to_monomial(v)])
                            KEEP_GOING = False
                            break
                    KEEP_GOING = False

    # print(hdict)
    for _ in keys_to_pop:
        hdict.pop(_)
    # print(hdict)
    

    return hdict,return_list


# given frob(g) = h, return a minimal set of elements of form [u,g] such that h = sum(x^u * g)
def to_uvg(h):

    print("beginning to_uvg")
    hdict, return_list = to_uvg_helper(h)
    # return_list = []
    # hdict = h.dict()
    hdict_keys = hdict.keys()

    degrees = set([degree_vector(a) for a in hdict_keys])

    
    # I think this might be badly recursive but should ultimately save matrix multiplications

    # first stratify by degrees
    for deg in degrees:
        mon_keys = set([a for a in hdict_keys if degree_vector(a) == deg])

        # u_checked = set()
        while len(mon_keys)>0:
            # search through all valid g options to find all valid u options
            # might be good to remember which values of u we've checked to avoid duplicates

            div_count_record = -1
            u_record = 0
            # u_checked = set()

            for v in Pn_minus_1_pts:


                for mon in mon_keys:
                    if all([mon[i]>=v[i] for i in range(n)]):
                        u = [mon[i]-v[i] for i in range(n)]
                        # if not tuple(u) in u_checked:
                        
                        # now check how many elements in mon_keys are divisible by u
                        div_count = len([0 for a in mon_keys if all([a[i]>=u[i] for i in range(n)]) and degree_vector([a[i]-u[i] for i in range(n)])>0])
                        if div_count > div_count_record:
                            div_count_record = div_count
                            u_record = u
                        # else:
                            # u_checked.add(tuple(u))

            # now append the tuple [u,g] and remove the corresponding monomials from mon_keys
            mons = [a for a in mon_keys if all([a[i]>=u_record[i] for i in range(n)])]
            mon_keys -= set(mons)
            g = sum([hdict[a] * vector_to_monomial([a[i]-u_record[i] for i in range(n)]) for a in mons])
            return_list.append([u_record,g])



    print("end to_uvg")
    print(return_list)
    return return_list



def reduce_griffiths_dwork(u,g):
    g_vec = vector(matrix(to_pn_minus_1_basis(g)))
    # todo: speed up!
    for v in P1_pts:
        if (u not in P1_pts) and all([u[i]>=v[i] for i in range(n)]):

            if not tuple(v) in Ruv_u_dict.keys():
                compute_Ruv(v)
            Ruv_const_dict_tuple_v = Ruv_const_dict[tuple(v)]
            Ruv_u_dict_tuple_v = Ruv_u_dict[tuple(v)]


            while (u not in P1_pts) and all([u[i]>=v[i] for i in range(n)]):
                u = [u[i]-v[i] for i in range(n)]
                # left g -> g A
                g_vec *=  (Ruv_const_dict_tuple_v + sum([u[i] * Ruv_u_dict_tuple_v[i] for i in range(n)]))
    g = from_pn_minus_1_basis(vector(g_vec))
    
    return u,g