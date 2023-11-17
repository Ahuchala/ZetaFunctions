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
