loadPackage "Schubert2"

k = 3
n = 4

d = 0




G = flagBundle {k,n-k} -- secretly Gr(k,n)
for d from 1 to 25 list (
    X = sectionZeroLocus(OO_G(d)); -- quadric and a cubic
    OmG = cotangentBundle G;
    OmX = cotangentBundle X;

    ls = for i from 0 to floor(dim X / 2) list
        abs(chi exteriorPower(i, OmX) - chi exteriorPower(i, OmG));

    print (join(ls,reverse ls))#0
)