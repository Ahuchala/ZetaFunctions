loadPackage "Schubert2";
k = 1
n = 5

d = 5

G = flagBundle {k,n-k} -- secretly Gr(k,n)

X = sectionZeroLocus(OO_G(d)); -- smooth hyperplane of degree d


OmG = cotangentBundle G;
OmX = cotangentBundle X;

ls = for i from 0 to floor(dim X / 2) list
    abs(chi exteriorPower(i, OmX) - chi exteriorPower(i, OmG));
print join(ls,reverse ls); -- primitive hodge numbers


Tx = tangentBundle X;
chern Tx -- Chern class of tangent bundle
chern exteriorPower(3, Tx) -- Chern class of Det(Tangent bundle)