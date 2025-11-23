loadPackage "Schubert2"

k = 2
n = 5

d = [2]
G = flagBundle {k,n-k} -- secretly Gr(k,n)

X = sectionZeroLocus(sum for i from 1 to #d list OO_G(i)); -- quadric and a cubic
OmG = cotangentBundle G;
OmX = cotangentBundle X;

