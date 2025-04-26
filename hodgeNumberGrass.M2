loadPackage "Schubert2"

k = 2
n = 5

d = [2,3]




G = flagBundle {k,n-k} -- secretly Gr(k,n)
-- for d from 1 to 25 list (
X = sectionZeroLocus(sum for i from 1 to #d list OO_G(i)); -- quadric and a cubic
OmG = cotangentBundle G;
OmX = cotangentBundle X;

ls = for i from 0 to floor(dim X / 2) list
    abs(chi exteriorPower(i, OmX) - chi exteriorPower(i, OmG));

if (n % 2 ==1) then (
    print join(ls,reverse ls);
) else (
-- include the second half, but don't repeat middle element
print (join(ls, for i from 1 to #ls-1  list ls#(#ls-i-1)));
);