needsPackage "Resultants";
needsPackage "TensorComplexes"; -- for multiSubsets function
needsPackage "Divisor";
needsPackage "Schubert2"; -- for Hodge number verification
allowableThreads = 4;
n = 10;
k = 3;

-- K = QQ;
K = ZZ/101;
d = 1;
permList = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {1, 3, 0, 6, 9, 8, 2, 5, 4, 7}, {2, 0, 6, 1, 8, 7, 3, 9, 5, 4}, {3, 6, 1, 2, 7, 4, 0, 8, 9, 5}, {4, 8, 9, 5, 0, 3, 7, 6, 1, 2}, {5, 7, 8, 9, 6, 0, 4, 1, 2, 3}, {6, 2, 3, 0, 5, 9, 1, 4, 7, 8}, {7, 9, 5, 4, 3, 2, 8, 0, 6, 1}, {8, 5, 4, 7, 2, 1, 9, 3, 0, 6}, {9, 4, 7, 8, 1, 6, 5, 2, 3, 0}};
-- permList = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {1, 8, 2, 6, 9, 3, 5, 4, 0, 7}, {3, 6, 2, 7, 1, 9, 4, 0, 5, 8}, {4, 9, 2, 1, 5, 0, 8, 6, 7, 3}, {5, 3, 2, 9, 0, 4, 7, 8, 6, 1}, {6, 5, 2, 4, 8, 7, 9, 1, 3, 0}, {7, 4, 2, 0, 6, 8, 1, 3, 9, 5}, {8, 0, 2, 5, 7, 6, 3, 9, 1, 4}, {9, 7, 2, 8, 3, 1, 0, 5, 4, 6}};

-- permList = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {1, 3, 0, 6, 9, 8, 2, 5, 4, 7}, {2, 0, 6, 1, 8, 7, 3, 9, 5, 4}, {3, 6, 1, 2, 7, 4, 0, 8, 9, 5}, {4, 8, 9, 5, 0, 3, 7, 6, 1, 2}, {5, 7, 8, 9, 6, 0, 4, 1, 2, 3}, {6, 2, 3, 0, 5, 9, 1, 4, 7, 8}, {7, 9, 5, 4, 3, 2, 8, 0, 6, 1}, {8, 5, 4, 7, 2, 1, 9, 3, 0, 6}, {9, 4, 7, 8, 1, 6, 5, 2, 3, 0}};


ANTISYMMETRIZE = false;
-- access generators like p_(0,2)

-- R = K[p_(0,0) .. p_(n-1,n-1)]

pluckerIdeal = Grassmannian(k-1,n-1,CoefficientRing=>K);
R = ring pluckerIdeal;




gensR = gens R;
numGens = #gensR;


countInversions = (perm) -> (
    lenPerm = #perm;
    ans = 0;
    for i from 0 to lenPerm-1 do (
  for j from i+1 to lenPerm-1 do (
if perm#j < perm#i then (
    ans += 1;
);
  );
    );
    return ans;
);

-- (1,2,3) -> false, (2,2,3) -> true
hasDuplicates = (ls) -> (
    return unique ls != ls;
);

-- define all p_I, not just for I increasing (i.e. the k wedges of x_i)
for inds in multiSubsets(n,k) do (
	-- set all p_ii to zero
	if (hasDuplicates(inds)) then (
  for perm in permutations(inds) do (
p_(toSequence perm) = 0;
  );
    ) else for perm in permutations(inds) do (
		-- e.g. set p_(i,j,k) = -p_(j,i,k)
  p_(toSequence perm) = (-1)^(countInversions(perm)) * p_(toSequence inds);
    );
);


-- f = sum(apply(gens R, i->random(0,100)*(i)^d))

randomIntPoly = (d, R) -> 
(
	g := random(d, R);
	monoms := flatten entries monomials g;
	sum apply(monoms, m -> random(100) * m)
);

-- Usage:
-- f = randomIntPoly(d, R);
-- f = -231*p_(0,2,6)  -231*p_(4,5,8) +231*p_(0,1,2)  +231*p_(0,1,3)  +231*p_(1,3,6)  +231*p_(2,3,6)  +231*p_(4,7,9)  +231*p_(4,8,9)  +231*p_(5,7,8)  +231*p_(5,7,9)
-- f = 230*p_(0,1,2)+231*p_(0,1,3)-59*p_(0,2,3)-59*p_(1,2,3)+148*p_(0,1,4)-16*p_(0,2,4)+150*p_(1,2,4)-204*p_(0,3,4)+184*p_(1,3,4)+168*p_(2,3,4)+199*p_(0,1,5)-184*p_(0,2,5)+148*p_(1,2,5)-168*p_(0,3,5)-40*p_(1,3,5)+150*p_(2,3,5)-204*p_(0,4,5)-148*p_(1,4,5)-150*p_(2,4,5)-168*p_(3,4,5)+59*p_(0,1,6)-231*p_(0,2,6)-59*p_(1,2,6)+59*p_(0,3,6)+231*p_(1,3,6)+231*p_(2,3,6)+57*p_(0,4,6)+148*p_(1,4,6)-40*p_(2,4,6)-199*p_(3,4,6)-204*p_(0,5,6)-57*p_(1,5,6)+148*p_(2,5,6)-16*p_(3,5,6)+57*p_(4,5,6)-40*p_(0,1,7)-199*p_(0,2,7)+168*p_(1,2,7)+57*p_(0,3,7)+16*p_(1,3,7)+204*p_(2,3,7)-57*p_(0,4,7)+150*p_(1,4,7)+168*p_(2,4,7)+204*p_(3,4,7)+199*p_(0,5,7)-40*p_(1,5,7)+184*p_(2,5,7)+16*p_(3,5,7)+59*p_(4,5,7)+148*p_(0,6,7)-150*p_(1,6,7)-184*p_(2,6,7)+148*p_(3,6,7)-148*p_(4,6,7)-148*p_(5,6,7)+184*p_(0,1,8)-148*p_(0,2,8)+204*p_(1,2,8)-148*p_(0,3,8)+199*p_(1,3,8)-57*p_(2,3,8)+148*p_(0,4,8)+184*p_(1,4,8)+16*p_(2,4,8)+199*p_(3,4,8)-184*p_(0,5,8)-199*p_(1,5,8)-148*p_(2,5,8)+40*p_(3,5,8)-231*p_(4,5,8)+150*p_(0,6,8)-168*p_(1,6,8)-16*p_(2,6,8)-40*p_(3,6,8)+40*p_(4,6,8)+16*p_(5,6,8)+148*p_(0,7,8)+168*p_(1,7,8)+204*p_(2,7,8)-57*p_(3,7,8)-59*p_(4,7,8)+231*p_(5,7,8)+150*p_(6,7,8)+16*p_(0,1,9)+40*p_(0,2,9)-57*p_(1,2,9)-150*p_(0,3,9)+148*p_(1,3,9)+148*p_(2,3,9)-16*p_(0,4,9)-148*p_(1,4,9)+40*p_(2,4,9)-184*p_(3,4,9)-168*p_(0,5,9)+57*p_(1,5,9)-148*p_(2,5,9)-150*p_(3,5,9)+59*p_(4,5,9)+168*p_(0,6,9)-204*p_(1,6,9)-199*p_(2,6,9)+184*p_(3,6,9)+199*p_(4,6,9)+204*p_(5,6,9)-40*p_(0,7,9)+16*p_(1,7,9)+199*p_(2,7,9)+148*p_(3,7,9)+231*p_(4,7,9)+231*p_(5,7,9)+184*p_(6,7,9)+150*p_(0,8,9)+204*p_(1,8,9)-57*p_(2,8,9)+148*p_(3,8,9)+231*p_(4,8,9)-59*p_(5,8,9)+168*p_(6,8,9)-59*p_(7,8,9);

scan(11, i -> a_i = random(-100, 100))
f = -a_0*(p_(0,2,4) + p_(0,4,9) + p_(2,6,8) + p_(3,5,6)) + a_0*(p_(0,1,9) + p_(1,3,7) + p_(1,7,9) + p_(2,4,8) + p_(3,5,7) + p_(5,6,8)) - a_1*(p_(0,1,7) + p_(0,7,9) + p_(1,3,5) + p_(1,5,7) + p_(2,4,6) + p_(3,6,8)) + a_1*(p_(0,2,9) + p_(2,4,9) + p_(3,5,8) + p_(4,6,8)) - a_2*(p_(0,4,7) + p_(1,2,9) + p_(1,5,6) + p_(2,3,8) + p_(2,8,9) + p_(3,7,8)) + a_2*(p_(0,3,7) + p_(0,4,6) + p_(1,5,9) + p_(4,5,6)) - a_3*(p_(0,2,3) + p_(1,2,3) + p_(1,2,6) + p_(4,7,8) + p_(5,8,9) + p_(7,8,9)) + a_3*(p_(0,1,6) + p_(0,3,6) + p_(4,5,7) + p_(4,5,9)) - a_4*(p_(0,2,8) + p_(0,3,8) + p_(1,4,5) + p_(1,4,9) + p_(2,5,8) + p_(2,5,9) + p_(4,6,7) + p_(5,6,7)) + a_4*(p_(0,1,4) + p_(0,4,8) + p_(0,6,7) + p_(0,7,8) + p_(1,2,5) + p_(1,3,9) + p_(1,4,6) + p_(2,3,9) + p_(2,5,6) + p_(3,6,7) + p_(3,7,9) + p_(3,8,9)) - a_5*(p_(0,3,9) + p_(1,6,7) + p_(2,4,5) + p_(3,5,9)) + a_5*(p_(0,6,8) + p_(0,8,9) + p_(1,2,4) + p_(1,4,7) + p_(2,3,5) + p_(6,7,8)) - a_6*(p_(0,3,5) + p_(0,5,9) + p_(1,6,8) + p_(3,4,5)) + a_6*(p_(0,6,9) + p_(1,2,7) + p_(1,7,8) + p_(2,3,4) + p_(2,4,7) + p_(6,8,9)) - a_7*(p_(0,2,5) + p_(0,5,8) + p_(2,6,7) + p_(3,4,9)) + a_7*(p_(0,1,8) + p_(1,3,4) + p_(1,4,8) + p_(2,5,7) + p_(3,6,9) + p_(6,7,9)) - a_8*(p_(0,2,7) + p_(1,5,8) + p_(2,6,9) + p_(3,4,6)) + a_8*(p_(0,1,5) + p_(0,5,7) + p_(1,3,8) + p_(2,7,9) + p_(3,4,8) + p_(4,6,9)) - a_9*(p_(0,3,4) + p_(0,4,5) + p_(0,5,6) + p_(1,6,9)) + a_9*(p_(1,2,8) + p_(1,8,9) + p_(2,3,7) + p_(2,7,8) + p_(3,4,7) + p_(5,6,9)) - a_10*(p_(0,2,6) + p_(4,5,8)) + a_10*(p_(0,1,2) + p_(0,1,3) + p_(1,3,6) + p_(2,3,6) + p_(4,7,9) + p_(4,8,9) + p_(5,7,8) + p_(5,7,9))

inversePerm = (perm) -> (
	n := #perm;
	inv := new MutableList from toList(0..n-1);
	for i from 0 to n-1 do (
		inv#(perm#i) = i;
	);
	toList inv
);

countInversions = (perm) -> (
    ans := 0;
    for i from 0 to #perm-2 do (
        for j from i+1 to #perm-1 do (
            if perm#j < perm#i then ans = ans + 1;
        );
    );
    ans
);

-- Build lookup table once
varLookup := hashTable apply(gens R, v -> (
    toList (baseName v)#1 => v
));

applyPermutationToPoly = (f, perm) -> (
    subList := apply(gens R, v -> (
        indices := toList (baseName v)#1;
        newIndices := apply(indices, i -> perm#i);
        sortedIndices := sort newIndices;
        sign := (-1)^(countInversions(newIndices));
        v => sign * varLookup#sortedIndices
    ));
    sub(f, subList)
);

sumOverPermutations = (f, permList) -> (
  if ANTISYMMETRIZE then (
    return sum apply(permList, perm -> (-1)^(countInversions(perm)) * applyPermutationToPoly(f, perm))
  );
  if not ANTISYMMETRIZE then (
    return sum apply(permList, perm -> applyPermutationToPoly(f, perm))
  );
);


-- permList ={
-- 	{0,1,2,3,4,5,6,7,8,9},
-- 	{1,2,0,3,4,5,6,7,8,9},
-- 	{2,0,1,3,4,5,6,7,8,9}
-- };
-- symF = sumOverPermutations(f, permList);
symF = f;
f = symF;
-- f = p_(0,2,5) + p_(1,3,6) + p_(2,4,7) + p_(3,0,8) + p_(4,1,9) + p_(0,9,7) + p_(1,5,8) + p_(2,6,9) + p_(3,7,5) + p_(4,8,6);
-- f = -192*p_(0,1,2)-195*p_(0,1,3)+112*p_(0,2,3)+122*p_(1,2,3)-9*p_(0,1,4)-47*p_(0,2,4)-115*p_(1,2,4)-112*p_(0,3,4)-37*p_(1,3,4)+47*p_(2,3,4)-195*p_(0,1,5)-195*p_(0,2,5)-37*p_(1,2,5)-37*p_(0,3,5)-192*p_(1,3,5)+9*p_(2,3,5)+34*p_(0,4,5)+145*p_(1,4,5)-18*p_(2,4,5)-149*p_(3,4,5)-122*p_(0,1,6)+37*p_(0,2,6)+9*p_(1,2,6)-149*p_(0,3,6)+18*p_(1,3,6)+34*p_(2,3,6)+115*p_(0,4,6)+115*p_(1,4,6)-122*p_(2,4,6)+145*p_(3,4,6)+192*p_(0,5,6)+112*p_(1,5,6)-195*p_(2,5,6)-145*p_(3,5,6)+149*p_(4,5,6)-37*p_(0,1,7)-145*p_(0,2,7)+149*p_(1,2,7)-192*p_(0,3,7)-195*p_(1,3,7)+18*p_(2,3,7)+122*p_(0,4,7)+192*p_(1,4,7)-34*p_(2,4,7)+195*p_(3,4,7)+9*p_(0,5,7)+112*p_(1,5,7)+34*p_(2,5,7)+122*p_(3,5,7)+18*p_(4,5,7)-115*p_(0,6,7)-47*p_(1,6,7)-18*p_(2,6,7)-34*p_(3,6,7)-9*p_(4,6,7)-47*p_(5,6,7)+18*p_(0,1,8)+112*p_(0,2,8)-145*p_(1,2,8)+47*p_(0,3,8)+34*p_(1,3,8)+47*p_(2,3,8)+47*p_(0,4,8)-149*p_(1,4,8)-112*p_(2,4,8)-112*p_(3,4,8)+122*p_(0,5,8)+149*p_(1,5,8)-192*p_(2,5,8)-115*p_(3,5,8)-145*p_(4,5,8)+9*p_(0,6,8)-34*p_(1,6,8)+195*p_(2,6,8)-18*p_(3,6,8)+192*p_(4,6,8)+37*p_(5,6,8)+115*p_(0,7,8)+145*p_(1,7,8)-149*p_(2,7,8)-9*p_(3,7,8)+37*p_(4,7,8)+115*p_(5,7,8)-122*p_(6,7,8)+34*p_(0,1,9)+149*p_(0,2,9)-115*p_(1,2,9)+145*p_(0,3,9)-9*p_(1,3,9)+115*p_(2,3,9)+18*p_(0,4,9)-122*p_(1,4,9)-9*p_(2,4,9)-192*p_(3,4,9)-18*p_(0,5,9)-47*p_(1,5,9)-122*p_(2,5,9)-115*p_(3,5,9)-34*p_(4,5,9)-145*p_(0,6,9)-47*p_(1,6,9)+192*p_(2,6,9)+149*p_(3,6,9)+37*p_(4,6,9)+112*p_(5,6,9)-149*p_(0,7,9)-112*p_(1,7,9)+145*p_(2,7,9)-37*p_(3,7,9)+195*p_(4,7,9)+47*p_(5,7,9)-112*p_(6,7,9)-34*p_(0,8,9)-18*p_(1,8,9)+37*p_(2,8,9)+122*p_(3,8,9)+195*p_(4,8,9)+9*p_(5,8,9)+195*p_(6,8,9)+192*p_(7,8,9)

-- d = (degree f)#0;
-- degf = d;


-- this needs to have pluckerIdeal inside too
-- print "Smoothness check"
-- assert isSmooth(pluckerIdeal+ideal f,IsGraded=>true)

print "Divisibility check"

-- check Fern and my vanishings
for t from 1 to min(n-1,n //d) do (
    print(t);
    assert not (n % (t*d)==0 and k*d*t % n == 0)
);


assert (gcd(k,n//d) == 1);
assert (gcd(n-k,n//d) == 1); -- probably redundant?
assert (n<5 or n % 2 == 0 or (k!=2 and k!= n-2) or (n+1)//2 % d != 0);




-- e.g. p_(0,1) -> (0,1)
getIndex = (mon) -> (
    return (baseName mon)#1;
);

-- e.g. 5*p_(0,1)*p_(1,2)^2 -> ({1,0,2,0,0,...,0})
getMonomials = (polynomial) -> (
    return exponents (polynomial);
);

-- e.g. 5* p_(0,1)*p_(1,2)^2 -> (| p_(0,1)p_(1,2)^2 |, {3} | 5 |)
getCoefficients = (polynomial) -> (
    return coefficients(polynomial);
);

exponentToMonomial = (ls) -> (
    return product (for i from 0 to (numGens-1) list (gensR#i) ^ (ls#i));
);

-- this only works for p_I^j, not products
differentiateMonomial = (i,j,mon) -> (
    monomialIndex = (exponents (mon))#0;
    -- silly way to find the single exponent that shows up
    exponent = sum(monomialIndex);
    -- select the entry 
    ind = position((exponents mon)#0,lambda->lambda==exponent);
    monomial = gensR#ind;
    ls = getIndex(monomial);
    if (not member(j,ls)) then (
  return 0;
    );
    if (i!=j) and( member(i,ls)) then (
  return 0;
    );

    -- set j to i (replace copies then edits a list)
    return exponent * p_(replace(position(ls, lambda->lambda==j),i,ls)) * p_ls^(exponent-1);
);


differentiatePolynomial = (i,j,polynomial) -> (
    ans = 0;
    mons = terms(polynomial);
    for monomial in mons do (
  exponentList = ((listForm monomial)#0)#0;
  coeff = ((listForm monomial)#0)#1;
  
  for lsIndex from 0 to (numGens -1) do (
if ((exponentList#lsIndex)>0) then (
    mon = (gensR#lsIndex) ^(exponentList#lsIndex);
    -- this is an annoying way to implement the product rule
    monomialOverMon = new MutableList from exponentList;
    monomialOverMon#lsIndex = 0;
    monomialOverMon = exponentToMonomial(toList monomialOverMon);
    ans +=  coeff * differentiateMonomial(i,j,mon) * monomialOverMon;
);
  );
    );
    return ans;

);




-- Jacobian ideal I
I = ideal flatten (flatten for i from 0 to n-1 list differentiatePolynomial(i,i,f), 
    flatten flatten for i from 0 to n-1 list (
	for j from i+1 to n-1 list (
		{differentiatePolynomial(i,j,f),
		differentiatePolynomial(j,i,f)}
		-- differentiatePolynomial(i,i,f)-differentiatePolynomial(j,j,f)}
	)
));

-- I += ideal flatten for i from 0 to n-1 list differentiatePolynomial(i,i,f)

I = trim I;

gensI = gens I;

if (# (entries gens I)#0 != n^2) then error (
    print "Warning: Jacobian ideal does not have expected number of generators";
    exit 1;
) ;



J = R /  (I + pluckerIdeal);
-- S = R/(pluckerIdeal + J);
-- S = R/(pluckerIdeal + J + antisymmetrize_ideal);

-- Hodge numbers of primitive cohomology
-- (R_f)_{(p+1)d-n} = H^{N-1-p,p}
for i from 0 to n+1 list hilbertFunction((i+1)*d - n,J)


for i from 1 to n-1 do if i == k*(n-k)/2 then print concatenate("Warning: nontrivial cokernel contribution for i =",toString i) else continue


-- e.g. 3/3 == 1
for i from 1 to n-1 do if (i==((2*n-1-d)/3) or i==((4*n-9-d)/3)) then print concatenate("Potential error with i =",toString i) else continue


print "Expected Hodge numbers:"


-- k = 3
-- n = 7

-- d = toList(d)




G = flagBundle {k,n-k}; -- secretly Gr(k,n)
-- for d from 1 to 25 list (
-- X = sectionZeroLocus(sum for i from 1 to #d list OO_G(i)); -- quadric and a cubic
X = sectionZeroLocus(OO_G(d)); 

OmG = cotangentBundle G;
OmX = cotangentBundle X;

ls = for i from 0 to floor(dim X / 2) list
    abs(chi exteriorPower(i, OmX) - chi exteriorPower(i, OmG));

if (k*(n-k) % 2 ==0) then (
    print join(ls,reverse ls);
) else (
-- include the second half, but don't repeat middle element
print (join(ls, for i from 1 to #ls-1  list ls#(#ls-i-1)));
);

-- print "Smoothness check"
-- assert isSmooth(pluckerIdeal+ideal f,IsGraded=>true)



 getBasisInJ = (deg, J) -> (
    B := basis(deg, J);
    if B == 0 then return {};
    flatten entries B
);

-- Function to get the matrix representation of a permutation acting on a basis in J
permutationMatrixInJ = (perm, basisElems, J) -> (
    -- basisElems should be a list of elements in J
    m := #basisElems;
    
    -- Lift to R, apply permutation, reduce back to J
    images := apply(basisElems, b -> (
        bLifted := lift(b, J);
        imageInR := applyPermutationToPoly(bLifted, perm);
        sub(imageInR, J)
    ));
    
    -- Express each image as a linear combination of basis elements
    coeffMatrix := {};
    for img in images do (
        imgLifted := lift(img, J);
        coeffs := apply(basisElems, b -> (
            bLifted := lift(b, J);
            coefficient(bLifted, imgLifted)
        ));
        coeffMatrix = append(coeffMatrix, coeffs);
    );
    
    transpose matrix coeffMatrix
);

-- Get basis in J at a given degree
allBases = for i from 0 to n+1 list 
(
	deg := (i+1)*d - n;
	if deg >= 0 then basis(deg, J) else 0
);

-- Then extract as needed
getBasisInJ = (deg, J) -> 
(
	i := (deg + n - d) // d;
	if i >= 0 and i < #allBases and allBases#i != 0 then
		flatten entries allBases#i
	else
		{}
);





-- Function to create direct sum of permutation matrices over all relevant degrees
permutationMatrixAllDegrees = (perm, J, n, d) -> (
  matrices := {};
  for i from 0 to n+1 do (
    deg := (i+1)*d - n;
    if deg >= 0 then (
      B := getBasisInJ(deg, J);
      if #B > 0 then (
        M := permutationMatrixInJ(perm, B, J);
        if ANTISYMMETRIZE then (
          M *= (-1)^(countInversions(perm) * (i*d-n)); -- Adjust sign based on action on f
        );
        matrices = append(matrices, M);
      );
    );
  );

  if #matrices == 0 then (
    return matrix{{0}};
  );

  directSum matrices
);

-- Collect traces
traceList := {};

for perm in permList do (
  bigM := permutationMatrixAllDegrees(perm, J, n, d);
   -- Adjust sign based on action on Omega
   if k % 2 == 1 then (
      sign := (-1)^(countInversions(perm));
      bigM *= sign;
   );
  traceList = append(traceList, trace bigM);
);

print permList;
print "Traces (in order of permList):";
print traceList;
print (toString symF);

print eigenvalues permutationMatrixAllDegrees(permList#3,J,n,d)
for i from 0 to #permList-1 list symF - (-1)^(countInversions permList#i)*applyPermutationToPoly(symF,permList#i)

needsPackage "RationalPoints2"
rationalPoints(ideal f, Amount => true)