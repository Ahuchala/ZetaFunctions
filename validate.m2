needsPackage "Resultants";
needsPackage "TensorComplexes"; -- for multiSubsets function
needsPackage "Divisor";
needsPackage "Schubert2"; -- for Hodge number verification
allowableThreads = 4;
n = 10;
k = 3;

K = QQ;
d = 1;
permList = {
  {0,1,2,3,4,5,6,7,8,9},
  {0,6,3,2,4,9,1,8,7,5},
  {4,1,2,3,0,5,6,8,7,9},
  {4,6,3,2,0,9,1,7,8,5},
  {7,1,2,3,8,5,6,4,0,9},
  {7,6,3,2,8,9,1,0,4,5},
  {8,1,2,3,7,5,6,0,4,9},
  {8,6,3,2,7,9,1,4,0,5}
};
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


f = sum(apply(gens R, i->random(0,100)*(i)^d))

randomIntPoly = (d, R) -> 
(
	g := random(d, R);
	monoms := flatten entries monomials g;
	sum apply(monoms, m -> random(100) * m)
);

-- Usage:
f = randomIntPoly(d, R);



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
    sum apply(permList, perm -> applyPermutationToPoly(f, perm))
);


-- permList ={
-- 	{0,1,2,3,4,5,6,7,8,9},
-- 	{1,2,0,3,4,5,6,7,8,9},
-- 	{2,0,1,3,4,5,6,7,8,9}
-- };
symF = sumOverPermutations(f, permList);
f = symF;

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


needsPackage "Resultants";
needsPackage "TensorComplexes";
needsPackage "Divisor";
needsPackage "Schubert2";
allowableThreads = 4;
n = 10;
k = 3;
K = QQ;
d = 1;

-- Your permutation lists from SageMath
allPermLists = {
    {{0,1,2,3,4,5,6,7,8,9}, {2,1,0,3,4,5,6,7,8,9}},
    {{0,1,2,3,4,5,6,7,8,9}, {0,3,2,1,4,5,6,7,9,8}},
    {{0,1,2,3,4,5,6,7,8,9}, {1,0,8,5,9,3,7,6,2,4}},
    {{0,1,2,3,4,5,6,7,8,9}, {3,1,5,0,4,2,9,7,8,6}},
    {{0,1,2,3,4,5,6,7,8,9}, {8,7,2,5,9,3,6,1,0,4}}
    -- ... add all your permutation lists here
};

testPermList = (permList) -> (
    pluckerIdeal := Grassmannian(k-1,n-1,CoefficientRing=>K);
    R := ring pluckerIdeal;
    gensR := gens R;
    numGens := #gensR;
    
    -- Setup coordinate definitions
    for inds in multiSubsets(n,k) do (
        if (unique inds != inds) then (
            for perm in permutations(inds) do p_(toSequence perm) = 0;
        ) else for perm in permutations(inds) do (
            p_(toSequence perm) = (-1)^(countInversions(perm)) * p_(toSequence inds);
        );
    );
    
    -- Generate random symmetric polynomial
    f := randomIntPoly(d, R);
    f = sumOverPermutations(f, permList);
    
    -- Compute Jacobian ideal
    I := ideal flatten (flatten for i from 0 to n-1 list differentiatePolynomial(i,i,f), 
        flatten flatten for i from 0 to n-1 list (
            for j from i+1 to n-1 list {
                differentiatePolynomial(i,j,f),
                differentiatePolynomial(j,i,f)
            }
        ));
    I = trim I;
    
    return (# (entries gens I)#0 == n^2)
);

-- Test all permutation lists
validIndices = {};
for idx from 0 to #allPermLists-1 do (
    print("Testing permutation list " | toString idx | "...");
    if testPermList(allPermLists#idx) then (
        validIndices = append(validIndices, idx);
        print("  ✓ VALID");
    ) else (
        print("  ✗ Invalid number of generators");
    );
);

print("\nValid permutation list indices: " | toString validIndices);