loadPackage "Resultants";
loadPackage "TensorComplexes"; -- for multiSubsets function
loadPackage "Divisor";

-- q = 991
-- p = 7


n = 5;
k = 2;
-- K = ZZ/q
K = QQ;

prec = 1;

-- access generators like p_(0,2)

-- R = K[p_(0,0) .. p_(n-1,n-1)]

pluckerIdeal = Grassmannian(k-1,n-1,CoefficientRing=>K);
R = ring pluckerIdeal;




gensR = gens R;
numGens = #gensR;

-- f = -46*p_(0,1)^3-39*p_(0,2)^3+p_(0,1)*p_(0,2)*p_(1,2)+35*p_(1,2)^3+12*p_(0,3)^3+36*p_(1,3)^3+43*p_(2,3)^3+p_(0,4)^3-5*p_(1,4)^3+4*p_(2,4)^3-24*p_(3,4)^3

--       3 4         4       7 4       1 4       5 4       5 4       7 4       9 4      
-- o12 = -p      + 2p      + -p      + -p      + -p      + -p      + -p      + -p      +
--       7 0,1,2     0,1,3   5 0,2,3   3 1,2,3   9 0,1,4   2 0,2,4   8 1,2,4   5 0,3,4  
--       ------------------------------------------------------------------------------------
--       3 4        7 4
--       -p      + --p
--       7 1,3,4   10 2,3,4



-- f = -3/7*p_(0,1,2)^4 + 2*p_(0,1,3)^4 - 7/5*p_(0,2,3)^4 - 1/3*p_(1,2,3)^4 - 5/9*p_(0,1,4)^4 - 5/2*p_(0,2,4)^4 - 7/8*p_(1,2,4)^4 - 9/5*p_(0,3,4)^4 - 3/7*p_(1,3,4)^4 + 7/10*p_(2,3,4)^4



f = sum(apply(gensR, i->random(K)*(i)^2))
-- f += random(QQ)*gensR#0 * gensR#1 * gensR#2^2
-- hypersurface V(f)
-- f = random(2,R)
-- f = p_(0,1)^2 + 2*p_(0,2)^2 +4*p_(0,3)^2 + 5*p_(0,4)^2 + 6*p_(1,2)^2+11*p_(1,3)^2+75*p_(1,4)^2+13*p_(2,3)^2+43*p_(2,4)^2+8*p_(3,4)^2
-- f = x01^2+2*x02^2+4*x03^2+5*x04^2+6*x12^2+11*x13^2+75*x14^2+13*x23^2+43*x24^2+8*x34^2

d = (degree f)#0;
degf = d;


print "Smoothness check"
assert isSmooth(ideal f,IsGraded=>true)

print "Divisibility check"

-- check Fern and my vanishings
for t from 1 to n //d do (
    assert not (n % (t*d)==0 and k*d*t % n == 0)
);


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

-- I = trim I;
gensI = gens I;
-- I = trim (I + ideal f);



J = R /  (I + pluckerIdeal);
-- S = R/(pluckerIdeal + J);
-- S = R/(pluckerIdeal + J + antisymmetrize_ideal);

-- Hodge numbers of primitive cohomology
-- (R_f)_{(p+1)d-n} = H^{N-1-p,p}
for i from 0 to n-1 list hilbertFunction((i+1)*d - n,J)


for i from 1 to n-1 do if i == k*(n-k)/2 then print concatenate("Warning: nontrivial cokernel contribution for i =",toString i) else continue


-- e.g. 3/3 == 1
for i from 1 to n-1 do if (i==((2*n-1-d)/3) or i==((4*n-9-d)/3)) then print concatenate("Potential error with i =",toString i) else continue
