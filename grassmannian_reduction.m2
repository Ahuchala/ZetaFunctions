q = 7
-- p = 7


n = 5
k = 2
-- K = ZZ/p
K = QQ

prec = 2

-- access generators like p_(0,2)

-- R = K[p_(0,0) .. p_(n-1,n-1)]

pluckerIdeal = Grassmannian(k-1,n-1,CoefficientRing=>K)
R = ring pluckerIdeal




gensR = gens R
numGens = #gensR

-- hypersurface V(f)
f = random(2,R)
-- f = x01^2+2*x02^2+4*x03^2+5*x04^2+6*x12^2+11*x13^2+75*x14^2+13*x23^2+43*x24^2+8*x34^2

d = (degree f)#0;
degf = d;
-- f = sum(apply(gensR, i->i^2))

for i from 0 to n-1 do (
	for j from i+1 to n-1 do (
		p_(j,i) = -p_(i,j)
	)
)
for i from 0 to n-1 do (
	p_(i,i) = 0
)



differentiatePolynomial = (i,j,f) -> (
	-- M,C -> monomial, coefficient matrices
	return sum for k from 0 to n-1 list (
		if k == j or k == i then 0 else p_(k,i) * diff(p_(k,j), f)
	)
);



-- Jacobian ideal I
I = ideal flatten flatten for i from 0 to n-1 list (
	for j from i+1 to n-1 list (
		{differentiatePolynomial(i,j,f),
		differentiatePolynomial(j,i,f), 
		differentiatePolynomial(i,i,f)-differentiatePolynomial(j,j,f)}
	)
)

-- I = trim I;
gensI = gens I;
-- I = trim (I + ideal f);


for i from 0 to n-1 list (
	for j from i+1 to n-1 list (
		{differentiatePolynomial(i,j,f),
		differentiatePolynomial(j,i,f), 
		differentiatePolynomial(i,i,f)-differentiatePolynomial(j,j,f)}
	)
)

-- for i from 0 to n-1 list (
-- 	for j from i+1 to n-1 list (
--         {(mingens I) #? (differentiatePolynomial(i,j,f)),
--         (mingens I) #? (differentiatePolynomial(j,i,f)),
--         (mingens I) #? (differentiatePolynomial(i,i,f)-differentiatePolynomial(j,j,f))}
--     )
-- )





J = R / (pluckerIdeal + I);
-- S = R/(pluckerIdeal + J);
-- S = R/(pluckerIdeal + J + antisymmetrize_ideal);

-- Hodge numbers of primitive cohomology
for i from 1 to n-1 list hilbertFunction((i+1)*d - n,J)

for i from 1 to n-1 do if i == k*(n-k)/2 then print concatenate("Warning: nontrivial cokernel contribution for i =",toString i) else continue


-- a in b -> b#?a
for i from 1 to n-1 do if ZZ#? ((2*n-1-d)/3) or ZZ#? ((4*n-9-d)/3) then print concatenate("Potential error with i =",toString i) else continue

use R
g = sum(apply(gensR, i-> i^8))

-- pth power on monomials, identity on coefficients
sigma = (g) -> (
	(mons,coeff) = coefficients(g);
	return (flatten entries (matrix {apply(flatten entries mons,i->i^q)} * coeff))#0;
);

frobenius = (b,prec) -> (
	degreeg = (degree b)#0;
	d = (degreeg) // degf;
	summer = 0;
	sigmag = sigma(b);
	fj = 1;
	for j from 0 to prec-1 do (
		const = binomial(-d,j) * binomial(d+prec-1,d+j);
		summer += const * sigmag * fj;
		fj *= f;
	);
	return summer;
);

-- In: polynomial g representing form g/f^d * omega
-- Out: equivalent elements of B modulo griffiths-dwork relations
reducePolynomial = (g) -> (
-- idea: naive griffiths-dwork reduction taking g in the jacobian ideal, ie
-- g = sum a_ij D_ij f + sum b_ij (D_ii f - D_jj f)
--   == sum D_ij a_ij + sum (D_ii - D_jj) b_ij
    ans = 0;
    degreeg = (degree g)#0;
    for i from 0 to (degreeg) // degf do (
        h = g % I;
        g -= h;
        ans += h;
        ls = flatten entries (g // gensI);
        g = 0;
        ind = 0;
        for i from 0 to n-1 do (
            for j from i+1 to n-1 do (
                g += differentiatePolynomial(i,j,ls#ind);
                g += differentiatePolynomial(j,i,ls#(ind+1));
                g += differentiatePolynomial(i,i,ls#(ind+2))-differentiatePolynomial(j,j,ls#(ind+2));
                ind += 3;
            );
        );
    );
	return ans;
);


frobOnCohom = (g) -> (
	return substitute(frobenius(g * (product gens R),prec) / (product gens R),R);
);

-- substitute((flatten entries basis(1,J))#0,R)
-- h = frobOnCohom(substitute((flatten entries basis(1,J))#0,R))
-- reducePolynomial(h)