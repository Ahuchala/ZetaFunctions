-- p = 7
n = 5
k = 2
-- K = ZZ/p
K = QQ

-- access generators like p_(0,2)

-- R = K[p_(0,0) .. p_(n-1,n-1)]

pluckerIdeal = Grassmannian(k-1,n-1,CoefficientRing=>K)
R = ring pluckerIdeal




gensR = gens R
numGens = #gensR

-- hypersurface V(f)
f = random(2,R)
-- f = x01^2+2*x02^2+4*x03^2+5*x04^2+6*x12^2+11*x13^2+75*x14^2+13*x23^2+43*x24^2+8*x34^2

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



-- Jacobian ideal J
jacobianIdeal = ideal flatten flatten for i from 0 to n-1 list (
	for j from i+1 to n-1 list (
		{differentiatePolynomial(i,j,f),
		differentiatePolynomial(j,i,f), 
		differentiatePolynomial(i,i,f)-differentiatePolynomial(j,j,f)}
	)
)

jacobianIdeal = trim (jacobianIdeal + ideal f);



J = R / (pluckerIdeal + jacobianIdeal);
-- S = R/(pluckerIdeal + J);
-- S = R/(pluckerIdeal + J + antisymmetrize_ideal);

for i from 0 to 10 list hilbertFunction(i,J)

