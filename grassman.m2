-- p = 7
n = 5
k = 2
-- K = ZZ/p
K = QQ

-- access generators like p_(0,2)
-- I=Grassmannian(k-1,n-1,CoefficientRing=>F); 
-- R=ring I;
-- R = K[p_(0,0) .. p_(n-1,n-1)]


-- antisymmetrize_ideal = ideal flatten for i from 0 to n-1 list (
-- 	for j from i to n-1 list (
-- 		p_(i,j) + p_(j,i)
-- 	)
-- )

plucker_ideal = Grassmannian(k-1,n-1,CoefficientRing=>K)
R = ring plucker_ideal




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


-- act on x_rs * x_hk with D^i_j = x_i d/dx_j
-- from bottom of page 9 of Fatighenti + Mongardi
-- differentiateMonomial = (vi,vj,vr,vs,vh,vk) -> (
-- 	ans = 0;
-- 	if vj == vr then (
-- 		ans += p_(vi,vs) * p_(vh,vk);
-- 	);
-- 	if vj == vs then (
-- 		ans += p_(vr,vi) * p_(vh,vk);
-- 	);
-- 	if vj == vh then (
-- 		ans += p_(vr,vs) * p_(vi,vk);
-- 	);
-- 	if vj == vk then (
-- 		ans += p_(vr,vs) * p_(vh,vi);
-- 	);
-- 	return ans;
-- );


differentiatePolynomial = (i,j,f) -> (
	-- M,C -> monomial, coefficient matrices
	return sum for k from 0 to n-1 list (
		if k == j or k == i then 0 else p_(k,i) * diff(p_(k,j), f)
	)
);



-- for i from 0 to n-1 do (
-- 	for j from i+1 to n-1 do (
-- 		-- act on f with D^i_j = x_i d/dx_j
-- 		-- this will be zero except on terms with j in the index
-- 		-- ie a x_{jl} -> a x_{il}
-- 		ans = 0;
-- 		for l from j+1 to n-1 do (
-- 			ans += coefficient(f,p_(j,l)) * p_(i,l);
-- 		)
-- 		-- ls = append(ls,ans);
-- 	)
-- )

-- Jacobian ideal J
J = ideal flatten flatten for i from 0 to n-1 list (
	for j from 0 to n-1 list (
		if i == j then {0} else (
			{differentiatePolynomial(i,j,f), differentiatePolynomial(i,i,f)-differentiatePolynomial(j,j,f)}
		)
	)
)

J = trim J;



S = R / J;
-- S = R/(plucker_ideal + J);
-- S = R/(plucker_ideal + J + antisymmetrize_ideal);

for i from 0 to 10 list hilbertFunction(i,S)

