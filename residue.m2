k = 2
n = 5

lsX = x_{1,1}..x_{n,n}

R = QQ[lsX]

-- degree of f in p_ij
d = 2


ls = subsets((for i from 1 to n list i),k)
lsMinusOne = subsets((for i from 1 to n list i-1),k)
mat = transpose matrix for j from 1 to k list for i from 1 to n list x_{i,j}

--initialize plucker coordinates p_I
for a from 0 to binomial(n,k)-1 do p_(ls#a) = det mat^(lsMinusOne#a)

for i from 1 to n do (
	for j from i+1 to n do (
		p_{j,i} = -p_{i,j}
	)
)
for i from 1 to n do (
	p_{i,i} = 0
)


-- f = p_{1,2}*p_{2,3} + p_{2,4} * p_{1,4}+p_{1,3}^2 +p_{1,5}^2
f = p_{1,2}^2 + 2*p_{1,3}^2 +4*p_{1,4}^2 + 5*p_{1,5}^2 + 6*p_{2,3}^2+11*p_{2,4}^2+75*p_{2,5}^2+13*p_{3,4}^2+43*p_{3,5}^2+8*p_{4,5}^2

Eij = (i,j,g) -> (
    return sum(
        for l from 1 to k list x_{j,l}*diff(x_{i,l},g)
    );
);

I = ideal flatten(for i from 1 to n list for j from 1 to n list Eij(i,j,f))
-- I = ideal flatten(for i from 1 to n list for j from 1 to n list (if i!= j then Eij(i,j,f) else 0))
-- I += ideal flatten(for i from 1 to n list for j from 1 to n list (if i!= j then Eij(i,i,f)-Eij(j,j,f) else 0))
-- I += ideal f
P = ideal flatten(for i from 1 to n list for j from i+1 to n list p_{i,j})
-- need to introduce plucker relations
-- pluckerIdeal = p_{1,2}*p_{3,4} - p_{1,3}*p_{2,4}+p_{1,4}*p_{2,3}
J = P/I;
-- toString hilbertSeries(J)

-- Eij(1,2,f) - (4*p_{1,3}*p_{2,3} + 8*p_{1,4}*p_{2,4} + 10 * p_{1,5}*p_{2,5})
Eij(4,4,f)-Eij(5,5,f)-(8*p_{1,4}^2-10*p_{1,5}^2+22*p_{2,4}^2-150*p_{2,5}^2+26*p_{3,4}^2-86*p_{3,5}^2)