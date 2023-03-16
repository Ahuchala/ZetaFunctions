using Nemo
p = 7
prec = 3

S = residue_ring(ZZ,p^prec)
R, (x,y,z) = S["x","y","z"]

Rgens = gens(R)
n = size(Rgens)[1]

weight = [1,1,1]
f = x*y*z + 2*x^2 * y + 3*x^2 * 4*z+ x*y^2 + 5*y^2 * z + 6*x^3 + 7*y^3 + 8*z^3
d = total_degree(f)

function monomial_to_vector(g)
	return collect(exponent_vectors(g))[1]
end

function polynomial_to_vector(g)
	return collect(exponent_vectors(g))
end

function vector_to_monomial(v::Vector{<:Integer})
	@assert !isempty(v) && all(>=(0), v)
	@assert length(v) == n
	return Rgens.^v
end

fweight = sum(weight.*polynomial_to_vector(f)[1])
@assert all([fweight==sum(weight.*polynomial_to_vector(f)[i]) for i = 1:size(polynomial_to_vector(f),1)])
