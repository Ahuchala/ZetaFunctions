using Singular
using Polyhedra
# using LazySets
import GLPK
lib = DefaultLibrary{Int64}(GLPK.Optimizer)
# import Singular.lift_std_syz

p = 7
prec = 3

# use documentation from https://oscar-system.github.io/Singular.jl/latest/ideal/
S = residue_ring(ZZ, p^prec)
R, (x, y, z) = polynomial_ring(S, ["x", "y", "z"])

Rgens = gens(R)
n = size(Rgens)[1]

weight = [1,1,1]
f = x*y*z + 2*x^2 * y + 3*x^2 * 4*z+ x*y^2 + 5*y^2 * z + 6*x^3 + 7*y^3 + 8*z^3
d = total_degree(f)
fdegree = d

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


vertices = [[[Int(j == i) for j = 1:n] for i = 1:n] ; [[0 for i = 1:n]]]
for i = 1:n
	for j = 1:n
		vertices[i][j] *= fweight//weight[i]
	end
end



function scale_by_d(scale)
	return [v.*scale for v in vertices]
end

println(scale_by_d(1))

# v a point in point in P_dl
function affine_vector_to_monomial(v::Vector{<:Integer}, l::T1 where {T1<:Integer})
	return Rgens[1]^(((l*fdegree)-sum([v[i-1] * weight[i] for i =2:n]))÷ weight[1]) * prod([Rgens[i+1]^v[i] for i =1:n-1])
end

# this only works for weighted projective space
function degree(g)
    if g == R(0)
    	return -1
    elseif g == R(1)
        return 0
    end
#     shouldn't matter which one we pick
    mon = collect(monomials(g))[1]
    return sum(weight .* monomial_to_vector(mon)) ÷ fweight
    #     return (g.degree()//fdegree
end

df = [Rgens[i] * derivative(f,Rgens[i]) for i = 1:n]
I = Ideal(R,df)

# I would love to involve groebner basis stuff but I think Singular.jl is scuffed

# I would love not to compute the syzigies of I but I think there's a bug
# import Singular.lift_std
# groebner_I, T = lift_std(I)

# should have like matrix(groebner_I) = matrix(I) * I_T
# groebner_I, I_T, syz = lift_std_syz(I,complete_reduction=true)
# println(groebner_I)
# J = QuotientRing(R,groebner_I)

# println(J)


function compute_primitive_cohomology()
	cohomology_basis = []
	for scale = 1:n-1
		# P_int = list(LatticePolytope(scale_by_d(scale)).interior_points())
        # monomial_int = [I.reduce(affine_vector_to_monomial(_,scale)) for _ in P_int]

		P = polyhedron(vrep(scale_by_d(scale)),lib)
		println(P)
        # want to find (P_int + If)/(If)
	end
	return cohomology_basis
end

B = compute_primitive_cohomology()


