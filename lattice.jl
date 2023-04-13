using Singular
using Polyhedra
# using LazySets
import GLPK
lib = DefaultLibrary{Int64}(GLPK.Optimizer)
# import Singular.lift_std_syz

p = 7
prec = 3

# use documentation from https://oscar-system.github.io/Singular.jl/latest/ideal/
# S = residue_ring(ZZ, p^prec)

# R, (w, x, y, z) = polynomial_ring(S, ["w", "x", "y", "z"])
# weight = [1,1,1,1]
# f = w^4 + x^4 + y^4 + z^4

R, (x, y, z) = polynomial_ring(S, ["x", "y", "z"])
weight = [1,1,1]
# f = x^4 + y^4 + z^4
f = x*y*z + 2*x^2 * y + 3*x^2 * 4*z+ x*y^2 + 5*y^2 * z + 6*x^3 + 7*y^3 + 8*z^3

Rgens = gens(R)
n = size(Rgens)[1]


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


vertices = [[[Int(j == i) for j = 1:n-1] for i = 1:n-1] ; [[0 for i = 1:n-1]]]
for i = 1:n-1
	for j = 1:n-1
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
if n == 3
	J, (x,y,z) = QuotientRing(R,std(I))
else
	J, (w,x,y,z) = QuotientRing(R,std(I))
end


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
	cohomology_basis = spoly{n_Zn}[]
	for scale = 1:n-1
		# P_int = list(LatticePolytope(scale_by_d(scale)).interior_points())
        # monomial_int = [I.reduce(affine_vector_to_monomial(_,scale)) for _ in P_int]
		vertex_set = scale_by_d(scale)
		P = polyhedron(vrep(vertex_set),lib)

		# build a net of integer points and intersect with P?

		# probably should start at 0 but can rule those out as not lying on interior
		tuples = Base.Iterators.product([1:(maximum(maximum(vertex_set))-1) for i in 1:n-1]...) #max(vert)minus 1 is ok??
		net = polyhedron(vrep(vec([[i for i in aa] for aa in tuples])),lib)
		Net = [a for a in points(net)]

		# only keep the points in the interior
		Net = filter(a -> all([ininterior(a,h) for h in halfspaces(hrep(P))]), Net)
		integer_points = Net

		# integer_points = [point for point in points(intersect(P,polyhedron(vrep(vec([[i for i in aa] for aa in tuples])),lib)))]
		# integer_points = filter(a -> all([ininterior(a,blah) for blah in halfspaces(hrep(P))]), integer_points)
		# println(integer_points)
		# integer_points = [all([ininterior(a,blah) for blah in halfspaces(hrep(P))]) ? a : -100 for a in integer_points]
		# deleteat!(integer_points, integer_points .== -100); #this is hacky

		# need to check if each point added is linearly independent of the current basis
		primitive_ideal = Ideal(J)
		cohomology_basis_pure = spoly{n_Zn}[]
		for a in integer_points
			b = affine_vector_to_monomial(Int.(a),scale)
			if !contains(primitive_ideal, Ideal(J,b))
				push!(cohomology_basis_pure,b)
				primitive_ideal = Ideal(J,cohomology_basis_pure)
		# cohomology_basis = [cohomology_basis ; [affine_vector_to_monomial(Int.(a),scale) for a in integer_points]]
        # want to find (P_int + If)/(If)
			end
		end
		cohomology_basis = [cohomology_basis; cohomology_basis_pure]
	end
	return cohomology_basis
end

B = compute_primitive_cohomology()
println(B)
println("actual/expected cohomology basis size: ",size(B,1)," ", (-1)^(n)*Integer(((1-d)^(n)-1)/d+1))


