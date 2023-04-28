using Singular
using Polyhedra
using Base.Iterators
import GLPK
lib = DefaultLibrary{Int64}(GLPK.Optimizer)

p = 11
prec = 3

# use documentation from https://oscar-system.github.io/Singular.jl/latest/ideal/
zp = residue_ring(ZZ, p^prec)

# R, (w, x, y, z) = polynomial_ring(zp, ["w", "x", "y", "z"])
# weight = [1,1,1,1]
# f = w^4+2*w*x^3-2*x^4-x^3*y-x^2*y^2-y^4+w^3*z-x^3*z-2*w^2*y*z+2*w*x*y*z-x^2*y*z-w*y^2*z+2*x*y^2*z-2*y^3*z-w^2*z^2-2*w*x*z^2+x^2*z^2-2*w*y*z^2+x*y*z^2+y^2*z^2+2*w*z^3+2*x*z^3-2*y*z^3-2*z^4
# # f = w^4 + x^4 + y^4 + z^4


# TODO: check smoothness, nondegeneracy

R, (x, y, z) = polynomial_ring(zp, ["x", "y", "z"])
weight = [1,1,1]
# f = x^5 + y^5 + z^5
# f = x^4 + y^4 + z^4
f = x*y*z + 2*x^2 * y + 3*x^2 * 4*z+ x*y^2 + 5*y^2 * z + 6*x^3 + 7*y^3 + 8*z^3

Rgens = gens(R)
n = size(Rgens)[1]


S, (x,y,z,u1, u2, u3) = polynomial_ring(zp, ["x","y","z","u1", "u2", "u3"])
U = [u1,u2,u3]
if n == 4
	S, (w,x,y,z,u1, u2, u3, u4) = polynomial_ring(zp, ["w","x","y","z","u1", "u2", "u3", "u4"])
	U = [u1,u2,u3,u4]
end
Sgens = gens(S)


d = total_degree(f)
fdegree = d

df = [Rgens[i] * derivative(f,Rgens[i]) for i = 1:n]
I = Ideal(R,df)
I_std = std(I)
if n == 3
	J, (x,y,z) = QuotientRing(R, I_std)
else
	J, (w,x,y,z) = QuotientRing(R,I_std)
end

Jgens = gens(J)



# forgets coefficient
function monomial_to_vector(g)
	return collect(exponent_vectors(g))[1]
end

# remembers coefficient
function polynomial_to_vector(g)
	return [collect(coefficients(g)),collect(exponent_vectors(g))]
end

function vector_to_monomial(v::Vector{<:Integer})
	@assert !isempty(v) && all(>=(0), v)
	@assert length(v) == n
	return prod(Rgens.^v)
end

# input is like [[[expvector], coeff], [[expvector], coeff],...]
function vector_to_polynomial(v)
	if isempty(v)
		return R(0)
	end
	return sum([vector_to_monomial(a[1]) * a[2] for a in v])
end

# returns sum xi dg/dxi, using sum xi d(x^u)/dxi = sum(u) x^u
function toric_derivative(g)
	return [sum(collect(exponent_vectors(monomial))[1]) * monomial for monomial in monomials(g)]
end

# same as toric_derviative, but input is of the form [[[expvector], coeff], [[expvector], coeff],...]
function toric_derivative_vector(v)
	return [[a[1], sum(a[1])*a[2]] for a in v]
end


# change a polynomial g to coefficients in ring B
function change_ring(g, ringB)
	if g == 0
		return ringB(0)
	end
	coeffs = collect(coefficients(g))
	mons = collect(monomials(g))
	return sum([coeffs[i] * prod(gens(ringB)[1:n].^collect(exponent_vectors(mons[i]))[1][1:n]) for i =1:size(coeffs,1)])
end


# in: monomial in R
# out: monomial in J = R/I
function R_to_J(g)
	return change_ring(g,J)
end

function R_to_S(g)
	return change_ring(g,S)
end

function J_to_R(g)
	return change_ring(g,R)
end

# in: polynomial in R
# out: polynomial in J = R/I
function poly_R_to_J(g)
	return change_ring(g,J)
end

function poly_J_to_R(g)
	return change_ring(g,R)
end

function poly_R_to_S(g)
	return change_ring(g,S)
end

# f_J = poly_R_to_J(f)
# df_J = R_to_J.(df)
# I_J = Ideal(J,df_J)
# I_J_std = std(I_J)

f_S = poly_R_to_S(f)
df_S = R_to_S.(df)
I_S = Ideal(S,df_S)
I_S_std = std(I_S)



fweight = sum(weight.*polynomial_to_vector(f)[2][1])
@assert all([fweight==sum(weight.*polynomial_to_vector(f)[2][i]) for i = 1:size(polynomial_to_vector(f),1)])


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

function affine_vector_to_J_monomial(v::Vector{<:Integer}, l::T1 where {T1<:Integer})
	return Jgens[1]^(((l*fdegree)-sum([v[i-1] * weight[i] for i =2:n]))÷ weight[1]) * prod([Jgens[i+1]^v[i] for i =1:n-1])
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




function compute_primitive_cohomology()
	cohomology_basis = spoly{n_Zn}[]
	for scale = 1:n-1

		vertex_set = scale_by_d(scale)
		P = polyhedron(vrep(vertex_set),lib)

		# build a net of integer points and intersect with P

		# probably should start at 0 but can rule those out as not lying on interior
		tuples = Base.Iterators.product([1:(maximum(maximum(vertex_set))-1) for i in 1:n-1]...) #max(vert)minus 1 is ok??
		net = polyhedron(vrep(vec([[i for i in aa] for aa in tuples])),lib)
		Net = [a for a in points(net)]

		# only keep the points in the interior
		Net = filter(a -> all([ininterior(a,h) for h in halfspaces(hrep(P))]), Net)
		integer_points = Net

		# need to check if each point added is linearly independent of the current basis
		primitive_ideal = Ideal(J)
		cohomology_basis_pure = spoly{n_Zn}[]
		for a in integer_points
			b = affine_vector_to_monomial(Int.(a),scale)
			if !contains(primitive_ideal, Ideal(J,b))
				push!(cohomology_basis_pure,b)
				# want to find (P_int + If)/(If)
				primitive_ideal = Ideal(J,cohomology_basis_pure)
			end
		end
		cohomology_basis = [cohomology_basis; cohomology_basis_pure]
	end
	return cohomology_basis
end

B = compute_primitive_cohomology()
println(B)
len_B = size(B,1)
@assert len_B == (-1)^(n)*Integer(((1-d)^(n)-1)/d+1)

function sigma(g) 
	mons = collect(monomials(g)) # I really really hope these are always in the same order
	coeffs = collect(coefficients(g))
	return sum([coeffs[i] * mons[i]^p for i = 1:size(coeffs,1)])
end

# TODO: distinguish between this precision and prec
function frobenius(g)
	d = degree(g)
	summer = 0
	sigma_g = sigma(g)
	fj = R(1)
	for j = 0:prec-1
		summer +=  binomial(-d,j) * binomial(d+prec-1,d+j) * sigma_g*sigma(fj)
		fj *= f
	end
	return summer
end

function frobenius_on_cohom(i)
	return frobenius(B[i]) * p^(n-2)
end


# qr_dict = Dict{spoly{n_Zn},Vector{Vector{Vector{Any}}}}
# syntax: merge!(qr_dict, Dict(x => [x,y,z]))
P1 = nothing
qr_dict = Dict()

# don't think I need scale = 0??
for scale = 1:n#+1
	vertex_set = scale_by_d(scale)
	P = polyhedron(vrep(vertex_set),lib)
	tuples = Base.Iterators.product([0:(maximum(maximum(vertex_set))) for i in 1:n-1]...)
	net = polyhedron(vrep(vec([[i for i in aa] for aa in tuples])),lib)
	Net = [a for a in points(net)]
	Net = filter(a -> all([in(a,h) for h in halfspaces(hrep(P))]), Net)
	integer_points = Net
	monomial_int = [affine_vector_to_monomial(Int.(pt),scale) for pt in integer_points]
	if scale == 1
		global P1 = monomial_int
	end
	for monomial in monomial_int
		r = Singular.reduce(monomial, I_std)  # poly, ideal
		q = monomial - r
		q = Singular.lift(I,Ideal(R,q))[1][1] #Ideal, subideal

		# TODO: clean this up, should just be vectors of ints
		# r = collect(exponent_vectors(r))
		# q now looks like

		# 28-element Vector{Any}:
		#  (3, [3, 0, 0], 136)
		#  (2, [3, 0, 0], 185)
		#  (1, [3, 0, 0], 324)
		#  (3, [2, 1, 0], 173)
		#  (2, [2, 1, 0], 173)
		#  (1, [2, 1, 0], 98)
		#  (3, [1, 2, 0], 32)
		#  (2, [1, 2, 0], 256)
		#  (1, [1, 2, 0], 111)
		#  (3, [0, 3, 0], 329)
		#  (2, [0, 3, 0], 98)
		#  (1, [0, 3, 0], 259)
		#  (3, [2, 0, 1], 339)
		#  ⋮

		# and r like

		# 31*z^6

		# r = polynomial_to_vector(r)
		# now r looks like

# 		2-element Vector{Vector}:
		 # n_Zn[286, 231, 286, 116, 117]
		 # [[1, 2, 0], [0, 3, 0], [1, 1, 1], [0, 2, 1], [0, 0, 3]]


		# q = collect(q)
		q =  [[[a[2],a[3]] for a in q if a[1] == i] for i = 1 : n]
		# q = collect(flatten(q))
		q = [vector_to_polynomial(a) for a in q]
		# println(q)
		# return
		# now q looks like 

		# 3-element Vector{Vector{Vector{Any}}}:
		#  [[[3, 0, 0], 324], [[2, 1, 0], 98], [[1, 2, 0], 111], [[0, 3, 0], 259], [[2, 0, 1], 49], [[1, 1, 1], 200], [[0, 2, 1], 113], [[1, 0, 2], 166], [[0, 1, 2], 309], [[0, 0, 3], 169]]
		#  [[[3, 0, 0], 185], [[2, 1, 0], 173], [[1, 2, 0], 256], [[0, 3, 0], 98], [[2, 0, 1], 196], [[1, 1, 1], 325], [[0, 2, 1], 316], [[1, 0, 2], 171], [[0, 1, 2], 21], [[0, 0, 3], 55]]
		#  [[[3, 0, 0], 136], [[2, 1, 0], 173], [[1, 2, 0], 32], [[0, 3, 0], 329], [[2, 0, 1], 339], [[1, 1, 1], 33], [[0, 2, 1], 103], [[0, 0, 3], 13]]

		qr_dict[monomial] = [q,r]

		# merge!(qr_dict, Dict(monomial => [q,r]))

	end
	# println(qr_dict)
end


R_dict = Dict()


# In: u, h such that g = x^u * h
# Out: 
function reduce_monomial(u,h)
 
	# write h = \sum g_i x_i df/dx_i + \sum c_j mu_lj
	# then  x^u h \equiv 1/(const) * x^u (sum_i ui gi  + xi dgi/dxi) + \sum c_j mu_lj
	# println(h, parent(h))
	h_R = h
	ring = parent(h)
	u_coefficient = 1
	if parent(u[1]) == S
		ring = S
		ideal = I_S_std
		h = change_ring(h,S)
		h_R = change_ring(h,R)
		if change_ring(h_R,ring) != h
			# literally h / h_R
			u_coefficient = prod(Sgens[n+1:end].^(collect(exponent_vectors(h))[1][n+1:end]))
		end
	end


	if reduce(h,ideal) == h
		# I'm confused why this doesn't return x^u h
		return h
	end
	# r and q are backwards in sage :(
	q,r = qr_dict[h_R]

	ans = u_coefficient * change_ring(r + sum([Rgens[i] * derivative(q[i],Rgens[i]) for i = 1:n]),ring)

	mons = collect(monomials(ans))
	coeffs = collect(coefficients(ans))
	# println(q,r,u,h)
	if ans != ring(0)
		ans = sum([coeffs[i] * reduce_monomial(u,mons[i]) for i = 1 : size(coeffs,1)])
	end

	for i = 1:n
    	ri_mons = [u_coefficient*change_ring(a,ring) for a in collect(monomials(q[i]))]
    	ri_coeffs = collect(coefficients(q[i]))
    	for j = 1 :size(ri_coeffs,1)
    		ans += u[i] * ri_coeffs[j] * reduce_monomial(u,ri_mons[j])
        end
    end
    return ans
end

function compute_psi_basis()
	psi_basis = spoly{n_Zn}[R(1)]
	for scale = 1:n-1

		vertex_set = scale_by_d(scale)
		P = polyhedron(vrep(vertex_set),lib)

		# build a net of integer points and intersect with P

		# probably should start at 0 but can rule those out as not lying on interior
		tuples = Base.Iterators.product([0:(maximum(maximum(vertex_set))) for i in 1:n-1]...)
		net = polyhedron(vrep(vec([[i for i in aa] for aa in tuples])),lib)
		Net = [a for a in points(net)]

		# only keep the points in the interior
		Net = filter(a -> all([in(a,h) for h in halfspaces(hrep(P))]), Net)
		integer_points = Net

		primitive_ideal = Ideal(J)
		primitive_ideal_std = std(primitive_ideal)


		psi_basis_pure = spoly{n_Zn}[]
		for a in integer_points
			b = affine_vector_to_J_monomial(Int.(a),scale)
			b_R = change_ring(b,R)
			for monomial in map(R_to_J, monomials(qr_dict[b_R][2]))
				if !(monomial in psi_basis_pure)
					if reduce(monomial,primitive_ideal_std) != 0
						push!(psi_basis_pure,monomial)
						primitive_ideal = Ideal(J,psi_basis_pure)
						primitive_ideal_std = std(primitive_ideal)
					end
				end
			end
		end
		psi_basis = [psi_basis; psi_basis_pure]
	end
	return [change_ring(psi,R) for psi in psi_basis]
end


psi_basis = compute_psi_basis()
println(size(psi_basis,1))


sv_dict = Dict()


function R_to_S(g)
	return prod(Sgens[1:n].^(collect(exponent_vectors(g))[1]))
end

function J_to_S(g)
	return change_ring(g,S)
end

function poly_J_to_S(g)
	return change_ring(g,S)
end

# psi_basis = [R_to_S(a) for a in psi_basis]
# P1 = [R_to_S(a) for a in P1]

# f_S = poly_J_to_S(f)
# df_S = [Sgens[i] * derivative(f_S,Sgens[i]) for i = 1:n]
# I_S = Ideal(S,df_S)
# I_S_std = std(I_S)


# for mon in psi_basis
#     for v in P1
#     	println(mon,v)
#         sv_dict[(mon,v)] = reduce_monomial(U, v*mon)
#     end
# end


# assemble info from sv_dict into a dict of matrices

s_dict = Dict()
len_psi_basis = size(psi_basis,1)

function compute_s_dict_v(v)
	mat = [[S(0) for a = 1:len_psi_basis] for aa = 1:len_psi_basis]
	for i = 1:len_psi_basis
		sv_dict[psi_basis[i],v] = reduce_monomial(U, v*psi_basis[i])
		row_vec = sv_dict[psi_basis[i],v]
		println(row_vec)
		row_vec_monomials = collect(monomials(row_vec))
		row_vec_coeffs = collect(coefficients(row_vec))
		for j = 1:size(row_vec_coeffs,1)
			h = row_vec_monomials[j]
			h_R = change_ring(h,R)
			u_coefficient = prod(Sgens[n+1:end].^(collect(exponent_vectors(h))[1][n+1:end]))
			mon_index = findall(xx->xx==h_R,psi_basis)

			# can this be empty if X is smooth?
			# if !isempty(mon_index)
				mat[i][mon_index[1]] = mat[i][mon_index[1]] + u_coefficient * row_vec_coeffs[j]
			# end
		end
	end
	s_dict[v] = mat
end

for v in P1
	compute_s_dict_v(v)
end


function reduce_uv(monomial, coeff)
	u = monomial_to_vector(monomial)
	g = transpose([R(0) for a = 1:len_psi_basis])
# 	 A = [0 1 ; 2 0]
# 2×2 Matrix{Int64}:
#  0  1
#  2  0
	g[1] = coeff
	for v in P1
		v_vec = monomial_to_vector(v)
		if !(v in keys(s_dict))
			compute_s_dict_v(v)
		end
		mat = s_dict[v]
		while all(u.>=v_vec)
			for ii = 1:n
				u[ii] -= v_vec[ii]
			end
			m = 1
			if n == 3
				m = [[R(a(0,0,0,u[1],u[2],u[3])) for a in aa] for aa in mat]
			elseif n == 4
				m = [[R(a(0,0,0,0,u[1],u[2],u[3],u[4])) for a in aa] for aa in mat]
				# a(0,0,0,1,1,1)
			end
			# todo: learn how to write matrices properly
			m = reduce(hcat, m)
			g = g*m# not transpose()
		end
	end
	h = sum([g[ii] * psi_basis[ii] for ii = 1:len_psi_basis])
    return h
end

frob_matrix = [[R(0) for i =1:len_B] for j =1:len_B]


for i = 1:len_B
    ans = R(0)
    fro = frobenius_on_cohom(i)
    println(B[i])
    mons = collect(monomials(fro))
    coeffs = R.(collect(coefficients(fro)))
    for k = 1 : size(coeffs,1)
    	println(fro)
    	# todo: sort out this division nonsense
        h = reduce_uv(mons[k],coeffs[k]) ÷ R(factorial(big(degree(mons[k])-1)))
        for b in B
            ans += b * (h.monomial_coefficient(b)%p^prec )
        end
        for j =1:len_B
            frob_matrix[i][j] = ans.monomial_coefficient(R(B[j])) % p^prec
        end
    end
end
        
# frob_matrix = matrix(frob_matrix)
println(frob_matrix)
# print(frob_matrix.characteristic_polynomial())
# print(frob_matrix.characteristic_polynomial() % p^prec)