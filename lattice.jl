using Singular
using Polyhedra
using Base.Iterators
import GLPK
lib = DefaultLibrary{Int64}(GLPK.Optimizer)

DEBUG = true


p = 11
# p = 13
# n = 4
# d = 4
# len_B = (-1)^(n)*Integer(((1-d)^(n)-1)/d+1)

include("aux_functions.jl")
# include("compute_prec.jl")

# desired_prec = Int(ceil(log(2*fdegree/floor(len_B/2))/log(p) + (n)*floor(len_B/2)/2)) #???
# desired_prec = min(desired_prec,Int(ceil(log(2*binomial(len_B,(len_B+1)÷2)*p^((n/2)*((len_B+1)÷2)))/log(p))))
# arithmetic_prec = find_s(desired_prec)
# desired_prec = 3
# arithmetic_prec = 2

# frob_prec = arithmetic_prec
# prec = desired_prec

prec = 12
frob_prec = 8

# prec = n
# prec = n+1

# use documentation from https://oscar-system.github.io/Singular.jl/latest/ideal/
zp = residue_ring(ZZ, p^prec)

# R, (w, x, y, z) = polynomial_ring(zp, ["w", "x", "y", "z"])
# weight = [1,1,1,1]
# f = w^4+2*w*x^3-2*x^4-x^3*y-x^2*y^2-y^4+w^3*z-x^3*z-2*w^2*y*z+2*w*x*y*z-x^2*y*z-w*y^2*z+2*x*y^2*z-2*y^3*z-w^2*z^2-2*w*x*z^2+x^2*z^2-2*w*y*z^2+x*y*z^2+y^2*z^2+2*w*z^3+2*x*z^3-2*y*z^3-2*z^4
# f = w^4 + x^4 + y^4 + z^4


# TODO: check smoothness, nondegeneracy

R, (x, y, z) = polynomial_ring(zp, ["x", "y", "z"])
# weight = [1,3,1]
# f = y^2 - x^6 - x^3*z^3
weight = [1,1,1]
# f = x^5 + y^5 + z^5+3*x*y*z^3
f = x^4 + y^4 + z^4-x*y*z^2+4*x^2*z^2
# f = x^3 + y^3 + z^3
# f = x*y*z + 2*x^2 * y + 3*x^2 * 4*z+ x*y^2 + 5*y^2 * z + 6*x^3 + 7*y^3 + 8*z^3

Rgens = gens(R)
n = size(Rgens)[1]



S, (x,y,z,u1, u2, u3) = polynomial_ring(zp, ["x","y","z","u1", "u2", "u3"])
U = [u1,u2,u3]
if n == 4
	S, (w,x,y,z,u1, u2, u3, u4) = polynomial_ring(zp, ["w","x","y","z","u1", "u2", "u3", "u4"])
	U = [u1,u2,u3,u4]
end
Sgens = gens(S)

f_coeff = collect(coefficients(f))
f_mons = collect(monomials(f))
# f = sum([(Integer(BigInt(f_coeff[i]))%p)*f_mons[i] for i in 1:size(f_coeff,1)])

d = total_degree(f)
fdegree = d

# k = (-1)^(n)*Integer(((1-d)^(n)-1)/d+1)
# frob_prec = Int(ceil(log(2*binomial(k,(k+1)÷2)*p^((n/2)*((k+1)÷2)))/log(p)))
println("Using arithmetic with precision p^",prec," and frobenius precision ", frob_prec)

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
	if DEBUG
		@assert !isempty(v) && all(>=(0), v)
		@assert length(v) == n
	end
	return prod(Rgens.^v)
end

# input is like [[[expvector], coeff], [[expvector], coeff],...]
function vector_to_polynomial(v)
	if isempty(v)
		return R(0)
	end
	return sum([vector_to_monomial(a[1]) * a[2] for a in v])
end

function vector_to_monomial_J(v::Vector{<:Integer})
	if DEBUG
		@assert !isempty(v) && all(>=(0), v)
		@assert length(v) == n
	end
	return prod(Jgens.^v)
end

function vector_to_polynomial_J(v)
	if isempty(v)
		return J(0)
	end
	return sum([vector_to_monomial_J(a[1]) * a[2] for a in v])
end

# change a polynomial g to coefficients in ring B
function change_ring(g, ringB)
	if iszero(g)
		return ringB(0)
	end
	coeffs = collect(coefficients(g))
	mons = collect(monomials(g))
	return sum([coeffs[i] * prod(gens(ringB)[1:n].^collect(exponent_vectors(mons[i]))[1][1:n]) for i =1:size(coeffs,1)])
end


f_S = poly_R_to_S(f)
df_S = R_to_S.(df)
I_S = Ideal(S,df_S)
I_S_std = std(I_S)



fweight = sum(weight.*polynomial_to_vector(f)[2][1])

if DEBUG
	@assert all([fweight==sum(weight.*polynomial_to_vector(f)[2][i]) for i = 1:size(polynomial_to_vector(f),1)])
end

vertices = [[[Int(j == i) for j = 1:n-1] for i = 1:n-1] ; [[0 for i = 1:n-1]]]
for i = 1:n-1
	for j = 1:n-1
		vertices[i][j] *= fweight//weight[i+1]
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
    if g == 0
    	return -1
    elseif g == 1
        return 0
    end
#     shouldn't matter which one we pick
    mon = collect(monomials(g))[1]
    return sum(weight .* monomial_to_vector(mon)) ÷ fweight
end



function compute_primitive_cohomology()
	cohomology_basis = spoly{n_Zn}[]
	for scale = 1:n

		vertex_set = scale_by_d(scale)
		m_v_s = maximum(maximum.(vertex_set))-1
		P = polyhedron(vrep(vertex_set),lib)

		# build a net of integer points and intersect with P

		# probably should start at 0 but can rule those out as not lying on interior
		tuples = Base.Iterators.product([1:m_v_s for i in 1:n-1]...) #max(vert)minus 1 is ok??
		integer_points = polyhedron(vrep(vec([[i for i in aa] for aa in tuples])),lib)
		integer_points = [a for a in points(integer_points)]

		# only keep the points in the interior
		integer_points = filter(a -> all([ininterior(a,h) for h in halfspaces(hrep(P))]), integer_points)

		# need to check if each point added is linearly independent of the current basis
		primitive_ideal = Ideal(J) #TODO: groebner?
		primitive_ideal_std = std(primitive_ideal)
		cohomology_basis_pure = spoly{n_Zn}[]
		for a in integer_points
			b = affine_vector_to_J_monomial(Int.(a),scale)
			if Singular.reduce(b,primitive_ideal_std) != 0
			# if !contains(primitive_ideal,Ideal(J,b))
				push!(cohomology_basis_pure,change_ring(b,R))
				# want to find (P_int + If)/(If)
				primitive_ideal = Ideal(J,cohomology_basis_pure)
				primitive_ideal_std = std(primitive_ideal)
			end
		end
		cohomology_basis = [cohomology_basis; cohomology_basis_pure]
	end
	return cohomology_basis
end

B = compute_primitive_cohomology()
println(B)
len_B = size(B,1)

if DEBUG
	if weight != [1 for i = 1:n]
		println("Warning: basis length not checked for weight ", weight)
	else
		@assert size(B,1) == (-1)^(n)*Integer(((1-d)^(n)-1)/d+1)
	end
end

max_degree_B = maximum([degree(a) for a in B])

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
	for j = 0:frob_prec-1
		if DEBUG
			@assert binomial(-d,j) * binomial(d+frob_prec-1,d+j) == big(binomial(-d,j)) * big(binomial(d+frob_prec-1,d+j))
		end
		summer +=  binomial(-d,j) * binomial(d+frob_prec-1,d+j) * sigma_g*sigma(fj)
		fj *= f
	end
	return summer
end

function frobenius_on_cohom(i)
	return frobenius(B[i]) * p^(n-2)
end


# qr_dict = Dict{spoly{n_Zn},Vector{Vector{Vector{Any}}}}
P1 = nothing
# TODO: sizehint?
qr_dict = Dict()

# don't think I need scale = 0??
for scale = 1:n#+1
	vertex_set = scale_by_d(scale)
	P = polyhedron(vrep(vertex_set),lib)
	tuples = Base.Iterators.product([0:(maximum(maximum.(vertex_set))) for i in 1:n-1]...)
	net = polyhedron(vrep(vec([[i for i in aa] for aa in tuples])),lib)
	Net = [Integer.(a) for a in points(net)]
	Net = filter(a -> all([in(a,h) for h in halfspaces(hrep(P))]), Net)
	# integer_points = Net
	monomial_int = [affine_vector_to_monomial(pt,scale) for pt in Net]
	if scale == 1
		global P1 = monomial_int
	end
	for monomial in monomial_int
		rem = Singular.reduce(monomial, I_std)  # poly, ideal
		quo = monomial - rem
		# this feels bound to cause a bug
		# TODO: lift groebner??
		quo = Singular.lift(I,Ideal(R,quo))[1][1] #Ideal, subideal
		quo =  [[[a[2],a[3]] for a in quo if a[1] == i] for i = 1 : n]
		# todo: not vectorized
		quo = [vector_to_polynomial(a) for a in quo]
		# now q looks like 

		# 3-element Vector{Vector{Vector{Any}}}:
		#  [[[3, 0, 0], 324], [[2, 1, 0], 98], [[1, 2, 0], 111], [[0, 3, 0], 259], [[2, 0, 1], 49], [[1, 1, 1], 200], [[0, 2, 1], 113], [[1, 0, 2], 166], [[0, 1, 2], 309], [[0, 0, 3], 169]]
		#  [[[3, 0, 0], 185], [[2, 1, 0], 173], [[1, 2, 0], 256], [[0, 3, 0], 98], [[2, 0, 1], 196], [[1, 1, 1], 325], [[0, 2, 1], 316], [[1, 0, 2], 171], [[0, 1, 2], 21], [[0, 0, 3], 55]]
		#  [[[3, 0, 0], 136], [[2, 1, 0], 173], [[1, 2, 0], 32], [[0, 3, 0], 329], [[2, 0, 1], 339], [[1, 1, 1], 33], [[0, 2, 1], 103], [[0, 0, 3], 13]]

		qr_dict[monomial] = [quo,rem]

	end
end

R_dict = Dict()


function compute_psi_basis()
	psi_basis = spoly{n_Zn}[R(1)]
	for scale = 1:n

		vertex_set = scale_by_d(scale)
		P = polyhedron(vrep(vertex_set),lib)

		# build a net of integer points and intersect with P

		# probably should start at 0 but can rule those out as not lying on interior
		tuples = Base.Iterators.product([0:(maximum(maximum.(vertex_set))) for i in 1:n-1]...)
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
					if Singular.reduce(monomial,primitive_ideal_std) != 0
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
len_psi_basis = size(psi_basis,1)
println(len_psi_basis)


# this has been rewritten and is probably correct

# In: u, h such that g = x^u * h
# Out: m where R(m) \in psi_basis
function reduce_monomial(u,h)

	# write h = \sum g_i x_i df/dx_i + \sum c_j mu_lj
	# then  x^u h \equiv x^u [1/(const) (sum_i ui gi  + xi dgi/dxi) + \sum c_j mu_lj ]
	

	# We'll treat h = h_S as an element of S, and use h_R = R(h)


	# We can assume parent(u) == S
	if DEBUG
		@assert parent(u[1]) == S
		@assert parent(h) == S
	end

	h_R = change_ring(h,R)

	# since h is a monomial, we can consider only the U term in front
	u_coefficient = prod(Sgens[n+1:end].^(collect(exponent_vectors(h))[1][n+1:end]))


	if DEBUG
		@assert change_ring(h_R,S) * u_coefficient == h
	end

	# if h_R is in psi_basis, then g = x^u h is already reduced
	if Singular.reduce(h_R,I_std) == h_R
		return h
	end

	# we write h_R = sum_i gi xi df/dxi + (remainder)
	#              = q + r
	q,r = qr_dict[h_R]
	q_S = sum([u[i]*change_ring(q[i],S) for i = 1:n])

	ans = u_coefficient * (change_ring(r + toric_derivative_vector(q),S) + q_S)


    # now we need to reduce all monomials in ans
    ans_mons = collect(monomials(ans))
    ans_coeffs = collect(coefficients(ans))
    num_coeffs = size(ans_coeffs,1)

    return sum([ans_coeffs[i]*reduce_monomial(u,ans_mons[i]) for i = 1:num_coeffs])
end



s_dict = Dict()



# I think this calculation works
# I no longer think this works
function compute_s_dict_v(v)
	mat = [[S(0) for a = 1:len_psi_basis] for aa = 1:len_psi_basis]
	for i = 1:len_psi_basis
		row_vec = reduce_monomial(U, change_ring(v*psi_basis[i],S))
		row_vec_monomials = collect(monomials(row_vec))
		row_vec_coeffs = collect(coefficients(row_vec))
		for j = 1:size(row_vec_coeffs,1)
			h = row_vec_monomials[j]
			h_R = change_ring(h,R)
			u_coefficient = prod(Sgens[n+1:end].^(collect(exponent_vectors(h))[1][n+1:end]))
			mon_index = findall(xx->xx==h_R,psi_basis)

			if DEBUG
				if isempty(mon_index)
					println(h_R," not contained in ", psi_basis)
				end
			end
			
			mat[i][mon_index[1]] = mat[i][mon_index[1]] + u_coefficient * row_vec_coeffs[j]
		end
	end
	s_dict[v] = mat
end

# this has bugs!
function reduce_uv(monomial, coeff)
	u = monomial_to_vector(monomial)
	g = transpose([R(0) for a = 1:len_psi_basis])
	g[1] = coeff
	for v in P1

		v_vec = monomial_to_vector(v)
		if all(u.>=v_vec)
			if !haskey(s_dict,v)
				compute_s_dict_v(v)
			end
			# mat = s_dict[v]
			mat = transpose(reduce(hcat,s_dict[v]))
			while all(u.>=v_vec)

				u .-= v_vec

				m = 1
				if n == 3
					m = [evaluate(a,[0,0,0,u[1],u[2],u[3]]) for a in mat]
					# m = [[evaluate(a,[0,0,0,u[1],u[2],u[3]]) for a in aa] for aa in mat]
				elseif n == 4
					m = [evaluate(a,[0,0,0,0,u[1],u[2],u[3],u[4]]) for a in mat]
					# m = [[evaluate(a,[0,0,0,0,u[1],u[2],u[3],u[4]]) for a in aa] for aa in mat]
				else
					println("error: not implemented for n > 4")
				end
				# todo: learn how to write matrices properly

				# m = transpose(reduce(hcat, m))
				g = g*m
			end
		end
	end
	# why are these different!?!?

	# return sum(g .* psi_basis)
	return sum([g[ii] * psi_basis[ii] for ii = 1:len_psi_basis])

end

frob_matrix = [[zp(0) for i = 1:len_B] for j = 1:len_B]

# println(len_B)



graded_I_B = []
graded_I_B_std = []
for i = 1 : max_degree_B
	global graded_I_B = [graded_I_B ; Ideal(J,filter(a -> (degree(a) == i),B))]
	global graded_I_B_std = [graded_I_B_std; std(Ideal(J,filter(a -> (degree(a) == i),B)))]
end

cache_reduction = Dict()

# this is broken
# todo: reduce everything in psi_basis
# function reduce_to_basis(g)
# 	# println(g)
# 	# g  = change_ring(g,J)

# 	g_mons = collect(monomials(g))
# 	g_coeffs = collect(coefficients(g))
# 	ans = R(0)
# 	for j = 1 : size(g_coeffs,1)
# 		mon = g_mons[j]
# 		if !haskey(cache_reduction,mon)
# 			rem = Singular.reduce(mon, I_std)  # poly, ideal
# 			if rem == mon
# 				cache_reduction[mon] = mon
# 			else
# 				quo = mon - rem
# 				quo = Singular.lift(I,Ideal(R,quo))[1][1] #Ideal, subideal
# 				quo =  [[[a[2],a[3]] for a in quo if a[1] == i] for i = 1 : n]
# 				quo = [vector_to_polynomial(a) for a in quo]
# 				dict_ans = rem + toric_derivative_vector(quo)

# 				if degree(dict_ans) > 0
# 					dict_mons = collect(monomials(dict_ans))
# 					dict_coeffs = collect(coefficients(dict_ans))
# 					s = size(dict_coeffs,1)
# 					dict_ans = sum([dict_coeffs[i]*reduce_to_basis(dict_mons[i]) for i = 1:s])
# 				end

# 				cache_reduction[mon] = dict_ans
# 			end
# 		end
# 		ans += g_coeffs[j] * cache_reduction[mon]
# 	end
# 	return change_ring(ans,R)
# end

function final_reduce_to_basis(g)
	g  = change_ring(g,J)
	g_mons = collect(monomials(g))
	g_coeffs = collect(coefficients(g))
	ans = J(0)
	for i = 1 : size(g_coeffs,1)
		mon = g_mons[i]
		if !haskey(cache_reduction,mon)
			dict_ans = J(0)
			monomial_to_reduce = mon
			for i = max_degree_B:-1:1
				if !iszero(monomial_to_reduce)
					monomial_to_reduce = change_ring(monomial_to_reduce,J)
					ideal = graded_I_B[i]
					ideal_std = graded_I_B_std[i]
					rem = Singular.reduce(monomial_to_reduce,ideal_std)
					quo = Singular.lift(ideal,Ideal(J,monomial_to_reduce-rem))[1]
					s = size(gens(ideal),1)
					quo =  [[[a[2],a[3]] for a in quo[1] if a[1] == i] for i = 1 : s]
					quo = [vector_to_polynomial(a) for a in quo]
					
					dict_ans += sum([change_ring(quo[i],J)*gens(ideal)[i] for i = 1 : s])
					monomial_to_reduce = rem
				end
			end
			cache_reduction[mon] = dict_ans
		end
		ans += g_coeffs[i] * cache_reduction[mon]
	end
	return change_ring(ans,R)
end

for i = 1:len_B
    fro = frobenius_on_cohom(i)
    println(B[i])
    mons = collect(monomials(fro))
    coeffs = R.(collect(coefficients(fro)))
    for k = 1 : size(coeffs,1)
    	# println(fro)
    	# todo: sort out this division nonsense
    	# could and should also cache this factorial thing

    	h = reduce_uv(mons[k],coeffs[k])
    	if h != 0
    		# h = reduce_to_basis(h)
    		# ans_mons = collect(monomials(h))
	    	# ans_coeffs = collect(coefficients(h))
	    	# h = sum([factorial(max(0,degree(ans_mons[j])-1)) * ans_mons[j]*ans_coeffs[j] for j=1:size(ans_coeffs,1)])

    		h = final_reduce_to_basis(h)
	        ans_mons = collect(monomials(h))
	    	ans_coeffs = collect(coefficients(h))
	    	denom = zp(factorial(big(degree(mons[k])-1)))

	    	ans_coeffs = [numerator(a // denom)*zp(invmod(BigInt(denominator(a // denom)),p^prec)) for a in ans_coeffs]

	    	# now add to frob matrix
	    	for j = 1:size(ans_coeffs,1)
	    		mon_index = findall(xx->xx==ans_mons[j],B)
	    		if isempty(mon_index) #!?
	    			println("error: ", ans_mons[j])
	    		else
	    			frob_matrix[i][mon_index[1]] = frob_matrix[i][mon_index[1]]+ans_coeffs[j]
	    		end
	    	end
	    end
    end
end
        

println(frob_matrix)

