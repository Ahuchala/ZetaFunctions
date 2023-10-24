using Singular
using Polyhedra
using Base.Iterators
import GLPK
lib = DefaultLibrary{Int64}(GLPK.Optimizer)

DEBUG = true


# p = 389
p = 7
# n = 4
# d = 4
# len_B = (-1)^(n)*Integer(((1-d)^(n)-1)/d+1)


# include("compute_prec.jl")

# desired_prec = Int(ceil(log(2*fdegree/floor(len_B/2))/log(p) + (n)*floor(len_B/2)/2)) #???
# desired_prec = min(desired_prec,Int(ceil(log(2*binomial(len_B,(len_B+1)÷2)*p^((n/2)*((len_B+1)÷2)))/log(p))))
# arithmetic_prec = find_s(desired_prec)
# desired_prec = 3
# arithmetic_prec = 2

# frob_prec = arithmetic_prec
# prec = desired_prec

# prec = 4
# prec = -1
frob_prec = 3
# prec = frob_prec+3
prec = 7


# prec = n
# prec = n+1

# use documentation from https://oscar-system.github.io/Singular.jl/latest/ideal/
# zp = residue_ring(ZZ, p^prec)

# residue_ring appears bugged, manually truncating precision with
modular_arithmetic = false
zp = QQ


if zp != QQ
	if big(p)^prec > -1 + 2^63
		zp = residue_ring(ZZ, big(p)^prec)
		println("warning: using big ints for large modulus ",big(p)^prec)
	end
end

# R, (w, x, y, z) = polynomial_ring(zp, ["w", "x", "y", "z"])
# weight = [1,1,1,1]
# f = w^4+2*w*x^3-2*x^4-x^3*y-x^2*y^2-y^4+w^3*z-x^3*z-2*w^2*y*z+2*w*x*y*z-x^2*y*z-w*y^2*z+2*x*y^2*z-2*y^3*z-w^2*z^2-2*w*x*z^2+x^2*z^2-2*w*y*z^2+x*y*z^2+y^2*z^2+2*w*z^3+2*x*z^3-2*y*z^3-2*z^4
# f = w^4 + x^4 + y^4 + z^4


# TODO: check smoothness, nondegeneracy

R, (x, y, z) = polynomial_ring(zp, ["x", "y", "z"])


# latest todo
# zp = ZZ
# R_ZZ, (x,y,z) = polynomial_ring(zp, ["x", "y", "z"])
# pn_ideal = Ideal(R_ZZ,R_ZZ(p^prec))
# R,(x, y, z) = QuotientRing(R_ZZ,pn_ideal)
# zp = QQ


weight = [1,3,1]
# f = y^2 - x^6 - x^3*z^3-z^6
f = y^2 - (-2*x^6-x^5*z+3*x^4*z^2+x^3*z^3-2*x^2*z^4+x*z^5+3*z^6)

# weight = [1,1,1]
# f = x^5 + y^5 + z^5+3*x*y*z^3

# f = x^4 + y^4 + z^4-x*y*z^2+4*x^2*z^2
# f = x^3 + y^3 + z^3
# f = x*y*z + 2*x^2 * y + 3*x^2 * 4*z+ x*y^2 + 5*y^2 * z + 6*x^3 + 7*y^3 + 8*z^3

Rgens = gens(R)
n = size(Rgens)[1]




f_coeff = collect(coefficients(f))
f_mons = collect(monomials(f))
# f = sum([(Integer(BigInt(f_coeff[i]))%p)*f_mons[i] for i in 1:size(f_coeff,1)])

d = total_degree(f)
fdegree = d
println(binomial(d*(n-1),n-1), ' ', binomial(d*(n-1)+n-1,n-1))

# k = (-1)^(n)*Integer(((1-d)^(n)-1)/d+1)
# frob_prec = Int(ceil(log(2*binomial(k,(k+1)÷2)*p^((n/2)*((k+1)÷2)))/log(p)))
println("Using arithmetic with precision p^",prec," and frobenius precision ", frob_prec)

include("aux_functions.jl")


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
	if DEBUG
		g_terms = collect(terms(g))
		@assert all([g_terms[i] == coeffs[i]*mons[i] for i = 1:size(coeffs,1)])
	end
	return sum([coeffs[i] * prod(gens(ringB)[1:n].^collect(exponent_vectors(mons[i]))[1][1:n]) for i =1:size(coeffs,1)])
end


# f_S = poly_R_to_S(f)
# df_S = R_to_S.(df)
# I_S = Ideal(S,df_S)
# I_S_std = std(I_S)



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
    if DEBUG
    	# @assert sum(weight .* monomial_to_vector(mon)) % fweight == 0
    	@assert all([sum(weight .* monomial_to_vector(mon)) == sum(weight .* monomial_to_vector(collect(monomials(g))[i])) for i = 2:size(collect(monomials(g)),1)])
    end
    return sum(weight .* monomial_to_vector(mon)) ÷ fweight
end


function polynomial_degree(g)
	return maximum([degree(mon) for mon in monomials(g)])
end

function compute_primitive_cohomology()
	cohomology_basis = typeof(Rgens[1])[]
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
		cohomology_basis_pure = typeof(Rgens[1])[]
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

# B = [z^8, x^2*y*z, x*y^2*z, x*y*z^2, y*z^7, x*z^7]	
# B= [x*y*z,z^6]
# B = [change_ring(blah,R) for blah in B]

# this is the notation of the paper
B = compute_primitive_cohomology()

# this makes it agree with the basis used by Singular lifts
B = Set(collect(flatten([collect(monomials(Singular.reduce(blah,I_std))) for blah in B])))
B = [blah for blah in B]


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

	if DEBUG
		g_terms = collect(terms(g))
		@assert all([g_terms[i] == coeffs[i]*mons[i] for i = 1:size(coeffs,1)])
	end
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

function divide_by_x1xn(g)
	coeffs = collect(coefficients(g))
	mons = collect(monomials(g))
	return sum([vector_to_monomial(monomial_to_vector(mons[i]).-1) * coeffs[i] for i = 1:size(coeffs,1)])
end

function frobenius_on_cohom(i)
	# return frobenius(B[i]) * p
	# return divide_by_x1xn(frobenius(B[i]) * p^(n-2))
	return frobenius(B[i]) * p^(n-2)
	# return frobenius(divide_by_x1xn(B[i]))* prod(Rgens)*p
	# this is what I think it should be...
	# return divide_by_x1xn(frobenius(B[i])) * p
end


vertex_set = scale_by_d(1)
P = polyhedron(vrep(vertex_set),lib)
tuples = Base.Iterators.product([0:(maximum(maximum.(vertex_set))) for i in 1:n-1]...)
net = polyhedron(vrep(vec([[i for i in aa] for aa in tuples])),lib)
Net = [Integer.(a) for a in points(net)]
Net = filter(a -> all([in(a,h) for h in halfspaces(hrep(P))]), Net)
# integer_points = Net
monomial_int = [affine_vector_to_monomial(pt,1) for pt in Net]
P1 = monomial_int
size_P1 = size(P1,1)

vertex_set = scale_by_d(n-1)
P = polyhedron(vrep(vertex_set),lib)
tuples = Base.Iterators.product([0:(maximum(maximum.(vertex_set))) for i in 1:n-1]...)
net = polyhedron(vrep(vec([[i for i in aa] for aa in tuples])),lib)
Net = [Integer.(a) for a in points(net)]
Net = filter(a -> all([in(a,h) for h in halfspaces(hrep(P))]), Net)
# integer_points = Net
monomial_int = [affine_vector_to_monomial(pt,n-1) for pt in Net]
Pn = monomial_int
size_Pn = size(Pn,1)

# TODO: sizehint?
qr_dict = Dict()

scale = n
vertex_set = scale_by_d(scale)
P = polyhedron(vrep(vertex_set),lib)
tuples = Base.Iterators.product([0:(maximum(maximum.(vertex_set))) for i in 1:n-1]...)
net = polyhedron(vrep(vec([[i for i in aa] for aa in tuples])),lib)
Net = [Integer.(a) for a in points(net)]
Net = filter(a -> all([in(a,h) for h in halfspaces(hrep(P))]), Net)
monomial_int = [affine_vector_to_monomial(pt,scale) for pt in Net]
for monomial in monomial_int

	# this appears to have bugs in positive characteristic

	rem = Singular.reduce(monomial, I_std)  # poly, ideal
	quo = monomial - rem
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

	# if modular_arithmetic
	# 	for i = 1 : n
	# 		quo_mons = collect(monomials(quo[i]))
	# 		quo_coeffs = collect(coefficients(quo[i]))
	# 		size_quo_coeffs = size(quo_coeffs,1)
	# 		if size_quo_coeffs > 0
	# 			quo_coeffs = [(mod(BigInt(numerator(a)), p^prec)) * BigInt(invmod(BigInt(denominator(a)),p^prec)) for a in quo_coeffs]
	# 			# println(quo)
	# 			# println(sum(quo_coeffs .* quo_mons))
	# 			quo[i] = sum(quo_coeffs .*  quo_mons)
	# 		end
	# 	end
	# 	rem_mons = collect(monomials(rem))
	# 	rem_coeffs = collect(coefficients(rem))
	# 	size_rem_coeffs = size(rem_coeffs,1)
	# 	if size_rem_coeffs > 0
	# 		rem_coeffs = [(mod(BigInt(numerator(a)), p^prec)) * BigInt(invmod(BigInt(denominator(a)),p^prec)) for a in rem_coeffs]
	# 		rem = sum(rem_coeffs .*  rem_mons)
	# 	end
	# end



	# 	println(quo)
	# 	quo = [(mod(BigInt(numerator(a)), p^prec)) * BigInt(invmod(BigInt(denominator(a)),p^prec)) for a in quo]
	# 	rem = (mod(BigInt(numerator(rem)), p^prec)) * BigInt(invmod(BigInt(denominator(rem)),p^prec))
	# end

	qr_dict[monomial] = [quo,rem]
end

# end

# keys are nu
# out is assembled like R_const + sum(mu .* R_mu) + j*S
R_mu_dict = Dict()
R_const_dict = Dict()
S_dict = Dict()


function compute_R_S(v)
	if !haskey(R_mu_dict,[v]) 

		mon_v = vector_to_monomial(v)
		mat_R_const = Array{typeof(zp(1))}(undef,size_Pn,size_Pn)
		mat_R_const .= ZERO
		mat_R_mu = Array{typeof(zp(1))}(undef,size_Pn,size_Pn,n)
		mat_R_mu .= ZERO
		mat_S = Array{typeof(zp(1))}(undef,size_Pn,size_Pn)
		mat_S .= ZERO

		for ind = 1 : size(Pn,1)

			vg = Pn[ind] * mon_v
			pi_n = qr_dict[vg][1]

			mat_R_const_poly = toric_derivative_vector(pi_n)
			mat_S_poly = sum(v.* pi_n)

			if mat_R_const_poly != 0

				mat_R_const_poly_mons = collect(monomials(mat_R_const_poly))
				mat_R_const_poly_coeffs = collect(coefficients(mat_R_const_poly))
				for i = 1 : size(mat_R_const_poly_coeffs,1)
					mon_ind = findall(xx->xx==mat_R_const_poly_mons[i],Pn)[1]
					mat_R_const[mon_ind,ind] = mat_R_const_poly_coeffs[i]
				end
			end

			for j = 1 : n
				if pi_n[j] != 0
					mat_R_mu_poly_mons = collect(monomials(pi_n[j]))
					mat_R_mu_poly_coeffs = collect(coefficients(pi_n[j]))
					for i = 1 : size(mat_R_mu_poly_mons,1)
						mon_ind = findall(xx->xx==mat_R_mu_poly_mons[i],Pn)[1]
						mat_R_mu[mon_ind,ind,j] = mat_R_mu_poly_coeffs[i]
					end
				end
			end

				
			if mat_S_poly != 0
				mat_S_poly_mons = collect(monomials(mat_S_poly))
				mat_S_poly_coeffs = collect(coefficients(mat_S_poly))

				for i = 1 : size(mat_S_poly_coeffs,1)
					mon_ind = findall(xx->xx==mat_S_poly_mons[i],Pn)[1]
					mat_S[mon_ind,ind] = mat_S_poly_coeffs[i]
				end
			end

		end

		# if modular_arithmetic
		# 	# 
		# 	R_mu_dict[v] = [(mod(BigInt(numerator(a)), p^prec)) * BigInt(invmod(BigInt(denominator(a)),p^prec)) for a in mat_R_mu]
		# 	R_const_dict[v] = [(mod(BigInt(numerator(a)), p^prec)) * BigInt(invmod(BigInt(denominator(a)),p^prec)) for a in mat_R_const]
		# 	S_dict[v] = [(mod(BigInt(numerator(a)), p^prec)) * BigInt(invmod(BigInt(denominator(a)),p^prec)) for a in mat_S]
		# else
			R_mu_dict[v] = mat_R_mu
			R_const_dict[v] = mat_R_const
			S_dict[v] = mat_S
		# end
	end

	return [R_const_dict[v], R_mu_dict[v], S_dict[v]]
end

	

# this has bugs?
function reduce_uv(monomial)
	# return monomial * coeff
	g = 0
	u = monomial_to_vector(monomial)

	for pn in Pn
		# if pn divides monomial
		if all(monomial_to_vector(pn) .<= u)
			g = pn
		end
	end
	if DEBUG
		if g == 0
			println("warning: no monomial divisors found for ", monomial)
		end
	end

	monomial = monomial / g

	u = monomial_to_vector(monomial)
	g = transpose(to_Pn_basis(g))
	# g = [R(0) for a = 1:size_Pn]
	# g[1] = coeff
	for v in P1

		v_vec = monomial_to_vector(v)
		if all(u.>=v_vec)
			mat_R_const, mat_R_mu ,mat_S = compute_R_S(v_vec)

			while all(u.>=v_vec) #&& degree(vector_to_monomial(u))> fdegree

				u .-= v_vec

				# m = 1
				# if n == 3
				# 	m = [evaluate(a,[0,0,0,u[1],u[2],u[3]]) for a in mat]
				# elseif n == 4
				# 	m = [evaluate(a,[0,0,0,0,u[1],u[2],u[3],u[4]]) for a in mat]
				# else
				# 	println("error: not implemented for n > 4")
				# end
				# m = R_mat
				g = (mat_R_const + sum(u[i] * mat_R_mu[:,:,i] for i = 1:n))*g
				# g = transpose(m)*g
			end
		end
	end
	return vector_to_monomial(u)*sum(g .* Pn)
end


# println(len_B)



graded_B = []
graded_B_std = []
graded_B_J = []
graded_B_J_std = []

B_J = [change_ring(blah,J) for blah in B]

# psi_basis_J = [change_ring(blah,J) for blah in psi_basis]
# max_degree_psi = maximum([degree(a) for a in psi_basis])

# for i = 0 : max_degree_psi
# 	global graded_psi = [graded_psi ; Ideal(J,filter(a -> (degree(a) == i),psi_basis_J))]
# 	global graded_psi_std = [graded_psi_std; std(Ideal(J,filter(a -> (degree(a) == i),psi_basis_J)))]
# end

for i = 1 : max_degree_B
	global graded_B = [graded_B ; Ideal(R,filter(a -> (degree(a) == i),B))]
	global graded_B_std = [graded_B_std; std(Ideal(R,filter(a -> (degree(a) == i),B)))]
end

for i = 1 : max_degree_B
	global graded_B_J = [graded_B_J ; Ideal(J,filter(a -> (degree(a) == i),B_J))]
	global graded_B_J_std = [graded_B_J_std; std(Ideal(J,filter(a -> (degree(a) == i),B_J)))]
end

global reduction_matrix = Array{typeof(zp(1))}(undef,size_Pn,len_B)
reduction_matrix .= ZERO
# global reduction_matrix = [[zp(0) for a = 1:len_B] for aa = 1:size_Pn]

cache_reduction = Dict()



function compute_reduction_matrix()
	for ind = 1:size_Pn
		g = Pn[ind]
		# psi = Pn[ind]
		# psi = change_ring(Pn[ind],J) #??
		# if psi == 1
			# if DEBUG
				# @assert(ind == 1)
				# assert psi_basis is ordered by increasing degree
			# end
		# else

			# deg = degree(psi)
			# ideal = graded_B[deg]
			# ideal_std = graded_B_std[deg]

			# rem = Singular.reduce(psi,ideal_std)
			# quo,r = Singular.lift(ideal,Ideal(R,psi-rem))
			# quo = quo[1]
			# s = size(gens(ideal),1)
			# quo =  [[[a[2],a[3]] for a in quo if a[1] == j] for j = 1 : s]
			# quo = [vector_to_polynomial(a) for a in quo]

			# ans = sum([change_ring(quo[j],J)*gens(ideal)[j] for j = 1 : s])
		ans = reduce_polynomial(g)

		ans_coeffs = collect(coefficients(ans))
		ans_mons = collect(monomials(ans))

		for j = 1:size(ans_coeffs,1)
			if ans_mons[j] in B
	    		mon_index = findall(xx->xx==ans_mons[j],B)[1]
	    		reduction_matrix[ind,mon_index] = ans_coeffs[j]
	    	end
    	end

	end
	# end


end

compute_reduction_matrix()
# reduction_matrix = transpose(reduce(hcat,reduction_matrix))

function to_Pn_basis(h)
	ans = [zp(0) for blah = 1: size_Pn]

	ans_coeffs = collect(coefficients(h))
	ans_mons = collect(monomials(h))

	for j = 1:size(ans_coeffs,1)
		mon_index = findall(xx->xx==ans_mons[j],Pn)[1]
		ans[mon_index] = ans_coeffs[j]
	end
	return transpose(ans)
end

function from_B_basis(v)
	return sum([B[j] * v[j] for j = 1:len_B])
end



frob_matrix = [[BigInt(0) for i = 1:len_B] for j = 1:len_B]


for i = 1:len_B
    fro = frobenius_on_cohom(i)
    println(B[i])
    mons = collect(monomials(fro))
    coeffs = R.(collect(coefficients(fro)))
    if DEBUG
    	fro_terms = collect(terms(fro))
    	@assert all([fro_terms[i] == mons[i]*coeffs[i] for i = 1:size(coeffs,1)])
    end
    for k = 1 : size(coeffs,1)

    	# could and should also cache this factorial thing

    	h = coeffs[k] * reduce_uv(mons[k])
    	# h 

    	h = from_B_basis(to_Pn_basis(h)*reduction_matrix)
    	if h != 0
    		denom = factorial(big(degree(mons[k]))-1)

	        ans_mons = collect(monomials(h))
	        ans_coeffs = 1
	        if zp == QQ
	        	ans_coeffs = BigInt.(numerator.(collect(coefficients(h)))).//BigInt.(denominator.(collect(coefficients(h))))
	        else
	    		ans_coeffs = BigInt.(collect(coefficients(h))).//1
	    	end
	    	ans_coeffs = [factorial(max(0,degree(ans_mons[j]))-1) //denom * ans_coeffs[j] for j=1:size(ans_coeffs,1)]


	    	if DEBUG
	    		ans_terms = collect(terms(h))
	    		# @assert all([ans_terms[i] == ans_coeffs[i]*ans_mons[i] for i = 1:size(ans_coeffs,1)])
	    	end

	    	# h = sum([factorial(big(degree(ans_mons[j]))-1)//factorial( big(degree(mons[k])-1)) * ans_mons[j]*ans_coeffs[j] for j=1:size(ans_coeffs,1)])

	    	# h = sum([factorial(max(0,degree(ans_mons[j]))-1) //denom * ans_mons[j]*ans_coeffs[j] for j=1:size(ans_coeffs,1)])
	    	if DEBUG
	    		ans_terms = collect(terms(h))
	    		# @assert all([ans_terms[i] == ans_coeffs[i]*ans_mons[i] for i = 1:size(ans_coeffs,1)])
	    	end
	    	# this is the right denom
	    	# ans_mons = collect(monomials(h))
	    	# ans_coeffs = collect(coefficients(h))

	    	# if zp == QQ
	    		# ans_coeffs //= denom
	    		ans_coeffs = BigInt.(numerator.(ans_coeffs)) .* BigInt.([invmod(BigInt(denominator(a)),big(p)^prec) for a in ans_coeffs])
	    		
	    	# else
	    		# ans_coeffs = BigInt.(ans_coeffs)
	    		# ans_coeffs //= denom
	    		# ans_coeffs = BigInt.(numerator.(ans_coeffs)) .* BigInt.([invmod(BigInt(denominator(a)),big(p)^prec) for a in ans_coeffs])
	    		
	    	# end

	    	# ans_coeffs = [zp(a) for a in ans_coeffs]

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
        
for i = 1:len_B
	for j = 1:len_B
		# if zp == QQ
			frob_matrix[i][j] =(mod(BigInt(numerator(frob_matrix[i][j])), p^prec)) * BigInt(invmod(BigInt(denominator(frob_matrix[i][j])),p^prec))
		# else
			# frob_matrix[i][j] =(mod(BigInt(numerator(frob_matrix[i][j])), p^prec)) * BigInt(invmod(BigInt(denominator(frob_matrix[i][j])),p^prec))
		# end
	end
end

println(frob_matrix)