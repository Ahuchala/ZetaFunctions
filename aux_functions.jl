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

# return dg/dxj
function my_derivative(g,j)
	ring_gens = gens(parent(g))
	local mons = collect(monomials(g))
	coeffs = collect(coefficients(g))
	ans = parent(g)(0)
	for i = 1:length(coeffs)
		u = monomial_to_vector(mons[i])
		# for j = 1:n
			if u[j] > 0
				u[j] -= 1
				ans += coeffs[i] * (u[j]+1)*prod(ring_gens.^u)
				u[j] += 1
			end
		# end
	end
	return ans
end


# returns sum xi dg/dxi, using sum xi d(x^u)/dxi = sum(u) x^u
function toric_derivative(g)
	return sum([sum(collect(exponent_vectors(monomial))[1]) * monomial for monomial in monomials(g)])
end

# same as toric_derviative, but input is of the form [[[expvector], coeff], [[expvector], coeff],...]
function toric_derivative_vector(v)
	return [[a[1], sum(a[1])*a[2]] for a in v]
end