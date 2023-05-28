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
# function toric_derivative(g)
# 	return sum([sum(collect(exponent_vectors(monomial))[1]) * monomial for monomial in monomials(g)])
# end

# same as toric_derviative, but input is of the form g1,...,gn
# and returns sum xi dg_i dxi
function toric_derivative_vector(g)
	coeff_vec_g = [collect(coefficients(g[i])) for i = 1:n]
	mon_vec_g = [collect(monomials(g[i])) for i = 1:n]
	exp_vec_g = [collect(exponent_vectors(g[i])) for i = 1:n]
	ans = 0
	for i = 1:n
		num_coeffs = size(coeff_vec_g[i],1)
		for j = 1:num_coeffs
			ans += coeff_vec_g[i][j]*mon_vec_g[i][j]*exp_vec_g[i][j][i]
		end
	end
	return ans
end