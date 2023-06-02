
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