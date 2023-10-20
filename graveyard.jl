
I_f = Ideal(R,[derivative(f,Rgens[i]) for i = 1:n])
I_f_std = std(I_f)

shifted_B = [divide_by_x1xn(b) for b in B]
reduction_dict = Dict()

function reduce_to_basis(g)
	summer = 0

# maybe could make this a set
	term_list = collect(terms(g))
	while !isempty(term_list)
		term = pop!(term_list)
		# println(term)
		coeff = collect(coefficients(term))[1]
		mon = collect(monomials(term))[1]
		if !(mon in shifted_B)
			if !haskey(reduction_dict,mon)
				# reduce it
				rem = Singular.reduce(mon, I_f_std)  # poly, ideal
				quo = mon - rem
				quo = Singular.lift(I_f,Ideal(R,quo))[1][1] #Ideal, subideal
				quo =  [[[a[2],a[3]] for a in quo if a[1] == i] for i = 1 : n]
				quo = [vector_to_polynomial(a) for a in quo]

				reduction_dict[mon] = rem + sum([derivative(quo[i],Rgens[i]) for i = 1:n])

			end
			result = coeff * reduction_dict[mon]
			term_list = [term_list; collect(terms(result))]
		else
			summer += term
		end
	end
	return summer
end
