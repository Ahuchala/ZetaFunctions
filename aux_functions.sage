
fdeg = fdegree
prod_rgens = prod(Rgens)
def degree(g):

	gdeg = sum([weights[i] * monomial_to_vector(g.monomials()[0])[i] for i in range(n)]) // fdeg
	# if DEBUG:
		# assert(all([sum([weights[i]*monomial_to_vector(a)[i] for i in range(n)])==gdeg*fdeg for a in g.monomials()]))
	return gdeg

def degree_vector(v):
	return sum([weights[i] * v[i] for i in range(n)]) // fdeg


# return the degree of g * xyz, etc
def affine_degree(g):
	return degree(prod_rgens*g)

