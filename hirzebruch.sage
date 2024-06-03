from sage.rings.polynomial.polydict import ETuple


# based on code from https://pbelmans.ncag.info/blog/2018/01/27/hodge-numbers-of-complete-intersections/


def h_pq(d,n):
	r = len(d)
	Poly_ring.<y,z> = PowerSeriesRing(ZZ,n+1-r)
	xi = 1/((1-z)*(1-z*y)) * prod(((1+z*y)^d[i] - (1-z)^d[i])/((1+z*y)^d[i]+y*(1-z)^d[i]) for i in range(r))/z^(r)

	# f = 1/((1+x)*(1+y))*prod(-1+((1+x)^d[i]-(1+y)^d[i])/(x*(1+y)^d[i]-y*(1+x)^d[i]) for i in range(len_d)) + 1/(1-x*y)
	f_dict = xi.dict()
	ls = [f_dict[ETuple([i,n-r-i])]  if ETuple([i,n-r-i]) in f_dict else 0 for i in range(n + 1-r)]
	for i in range(n+1-r):
		if 2*i != n-r:
			ls[i] = (-1)^(n-r-i)*(ls[i] -(-1)^i)#??
		else:
			ls[i] *= (-1)^i
	return ls

# in: a list like [1,20,1]
# out: the hodge polygon slopes including multiplicity, e.g. [0, 1,1,1,1...,1,2]
def h_pq_to_hodge_polygon(ls):
	return_ls = []
	slope = 0
	for i in range(len(ls)):
		for j in range(ls[i]):
			return_ls.append(slope)
		slope += 1
	return return_ls