# %runfile mat_mul.pyx


if USE_CYTHON:
	load("mat_mul.pyx")

# returns A*v
def multiply_matrix_vector(A,v):
    n = len(A[0])
    return vector(sage_mat_vec_mul(n,A,v))

def mat_mul(A,B,k,g_vec):
	if all(vec == 0 for vec in g_vec): #g_vec is zero
		return g_vec

	# maybe escape if j is small
	if USE_CYTHON and k > 5:
		return vector(sage_iterated_mat_vec_mul(len(A[0]),k,A,B,g_vec,p^(prec+arithmetic_precision_increase)))
	
	for j in range(1,k+1):
		g_vec *= A - j *B
	return g_vec

	# for j in range(1,k+1):
	# 	# g_vec_S *= A_S - j*B_S
	# 	if USE_CYTHON:
	# 		g_vec = multiply_matrix_vector((A-j*B).transpose(),g_vec)
	# 		print(g_vec)
	# 		g_vec = [_ % p^(prec+2) for  _ in g_vec]
	# 	else:
	# 		g_vec *= A - j *B
	# 	# todo: only do remainder when needed, and find the correct precision bound
	# 	# g_vec = np.remainder(g_vec,p^(prec+2))
	# return g_vec

