# S = GF(p^(prec+2))

def mat_mul(A,B,k,g_vec):
	# maybe escape if j is small

	# g_vec_S = vector(S,g_vec)
	# A_S = matrix(S,A)
	# B_S = matrix(S,B)
	# A_np = np.matrix(A,dtype=np.int32)
	# B_np = np.matrix(B,dtype=np.int32)

	for j in range(1,k+1):
		# g_vec_S *= A_S - j*B_S
		g_vec *= A - j *B
		# todo: only do remainder when needed, and find the correct precision bound
		# g_vec = np.remainder(g_vec,p^(prec+2))
	return g_vec
	# print(g_vec_S)
	# return [ZZ(_) for _ in g_vec_S]

# import cython
# cimport numpy as np
# cimport cython

# # returns \prod_{i=1}^k (A+jB) g_vec
# def np_mat_mul(A,B,k,g_vec):
# 	# maybe escape if j is small

# 	A_np = np.matrix(A,dtype=np.int32)
# 	B_np = np.matrix(B,dtype=np.int32)

# 	for j in range(1,k+1):
# 		g_vec *= A_np-j*B_np
# 		# todo: only do remainder when needed, and find the correct precision bound
# 		g_vec = np.remainder(g_vec,p^(prec+2))
# 	return g_vec

# def mat_mul(A,B,k,g_vec):
# 	# maybe escape if j is small

# 	A_c = initialize_cython_matrix(A)
# 	B_c = initialize_cython_matrix(B)

# 	return cmat_mul(A_c,B_c,k,g_vec)
# 	# for j in range(1,k+1):
# 		# g_vec *= A_np-j*B_np
# 		# todo: only do remainder when needed, and find the correct precision bound
# 		# g_vec = np.remainder(g_vec,p^(prec+2))
# 	# return g_vec

# cdef np.ndarray[np.int32, ndim=2] cmat_mul(np.ndarray[np.int32_t, ndim=2] A, np.ndarray[np.int32, ndim=2] B, int k, np.ndarray[np.int32, ndim=1] g_vec):
# 	# maybe escape if j is small

# 	for j in range(1,k+1):
# 		g_vec *= A-j*B
# 		# todo: only do remainder when needed, and find the correct precision bound
# 		g_vec = np.remainder(g_vec,p^(prec+2))
# 	return g_vec


# def initialize_cython_matrix(A):
# 	# given a sage matrix A, return a cython/numpy array A_c
# 	A_c = np.matrix(A,dtype=np.int32)
# 	return A_c