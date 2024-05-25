load("mat_mul.sage")
A = matrix([(0, 2, 1, 1, -1, 0, 34), (1, 0, 6, 1, -1, -2, -1), (1, 1, 0, -1, -3, 23, 15), (5, 4, 0, 1, -6, 4, 2), (1, -20, -1, -1, 1, 0, 3), (61, 1, 3, 2, -1, -1, -5), (-2, -3, -1, -2, 0, 0, 2)])
v = vector([3, 29, -1, 1, 0, 0, 2])
assert multiply_matrix_vector(A,v) == vector((126, -4, 61, 136, -571, 201, -90))

A = random_matrix(ZZ,7,7)
v = vector(random_matrix(ZZ,7,1))
assert multiply_matrix_vector(A,v) == A*v