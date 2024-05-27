# DTYPE long long



#compile with:     %runfile mat_mul.pyx

# cdef extern from "mat_mul.c":
#     int** c_mat_mul(A,B)

# cdef extern from "mat_mul.c":
#     int* sum_of_two_vectors(int n, int * vec1, int * vec2)

cdef extern from "mat_mul.c":
    int* c_mat_vec_mul(int n, int * A, int * v)
    long long* c_iterated_mat_vec_mul(int n, int r, long long * A, long long * B, long long * v, long long p_prec)



from libc.stdlib cimport malloc, free

# def my_bridge_function(A,B):
#     return c_mat_mul(A,B) # This is the C function from mat_mul.c

# def sage_sum_of_vectors(n, list1, list2):
#     cdef int * vec1 = <int *> malloc(n*sizeof(int))
#     cdef int * vec2 = <int *> malloc(n*sizeof(int))

#     # Fill the vectors
#     for i in range(n):
#         vec1[i] = list1[i]
#         vec2[i] = list2[i]

#     # Call the C function
#     cdef int * vec3 = sum_of_two_vectors(n,vec1,vec2)

#     # Save the answer in a Python object
#     answer = [vec3[i] for i in range(n)]

#     free(vec1)
#     free(vec2)
#     free(vec3)

#     return answer

# important: A is passed into c as a vector
def sage_mat_vec_mul(n, matrix_A, vec_v):
    cdef int * A = <int *> malloc(n*n*sizeof(int))
    cdef int * v = <int *> malloc(n*sizeof(int))

    # Fill the vectors
    for i in range(n):
        for j in range(n):
            A[i+n*j] = matrix_A[i][j]
        v[i] = vec_v[i]

    # Call the C function
    cdef int * vec3 = c_mat_vec_mul(n,A,v)

    # Save the answer in a Python object
    answer = [vec3[i] for i in range(n)]

    free(A)
    free(v)
    free(vec3)

    return answer

# returns prod_{k=1}^{r+1}v(A+jB) mod p^prec
def sage_iterated_mat_vec_mul(n, r, matrix_A, matrix_B, vec_v, p_prec):
    cdef long long * A = <long long *> malloc(n*n*sizeof(long long))
    cdef long long * B = <long long *> malloc(n*n*sizeof(long long))
    cdef long long * v = <long long *> malloc(n*sizeof(long long))

    # Fill the vectors
    # note the transposes
    for i in range(n):
        for j in range(n):
            A[i+n*j] = matrix_A[j][i]
            B[i+n*j] = matrix_B[j][i]
        v[i] = vec_v[i]

    # Call the C function
    cdef long long * vec3 = c_iterated_mat_vec_mul(n,r,A,B,v, p_prec)

    # Save the answer in a Python object
    answer = [vec3[i] for i in range(n)]

    free(A)
    free(B)
    free(v)
    free(vec3)

    return answer