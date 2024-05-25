// int** c_mat_mul(int** a, int** b) {
//     size_t dim_a = sizeof(a) / sizeof(a[0]);
//     int** c = {0}; // sets to zero
//     int i, j, k;
//     for(i=0;i<dim_a;i++){
//         for(j=0;j<dim_a;j++){
//             for(k=0;k<dim_a;k++){
//                 c[i][j]+=a[i][k]*b[k][j];
//             }
//         }
//     }
//     return c;
// }

// int * sum_of_two_vectors(int n, int * vec1, int * vec2){
//   /*
//      INPUT : two arrays vec1,vec2 of n integers
//      OUTPUT: an array of size n equal to vec1+vec2
//   */
//   int * sum = (int *) malloc(n*sizeof(int));
//   int i;

//   for(i=0;i<n;i++)
//     sum[i] = vec1[i] + vec2[i];
//   return sum;
// }

int * c_mat_vec_mul(int n, int * A, int * v){
  /*
     INPUT : a vectorized matrix A and a vector v of integers of dimension n
     OUTPUT: a vector A*v
  */

    int * return_vec = (int *) malloc(n*sizeof(int));
    int i,j;
    for(i = 0; i < n; i++) {
        return_vec[i] = 0;
    }

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            return_vec[i] += A[i+n*j] * v[j];
            }
        }
    return return_vec;
}

// returns prod_{k=1}^{r}(A+jB)v
long * c_iterated_mat_vec_mul(int n, int r, long * A, long * B, long * v, long p_prec){
//  n = dim(A)

    long * old_vec = (long *) malloc(n*sizeof(long));
    long * new_vec = (long *) malloc(n*sizeof(long));
    int i,j, k;
    for(i = 0; i < n; i++) {
        old_vec[i] = v[i];
        new_vec[i] = 0;
    }

    for (k = 1; k < r+1; k++) {
        for(i = 0; i < n; i++) {
            for(j = 0; j < n; j++) {
                new_vec[i] += A[i+n*j] * old_vec[j] - k*B[i+n*j] * old_vec[j];
            }
        }
        for(i = 0; i < n; i++) {
            old_vec[i] = new_vec[i] % p_prec;
            new_vec[i] = 0;
        }
    }
    free(new_vec);
    return old_vec;
}