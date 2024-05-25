int main() {
    int matA[2][2]={0,1,2,3};
    int matB[2][2]={0,1,2,3};
    int matC[2][2];
    int i, j, k;
    for(i=0;i<M;i++){
        for(j=0;j<K;j++){
            matC[i][j]=0;
            for(k=0;k<N;k++){
                matC[i][j]+=matA[i][k]*matB[k][j];
            }
        }
    }
}