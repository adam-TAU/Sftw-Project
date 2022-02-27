#include <stdio.h>
void *calloc(size_t nitems, size_t size);

int numOfVectors = 5;

void printMatrix(double** mat, int numOfRows, int numOfCols) {
    /*prints a matrix*/
    int i, j;
    for (i = 0; i < numOfRows; i++) {
        for (j = 0; j < numOfCols; j++) {
            if ((mat[i][j]<0)&&(mat[i][j]>-0.00005)){
                mat[i][j] = 0;
            }
            printf("%.4f", mat[i][j]); /*format the floats precision to 4 digits*/
            if (j < numOfCols - 1) {
                printf(",");
            }
        }
        if (i < numOfRows - 1) {
            printf("\n");
        }
    }
}

void printJacobi(double **A, double **V) {
    /*gets A matrix (for eigenvalues) and V matrix (for eigenvectors) 
    and prints them according to instructions*/
    int i,j;
    for (i = 0; i < numOfVectors; i++) {
        if ((A[i][i]<0)&&(A[i][i]>-0.00005)){
                A[i][i] = 0;
        }
        printf("%.4f", A[i][i]); /*eigenvalues, Format to 4 digits*/
            if (i < numOfVectors - 1) {
                printf(",");
            }
    }
    printf("\n");
    for (i = 0; i < numOfVectors; i++) {
        for (j = 0; j < numOfVectors; j++) {
            if ((V[j][i]<0)&&(V[j][i]>-0.00005)){
                V[j][i] = 0;
            }
            printf("%.4f", V[j][i]); /*Transpose V, Format to 4 digits*/
            if (j < numOfVectors - 1) {
                printf(",");
            }
        }
        if ( i < numOfVectors - 1) {
            printf("\n");
        }
    }
}


int main() {
    int i,j;
    double **mat1, **mat2;
    mat1 = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    mat2 = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    for (i = 0; i < numOfVectors; i++) {
        mat1[i] = (double *)calloc(numOfVectors, sizeof(double));
        mat2[i] = (double *)calloc(numOfVectors, sizeof(double));
    }

    for (i=0;i<numOfVectors;i++){
        for (j=0;j<numOfVectors;j++){
            if (i==0)
                mat1[i][j] = mat2[i][j] = 0;
            else if (i==1)
                mat1[i][j] = mat2[i][j] = 0.00005;
            else if (i==2)
                mat1[i][j] = mat2[i][j] = -0.00005;
            else if (i==3)
                mat1[i][j] = mat2[i][j] = -0.000049;
            else if (i==4)
                mat1[i][j] = mat2[i][j] = -0.000051;
        }
    }

    /*printMatrix(mat1,numOfVectors,numOfVectors);*/
    printJacobi(mat1,mat2);
    return 0;
}
