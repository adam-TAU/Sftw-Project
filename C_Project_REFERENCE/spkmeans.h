#include <stdio.h>
#ifndef SPKMEANS_H_
#define SPKMEANS_H_

typedef struct eigenVector {
    double eigenVal;
    int columnIndex;
} eigenVector;  

int k, dimension, numOfVectors, changes, max_iter;
float rawK, rawMaxIter;
double *eigenVals, *eigenGaps;
double **vectors, **centroids, **wam, **ddg, **lnorm, **V, **U;
int **clusters, *clustersSizes;
char *goal;
eigenVector *eigenVectors;

void errorAssert(int cond, int isInputError);
int calcDimension(char buffer[]);
void readFile(FILE *file);
void assignUToVectors(void); 
void initCentroids(void); 
double distance(double *vector1, double *vector2);
int closestCentroid(double *vector);
void assignVectorToCluster(void); 
double* calcCentroidForCluster(int clusterInd);
void updateCentroidValue(void);
void deepClone(double **a, double** b);
void printMatrix(double** mat, int numOfRows, int numOfCols); 
double** matrixMultiplication(double** a, double** b);
void squareMatrixTranspose(double **matrix, int numOfRows);
double calcWeightsForAdjacencyMatrix(double *vector1, double *vector2);
double** weightedAdjacencyMatrix(void);
double** diagonalDegreeMatrix(int calcWam, int toPrint);
double** laplacianNorm(void);
int* maxOffDiagonalValue(double** mat);
double calcTheta(double **matrix, int i, int j);
double calcT(double theta);
double calcC(double t);
double calcS(double t, double c);
void createRotationMatrixP(double** P, int maxRow, int maxCol, double c, double s);
void updateAPrime(double** A, double** APrime, int i, int j, double c, double s);
double calcOffSquared(double** mat);
int checkConvergence(double** A, double** APrime);
void printJacobi(double **A, double **V); 
double** jacobi(double **A, int toPrint);
int compareEigenVectors(const void *a, const void *b); 
void sortEigenVectorsAndValues(void); 
int eigengapHeuristic(void);
void normalizeUMatrix(void); 
void createUMatrix(void);
void free2DDoubleArray(double ** arr, int numOfElements);
void free2DIntArray(int ** arr, int numOfElements);
void freeMemory(void);

#endif