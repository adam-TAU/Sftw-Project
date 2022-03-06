#ifndef GRAPH_H
#define GRAPH_H


#include <stdlib.h>
#include <math.h>
#include "matrix.h"


/* Find the Euclidean norm of two vectors of the same dim x 1 */
double graph_euclidena_norm(double v1[], double v2[], int dim);


/* Create an Identity Matrix, with the dimensions of dim x dim */
matrix_t graph_identity_matrix(size_t dim);


/* Calculate and return the adjacent matrix of the given matrix <mat> */
matrix_t graph_adjacent_matrix(double input[][], int dim);


/* Calculate and return the diagonal degree matrix of the given matrix <mat> */
matrix_t graph_diagonal_degree_matrix(matrix_t mat);


/* Calculate and return the normalized graph Laplacian matrix of the given matrix <mat> */
matrix_t graph_normalized_laplacian(double input[][], int dim);



#endif
