#ifndef GRAPH_H
#define GRAPH_H


#include <stdlib.h>
#include <math.h>
#include "matrix.h"

/* Define a structure that will hold a vector's coordinates. This will also be used int the spkmeans.c mechanism */
typedef struct {
    double *data;
} dpoint_t;


/* Find the Euclidean norm of two vectors of the same dim x 1 */
double euclidean_norm(dpoint_t v1, dpoint_t v2, size_t dim);


/* Calculate and return the adjacent matrix of the given matrix <mat> */
matrix_t graph_adjacent_matrix(dpoint_t input[], size_t dim);


/* Calculate and return the diagonal degree matrix of the given matrix <mat> */
matrix_t graph_diagonal_degree_matrix(matrix_t mat);


/* Calculate and return the normalized graph Laplacian matrix of the given matrix <mat> */
matrix_t graph_normalized_laplacian(dpoint_t input[], size_t dim);



#endif
