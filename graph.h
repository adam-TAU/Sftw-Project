#ifndef GRAPH_H
#define GRAPH_H


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"


/* Find the Euclidean norm of two vectors of the same dim x 1 */
double euclidean_norm(dpoint_t v1, dpoint_t v2, size_t dim);


/* Calculate and return the adjacent matrix of the given matrix <mat> */
matrix_t graph_adjacent_matrix(dpoint_t input[], size_t num_data, size_t dim);


/* Calculate and return the diagonal degree matrix of the given matrix <mat>.
 * In such case that the boolean is_sqrt equals <true>, the returned matrix is D ^ (-1/2) */
matrix_t graph_diagonal_degree_matrix(matrix_t mat, bool is_sqrt);


/* Calculate and return the normalized graph Laplacian matrix of the given matrix <mat> */
matrix_t graph_normalized_laplacian(dpoint_t input[], size_t num_data, size_t dim);



#endif
