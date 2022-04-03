#ifndef GRAPH_H
#define GRAPH_H


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"


/* Find the Euclidean norm of two vectors of the same dim x 1 */
double euclidean_norm(dpoint_t v1, dpoint_t v2, size_t dim);


/* Calculate and return the weighted adjacency matrix of the given list of datapoints.
   
   In case of allocation failure, the return value has a `data` field of `NULL`. */
int graph_adjacent_matrix(dpoint_t input[], size_t num_data, size_t dim, matrix_t* output);


/* Calculate and return the diagonal degree matrix of the given matrix <mat>.
   In such case that the boolean is_sqrt equals <true>, the returned matrix is D ^ (-1/2).
   
   In case of allocation failure, the return value has a `data` field of `NULL`. */
int graph_diagonal_degree_matrix(matrix_t mat, bool is_sqrt, matrix_t* output);


/* Calculate the normalized graph Laplacian matrix of the given list of datapoints,
   and store the result in `output`. On success, returns 0.
   
   In case of any allocation failure, the return value is `BAD_ALLOC`.*/
int graph_normalized_laplacian(dpoint_t input[], size_t num_data, size_t dim, matrix_t *output);



#endif
