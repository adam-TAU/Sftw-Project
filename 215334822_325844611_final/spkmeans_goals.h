#ifndef SPKMEANS_GOALS_H
#define SPKMEANS_GOALS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "graph.h"
#include "eigen.h"


extern void assert_other(bool condition);

/* A function to print the weighted adjacency matrix out of the given vectors. Store the output into <output> */
int build_weighted_adjacency_matrix(matrix_t *output);

/* A function to print the diagonal degree matrix of the given vectors. Store the output into <output> */
int build_diagonal_degree_matrix(matrix_t *output);

/* A function to print the normalized graph laplacian matrix of the given vectors. Store the output into <output> */
int build_normalized_laplacian(matrix_t *output);


/* A function to print the eigen values and eigen vectors of the given input. Store the output into <output> */
int build_jacobi_output(matrix_t *output);

/* A function to perform the whole algorithm of the spectral clustering algorithm. It returns the spectral datapoints for kmeans++ as a matrix.
 * Store the output into <output> */
int build_T_of_spectral_kmeans(size_t K, matrix_t *output);






#endif
