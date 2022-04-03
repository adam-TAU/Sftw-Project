#ifndef SPKMEANS_GOALS_H
#define SPKMEANS_GOALS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "graph.h"
#include "eigen.h"


extern void assert_other(bool condition);

/* A function to print the weighted adjacency matrix out of the given vectors */
int print_weighted_adjacency_matrix(void);

/* A function to print the diagonal degree matrix of the given vectors */
int print_diagonal_degree_matrix(void);

/* A function to print the normalized graph laplacian matrix of the given vectors */
int print_normalized_laplacian(void);


/* A function to print the eigen values and eigen vectors of the given input. 
 * Pre-Conditions:
 * 		<mat> must be a "real" matrix
 * 		<mat> must be a "symmetric" matrix */
int print_jacobi_output(void);

/* A function to perform the whole algorithm of the spectral clustering algorithm. It returns the spectral datapoints for kmeans++ as a matrix */
int get_T_of_spectral_kmeans(size_t K, matrix_t *output);






#endif
