#ifndef SPKMEANS_GOALS_H
#define SPKMEANS_GOALS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "graph.h"
#include "eigen.h"


/* A function to print the weighted adjacency matrix out of the given vectors */
void print_weighted_adjacency_matrix(void);

/* A function to print the diagonal degree matrix of the given vectors */
void print_diagonal_degree_matrix(void);

/* A function to print the normalized graph laplacian matrix of the given vectors */
void print_normalized_laplacian(void);


/* A function to print the eigen values and eigen vectors of the given input. 
 * Pre-Conditions:
 * 		<mat> must be a "real" matrix
 * 		<mat> must be a "symmetric" matrix */
void print_jacobi_output(void);

/* A function to perform the whole algorithm of the the spectral clustering algorithm. It prints out the centroids induced from that algorithm */
matrix_t get_T_of_spectral_kmeans(size_t K);






#endif
