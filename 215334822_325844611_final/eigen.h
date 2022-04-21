#ifndef EIGEN_H
#define EIGEN_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"

/* Configuring convergence conditions to the Jacobi algorithm */
#define max_jacobi_iterations 100 /* must be 100 */
#define epsilon 1e-5

/* Defining exception types */
#define ASYMMETRIC 3
#define NOT_REAL 4
#define HEURISTIC_PICKED_1 5


/* Define a structure that will hold an eigen value's <value> and <ind at matrix> */
typedef struct eigen_value {
	double value;
	size_t col;
} eigen;

/* Define a structure that will hold the output of the Jacobi algorithm. That includes the eigen values as well as the eigen_vectors matrix
 * If K is desired in the output, you can just check that throught the amount of cols of K_eigen_vectors.
 * K == K_eigen_vectors.cols */
typedef struct jacobi_output {
	matrix_t K_eigen_vectors;
	eigen* eigen_values;
	int signal;
} jacobi_output;


/* Free a jacobi output - safely! */
void eigen_free_jacobi_safe(jacobi_output *out);

/* Print the output of the jacobi algorithm */
void eigen_print_jacobi(jacobi_output out);

/* Return the output of the Jacobi algorithm when applied to the matrix <mat>.
 * Pre-Conditions:
 * 		<mat> must be a "real" matrix
 * 		<mat> must be a "symmetric" matrix.
 *		<K> must be lower or equal to <mat.rows>
 *
 * If K == 0: using the heuristic gap, determine a new positive K, and return the first (while sorted) K eigen values along with their eigen vectors
 * If 0 < K < mat.rows: return the first (while sorted) K eigen values along with thier eigen vectors
 * If K = mat.rows: return all of the K eigen values along with their eigen vectors. This part doesn't necessarily return the eigen values sorted */
int eigen_jacobi(matrix_t mat, size_t K, jacobi_output* output);

/* Define a "compare" function between two "eigen"-s. */
int eigen_compare(const void* eigen1, const void* eigen2);

/* Return the sign of a value.
 * (val >= 0) <=> ret == 1 */
int sign(double val);





#endif
