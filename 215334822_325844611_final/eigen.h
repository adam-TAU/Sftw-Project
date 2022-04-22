#ifndef EIGEN_H
#define EIGEN_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"

/* Configuring convergence conditions to the Jacobi algorithm */
#define max_jacobi_iterations 100 /* must be 100 */
#define epsilon 1e-5


/* Define a structure that will hold an eigen value's <value> and <ind at matrix> */
typedef struct eigen_value_t {
	double value;
	size_t col;
} eigen_t;

/* Define a structure that will hold the output of the Jacobi algorithm. That includes the eigen values as well as the eigen_vectors matrix
 * If K is desired in the output, you can just check that throught the amount of cols of eigen_vectors.
 * K == eigen_vectors.cols */
typedef struct jacobi_t {
	matrix_t eigen_vectors;
	eigen_t* eigen_values;
	int signal;
} jacobi_t;


/* Return the output of the Jacobi algorithm when applied to the matrix <mat>.
 * Pre-Conditions:
 * 		<mat> must be a "real" matrix
 * 		<mat> must be a "symmetric" matrix.
 *		<K> must be lower or equal to <mat.rows>
 *
 * If K == 0: using the heuristic gap, determine a new positive K, and return the first (while sorted) K eigen values along with their eigen vectors
 * If 0 < K < mat.rows: return the first (while sorted) K eigen values along with thier eigen vectors
 * If K = mat.rows: return all of the K eigen values along with their eigen vectors. This part doesn't necessarily return the eigen values sorted */
int eigen_jacobi(matrix_t mat, size_t K, jacobi_t* output);

/* Cast a variable of type <jacobi_output> into <matrix_t>, for printing purposes only!
 * If we print that matrix, we will get the desired printage of a jacobi output.
 * The outputted matrix is stored into the <output> argument.
 * This conversion frees the contents of <origin>. */
int eigen_jacobi_to_mat(jacobi_t origin, matrix_t *output);

/* Define a "compare" function between two "eigen"-s. */
int eigen_compare(const void* eigen1, const void* eigen2);

/* Return the sign of a value.
 * (val >= 0) <=> ret == 1 */
int sign(double val);





#endif
