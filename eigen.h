#ifndef EIGEN_H
#define EIGEN_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"

/* Configuring convergence conditions to the Jacobi algorithm */
#define max_jacobi_iterations 100
#define epsilon pow(10, -15) 

/* Defining exception types */
#define ASYMMETRIC 3
#define NOT_REAL 4


/* Define a structre the will hold the indices in the matrix of a value.
 * Inconvient in terms of cache locality, therefore why we chose to have this only for a specific use in this file */
typedef struct matrix_ind {
	size_t i;
	size_t j;
} matrix_ind;


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
} jacobi_output;


/* Free a jacobi output - safely! */
void eigen_free_jacobi_safe(jacobi_output out);

/* Print the output of the jacobi algorithm */
void eigen_print_jacobi(jacobi_output out);

/* Return the output of the Jacobi algorithm when applied to the matrix <mat>
 * Pre-Conditions:
 * 		<mat> must be a "real" matrix
 * 		<mat> must be a "symmetric" matrix.
 * If sort equals <true>, then the first K eigen vectors and values will be returned. Else, they would be random (unless K is equal to the number of eigen values) */
jacobi_output eigen_jacobi(matrix_t mat, size_t K, bool sort);

/* Given the diagonal matrix <mat>, pull out its eigen values - sort them, determine K (the amount of eigen vectors we want), and form an eigen-vectors matrix.
 * This part might use the eigen heuristic gap if it was given an invalid K as an argument.
 * If sort equals <true>, then the first K eigen vectors and values will be returned. Else, they would be random (unless K is equal to the number of eigen values).
 * If K == 0: pick K by the heuristic gap.
 * If 0 < K <= mat_vectors.cols: pick the first K eigen values/vectors.
		If sort: return them sorted.
		else: return them unsorted (unpractical, doesn't make sense, and isn't reached by this program. But for the completency of this module, we kept this part)
 * If K > mat_vectors.cols: return all of the eigen values/vectors:
 		If sort: return them sorted.
 		else: return them unsorted. */
jacobi_output eigen_format_eigen_vectors(matrix_t mat_vectors, matrix_t mat_eigens, size_t K, bool sort);

/* In case K wasn't given as an input, then try to determine it using the eigen heuristic gap.
 * Pre-Conditions:
 *		the given array of eigen values must be sorted by vaule (remember eigen is a struct) */
size_t eigen_heuristic_gap(eigen* sorted_eigen_values, size_t rows);

/* Define a "compare" functino between two "eigen"-s. */
int compare(const void* eigen1, const void* eigen2);

/* Given an already diagonal matrix, extrat its eigen values.
 * If sort equals <true>, sort the eigen values.. */
eigen* eigen_extract_eigen_values(matrix_t mat, bool sort);

/* In the jacobi algorithm, this is the function that transforms A_tag (the next matrix in the recurssive algorithm), through the current A matrix */
void eigen_update_jacobi_A_tag(matrix_t A_tag, matrix_t A);

/* Given a symmetric matrix called <mat>, find the competent rotation matrix for it */
matrix_t eigen_build_rotation_matrix(matrix_t mat);

/* Given two matrices, determine the distance between their sum of squared off-diagonals */
double eigen_distance_of_squared_off(matrix_t mat1, matrix_t mat2);

/* Given a matrix, return the sum of squared off-diagonals. */
double eigen_sum_squared_off(matrix_t mat);

/* Return the index of an off diagonal inside the one matrix's array that holds the largest absolute value */  
matrix_ind eigen_ind_of_largest_offdiagonal(matrix_t mat); 

/* Calculae the values of 'c' and 's' of the competent rotation matrix */
void eigen_calc_c_s(double* c, double *s, matrix_t mat, matrix_ind loc);


/* Return the sign of a value.
 * (val >= 0) <-> ret == 1 */
int sign(double val);



#endif
