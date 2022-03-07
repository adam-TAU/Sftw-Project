#ifndef EIGEN_H
#define EIGEN_H

#include <stdlib.h>
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


/* Return the output of the Jacobi algorithm when applied to the matrix <mat>
 * Pre-Conditions:
 *		<mat> must be a symmetric matrix
 * 		<mat> must be a "real" matrix */
matrix_t eigen_jacobi(matrix_t mat);

/* Given the diagonal matrix <mat>, pull out its eigen values, and determine the index of the most distanced gap between the ordered eigen values */
size_t eigen_heuristic_gap(matrix_t mat);

/* Given an already diagonal matrix, extrat its eigen values and sort them by order. */
eigen* eigen_extract_eigen_values(matrix_t mat);

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


/* Define a "compare" functino between two "eigen"-s. */
int compare(const void* eigen1, const void* eigen2);


#endif
