#ifndef MATRIX_H
#define MATRIX_H
#include <stdlib.h>

#define DIM_MISMATCH 1
#define BAD_ALLOC 2

typedef struct matrix {
	double *data;
	size_t rows;
	size_t cols;
	size_t len;
} matrix_t;


/* Define a structre the will hold the indices in the matrix of a value.
 * Inconvient in terms of cache locality, therefore why we chose to have this only for a specific use in this file */
typedef struct matrix_ind {
	size_t i;
	size_t j;
} matrix_ind;

/* Define a structure that will hold a vector's coordinates. This will also be used in the spkmeans.c mechanism */
typedef struct {
	double *data;
    size_t current_set;
} dpoint_t;


#define bool int
#define true 1
#define false 0

/*
 * BASIC METHODS
 */

/* Creates a new matrix with the given dimensions. The created matrix must be
   freed with `matrix_free`. The created matrix is zero-initialized.

   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
int matrix_new(size_t rows, size_t cols, matrix_t* output);

/* Clones the given matrix into a newly allocated matrix.

   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
int matrix_clone(matrix_t mat, matrix_t* output);

/* Copies data from one matrix into another.
   If the matrices are different in size, no change is made and
   DIM_MISMATCH is returned. */
int matrix_copy(matrix_t dest, matrix_t src);

/* Builds a matrix out of an already existing dataset of <dpoint_t>s.

   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
int matrix_build_from_dpoints(dpoint_t* vectors, size_t num_vectors, size_t dim, matrix_t* output);

/* Gets the desired element from the given matrix. */
double matrix_get(matrix_t mat, size_t i, size_t j);

/* Set the desired element from the given matrix. */
void matrix_set(matrix_t mat, size_t i, size_t j, double val);

/* Prints the given matrix's rows. */
void matrix_print_rows(matrix_t mat);

/* Frees a given matrix that was allocated using any matrix method. */
void matrix_free(matrix_t mat);

/* Calls matrix_free on `mat`, but only if `mat.data` isn't `NULL`. */
void matrix_free_safe(matrix_t mat);

/* Create an Identity Matrix, with the dimensions of dim x dim.

   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
int matrix_identity(size_t dim, matrix_t* output);

/* Given a matrix, return the sum of squared off-diagonals. */
double matrix_sum_squared_off(matrix_t mat);

/* Return the index of an off diagonal inside the one matrix's array that holds the largest absolute value.
 * Pre-condition: <sym_mat> must be a symmetric matrix!
 * Notice: due to the pre-condition, only the higher half of the matrix will be scanned for the largest off diagonal value.
 * Therefore, the index returned will and must be out of the higher half of the matrix. */  
matrix_ind matrix_ind_of_largest_offdiagonal(matrix_t sym_mat);


#endif /* MATRIX_H */
