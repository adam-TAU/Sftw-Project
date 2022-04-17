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

/* Swaps the data contained in the given matrices, efficiently (by swapping pointers). */
void matrix_swap(matrix_t *mat1, matrix_t *mat2);

/* Copies data from one matrix into another.
   If the matrices are different in size, no change is made and
   DIM_MISMATCH is returned. */
int matrix_copy(matrix_t dest, matrix_t src);

/* Builds a matrix out of an already existing dataset of <dpoint_t>s.

   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
int matrix_build_from_dpoints(dpoint_t* vectors, size_t num_vectors, size_t dim, matrix_t* output);

/* Calculates the index of the desired element for use with the matrix's
   inner `data` field. This is mainly used for optimization.

   This function also performs a bounds check on the given values: in case
   of failure, an error is printed to the console and the returned value
   is the maximum for `size_t`. */
size_t matrix_calc_index(matrix_t mat, size_t i, size_t j);

/* Gets the desired element from the given matrix.

   This function calls `matrix_calc_index`. */
double matrix_get(matrix_t mat, size_t i, size_t j);

/* Set the desired element from the given matrix.

   This function calls `matrix_calc_index`. */
void matrix_set(matrix_t mat, size_t i, size_t j, double val);

/* Pretty-prints the given matrix. */
void matrix_petty_print(matrix_t mat);

/* Prints the given matrix's rows. */
void matrix_print_rows(matrix_t mat);

/* Prints the given matrix's columns. */
void matrix_print_cols(matrix_t mat);

/* Frees a given matrix that was allocated using any matrix method. */
void matrix_free(matrix_t mat);

/* Calls matrix_free on `mat`, but only if `mat.data` isn't `NULL`. */
void matrix_free_safe(matrix_t mat);

/*
 * BASIC ARITHMETIC METHODS
 */

/* Adds the second matrix into the first. The return value is `0` only if
   the dimensions match correctly, and `DIM_MISMATCH` otherwise.

   In case of failure, the first matrix isn't modified. */
int matrix_add_assign(matrix_t self, matrix_t other);

/* Adds the given matrices and stores the newly allocated matrix in `output`.
   The return value is `0` if and only if the operation was successful.
   The output matrix must be unrelated to the input matrices.

Errors:

In case of dimension mismatch, the return value is `DIM_MISMATCH`.
In case of allocation failure, the return value is `BAD_ALLOC`. */
int matrix_add(matrix_t mat1, matrix_t mat2, matrix_t *output);

/* Multiplies the given scalar into the given matrix. */
void matrix_mul_scalar_assign(matrix_t mat, double scalar);

/* Multiplies the given scalar and matrix, and returns the allocated result.

   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
int matrix_mul_scalar(matrix_t mat, double scalar, matrix_t* output);

/*
 * MATRIX MULTIPLICATION
 */

/* Multiplies the given matrices to produce a new matrix, stored in `output`.
Precondition: `mat1.cols == mat2.rows`.
The output matrix must be unrelated to the input matrices.

Errors:

In case of dimension mismatch, the return value is `DIM_MISMATCH`.
In case of allocation failure, the return value is `BAD_ALLOC`. */
int matrix_mul(matrix_t mat1, matrix_t mat2, matrix_t *output);


/* Calculates `mat1` * `mat2`, and stores the result in `output`.
   Precondition: `mat1.cols == mat2.rows`, `mat1.rows == output.rows`, and `mat2.cols == output.rows`.
   Additionally, `output` MUST point to a different memory location than either `mat1` or `mat2`.

   This function performs no allocations - therefore, the only possible error value is `DIM_MISMATCH`. */
int matrix_mul_buffer(matrix_t mat1, matrix_t mat2, matrix_t output);


/* Create an Identity Matrix, with the dimensions of dim x dim.

   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
int matrix_identity(size_t dim, matrix_t* output);


/* Overwrites the given matrix so that it becomes an identity matrix.
   The input matrix must be of square dimensions. */
int matrix_set_identity(matrix_t mat);

/* Given a matrix named <mat>, return <mat ^ t>: meaning its transposed matrix.

   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
int matrix_transpose(matrix_t mat, matrix_t* output);

/* Given a matrix, return the sum of squared off-diagonals. */
double matrix_sum_squared_off(matrix_t mat);

/* Return the index of an off diagonal inside the one matrix's array that holds the largest absolute value.
 * Pre-condition: <sym_mat> must be a symmetric matrix!
 * Notice: due to the pre-condition, only the higher half of the matrix will be scanned for the largest off diagonal value.
 * Therefore, the index returned will and must be out of the higher half of the matrix. */  
matrix_ind matrix_ind_of_largest_offdiagonal(matrix_t sym_mat);


#endif /* MATRIX_H */
