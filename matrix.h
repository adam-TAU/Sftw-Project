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


/* Define a structure that will hold a vector's coordinates. This will also be used in the spkmeans.c mechanism */
typedef struct {
    double *data;
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
matrix_t matrix_new(size_t rows, size_t cols);

/* Clones the given matrix into a newly allocated matrix.

   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
matrix_t matrix_clone(matrix_t mat);

/* Swaps the data contained in both matrix instances.
   This doesn't allocate, and is useful for loops. */
void matrix_swap(matrix_t *mat1, matrix_t *mat2);

/* Builds a matrix out of an already existing dataset of <dpoint_t>s.
   
   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
matrix_t matrix_build(dpoint_t* vectors, size_t num_vectors, size_t dim);

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
matrix_t matrix_mul_scalar(matrix_t mat, double scalar);

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


/* Multiplies the given matrices to produce a new matrix, and stores it into <mat1>.
   Precondition: `mat1.cols == mat1.rows == mat2.rows == mat2.cols`.

   Errors:

   In case of dimension mismatch, the return value is `DIM_MISMATCH`.
   In case of allocation failure, the return value is `BAD_ALLOC`. */
int matrix_mul_assign_to_first(matrix_t *mat1, matrix_t mat2);


/* Multiplies the given matrices to produce a new matrix, and stores it into <mat2>.
   Precondition: `mat1.cols == mat1.rows == mat2.rows == mat2.cols`.

   Errors:

   In case of dimension mismatch, the return value is `DIM_MISMATCH`.
   In case of allocation failure, the return value is `BAD_ALLOC`. */
int matrix_mul_assign_to_second(matrix_t mat1, matrix_t *mat2);



/* Create an Identity Matrix, with the dimensions of dim x dim.

   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
matrix_t matrix_identity(size_t dim);

/* Given a matrix named <mat>, return <mat ^ t>: meaning its transposed matrix.

   In case of allocation failure, the output matrix has a `data` field of
   `NULL`. */
matrix_t matrix_transpose(matrix_t mat);

#endif /* MATRIX_H */
