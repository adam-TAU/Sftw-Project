#include "matrix.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

int matrix_new(size_t rows, size_t cols, matrix_t *output) {
	size_t len = rows * cols;
	double *data = calloc(len, sizeof(double));
	if (NULL == data) return BAD_ALLOC;

	output->data = data;
	output->rows = rows;
	output->cols = cols;
	output->len = len;

	return 0;
}

int matrix_clone(matrix_t mat, matrix_t *output) {

	if (matrix_new(mat.rows, mat.cols, output)) return BAD_ALLOC;

	matrix_copy(*output, mat); /* no need to check error code - dims are equal */

	return 0;
}

void matrix_swap(matrix_t *mat1, matrix_t *mat2) {
    matrix_t temp = *mat1;
    *mat1 = *mat2;
    *mat2 = temp;
}

int matrix_copy(matrix_t dest, matrix_t src) {
    if(!(dest.rows == src.rows && dest.cols == src.cols)) return DIM_MISMATCH;
    memcpy(dest.data, src.data, sizeof(double) * src.len);
    return 0;
}


int matrix_build_from_dpoints(dpoint_t* vectors, size_t num_vectors, size_t dim, matrix_t *output) {
	size_t i, j;
	if (matrix_new(num_vectors, dim, output)) return BAD_ALLOC;

	for (i = 0; i < output->rows; i++) {
		for (j = 0; j < output->cols; j++) {
			matrix_set(*output, i, j, vectors[i].data[j]); 
		}
	}

	return 0;
}


size_t matrix_calc_index(matrix_t mat, size_t i, size_t j) {
	if(i >= mat.rows) {
		printf("invalid index for matrix: the number of rows is %lu but the "
				"index is %lu\n",
				mat.rows, i);
		return (size_t)-1;
	}

	if(j >= mat.cols) {
		printf("invalid index for matrix: the number of cols is %lu but the "
				"index is %lu\n",
				mat.cols, j);
		return (size_t)-1;
	}

	return i * mat.cols + j;
}

double matrix_get(matrix_t mat, size_t i, size_t j) {
	size_t index = matrix_calc_index(mat, i, j);
	return mat.data[index];
}


void matrix_set(matrix_t mat, size_t i, size_t j, double val) {
	size_t index = matrix_calc_index(mat, i, j);
	mat.data[index] = val;
}


void matrix_petty_print(matrix_t mat) {
	size_t idx, row, col;
	printf("matrix(%lu x %lu):\n", mat.rows, mat.cols);

	idx = 0;
	for(row = 0; row < mat.rows; row++) {
		for(col = 0; col < mat.cols; col++) {
			printf("%.4f ", mat.data[idx]);
			idx++;
		}
		puts("");
	}
}


void matrix_print_rows(matrix_t mat) {
	size_t idx, row, col;

	idx = 0;
	for(row = 0; row < mat.rows; row++) {
		for(col = 0; col < mat.cols; col++) {
			printf("%.4f", mat.data[idx]);
			if (col < mat.cols - 1) printf(",");
			idx++;
		}
		puts("");
	}
}


void matrix_print_cols(matrix_t mat) {
	size_t i, j;

	for(j = 0; j < mat.cols; j++) {
		for(i = 0; i < mat.rows; i++) {
			printf("%.4f", matrix_get(mat, i, j));
			if (i < mat.rows - 1) printf(",");
		}
		puts("");
	}
}


void matrix_free(matrix_t mat) {
	free(mat.data);
}

void matrix_free_safe(matrix_t mat) {
	if(NULL != mat.data) {
		matrix_free(mat);
	}
}

int matrix_add_assign(matrix_t self, matrix_t other) {
	size_t idx;
	if(!(self.rows == other.rows && self.cols == other.cols)) {
		return DIM_MISMATCH;
	}

	for(idx = 0; idx < self.len; idx++) {
		self.data[idx] += other.data[idx];
	}

	return 0;
}

int matrix_add(matrix_t mat1, matrix_t mat2, matrix_t *output) {
	size_t idx;
	if(!(mat1.rows == mat2.rows && mat1.cols == mat2.cols)) {
		return DIM_MISMATCH;
	}

	if (matrix_new(mat1.rows, mat1.cols, output)) return BAD_ALLOC;

	for(idx = 0; idx < mat1.len; idx++) {
		output->data[idx] = mat1.data[idx] + mat2.data[idx];
	}

	return 0;
}

void matrix_mul_scalar_assign(matrix_t mat, double scalar) {
	size_t idx;
	for(idx = 0; idx < mat.len; idx++) {
		mat.data[idx] *= scalar;
	}
}

int matrix_mul_scalar(matrix_t mat, double scalar, matrix_t* output) {
	size_t idx;
	if (matrix_new(mat.rows, mat.cols, output)) return BAD_ALLOC;

	for(idx = 0; idx < mat.len; idx++) {
		output->data[idx] = mat.data[idx] * scalar;
	}

	return 0;
}

int matrix_mul(matrix_t mat1, matrix_t mat2, matrix_t *output) {
	if (matrix_new(mat1.rows, mat2.cols, output)) return BAD_ALLOC;
	return matrix_mul_buffer(mat1, mat2, *output);
}


int matrix_mul_buffer(matrix_t mat1, matrix_t mat2, matrix_t output) {
    size_t idx, i, j, k;

    if(!(mat1.cols == mat2.rows && mat1.rows == output.rows && mat2.cols == output.rows)) return DIM_MISMATCH;
    
    idx = 0;
	for(i = 0; i < output.rows; i++) {
		for(j = 0; j < output.cols; j++) {
			double sum = 0;
			for(k = 0; k < mat1.cols; k++) {
				sum += matrix_get(mat1, i, k) * matrix_get(mat2, k, j);
			}

			output.data[idx] = sum;

			idx++;
		}
	}

    return 0;
}



int matrix_identity(size_t dim, matrix_t* output) {
	size_t i;
	if (matrix_new(dim, dim, output)) return BAD_ALLOC;

	for (i = 0; i < dim; i++) {
		matrix_set(*output, i, i, 1);
	}

	return 0;
}

int matrix_set_identity(matrix_t mat) {
    size_t i, j, idx;
    if(mat.rows != mat.cols) return DIM_MISMATCH;

    idx = 0;
    for(i = 0; i < mat.rows; i++) {
        for(j = 0; j < mat.cols; j++) {
            mat.data[idx] = (i == j) ? 1 : 0;
            idx++;
        }
    }

    return 0;
}



int matrix_transpose(matrix_t mat, matrix_t* output) {
	size_t i, j;
	if (matrix_new(mat.cols, mat.rows, output)) return BAD_ALLOC;

	for (i = 0; i < output->rows; i++) {
		for (j = 0; j < output->cols; j++) {
			matrix_set(*output, i, j, matrix_get(mat, j, i));
		}
	}

	return 0;
}




double matrix_sum_squared_off(matrix_t mat) {
	double sum = 0.0;
	size_t i, j;

	for (i = 0; i < mat.rows; i++) {
		for (j = 0; j < mat.cols; j++) {
			sum += (i != j) ? pow( matrix_get(mat, i, j), 2 ) : 0;
		}
	}

	return sum;
}



matrix_ind matrix_ind_of_largest_offdiagonal(matrix_t sym_mat) {
	matrix_ind output;
	size_t i, j;
	double current_max = -1;

	for (i = 0; i < sym_mat.rows; i++) {
		for (j = i + 1; j < sym_mat.cols; j++) {
			double tmp = fabs(matrix_get(sym_mat, i, j));
			if (tmp > current_max) {
				output.i = i;
				output.j = j;
				current_max = tmp;
			}
		}
	}
	return output;
}


bool matrix_is_diagonal(matrix_t mat) {
	size_t i, j;
	
	for (i = 0; i < mat.rows; i++) {
		for (j = 0; j < mat.cols; j++) {
			if (i != j) {
				if (0 != matrix_get(mat, i, j)) return false;
			}
		}
	}
	
	return true;
}

