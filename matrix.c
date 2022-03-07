#include "matrix.h"
#include <stdio.h>

matrix_t matrix_new(size_t rows, size_t cols) {
    matrix_t output;
    size_t len = rows * cols;
    double *data = calloc(len, sizeof(double));

    output.data = data;
    output.rows = rows;
    output.cols = cols;
    output.len = len;

    return output;
}

matrix_t matrix_clone(matrix_t mat) {
    matrix_t cloned = matrix_new(mat.rows, mat.cols);
    size_t idx;
    if(NULL == cloned.data) {
        return cloned; /* return early to avoid accessing `data` */
    }

    for(idx = 0; idx < mat.len; idx++) {
        cloned.data[idx] = mat.data[idx];
    }

    return cloned;
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


void matrix_print(matrix_t mat) {
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

    *output = matrix_new(mat1.rows, mat1.cols);
    if(NULL == output->data) {
        return BAD_ALLOC;
    }

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

matrix_t matrix_mul_scalar(matrix_t mat, double scalar) {
    size_t idx;
    matrix_t output = matrix_new(mat.rows, mat.cols);
    if(NULL == output.data) {
        return output;
    }

    for(idx = 0; idx < mat.len; idx++) {
        output.data[idx] = mat.data[idx] * scalar;
    }

    return output;
}

int matrix_mul(matrix_t mat1, matrix_t mat2, matrix_t *output) {
    size_t idx, i;
    if(mat1.cols != mat2.rows) {
        return DIM_MISMATCH;
    }

    *output = matrix_new(mat1.rows, mat2.cols);
    if(NULL == output->data) {
        return BAD_ALLOC;
    }

    idx = 0;
    for(i = 0; i < output->rows; i++) {
        size_t j;
        for(j = 0; j < output->cols; j++) {
            size_t k;
            double sum = 0;
            for(k = 0; k < mat1.cols; k++) {
                sum += matrix_get(mat1, i, k) * matrix_get(mat2, k, j);
            }

            output->data[idx] = sum;

            idx++;
        }
    }

    return 0;
}



int matrix_mul_assign(matrix_t mat1, matrix_t mat2) {
	int signal; 
    matrix_t* output;
    output = malloc(sizeof(matrix_t));
	signal = matrix_mul(mat1, mat2, output);
	
	if (output->len == mat1.len) {
		matrix_free(mat1);
		mat1.data = output->data;
	} else {
		signal = 1;
	}
	
	matrix_free_safe(*output);
	free(output);
	return signal;
}


matrix_t matrix_identity_matrix(size_t dim) {
	size_t i;
	matrix_t output;
	output = matrix_new(dim, dim);
 
 	for (i = 0; i < dim; i++) {
		matrix_set(output, i, i, 1);
 	}
 	
 	return output;
}



matrix_t matrix_transpose(matrix_t mat) {
	size_t i, j;
	matrix_t output;
	
	output = matrix_new(mat.cols, mat.rows);
	for (i = 0; i < output.rows; i++) {
		for (j = 0; j < output.cols; j++) {
			matrix_set(output, i, j, matrix_get(mat, j, i));
		}
	}
	
	return output;
}










