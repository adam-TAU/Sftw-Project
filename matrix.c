#include "matrix.h"
#include <stdio.h>

matrix_t matrix_new(size_t rows, size_t cols) {
    size_t len = rows * cols;
    double *data = calloc(len, sizeof(double));

    matrix_t output = {.data = data, .rows = rows, .cols = cols, .len = len};
    return output;
}

matrix_t matrix_clone(matrix_t mat) {
    matrix_t cloned = matrix_new(mat.rows, mat.cols);
    if(NULL == cloned.data) {
        return cloned; // return early to avoid accessing `data`
    }

    for(size_t idx = 0; idx < mat.len; idx++) {
        cloned.data[idx] = mat.data[idx];
    }

    return cloned;
}

size_t matrix_calc_index(matrix_t mat, size_t i, size_t j) {
    if(i >= mat.rows) {
        printf("invalid index for matrix: the number of rows is %zu but the "
               "index is %zu\n",
               mat.rows, i);
        return (size_t)-1;
    }

    if(j >= mat.cols) {
        printf("invalid index for matrix: the number of cols is %zu but the "
               "index is %zu\n",
               mat.cols, j);
        return (size_t)-1;
    }

    return i * mat.cols + j;
}

double matrix_get(matrix_t mat, size_t i, size_t j) {
    size_t index = matrix_calc_index(mat, i, j);
    return mat.data[index];
}

void matrix_free(matrix_t mat) {
    free(mat.data);
}

int matrix_add_assign(matrix_t self, matrix_t other) {
    if(!(self.rows == other.rows && self.cols == other.cols)) {
        return DIM_MISMATCH;
    }

    for(size_t idx = 0; idx < self.len; idx++) {
        self.data[idx] += other.data[idx];
    }

    return 0;
}

int matrix_add(matrix_t mat1, matrix_t mat2, matrix_t *output) {
    if(!(mat1.rows == mat2.rows && mat1.cols == mat2.cols)) {
        return DIM_MISMATCH;
    }

    *output = matrix_new(mat1.rows, mat1.cols);
    if(NULL == output->data) {
        return BAD_ALLOC;
    }

    for(size_t idx = 0; idx < mat1.len; idx++) {
        output->data[idx] = mat1.data[idx] + mat2.data[idx];
    }

    return 0;
}

void matrix_mul_scalar_assign(matrix_t mat, double scalar) {
    for(size_t idx = 0; idx < mat.len; idx++) {
        mat.data[idx] *= scalar;
    }
}

matrix_t matrix_mul_scalar(matrix_t mat, double scalar) {
    matrix_t output = matrix_new(mat.rows, mat.cols);
    if(NULL == output.data) {
        return output;
    }

    for(size_t idx = 0; idx < mat.len; idx++) {
        output.data[idx] = mat.data[idx] * scalar;
    }

    return output;
}

int matrix_mul(matrix_t mat1, matrix_t mat2, matrix_t *output) {
    if(mat1.cols != mat2.rows) {
        return DIM_MISMATCH;
    }

    *output = matrix_new(mat1.rows, mat2.cols);
    if(NULL == output->data) {
        return BAD_ALLOC;
    }

    size_t idx = 0;
    for(size_t i = 0; i < output->rows; i++) {
        for(size_t j = 0; j < output->cols; j++) {
            double sum = 0;
            for(size_t k = 0; k < mat1.cols; k++) {
                sum += matrix_get(mat1, i, k) * matrix_get(mat2, k, j);
            }

            output->data[idx] = sum;

            idx++;
        }
    }

    return 0;
}
