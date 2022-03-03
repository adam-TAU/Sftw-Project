#include "bool.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define EPSILON 0.000001

static bool eq_loose(const double *data1, const double *data2, size_t count) {
    size_t i;
    for(i = 0; i < count; i++) {
        if(fabs(data1[i] - data2[i]) >= EPSILON) {
            return false;
        }
    }

    return true;
}

static bool simple_add(void) {
    double data1[6] = {1, 2, 3, 4, 5, 6};
    double data2[6] = {2, 3, 4, 5, 6, 7};
    double data3[6] = {3, 5, 7, 9, 11, 13};
    int err;
    matrix_t m1, m2, m3;
    m1 = matrix_new(2, 3);
    m2 = matrix_new(2, 3);

    if(NULL == m1.data || NULL == m2.data) {
        memset(&m3, 0, sizeof(matrix_t));
        goto error;
    }

    memcpy(m1.data, data1, m1.len * sizeof(double));
    memcpy(m2.data, data2, m2.len * sizeof(double));

    if(0 != (err = matrix_add(m1, m2, &m3))) {
        puts("simple_add failed! error while adding");
        goto error;
    }

    if(!eq_loose(data3, m3.data, m3.len)) {
        puts("simple_add failed! unexpected result:");
        matrix_print(m3);
        goto error;
    }

    return true;

error:
    matrix_free_safe(m1);
    matrix_free_safe(m2);
    matrix_free_safe(m3);
    return false;
}

static bool simple_mul(void) {
    double data1[6] = {1, 2, 3, 4, 5, 6};
    double data2[6] = {2, 3, 4, 5, 6, 7};
    double data3[4] = {28, 34, 64, 79};
    matrix_t m1, m2, m3;
    int err;

    m1 = matrix_new(2, 3);
    m2 = matrix_new(3, 2);

    if(NULL == m1.data || NULL == m2.data) {
        memset(&m3, 0, sizeof(matrix_t));
        goto error;
    }

    memcpy(m1.data, data1, m1.len * sizeof(double));
    memcpy(m2.data, data2, m2.len * sizeof(double));

    if(0 != (err = matrix_mul(m1, m2, &m3))) {
        puts("simple_mul failed! error while multiplying");
        goto error;
    }

    if(!eq_loose(data3, m3.data, m3.len)) {
        puts("simple_mul failed! unexpected result:");
        matrix_print(m3);
        goto error;
    }

    return true;

error:
    matrix_free_safe(m1);
    matrix_free_safe(m2);
    matrix_free_safe(m3);
    return false;
}

int main(void) {
    size_t i, test_len = 2;
    int (*tests[])(void) = {simple_add, simple_mul};

    int passed_tests = 0;
    for(i = 0; i < test_len; i++) {
        if(tests[i]()) {
            passed_tests++;
        }
        puts("------");
    }

    printf("\nPassed tests: %d/%lu\n", passed_tests, test_len);

    return 0;
}
