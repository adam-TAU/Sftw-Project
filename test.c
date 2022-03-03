#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define EPSILON 0.000001

int eq_loose(const double *data1, const double *data2, size_t count) {
    size_t i;
    for(i = 0; i < count; i++) {
        if(fabs(data1[i] - data2[i]) >= EPSILON) {
            return 0;
        }
    }

    return 1;
}

int simple_add() {
    double data1[6] = {1, 2, 3, 4, 5, 6};
    double data2[6] = {2, 3, 4, 5, 6, 7};
    double data3[6] = {3, 5, 7, 9, 11, 13};
    int err;
    matrix_t m1, m2, m3;
    m1 = matrix_new(2, 3);
    m2 = matrix_new(2, 3);

    memcpy(m1.data, data1, 6 * sizeof(double));
    memcpy(m2.data, data2, 6 * sizeof(double));

    if(0 != (err = matrix_add(m1, m2, &m3))) {
        puts("simple_add failed! error while adding");
        return 0;
    }

    if(!eq_loose(data3, m3.data, m3.len)) {
        puts("simple_add failed! unexpected result:");
        matrix_print(m3);
        return 0;
    }

    return 1;
}

int main() {
    size_t i, test_len = 1;
    int (*tests[])(void) = {simple_add};

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
