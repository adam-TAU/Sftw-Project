#ifndef SPKMEANS_H
#define SPKMEANS_H


#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define EPSILON 0.001
#define bool int
#define true 1
#define false 0

typedef struct {
    double *data;
} dpoint_t;

typedef struct {
    dpoint_t current_centroid;
    dpoint_t sum;
    int count;
} set_t;




void assert_other(bool condition);
void initialize_sets(int*);
void assign_to_closest(dpoint_t dpoint);
double sqdist(dpoint_t p1, dpoint_t p2);
void add_to_set(set_t *set, dpoint_t dpoint);
int update_centroid(set_t *set);
void init_datapoint(dpoint_t *dpoint);
void free_datapoint(dpoint_t);
void free_program(void);
void converge(int max_iter);
void parse_args(void);








#endif
