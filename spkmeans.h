#ifndef SPKMEANS_H
#define SPKMEANS_H



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"
#include "graph.h"
#include "eigen.h"



#define MAX_ITER 300
#define EPSILON 0.001
#define bool int
#define true 1
#define false 0


typedef int make_iso_compilers_happy;

typedef struct {
    dpoint_t current_centroid;
    dpoint_t sum;
    int count;
} set_t;




void print_weighted_adjacency_matrix(dpoint_t vectors[]);
void print_diagonal_degree_matrix(dpoint_t vectors[]);
void print_normalized_laplacian(dpoint_t vectors[]);


/* A function to print the eigen values and eigen vectors of the given input. 
 * Pre-Conditions:
 * 		<mat> must be a "real" matrix
 * 		<mat> must be a "symmetric" matrix */
void print_jacobi_output(dpoint_t vectors[]);


void print_spectral_kmeans(dpoint_t vectors[], size_t K);





#endif
