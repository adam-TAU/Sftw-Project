#ifndef SPKMEANS_H
#define SPKMEANS_H


#include "spkmeans_goals.h"
#include <string.h>



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


/* A function to pass data about the vectors into the main mechanism of this file, produced by the spkmeansmodule.c interface.
 * That includes:
 * 		1. vectors (the vectors we want to process - corresponding to the wanted goal)
 *		2. vectors_dim (the dimension of each given vector)
 * 		3. vectors_amount (the amount of vectors given)
 *		4. goal (what we wish to do with the matrix formed from the given vectors)
 * 		5. (Optional) output (a pointer to a variable of type matrix_t that the caller wishes to store the T matrix at (relevant to goal == "spk" only) */
void spkmeans_pass_vectors_info(dpoint_t *vectors_from_py, size_t vectors_dim_from_py, size_t vectors_amount_from_py, size_t K_from_py, char* goal_from_py, matrix_t* output);


/* A function to pass data about the datapoints and the crucial fectors when clustering datasets using the K-means algorithm.
 * That includes:
 *		1. datapoints
 * 		2. K (amount of clusters)
 *		3. dim (dimension of each datapoint)
 * 		4. num_data (the amount of datapoints we have to cluster) */
void spkmeans_pass_kmeans_info(dpoint_t *datapoints_from_py, size_t K_from_py, size_t dim_from_py, size_t num_data_from_py);

#endif
