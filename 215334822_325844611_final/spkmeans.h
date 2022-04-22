#ifndef SPKMEANS_H
#define SPKMEANS_H

#include "spkmeans_goals.h"
#include <string.h>

#define MAX_ITER 100
#define EPSILON 0.00001

typedef int make_iso_compilers_happy;

typedef struct {
    dpoint_t current_centroid;
    dpoint_t sum;
    int count;
} set_t;

/***************************** EXTERNAL VARIABLES *************************/
extern char *goal;
extern size_t K;

extern size_t dim;
extern size_t num_data;
extern dpoint_t *datapoints;
extern set_t *sets;
/*************************************************************************/

/**************************** MECHANISM'S INTERFACES
 * ************************************/
/* A function to pass data about the vectors into the main mechanism of this
 *file, produced by the spkmeansmodule.c interface. That includes:
 * 		1. the filename of the datapoints (C has a parsing function that handles
 *parsing of datapoints from a file)
 * 		2. output (a pointer to a variable of type matrix_t that the caller
 *wishes to store the T matrix at) a. in case output == NULL, print the matrix
 *outputted by the wanted goal b. in case output != NULL, store the matrix
 *outputted by the wanted goal into the output variable */
int spkmeans_pass_goal_info_and_run(char *infile, matrix_t *output);

/* A function that passes the initial centroids indices into the Kmeans
 * mechanism */
void spkmeans_pass_kmeans_info_and_run(size_t *initial_centroids_indices);
/*************************************************************************/

/**************************** AUXILIARY FUNCTIONS
 * ************************************/
/* A function used to assert a certain condition and print out "An error has
 * occurred" in case the condition isn't met */
void assert_other(bool condition);

/* A function used to initialize a datapoint (i.e., allocate an array to it) */
void init_datapoint(dpoint_t *dpoint);

/* A function used to free a datapoint */
void free_datapoint(dpoint_t);

/* A function used to free the program of the kmeans algorithm */
void free_program(void);
/*************************************************************************/

#endif
