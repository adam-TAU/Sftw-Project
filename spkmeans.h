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
extern char* goal;
extern size_t K;

extern size_t dim;
extern size_t num_data;
extern dpoint_t *datapoints;
/*************************************************************************/







/*********************** EXTERNAL FUNCTIONS ******************************/
void init_datapoint(dpoint_t *dpoint);
void free_program(void);
/*************************************************************************/



/**************************** API ************************************/
/* A function to pass data about the vectors into the main mechanism of this file, produced by the spkmeansmodule.c interface.
 * That includes:
 * 		1. vectors (the vectors we want to process - corresponding to the wanted goal)
 *		2. vectors_dim (the dimension of each given vector)
 * 		3. vectors_amount (the amount of vectors given)
 *		4. goal (what we wish to do with the matrix formed from the given vectors)
 * 		5. (Optional) output (a pointer to a variable of type matrix_t that the caller wishes to store the T matrix at (relevant to goal == "spk" only) */
int spkmeans_pass_goal_info_and_run(char *infile, matrix_t *output);


/* A function that passes the initial centroids indices into the Kmeans mechanism */
void spkmeans_pass_kmeans_info_and_run(size_t *initial_centroids_indices);
/*************************************************************************/


#endif
