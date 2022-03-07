#include "graph.h"
#include "eigen.h"
#include <stdio.h>
#include <stdlib.h>

typedef int make_iso_compilers_happy;


int main(int argv, char* args[]) {
	if (argv > 0) printf("%s\n", args[0]);
	printf("success: %d\n", argv);
	return 0;
}




/*****************************************************************************/


void assert_other(bool condition) {
    if(!condition) {
        puts("An Error Has Occurred");
        free_program();
        exit(1);
    }
}


/*****************************************************************************/

double** datapoints_arg = NULL;
int* initial_centroids_indices = NULL;

int dim;
int num_data;
dpoint_t *datapoints = NULL;

int K;
set_t *sets = NULL;

double epsilon;

double** centroids_c = NULL;

/*****************************************************************************/



/* This function parses the arguments retrieved from Python, into our
global arguments that serve the k-means algorithm) */
void parse_args(void) {
        /* actual parsing */
        #undef EPSILON
        #define EPSILON epsilon

        int i;
        datapoints = calloc(num_data, sizeof(*datapoints));
        assert_other(NULL != datapoints);
        for(i = 0; i < num_data; i++) {
                datapoints[i].data = datapoints_arg[i];
        }
        initialize_sets(initial_centroids_indices);
}





/* This function was created to reduce code duplication. It is the engine of the algorithm,
 * providing the structurization of the clusters until convergence
 */
void converge(int max_iter) {
        int iter, i, updated_centroids;
    for(iter = 0; iter < max_iter; iter++) {
        for(i = 0; i < num_data; i++) {
            assign_to_closest(datapoints[i]);
        }

        updated_centroids = 0;
        for(i = 0; i < K; i++) {
            updated_centroids += update_centroid(&sets[i]);
        }

        if(updated_centroids == 0) { /* Convergence */
            break;
        }
    }
}




/* Assigns the given datapoint to the closest set that it can find, using the
 * sqdist function. */
void assign_to_closest(dpoint_t dpoint) {
    int i, min_idx = 0;
    double min_dist = -1.0;

    for(i = 0; i < K; i++) {
        double dist = sqdist(sets[i].current_centroid, dpoint);

        if((min_dist < 0.0) || (dist < min_dist)) {
            min_idx = i;
            min_dist = dist;
        }
    }

    add_to_set(&sets[min_idx], dpoint);
}

/* Updates the centroid of the given set using its stored `sum` and `count`
 * properties, while also resetting them to 0 for the next iteration. */
int update_centroid(set_t *set) {
    double dist;
    int i;

    for(i = 0; i < dim; i++) {
        set->sum.data[i] /= (double)set->count;
    }

    dist = sqrt(sqdist(set->sum, set->current_centroid));

    for(i = 0; i < dim; i++) {
        set->current_centroid.data[i] = set->sum.data[i];
        set->sum.data[i] = 0.0;
    }

    set->count = 0;

    return (dist >= EPSILON) ? 1 /* If this set's centroid changed, return 1 */
                             : 0;
}

/* Calculates the squared distance between two given datapoints. */
double sqdist(dpoint_t p1, dpoint_t p2) {
    double dot = 0;
    int i;

    for(i = 0; i < dim; i++) {
        double temp = p1.data[i] - p2.data[i];
        dot += temp * temp;
    }

    return dot;
}

/* Adds the given datapoint to the provided set, taking into account both the
 * `sum` and `count` properties. */
void add_to_set(set_t *set, dpoint_t dpoint) {
    int i;

    set->count += 1;
    for(i = 0; i < dim; i++) {
        set->sum.data[i] += dpoint.data[i];
    }
}



/* Initializes all of the sets, both allocating memory for the `sum` and
 * `current_centroid` properties and copying the data from a given index list
 *  that indicates what datapoint we should initialize the current centroid with.
 *  In the case of a regular kmeans.c run in HW1, we would power this with the plain [0,1,...,K-1] index list.
 *  Though, with the kmeans++, we will power it with the observation indices list */
void initialize_sets(int* indices) {
    int i, j;

    sets = calloc(K, sizeof(*sets));
    assert_other(NULL != sets);

    for(i = 0; i < K; i++) {
        /* count is already zero. We just need to allocate the centroid and sum
           datapoints. */
        init_datapoint(&sets[i].sum);
        init_datapoint(&sets[i].current_centroid);

        /* Copy initial current_centroid from i-th datapoint */
        for(j = 0; j < dim; j++) {
            sets[i].current_centroid.data[j] = datapoints[indices[i]].data[j];
        }
    }
}




/* Initializes a single datapoint - allocates enough space for it and sets all
 * the values to zero. */
void init_datapoint(dpoint_t *dpoint) {
    assert_other(dim > 0);

    dpoint->data = calloc(dim, sizeof(*dpoint->data));
    assert_other(NULL != dpoint->data);
}

/* Frees the given datapoint. If it's already been freed or not yet allocated,
 * this function safely does nothing. */
void free_datapoint(dpoint_t dpoint) {
    if(NULL != dpoint.data) {
        free(dpoint.data);
    }
}

/* Frees all of the memory allocated by the program. If a certain variable
 * hasn't been allocated yet, this function does not attempt to free it. */
void free_program() {
    int i;

    if(NULL != datapoints) {
        for(i = 0; i < num_data; i++) {
            free_datapoint(datapoints[i]);
        }
        free(datapoints);
    }

    if(NULL != sets) {
        for(i = 0; i < K; i++) {
            free_datapoint(sets[i].current_centroid);
            free_datapoint(sets[i].sum);
        }
        free(sets);
    }

    /* free-ing the Python-given datapoints */
    if(NULL != datapoints_arg) {
        /* the actual vectors of the datapoints have been linked
         * directly to datapoints[i].data (line 107, function parse_args()) -> hence why we don't free
         * each vector of datapoints_arg individually */
        free(datapoints_arg);
    }

    /* free-ing the Python-given centroids-indices */
    if(NULL != initial_centroids_indices) free(initial_centroids_indices);

    /* free-ing the centroids which were used to build values to return to Python
     * We've already free-d each particular centroid with in the free_datapoint(sets[i].current_centroid) */
    if(NULL != centroids_c) { 
        free(centroids_c);
    }
}
