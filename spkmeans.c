#include "spkmeans.h"


static void assert_input(bool condition);
static void assert_other(bool condition);
static void collect_data(const char *filename);
static void initialize_sets(void);
static void get_num_and_dim(FILE *file);
static void parse_datapoint(FILE *file, dpoint_t *dpoint);
static void assign_to_closest(dpoint_t dpoint);
static double sqdist(dpoint_t p1, dpoint_t p2);
static void add_to_set(set_t *set, dpoint_t dpoint);
static int update_centroid(set_t *set);
static void parse_args(int argc, char **argv, char **infile);
static void init_datapoint(dpoint_t *dpoint);
static void free_datapoint(dpoint_t);
static void free_program(void);

/*****************************************************************************/

static void assert_input(bool condition) {
    if(!condition) {
        puts("Invalid Input!");
        free_program();
        exit(1);
    }
}

static void assert_other(bool condition) {
    if(!condition) {
        puts("An Error Has Occurred");
        free_program();
        exit(1);
    }
}



/*****************************************************************************/

size_t K = 0;
size_t dim = 0;
size_t num_data = 0;
dpoint_t *datapoints = NULL;

char* goal;
set_t *sets = NULL;


/************************* INTERFACE FOR GOALS *******************************/

void print_weighted_adjacency_matrix(dpoint_t vectors[]) {
	matrix_t output;
	output = graph_adjacent_matrix(vectors, dim);
	matrix_print_rows(output);
	matrix_free_safe(output);
}

void print_diagonal_degree_matrix(dpoint_t vectors[]) {
	matrix_t WAM, output;
	WAM = graph_adjacent_matrix(vectors, dim);
	output = graph_diagonal_degree_matrix(WAM);
	
	matrix_print_rows(output);
	matrix_free_safe(output);
	matrix_free_safe(WAM);
}

void print_normalized_laplacian(dpoint_t vectors[]) {
	matrix_t output;
	output = graph_normalized_laplacian(vectors, dim);
	matrix_print_rows(output);
	matrix_free_safe(output);
}

void print_jacobi_output(dpoint_t vectors[]) {
	size_t i, j;
	matrix_t jacobi_input;
	jacobi_output output;
	
	assert_other(num_data == dim);
	jacobi_input = matrix_new(num_data, dim);
	
	for (i = 0; i < num_data; i++) {
		for (j = 0; j < dim; j++) {
			matrix_set(jacobi_input, i, j, vectors[i].data[j]); 
		}
	}

	output = eigen_jacobi(jacobi_input);
	eigen_print_jacobi(output);
	
	matrix_free_safe(jacobi_input);
	eigen_free_jacobi_safe(output);
}

void print_spectral_kmeans(dpoint_t vectors[], size_t K) {
	printf("nothing %li %f\n", K, vectors[0].data[0]);
}


/*****************************************************************************/

int main(int argc, char **argv) {
    char *infile;
    size_t i, iter, updated_centroids;

    parse_args(argc, argv, &infile);

    collect_data(infile);
    printf("dim = %li, N = %li\n", dim, num_data);

    initialize_sets();

    for(iter = 0; iter < MAX_ITER; iter++) {
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

    free_program();
    return 0;
}

/* Assigns the given datapoint to the closest set that it can find, using the
 * sqdist function. */
static void assign_to_closest(dpoint_t dpoint) {
    size_t i, min_idx = 0;
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
static int update_centroid(set_t *set) {
    double dist;
    size_t i;

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
static double sqdist(dpoint_t p1, dpoint_t p2) {
    double dot = 0;
    size_t i;

    for(i = 0; i < dim; i++) {
        double temp = p1.data[i] - p2.data[i];
        dot += temp * temp;
    }

    return dot;
}

/* Adds the given datapoint to the provided set, taking into account both the
 * `sum` and `count` properties. */
static void add_to_set(set_t *set, dpoint_t dpoint) {
    size_t i;

    set->count += 1;
    for(i = 0; i < dim; i++) {
        set->sum.data[i] += dpoint.data[i];
    }
}

/* Initializes all of the sets, both allocating memory for the `sum` and
 * `current_centroid` properties and copying the data from the relevant
 * datapoint. */
static void initialize_sets() {
    size_t i, j;

    sets = calloc(K, sizeof(*sets));
    assert_other(NULL != sets);

    for(i = 0; i < K; i++) {
        /* count is already zero. We just need to allocate the centroid and sum
           datapoints. */
        init_datapoint(&sets[i].sum);
        init_datapoint(&sets[i].current_centroid);

        /* Copy initial current_centroid from i-th datapoint */
        for(j = 0; j < dim; j++) {
            sets[i].current_centroid.data[j] = datapoints[i].data[j];
        }
    }
}

/* Given an input filename, gathers all of the datapoints stored in that file,
 * while also figuring out what `dim` and `num_data` are supposed to be. */
static void collect_data(const char *filename) {
    FILE *input;
    size_t i;

    input = fopen(filename, "r");
    assert_input(NULL != input);
    get_num_and_dim(input);

    datapoints = calloc(num_data, sizeof(*datapoints));
    assert_other(NULL != datapoints);

    for(i = 0; i < num_data; i++) {
        parse_datapoint(input, &datapoints[i]);
    }

    fclose(input);
}

/* Parses a single datapoint from the given file, assuming that `dim` has
 * already been figured out. */
static void parse_datapoint(FILE *file, dpoint_t *dpoint) {
    size_t i;

    init_datapoint(dpoint);

    for(i = 0; i < dim; i++) {
        /* The following ',' is okay, because even if it isn't found parsing
           will be successful. */
        fscanf(file, "%lf,", &dpoint->data[i]);
    }

    /* Get rid of extra whitespace. */
    fscanf(file, "\n");
}

/* Determines `num_data` and `dim` from the current file by inspecting line
 * structure and amount. */
static void get_num_and_dim(FILE *file) {
    int c;

    dim = 1; /* Starting with 1 because the amount of numbers is always 1 more
                than the amount of commas. */
    num_data = 0;

    rewind(file);
    while(EOF != (c = fgetc(file))) {
        if(c == '\n') {
            num_data++;
        } else if(c == ',' && num_data == 0) {
            dim++;
        }
    }
    rewind(file);
}

/* Parses the arguments given to the program into K, MAX_ITER, input_file and
 * output_file. */
static void parse_args(int argc, char **argv, char **infile) {

    assert_input(argc == 2);

    goal = argv[1];
    assert_input(! ( strcmp(goal, "wam") && strcmp(goal, "ddg") && strcmp(goal, "lnorm") && strcmp(goal, "jacobi") ) );

    *infile = argv[2];
}

/* Initializes a single datapoint - allocates enough space for it and sets all
 * the values to zero. */
static void init_datapoint(dpoint_t *dpoint) {
    assert_other(dim > 0);

    dpoint->data = calloc(dim, sizeof(*dpoint->data));
    assert_other(NULL != dpoint->data);
}

/* Frees the given datapoint. If it's already been freed or not yet allocated,
 * this function safely does nothing. */
static void free_datapoint(dpoint_t dpoint) {
    if(NULL != dpoint.data) {
        free(dpoint.data);
    }
}

/* Frees all of the memory allocated by the program. If a certain variable
 * hasn't been allocated yet, this function does not attempt to free it. */
static void free_program() {
    size_t i;

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
}

