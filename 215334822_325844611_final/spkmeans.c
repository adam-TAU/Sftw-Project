#include "spkmeans.h"

/******************************************************************************/

static int handle_goal(matrix_t *output);
static void kmeans(size_t *initial_centroids_indices);

static const char *get_filename_ext(const char *filename);
static void collect_data(const char *filename);
static void initialize_sets(size_t *initial_centroids_indices);
static void get_num_and_dim(FILE *file);
static void parse_datapoint(FILE *file, dpoint_t *dpoint);
static void assign_to_closest(dpoint_t *dpoint);
static double sqdist(dpoint_t p1, dpoint_t p2);
static void add_to_set(set_t *set, dpoint_t dpoint);
static int update_centroid(set_t *set);
static void parse_args(int argc, char **argv, char **infile);

/**************************** AUXILIARY FUNCTIONS
 * *********************************/
void assert_input(bool condition) {
    if(!condition) {
        printf("Invalid Input!");
        free_program();
        exit(1);
    }
}

void assert_other(bool condition) {
    if(!condition) {
        printf("An Error Has Occurred");
        free_program();
        exit(1);
    }
}
/*****************************************************************************/

/*************************** VARIABLES
 * *****************************************/
char *goal = "no goal yet";
size_t K = 0;

size_t dim = 0;
size_t num_data = 0;
dpoint_t *datapoints = NULL;

set_t *sets = NULL;
/*****************************************************************************/

/********************** USED BY THE CPython INTERFACE
 * ******************************/
int spkmeans_pass_goal_info_and_run(char *infile, matrix_t *output) {
    collect_data(infile);
    return handle_goal(output);
}

void spkmeans_pass_kmeans_info_and_run(size_t *initial_centroids_indices) {
    kmeans(initial_centroids_indices);

    if(initial_centroids_indices != NULL) {
        free(initial_centroids_indices);
    }
}
/*****************************************************************************/

/*************************** 2 SEPARATE MAIN MECHANISMS
 * ************************************/

static int handle_goal(matrix_t *output) {
    int signal;

    /* Build the output corresponding to the wanted goal */
    if(strcmp(goal, "wam") == 0) {
        if((signal = build_weighted_adjacency_matrix(output)))
            goto error;
    }
    if(strcmp(goal, "ddg") == 0) {
        if((signal = build_diagonal_degree_matrix(output)))
            goto error;
    }
    if(strcmp(goal, "lnorm") == 0) {
        if((signal = build_normalized_laplacian(output)))
            goto error;
    }
    if(strcmp(goal, "jacobi") == 0) {
        if((signal = build_jacobi_output(output)))
            goto error;
    }
    if(strcmp(goal, "spk") == 0)
    { /* available only for the CPython interface */
        if((signal = build_T_of_spectral_kmeans(K, output)))
            goto error;
    }

    return 0;

error:

    /* This mechanism is the last endpoint which uses the datapoints, hence we
     * can free the resources responsively */
    if(NULL != datapoints) {
        size_t i;
        for(i = 0; i < num_data; i++) {
            free_datapoint(datapoints[i]);
        }
        free(datapoints);
    }

    return signal;
}

static void kmeans(size_t *initial_centroids_indices) {
    size_t i, iter, updated_centroids;

    initialize_sets(initial_centroids_indices);

    for(iter = 0; iter < MAX_ITER; iter++) {
        for(i = 0; i < num_data; i++) {
            assign_to_closest(&datapoints[i]);
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
/*****************************************************************************/

/****************************** MAIN FUNCTION
 * ************************************/
int main(int argc, char **argv) {
    char *infile;
    int signal;
    matrix_t output;

    /* Parse args and collect data from file */
    parse_args(argc, argv, &infile);
    collect_data(infile);

    /* Power the wanted goal */
    if((signal = handle_goal(&output))) {
        assert_other(false);
    }

    /* Print and free */
    matrix_print_rows(output);
    matrix_free_safe(output);

    return 0;
}
/*****************************************************************************/

/***************************** KMEANS++ MECHANISM **************************/
/* Assigns the given datapoint to the closest set that it can find, using the
 * sqdist function. */
static void assign_to_closest(dpoint_t *dpoint) {
    size_t i, min_idx = 0;
    double min_dist = -1.0;

    for(i = 0; i < K; i++) {
        double dist = sqdist(sets[i].current_centroid, *dpoint);

        if((min_dist < 0.0) || (dist < min_dist)) {
            min_idx = i;
            min_dist = dist;
        }
    }

    add_to_set(&sets[min_idx], *dpoint);
    dpoint->current_set = min_idx;
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
static void initialize_sets(size_t *initial_centroids_indices) {
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
            sets[i].current_centroid.data[j] =
                datapoints[initial_centroids_indices[i]].data[j];
        }
    }
}

/* Given a file name, it returns the file extension of the file */
const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename)
        return "";
    return dot + 1;
}

/* Given an input filename, gathers all of the datapoints stored in that file,
 * while also figuring out what `dim` and `num_data` are supposed to be. */
static void collect_data(const char *filename) {
    FILE *input;
    size_t i;

    /* Asserting that the file extension is either .csv or .txt */
    const char *file_ext = get_filename_ext(filename);
    assert_input((strcmp(file_ext, "csv") == 0) ||
                 (strcmp(file_ext, "txt") == 0));

    /* Extracting the data from the input file */
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

    assert_input(argc == 3);

    goal = argv[1];
    assert_input(!(strcmp(goal, "wam") && strcmp(goal, "ddg") &&
                   strcmp(goal, "lnorm") && strcmp(goal, "jacobi")));

    *infile = argv[2];
}

/* Initializes a single datapoint - allocates enough space for it and sets all
 * the values to zero. */
void init_datapoint(dpoint_t *dpoint) {
    assert_other(dim > 0);

    dpoint->data = calloc(dim, sizeof(*dpoint->data));
    assert_other(NULL != dpoint->data);

    dpoint->current_set = (size_t)-1;
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
    size_t i = 0;

    if(NULL != sets) {
        for(i = 0; i < K; i++) {
            free_datapoint(sets[i].current_centroid);
            free_datapoint(sets[i].sum);
        }
        free(sets);
    }

    if(NULL != datapoints) {
        for(i = 0; i < num_data; i++) {
            free_datapoint(datapoints[i]);
        }
        free(datapoints);
    }
}
/*****************************************************************************/
