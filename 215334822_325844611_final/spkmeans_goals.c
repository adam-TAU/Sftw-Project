#include "spkmeans_goals.h"

/************************** ADD ERROR HANDLING *******************************/

/* Global variables crucial for preforming the goals. They're defined by the
 * spkmeans.c file, no matter the situation, hence the 'extern' type */
extern size_t dim;
extern size_t num_data;
extern dpoint_t *datapoints;

extern void assert_other(int condition);
extern void assert_input(int condition);

/************************* INTERFACE FOR GOALS *******************************/

int build_weighted_adjacency_matrix(matrix_t *output) {

    /* Find the WAM matrix */
    if(graph_adjacent_matrix(datapoints, num_data, dim, output)) {
        return BAD_ALLOC;
    }

    return 0;
}

int build_diagonal_degree_matrix(matrix_t *output) {

    /* Find the DDG matrix */
    if(graph_diagonal_degree_matrix(datapoints, num_data, dim, false, output)) {
        return BAD_ALLOC;
    }

    return 0;
}

int build_normalized_laplacian(matrix_t *output) {

    /* Find the LNORM matrix */
    if(graph_normalized_laplacian(datapoints, num_data, dim, output)) {
        return BAD_ALLOC;
    }

    return 0;
}

int build_jacobi_output(matrix_t *output) {
    matrix_t jacobi_input;
    jacobi_t jacobi_res;

    /* making sure that the given vectors' dataset represents a symmetric matrix
     * (else jacobi isn't feasible) */
    assert_input(num_data == dim);

    /* Converting the input into a matrix and sending it into the jacobi
     * algorithm */
    if(matrix_build_from_dpoints(datapoints, num_data, dim, &jacobi_input))
        goto error;

    /* Extracting all of the eigen values (num_data eigen values) */
    if(eigen_jacobi(jacobi_input, num_data, &jacobi_res))
        goto error;

    /* Converting the output format from a <jacobi_res> into a <matrix_t> */
    if(eigen_jacobi_to_mat(jacobi_res, output))
        goto error;

    /* Free-ing the matrix that was created as the jacobi's algorithm's input */
    matrix_free(jacobi_input);

    return 0;

error:
    matrix_free_safe(jacobi_input);
    if(NULL != jacobi_res.eigen_values)
        free(jacobi_res.eigen_values);
    matrix_free_safe(jacobi_res.eigen_vectors);
    return BAD_ALLOC;
}

int build_T_of_spectral_kmeans(size_t K, matrix_t *output) {
    matrix_t L_norm;
    jacobi_t jacobi_res;
    size_t i, j;

    /* Finding the graph normalized laplacian matrix */
    if(graph_normalized_laplacian(datapoints, num_data, dim, &L_norm))
        goto error;

    /* Applying the jacbobi algorithm upon the graph normalized laplacian
     * matrix. This extracts the first k eigen values and their corresponding
     * eigen vectors, sortedly */
    if(eigen_jacobi(L_norm, K, &jacobi_res))
        goto error;

    /* Creating the T matrix */
    if(matrix_new(jacobi_res.eigen_vectors.rows, jacobi_res.eigen_vectors.cols,
                  output))
        goto error;

    /* Building the T matrix */
    for(i = 0; i < output->rows; i++) {
        double sum_squared_of_rows = 0;
        double norm_of_row;

        /* Calculating the sum of squared of the row */
        for(j = 0; j < output->cols; j++) {
            sum_squared_of_rows +=
                pow(matrix_get(jacobi_res.eigen_vectors, i, j), 2);
        }

        /* Calculating the norm of the row with the sum of squared of the row */
        if(0 == (norm_of_row = pow(sum_squared_of_rows, 0.5)))
        { /* if the sum of the row equals 0, keep the row as is */
            norm_of_row = 1;
        }

        /* Normalize the jacobi output into the new matrix: T */
        for(j = 0; j < output->cols; j++) {
            matrix_set(*output, i, j,
                       matrix_get(jacobi_res.eigen_vectors, i, j) /
                           norm_of_row);
        }
    }

    /* Free-ing and Returning */
    matrix_free(L_norm);
    free(jacobi_res.eigen_values);
    matrix_free(jacobi_res.eigen_vectors);

    return 0;

error:
    matrix_free_safe(L_norm);
    if(NULL != jacobi_res.eigen_values)
        free(jacobi_res.eigen_values);
    matrix_free_safe(jacobi_res.eigen_vectors);
    matrix_free_safe(*output);
    return BAD_ALLOC;
}

/*****************************************************************************/
