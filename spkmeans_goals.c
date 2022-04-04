#include "spkmeans_goals.h"


/************************** ADD ERROR HANDLING *******************************/

/* Global variables crucial for preforming the goals. They're defined by the spkmeans.c file, no matter the situation, hence the 'extern' type */
extern size_t dim;
extern size_t num_data;
extern dpoint_t *datapoints;

extern void assert_other(int condition);

/************************* INTERFACE FOR GOALS *******************************/

int print_weighted_adjacency_matrix() {
	matrix_t output;
	if (0 != graph_adjacent_matrix(datapoints, num_data, dim, &output)) {
		return BAD_ALLOC;
	}

	/* Printing and free-ing */
	matrix_print_rows(output);
	matrix_free(output);
	return 0;
}

int print_diagonal_degree_matrix() {
	matrix_t WAM, output;
	
	if (0 != graph_adjacent_matrix(datapoints, num_data, dim, &WAM)) goto error;
	
	if (0 != graph_diagonal_degree_matrix(WAM, true, &output)) goto error;

	/* Printing and free-ing */
	matrix_print_rows(output);
	return 0;

error:
	matrix_free_safe(output);
	matrix_free_safe(WAM);
	return BAD_ALLOC;
}

int print_normalized_laplacian() {
	matrix_t output;
	if (0 != graph_normalized_laplacian(datapoints, num_data, dim, &output)) {
		return BAD_ALLOC;
	}

	/* Printing and free-ing */
	matrix_print_rows(output);
	matrix_free(output);
	return 0;
}

int print_jacobi_output() {
	matrix_t jacobi_input;
	jacobi_output output;

	/* making sure that the given vectors' dataset represents a symmetric matrix (else jacobi isn't feasible) */
	assert_other(num_data == dim);
	
	/* Converting the input into a matrix and sending it into the jacobi algorithm */
	if (0 != matrix_build_from_dpoints(datapoints, num_data, dim, &jacobi_input)) goto error;
	
	if (0 != eigen_jacobi(jacobi_input, num_data + 1, &output)) goto error;

	/* Printing and free-ing */
	eigen_print_jacobi(output);
	matrix_free(jacobi_input);
	free(output.eigen_values);
	matrix_free(output.K_eigen_vectors);
	return 0;

error:
	matrix_free_safe(jacobi_input);
	if (NULL != output.eigen_values) free(output.eigen_values);
	matrix_free_safe(output.K_eigen_vectors);
	return BAD_ALLOC;
}

int get_T_of_spectral_kmeans(size_t K, matrix_t* output) {
	matrix_t L_norm;
	jacobi_output jacobi_res;
	size_t i, j;


	if (0 != graph_normalized_laplacian(datapoints, num_data, dim, &L_norm)) goto error;
	
	if (0 != eigen_jacobi(L_norm, K, &jacobi_res)) goto error;

	if (0 != matrix_new(jacobi_res.K_eigen_vectors.rows, jacobi_res.K_eigen_vectors.cols, output)) goto error;

	for (i = 0; i < output->rows; i++) {
		double norm_of_row = 0;
		for (j = 0; j < output->cols; j++) {
			norm_of_row += pow( matrix_get(jacobi_res.K_eigen_vectors, i, j), 2 );
		}
		norm_of_row = pow( norm_of_row, 0.5 );

		for (j = 0; j < output->cols; j++) {
			matrix_set(*output, i, j, matrix_get(jacobi_res.K_eigen_vectors, i, j) / norm_of_row );
		}
	}

	/* Free-ing and Returning */
	free(jacobi_res.eigen_values);
	matrix_free(jacobi_res.K_eigen_vectors);
	matrix_free(L_norm);
	return 0;

error:
	matrix_free_safe(L_norm);
	matrix_free_safe(*output);
	if (NULL != jacobi_res.eigen_values) free(jacobi_res.eigen_values);
	matrix_free_safe(jacobi_res.K_eigen_vectors);
	return BAD_ALLOC;

}

/*****************************************************************************/
