#include "spkmeans_goals.h"


/* Global variables crucial for preforming the goals. They're defined by the spkmeans.c file, no matter the situation, hence the 'extern' type */
extern size_t dim;
extern size_t num_data;
extern dpoint_t *datapoints;

extern void assert_other(int condition);

/************************* INTERFACE FOR GOALS *******************************/

void print_weighted_adjacency_matrix() {
	matrix_t output = graph_adjacent_matrix(datapoints, num_data, dim);
    assert_other(NULL != output.data);
	
	/* Printing and free-ing */
	matrix_print_rows(output);
	matrix_free(output);
}

void print_diagonal_degree_matrix() {
	matrix_t WAM, output;
	WAM = graph_adjacent_matrix(datapoints, num_data, dim);
    assert_other(NULL != WAM.data);

	output = graph_diagonal_degree_matrix(WAM, false);
    assert_other(NULL != output.data);
	
	/* Printing and free-ing */
	matrix_print_rows(output);
	matrix_free(output);
	matrix_free(WAM);
}

void print_normalized_laplacian() {
	matrix_t output;
	int signal = graph_normalized_laplacian(datapoints, num_data, dim, &output);
    assert_other(signal == 0);
	
	/* Printing and free-ing */
	matrix_print_rows(output);
	matrix_free(output);
}

void print_jacobi_output() {
	matrix_t jacobi_input;
	jacobi_output output;
	
	/* making sure that the given vectors' dataset represents a symmetric matrix (else jacobi isn't feasible) */
	assert_other(num_data == dim);
	jacobi_input = matrix_build(datapoints, num_data, dim);
    assert_other(jacobi_input.data != NULL);

	output = eigen_jacobi(jacobi_input, num_data + 1);
    assert_other(output.signal == 0);
	
	/* Printing and free-ing */
	eigen_print_jacobi(output);
	matrix_free(jacobi_input);
	free(output.eigen_values);
    matrix_free(output.K_eigen_vectors);
}

matrix_t get_T_of_spectral_kmeans(size_t K) {
	matrix_t L_norm, T_points;
	jacobi_output jacobi_out;
	size_t i, j;
    int signal;
	
	signal = graph_normalized_laplacian(datapoints, num_data, dim, &L_norm);
    assert_other(signal == 0);

	jacobi_out = eigen_jacobi(L_norm, K);
    assert_other(jacobi_out.signal == 0);

	
	T_points = matrix_new(jacobi_out.K_eigen_vectors.rows, jacobi_out.K_eigen_vectors.cols);
    assert_other(T_points.data != NULL);

	for (i = 0; i < T_points.rows; i++) {
		double norm_of_row = 0;
		for (j = 0; j < T_points.cols; j++) {
			norm_of_row += pow( matrix_get(T_points, i, j), 2 );
		}
		norm_of_row = pow( norm_of_row, 0.5 );
		
		for (j = 0; j < T_points.cols; j++) {
			matrix_set(T_points, i, j, matrix_get(T_points, i, j) / norm_of_row );
		}
	}
	
	/* Free-ing and Returning */
	free(jacobi_out.eigen_values);
    matrix_free(jacobi_out.K_eigen_vectors);
	matrix_free(L_norm);
	return T_points;
}



/*****************************************************************************/
