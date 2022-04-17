#include "graph.h"
#include "matrix.h"



double euclidean_norm(dpoint_t v1, dpoint_t v2, size_t dim) {
	size_t i;
	double sum = 0.0;

	for (i = 0; i < dim; i++) {
		sum += pow((v1.data[i] - v2.data[i]), 2);
	}	

	return sqrt(sum);
}



int graph_adjacent_matrix(dpoint_t input[], size_t num_data, size_t dim, matrix_t* output) {
	size_t i, j;
	if (matrix_new(num_data, num_data, output) != 0) {
		return BAD_ALLOC;
	}


	for (i = 0; i < num_data; i++) {
		/* starting with j = i + 1 ensures none of the diagonal values are nonzero */
		for (j = i + 1; j < num_data; j++) {
			double tmp = exp( euclidean_norm( input[i], input[j], dim ) * (-0.5) );
			matrix_set(*output, i, j, tmp);
			matrix_set(*output, j, i, tmp);
		}
	}
	return 0;
}



int graph_diagonal_degree_matrix(matrix_t mat, bool is_sqrt, matrix_t* output) {
	size_t i, j;
	size_t dimension = mat.rows;

	if(matrix_new(mat.rows, mat.cols, output) != 0) {
		return BAD_ALLOC;
	}

	for (i = 0; i < dimension; i++) {
		double sum = 0.0, result;

		for (j = 0; j < dimension; j++) {
			sum += matrix_get(mat, i, j);
		}

		result = is_sqrt ? (1 / sqrt(sum)) : sum;
		matrix_set(*output, i, i, result);
	}

	return 0;
}




int graph_normalized_laplacian(dpoint_t input[], size_t num_data, size_t dim, matrix_t *output) {
	matrix_t W;
	double* D_sqrt;
	size_t i, j;

	/* in case of an error */
	W.data = NULL;
	D_sqrt = NULL;
	output->data = NULL;
	
	/* Build the WAM matrix */
	if (graph_adjacent_matrix(input, num_data, dim, &W)) goto error;
	
	/* Create an array that contains at index <i> the sum of row <i> in matrix WAM <W> */
	D_sqrt = (double*) calloc(W.rows, sizeof(double));
	for (i = 0; i < W.rows; i++) {
		double sum = 0.0;
		for (j = 0; j < W.rows; j++) {
			sum += matrix_get(W, i, j);
		}
		
		D_sqrt[i] = 1 / sqrt(sum);
	}

	/* Build the output */
	matrix_new(W.rows, W.rows, output);
	for (i = 0; i < output->rows; i++) {
		for (j = 0; j < output->cols; j++) {
			/* Lnorm = I - D_sqrt * WAM * D_sqrt */
			double mult_val;
			mult_val = D_sqrt[i] * matrix_get(W, i, j) * D_sqrt[j];
			matrix_set(*output, i, j, (i == j) * 1 - mult_val); 
		}
	}

	/* Free-ing */
	if (D_sqrt != NULL) free(D_sqrt);
	matrix_free(W);
	return 0;

error:
	/* Free-ing */
	if (D_sqrt != NULL) free(D_sqrt);
	matrix_free_safe(W);
	return BAD_ALLOC;
}


