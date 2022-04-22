#include "graph.h"
#include "matrix.h"


/* Find the Euclidean norm of two vectors of the same dim x 1 */
static double euclidean_norm(dpoint_t v1, dpoint_t v2, size_t dim) {
	size_t i;
	double sum = 0.0;

	for (i = 0; i < dim; i++) {
		sum += pow((v1.data[i] - v2.data[i]), 2);
	}	

	return sqrt(sum);
}



int graph_adjacent_matrix(dpoint_t input[], size_t num_data, size_t dim, matrix_t* output) {
	size_t i, j;

	/* Creating the output matrix */
	if (matrix_new(num_data, num_data, output) != 0) goto error;

	/* Building the output matrix */
	for (i = 0; i < num_data; i++) {
		/* starting with j = i + 1 ensures none of the diagonal values are nonzero */
		for (j = i + 1; j < num_data; j++) {
			double tmp = exp( euclidean_norm( input[i], input[j], dim ) * (-0.5) );
			matrix_set(*output, i, j, tmp);
			matrix_set(*output, j, i, tmp);
		}
	}
	return 0;

error:
	/* Free-ing */
	matrix_free_safe(*output);
	return BAD_ALLOC;
}



int graph_diagonal_degree_matrix(dpoint_t input[], size_t num_data, size_t dim, bool is_sqrt, matrix_t* output) {
	size_t i, j;
	matrix_t W;

	/* in case of an error */
	W.data = NULL;
	output->data = NULL;

	/* Build the WAM matrix */
	if (graph_adjacent_matrix(input, num_data, dim, &W)) goto error;

	/* Create the output matrix */
	if (matrix_new(W.rows, W.cols, output)) goto error;

	/* Build the output matrix */
	for (i = 0; i < W.rows; i++) {
		double sum = 0.0;

		/* Calculate row-specific sum */
		for (j = 0; j < W.cols; j++) {
			sum += matrix_get(W, i, j);
		}

		/* Insert the value into the diagonal degree matrix. Apply a square root in case we want D^(-1/2) */
		matrix_set(*output, i, i, (is_sqrt ? (1 / sqrt(sum)) : sum));
	}

	/* Free-ing */
	matrix_free_safe(W);

	return 0;

error:
	/* Free-ing */
	matrix_free_safe(*output);
	return BAD_ALLOC;
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

	/* Build D_sqrt manually. There's no need to create a matrix of size n^2 just to know it's diagonal (which is n values) */
	D_sqrt = (double*) calloc(W.rows, sizeof(double));
	for (i = 0; i < W.rows; i++) {
		double sum = 0.0;
		for (j = 0; j < W.rows; j++) {
			sum += matrix_get(W, i, j);
		}

		D_sqrt[i] = 1 / sqrt(sum);
	}

	/* Creating the output matrix */
	if (matrix_new(W.rows, W.rows, output)) goto error;

	/* Building the output matrix */
	for (i = 0; i < output->rows; i++) {
		for (j = 0; j < output->cols; j++) {
			/* L_norm = I - D_sqrt * WAM * D_sqrt */
			double mult_val;

			/* Matrix multiplication avoided by in-place editing of the output matrix.
			 * This is avoided since D_sqrt is always a diagonal matrix, and thus we represent D_sqrt's diagonal as an array */
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
	matrix_free_safe(*output);
	return BAD_ALLOC;
}


