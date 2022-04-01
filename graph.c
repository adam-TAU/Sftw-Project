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



matrix_t graph_adjacent_matrix(dpoint_t input[], size_t num_data, size_t dim) {
	size_t i, j;
	matrix_t output = matrix_new(num_data, num_data);

    if(NULL == output.data) return output;

	for (i = 0; i < num_data; i++) {
        /* starting with j = i + 1 ensures none of the diagonal values are nonzero */
		for (j = i + 1; j < num_data; j++) {
			double tmp = exp( euclidean_norm( input[i], input[j], dim ) * (-0.5) );
			matrix_set(output, i, j, tmp);
			matrix_set(output, j, i, tmp);
		}
	}
	
	return output;
}



matrix_t graph_diagonal_degree_matrix(matrix_t mat, bool is_sqrt) {
	size_t i, j;
	size_t dimension = mat.rows;
	matrix_t output = matrix_new(mat.rows, mat.cols);
	
	if(NULL == output.data) return output;
    
	for (i = 0; i < dimension; i++) {
		double sum = 0.0, result;
		
		for (j = 0; j < dimension; j++) {
			sum += matrix_get(mat, i, j);
		}
		
		result = is_sqrt ? (1 / sqrt(sum)) : sum;
		matrix_set(output, i, i, result);
	}
	
	return output;
}





int graph_normalized_laplacian(dpoint_t input[], size_t num_data, size_t dim, matrix_t *output) {
	matrix_t D, W, I, MULT;

    /* in case of an error */
    D.data = NULL;
    W.data = NULL;
    I.data = NULL;
    MULT.data = NULL;
	
	I = matrix_identity(num_data);
    if(NULL == I.data) goto error;

	W = graph_adjacent_matrix(input, num_data, dim);
    if(NULL == W.data) goto error;

	D = graph_diagonal_degree_matrix(W, true);
    if(NULL == D.data) goto error;
	
	if(0 != matrix_mul(D, W, &MULT)) goto error;
	if(0 != matrix_mul_assign_to_first(&MULT, D)) goto error;
	matrix_mul_scalar_assign(MULT, -1);
	
	*output = I;
	matrix_add_assign(*output, MULT);

    matrix_free(MULT);
    matrix_free(W);
    matrix_free(D);
    /* do not free `I`, because it's the output */
	
	return 0;

error:
    matrix_free_safe(MULT);
	matrix_free_safe(W);
	matrix_free_safe(D);
    matrix_free_safe(I);
    return BAD_ALLOC;
}
