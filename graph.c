#include "graph.h"



double euclidean_norm(dpoint_t v1, dpoint_t v2, size_t dim) {
	size_t i;
	double sum = 0.0;
	
	for (i = 0; i < dim; i++) {
		dim += pow((v1.data[i] - v2.data[i]), 2);
	}	
	
	return pow(sum, 0.5);
}



matrix_t graph_adjacent_matrix(dpoint_t input[], size_t dim) {
	size_t i, j;
	matrix_t output;
	
	output = matrix_new(dim, dim);
	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			double tmp;
			tmp = exp( euclidean_norm( input[i], input[j], dim ) / (-2));
			matrix_set(output, i, j, tmp);
		}
	}
	
	return output;
}



matrix_t graph_diagonal_degree_matrix(matrix_t mat) {
	size_t i, j;
	size_t dim;
	matrix_t output;
	
	dim = mat.rows;
	output = matrix_new(mat.rows, mat.cols);
	for (i = 0; i < dim; i++) {
		double sum = 0.0, result;
		
		for (j = 0; j < dim; j++) {
			sum += matrix_get(mat, i, j);
		}
		
		result = (1 / sqrt(sum));
		matrix_set(output, i, i, result);
	}
	
	return output;
}





matrix_t graph_normalized_laplacian(dpoint_t input[], size_t dim) {
	matrix_t D, W, I, L_norm, MULT;
	
	I = matrix_identity_matrix(dim);
	W = graph_adjacent_matrix(input, dim);
	D = graph_diagonal_degree_matrix(W);
	
	matrix_mul(D, W, &MULT);
	matrix_mul_assign(MULT, D);
	matrix_mul_scalar_assign(MULT, -1);
	
	L_norm = I;
	matrix_add(L_norm, MULT, &L_norm);
	
	matrix_free_safe(W);
	matrix_free_safe(D);
	matrix_free_safe(D);
	
	return L_norm;
}
