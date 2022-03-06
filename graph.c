#include "graph.h"



double graph_euclidena_norm(double v1[], double v2[], int dim) {
	int i;
	double sum = 0.0;
	
	for (i = 0; i < dim; i++) {
		dim += pow((v1[i] - v2[i]), 2);
	}	
	
	return pow(sum, 0.5);
}




matrix_t graph_identity_matrix(size_t dim) {
	int i;
	matrix_t output;
	output = matrix_new(rows, cols);
 
 	for (i = 0; i < dim; i++) {
		matrix_set(output, i, i, 1);
 	}
}



matrix_t graph_adjacent_matrix(double input[][], int dim) {
	int i, j;
	matrix_t output;
	
	output = matrix_new(dim, dim);
	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			double tmp;
			tmp = exp( graph_euclidean_norm( input[i], input[j] ) / (-2));
			matrix_set(output, tmp);
		}
	}
}



matrix_t graph_diagonal_degree_matrix(matrix_t mat) {
	int i, j;
	size_t dim;
	matrix_t output;
	
	dim = mat.rows;
	output = matrix_new(mat.rows, mat.cols);
	for (i = 0; i < dim; i++) {
		double sum, result;
		
		for (j = 0; j < dim; j++) {
			sum += matrix_get(mat, i, j);
		}
		
		result = (1 / sqrt(sum));
		matrix_set(output, i, i, result);
	}
}





matrix_t graph_normalized_laplacian(double input[][], int dim) {
	matrix_t D, W, I, L_norm, MULT;
	
	I = graph_identity_matrix(dim);
	W = graph_adjacent_matrix(double input[][], dim);
	D = graph_diagonal_degree_matrix(W);
	
	matrix_mul(D, W, &MULT);
	matrix_mul_assign(MULT, D);
	matrix_mul_scalar_assign(MULT, -1);
	
	L_norm = I;
	matrix_add(L_norm, MULT, &L_norm);
	
	return L_norm;
}
