#include "eigen.h"




void eigen_free_jacobi_safe(jacobi_output out) {
	if (out.eigen_values != NULL) free(out.eigen_values);
	matrix_free_safe(out.K_eigen_vectors);
}


void eigen_print_jacobi(jacobi_output out) {
	size_t i;
	
	for (i = 0; i < out.K_eigen_vectors.cols; i++) {
		printf("%.4f", out.eigen_values[i].value);
		if (i < out.K_eigen_vectors.cols - 1) printf(", ");
	}
	
	puts("");
	matrix_print_cols(out.K_eigen_vectors);
}



jacobi_output eigen_jacobi(matrix_t mat) {
	matrix_t prev, next, P_rotation, P_rotation_t, P_multiplication;
	jacobi_output jacobi_output;

	P_multiplication = matrix_identity_matrix(mat.rows);
	next = matrix_clone(mat);
	
	do {
		matrix_free_safe(prev);
		prev = next;
		
		P_rotation = eigen_build_rotation_matrix(prev);
		P_rotation_t = matrix_transpose(P_rotation);
		matrix_mul_assign(P_multiplication, P_rotation);

		matrix_mul(P_rotation_t, prev, &next);
		matrix_mul_assign(next, P_rotation);
		
		matrix_free_safe(P_rotation);
	} while (eigen_distance_of_squared_off(prev, next) <= epsilon);
	
	/* Simply sort the the eigen vectors and eigen values */
	jacobi_output = eigen_format_eigen_vectors(P_multiplication, next, P_multiplication.cols);
	
	matrix_free_safe(prev);
	matrix_free_safe(next);
	matrix_free_safe(P_rotation);
	matrix_free_safe(P_multiplication);
	
	return jacobi_output;
}



jacobi_output eigen_format_eigen_vectors(matrix_t mat_vectors, matrix_t mat_eigens, size_t K) {
	matrix_t U_eigen_vectors;
	size_t i, j;
	eigen* sorted_eigen_values;
	size_t eigen_col;
	jacobi_output result;
	
	/* find K */
	sorted_eigen_values = eigen_extract_eigen_values(mat_eigens);
	if (K == 0 || K > U_eigen_vectors.cols) { /* this tells us that K wasn't determined through the command line interface. also a safety mechanism */
		K = eigen_heuristic_gap(sorted_eigen_values, mat_eigens.rows);
	}
	
	/* Form a matrix with the K-first eigen values */
	U_eigen_vectors = matrix_new(mat_vectors.rows, K);
	for (j = 0; j < K; j++) {
		eigen_col = sorted_eigen_values[j].col;
		
		for (i = 0; i < U_eigen_vectors.rows; i++) {
			matrix_set(U_eigen_vectors, i, j, matrix_get(mat_vectors, i, eigen_col) ); 
		}
	}
	
	result.K_eigen_vectors = U_eigen_vectors;
	result.eigen_values = sorted_eigen_values;
	return result;
}


size_t eigen_heuristic_gap(eigen* sorted_eigen_values, size_t rows) {
	size_t i, K = rows;
	double tmp;
	double max_gap = -1;
	
	for (i = 0; i < floor( rows/2 ); i++) {
		tmp = fabs( sorted_eigen_values[i].value - sorted_eigen_values[i+1].value );
		
		if (tmp > max_gap) {
			K = i + 1;
		}
	}
	
	return K;
}



int eigen_compare(const void* eigen1, const void* eigen2) {
	if ( ((eigen*)eigen1)->value < ((eigen*)eigen2)->value ) return -1;
	else if ( ((eigen*)eigen1)->value > ((eigen*)eigen2)->value ) return 1;
	
	return 0;
}



eigen* eigen_extract_eigen_values(matrix_t mat) {
	size_t i;
	eigen* eigen_values;
	
	eigen_values = malloc(mat.rows * sizeof(eigen));
	for (i = 0; i < mat.cols; i++) {
		eigen_values[i].value = matrix_get(mat, i, i);
		eigen_values[i].col = i;
	}
	
	qsort(eigen_values, mat.rows, sizeof(eigen), eigen_compare);
	return eigen_values;
}



matrix_t eigen_build_rotation_matrix(matrix_t mat) {
	matrix_ind loc;
	double *c, *s;
	matrix_t P_rotation;
	
	
	c = malloc(sizeof(double));
	s = malloc(sizeof(double));
	
	loc = eigen_ind_of_largest_offdiagonal(mat);
	P_rotation = matrix_identity_matrix(mat.rows);
	eigen_calc_c_s(c, s, mat, loc);
	
	
	matrix_set(P_rotation, loc.i, loc.i, *c); 
	matrix_set(P_rotation, loc.j, loc.j, *c);
	
	if (loc.i < loc.j) {
		matrix_set(P_rotation, loc.i, loc.j, *s);  
		matrix_set(P_rotation, loc.j, loc.i, -(*s));
	}
	
	if (loc.i > loc.j) {
		matrix_set(P_rotation, loc.j, loc.i, *s);  
		matrix_set(P_rotation, loc.i, loc.j, -(*s));
	}
	
	return P_rotation;
}




double eigen_distance_of_squared_off(matrix_t mat1, matrix_t mat2) {
	return eigen_sum_squared_off(mat1) - eigen_sum_squared_off(mat2);
}



double eigen_sum_squared_off(matrix_t mat) {
	double sum = 0.0;
	size_t i, j;
	
	for (i = 0; i < mat.rows; i++) {
		for (j = 0; j < mat.cols; j++) {
			sum += (i != j) ? matrix_get(mat, i, j) : 0;
		}
	}
	
	return sum;
}




matrix_ind eigen_ind_of_largest_offdiagonal(matrix_t mat) {
	matrix_ind output;
	size_t i, j;
	double current_max = -1;
	
	for (i = 0; i < mat.rows; i++) {
		for (j = 0; j > mat.cols; j++) {
			if (fabs(matrix_get(mat, i, j)) > current_max) {
				output.i = i;
				output.j = j;
			}
		}
	}
	
	return output;
}




void eigen_calc_c_s(double* c, double *s, matrix_t mat, matrix_ind loc) {
	double theta, tmp;
	
	theta = (matrix_get(mat, loc.j, loc.j) - matrix_get(mat, loc.i, loc.i)) / (2 * matrix_get(mat, loc.i, loc.j));
	tmp = sign(theta) / ( fabs(theta) + sqrt( pow(theta, 2) + 1 ) );
	*c = 1 / sqrt( pow(tmp, 2) + 1 );
	*s = tmp * (*c);
}



int sign(double val) {
	if (val >= 0) return 1;
	return 0;
}









