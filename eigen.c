#include "eigen.h"









matrix_t eigen_jacobi(matrix_t mat, int K) {
	matrix_t prev, next, P_rotation, P_rotation_t, P_multiplication, U_output;
	
	P_multiplication = matrix_identity(mat.rows);
	next = matrix_clone(mat);
	
	do {
		matrix_free_safe(prev);
		prev = next;
		
		P_rotation = eigen_build_rotation_matrix(prev);
		P_rotation_t = matrix_transpose(P_rotation);
		matrix_mul_assign(P_multiplication, P_rotation);

		matrix_mul(P_rotation_t, prev, &next);
		matrix_mul_assign(next, P_rotation);
				
	} while (eigen_distance_of_squared_off(prev, next) <= epsilon);
	
	matrix_free_safe(prev);
	U_output = eigen_calc_eigen_vectors(next, K);
	return U_output;
}



matrix_t eigen_calc_eigen_vectors(matrix_t mat, int K) {
	size_t K; matrix_t U_output;
	size_t i, j;
	double max_gap = -1;
	eigen* sorted_eigen_values;
	size_t eigen_col
	
	/* find K */
	sorted_eigen_values = eigen_extract_eigen_values(mat);
	if (K < 0) {
		K = eigen_heuristic_gap(sorted_eigen_values);
	}
	
	/* Form a matrix with the K-first eigen values */
	U_output = matrix_new(mat.rows, K);
	for (j = 0; j < K; j++) {
		eigen_col = sorted_eigen_values[j].col;
		
		for (i = 0; i < U_output.rows; i++) {
			matrix_set(U_output, i, j, matrix_get(mat, i, eigen_col) ); 
		}
	}
	
	return U_output;
}


size_t eigen_heuristic_gap(eigen* sorted_eigen_values) {
	for (i = 0; i < floor( mat.rows/2 ); i++) {
		double tmp;
		tmp = fabs( sorted_eigen_values[i].value - sorted_eigen_values[i+1].value )
		
		if (tmp > max_gap) {
			K = i + 1;
		}
	}
}



eigen* eigen_extract_eigen_values(matrix_t mat) {
	int i;
	eigen* eigen_values;
	
	eigen_values = malloc(mat.rows * sizeof(eigen));
	for (i = 0; i < mat.cols; i++) {
		eigen_values[i].data = matrix_get(mat, i, i);
		eigen_values[i].col = i;
	}
	
	qsort(eigen_values, mat.rows, sizeof(eigen), compare);
	return eigen_values;
}



matrix_t eigen_build_rotation_matrix(matrix_t mat) {
	matrix_ind loc;
	double *c, *s;
	matrix_t P_rotation;
	
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
	double sum;
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



int eigen_compare(const void* eigen1, const void* eigen2) {
	if ( ((eigen*)eigen1)->value < ((eigen*)eigen2)->value ) return -1;
	else if ( ((eigen*)eigen1)->value > ((eigen*)eigen2)->value ) return 1;
	
	return 0;
}





