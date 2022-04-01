#include "eigen.h"
#include "matrix.h"





void eigen_print_jacobi(jacobi_output out) {
	size_t i;
	
	for (i = 0; i < out.K_eigen_vectors.cols; i++) {
		printf("%.4f", out.eigen_values[i].value);
		if (i < out.K_eigen_vectors.cols - 1) printf(",");
	}
	
	puts("");
	matrix_print_cols(out.K_eigen_vectors);
}



jacobi_output eigen_jacobi(matrix_t mat, size_t K) {
	size_t iterations = 0;
	matrix_t A, A_tag, P, V;
	jacobi_output jacobi_output;
    matrix_ind loc;
    double s, c;

	V = matrix_identity(mat.rows);
	A_tag = matrix_clone(mat);
	A = matrix_clone(mat);
    P.data = NULL;

    if(NULL == V.data || NULL == A_tag.data || NULL == A.data) goto error;
	
	do {
		iterations++;
	
		matrix_swap(&A, &A_tag);
        loc = eigen_ind_of_largest_offdiagonal(A);
        eigen_calc_c_s(&c, &s, A, loc);

		P = eigen_build_rotation_matrix(A, loc, c, s);
        if(NULL == P.data) goto error;

		matrix_mul_assign_to_first(&V, P);
        if(NULL == V.data) goto error;
				
		eigen_update_jacobi_A_tag(A_tag, A, loc, c, s);
		matrix_free(P);
		
	} while (( eigen_distance_of_squared_offdiagonals(A, A_tag) > epsilon ) && ( iterations <= 100 ));
	
	/* Extract the eigen values and eigen vectors and insert them into an output format */
	jacobi_output = eigen_format_eigen_vectors(V, A_tag, K);
    if(jacobi_output.signal != 0) goto error;
	
	/* In this case we had to create another matrix to hold the wanted eigen vectors */
	if (K <= V.cols) { 
		matrix_free(V);
	} 
	
	matrix_free(A);
	matrix_free(A_tag);
	return jacobi_output;

error:
    matrix_free_safe(A);
    matrix_free_safe(A_tag);
    matrix_free_safe(P);
    matrix_free_safe(V);
    jacobi_output.signal = BAD_ALLOC;
    return jacobi_output;
}



jacobi_output eigen_format_eigen_vectors(matrix_t mat_vectors, matrix_t mat_eigens, size_t K) {
	size_t i, j;
	eigen* sorted_eigen_values = NULL;
	matrix_t U_eigen_vectors;
	jacobi_output result;
	
	
	if (K <= mat_vectors.cols) {	
		/* find K */
		sorted_eigen_values = eigen_extract_eigen_values(mat_eigens, true);
        if(NULL == sorted_eigen_values) goto error;
		
		/* If K == 0, it means the CMD asked us to use the heuristic gap to determine K */
		if (K == 0) { 
			K = eigen_heuristic_gap(sorted_eigen_values, mat_eigens.rows);
		}
		
		/* Form a matrix with the K-first eigen values */
		U_eigen_vectors = matrix_new(mat_vectors.rows, K);
        if(NULL == U_eigen_vectors.data) goto error;

		for (j = 0; j < K; j++) {
			size_t eigen_col = sorted_eigen_values[j].col;
			
			for (i = 0; i < U_eigen_vectors.rows; i++) {
				matrix_set(U_eigen_vectors, i, j, matrix_get(mat_vectors, i, eigen_col) ); 
			}
		}
		
		/* Format the output */
		result.K_eigen_vectors = U_eigen_vectors;
		result.eigen_values = sorted_eigen_values;
		
	} else {
		/* If K > mat_vectors.cols, it means that the jacobi algorithm was powered alone.
		 * That i, since K > mat_vectors.cols is prohibited by the python CMD interface.
		 * Moreover, it's since we use this case as an indicator to when jacobi was powered without any future spectral clustering use. 
		 * In such case, a jacobi algorithm alone isn't due to any specification of K, and we will return all of the eigen values/vectors (unsorted) */
		result.K_eigen_vectors = mat_vectors;
		result.eigen_values = eigen_extract_eigen_values(mat_eigens, false);
        if(NULL == result.eigen_values) goto error;
	}
	
	return result;

error:
    if(sorted_eigen_values) {
        free(sorted_eigen_values);
    }
    result.signal = BAD_ALLOC;
    return result;
}


size_t eigen_heuristic_gap(eigen* sorted_eigen_values, size_t rows) {
	size_t i, K = rows;
	double tmp;
	double max_gap = -1;
	
	for (i = 0; i < rows/2; i++) {
		tmp = fabs( sorted_eigen_values[i].value - sorted_eigen_values[i+1].value );
		
		if (tmp > max_gap) {
			K = i + 1;
		}
	}
	
	return K;
}



int eigen_compare(const void* eigen1, const void* eigen2) {
    return ((eigen*)eigen1)->value - ((eigen*)eigen2)->value;
}



eigen* eigen_extract_eigen_values(matrix_t mat, bool sort) {
	size_t i;
	eigen* eigen_values;
	
	eigen_values = malloc(mat.rows * sizeof(eigen));
    if(NULL == eigen_values) return NULL;

	for (i = 0; i < mat.cols; i++) {
		eigen_values[i].value = matrix_get(mat, i, i);
		eigen_values[i].col = i;
	}
	
	if (sort) {
		qsort(eigen_values, mat.rows, sizeof(eigen), eigen_compare);
	}
	
	return eigen_values;
}


void eigen_update_jacobi_A_tag(matrix_t A_tag, matrix_t A, matrix_ind loc, double c, double s) {
	double c2, s2, Aii, Ajj, Aij;
	size_t i, j, r;
	
	i = loc.i;
	j = loc.j;
	
	for (r = 0; r < A.rows; r++) {
		double a_ri, a_rj;
		
		a_ri = matrix_get(A, r, i);
		a_rj = matrix_get(A, r, j);
		
		if (r != i && r != j) {
			matrix_set(A_tag, r, i, c * a_ri - s * a_rj);
			matrix_set(A_tag, i, r, c * a_ri - s * a_rj);
			matrix_set(A_tag, r, j, c * a_rj + s * a_ri);
			matrix_set(A_tag, j, r, c * a_rj + s * a_ri);
		}
	}
	
	matrix_set(A_tag, i, j, 0);
	matrix_set(A_tag, j, i, 0);

    c2 = c*c;
    s2 = s*s;
    Aii = matrix_get(A, i, i);
    Ajj = matrix_get(A, j, j);
    Aij = matrix_get(A, i, j);
	matrix_set(A_tag, i, i, c2 * Aii + s2 * Ajj - 2 * c * s * Aij);
	matrix_set(A_tag, j, j, s2 * Aii + c2 * Ajj + 2 * c * s * Aij);
}



matrix_t eigen_build_rotation_matrix(matrix_t mat, matrix_ind loc, double c, double s) {
    matrix_t P_rotation;
	size_t i, j;
	
	i = loc.i;
	j = loc.j;
	
	P_rotation = matrix_identity(mat.rows);
    if(NULL == P_rotation.data) return P_rotation;
	
	matrix_set(P_rotation, i, i, c); 
	matrix_set(P_rotation, j, j, c);

    s *= (i < j) ? 1 : -1;
	matrix_set(P_rotation, i, j, s);  
	matrix_set(P_rotation, j, i, -s);

	
	return P_rotation;
}




double eigen_distance_of_squared_offdiagonals(matrix_t mat1, matrix_t mat2) {
	return eigen_sum_squared_off(mat1) - eigen_sum_squared_off(mat2);
}



double eigen_sum_squared_off(matrix_t mat) {
	double sum = 0.0;
	size_t i, j;
	
	for (i = 0; i < mat.rows; i++) {
		for (j = 0; j < mat.cols; j++) {
			sum += (i != j) ? pow( matrix_get(mat, i, j), 2 ) : 0;
		}
	}
	
	return sum;
}




matrix_ind eigen_ind_of_largest_offdiagonal(matrix_t mat) {
	matrix_ind output;
	size_t i, j;
	double current_max = -1;
	
	for (i = 0; i < mat.rows; i++) {
		for (j = i + 1; j < mat.cols; j++) { /* The jacobi algorithm is due to get a symmetric matrix, therefore we check only half of the off_diagonals */
			double tmp = fabs(matrix_get(mat, i, j));
			if (tmp > current_max) {
				output.i = i;
				output.j = j;
				current_max = tmp;
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
	return -1;
}
