#include "eigen.h"
#include "matrix.h"
#include <math.h>



/********************************************* STATIC FUNCTION DECLARATIONS (JACOBI's ALGORITHM) **************************************************************/
/* Given the diagonal matrix <mat_of_eigens>, pull the eigen values out of its diagonal,
 * sort them if needed, and determine how many of them we need to store in the <output> argument.
 * The determination of the amount of eigen values, is done by the value of <K>. If K==0, then we use the heuristic gap to determine a new K.
 * Once K is determined, we store the <K>-first eigen vectors into the <output> argument.
 * Most of this function's work is to simply format the output of the jacobi algorithm - And when needed, apply the heuristic gap. */
static int jacobi_format_output(matrix_t mat_vectors, matrix_t mat_of_eigens, size_t K, jacobi_t* output);

/* In case K==0 was given as input, we try to determine a new valid K using the eigen heuristic gap.
 * <eigen_values_amount> is the amount of eigen values. We need that for the heuristic's algorithm.
 * Pre-Conditions:
 *		the given array of eigen values must be sorted by value (remember: eigen is a struct) */
static size_t jacobi_eigen_heuristic(eigen_t* sorted_eigen_values, size_t eigen_values_amount);

/* Given an the last stage of A_tag in the jacobi algorithm, extract its eigen values.
   If sort equals <true>, sort the eigen values.
   An array of eigen values would be assigned to the output argument. */
static int jacobi_extract_eigen_values(matrix_t mat, bool sort, eigen_t** output);

/* In the jacobi algorithm, this is the function that transforms A_tag (the next matrix in the recursive algorithm), through the current A matrix */
static void jacobi_update_A_tag(matrix_t A_tag, matrix_t A, matrix_ind_t loc, double c, double s);

/* Given the <c> and <s> and <i,j> (in <loc>), which we are supposed to build a rotation matrix upon, simply
 * apply the changes *in-place*, that would have occurred due to a right-hand multiplication in the rotation matrix.
 * This changes will be applied to the matrix <V>.
 * Pre-condition: 
 - matrix_ind_t <loc> must be an index of the the largest off diagonal value, in the upper half of the current jacobi matrix.
 This means, that <loc> must point to an index that satisfies i < j. */
static void jacobi_apply_rotation(matrix_t V, matrix_ind_t loc, double c, double s);

/* Calculate the rotation matrix using the given data, and store the result in the pre-allocated `output`. */
int eigen_build_rotation_matrix(matrix_ind_t loc, double c, double s, matrix_t output);

/* Given two matrices, determine the distance between their sum of squared off-diagonals */
static double jacobi_distance_of_squared_offdiagonals(matrix_t mat1, matrix_t mat2);

/* Calculae the values of 'c' and 's' of the desired rotation matrix */
static void jacobi_calc_c_s(double* c, double *s, matrix_t current_jacobi_mat, matrix_ind_t loc);
/*****************************************************************************************************************************************/







/********************************************* GLOBAL FUNCTIONS OF THE EIGEN MODULE **************************************************************/
void eigen_print_jacobi(jacobi_t out) {
	size_t i;

	for (i = 0; i < out.eigen_vectors.cols; i++) {

		/* If an eigen value will be printed as -0.0000, we'll replace the printage with 0.0000 */
		if ( (out.eigen_values[i].value > -0.0001) && (out.eigen_values[i].value < 0) ) {
			printf("%.4f", (double)0);
		} else {
			printf("%.4f", out.eigen_values[i].value);
		}

		/* Segregation of eigen values */
		if (i < out.eigen_vectors.cols - 1) printf(",");
	}

	puts("");
	matrix_print_rows(out.eigen_vectors);
}



int eigen_jacobi(matrix_t mat, size_t K, jacobi_t* output) {
	size_t iterations;
	matrix_t A, A_tag, V;
	matrix_ind_t loc;
	double s, c;

	A.data = NULL;
	A_tag.data = NULL;
	V.data = NULL;

	if (matrix_identity(mat.rows, &V)) goto error;
	if (matrix_clone(mat, &A_tag)) goto error;
	if (matrix_clone(mat, &A)) goto error;

	for(iterations = 0; iterations < max_jacobi_iterations; iterations++) {
		matrix_copy(A, A_tag); /* matrices are created with equal dims - no error check */

		loc = matrix_ind_of_largest_offdiagonal(A); 
		if (0 == matrix_get(A, loc.i, loc.j)) break;/* stop the algorithm if the matrix of the last iteration is diagonal (the next step will result in nan-s) */

		jacobi_calc_c_s(&c, &s, A, loc);
		jacobi_apply_rotation(V, loc, c, s); /* in-place multiplication of the rotation matrix of the current iteration and V (the output eigen vector matrix) */
		jacobi_update_A_tag(A_tag, A, loc, c, s); /* updating the current matrix of the jacobi algorith into the new matrix of the next iteration */

		if (jacobi_distance_of_squared_offdiagonals(A, A_tag) <= epsilon) break;

	}

	/* Extract the eigen values and eigen vectors and insert them into an output format */
	if (jacobi_format_output(V, A_tag, K, output)) goto error;

	/* If  K < num_data, then goal = spk, and we have to sort the eigen values, 
	 * therefore creating another matrix for the eigen vectors - hence, we don't need V anymore */
	if (K < V.cols) { 
		matrix_free(V);
	} 

	/* Free-ing */
	matrix_free(A);
	matrix_free(A_tag);
	return 0;

error:
	matrix_free_safe(A);
	matrix_free_safe(A_tag);
	matrix_free_safe(V);

	return BAD_ALLOC;
}



int eigen_jacobi_to_mat(jacobi_t origin, matrix_t *output) {
	size_t i, j;

	/* Creating the output matrix */
	if (matrix_new(origin.eigen_vectors.rows + 1, origin.eigen_vectors.cols, output)) {
		return BAD_ALLOC;
	}


	/* Building the first row of the matrix (the eigen values) */
	for (j = 0; j < output->cols; j++) {
		double value;

		/* Building the eigen value that will be inserted into the output matrix */
		if ( (origin.eigen_values[j].value > -0.0001) && (origin.eigen_values[j].value < 0) ) { /* If we need to round up an eigen value */
			value = 0.0;
		} else {
			value = origin.eigen_values[j].value;
		}

		/* Inserting the eigen value into the output matrix */
		matrix_set(*output, 0, j, value);
	}


	/* Building the rest of the rows of the matrix (the eigen vectors) */
	for (i = 1; i < output->rows; i++) {
		for (j = 0; j < output->cols; j++) {
			matrix_set(*output, i, j, matrix_get(origin.eigen_vectors, i - 1, j));
		}
	}

	/* Free the original output format */
	free(origin.eigen_values);
	matrix_free(origin.eigen_vectors);

	return 0;
}



int eigen_compare(const void* eigen1, const void* eigen2) {
	double eigen1_val = ((eigen_t*)eigen1)->value;
	double eigen2_val = ((eigen_t*)eigen2)->value;

	return (eigen1_val > eigen2_val) ? 1 : ( (eigen1_val < eigen2_val) ? -1 : 0 );
}



int sign(double val) {
	if (val >= 0) return 1;
	return -1;
}
/*****************************************************************************************************************************************/
















/********************************************* STATIC FUNCTION DEFINITIONS (RELATED TO JACOBI's ALGORITHM) **************************************************************/
static int jacobi_format_output(matrix_t mat_vectors, matrix_t mat_of_eigens, size_t K, jacobi_t* output) {
	size_t i, j;
	eigen_t* sorted_eigen_values = NULL;
	matrix_t eigen_vectors;


	if (K < mat_vectors.cols) {	/* The goal was spk, since K == mat_vectors.cols is prohibited by the Python CMD interface */
		/* sort the eigen values */
		if (jacobi_extract_eigen_values(mat_of_eigens, true, &sorted_eigen_values)) goto error;

		/* If K == 0, it means the CMD asked us to use the heuristic gap to determine K */
		if (K == 0) { 
			K = jacobi_eigen_heuristic(sorted_eigen_values, mat_of_eigens.rows);
		}

		/* Form a matrix with the K-first eigen values */
		if (matrix_new(mat_vectors.rows, K, &eigen_vectors)) goto error;

		/* For each eigen value out of the first K, find its corresponding column in V, and copy it into <eigen_vectors> */
		for (j = 0; j < K; j++) {
			size_t eigen_col = sorted_eigen_values[j].col;

			for (i = 0; i < eigen_vectors.rows; i++) {
				matrix_set(eigen_vectors, i, j, matrix_get(mat_vectors, i, eigen_col) ); 
			}
		}

		/* Format the output */
		output->eigen_vectors = eigen_vectors;
		output->eigen_values = sorted_eigen_values;

	} else if (K == mat_vectors.cols) {
		/* If K = mat_vectors.cols, it means that the jacobi algorithm was powered alone.
		 * That is, since K = mat_vectors.cols is prohibited by the python CMD interface.
		 * Moreover, it's since we use this case as an indicator to when jacobi was powered without any future spectral clustering use. 
		 * In such case, a jacobi algorithm alone isn't due to any specification of K, and we will return all of the eigen values/vectors (unsorted) */

		/* Extracting all of the eigen values, without sorting */
		if (jacobi_extract_eigen_values(mat_of_eigens, false, &output->eigen_values)) goto error;

		/* Extracting all of the eigen vectors */
		output->eigen_vectors = mat_vectors;
	}

	return 0;

error:
	if(sorted_eigen_values) {
		free(sorted_eigen_values);
	}
	return BAD_ALLOC;
}


static size_t jacobi_eigen_heuristic(eigen_t* sorted_eigen_values, size_t size) {
	size_t i, K = size;
	double tmp;
	double max_gap = -1;

	for (i = 0; i < size/2; i++) {
		tmp = fabs( sorted_eigen_values[i].value - sorted_eigen_values[i+1].value );

		if (tmp > max_gap) {
			max_gap = tmp;
			K = i + 1;
		}
	}

	return K;
}


static int jacobi_extract_eigen_values(matrix_t mat, bool sort, eigen_t** output) {
	size_t i;

	*output = malloc(mat.rows * sizeof(eigen_t));
	if(NULL == *output) return BAD_ALLOC;

	for (i = 0; i < mat.cols; i++) {
		(*output)[i].value = matrix_get(mat, i, i);
		(*output)[i].col = i;
	}

	if (sort) {
		qsort(*output, mat.cols, sizeof(eigen_t), eigen_compare);
	}

	return 0;	
}


static void jacobi_update_A_tag(matrix_t A_tag, matrix_t A, matrix_ind_t loc, double c, double s) {
	double c2, s2, Aii, Ajj, Aij;
	size_t i, j, r;

	i = loc.i;
	j = loc.j;

	for (r = 0; r < A.rows; r++) {	
		if (r != i && r != j) {
			double a_ri, a_rj;

			a_ri = matrix_get(A, r, i);
			a_rj = matrix_get(A, r, j);

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


static void jacobi_apply_rotation(matrix_t V, matrix_ind_t loc, double c, double s) {
	size_t row;

	/* P is essentially an identity matrix with 4 values changed. No need for a robust matrix multiplication algorithm.
	 * We'll change just the values in V that are supposed to be changed due to such multiplication, accordingly */
	for (row = 0; row < V.rows; row++) {
		double row_i, row_j;
		row_i = c * matrix_get(V, row, loc.i) - s * matrix_get(V, row, loc.j);
		row_j = s * matrix_get(V, row, loc.i) + c * matrix_get(V, row, loc.j);
		matrix_set(V, row, loc.i, row_i);
		matrix_set(V, row, loc.j, row_j);
	}
}



static double jacobi_distance_of_squared_offdiagonals(matrix_t mat1, matrix_t mat2) {
	return matrix_sum_squared_off(mat1) - matrix_sum_squared_off(mat2);
}



static void jacobi_calc_c_s(double* c, double *s, matrix_t current_jacobi_mat, matrix_ind_t loc) {
	double theta, tmp;

	theta = (matrix_get(current_jacobi_mat, loc.j, loc.j) - matrix_get(current_jacobi_mat, loc.i, loc.i)) / (2 * matrix_get(current_jacobi_mat, loc.i, loc.j));
	tmp = sign(theta) / ( fabs(theta) + sqrt( pow(theta, 2) + 1 ) );
	*c = 1 / sqrt( pow(tmp, 2) + 1 );
	*s = tmp * (*c);
}
/*****************************************************************************************************************************************/



