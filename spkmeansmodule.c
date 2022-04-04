#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"





#define PY_ERROR -1




/**************************************************************************/
static PyObject* run_goal(PyObject *self, PyObject *args);
static PyObject* kmeans_fit(PyObject *self, PyObject *args);

static int matrixToList(const matrix_t mat, PyObject **output);
static int listToArray_D(PyObject *list, size_t length, double* output);
static int listToArray_L(PyObject *list, size_t length, size_t* output);
static int py_parse_args(PyObject*);

/**************************************************************************/


/**************************************************************************/
double** datapoints_arg = NULL;
size_t* initial_centroids_indices = NULL;
/**************************************************************************/



/************************* configuring the C API ****************************************************/
static PyObject* run_goal(PyObject *self, PyObject *args) {
	char* error_str;
	char* infile;
	matrix_t output;
	PyObject* T_points = NULL;
	int signal;

	/* Fetch the string of the infile */
	if(!PyArg_ParseTuple(args, "s", &infile)) {
		return NULL;
	}

	/* Perform the wanted goal's operation and return the result (if there's any) */
	if ( 0 != (signal = spkmeans_pass_goal_info_and_run(infile, &output)) ) { /* run failed */
		goto error_goal;

	} else if (strcmp(goal, "spk") == 0) { /* run succeeded. goal was "spk", therefore we have to return to python the output of the goal's mechanism */

		/* Initializing T_points as a list */
		if ( NULL == (T_points = PyList_New(K)) ) {
			PyErr_SetString(PyExc_RuntimeError, "Error: trying to create a new python list!");
			goto error_generic;
		}

		/* Buildin the T_points matrix into a python list of lists */
		if (0 != matrixToList(output, &T_points)) goto error_generic;

		return T_points;
	}

	Py_RETURN_NONE;

error_generic:
	Py_XDECREF(T_points);
	return NULL;

error_goal:
	if (signal == BAD_ALLOC) error_str = "Error: Bad Allocation!";
	if (signal == DIM_MISMATCH) error_str = "Error: Dimension Mismatch!";
	PyErr_SetString(PyExc_RuntimeError, error_str);
	return NULL;

}


static PyObject* kmeans_fit(PyObject *self, PyObject *args) {
	/* parsing the given lists as arrays (If an error has been captured
	 * a PyExc has been set, and we return NULL */
	if (0 != py_parse_args(args)) {
		PyErr_SetString(PyExc_RuntimeError, "Error: Argument parsing failed!");
		return NULL;
	}

	/* building the returned centroids' list */
	spkmeans_pass_kmeans_info_and_run(initial_centroids_indices);
	free_program();
	Py_RETURN_NONE;
}
/**************************************************************************/






/***************************** Generic C API Functions ***************************/

/* This parses the given Python arguments into C-represented Objects + manage Reference counts of Py args 
 * Returns 0 on success, and 1 on failure */
static int py_parse_args(PyObject *args) {
	size_t i;
	PyObject *datapoints_py = NULL;
	PyObject *initial_centroids_indices_py = NULL;

	/* Fetching Arguments from Python */
	if(!PyArg_ParseTuple(args, "iiidiOO", &K, &dim, &num_data, &datapoints_py, &initial_centroids_indices_py)) {
		return 1;
	}

	/* Parsing the datapoints */
	datapoints = calloc(num_data, sizeof(*datapoints));
	assert_other(NULL != datapoints);

	for (i = 0; i < num_data; i++) {
		PyObject *tmpItem;

		/* PyList_GetItem returns a borrowed reference - no need to Py_DECREF */
		if ( NULL == (tmpItem = PyList_GetItem(datapoints_py, i)) ) {
			PyErr_SetString(PyExc_RuntimeError, "Error: couldn't fetch an item from a python list!");
			goto error;
		}

		/* Appending the comprehensive array */
		init_datapoint(&datapoints[i]);
		if (0 != listToArray_D(tmpItem, dim, datapoints[i].data)) goto error;
	}

	initial_centroids_indices = calloc(K, sizeof(size_t));
	assert_other(NULL != initial_centroids_indices);
	if (0 != listToArray_L(initial_centroids_indices_py, K, initial_centroids_indices)) goto error;

	Py_XDECREF(datapoints_py);
	return 0;


error:
	/* if any of the CPython functions fail:
	 * An error of py_parse_args doesn't trigger the free_program, 
	 * since the program hasn't advanced enough, therefore we free datapoints_arg here */
	if (datapoints_arg != NULL) {
		for (i = 0; i < num_data; i++) {
			if (datapoints_arg[i] != NULL) free(datapoints_arg[i]);
		}
		free(datapoints_arg);
	}
	Py_XDECREF(datapoints_py);
	Py_XDECREF(initial_centroids_indices_py);
	return 1;

}


/* This builds a PyList out of an existing array
 * No need to care about reference counts - it's managed by fit_capi */
static int matrixToList(const matrix_t mat, PyObject **output) {
	PyObject *pyfloat = NULL;
	size_t i, j;

	/* Creating outer list */
	if ( NULL == (*output = PyList_New(mat.rows)) ) {
		PyErr_SetString(PyExc_RuntimeError, "Creation of a PyList has failed!");
		goto error;
	}

	for (i = 0; i < mat.rows; ++i) {
		PyObject* tmpList;

		/* Creating inner list */
		if ( NULL == (tmpList = PyList_New(mat.cols)) ) {
			PyErr_SetString(PyExc_RuntimeError, "Creation of a PyList has failed!");
			goto error;
		}

		/* Building inner list */
		for (j = 0; j < mat.cols; ++j) {

			if(NULL == (pyfloat = PyFloat_FromDouble( matrix_get(mat, i, j) ) ) ) {
				PyErr_SetString(PyExc_TypeError, "Double to Float conversion failed!");
				goto error;
			}
			if(0 != PyList_SetItem(tmpList, (Py_ssize_t)j, pyfloat)) {
				PyErr_SetString(PyExc_IndexError, "Index out of bounds or the index isn't an integer!");
				goto error;
			}
		}

		/* Inserting inner list */
		if(0 != PyList_SetItem(*output, (Py_ssize_t)i, tmpList)) {
			PyErr_SetString(PyExc_IndexError, "Index out of bounds or the index isn't an integer!");
			goto error;
		}
	}

	return 0;

error:
	Py_XDECREF(*output);
	Py_XDECREF(pyfloat);
	return PY_ERROR;
}




/* This parses a python Floats' List into a C Double's array
 * No need to worry about reference counts, it's managed by py_parse_args() */
static int listToArray_D(PyObject *list, size_t length, double* output) {
	size_t i;
	PyObject *pypoint = NULL;

	/* first check if the given PyObject is indeed a list */
	if (!PyList_Check(output)) {
		PyErr_SetString(PyExc_TypeError, "The passed argument isn't a list");
		goto error;
	}


	/* Insert the data into the array */
	for(i = 0; i < length; ++i) {
		/* PyList_GetItem returns a borrowed reference - no need to Py_DECREF */
		pypoint = PyList_GetItem(list, (Py_ssize_t)i);
		if (!PyFloat_Check(pypoint)) {
			PyErr_SetString(PyExc_TypeError, "Must pass an list of floats");
			goto error;
		}
		output[i] = (double) PyFloat_AsDouble(pypoint);
	}
	return 0;


error:
	/* If any of the CPython functions fail */
	Py_XDECREF(pypoint);
	return PY_ERROR;
}




/* This parses a python Integers' List into a C Long's array
 * No need to worry about reference counts, it's managed by py_parse_args() */
static int listToArray_L(PyObject *list, size_t length, size_t* output) {
	size_t i;
	PyObject *pypoint = NULL;

	/* first check if the given PyObject is indeed a list */
	if (!PyList_Check(output)) {
		PyErr_SetString(PyExc_TypeError, "The passed argument isn't a list");
		goto error;
	}


	/* Insert the data into the array */
	for(i = 0; i < length; ++i) {
		/* PyList_GetItem returns a borrowed reference - no need to Py_DECREF */
		pypoint = PyList_GetItem(list, (Py_ssize_t)i);
		if (!PyLong_Check(pypoint)) {
			PyErr_SetString(PyExc_TypeError, "Must pass an list of floats");
			goto error;
		}
		output[i] = (size_t) PyLong_AsLong(pypoint);
	}
	return 0;


error:
	/* If any of the CPython functions fail */
	Py_XDECREF(pypoint);
	return PY_ERROR;
}
/**************************************************************************/




/**************************************************************************/
static PyMethodDef capiMethods[] = {
	{"goal",
		(PyCFunction) run_goal,
		METH_VARARGS,
		PyDoc_STR("Perform the wanted operations on the given datapoints, corresponding to the determined 'goal'")
	},
	{"kmeans_fit",
		(PyCFunction) kmeans_fit,
		METH_VARARGS,
		PyDoc_STR("Given a set of datapoints, an array of the indices of the initial centroids (induced from kmeans++'s first step), perform the kmeans algorithm")
	},
	{NULL, NULL, 0, NULL}
};




static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"spkmeans",
	NULL,
	-1,
	capiMethods
};



PyMODINIT_FUNC
PyInit_mykmeanssp(void) {
	PyObject *m;
	m = PyModule_Create(&moduledef);
	if (!m) {
		return NULL;
	}
	return m;
}
/**************************************************************************/
