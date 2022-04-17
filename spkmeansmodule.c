#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"





#define PY_ERROR -1




/**************************************************************************/
static PyObject* run_goal(PyObject *self, PyObject *args);
static PyObject* kmeans_fit(PyObject *self, PyObject *args);

static int matrixToList(const matrix_t mat, PyObject **output);
static int extract_datapoint_indices(PyObject **output);
static int listToArray_D(PyObject *list, size_t length, double** output);
static int listToArray_L(PyObject *list, size_t length, size_t** output);
static int py_kmeans_parse_args(PyObject*);

/**************************************************************************/


/**************************************************************************/
static double** datapoints_arg = NULL;
static size_t* initial_centroids_indices = NULL;
/**************************************************************************/



/************************* configuring the C API ****************************************************/
static PyObject* run_goal(PyObject *self, PyObject *args) {
	char* infile;
	matrix_t output;
	PyObject* T_points = NULL;
	int signal;

    output.data = NULL; // in case of an error, `output`'s data field is freed if it's not null

	/* Fetch the string of the infile */
	if(!PyArg_ParseTuple(args, "lss", &K, &goal, &infile)) goto error;

	/* Perform the wanted goal's operation and return the result (if there's any) */
	if ( (signal = spkmeans_pass_goal_info_and_run(infile, &output)) ) { /* run failed */
		goto error;

	} else if (strcmp(goal, "spk") == 0) { /* run succeeded. goal was "spk", therefore we have to return to python the output of the goal's mechanism */
		
		/* Initializing T_points as a list */
		if ( NULL == (T_points = PyList_New(K)) ) goto error;

		/* Buildin the T_points matrix into a python list of lists */
		if (matrixToList(output, &T_points)) goto error;
		matrix_free_safe(output);
		return T_points;
	}

	Py_RETURN_NONE;

error:
	matrix_free_safe(output);
	Py_XDECREF(T_points);
	assert_other(false);
	Py_RETURN_NONE;
}


static PyObject* kmeans_fit(PyObject *self, PyObject *args) {
    PyObject *dpoint_assignments;

	/* parsing the given lists as arrays (If an error has been captured
	 * a PyExc has been set, and we return NULL */
    assert_other(0 == py_kmeans_parse_args(args));

	/* building the returned centroids' list */
	spkmeans_pass_kmeans_info_and_run(initial_centroids_indices);
    assert_other(0 == extract_datapoint_indices(&dpoint_assignments));

    free_datapoints(); // they weren't freed before because we needed to access their data
	return dpoint_assignments;
}
/**************************************************************************/






/***************************** Generic C API Functions ***************************/

/* This parses the given Python arguments into C-represented Objects + manage Reference counts of Py args 
 * Returns 0 on success, and 1 on failure */
static int py_kmeans_parse_args(PyObject *args) {
	size_t i;
	PyObject *datapoints_py = NULL;
	PyObject *initial_centroids_indices_py = NULL;
	int signal;
	
	/* Fetching Arguments from Python */
	if(!PyArg_ParseTuple(args, "OllOl", &datapoints_py, &num_data, &dim, &initial_centroids_indices_py, &K)) {
		signal = PY_ERROR;
	}

	/* Parsing the datapoints: creating the datapoints array */
	datapoints = calloc(num_data, sizeof(dpoint_t));
	if (NULL == datapoints) {
		signal = BAD_ALLOC;
		goto error;
	}

	/* Building the datapoints array: inserting inner lists */
	for (i = 0; i < num_data; i++) {
		PyObject *tmpItem;

		/* Extracting inner object of the given list */
		if ( NULL == (tmpItem = PyList_GetItem(datapoints_py, i)) ) {
			signal = PY_ERROR;
			goto error;
		}

		/* Parsing the list into an array of doubles and storing it as a datapoints in the datapoints array */
		if ( (signal = listToArray_D(tmpItem, dim, &datapoints[i].data)) ) goto error;
	}

	/* Parsing the initial centroids indices array: extracting the list from python into the array */
	if ( (signal = listToArray_L(initial_centroids_indices_py, K, &initial_centroids_indices)) ) goto error;

	return 0;


error:
	/* if any of the CPython functions fail:
	 * An error of py_parse_args doesn't trigger the free_program, 
	 * since the program hasn't advanced enough, therefore we free datapoints + initial_centroids_indices here */
	if (datapoints != NULL) {
		for (i = 0; i < num_data; i++) {
			free_datapoint(datapoints[i]);
		}
		free(datapoints_arg);
	}
	
	Py_XDECREF(datapoints_py);
	Py_XDECREF(initial_centroids_indices_py);
	return signal;

}


/* This builds a PyList out of an existing matrix.
 * Creates an untracked reference. */
static int matrixToList(const matrix_t mat, PyObject **output) {
	PyObject *pyfloat = NULL;
	size_t i, j;

	/* Creating outer list */
	if ( NULL == (*output = PyList_New(mat.rows)) ) goto error;

	for (i = 0; i < mat.rows; ++i) {
		PyObject* tmpList;

		/* Creating inner list */
		if ( NULL == (tmpList = PyList_New(mat.cols)) ) goto error;

		/* Building inner list */
		for (j = 0; j < mat.cols; ++j) {

			if(NULL == (pyfloat = PyFloat_FromDouble( matrix_get(mat, i, j) ) ) ) goto error;
			if(PyList_SetItem(tmpList, (Py_ssize_t)j, pyfloat)) goto error;
		}

		/* Inserting inner list */
		if(PyList_SetItem(*output, (Py_ssize_t)i, tmpList)) goto error;
		
	}

	return 0;

error:
	Py_XDECREF(*output);
	return PY_ERROR;
}



/* Creates a PyList out of the indices of centroids of the datapoints. */
static int extract_datapoint_indices(PyObject **output) {
	PyObject *pylong = NULL;
	size_t i;

	/* Creating outer list */
	if ( NULL == (*output = PyList_New((Py_ssize_t)num_data)) ) goto error;

	for (i = 0; i < num_data; ++i) {
		if(NULL == (pylong = PyLong_FromSize_t( datapoints[i].current_set ) ) ) goto error;
		if(PyList_SetItem(*output, (Py_ssize_t)i, pylong)) goto error;	
	}

	return 0;

error:
	Py_XDECREF(*output);
	return PY_ERROR;
}




/* This parses a python Floats' List into a C Double's array
 * No need to worry about reference counts, it's managed by py_parse_args(). */
static int listToArray_D(PyObject *list, size_t length, double** output) {
	size_t i;
	PyObject *pypoint = NULL;
	int signal;
	
	/* first check if the given PyObject is indeed a list */
	if (!PyList_Check(list)) {
		signal = PY_ERROR;
		goto error;
	}
	
	/* Initialize the array */
	(*output) = calloc(length, sizeof(**output));
	if (NULL == (*output)) {
		signal = BAD_ALLOC;
		goto error;
	}


	/* Insert the data into the array */
	for(i = 0; i < length; ++i) {
		/* PyList_GetItem returns a borrowed reference - no need to Py_DECREF */
		pypoint = PyList_GetItem(list, (Py_ssize_t)i);
		if (!PyFloat_Check(pypoint)) {
			signal = PY_ERROR;
			goto error;
		}
		(*output)[i] = (double) PyFloat_AsDouble(pypoint);
	}
	
	return 0;

error:
	if ( NULL != (*output) ) free(*output);
	return signal;
}




/* This parses a python Integers' List into a C Long's array
 * No need to worry about reference counts, it's managed by py_parse_args(). */
static int listToArray_L(PyObject *list, size_t length, size_t** output) {
	size_t i;
	PyObject *pypoint = NULL;
	int signal;
	
	/* first check if the given PyObject is indeed a list */
	if (!PyList_Check(list)) {
		signal = PY_ERROR;
		goto error;
	}

	/* Initialize the array */
	(*output) = calloc(length, sizeof(**output));
	if (NULL == (*output)) {
		signal = BAD_ALLOC;
		goto error;
	}

	/* Insert the data into the array */
	for(i = 0; i < length; ++i) {
		/* PyList_GetItem returns a borrowed reference - no need to Py_DECREF */
		pypoint = PyList_GetItem(list, (Py_ssize_t)i);
		if (!PyLong_Check(pypoint)) {
			signal = PY_ERROR;
			goto error;
		}
		(*output)[i] = (size_t) PyLong_AsLong(pypoint);
	}
	
	return 0;


error:
	if ( NULL != (*output) ) free(*output);
	return signal;
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
PyInit_spkmeans(void) {
	PyObject *m;
	m = PyModule_Create(&moduledef);
	if (!m) {
		return NULL;
	}
	return m;
}
/**************************************************************************/
