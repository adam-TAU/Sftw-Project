#define PY_SSIZE_T_CLEAN
#include "spkmeans.h"
#include <Python.h>

#define PY_ERROR -1

/**************************************************************************/
static PyObject *run_goal(PyObject *self, PyObject *args);
static PyObject *kmeans_fit(PyObject *self, PyObject *args);

static int matrixToList(const matrix_t mat, PyObject **output);
static int listToArray_D(PyObject *list, size_t length, double **output);
static int listToArray_L(PyObject *list, size_t length, size_t **output);
static int py_kmeans_parse_args(PyObject *);

/**************************************************************************/

/**************************************************************************/
static double **datapoints_arg = NULL;
static size_t *initial_centroids_indices = NULL;
/**************************************************************************/

/************************* configuring the C API
 * ****************************************************/
static PyObject *run_goal(PyObject *self, PyObject *args) {
    char *infile;
    matrix_t output;
    PyObject *py_output = NULL;

    output.data = NULL;

    /* Fetch the string of the infile */
    if(!PyArg_ParseTuple(args, "lss", &K, &goal, &infile))
        goto error;

    /* Perform the wanted goal's operation and return the result (if there's
     * any) */
    if(spkmeans_pass_goal_info_and_run(infile, &output))
        goto error;

    /* Build the matrix that was created by the goal */
    if(matrixToList(output, &py_output))
        goto error;

    /* Free and return */
    matrix_free(output);
    free_program();
    
    return py_output;

error:
    matrix_free_safe(output);
    free_program();
    Py_XDECREF(py_output);
    assert_other(false);
    return NULL;
}

static PyObject *kmeans_fit(PyObject *self, PyObject *args) {
    PyObject *py_output;
    matrix_t centroids_mat;
    size_t i, j;

    centroids_mat.data = NULL; // in case of an error, `centroids_mat`'s data
                               // field is freed if it's not null

    /* parsing the given lists as arrays (If an error has been captured
     * a PyExc has been set, and we return NULL */
    assert_other(0 == py_kmeans_parse_args(args));

    /* building the returned centroids' list */
    spkmeans_pass_kmeans_info_and_run(initial_centroids_indices);

    /* Creating the matrix that will hold the centroids */
    if(matrix_new(K, dim, &centroids_mat))
        goto error;

    /* Building the matrix that will hold the centroids */
    for(i = 0; i < K; i++) {
        for(j = 0; j < dim; j++) {
            matrix_set(centroids_mat, i, j, sets[i].current_centroid.data[j]);
        }
    }

    /* Build the matrix that was created by the centroids */
    if(matrixToList(centroids_mat, &py_output))
        goto error;

    /* Free and return */
    matrix_free(centroids_mat);
    free_program();
    
    return py_output;

error:
    matrix_free_safe(centroids_mat);
    free_program();
    assert_other(false);
    return NULL;
}
/**************************************************************************/

/***************************** Generic C API Functions
 * ***************************/

/* This parses the given Python arguments into C-represented Objects + manage
 * Reference counts of Py args Returns 0 on success, and 1 on failure */
static int py_kmeans_parse_args(PyObject *args) {
    size_t i;
    PyObject *datapoints_py = NULL;
    PyObject *initial_centroids_indices_py = NULL;
    int signal;

    /* Fetching Arguments from Python */
    if(!PyArg_ParseTuple(args, "OllOl", &datapoints_py, &num_data, &dim,
                         &initial_centroids_indices_py, &K))
    {
        signal = PY_ERROR;
    }

    /* Parsing the datapoints: creating the datapoints array */
    datapoints = calloc(num_data, sizeof(dpoint_t));
    if(NULL == datapoints) {
        signal = BAD_ALLOC;
        goto error;
    }

    /* Building the datapoints array: inserting inner lists */
    for(i = 0; i < num_data; i++) {
        PyObject *tmpItem;

        /* Extracting inner object of the given list */
        if(NULL == (tmpItem = PyList_GetItem(datapoints_py, i))) {
            signal = PY_ERROR;
            goto error;
        }

        /* Parsing the list into an array of doubles and storing it as a
         * datapoints in the datapoints array */
        if((signal = listToArray_D(tmpItem, dim, &datapoints[i].data)))
            goto error;
    }

    /* Parsing the initial centroids indices array: extracting the list from
     * python into the array */
    if((signal = listToArray_L(initial_centroids_indices_py, K,
                               &initial_centroids_indices)))
        goto error;

    return 0;

error:
    /* if any of the CPython functions fail:
     * An error of py_parse_args doesn't trigger the free_program,
     * since the program hasn't advanced enough, therefore we free datapoints +
     * initial_centroids_indices here */
    if(datapoints != NULL) {
        for(i = 0; i < num_data; i++) {
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
    if(NULL == (*output = PyList_New(mat.rows)))
        goto error;

    for(i = 0; i < mat.rows; ++i) {
        PyObject *tmpList;

        /* Creating inner list */
        if(NULL == (tmpList = PyList_New(mat.cols)))
            goto error;

        /* Building inner list */
        for(j = 0; j < mat.cols; ++j) {

            if(NULL == (pyfloat = PyFloat_FromDouble(matrix_get(mat, i, j))))
                goto error;
            if(PyList_SetItem(tmpList, (Py_ssize_t)j, pyfloat))
                goto error;
        }

        /* Inserting inner list */
        if(PyList_SetItem(*output, (Py_ssize_t)i, tmpList))
            goto error;
    }

    return 0;

error:
    Py_XDECREF(*output);
    return PY_ERROR;
}

/* This parses a python Floats' List into a C Double's array
 * No need to worry about reference counts, it's managed by py_parse_args(). */
static int listToArray_D(PyObject *list, size_t length, double **output) {
    size_t i;
    PyObject *pypoint = NULL;
    int signal;

    /* first check if the given PyObject is indeed a list */
    if(!PyList_Check(list)) {
        signal = PY_ERROR;
        goto error;
    }

    /* Initialize the array */
    (*output) = calloc(length, sizeof(**output));
    if(NULL == (*output)) {
        signal = BAD_ALLOC;
        goto error;
    }

    /* Insert the data into the array */
    for(i = 0; i < length; ++i) {
        /* PyList_GetItem returns a borrowed reference - no need to Py_DECREF */
        pypoint = PyList_GetItem(list, (Py_ssize_t)i);
        if(!PyFloat_Check(pypoint)) {
            signal = PY_ERROR;
            goto error;
        }
        (*output)[i] = (double)PyFloat_AsDouble(pypoint);
    }

    return 0;

error:
    if(NULL != (*output))
        free(*output);
    return signal;
}

/* This parses a python Integers' List into a C Long's array
 * No need to worry about reference counts, it's managed by py_parse_args(). */
static int listToArray_L(PyObject *list, size_t length, size_t **output) {
    size_t i;
    PyObject *pypoint = NULL;
    int signal;

    /* first check if the given PyObject is indeed a list */
    if(!PyList_Check(list)) {
        signal = PY_ERROR;
        goto error;
    }

    /* Initialize the array */
    (*output) = calloc(length, sizeof(**output));
    if(NULL == (*output)) {
        signal = BAD_ALLOC;
        goto error;
    }

    /* Insert the data into the array */
    for(i = 0; i < length; ++i) {
        /* PyList_GetItem returns a borrowed reference - no need to Py_DECREF */
        pypoint = PyList_GetItem(list, (Py_ssize_t)i);
        if(!PyLong_Check(pypoint)) {
            signal = PY_ERROR;
            goto error;
        }
        (*output)[i] = (size_t)PyLong_AsLong(pypoint);
    }

    return 0;

error:
    if(NULL != (*output))
        free(*output);
    return signal;
}
/**************************************************************************/

/**************************************************************************/
static PyMethodDef capiMethods[] = {
    {"goal", (PyCFunction)run_goal, METH_VARARGS,
     PyDoc_STR("Perform the wanted operations on the given datapoints, "
               "corresponding to the determined 'goal'")},
    {"kmeans_fit", (PyCFunction)kmeans_fit, METH_VARARGS,
     PyDoc_STR("Given a set of datapoints, an array of the indices of the "
               "initial centroids (induced from kmeans++'s first step), "
               "perform the kmeans algorithm")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {PyModuleDef_HEAD_INIT, "spkmeans", NULL,
                                       -1, capiMethods};

PyMODINIT_FUNC PyInit_spkmeans(void) {
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if(!m) {
        return NULL;
    }
    return m;
}
/**************************************************************************/
