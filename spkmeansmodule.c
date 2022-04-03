#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"



/**************************************************************************/
static PyObject* goal(PyObject *self, PyObject *args);
static PyObject* kmeans_fit(PyObject *self, PyObject *args);

/**************************************************************************/



static PyObject* goal(PyObject *self, PyObject *args) {
	char* infile;
	matrix_t output;
	
	/* Fetch the string of the infile */
    if(!PyArg_ParseTuple(args, "s", &infile)) {
		return NULL;
    }
	 
	 /* Perform the wanted goal's operation and return the result (if there's any) */
	spkmeans_pass_goal_info_and_run(infile, &output);
	
	if (output.data == NULL) {
		Py_RETURN_NONE;
	} else {
		/* Convert the output matrix into a list of lists */
		return NULL;
	}
}


/************************* configuring the C API ****************************************************/

static PyObject* kmeans_fit(PyObject *self, PyObject *args) {
    int i;
	PyObject *centroids_py;
	PyObject *result;
	
    /* parsing the given lists as arrays (If an error has been captured
     * a PyExc has been set, and we return NULL */
	if (1 == py_parse_args(args)) return NULL;

    /* building the returned centroids' list */
    centroids_c = fit_c();
    centroids_py = PyList_New(K);
    for(i = 0; i < K; ++i) {
    		PyObject *tmpList = arrayToList_D(centroids_c[i], dim);
    		if(NULL == tmpList) goto failed;
    		
    		/* Heads up! PyList_SetItem steals a reference, so Py_DECREF(centroids_py) 
    		 * would Py_DECREF(tmpList) to 0 */
            if(0 != PyList_SetItem(centroids_py, i, tmpList)) goto failed;
    }
    result = Py_BuildValue("O", centroids_py);
	Py_DECREF(result);
	
    /* free-ing the program and returning the centroids */
    free_program();
    return result;
    
    
    failed:
    /* if parsing args went corretly tho we got an error */
    Py_DECREF(centroids_py);
    free_program();
    return NULL;
}


/***************************** Generic C API Functions ***************************/

/* This parses the given Python arguments into C-represented Objects + manage Reference counts of Py args 
 * Returns 0 on success, and 1 on failure */
static int py_parse_args(PyObject *args) {
    int i;
    PyObject *datapoints_py;
    PyObject *initial_centroids_py;
    
    /* Fetching Arguments from Python */
    if(!PyArg_ParseTuple(args, "iiidiOO", &K, &dim, &num_data, &datapoints_py, &initial_centroids_py)) {
		return 1;
    }
    
    /* Parsing the fetched Arguments into C-represented Objects */
    datapoints_arg = (double**)calloc(num_data, sizeof(double*));
    assert_other(NULL != datapoints_arg);
    
    for (i = 0; i < num_data; i++) {
    		/* PyList_GetItem returns a borrowed reference - no need to Py_DECREF */
    		PyObject *tmpItem = PyList_GetItem(datapoints_py, i);
    		if (NULL == tmpItem) goto failed;
    		
    		double *tmpArray = listToArray_D(tmpItem, dim);
    		if (NULL == (datapoints_arg[i] = tmpArray)) goto failed;
    }
    
    int *tmpArray = listToArray_I(initial_centroids_py, K);
    if (NULL == (initial_centroids_indices = tmpArray)) goto failed;
    Py_DECREF(datapoints_py);
   	return 0;
   	
   	
   	failed:
   	/* if any of the CPython functions fail:
   	 * An error of py_parse_args doesn't trigger the free_program, 
   	 * since the program hasn't advanced enough, therefore we free datapoints_arg here */
   	if (datapoints_arg != NULL) {
		for (i = 0; i < num_data; i++) {
			if (datapoints_arg[i] != NULL) free(datapoints_arg[i]);
		}
		free(datapoints_arg);
	}
	Py_DECREF(datapoints_py);
	return 1;
    
}


/* This builds a PyList out of an existing array
 * No need to care about reference counts - it's managed by fit_capi */
static PyObject* arrayToList_D(const double *const array, int length) {
        PyObject *list, *pyfloat;
        int i;

        list = PyList_New(length);
        for (i = 0; i < length; ++i) {
        		if(NULL == (pyfloat = PyFloat_FromDouble(array[i]))) {
        			PyErr_SetString(PyExc_TypeError, "Double to Float conversion failed!");
        			return NULL;
        		}
                if(0 != PyList_SetItem(list, (Py_ssize_t)i, pyfloat)) {
              		PyErr_SetString(PyExc_IndexError, "Index out of bounds or the index isn't an integer!");
              		return NULL;
              	}
        }
        return list;
}




/* This parses a python Floats' List into a C Double's array
 * No need to worry about reference counts, it's managed by py_parse_args() */
static double* listToArray_D(PyObject *list, int length) {
        int i;
        double* result;
        PyObject *pypoint;

        /* first check if the given PyObject is indeed a list */
        if (!PyList_Check(list)) {
                PyErr_SetString(PyExc_TypeError, "The passed argument isn't a list");
                return NULL;
        }

        /* casting the list of float to an array of doubles */
        result = (double*)calloc(length, sizeof(double));
        assert_other(NULL != result); 
        
        for(i = 0; i < length; ++i) {
        		/* PyList_GetItem returns a borrowed reference - no need to Py_DECREF */
                pypoint = PyList_GetItem(list, (Py_ssize_t)i);
                if (!PyFloat_Check(pypoint)) {
                        PyErr_SetString(PyExc_TypeError, "Must pass an list of floats");
                        goto failed;
                }
                result[i] = (double) PyFloat_AsDouble(pypoint);
        }
        return result;
        
        
        failed:
        /* If any of the CPython functions fail */
        if (result != NULL) free(result);
        return NULL;
}




/* This parses a python Integer's List into a C Integer's array
 * No need to worry about reference counts, it's managed by py_parse_args() */
static int* listToArray_I(PyObject *list, int length) {
        int i;
        int* result;
        PyObject *pypoint;

        /* first check if the given PyObject is indeed a list */
        if (!PyList_Check(list)) {
                PyErr_SetString(PyExc_TypeError, "The passed argument isn't a list!");
                return NULL;
        }

        /* casting the list of integers to an array of integers */
        result = (int*)calloc(length, sizeof(int));
        assert_other(NULL != result);
        
        for(i = 0; i < length; ++i) {
        		/* PyList_GetItem returns a borrowed reference - no need to Py_DECREF */
                pypoint = PyList_GetItem(list, (Py_ssize_t)i);
                if (!PyLong_Check(pypoint)) {
                        PyErr_SetString(PyExc_TypeError, "Must pass an list of floats!");
                        goto failed;
                }
                result[i] = (int) PyLong_AsLong(pypoint);
        }
        return result;
        
        failed:
        /* If any of the CPython functions fail */
        if (result != NULL) free(result);
        return NULL;
}



/**************************************************************************/

/**************************************************************************/


static PyMethodDef capiMethods[] = {
        {"goal",
                (PyCFunction) goal,
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
