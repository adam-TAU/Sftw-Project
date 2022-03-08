#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"



/**************************************************************************/

static PyObject* goal(PyObject *self, PyObject *args);
static PyObject* kmeans(PyObject *self, PyObject *args);


/**************************************************************************/



static PyObject* goal(PyObject *self, PyObject *args) {
	matrix_t output; 
	spkmeans_pass_goal_info_and_run(infile, &output);
	
	if (output.data == NULL) {
		Py_RETURN_NONE;
	} else {
		/* Convert the output matrix into a list of lists */
	}
	
}




static PyObject* kmeans_fit(PyObject *self, PyObject *args) {
	/* Do the same thing as the HW2 does (added more abstract layers this time, so a different implementation might be due (differs in a maximum of one line of code) */
	Py_RETURN_NONE;
}


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
