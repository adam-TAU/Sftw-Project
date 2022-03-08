#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"



/**************************************************************************/

static PyObject* goal(PyObject *self, PyObject *args);
static PyObject* kmeans(PyObject *self, PyObject *args);


/**************************************************************************/



static PyObject* goal(PyObject *self, PyObject *args) {
	Py_RETURN_NONE;
}




static PyObject* kmeans(PyObject *self, PyObject *args) {
	Py_RETURN_NONE;
}


/**************************************************************************/


static PyMethodDef capiMethods[] = {
        {"goal",
                (PyCFunction) goal,
                METH_VARARGS,
                PyDoc_STR("Perform the wanted operations on the given datapoints, corresponding to the determined 'goal'")
        },
        {"kmeans",
                (PyCFunction) kmeans,
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
