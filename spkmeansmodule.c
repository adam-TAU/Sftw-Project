#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <spkmeans.h>











/**************************************************************************/


static PyMethodDef capiMethods[] = {
        {"fit",
                (PyCFunction) fit_capi,
                METH_VARARGS,
                PyDoc_STR("TBD")
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
