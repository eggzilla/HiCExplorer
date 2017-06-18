#include <Python.h>
#include <iostream>
#include <vector>
#include "powerLawNoiseReduction.h"



static PyObject* powerLawNoiseReduction(PyObject* self, PyObject* args) {

    // instances, features, data : list
    // window size, threshold variance, threshold_abs_mean

    uint32_t windowSize;
    float thresholdVariance;
    float thresholdAbsMean;
    PyObject* instancesListObj, *featuresListObj, *dataListObj;
    uint32_t instancesListLength;
    uint32_t matrixSize;
    
    if (!PyArg_ParseTuple(args, "O!O!O!IffII", 
                            &PyList_Type, &instancesListObj, 
                            &PyList_Type, &featuresListObj,
                            &PyList_Type, &dataListObj,
                            &windowSize, &thresholdVariance,
                            &thresholdAbsMean,
                            &instancesListLength,
                            &matrixSize))
        return NULL;
    
    PowerLawNoiseReduction* powerLawNoiseReduction = new PowerLawNoiseReduction(instancesListLength, instancesListLength, matrixSize);
    powerLawNoiseReduction->parsePythonToCpp(instancesListObj, featuresListObj, dataListObj, instancesListLength);
    return Py_BuildValue("i", static_cast<int>(hashValue));
}

static PyMethodDef noiseReduction[] = {
    {"c_powerLawNoiseReduction", powerLawNoiseReduction, METH_VARARGS, "Calculate a hash value."},
    {NULL, NULL, 0, NULL}
};

// definition of the module for python
PyMODINIT_FUNC
init_c_hash_integer(void)
{
    (void) Py_InitModule("_c_noise_reduction", noiseReduction);
}

