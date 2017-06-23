
#include "powerLawNoiseReduction.h"
#include "typeDefinitions.h"


static PyObject* powerLawNoiseReduction(PyObject* self, PyObject* args) {

    // instances, features, data : list
    // window size, threshold variance, threshold_abs_mean

    uint32_t windowSize;
    float thresholdVariance;
    float thresholdAbsMean;
    float power;
    PyObject* instancesListObj, *featuresListObj, *dataListObj;
    uint32_t instancesListLength;
    uint32_t matrixSize;
    uint32_t threads;
    uint32_t iterations;
    
    if (!PyArg_ParseTuple(args, "O!O!O!IffIIfII", 
                            &PyList_Type, &instancesListObj, 
                            &PyList_Type, &featuresListObj,
                            &PyList_Type, &dataListObj,
                            &windowSize, &thresholdVariance,
                            &thresholdAbsMean,
                            &instancesListLength,
                            &matrixSize,
                            &power, &threads, &iterations))
        return NULL;
    
    PowerLawNoiseReduction* powerLawNoiseReduction = new PowerLawNoiseReduction(instancesListLength, matrixSize,
                                                            windowSize, thresholdVariance, thresholdAbsMean, threads);
    powerLawNoiseReduction->parsePythonToCpp(instancesListObj, featuresListObj, dataListObj);
    // powerLawNoiseReduction->computeGenomicMean();
    powerLawNoiseReduction->correctInteractions(power, iterations);
    PyObject* returnList = powerLawNoiseReduction->parseCppToPython();
    delete powerLawNoiseReduction;
    return returnList;
}

static PyMethodDef noiseReduction[] = {
    {"c_powerLawNoiseReduction", powerLawNoiseReduction, METH_VARARGS, "Calculate a hash value."},
    {NULL, NULL, 0, NULL}
};

// definition of the module for python
PyMODINIT_FUNC
init_c_noise_reduction(void)
{
    (void) Py_InitModule("_c_noise_reduction", noiseReduction);
}

