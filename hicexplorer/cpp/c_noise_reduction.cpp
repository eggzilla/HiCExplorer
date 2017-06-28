
#include "powerLawNoiseReduction.h"
#include "typeDefinitions.h"


static PyObject* powerLawNoiseReduction(PyObject* self, PyObject* args) {

    // instances, features, data : list
    // window size, threshold variance, threshold_abs_mean

    uint32_t windowSize;
    float thresholdVariance;
    float power;
    PyObject* instancesListObj, *featuresListObj, *dataListObj;
    uint32_t instancesListLength;
    uint32_t matrixSize;
    uint32_t threads;
    // uint32_t iterations;
    uint32_t removeLowInteractionCount;
    
    if (!PyArg_ParseTuple(args, "O!O!O!IfIIfII", 
                            &PyList_Type, &instancesListObj, 
                            &PyList_Type, &featuresListObj,
                            &PyList_Type, &dataListObj,
                            &windowSize, &thresholdVariance,
                            &instancesListLength,
                            &matrixSize,
                            &power, &threads,
                            &removeLowInteractionCount))
        return NULL;
    
    PowerLawNoiseReduction* powerLawNoiseReduction = new PowerLawNoiseReduction(instancesListLength, matrixSize,
                                                            windowSize, thresholdVariance, threads,
                                                            removeLowInteractionCount);
    powerLawNoiseReduction->parsePythonToCpp(instancesListObj, featuresListObj, dataListObj);
    // powerLawNoiseReduction->computeGenomicMean();
    powerLawNoiseReduction->correctInteractions(power);
    PyObject* returnList = powerLawNoiseReduction->parseCppToPython();
    delete powerLawNoiseReduction;
    return returnList;
}

static PyObject* powerLawNoiseReduction_h5(PyObject* self, PyObject* args) {

    // instances, features, data : list
    // window size, threshold variance, threshold_abs_mean

    char* matrixPath;
    uint32_t windowSize;
    float thresholdVariance;
    float power;
    uint32_t threads;
    uint32_t removeLowInteractionCount;
    char* matrixPathOutput;
    
    if (!PyArg_ParseTuple(args, "sIffIIs", 
                            &matrixPath,
                            &windowSize, &thresholdVariance,
                            &power, &threads,
                            &removeLowInteractionCount,
                            &matrixPathOutput))
        return NULL;
    
    PowerLawNoiseReduction* powerLawNoiseReduction = new PowerLawNoiseReduction(matrixPath,
                                                            windowSize, thresholdVariance, threads,
                                                            removeLowInteractionCount,
                                                            matrixPathOutput);
    // powerLawNoiseReduction->parsePythonToCpp(instancesListObj, featuresListObj, dataListObj);
    // powerLawNoiseReduction->computeGenomicMean();
    powerLawNoiseReduction->correctInteractions(power);
    powerLawNoiseReduction->writeH5();
    // PyObject* returnList = powerLawNoiseReduction->parseCppToPython();
    delete powerLawNoiseReduction;
    // return returnList;
    Py_RETURN_NONE;
}

static PyMethodDef noiseReduction[] = {
    {"c_powerLawNoiseReduction", powerLawNoiseReduction, METH_VARARGS, "Calculate a hash value."},
    {"c_powerLawNoiseReduction_h5", powerLawNoiseReduction_h5, METH_VARARGS, "Calculate a hash value."},
    
    {NULL, NULL, 0, NULL}
};

// definition of the module for python
PyMODINIT_FUNC
init_c_noise_reduction(void)
{
    (void) Py_InitModule("_c_noise_reduction", noiseReduction);
}

