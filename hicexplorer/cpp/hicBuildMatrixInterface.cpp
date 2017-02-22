#include <Python.h>
#include "hicBuildMatrix.h"

static PyObject* createObject(PyObject* self, PyObject* args) {
    char* samFileOne;
    char* samFileTwo;
    char* outputBam;
    size_t binSize;
    char* restrictionCutFile;
    size_t minDistance;
    size_t maxDistance;
    char* restrictionSequence;
    char* outFileName;
    char* region;
    int removeSelfLigation;
    int removeSelfCircles;
    size_t minMappingQuality;
    int skipDuplicationCheck;

    if (!PyArg_ParseTuple(args, "sssskskksssiiki", 
                            &samFileOne, &samFileTwo, &outputBam, 
                            &binSize, &restrictionCutFile, &minDistance,
                            &maxDistance, &restrictionSequence, &outFileName,
                            &region, &removeSelfLigation, &removeSelfCircles,
                            &minMappingQuality, &skipDuplicationCheck))
        return NULL;
    
    bool removeSelfLigationBool = (removeSelfLigation == 1) ? true : false;
    bool removeSelfCirclesBool = (removeSelfCircles == 1) ? true : false;
    bool skipDuplicationCheckBool = (skipDuplicationCheck == 1) ? true : false;
    HicBuildMatrix* hicBuildMatrix = HicBuildMatrix(samFileOne, samFileTwo, outputBam, 
                                                    binSize, restrictionCutFile, minDistance,
                                                    maxDistance, restrictionSequence, outFileName,
                                                    region, removeSelfLigationBool, removeSelfCirclesBool,
                                                    minMappingQuality, skipDuplicationCheckBool);
    size_t addressHicBuildObject = reinterpret_cast<size_t>(hicBuildMatrix);
    PyObject * pointerToPython = Py_BuildValue("k", addressHicBuildObject);
    
    return pointerToPython;
}

static PyObject* buildMatrix(PyObject* self, PyObject* args) {

}

// definition of avaible functions for python and which function parsing fucntion in c++ should be called.
static PyMethodDef hicBuildMatrixFunctions[] = {
    {"create_object", createObject, METH_VARARGS, "Create a hicBuildMatrix cpp object."},
    {"buildMatrix", buildMatrix, METH_VARARGS, "Build the hic matrix."},
    
    {NULL, NULL, 0, NULL}
};

// definition of the module for python
PyMODINIT_FUNC
init_hicBuildMatrix(void)
{
    (void) Py_InitModule("_hicBuildMatrix", hicBuildMatrixFunctions);
}