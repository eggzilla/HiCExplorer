#ifndef POWER_LAW_NOISE_REDUCTION_H
#define POWER_LAW_NOISE_REDUCTION_H

#include "typeDefinitions.h"
#include "h5Interface.h"

class PowerLawNoiseReduction {
    private:
        std::unordered_map<uint32_t, std::vector<matrixElement>* >* mGenomicDistance;
        std::vector<double>* mGenomicDistanceMean;
        std::vector<double>* mGenomicDistanceTotalInteractionCount;
        
        uint32_t mElementCount;
        uint32_t mMatrixSize;
        uint32_t mWindowSize;
        float mThresholdVariance;
        uint32_t mNumberOfCores;
        uint32_t mMaxElement;
        uint32_t mMinElement;
        uint32_t mRemoveLowInteractionCount;

        // char* mMatrixPath;

        H5Interface* h5Interface = NULL;
    
    public:
        PowerLawNoiseReduction();
        PowerLawNoiseReduction(uint32_t pElementCount, uint32_t pMatrixSize, uint32_t pWindowSize,
                                float pThresholdVariance, uint32_t pNumberOfCores,
                                uint32_t pRemoveLowInteractionCount);
        PowerLawNoiseReduction(char*pMatrixPath, uint32_t pWindowSize,
                                float pThresholdVariance, uint32_t pNumberOfCores,
                                uint32_t pRemoveLowInteractionCount, char* pMatrixPathOutput);
        ~PowerLawNoiseReduction();

        void parsePythonToCpp(PyObject * pInstancesListObj, PyObject * pFeaturesListObj, 
                                PyObject * pDataListObj);
        PyObject* parseCppToPython();
        void writeH5();
        void computeGenomicMean();
        void correctInteractions(float pPower);
};
#endif // TYPE_DEFINTIONS_BASIC_H