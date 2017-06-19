#ifndef POWER_LAW_NOISE_REDUCTION_H
#define POWER_LAW_NOISE_REDUCTION_H

#include "typeDefinitions.h"


class PowerLawNoiseReduction {
    private:
        std::unordered_map<uint32_t, std::vector<matrixElement>* >* mGenomicDistance;
        std::vector<double>* mGenomicDistanceMean;
        uint32_t mElementCount;
        uint32_t mMatrixSize;
        uint32_t mWindowSize;
        float mThresholdVariance;
        float mThresholdAbsMean;
        uint32_t mNumberOfCores;
    
    public:
        PowerLawNoiseReduction();
        PowerLawNoiseReduction(uint32_t pElementCount, uint32_t pMatrixSize, uint32_t pWindowSize,
                                float pThresholdVariance, float pThresholdAbsMean, uint32_t pNumberOfCores);
        ~PowerLawNoiseReduction();

        void parsePythonToCpp(PyObject * pInstancesListObj, PyObject * pFeaturesListObj, 
                                PyObject * pDataListObj);
        PyObject* parseCppToPython();
        void computeGenomicMean();
        void correctInteractions(float pPower, uint32_t pIterations);
};
#endif // TYPE_DEFINTIONS_BASIC_H