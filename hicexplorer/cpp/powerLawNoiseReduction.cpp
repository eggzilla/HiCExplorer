// std::unordered_map<uint32_t, matrixElement> genomicDistance;
//         std::unordered_map<uint32_t, double> genomicDistanceMean;
    
//     public:
//         PowerLawNoiseReduction();
//         ~PowerLawNoiseReduction();

//         void parsePythonToCpp();
//         PyObject* parseCppToPython();
//         void computeGenomicMean();
//         void correctDistance(uint32_t pDistance);

#import "powerLawNoiseReduction.h"

PowerLawNoiseReduction::PowerLawNoiseReduction(uint32_t pElementCount, size_t pMatrixSize, 
                                                uint32_t pWindowSize, float pThresholdVariance, 
                                                float pThresholdAbsMean, uint32_t pNumberOfCores) {
    mGenomicDistance = new std::unordered_map<uint32_t, std::vector<matrixElement>*>(pMatrixSize);
    for (uint32_t i = 0; i < pMatrixSize; ++i) {
        mGenomicDistance->operator[](i) = new std::vector<matrixElement>();
    }
    mGenomicDistanceMean = new std::vector<double>(pMatrixSize, 0.0);
    mElementCount = pElementCount;
    mMatrixSize = pMatrixSize;
    mWindowSize = pWindowSize;
    mThresholdVariance = pThresholdVariance;
    mThresholdAbsMean = pThresholdAbsMean;
    mNumberOfCores = pNumberOfCores;
}

~PowerLawNoiseReduction::PowerLawNoiseReduction() {
    for (uint32_t i = 0; i < pMatrixSize; ++i) {
        delete mGenomicDistance->operator[](i);
    }
    delete mGenomicDistance;
    delete mGenomicDistanceMean;
}

void PowerLawNoiseReduction::parsePythonToCpp(PyObject * pInstancesListObj, PyObject * pFeaturesListObj, PyObject * pDataListObj) {
    PyObject * pyObject_instance;
    PyObject * pyObject_feature;
    PyObject * pyObject_data;
    uint32_t instance;
    uint32_t feature;
    uint32_t data;
    uint32_t genomicDistance = 0;

    for(uint32_t i = 0; i < mElementCount; ++i) {
        pyObject_instance = PyList_GetItem(pInstancesListObj, i);
        pyObject_feature = PyList_GetItem(pFeaturesListObj, i);
        pyObject_data = PyList_GetItem(pDataListObj, i); 
        
        PyArg_Parse(instanceSize_tObj, "I", &instance);
        PyArg_Parse(featureSize_tObj, "I", &feature);
        PyArg_Parse(dataSize_tObj, "I", &data);

        matrixElement element;
        element.x = instance;
        element.y = feature;
        element.data = data;

        genomicDistance = abs(instance - feature);
        (mGenomicDistance->operator[](genomicDistance))->push_back(element);
    }
}

PyObject* PowerLawNoiseReduction::parseCppToPython() {

}

void PowerLawNoiseReduction::computeGenomicMean() {

    #pragma omp parallel for num_threads(mNumberOfCores)
    for (auto it = mGenomicDistance->begin(); it != mGenomicDistance->end(); ++it) {
        double interactionCount = 0;
        for (auto itVector = (it->second)->begin(); itVector != (it->second)->end(); ++itVector)  {
            interactionCount += (*itVector).data;
        }
        mGenomicDistanceMean->operator[]((it->first)) = interactionCount;
    }
}

PowerLawNoiseReduction::correctIterations(uint32_t pDistance, float pPower) {

    #pragma omp parallel for num_threads(mNumberOfCores)
    for (uint32_t i = 0; i < mGenomicDistanceMean->size() - mWindowSize; i += mWindowSize) {
        
        std::vector<double> normalizedInteractionCount(mWindowSize);
        // normalize values
        // x_i - min(x) / max(x) - min(x)
        double maxElement = *(std::max_element(mGenomicDistanceMean->begin() + i, mGenomicDistanceMean->begin() + i + mWindowSize));
        double minElement = *(std::min_element(mGenomicDistanceMean->begin() + i, mGenomicDistanceMean->begin() + i + mWindowSize));
        double meanAbsolute = 0;
        for (uint32_t j = 0; j < mWindowSize; ++j) {
            normalizedInteractionCount[j] = (mGenomicDistanceMean->operator[](i+j) - minElement) / (maxElement - minElement);
            meanAbsolute += mGenomicDistanceMean->operator[](i+j);
        }
        meanAbsolute /= mWindowSize;

        // compute variance on normalized data
        meanNormalized = 0;
        for (uint32_t j = 0; j < mWindowSize; ++j) {
            meanNormalized += normalizedInteractionCount[j];
        }
        meanNormalized /= mWindowSize;
        double varianceNormalized = 0;
        for (uint32_t j = 0; j < mWindowSize; ++j) {
            varianceNormalized += pow(normalizedInteractionCount[j] - meanNormalized, 2);
        } 
        varianceNormalized /= mWindowSize;
        
        // skip if variance is too large and mean of the absolute interactions are not too low
        if (varianceNormalized > mThresholdVariance && meanAbsolute > mThresholdAbsMean) {
            continue;
        }
        maxElement = 0;
        minElement = std::numeric_limits<double>::max();
        // correct values based on power law
        for (uint32_t j = 0; j < mWindowSize; ++j) {
            double meanInteractionCountDifference = abs(meanAbsolute - mGenomicDistanceMean->operator[](i +j));
            uint32_t numberOfInteractionAreas = (mGenomicDistance->operator[](i+j))->size();
            std::vector<double> changeValues (numberOfInteractionAreas, 0.0);
            if (meanAbsolute < mGenomicDistanceMean->operator[](i+j)) {
                // the rich-getting-richer
                for (uint32_t k = 0; k < changeValues.size(); ++k) {
                    changeValues[k] = ((1 - (mGenomicDistance->operator[](i+j)->operator[](k)) / mGenomicDistanceMean->operator[](i+j)));
                    changeValues[k] = pow(changeValues[k], power);
                    changeValues[k] = changeValues[k] * meanInteractionCountDifference;
                    
                    if (changeValues[k] > maxElement) {
                        maxElement = changeValues[k];
                    } else if (changeValues[k] < minElement) {
                        minElement = changeValues[k];
                    }
                }
                
            } else {
                // the poor-getting-poorer
                for (uint32_t k = 0; k < changeValues.size(); ++k) {
                    changeValues[k] = mGenomicDistance->operator[](i+j)->operator[](k) / mGenomicDistanceMean->operator[](i+j);
                    changeValues[k] = pow(changeValues[k], power);
                    changeValues[k] = changeValues[k] * meanInteractionCountDifference;
                    
                    if (changeValues[k] > maxElement) {
                        maxElement = changeValues[k];
                    } else if (changeValues[k] < minElement) {
                        minElement = changeValues[k];
                    }
                }
                
            }

            // normalize values to range [0, meanInteractionCountDifference]
            for (uint32_t k = 0; k < mWindowSize; ++k) {
                changeValues[k] = meanInteractionCountDifference * ((changeValues[k] - minElement) / (maxElement - minElement));
            }
            if (meanAbsolute < mGenomicDistanceMean->operator[](i+j)) {
                for(uint32_t k = 0; k < (mGenomicDistance->operator[](i+j))->size(); ; ++k) {
                    (mGenomicDistance->operator[](i+j))->operator[](k).data += changeValues[k];
                }
            } else {
                for(uint32_t k = 0; k < (mGenomicDistance->operator[](i+j))->size(); ; ++k) {
                    (mGenomicDistance->operator[](i+j))->operator[](k).data -= changeValues[k];
                    if (mGenomicDistance->operator[](i+j))->operator[](k).data < 0) {
                        mGenomicDistance->operator[](i+j))->operator[](k).data = 0;
                    }
                }
            }
        }
    }
}