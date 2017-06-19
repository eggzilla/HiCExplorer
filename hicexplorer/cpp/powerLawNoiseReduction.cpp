// std::unordered_map<uint32_t, matrixElement> genomicDistance;
//         std::unordered_map<uint32_t, double> genomicDistanceMean;
    
//     public:
//         PowerLawNoiseReduction();
//         ~PowerLawNoiseReduction();

//         void parsePythonToCpp();
//         PyObject* parseCppToPython();
//         void computeGenomicMean();
//         void correctDistance(uint32_t pDistance);

#include "powerLawNoiseReduction.h"
PowerLawNoiseReduction::PowerLawNoiseReduction(){};
PowerLawNoiseReduction::PowerLawNoiseReduction(uint32_t pElementCount, uint32_t pMatrixSize, 
                                                uint32_t pWindowSize, float pThresholdVariance, 
                                                float pThresholdAbsMean, uint32_t pNumberOfCores) {
    mGenomicDistance = new std::unordered_map<uint32_t, std::vector<matrixElement>*>();
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

PowerLawNoiseReduction::~PowerLawNoiseReduction() {
    for (uint32_t i = 0; i < mMatrixSize; ++i) {
        delete mGenomicDistance->operator[](i);
    }
    delete mGenomicDistance;
    delete mGenomicDistanceMean;
}

void PowerLawNoiseReduction::parsePythonToCpp(PyObject * pInstancesListObj, PyObject * pFeaturesListObj, PyObject * pDataListObj) {
    // PyObject * pyObject_instance;
    // PyObject * pyObject_feature;
    // PyObject * pyObject_data;
    uint32_t instance;
    uint32_t feature;
    uint32_t data;
    uint32_t genomicDistance = 0;

    #pragma omp parallel for num_threads(mNumberOfCores)    
    for(uint32_t i = 0; i < mElementCount; ++i) {
        PyObject * pyObject_instance = PyList_GetItem(pInstancesListObj, i);
        PyObject * pyObject_feature = PyList_GetItem(pFeaturesListObj, i);
        PyObject * pyObject_data = PyList_GetItem(pDataListObj, i); 
        
        PyArg_Parse(pyObject_instance, "I", &instance);
        PyArg_Parse(pyObject_feature, "I", &feature);
        PyArg_Parse(pyObject_data, "I", &data);

        matrixElement element;
        element.x = instance;
        element.y = feature;
        element.data = data;

        genomicDistance = abs(instance - feature);
        #pragma omp critical
        (mGenomicDistance->operator[](genomicDistance))->push_back(element);
    }
    // delete pInstancesListObj;
    // delete pFeaturesListObj;
    // delete pDataListObj;
}

PyObject* PowerLawNoiseReduction::parseCppToPython() {

    PyObject* instances = PyList_New(mElementCount);
    PyObject* features = PyList_New(mElementCount);
    PyObject* data = PyList_New(mElementCount);
    uint32_t counter = 0;
    for (auto it = mGenomicDistance->begin(); it != mGenomicDistance->end(); ++it) {
        for (auto itVector = (it->second)->begin(); itVector != (it->second)->end(); ++itVector) {
            PyObject* pyObject_instance = Py_BuildValue("i", static_cast<int>((*itVector).x));
            PyObject* pyObject_feature = Py_BuildValue("i", static_cast<int>((*itVector).y));
            PyObject* pyObject_data = Py_BuildValue("i", static_cast<int>((*itVector).data));
            PyList_SetItem(instances, counter, pyObject_instance);
            PyList_SetItem(features, counter, pyObject_feature);
            PyList_SetItem(data, counter, pyObject_data);
            ++counter;            
        }
    } 
    PyObject* returnList = PyList_New(3);
    PyList_SetItem(returnList, 0, instances);
    PyList_SetItem(returnList, 1, features);
    PyList_SetItem(returnList, 2, data);

    return returnList;
}

void PowerLawNoiseReduction::computeGenomicMean() {

    #pragma omp parallel for num_threads(mNumberOfCores)
    for (uint32_t i = 0; i < mGenomicDistance->size(); ++i) {
        auto it = mGenomicDistance->begin();
        std::advance(it, i);
        double interactionCount = 0;
        for (auto itVector = (it->second)->begin(); itVector != (it->second)->end(); ++itVector)  {
            interactionCount += (*itVector).data;
        }
        mGenomicDistanceMean->operator[]((it->first)) = interactionCount;
    }    
}

void PowerLawNoiseReduction::correctInteractions(float pPower, uint32_t pIterations) {

    for (uint32_t m = 0; m < pIterations; ++m) {
        std::cout << "Iteration: " << m << "/" << pIterations << std::endl;
        #pragma omp parallel for num_threads(mNumberOfCores)
        for (uint32_t i = 0; i < mGenomicDistanceMean->size() - mWindowSize; ++i) {
            
            std::vector<double> normalizedInteractionCount(mWindowSize, 0.0);
            // normalize values
            // x_i - min(x) / max(x) - min(x)
            double maxElement = *(std::max_element(mGenomicDistanceMean->begin() + i, mGenomicDistanceMean->begin() + i + mWindowSize));
            double minElement = *(std::min_element(mGenomicDistanceMean->begin() + i, mGenomicDistanceMean->begin() + i + mWindowSize));
            double meanAbsolute = 0.0;
            for (uint32_t j = 0; j < mWindowSize; ++j) {
                normalizedInteractionCount[j] = (mGenomicDistanceMean->operator[](i+j) - minElement) / (maxElement - minElement);
                meanAbsolute += mGenomicDistanceMean->operator[](i+j);
            }
            meanAbsolute /= mWindowSize;

            // compute variance on normalized data
            double meanNormalized = 0.0;
            for (uint32_t j = 0; j < mWindowSize; ++j) {
                meanNormalized += normalizedInteractionCount[j];
            }
            meanNormalized /= mWindowSize;
            double varianceNormalized = 0.0;
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
                
                double meanInteractionCountDifference = abs(meanAbsolute - mGenomicDistanceMean->operator[](i + j));
                
                uint32_t numberOfInteractionAreas = (mGenomicDistance->operator[](i+j))->size();
                if (numberOfInteractionAreas == 0) {
                    continue;
                }
                
                std::vector<double> changeValues (numberOfInteractionAreas, 0.0);
                
                if (meanAbsolute < mGenomicDistanceMean->operator[](i+j)) {
                    // the rich-getting-richer
                    
                    for (uint32_t k = 0; k < changeValues.size(); ++k) {
                        changeValues[k] = mGenomicDistance->operator[](i+j)->operator[](k).data / mGenomicDistanceMean->operator[](i+j);
                        changeValues[k] = 1 - changeValues[k];
                        changeValues[k] = pow(changeValues[k], pPower);
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
                        
                        changeValues[k] = mGenomicDistance->operator[](i+j)->operator[](k).data / mGenomicDistanceMean->operator[](i+j);
                        changeValues[k] = pow(changeValues[k], pPower);
                        changeValues[k] = changeValues[k] * meanInteractionCountDifference;
                        
                        if (changeValues[k] > maxElement) {
                            maxElement = changeValues[k];
                        } else if (changeValues[k] < minElement) {
                            minElement = changeValues[k];
                        }
                    }
                    
                }

                // normalize values to range [0, meanInteractionCountDifference]
                for (uint32_t k = 0; k < changeValues.size(); ++k) {
                    changeValues[k] = meanInteractionCountDifference * ((changeValues[k] - minElement) / (maxElement - minElement));
                    if (changeValues[k] != changeValues[k]) {
                        // according to IEEE f != f is true if f == nan
                        // https://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c
                        changeValues[k] = 0;
                    }
                }
                
                if (meanAbsolute < mGenomicDistanceMean->operator[](i+j)) {
                    for(uint32_t k = 0; k < (mGenomicDistance->operator[](i+j))->size(); ++k) {
                        (mGenomicDistance->operator[](i+j))->operator[](k).data += changeValues[k];
                        if ((mGenomicDistance->operator[](i+j))->operator[](k).data < 10) {
                            (mGenomicDistance->operator[](i+j))->operator[](k).data = 0;
                        }
                    }
                } 
                else {
                    for(uint32_t k = 0; k < (mGenomicDistance->operator[](i+j))->size(); ++k) {
                        (mGenomicDistance->operator[](i+j))->operator[](k).data -= changeValues[k];
                        if ((mGenomicDistance->operator[](i+j))->operator[](k).data < 10) {
                            mGenomicDistance->operator[](i+j)->operator[](k).data = 0;
                        }
                    }
                }
            }
        }
        this->computeGenomicMean();
    }
}