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
inline double fastPow(double a, double b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}
PowerLawNoiseReduction::PowerLawNoiseReduction(){};
PowerLawNoiseReduction::PowerLawNoiseReduction(uint32_t pElementCount, uint32_t pMatrixSize, 
                                                uint32_t pWindowSize, float pThresholdVariance, 
                                                uint32_t pNumberOfCores, uint32_t pRemoveLowInteractionCount) {
    mGenomicDistance = new std::unordered_map<uint32_t, std::vector<matrixElement>*>();
    for (uint32_t i = 0; i < pMatrixSize; ++i) {
        mGenomicDistance->operator[](i) = new std::vector<matrixElement>();
    }
    mGenomicDistanceMean = new std::vector<double>(pMatrixSize, 0.0);
    mGenomicDistanceTotalInteractionCount = new std::vector<double>(pMatrixSize, 0.0);
    mElementCount = pElementCount;
    mMatrixSize = pMatrixSize;
    mWindowSize = pWindowSize;
    mThresholdVariance = pThresholdVariance;
    mNumberOfCores = pNumberOfCores;
    mMaxElement = 0;
    mMinElement = std::numeric_limits<uint32_t>::max();
    mRemoveLowInteractionCount = pRemoveLowInteractionCount;
    // readWriteH5 = NULL;
}

PowerLawNoiseReduction::PowerLawNoiseReduction(char* pMatrixPath, uint32_t pWindowSize,
                                float pThresholdVariance, uint32_t pNumberOfCores,
                                uint32_t pRemoveLowInteractionCount) {

    // constructor to be used in case that the h5 file is read in c++ and not in python
    // be careful, a few variables are not initalized because of missing informaiton stored in h5!
    // mElementCount = pElementCount;

    h5Interface = new H5Interface(pMatrixPath);
    mGenomicDistance = h5Interface->readMatrix(pRemoveLowInteractionCount);
    mMatrixSize = mGenomicDistance->size();
    mGenomicDistanceMean = new std::vector<double>(mMatrixSize, 0.0);
    mGenomicDistanceTotalInteractionCount = new std::vector<double>(mMatrixSize, 0.0);

    mWindowSize = pWindowSize;
    mThresholdVariance = pThresholdVariance;
    mNumberOfCores = pNumberOfCores;
    mMaxElement = 0;
    mMinElement = std::numeric_limits<uint32_t>::max();
    mRemoveLowInteractionCount = pRemoveLowInteractionCount;
}

PowerLawNoiseReduction::~PowerLawNoiseReduction() {
        std::cout << __LINE__ << std::endl;

    #pragma omp parallel for num_threads(mNumberOfCores)
    for (uint32_t i = 0; i < mMatrixSize; ++i) {
        delete mGenomicDistance->operator[](i);
    }
        std::cout << __LINE__ << std::endl;
    
    delete mGenomicDistance;
        std::cout << __LINE__ << std::endl;
    
    delete mGenomicDistanceMean;
        std::cout << __LINE__ << std::endl;
    
    delete mGenomicDistanceTotalInteractionCount;
        std::cout << __LINE__ << std::endl;
    
    if (h5Interface != NULL) {
        std::cout << __LINE__ << std::endl;
        
        delete h5Interface;
    }
        std::cout << __LINE__ << std::endl;

}

void PowerLawNoiseReduction::parsePythonToCpp(PyObject * pInstancesListObj, PyObject * pFeaturesListObj, PyObject * pDataListObj) {
    
    std::cout << "Parse to c++ datastructure..." << std::endl;
    
    uint32_t instance;
    uint32_t feature;
    uint32_t data;
    uint32_t genomicDistance = 0;
   
    for(uint32_t i = 0; i < mElementCount; ++i) {

        PyObject * pyObject_instance = PyList_GetItem(pInstancesListObj, i);
        PyObject * pyObject_feature = PyList_GetItem(pFeaturesListObj, i);
        PyObject * pyObject_data = PyList_GetItem(pDataListObj, i); 
        
        PyArg_Parse(pyObject_instance, "I", &instance);
        PyArg_Parse(pyObject_feature, "I", &feature);
        PyArg_Parse(pyObject_data, "I", &data);

        if (data <= mRemoveLowInteractionCount) {
            continue;
        }
        if (data > mMaxElement) {
            mMaxElement = data;
        }
        matrixElement element;
        element.x = instance;
        element.y = feature;
        element.data = data;

        genomicDistance = abs(instance - feature);
        (mGenomicDistance->operator[](genomicDistance))->push_back(element);
    }
    std::cout << "Parse to c++ datastructure...Done!" << std::endl;
    
}

PyObject* PowerLawNoiseReduction::parseCppToPython() {

    PyObject* instances = PyList_New(0);
    PyObject* features = PyList_New(0);
    PyObject* data = PyList_New(0);
    uint32_t counter = 0;
    for (auto it = mGenomicDistance->begin(); it != mGenomicDistance->end(); ++it) {
        for (auto itVector = (it->second)->begin(); itVector != (it->second)->end(); ++itVector) {
            if ((*itVector).data <= mRemoveLowInteractionCount) {
                continue;
            }
            PyObject* pyObject_instance = Py_BuildValue("i", static_cast<int>((*itVector).x));
            PyObject* pyObject_feature = Py_BuildValue("i", static_cast<int>((*itVector).y));
            PyObject* pyObject_data = Py_BuildValue("i", static_cast<int>((*itVector).data));
            PyList_Append(instances, pyObject_instance);
            PyList_Append(features, pyObject_feature);
            PyList_Append(data, pyObject_data);
            // ++counter;            
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
        mGenomicDistanceMean->operator[]((it->first)) = interactionCount / (it->second)->size();
        mGenomicDistanceTotalInteractionCount->operator[]((it->first)) = interactionCount;
    }    
    std::cout << "Genomic mean is computed..." << std::endl;
}

void PowerLawNoiseReduction::correctInteractions(float pPower) {
    std::cout << "Start power law correction..." << std::endl;

    float* maxValue = new float [mNumberOfCores];
    float* minValue = new float [mNumberOfCores];
    
    int epsilon = 10; 
   
    this->computeGenomicMean();
    for (uint32_t i = 0; i < mNumberOfCores; ++i) {
        maxValue[i] = 0;
    }
    for (uint32_t i = 0; i < mGenomicDistanceMean->size() - mWindowSize; ++i) {
        
        std::vector<float> normalizedInteractionCount(mWindowSize, 0.0);
        
        float maxElement = *(std::max_element(mGenomicDistanceMean->begin() + i, mGenomicDistanceMean->begin() + i + mWindowSize));
        float minElement = *(std::min_element(mGenomicDistanceMean->begin() + i, mGenomicDistanceMean->begin() + i + mWindowSize));
        float meanAbsolute = 0.0;
        for (uint32_t j = 0; j < mWindowSize; ++j) {
            normalizedInteractionCount[j] = (mGenomicDistanceMean->operator[](i+j) - minElement) / (maxElement - minElement);
            meanAbsolute += mGenomicDistanceMean->operator[](i+j);
        }
        meanAbsolute /= mWindowSize;

        // compute variance on normalized data
        float meanNormalized = 0.0;
        for (uint32_t j = 0; j < mWindowSize; ++j) {
            meanNormalized += normalizedInteractionCount[j];
        }
        meanNormalized /= mWindowSize;
        float varianceNormalized = 0.0;
        for (uint32_t j = 0; j < mWindowSize; ++j) {
            varianceNormalized += pow((float) (normalizedInteractionCount[j] - meanNormalized), 2);
        } 
        varianceNormalized /= mWindowSize;
        
        // skip if variance is too large and mean of the absolute interactions are not too low
        if (varianceNormalized > mThresholdVariance) {
            continue;
        }
        
        maxElement = 0;
        minElement = std::numeric_limits<float>::max();
        // correct values based on power law
    
        #pragma omp parallel for num_threads(mNumberOfCores)    
        for (uint32_t j = 0; j < mWindowSize; ++j) {
            float  meanInteractionCountDifference = abs(meanAbsolute - mGenomicDistanceMean->operator[](i + j));
            
            // to prevent / handle overflows
            // if (meanInteractionCountDifference < 0.0001) {
            //     continue;
            // }
            // if (meanInteractionCountDifference > mGenomicDistanceMean->operator[](i + j)) {
            //     continue;
            // }

            float numberOfInteractionAreas = (mGenomicDistance->operator[](i+j))->size();
            if (numberOfInteractionAreas == 0) {
                continue;
            }
            
            std::vector<float> changeValues (numberOfInteractionAreas, 0.0);
                            
            if (mGenomicDistanceMean->operator[](i+j)-epsilon < mGenomicDistance->operator[](i+j)->operator[](0).data && 
                            mGenomicDistance->operator[](i+j)->operator[](0).data < epsilon + mGenomicDistanceMean->operator[](i+j)){
                float proportion = 1 / numberOfInteractionAreas;
                float amount = proportion * meanInteractionCountDifference;
                if (amount < 1) {
                    amount = 1;
                }
                if (meanAbsolute < mGenomicDistanceMean->operator[](i+j)) {
                    for (uint32_t k = 0; k < numberOfInteractionAreas; ++k) {
                        mGenomicDistance->operator[](i+j)->operator[](k).data += amount;
                        if (mGenomicDistance->operator[](i+j)->operator[](k).data < 0) {
                            mGenomicDistance->operator[](i+j)->operator[](k).data = 0;
                        }
                    }
                } else {
                    for (uint32_t k = 0; k < numberOfInteractionAreas; ++k) {
                        mGenomicDistance->operator[](i+j)->operator[](k).data -= amount;
                        if (mGenomicDistance->operator[](i+j)->operator[](k).data < 0) {
                            mGenomicDistance->operator[](i+j)->operator[](k).data = 0;
                        }
                    }
                }
                continue;
            }
            // applying power law
            if (meanAbsolute < mGenomicDistanceMean->operator[](i+j)) {
                // the rich-getting-richer
                
                for (uint32_t k = 0; k < changeValues.size(); ++k) {
                    
                    changeValues[k] = (float) mGenomicDistance->operator[](i+j)->operator[](k).data / (float) mGenomicDistanceTotalInteractionCount->operator[](i+j);
                    changeValues[k] = 1 - changeValues[k];
                    changeValues[k] = pow(changeValues[k], pPower);
                    changeValues[k] = changeValues[k] * meanInteractionCountDifference;
                    // catch overflows
                    if (changeValues[k] < -0.00001) {
                        // printf("%lf", changeValues[k]);
                        changeValues[k] = 0;
                        continue;
                    }
                    if (changeValues[k] > maxElement) {
                        maxElement = changeValues[k];
                    } else if (changeValues[k] < minElement) {
                        minElement = changeValues[k];
                    }
                }
                
            } else {
                // the poor-getting-poorer
                
                for (uint32_t k = 0; k < changeValues.size(); ++k) {
                    
                    changeValues[k] = mGenomicDistance->operator[](i+j)->operator[](k).data / mGenomicDistanceTotalInteractionCount->operator[](i+j);
                    changeValues[k] = pow(changeValues[k], pPower);
                    changeValues[k] = changeValues[k] * meanInteractionCountDifference;
                    // catch / handle overflows
                    if (changeValues[k] < 0.001) {
                        changeValues[k] = 0;
                        continue;
                    }
                    if (changeValues[k] > maxElement) {
                        maxElement = changeValues[k];
                    } else if (changeValues[k] < minElement) {
                        minElement = changeValues[k];
                    }
                }
                
            }

            // normalize values to range [0, meanInteractionCountDifference]
            double differenceRange = maxElement - minElement;
            
            for (uint32_t k = 0; k < changeValues.size(); ++k) {
                
                if (changeValues[k] > meanInteractionCountDifference) {
                    changeValues[k] = 0;                        
                    continue;
                }
                
                changeValues[k] = changeValues[k] - minElement;
                changeValues[k] /= differenceRange;
                changeValues[k] = floor(changeValues[k] * meanInteractionCountDifference);
                
            }
            
            
            if (meanAbsolute < mGenomicDistanceMean->operator[](i+j)) {
                for(uint32_t k = 0; k < (mGenomicDistance->operator[](i+j))->size(); ++k) {
                    
                    (mGenomicDistance->operator[](i+j))->operator[](k).data += static_cast<int32_t>(changeValues[k]);
                    if ((mGenomicDistance->operator[](i+j))->operator[](k).data < 0) {
                        (mGenomicDistance->operator[](i+j))->operator[](k).data = 0;
                    }
                    if ((mGenomicDistance->operator[](i+j))->operator[](k).data > maxValue[omp_get_thread_num()]) {
                        maxValue[omp_get_thread_num()] = (mGenomicDistance->operator[](i+j))->operator[](k).data;
                    }
                }
            } 
            else {
                for(uint32_t k = 0; k < (mGenomicDistance->operator[](i+j))->size(); ++k) {
                    
                    (mGenomicDistance->operator[](i+j))->operator[](k).data -= static_cast<int32_t>(changeValues[k]);
                    if ((mGenomicDistance->operator[](i+j))->operator[](k).data < 0) {
                        mGenomicDistance->operator[](i+j)->operator[](k).data = 0;
                    }
                    if ((mGenomicDistance->operator[](i+j))->operator[](k).data > maxValue[omp_get_thread_num()]) {
                        maxValue[omp_get_thread_num()] = (mGenomicDistance->operator[](i+j))->operator[](k).data;
                    }
                }
            }
        }
    } 
    

    uint32_t maximalValue = 0;
    for (uint32_t i = 0; i < mNumberOfCores; ++i) {
        if (maxValue[i] > maximalValue) {
            maximalValue = maxValue[i];
        }
    }
    delete [] maxValue;

    uint32_t minimalValue = 0;
    for (uint32_t i = 0; i < mNumberOfCores; ++i) {
        if (minValue[i] > minimalValue) {
            minimalValue = minValue[i];
        }
    }
    delete [] minValue;

    if (maximalValue < mMaxElement) {
        maximalValue = mMaxElement;
    }

    int max = 0; 
    #pragma omp parallel num_threads(mNumberOfCores)
    for (uint32_t i = 0; i < mGenomicDistance->size(); ++i) {
        auto it = mGenomicDistance->begin();
        std::advance(it, i);
        for (auto itVector = (it->second)->begin(); itVector != (it->second)->end(); ++itVector)  {
            (*itVector).data = static_cast<int32_t>((float) mMaxElement * (((float)(*itVector).data) /  (float) maximalValue));
            if ((*itVector).data > max) {
                max = (*itVector).data;
            }
        }
    }  
    std::cout << "Start powerlaw correction...Done!" << std::endl;
    
}