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
                                                float pThresholdAbsMean, uint32_t pNumberOfCores) {
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
    mThresholdAbsMean = pThresholdAbsMean;
    mNumberOfCores = pNumberOfCores;
    mMaxElement = 0;
    mMinElement = std::numeric_limits<uint32_t>::max();
}

PowerLawNoiseReduction::~PowerLawNoiseReduction() {
    for (uint32_t i = 0; i < mMatrixSize; ++i) {
        delete mGenomicDistance->operator[](i);
    }
    delete mGenomicDistance;
    delete mGenomicDistanceMean;
    delete mGenomicDistanceTotalInteractionCount;
}

void PowerLawNoiseReduction::parsePythonToCpp(PyObject * pInstancesListObj, PyObject * pFeaturesListObj, PyObject * pDataListObj) {
    
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
        mGenomicDistanceMean->operator[]((it->first)) = interactionCount / (it->second)->size();
        mGenomicDistanceTotalInteractionCount->operator[]((it->first)) = interactionCount;
    }    
}

void PowerLawNoiseReduction::correctInteractions(float pPower, uint32_t pIterations) {

    float* maxValue = new float [mNumberOfCores];
    float* minValue = new float [mNumberOfCores];
    
    int epsilon = 10;
    for (uint32_t m = 0; m < pIterations; ++m) {
        std::cout << "Iteration: " << m << "/" << pIterations << std::endl;
        this->computeGenomicMean();
        for (uint32_t i = 0; i < mNumberOfCores; ++i) {
            maxValue[i] = 0;
        }
        #pragma omp parallel for num_threads(mNumberOfCores)
        for (uint32_t i = 0; i < mGenomicDistanceMean->size() - mWindowSize; i = i + mWindowSize) {
            
            std::vector<float> normalizedInteractionCount(mWindowSize, 0.0);
            // if (i == 1) {
            //     break;
            // }
            // normalize values
            // x_i - min(x) / max(x) - min(x)
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
            for (uint32_t j = 0; j < mWindowSize; ++j) {
                // if (j == 1) break;
                float  meanInteractionCountDifference = abs(meanAbsolute - mGenomicDistanceMean->operator[](i + j));
                
                float numberOfInteractionAreas = (mGenomicDistance->operator[](i+j))->size();
                if (numberOfInteractionAreas == 0) {
                    continue;
                }
                
                std::vector<float> changeValues (numberOfInteractionAreas, 0.0);
                // ~same interaction count at all positions: +/- epsilon
                // std::cout << "mean: " << mGenomicDistanceMean->operator[](i+j) << " data: ";
                // for (uint16_t k = 0; k < numberOfInteractionAreas; ++k) {
                //     std::cout << mGenomicDistance->operator[](i+j)->operator[](k).data << ", ";
                // }
                // std::cout << std::endl;
                // std::cout << std::endl;
                
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
                // applying power laws
                if (meanAbsolute < mGenomicDistanceMean->operator[](i+j)) {
                    // the rich-getting-richer
                    
                    for (uint32_t k = 0; k < changeValues.size(); ++k) {
                        if (mGenomicDistance->operator[](i+j)->operator[](k).data > 700) {
                            std::cout << "data: " << mGenomicDistance->operator[](i+j)->operator[](k).data;
                            std::cout << " mean: " << mGenomicDistanceMean->operator[](i+j);
                            changeValues[k] = (float) mGenomicDistance->operator[](i+j)->operator[](k).data / (float) mGenomicDistanceTotalInteractionCount->operator[](i+j);
                            std::cout << "1- x: " << 1 - changeValues[k];
                            changeValues[k] = 1 - changeValues[k];
                            
                            std::cout << " pow: " << pow(changeValues[k], pPower);
                            changeValues[k] = pow(changeValues[k], pPower);
                            
                            std::cout << " amount: " << changeValues[k] * meanInteractionCountDifference;
                            changeValues[k] = changeValues[k] * meanInteractionCountDifference;
                            std::cout << std::endl;
                        } else {
                        changeValues[k] = (float) mGenomicDistance->operator[](i+j)->operator[](k).data / (float) mGenomicDistanceTotalInteractionCount->operator[](i+j);
                        changeValues[k] = 1 - changeValues[k];
                        changeValues[k] = pow(changeValues[k], pPower);
                        changeValues[k] = changeValues[k] * meanInteractionCountDifference;
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
                        
                        if (mGenomicDistance->operator[](i+j)->operator[](k).data > 700) {
                            std::cout << "data: " << mGenomicDistance->operator[](i+j)->operator[](k).data;
                            std::cout << " mean: " << mGenomicDistanceMean->operator[](i+j);
                            changeValues[k] = mGenomicDistance->operator[](i+j)->operator[](k).data / mGenomicDistanceTotalInteractionCount->operator[](i+j);
                            
                            std::cout << " pow: " << pow(changeValues[k], pPower);
                            changeValues[k] = pow(changeValues[k], pPower);
                            
                            std::cout << " amount: " << changeValues[k] * meanInteractionCountDifference;
                            changeValues[k] = changeValues[k] * meanInteractionCountDifference;
                            std::cout << std::endl;
                        } else {
                            changeValues[k] = mGenomicDistance->operator[](i+j)->operator[](k).data / mGenomicDistanceTotalInteractionCount->operator[](i+j);
                            changeValues[k] = pow(changeValues[k], pPower);
                            changeValues[k] = changeValues[k] * meanInteractionCountDifference;
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
                    // std::cout << "uncorrected: " << changeValues[k] << std::endl;
                    // std::cout << "MeanInteractionCOuntDifference: " << meanInteractionCountDifference;
                    // std::cout << " minElement: " << minElement;
                    // std::cout << " maxElement: " << maxElement << std::endl;
                    
                    changeValues[k] = changeValues[k] - minElement;
                    changeValues[k] /= differenceRange;
                    changeValues[k] = floor(changeValues[k] * meanInteractionCountDifference);
                    if (changeValues[k] < pow(2, 10) * -1 || changeValues[k] > pow(2, 10)) {
                        // std::cout << "Overflow" << std::endl;
                        printf("%.2Lf\n", changeValues[k]);
                        printf("%lf\n", pow(2, 30) * -1);
                    }
                    
                    
                    // changeValues[k] = meanInteractionCountDifference * ((changeValues[k] - minElement) / differenceRange);
                    // if (changeValues[k] != changeValues[k]) {
                    //     // according to IEEE f != f is true if f == nan
                    //     // https://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c
                    //     changeValues[k] = 0;
                    // }
                    // std::cout << "corrected: " << changeValues[k] << std::endl;
                }
                
                
                if (meanAbsolute < mGenomicDistanceMean->operator[](i+j)) {
                    for(uint32_t k = 0; k < (mGenomicDistance->operator[](i+j))->size(); ++k) {
                        // std::cout << "RIchorgData: " << (mGenomicDistance->operator[](i+j))->operator[](k).data << std::endl;
                        if ((mGenomicDistance->operator[](i+j))->operator[](k).data > 800) {
                            std::cout << "More than 800: " << (mGenomicDistance->operator[](i+j))->operator[](k).data;
                            std::cout << " Correction factor: " << static_cast<float>(changeValues[k]) << std::endl;
                        }
                        (mGenomicDistance->operator[](i+j))->operator[](k).data += static_cast<int32_t>(changeValues[k]);
                        if ((mGenomicDistance->operator[](i+j))->operator[](k).data < 0) {
                            (mGenomicDistance->operator[](i+j))->operator[](k).data = 0;
                        }
                        if ((mGenomicDistance->operator[](i+j))->operator[](k).data > maxValue[omp_get_thread_num()]) {
                            maxValue[omp_get_thread_num()] = (mGenomicDistance->operator[](i+j))->operator[](k).data;
                        }
                        if ((mGenomicDistance->operator[](i+j))->operator[](k).data < minValue[omp_get_thread_num()]) {
                            minValue[omp_get_thread_num()] = (mGenomicDistance->operator[](i+j))->operator[](k).data;
                        }
                        // std::cout << "RichcorrectedData: " << (mGenomicDistance->operator[](i+j))->operator[](k).data << std::endl;
                    }
                } 
                else {
                    for(uint32_t k = 0; k < (mGenomicDistance->operator[](i+j))->size(); ++k) {
                        // std::cout << "PoororgData: " << (mGenomicDistance->operator[](i+j))->operator[](k).data << std::endl;
                        
                        (mGenomicDistance->operator[](i+j))->operator[](k).data -= static_cast<int32_t>(changeValues[k]);
                        if ((mGenomicDistance->operator[](i+j))->operator[](k).data < 0) {
                            mGenomicDistance->operator[](i+j)->operator[](k).data = 0;
                        }
                        if ((mGenomicDistance->operator[](i+j))->operator[](k).data > maxValue[omp_get_thread_num()]) {
                            maxValue[omp_get_thread_num()] = (mGenomicDistance->operator[](i+j))->operator[](k).data;
                        }
                        if ((mGenomicDistance->operator[](i+j))->operator[](k).data < minValue[omp_get_thread_num()]) {
                            minValue[omp_get_thread_num()] = (mGenomicDistance->operator[](i+j))->operator[](k).data;
                        }
                        // std::cout << "PoorcorrectedData: " << (mGenomicDistance->operator[](i+j))->operator[](k).data << std::endl;
                        
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
    int max = 0; 
    // #pragma omp parallel num_threads(mNumberOfCores)
    for (uint32_t i = 0; i < mGenomicDistance->size(); ++i) {
        auto it = mGenomicDistance->begin();
        std::advance(it, i);
        // double interactionCount = 0;
        for (auto itVector = (it->second)->begin(); itVector != (it->second)->end(); ++itVector)  {
            // std::cout << "data: " << (*itVector).data << std::endl;
            // std::cout << (double) mMaxElement * (((double)(*itVector).data) /  maxValue) << std::endl;
            (*itVector).data = static_cast<int32_t>((float) mMaxElement * (((float)(*itVector).data) /  maximalValue));
            if ((*itVector).data > max) {
                max = (*itVector).data;
            }
        }
        // mGenomicDistanceMean->operator[]((it->first)) = interactionCount;
    }    
    double foo = 0.025;
    std::cout << "foo: " << foo << std::endl;;
    double exponent = 1.2;
    std::cout << "pow(800, -1.2): ";
    printf("%lf", std::pow(foo, -exponent));
    std::cout<< " pow(800, 1.2): ";
    printf("%lf", std::pow(foo, exponent));
    std::cout << std::endl; 

    std::cout << "mMaxElement: " << mMaxElement << std::endl;

    std::cout << "maxValue: " << maximalValue << std::endl;
    std::cout << "minValue: " << minimalValue << std::endl;
    
    std::cout << "maxNorm: " << max << std::endl;
    
}