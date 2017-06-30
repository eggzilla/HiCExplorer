#ifndef H5_INTERFACE_H
#define H5_INTERFACE_H

#include "typeDefinitions.h"

#include <H5Cpp.h>
#include <hdf5.h>
// #include <H5File.h>
// #include <H5DataSet.h>

class H5Interface {
    private:
        char* mMatrixPath = NULL;
        char* mMatrixPathOutput = NULL;
    
    public:
        H5Interface(char* pMatrixPath, char* pMatrixPathOutput);
        ~H5Interface();
        std::unordered_map<uint32_t, std::vector<matrixElement>*>* readMatrix(uint32_t pRemoveLowInteractionCount);
        void writeMatrix(std::unordered_map<uint32_t, std::vector<matrixElement>*>* pGenomicDistanceMap, 
                                uint32_t pRemoveLowInteractionCount);
        void createAttribute(H5::H5Object* pDataset, std::string pName, std::string pValue);

};

#endif // H5_INTERFACE_H
