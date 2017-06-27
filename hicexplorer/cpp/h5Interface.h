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
    
    public:
        H5Interface(char* pMatrixPath);
        ~H5Interface();
        std::unordered_map<uint32_t, std::vector<matrixElement>*>* readMatrix();
        bool writeMatrix();

};

#endif // H5_INTERFACE_H
