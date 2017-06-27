#include "h5Interface.h"

H5Interface::H5Interface(char* pMatrixPath) {
    mMatrixPath = pMatrixPath;
}

H5Interface::~H5Interface() {
    // delete mMatrixPath;
}

std::unordered_map<uint32_t, std::vector<matrixElement>*>* H5Interface::readMatrix() {
    
    
    uint64_t test = 0;

    H5::H5File* file = new H5::H5File( mMatrixPath, H5F_ACC_RDONLY );
    H5::Group* group = new H5::Group (file->openGroup("matrix"));
    std::cout << __LINE__ << std::endl;

    // data vector
    H5::DataSet* dataH5 = new H5::DataSet(group->openDataSet("data"));
    hid_t id = dataH5->getId();
    hid_t dspace = H5Dget_space(id);
    hsize_t* dims = new hsize_t[1];
    hsize_t* maxdims = new hsize_t[1];
    int ndims = H5Sget_simple_extent_dims(dspace, dims, maxdims);
    // std::cout << "ndims: " << ndims << " dims: " << dims[0] << " maxdims: " << maxdims[0]<< std::endl;
    // how big and which datatype?
    size_t sizeData = dims[0];
    int64_t* data = new int64_t [sizeData];
    dataH5->read(data, H5::PredType::NATIVE_INT64);


    // indices vector
    H5::DataSet *indicesH5 = new H5::DataSet(group->openDataSet("indices"));
    id = indicesH5->getId();
    dspace = H5Dget_space(id);
    dims[0] = 0;
    maxdims[0] = 0;
    ndims = H5Sget_simple_extent_dims(dspace, dims, maxdims);
    size_t sizeIndices = dims[0];
   
    int32_t* indices = new int32_t [sizeIndices];
    indicesH5->read(indices, H5::PredType::NATIVE_INT32);
    
    std::cout << __LINE__ << std::endl;
    
    // index location of next line
    H5::DataSet *indptrH5 = new H5::DataSet(group->openDataSet("indptr"));
    id = indptrH5->getId();
    dspace = H5Dget_space(id);
    dims[0] = 0;
    maxdims[0] = 0;
    ndims = H5Sget_simple_extent_dims(dspace, dims, maxdims);
    size_t sizeIndptr = dims[0];
   
    int32_t* indptr = new int32_t [sizeIndptr];
    indptrH5->read(indptr, H5::PredType::NATIVE_INT32);
    std::cout << __LINE__ << std::endl;

    std::unordered_map<uint32_t, std::vector<matrixElement>*>* genomicDistanceMap = new std::unordered_map<uint32_t, std::vector<matrixElement>*>();
    for (uint32_t i = 0; i < sizeIndptr; ++i) {
        genomicDistanceMap->operator[](i) = new std::vector<matrixElement>();
    }
    std::cout << __LINE__ << std::endl;
    uint32_t genomicDistance = 0;
    for(uint32_t i = 0; i < sizeIndptr - 1; ++i) {
        for (uint32_t j = indptr[i]; j < indptr[j+1]; ++j) {

            matrixElement element;
            element.x = i;
            element.y = indices[j];
            element.data = data[j];

            genomicDistance = abs(element.x -  element.y);
            (genomicDistanceMap->operator[](genomicDistance))->push_back(element);
        }
    }

    
    std::cout << __LINE__ << std::endl;

    delete [] data;
    delete [] indices;
    delete [] indptr;
    delete dataH5;
    delete indicesH5;
    delete indptrH5;
    delete group;
    delete file;

    return genomicDistanceMap;
}

bool H5Interface::writeMatrix() {
    return false;
}