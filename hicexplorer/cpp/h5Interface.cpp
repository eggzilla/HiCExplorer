#include "h5Interface.h"

H5Interface::H5Interface(char* pMatrixPath, char* pMatrixPathOutput) {
    mMatrixPath = pMatrixPath;
    mMatrixPathOutput = pMatrixPathOutput;
}

H5Interface::~H5Interface() {
    // delete mMatrixPath;
}

std::unordered_map<uint32_t, std::vector<matrixElement>*>* H5Interface::readMatrix(uint32_t pRemoveLowInteractionCount) {
    
    
    // uint64_t test = 0;

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
            if (data[j] <= pRemoveLowInteractionCount) {
                continue;
            }
            matrixElement element;
            element.x = i;
            element.y = indices[j];
            element.data = data[j];

            genomicDistance = abs(element.x -  element.y);
            (genomicDistanceMap->operator[](genomicDistance))->push_back(element);
        }
    }

    
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

void H5Interface::writeMatrix(std::unordered_map<uint32_t, std::vector<matrixElement>*>* pGenomicDistanceMap, 
                                uint32_t pRemoveLowInteractionCount) {

    // parse data from computation structure to csr storage structure with indptr, indices and data array                                
    std::vector<int32_t> indptrVector;
    std::vector<int32_t> indicesVector;
    std::vector<int64_t> dataVector;
    
    int32_t counter = 0; 
    int32_t xOld = -1;
    for (auto it = pGenomicDistanceMap->begin(); it != pGenomicDistanceMap->end(); ++it) {
        for (auto itVector = (it->second)->begin(); itVector != (it->second)->end(); ++itVector) {
            if ((*itVector).data <= pRemoveLowInteractionCount) {
                continue;
            }
            if (xOld != (*itVector).x) {
                indptrVector.push_back(counter);
                xOld = (*itVector).x;
            }
            indicesVector.push_back((*itVector).y);
            dataVector.push_back((*itVector).data);
            ++counter;          
        }
    } 
    indptrVector.push_back(counter);

    int32_t* indptr = &indptrVector[0];
    int32_t* indices = &indicesVector[0];
    int64_t* data = &dataVector[0];

    
    // create new h5 file
    H5::H5File* file = new H5::H5File(mMatrixPathOutput, H5F_ACC_TRUNC);
    // create group '/matrix'
    H5::Group* groupMatrix = new H5::Group(file->createGroup("/matrix"));
    // create attributes
    this->createAttribute(groupMatrix, "CLASS", "GROUP");
    this->createAttribute(groupMatrix, "TITLE", "");
    this->createAttribute(groupMatrix, "VERSION", "1.0");
    
    // create 3 dataspaces: /matrix/indptr, /matrix/indices, /matrix/data
    int RANK = 1;
    // indptr
    hsize_t dims[1];
    dims[0] = indptrVector.size();
    H5::DataSpace* dataspaceIndptr = new H5::DataSpace(RANK, dims);
    // create dataset
    H5::DataSet* datasetIndptr = new H5::DataSet(file->createDataSet("/matrix/indptr", H5::PredType::NATIVE_INT32,
                                                *dataspaceIndptr));
    datasetIndptr->write(indptr, H5::PredType::NATIVE_INT32);

    // Create attributes
    this->createAttribute(datasetIndptr, "CLASS", "CARRAY");
    this->createAttribute(datasetIndptr, "TITLE", "");
    this->createAttribute(datasetIndptr, "VERSION", "1.1");

    
    // indices
    dims[0] = indicesVector.size();
    H5::DataSpace* dataspaceIndices = new H5::DataSpace(RANK, dims);
    // create dataset
    H5::DataSet* datasetIndices = new H5::DataSet(file->createDataSet("/matrix/indices", H5::PredType::NATIVE_INT32,
                                                *dataspaceIndices));
    datasetIndices->write(indices, H5::PredType::NATIVE_INT32);
    
    // Create attributes
    this->createAttribute(datasetIndices, "CLASS", "CARRAY");
    this->createAttribute(datasetIndices, "TITLE", "");
    this->createAttribute(datasetIndices, "VERSION", "1.1");
    

    // data
    dims[0] = dataVector.size();
    H5::DataSpace* dataspaceData = new H5::DataSpace(RANK, dims);
    // create dataset
    H5::DataSet* datasetData = new H5::DataSet(file->createDataSet("/matrix/data", H5::PredType::NATIVE_INT64,
                                                *dataspaceData));
    datasetData->write(data, H5::PredType::NATIVE_INT64);

    // Create attributes
    this->createAttribute(datasetData, "CLASS", "CARRAY");
    this->createAttribute(datasetData, "TITLE", "");
    this->createAttribute(datasetData, "VERSION", "1.1");

    delete dataspaceIndices;
    delete datasetIndices;

    delete dataspaceIndptr;
    delete datasetIndptr;

    delete dataspaceData;
    delete datasetData;

    
    // copy old content of /intervals
    // create group intervals
    // create dataspaces for 
    // H5::H5File* fileOld = new H5::H5File( mMatrixPath, H5F_ACC_RDONLY );
    // H5::Group* groupIntervalsOld = new H5::Group (file->openGroup("intervals"));

    // H5::Group* groupIntervalsCopy = new H5::Group(groupIntervalsOld);
    // H5::Group* groupIntervalsNew = new H5::Group(file->createGroup(groupIntervalsNew->getId()));

    
    // delete groupIntervalsNew;
    // delete groupIntervalsCopy;
    // delete groupIntervalsOld;
    // delete fileOld;
    delete file;
    
    return;
}

void H5Interface::createAttribute(H5::H5Object* pDataset, std::string pName, std::string pValue) {
    
    H5::StrType str_type(0, H5T_VARIABLE);
    H5::DataSpace att_space(H5S_SCALAR);
    H5::Attribute att = pDataset->createAttribute( pName, str_type, att_space );
    att.write( str_type, pValue );
    return;
}
