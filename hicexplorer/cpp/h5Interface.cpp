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
    H5::H5File* fileOld = new H5::H5File( mMatrixPath, H5F_ACC_RDONLY );
    H5::Group* groupIntervalsOld = new H5::Group (fileOld->openGroup("/intervals"));


    H5::Group* groupIntervals = new H5::Group(file->createGroup("/intervals"));
    this->createAttribute(groupIntervals, "CLASS", "GROUP");
    this->createAttribute(groupIntervals, "TITLE", "");
    this->createAttribute(groupIntervals, "VERSION", "1.0");
    

    // get data set '/intervals/end_list'
    H5::DataSet* data_end_list = new H5::DataSet(groupIntervalsOld->openDataSet("end_list"));
    hid_t id = data_end_list->getId();
    hid_t dspace = H5Dget_space(id);
    hsize_t* maxdims = new hsize_t[1];
    int ndims = H5Sget_simple_extent_dims(dspace, dims, maxdims);
    size_t sizeData = dims[0];
    int64_t* end_list = new int64_t [sizeData];
    data_end_list->read(end_list, H5::PredType::NATIVE_INT64);

    H5::DataSpace* dataspace_end_list = new H5::DataSpace(RANK, dims);
    // create dataset
    H5::DataSet* dataset_end_list = new H5::DataSet(file->createDataSet("/intervals/end_list", H5::PredType::NATIVE_INT64,
                                                *dataspace_end_list));
    dataset_end_list->write(end_list, H5::PredType::NATIVE_INT64);
    // Create attributes
    this->createAttribute(dataset_end_list, "CLASS", "CARRAY");
    this->createAttribute(dataset_end_list, "TITLE", "");
    this->createAttribute(dataset_end_list, "VERSION", "1.1");
    delete end_list;
    delete dataset_end_list;
    delete dataspace_end_list;
    delete data_end_list;



    // get data set '/intervals/extra_list'
    H5::DataSet* data_extra_list = new H5::DataSet(groupIntervalsOld->openDataSet("extra_list"));
    id = data_extra_list->getId();
    dspace = H5Dget_space(id);
    // hsize_t* dims = new hsize_t[1];
    // maxdims = new hsize_t[1];
    ndims = H5Sget_simple_extent_dims(dspace, dims, maxdims);
    sizeData = dims[0];
    int64_t* extra_list = new int64_t [sizeData];
    data_end_list->read(extra_list, H5::PredType::NATIVE_INT64);

    H5::DataSpace* dataspace_extra_list = new H5::DataSpace(RANK, dims);
    // create dataset
    H5::DataSet* dataset_extra_list = new H5::DataSet(file->createDataSet("/intervals/extra_list", H5::PredType::NATIVE_INT64,
                                                *dataspace_extra_list));
    dataset_extra_list->write(end_list, H5::PredType::NATIVE_INT64);
    // Create attributes
    this->createAttribute(dataset_extra_list, "CLASS", "CARRAY");
    this->createAttribute(dataset_extra_list, "TITLE", "");
    this->createAttribute(dataset_extra_list, "VERSION", "1.1");
    delete extra_list;
    delete dataset_extra_list;
    delete dataspace_extra_list;
    delete data_extra_list;
    
    H5::DataSet* data_start_list = new H5::DataSet(groupIntervalsOld->openDataSet("start_list"));
    id = data_start_list->getId();
    dspace = H5Dget_space(id);
    ndims = H5Sget_simple_extent_dims(dspace, dims, maxdims);
    sizeData = dims[0];
    std::cout << "sizeData foo: " << sizeData << std::endl;
    int64_t* start_list = new int64_t [sizeData];
    data_start_list->read(extra_list, H5::PredType::NATIVE_INT64);


    H5::DataSpace* dataspace_start_list = new H5::DataSpace(RANK, dims);
    // create dataset
    H5::DataSet* dataset_start_list = new H5::DataSet(file->createDataSet("/intervals/start_list", H5::PredType::NATIVE_INT64,
                                                *dataspace_start_list));
    dataset_start_list->write(start_list, H5::PredType::NATIVE_INT64);
    // Create attributes
    this->createAttribute(dataset_start_list, "CLASS", "CARRAY");
    this->createAttribute(dataset_start_list, "TITLE", "");
    this->createAttribute(dataset_start_list, "VERSION", "1.1");
    delete start_list;
    delete dataspace_start_list;
    delete dataset_start_list;
    delete data_start_list;


    // get data set '/intervals/chr_list'
    H5::DataSet* data_chr_list = new H5::DataSet(groupIntervalsOld->openDataSet("chr_list"));
    id = data_chr_list->getId();
    dspace = H5Dget_space(id);
    ndims = H5Sget_simple_extent_dims(dspace, dims, maxdims);
    sizeData = dims[0];
    std::cout << "ndims" << ndims << std::endl;
    std::cout << "sizeData: " << sizeData << std::endl;
    // char* rdata = new char[sizeData*22];
    H5::DataType dtype = data_chr_list->getDataType();
    // std::cout << "dataType: " << dtype.getId() << std::endl;
    // // char** chr_list = new char* [sizeData];
    // char** chr_list = new char* [sizeData];
    
    // data_chr_list->read(chr_list, dtype);
    // std::cout << "foo" << std::endl;
    char *chr_list = NULL;
    H5::StrType vlst(0, H5T_VARIABLE);
    // Read and verify the dataset string as a string of chars.
    data_chr_list->read(&chr_list, vlst);
    // HDfree(chr_list); 

    H5::DataSpace* dataspace_chr_list = new H5::DataSpace(RANK, dims);
    // // create dataset
    H5::DataSet* dataset_chr_list = new H5::DataSet(file->createDataSet("/intervals/chr_list", dtype,
                                                *dataspace_chr_list));
    dataset_chr_list->write(chr_list, dtype);
    // // Create attributes
    this->createAttribute(dataset_chr_list, "CLASS", "CARRAY");
    this->createAttribute(dataset_chr_list, "TITLE", "");
    this->createAttribute(dataset_chr_list, "VERSION", "1.1");
    // delete chr_list;
    // delete dataset_chr_list;
    // delete dataspace_chr_list;
    // delete data_chr_list;

    delete maxdims;
    // delete ndims;
    delete fileOld;
    delete groupIntervalsOld;
    delete groupIntervals;
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
