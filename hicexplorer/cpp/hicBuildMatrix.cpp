#include "hicBuildMatrix.h"
HicBuildMatrix::HicBuildMatrix() {

}
HicBuildMatrix::HicBuildMatrix(char* pSamFileOne, char* pSamFileTwo,
                                char* pOutputBam, size_t pBinSize,
                                char* pRestrictionCutFile, size_t pMinDistance,
                                size_t pMaxDistance, char* pRestrictionSequence,
                                char* pOutFileName, char* pRegion,
                                bool pRemoveSelfLigation, bool pRemoveSelfCircles,
                                size_t pMinMappingQuality, bool pSkipDuplicationCheck,
                                size_t pThreads) {
    
    mSamFileOne = pSamFileOne;
    mSamFileTwo = pSamFileTwo;
    mOutputBam = pOutputBam; 
    mBinSize = pBinSize;
    mRestrictionCutFile = pRestrictionCutFile; 
    mMinDistance = pMinDistance;
    mMaxDistance = pMaxDistance; 
    mRestrictionSequence = pRestrictionSequence;
    mOutFileName = pOutFileName;
    mRegion = pRegion;
    mRemoveSelfLigation = pRemoveSelfLigation;
    mRemoveSelfCircles = pRemoveSelfCircles;
    mMinMappingQuality = pMinMappingQuality;
    mSkipDuplicationCheck = pSkipDuplicationCheck;
    if (pThreads <= 0) {
        mThreads = omp_get_max_threads();
    } else {
        mThreads = pThreads;    
    }

}


HicBuildMatrix::~HicBuildMatrix() {

}
 
void HicBuildMatrix::buildMatrix() {

    std::vector<seqan::BamStream*> streamsFileOne(mThreads);
    std::vector<seqan::BamStream*> streamsFileTwo(mThreads);

    std::vector<seqan::BamAlignmentRecord*> recordFileOne(mThreads);
    std::vector<seqan::BamAlignmentRecord*> recordFileTwo(mThreads);
    
    
    for (size_t i = 0; i < mThreads; ++i) {
        streamsFileOne[i] = new seqan::BamStream(mSamFileOne);
        streamsFileTwo[i] = new seqan::BamStream(mSamFileTwo);
        recordFileOne[i] = new seqan::BamAlignmentRecord();
        recordFileTwo[i] = new seqan::BamAlignmentRecord();
    }
    // TODO: get bam file length and divide it by mThreads
    size_t chunkSize = 0;

    omp_set_dynamic(0);
    #pragma omp parallel for num_threads(mThreads)
    for (size_t i = OMP_GET_THREAD_NUM() * chunkSize; 
            i <  OMP_GET_THREAD_NUM() * (chunkSize + 1) 
            && !atEnd(*(streamsFileOne[OMP_GET_THREAD_NUM()])) && !atEnd(*(streamsFileTwo[OMP_GET_THREAD_NUM()])); ++i) {
              
                readRecord(*(recordFileOne[OMP_GET_THREAD_NUM()]), *(streamsFileOne[OMP_GET_THREAD_NUM()]));
                readRecord(*(recordFileTwo[OMP_GET_THREAD_NUM()]), *(streamsFileTwo[OMP_GET_THREAD_NUM()]));

    }




    std::cout << mSamFileOne << std::endl;
 // Open input stream, BamStream can read SAM and BAM files.
    seqan::BamStream bamStreamIn(mSamFileOne);
    // Open output stream, "-" means stdin on if reading, else stdout.
    seqan::BamStream bamStreamOut("-", seqan::BamStream::WRITE);
    // Copy header.  The header is automatically written out before
    // the first record.
    bamStreamOut.header = bamStreamIn.header;

    seqan::BamAlignmentRecord record;
    while (!atEnd(streamsFileOne[]))
    {
        readRecord(record, bamStreamIn);
        writeRecord(bamStreamOut, record);
    }
}