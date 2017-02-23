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

std::vector<seqan::BamAlignmentRecord*> HicBuildMatrix::getSupplementaryAlignment(seqan::BamAlignmentRecord* pBamAlignmentRecord, seqan::BamStream* pBamStream) {

}

void HicBuildMatrix::buildMatrix() {
    seqan::BamStream fileOne(mSamFileOne);
    seqan::BamStream fileTwo(mSamFileTwo);
    
    seqan::BamAlignmentRecord* recordFileOne;
    seqan::BamAlignmentRecord* recordFileTwo;
    std::vector<seqan::BamAlignmentRecord*> recordsBufferOne(1e6);
    std::vector<seqan::BamAlignmentRecord*> recordsBufferTwo(1e6); 
    std::vector<seqan::BamAlignmentRecord*> mateOneSupplementaryVector();
    std::vector<seqan::BamAlignmentRecord*> mateTwoSupplementaryVector();
    
    size_t lineCount = 0;
    while (!atEnd(fileOne) && !atEnd(fileTwo)) {
        recordFileOne = new seqan::BamAlignmentRecord();
        recordFileTwo = new seqan::BamAlignmentRecord();
        
        
        readRecord(*(recordFileOne), fileOne);
        readRecord(*(recordFileTwo), fileTwo);

        while (*(recordFileOne).flag & 256 == 256) {
            if (!atEnd(fileOne)) {
                readRecord(*(recordFileOne), fileOne);
            } else {
                break;
            }
        }
        while (*(recordFileTwo).flag & 256 == 256) {
            if (!atEnd(fileTwo)) {
                readRecord(*(recordFileTwo), fileTwo);
            } else {
                break;
            }
        } 
        assert (*(recordFileOne).qName == *(recordFileTwo).qName);

        mateOneSupplementaryVector.clear();
        mateTwoSupplementaryVector.clear();

        mateOneSupplementaryVector = getSupplementaryAlignment();
        mateTwoSupplementaryVector = getSupplementaryAlignment();
        lineCount++;
        // read a million lines and process them in parallel
        if (lineCount % 1e6 == 0) {

        }
    }



}