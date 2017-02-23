#include <stdio.h>
#include <iostream>
#include <seqan/bam_io.h>
#include <assert.h>

#include <omp.h>

#ifndef HIC_BUILD_MATRIX
#define HIC_BUILD_MATRIX

class HicBuildMatrix {
    private:
        char* mSamFileOne;
        char* mSamFileTwo;
        char* mOutputBam;
        size_t mBinSize;
        char* mRestrictionCutFile;
        size_t mMinDistance;
        size_t mMaxDistance;
        char* mRestrictionSequence;
        char* mOutFileName;
        char* mRegion;
        bool mRemoveSelfLigation;
        bool mRemoveSelfCircles;
        size_t mMinMappingQuality;
        bool mSkipDuplicationCheck;
        size_t mThreads;
        
    public:
    HicBuildMatrix();
    HicBuildMatrix(char* pSamFileOne, char* pSamFileTwo,
                    char* pOutputBam, size_t pBinSize,
                    char* pRestrictionCutFile, size_t pMinDistance,
                    size_t pMaxDistance, char* pRestrictionSequence,
                    char* pOutFileName, char* pRegion,
                    bool pRemoveSelfLigation, bool pRemoveSelfCircles,
                    size_t pMinMappingQuality, bool pSkipDuplicationCheck, size_t pThreads);
    ~HicBuildMatrix();
    void buildMatrix();
    void intervalListToIntervalTree(int interval_list);
    void getBins(int bin_size, int chrom_size, int region);
    void bedToIntervalVector(int bed_file_handler);
    void getRfBins(int rf_cut_intervals, int min_distance,int max_distance);
    void getChromSizes(int bam_handle);
    void checkDanglingEnd(int read, int dangling_sequences);
    std::vector<seqan::BamAlignmentRecord*> getSupplementaryAlignment(seqan::BamAlignmentRecord* pBamAlignmentRecord, seqan::BamStream* pBamStream);
    void getCorrectMap(int primary, int supplement_list);
    void enlargeBins(int bin_intervals, int chrom_sizes);
};


#endif // HIC_BUILD_MATRIX