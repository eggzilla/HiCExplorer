#include <stdio.h>
#include <iostream>
#include <seqan/bam_io.h>

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
    void get_bins(int bin_size, int chrom_size, int region);
    void bed2interval_list(int bed_file_handler);
    void get_rf_bins(int rf_cut_intervals, int min_distance,int max_distance);
    void get_chrom_sizes(int bam_handle);
    void check_dangling_end(int read, int dangling_sequences);
    void get_supplementary_alignment(int read, int pysam_obj);
    void get_correct_map(int primary, int supplement_list);
    void enlarge_bins(int bin_intervals, int chrom_sizes);
};


#endif // HIC_BUILD_MATRIX