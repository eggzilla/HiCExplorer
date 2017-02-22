#include <stdio.h>
#include "hicBuildMatrix.h"
int main (char** args) {

    char* pSamFileOne = "test-data/hic.bam";
    char* pSamFileTwo ="test-data/hic.bam";
    char* pOutputBam ="test-data/hic.bam";
    size_t pBinSize = 1;
    char* pRestrictionCutFile = "test-data/hic.bam";
    size_t pMinDistance = 1 ;
    size_t pMaxDistance = 1;
    char* pRestrictionSequence = "GATC";
    char* pOutFileName = "bal.bam";
    char* pRegion = "chr2L";
    bool pRemoveSelfLigation = true;
     bool pRemoveSelfCircles = false;
    size_t pMinMappingQuality = 1;
    bool pSkipDuplicationCheck = false;

    HicBuildMatrix* hic = new HicBuildMatrix(pSamFileOne, pSamFileTwo,
                     pOutputBam,  pBinSize,
                    pRestrictionCutFile, pMinDistance,
                    pMaxDistance,  pRestrictionSequence,
                    pOutFileName,  pRegion,
                pRemoveSelfLigation,  pRemoveSelfCircles,
                    pMinMappingQuality,  pSkipDuplicationCheck);
    hic->buildMatrix();
    delete hic;
    return 0;
}