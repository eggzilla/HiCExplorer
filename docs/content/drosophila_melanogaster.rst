Hi-C analysis of drosophila melanogaster using HiC-Explorer
===========================================================

This Howto shows a typical HiC data analysis with HiC-Explorer. A drosophila melanogaster dataset is analyzed.
All needed and used data including intermediate files are stored at zenodo.

Prepare for analysis
--------------------

Please download R1.fastq and R2.fastq

Build Index
------------
In this step we build an index for the mapping.

bowtie
^^^^^^

.. code:: bash

    bowtie2-build dm6.fa dm6_index_bt2/dm6 --threads 8

http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml


hisat
^^^^^

.. code:: bash

    hisat2-build dm6.fa dm6_index_ht2/dm6 -p 8
    

https://ccb.jhu.edu/software/hisat2/manual.shtml#the-hisat2-build-indexer

star(?)
^^^^^^^


Mapping
-------

bowtie
hisat
star

Build HiC Matrix
----------------


Correct data
------------


Plot Matrix
-----------

Merge bins?


Find TADs
---------


Plot TADs
---------