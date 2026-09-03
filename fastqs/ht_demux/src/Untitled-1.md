
# Instructions to write a custom FASTQ parsing script

## Script objective

I want to write a python script to split and combine FASTQ reads based on 4 nucleotide (nt) long  barcodes.  Here are my files:

expRW086_pool_10_S10_L001_1.fastq.gz  expRW086_pool_11_S11_L001_2.fastq.gz  expRW086_pool_12_S12_L002_1.fastq.gz
expRW086_pool_03_S3_L001_1.fastq.gz  expRW086_pool_10_S10_L001_2.fastq.gz  expRW086_pool_11_S11_L002_1.fastq.gz  expRW086_pool_12_S12_L002_2.fastq.gz
expRW086_pool_03_S3_L001_2.fastq.gz  expRW086_pool_10_S10_L002_1.fastq.gz  expRW086_pool_11_S11_L002_2.fastq.gz  expRW086_pool_12_S12_L003_1.fastq.gz
expRW086_pool_03_S3_L002_1.fastq.gz  expRW086_pool_10_S10_L002_2.fastq.gz  expRW086_pool_11_S11_L003_1.fastq.gz  expRW086_pool_12_S12_L003_2.fastq.gz
expRW086_pool_03_S3_L002_2.fastq.gz  expRW086_pool_10_S10_L003_1.fastq.gz  expRW086_pool_11_S11_L003_2.fastq.gz  expRW086_pool_12_S12_L004_1.fastq.gz
expRW086_pool_03_S3_L003_1.fastq.gz  expRW086_pool_10_S10_L003_2.fastq.gz  expRW086_pool_11_S11_L004_1.fastq.gz  expRW086_pool_12_S12_L004_2.fastq.gz
expRW086_pool_03_S3_L003_2.fastq.gz  expRW086_pool_10_S10_L004_1.fastq.gz  expRW086_pool_11_S11_L004_2.fastq.gz
expRW086_pool_03_S3_L004_1.fastq.gz  expRW086_pool_10_S10_L004_2.fastq.gz  expRW086_pool_12_S12_L001_1.fastq.gz
expRW086_pool_03_S3_L004_2.fastq.gz  expRW086_pool_11_S11_L001_1.fastq.gz  expRW086_pool_12_S12_L001_2.fastq.gz

Files  with "pool_03" in the name contain only Control data. Files with "pool_10", "pool_11, and "pool_12" contain time point samples from 3 different enzymes samples "WT_SpCas9", "SpG", and "SpRY".

Each sample  contains a 65 nt forward read and a 10 nt reverse read. The reverse read contains only the P5 barcode no data.  Once it is used to split the 
R1 read into the right sample category is can be ignored.

Each pool was sequences across 4 lanes L001 to L004. These can be combined.

The P7 barcode appears at 0-index position [2:6] of the Forward read (_1) the P5 barcode appears at position [2:6] of the  10 nt reverse read (_2).

Here is the list of barcodes and how they map to samples:

sample,description,P5_sample_barcode,P7_sample_barcode
sample_01,WT_SpCas9,ATGC,AGTA
sample_02,SpG,ATGC,ACTT
sample_03,SpRY,GCAA,TCCA
control_sample,untreated_library,ATGC,TCGC

For the Control samples (pool_03)  each forward and reverse read should be checked to see if it matches its barcodes and R1 reads that pass should be saved to a file in the output dir "expRW086_control_1.fastq.gz"

For the experimental data in each file set,  "pool_10", "pool_11, and "pool_12", each read should be checked to see if its forward AND reverse barcodes match each of the three sample pairs. If reads pairs match the Forward reads, that match should be put in their correct file, for example expRW086_SpCas9_pool_10_1.fastq.gz , expRW086_SpCas9_pool_10_1.fastq.gz  or expRW086_SpY_pool_10_1.fastq.gz.  In the end the ouptut dir should have 10 files. 

## Processing script:

Write a Rust script that uses bio-rust for the fastq file reading and writing.  The script should have a CLI interface that takes an --input input directory and an --output output directory (created if not present) and --threads, the number of independent processing threads to process each pool task on  (up to 4).
it should do parallel file compression perhaps with a tool like pzip if that is not built into bio-rust.  Also provide the Crate.toml to define dependencies. The script should output how many reads were put into each pool/sample category successfully out of the total number of reads.

If you have any question about the script please ask before attempting to produce it.
