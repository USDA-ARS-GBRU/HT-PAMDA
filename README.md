# Reanalysis of HT-PAMDA

For original repo readme see `Original_README.md`

This fork makes modification needed to collect data and run the original HT-PAMDA, 5 years after publication.

## Data availability.

Data in the `fastqs` folder are very small.

In the `example_PAM_library` each files has this many lines:
```
  266792
  266792
  265048
  265048
  272956
  272956
  269088
  269088
```
and many of the reads are truncated for example 

```
+
AAAAAEEEEE
@NB500929:497:HVYT5BGXB:2:11101:2933:1389 2:N:0:CGCTCATT+AGGATAGG
TGTCGCCGGT
```

for `example_PAMDA_data` we see this many lines:
```
  452616
  452616
  452656
  452656
  453096
  453096
  451920
  451920
  487576
  487576
  488676
  488676
  488832
  488832
  489908
  489908
  451344
  451344
  451916
  451916
  451880
  451880
  453716
  453716
```

Again with many truncated reads.

IF we look in SRA we can find these samples. for example [expRW086_pool_03_S3_L004](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR11182586&display=download)

But when we download that test file `SRR1118258.fastq` it is different from the example file. it has 5966928 rather than 269088 lines, 22x less. The sampling seems somewhat random.

## Attempting to run the example data

#Setting up a conda env for the code. 

the repo has these requirements:

```
matplotlib==3.1.3
numpy==1.18.1
pandas==1.0.3
scipy==1.4.1
tqdm==4.46.0
seaborn==0.10.1
```

Copilot suggest Python 3.6 for this set of dependencies. I'll try that first. 

```
CONDA_SUBDIR=osx-64 conda create -n htpamdaenv -c conda-forge python=3.6
conda activate htpamdaenv
pip install -r requirements.txt 
```

# running examplecode 

```
python3 library_QC.py
```
this runs with this output

```(htpamdaenv) ➜  code git:(master) python3 library_QC.py
Begin library QC
100.0% of reads mapped from expRW086_pool_03_S3 (66698 reads)                                                                           
100.0% of reads mapped from expRW086_pool_03_S3 (68239 reads)                                                                           
100.0% of reads mapped from expRW086_pool_03_S3 (67272 reads)                                                                           
100.0% of reads mapped from expRW086_pool_03_S3 (66262 reads)                                                                           
fastq files: 100%|████████████████████████████████████████████████████████████████████████████████████████| 4/4 [00:04<00:00,  1.02s/it]
reads: : 0it [00:00, ?it/s]       
writing compressed CSV output
summarizing raw read counts
grouping counts by indicated PAM bases
PAM library SPACER1 max:min ratio: 2.7975
PAM library SPACER1 90:10 ratio: 1.4959
PAM library SPACER1 skewness: -0.0559
PAM library SPACER2 max:min ratio: 2.5042
PAM library SPACER2 90:10 ratio: 1.522
PAM library SPACER2 skewness: -0.061
Library QC complete
```



```
time python3 fastq2count.py
```

```
100.0% of reads mapped from expRW086_pool_10_S10 (112980 reads)                                                                         
100.0% of reads mapped from expRW086_pool_11_S11 (122208 reads)                                                                         
100.0% of reads mapped from expRW086_pool_11_S11 (122477 reads)                                                                         
100.0% of reads mapped from expRW086_pool_10_S10 (113274 reads)                                                                         
100.0% of reads mapped from expRW086_pool_12_S12 (112979 reads)                                                                         
100.0% of reads mapped from expRW086_pool_10_S10 (113164 reads)                                                                         
100.0% of reads mapped from expRW086_pool_12_S12 (112970 reads)                                                                         
100.0% of reads mapped from expRW086_pool_12_S12 (113429 reads)                                                                         
100.0% of reads mapped from expRW086_pool_11_S11 (122169 reads)                                                                         
100.0% of reads mapped from expRW086_pool_12_S12 (112835 reads)                                                                         
100.0% of reads mapped from expRW086_pool_11_S11 (121894 reads)                                                                         
100.0% of reads mapped from expRW086_pool_10_S10 (113154 reads)                                                                         
fastq files: 100%|██████████████████████████████████████████████████████████████████████████████████████| 12/12 [00:21<00:00,  1.83s/it]
reads: : 0it [00:00, ?it/s]       
writing compressed CSV output
summarizing raw read counts
python fastq2count.py  25.37s user 5.12s system 118% cpu 25.633 total
```

```
tyme python rawcount2normcount.py
```

```
(htpamdaenv) ➜  code git:(master) time python rawcount2normcount.py
grouping counts by indicated PAM bases
sample:: 100%|████████████████████████████████████████████████████████████████████████████████████████████| 4/4 [00:00<00:00,  8.63it/s]
determining most enriched PAMs per sample
samples: 100%|████████████████████████████████████████████████████████████████████████████████████████████| 3/3 [00:00<00:00,  9.80it/s]
normalizing each sample:
samples: 100%|████████████████████████████████████████████████████████████████████████████████████████████| 3/3 [00:01<00:00,  1.57it/s]
python rawcount2normcount.py  4.83s user 4.79s system 222% cpu 4.318 total
```

```
time python3 normcount2rate.py
```
output:
```
time python3 normcount2rate.py
calculating rate constants
samples: 100%|████████████████████████████████████████████████████████████████████████████████████████████| 3/3 [00:01<00:00,  2.95it/s]
appending rate constants
output to CSV
python3 normcount2rate.py  2.62s user 4.63s system 343% cpu 2.107 total
```

```
time python rate2heatmap.py
```
Output:
```
time python rate2heatmap.py
samples: 100%|████████████████████████████████████████████████████████████████████████████████████████████| 3/3 [00:01<00:00,  1.55it/s]
python rate2heatmap.py  5.33s user 7.19s system 375% cpu 3.334 total
```

