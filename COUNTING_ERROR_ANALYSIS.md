# HT-PAMDA Counting & Demultiplexing Error Analysis

This document provides a detailed technical analysis of four undisputable errors identified in the HT-PAMDA pipeline codebase (`code/PAMDA.py`).

Each error is analyzed below with code examples, mathematical impacts, and specific evaluations of whether and how each bug impacts the provided example dataset (`fastqs/example_PAMDA_data`).

---

## 1. P7 Barcode Length Variable Assignment Bug (`fastq2count`)

### Description
In `code/PAMDA.py` (`fastq2count`), line 203 assigns `P7_sample_BC_len` using `P5_sample_BCs[0]` instead of `P7_sample_BCs[0]`.

### Code Location
`code/PAMDA.py` lines 201–203:
```python
P5_sample_BCs = variant_ids['P5_sample_barcode'].tolist()
P5_sample_BC_len = len(P5_sample_BCs[0])
P7_sample_BCs = variant_ids['P7_sample_barcode'].tolist()
P7_sample_BC_len = len(P5_sample_BCs[0])  # BUG: Should be P7_sample_BCs[0]
```

### Impact & Example Dataset Context
- **General Impact:** If custom P5 and P7 barcodes have different lengths (e.g., 4-nt P5 barcode and 6-nt P7 barcode), `P7_sample_BC_len` is incorrectly set to 4 instead of 6. `find_BCs_and_spacer` then slices only 4 bases from Read 1 for `P7_sample_BC`. As a result, no read ever matches the 6-base P7 barcode in `variant_dict`, leading to **100% unmapped reads (`wrong_barcode`)**.
- **Impact on Provided Example Dataset:** In the provided example files (`barcode_csv/example_PAMDA_barcodes.csv`), all P5 barcodes (`ATGC`, `GCAA`) are 4 nt long, and all P7 barcodes (`AGTA`, `ACTT`, `TCCA`, `TCGC`) are also 4 nt long. Because both P5 and P7 barcodes happen to be 4 nt long, `len(P5_sample_BCs[0]) == len(P7_sample_BCs[0]) == 4`. Thus, **this bug has no impact on the provided example dataset**, but remains a critical flaw if users attempt to use barcode sets with unequal P5 and P7 lengths.

---

## 2. Timepoint Index Mismatch & `KeyError` on Control Timepoint Index 0 (`rawcount2normcount`)

### Description
There is a fundamental column indexing mismatch between how raw counts are saved and how `rawcount2normcount` attempts to look up control sample counts when `CONTROL_RAW_COUNT_CSV` is `None` and `CONTROL_SAMPLE_TIMEPOINT_FASTQ` is provided.

### Code Location
In `fastq2count` (lines 213, 290), count columns are written as `Raw_Counts_1`, `Raw_Counts_2`, ..., `Raw_Counts_N`.

In `rawcount2normcount` (lines 306–309):
```python
for index, row in df.iterrows():
    if row['Sample'] == control_sample:
        control_dict[row['Spacer']][row['PAM']] = row['Norm_Counts_' +
                                                      str(control_sample_timepoint_fastq)]
```

### Impact & Likelihood on Provided Dataset
- **Workflow Context in Example:** In the default example execution, `inputs.py` uses `CONTROL_RAW_COUNT_CSV = 'output/example_PAM_library/PAMDA_1_raw_counts.csv.gz'`, meaning the untreated control library was run separately via `library_QC.py`. When `CONTROL_RAW_COUNT_CSV` is provided, line 302 hardcodes `control_sample_timepoint_fastq = 1`, reading `Norm_Counts_1` from the separate control file.
- **Alternative Documented Setup (Multiplexed Control):** The pipeline documentation explicitly advertises that users can pool the control sample into a timepoint FASTQ by setting `CONTROL_RAW_COUNT_CSV = None` and specifying `CONTROL_SAMPLE_TIMEPOINT_FASTQ`.
- **Failure Analysis:** If a user follows this documented mode and sets `CONTROL_SAMPLE_TIMEPOINT_FASTQ = 0` (indicating the control is pooled in timepoint 0):
  1. The loop executes line 309: `row['Norm_Counts_0']`.
  2. Because `Norm_Counts_0` is created **after** this loop on line 313 (`df['Norm_Counts_0'] = norm_counts_0`), the script **100% crashes with a `KeyError: 'Norm_Counts_0'`**.
  3. Furthermore, even if a user sets `CONTROL_SAMPLE_TIMEPOINT_FASTQ = 1` to bypass the crash, `Norm_Counts_1` in `df` corresponds to experimental timepoint 1 ($t = 60\text{s}$), NOT $t = 0\text{s}$. Counts from $t = 60\text{s}$ are incorrectly used as $t_0$ baseline abundance.

---

## 3. Overwritten Top-N Enrichment Normalization (`rawcount2normcount`)

### Description
In `rawcount2normcount`, the code fits linear regressions to select the top $N$ most enriched (non-cleaved/uptrending) PAMs per sample, calculates correction factors (`uptrend_corrections`), but then **completely overwrites and discards** the correction step on the next line of code.

### Code Location
`code/PAMDA.py` lines 343–346:
```python
for index, row in df.iterrows():
    ...
    for i in range(len(timepoints)):
        df.loc[index, 'Norm_Counts_' + str(i)] = row['Norm_Counts_' + str(i)] / \
                                                 uptrend_corrections[sample_spacer][i]
        df.loc[index, 'Norm_Counts_' + str(i)] = row['Norm_Counts_' + str(i)] / \
                                                 row['Norm_Counts_0']
```

### Mathematical Impact
- **Nullification of Linear Fits:** Yes, this bug **completely nullifies** the linear fits and top-$N$ uptrend PAM selection.
- **Why it happens:** `row` is an immutable pandas `Series` snapshot returned by `df.iterrows()`. Line 344 assigns the enrichment-corrected value to `df.loc[index, 'Norm_Counts_' + str(i)]`. However, on line 345, `row['Norm_Counts_' + str(i)]` STILL holds the **original uncorrected value**. Line 345 computes `original / row['Norm_Counts_0']` and overwrites `df.loc[index, 'Norm_Counts_' + str(i)]`.
- **Result:** Line 344 is rendered completely dead code. All linear regression fitting, slope calculations, and top-$N$ enrichment determinations performed on lines 320–341 are discarded. The final normalized counts in `PAMDA_2_norm_counts.csv` and rate fits in `PAMDA_3_rates.csv` receive **zero enrichment adjustment**.

---

## 4. Misaligned Barcode Start Offset Slicing During Strand Flip (`find_BCs_and_spacer`)

### Description
In `find_BCs_and_spacer`, when a spacer is detected in Read 2 instead of Read 1, the code swaps the barcode assignments but incorrectly applies `P5_sample_BC_start` to Read 1 and `P7_sample_BC_start` to Read 2.

### Code Location
`code/PAMDA.py` lines 582–589:
```python
else:
    spacer_loc = read_sequenceR2.find(spacers[sp])
    spacer_loc_rc = reverse_complement(read_sequenceR1).find(spacers[sp])
    if (spacer_loc != -1 or spacer_loc_rc != -1):
        spacer = sp
        P5_sample_BC = read_sequenceR1[P5_sample_BC_start:
                                       P5_sample_BC_start + P5_sample_BC_len]
        P7_sample_BC = read_sequenceR2[P7_sample_BC_start:
                                       P7_sample_BC_start + P7_sample_BC_len]
```

### Impact & Example Dataset Context
- **Asymmetric Read Structure:** In the sequencing strategy used for HT-PAMDA, sequencing is asymmetric: Read 1 is 65 bp long, while Read 2 is only 10 bp long (containing only the inline barcode).
- **Why it never triggers in the example dataset:**
  1. The target spacer sequences (e.g. `GGGCACGGGCAGCTTGCCGG`) are 20 nucleotides long. A 20-nt spacer **cannot physically exist** on Read 2 because Read 2 is only 10 nt long! Thus, `read_sequenceR2.find(...)` is never true, and this `else` block is **never entered** during the analysis of this dataset.
  2. Additionally, in `inputs.py`, `P5_SAMPLE_BARCODE_START = 2` and `P7_SAMPLE_BARCODE_START = 2`. Because both offsets are identical (`2`), slicing at offset 2 on R1 vs R2 yields the same offset regardless.
- **When this bug matters:** If HT-PAMDA is applied to standard symmetric paired-end datasets (e.g., 2x75bp or 2x150bp) where spacers can appear on Read 2, OR if P5 and P7 barcodes have different starting offsets (`P5_sample_BC_start != P7_sample_BC_start`), this code extracts sequence from the wrong physical positions, causing **100% demultiplexing failure (`wrong_barcode`)** for all reads with spacers on Read 2.

---

## Summary Matrix

| Error | Error Description | Impact on General Datasets | Impact on Provided Example Dataset |
| :--- | :--- | :--- | :--- |
| **1. P7 Barcode Length** | Uses P5 barcode length for P7 barcode slicing | 100% unmapped reads if P5/P7 barcode lengths differ | **No impact** (both P5 and P7 barcodes happen to be 4 nt) |
| **2. Timepoint Indexing** | Reads `Norm_Counts_0` before creation when multiplexing control | `KeyError` crash or wrong timepoint used for $t_0$ baseline | **No impact in default mode**; **100% crash** if multiplexed control mode is used |
| **3. Enrichment Normalization** | Line 345 overwrites line 344 using uncorrected `row` value | **Nullifies linear fits and top-$N$ enrichment correction** | **Direct mathematical impact** (enrichment corrections are discarded) |
| **4. Strand Flip Offsets** | Swaps physical read offsets during R2 spacer match | Demultiplexing failure for R2 spacer reads if offsets differ | **No impact** (spacer cannot fit in 10-bp R2; offsets both equal 2) |
