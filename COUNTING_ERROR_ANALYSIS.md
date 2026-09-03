# HT-PAMDA Counting & Demultiplexing Error Analysis

This document provides a detailed technical analysis of four critical, undisputable errors identified in the HT-PAMDA pipeline codebase (`code/PAMDA.py`). These errors cause crashes, incorrect read demultiplexing, indexing offsets, and silent math/normalization failures.

---

## 1. P7 Barcode Length Variable Assignment Bug (`fastq2count`)

### Description
In `code/PAMDA.py`, function `fastq2count`, line 203 sets `P7_sample_BC_len` using `P5_sample_BCs[0]` instead of `P7_sample_BCs[0]`.

### Code Location
`code/PAMDA.py` lines 201–203:
```python
P5_sample_BCs = variant_ids['P5_sample_barcode'].tolist()
P5_sample_BC_len = len(P5_sample_BCs[0])
P7_sample_BCs = variant_ids['P7_sample_barcode'].tolist()
P7_sample_BC_len = len(P5_sample_BCs[0])  # BUG: Should be P7_sample_BCs[0]
```

### Impact
When P5 and P7 barcodes have different lengths (for instance, a 4-nucleotide P5 barcode and a 6-nucleotide P7 barcode), `P7_sample_BC_len` is incorrectly set to 4 instead of 6.
When `find_BCs_and_spacer` extracts `P7_sample_BC`, it slices only 4 bases instead of 6. Consequently, no read ever matches the 6-base P7 barcode in `variant_dict`, leading to **100% unmapped reads (`wrong_barcode`)** and complete failure to demultiplex samples.

---

## 2. Timepoint Index Mismatch & `KeyError` on Control Timepoint Index 0 (`rawcount2normcount`)

### Description
There is a fundamental column indexing mismatch between how raw counts are stored and how `rawcount2normcount` attempts to access control sample counts when `CONTROL_RAW_COUNT_CSV` is `None` and `CONTROL_SAMPLE_TIMEPOINT_FASTQ` is provided.

### Code Location
In `fastq2count` (lines 213, 290):
`store_all_data` stores count lists indexed by 0-based timepoints (`0, 1, ..., N-1`). These are exported to `PAMDA_1_raw_counts.csv.gz` under column headers:
`Raw_Counts_1`, `Raw_Counts_2`, ..., `Raw_Counts_N`.

In `rawcount2normcount` (lines 306–314):
```python
for index, row in df.iterrows():
    if row['Sample'] == control_sample:
        control_dict[row['Spacer']][row['PAM']] = row['Norm_Counts_' +
                                                      str(control_sample_timepoint_fastq)]
```

### Impact
If a user sets `CONTROL_RAW_COUNT_CSV = None` and specifies `CONTROL_SAMPLE_TIMEPOINT_FASTQ = 0` (indicating that the control sample is multiplexed in the first timepoint fastq), the loop looks up `row['Norm_Counts_0']`.
However, `Norm_Counts_0` has **not been created yet** (it is created on line 313: `df['Norm_Counts_0'] = norm_counts_0`). The loop immediately crashes with:
`KeyError: 'Norm_Counts_0'`

Additionally, because 0-based timepoint fastq indices map to 1-based `Raw_Counts` columns (`timepoint = 0` -> `Raw_Counts_1`), indexing across control and experimental timepoints is shifted by 1.

---

## 3. Overwritten Top-N Enrichment Normalization (`rawcount2normcount`)

### Description
In `rawcount2normcount`, the top-N enrichment correction factor (`uptrend_corrections`) is calculated across timepoints, but the code overwrites the result in the subsequent loop without applying the correction factor.

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

### Impact
`row` is an immutable pandas `Series` snapshot returned by `df.iterrows()`. Updating `df.loc[index, ...]` on line 344 does NOT update `row`.
On line 345, `row['Norm_Counts_' + str(i)] / row['Norm_Counts_0']` uses the **original, uncorrected** `row['Norm_Counts_' + str(i)]` value and overwrites `df.loc[index, 'Norm_Counts_' + str(i)]`.
As a result, line 344 has no effect whatsoever. The top-N enrichment normalization is **completely discarded**, and uncorrected values are used for rate constant fitting.

---

## 4. Misaligned Barcode Start Offset Slicing During Strand Flip (`find_BCs_and_spacer`)

### Description
In `find_BCs_and_spacer`, when the spacer sequence is found in Read 2 rather than Read 1, the code swaps the barcode assignments but incorrectly applies `P5_sample_BC_start` to Read 1 and `P7_sample_BC_start` to Read 2.

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

### Impact
Physically in the library architecture:
- Read 1 always contains the P7 barcode at offset `P7_sample_BC_start`.
- Read 2 always contains the P5 barcode at offset `P5_sample_BC_start`.

When the spacer is in Read 2, the above code slices Read 1 at offset `P5_sample_BC_start` expecting a P5 barcode, and Read 2 at offset `P7_sample_BC_start` expecting a P7 barcode.
If `P5_sample_BC_start != P7_sample_BC_start`, this extracts sequence from the wrong physical positions of the reads, causing **demultiplexing failures (`wrong_barcode`)** for all reads where the spacer is oriented on Read 2.
