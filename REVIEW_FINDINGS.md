# HT-PAMDA: Review Findings and Proposed Corrections

## Purpose and scope

This note describes four issues in the published HT-PAMDA analysis code (`code/PAMDA.py`) and one place where the documentation does not match the code. Each issue was confirmed by running the code. The proposed corrections are deliberately minimal: they change only the lines needed and leave the method and structure of the pipeline as published.

The corrections are on the `proposed-fixes` branch (commit `5472048`), which is based on upstream `master` at `f215a53`. All line numbers below refer to `code/PAMDA.py` at `f215a53`.

| # | Function | Summary | Affects the bundled example data? |
|---|---|---|---|
| 1 | `fastq2count` | P7 barcode length is taken from the P5 barcode list | No |
| 2 | `rawcount2normcount` | Top-N enrichment correction is computed but not applied | **Yes** |
| 3 | `normcount2rate` | Zero-count filtering is misaligned when `USE_TIMEPOINTS` omits timepoint 0 | No (default settings) |
| 4 | `normcount2rate` | Fallback value of 100 for non-positive rates when no positive rate exists | No |
| 5 | Documentation | `PAM1_INDEX_RANK` / `PAM2_INDEX_RANK` description does not match the code | No |

## Validation approach

- The pipeline was run end to end through rate fitting on the bundled example data, using the settings in `code/inputs.py` and `code/library_QC_inputs.py`. The control library is `fastqs/example_PAM_library`; the timepoints are `fastqs/example_PAMDA_data`; `PAM_LENGTH = 4`, giving 3 samples × 2 spacers × 256 PAMs = 1536 rate constants. Plotting libraries were stubbed out because plotting does not affect counts or rates.
- Issues that do not occur with the example settings were reproduced with small synthetic inputs passed directly to the affected functions.
- The original and corrected code were run on the same inputs and their outputs compared.

---

## Issue 1: P7 barcode length is taken from the P5 barcode list

**Location:** `fastq2count`, line 288

```python
P5_sample_BCs = variant_ids['P5_sample_barcode'].tolist()
P5_sample_BC_len = len(P5_sample_BCs[0])
P7_sample_BCs = variant_ids['P7_sample_barcode'].tolist()
P7_sample_BC_len = len(P5_sample_BCs[0])
```

**Effect:** `P7_sample_BC_len` is set to the length of the first **P5** barcode. If the P5 and P7 barcodes in `BARCODE_CSV` differ in length, every extracted P7 barcode has the wrong length and fails the `P7_sample_BC in P7_sample_BCs` test. The affected reads are counted as `wrong_barcode`. The bundled example data is unaffected because all of its P5 and P7 barcodes are 4 nt long.

**Validation:** 10 synthetic read pairs with a 4-nt P5 barcode (`ATGC`) and a 6-nt P7 barcode (`TCGCGA`):

| | Reads mapped |
|---|---|
| Original | 0.0% (0 of 10) |
| Corrected | 100.0% (10 of 10) |

**Correction:**
```diff
-    P7_sample_BC_len = len(P5_sample_BCs[0])
+    P7_sample_BC_len = len(P7_sample_BCs[0])
```

---

## Issue 2: The top-N enrichment correction is computed but not applied

**Location:** `rawcount2normcount`, lines 543–547

```python
for i in range(len(timepoints)):
    df.loc[index, 'Norm_Counts_' + str(i)] = row['Norm_Counts_' + str(i)] / \
                                             uptrend_corrections[sample_spacer][i]
    df.loc[index, 'Norm_Counts_' + str(i)] = row['Norm_Counts_' + str(i)] / \
                                             row['Norm_Counts_0']
```

**Effect:** The function computes a per-timepoint correction from the top-N most enriched PAMs (lines 497–526). This corrects for the rise in relative abundance of uncleaved PAMs as cleavable ones are depleted (`TOP_N_NORMALIZE`). `row` is the copy returned by `df.iterrows()`, so the second assignment recomputes from the **uncorrected** value and overwrites the first. The stored value is `Norm_Counts_i / Norm_Counts_0`, and the enrichment correction has no effect on the normalized counts or on the rate constants.

Without the correction, uncleaved PAMs appear to increase over time. In the example data, the median normalized abundance of non-NGG PAMs in the WT SpCas9 sample goes 1.00 → 1.04 → 1.06 → 1.11 across the four timepoints. The fit then gives a negative rate constant, and lines 646–648 replace every non-positive rate with the smallest positive rate in the run.

**Validation (example data, 1536 rate constants):**

| | Original | Corrected |
|---|---|---|
| Rates ≤ 0, replaced by the run minimum | 961 (63%) | 26 (1.7%) |
| Replaced, by sample (01 / 02 / 03) | 364 / 350 / 247 | 9 / 9 / 8 |
| WT SpCas9, non-NGG PAMs: median k | 1.8 × 10⁻⁷ (the replacement value) | 8.6 × 10⁻⁵ |
| WT SpCas9, NGG PAMs: median k | 0.0348 | 0.0358 |
| Spearman correlation, original vs corrected, per sample/spacer | | 0.80 – 0.95 |

The fastest rates (canonical NGG PAMs for WT SpCas9) barely change. The difference is concentrated in slow and intermediate PAMs, which make up most of the heatmap's range. In the original output, most of these PAMs share a single replacement value.

**Correction:** apply the correction to both the numerator and the t0 reference, so that the t0 value stays exactly 1.0 as before.

```diff
-            df.loc[index, 'Norm_Counts_' + str(i)] = row['Norm_Counts_' + str(i)] / \
-                                                     uptrend_corrections[sample_spacer][i]
-            df.loc[index, 'Norm_Counts_' + str(i)] = row['Norm_Counts_' + str(i)] / \
-                                                     row['Norm_Counts_0']
+            df.loc[index, 'Norm_Counts_' + str(i)] = (row['Norm_Counts_' + str(i)] /
+                                                      uptrend_corrections[sample_spacer][i]) / \
+                                                     (row['Norm_Counts_0'] /
+                                                      uptrend_corrections[sample_spacer][0])
```

Note on the design choice: dividing only by `row['Norm_Counts_0']` instead of the corrected t0 value gives the same rate constants (to within 0.03% on the example data), because the fitted amplitude `a` in `y = a·e^(−kt)` absorbs a per-sample constant. That variant would change the t0 column of `PAMDA_2_norm_counts.csv` to a value other than 1.0, so the form above was chosen to keep that output's meaning unchanged.

---

## Issue 3: Zero-count filtering is misaligned when `USE_TIMEPOINTS` omits timepoint 0

**Location:** `normcount2rate`, lines 604–613

```python
tps = [x[1] for x in timepoint_indices]
obs_raw = [0.00000000001] + [row['Raw_Counts_' + str(x[0])]
                             for x in timepoint_indices if str(x[0]) != '0']
obs_norm = [row['Norm_Counts_' + str(x[0])] for x in timepoint_indices]

zero_indices = [i for i, e in enumerate(obs_raw) if e == 0]
for zero_index in sorted(zero_indices)[::-1]:
    del tps[zero_index]
    del obs_norm[zero_index]
    del obs_raw[zero_index]
```

**Effect:** A placeholder for t0 is always prepended to `obs_raw`, even when timepoint 0 is not in `USE_TIMEPOINTS`. In that case `obs_raw` is one element longer than `tps` and `obs_norm`, and the zero-count indices are off by one:

- If a PAM has zero reads at the last selected timepoint, `del tps[zero_index]` raises `IndexError: list assignment index out of range`, and the run stops.
- If a PAM has zero reads at another timepoint, the following timepoint is removed instead. The zero-count observation stays in the fit and a valid observation is dropped, with no error.

The README documents subsets of timepoints as a supported option. With the default (`USE_TIMEPOINTS = None`), or any subset that includes 0, the lists line up and results are unaffected.

**Validation:** a synthetic three-PAM input generated with k = 0.002, with no zero counts, a zero count at t = 60, or a zero count at t = 1920:

| `USE_TIMEPOINTS` | Original | Corrected |
|---|---|---|
| `None` | 0.002, 0.002, 0.002 | 0.002, 0.002, 0.002 |
| `[0,1,2,3]` | 0.002, 0.002, 0.002 | 0.002, 0.002, 0.002 |
| `[1,2,3]` | `IndexError` | 0.002, 0.002, 0.002 |

**Correction:** build `obs_raw` in the same order as `timepoint_indices`, using the placeholder only where timepoint 0 is selected.

```diff
-        obs_raw = [0.00000000001] + [row['Raw_Counts_' + str(x[0])]
-                                     for x in timepoint_indices if str(x[0]) != '0']
+        obs_raw = [0.00000000001 if str(x[0]) == '0' else row['Raw_Counts_' + str(x[0])]
+                   for x in timepoint_indices]
```

---

## Issue 4: Fallback value of 100 for non-positive rates

**Location:** `normcount2rate`, line 646

```python
min_k = min([100 if ((isinstance(x, float) and x <= 0) or x == 'NaN') else x for x in ks])
df['Rate_Constant_k'] = [x if ((isinstance(x, float) and x > 0) or x == 'NaN') else min_k for x in ks]
```

**Effect:** Non-positive rates are meant to be replaced by the smallest positive rate in the run. When a run has **no** positive rate, the minimum is the sentinel value 100, so every non-positive rate is reported as 100. That is the fastest value in the output, which reverses the intended meaning. This case is unlikely with real data, but a run that has failed (for example, one with no detectable cleavage) would be reported as uniformly fast cleavage.

**Validation:** a synthetic input in which every PAM increases over time:

| | Reported rates |
|---|---|
| Original | 100, 100 |
| Corrected | NaN, NaN |

Whenever at least one positive rate exists, as in the example data, this correction does not change the output.

**Correction:** take the minimum over positive rates only, and fall back to `'NaN'` (the convention already used in this function) when there are none.

```diff
-    min_k = min([100 if ((isinstance(x, float) and x <= 0) or x == 'NaN') else x for x in ks])
+    positive_ks = [x for x in ks if isinstance(x, float) and x > 0]
+    min_k = min(positive_ks) if len(positive_ks) != 0 else 'NaN'
```

---

## Issue 5 (documentation): `PAM1_INDEX_RANK` and `PAM2_INDEX_RANK`

**Location:** `README.md` (parameter descriptions) and the matching docstring in `code/inputs.py`. The code is `rate2heatmap`, lines 702–747.

**Current documentation:** "The position of the list corresponds to the position in the first or second half of the PAM and the integer value at each position in the list corresponds to the priority of the position. A higher integer value indicates a higher priority." The example `[2,1]` is said to give the first position the higher priority.

**Behavior of the code:** the **order** of the list sets the priority, and the integer **values** identify positions. The position named by the first element varies most slowly along the axis. With `PAM1_INDEX_RANK = [2,1]`, the rows come out `AA, CA, GA, TA, AC, …`, so the **second** position has the higher priority, which is the reverse of the documented behavior. The built-in defaults follow the code's convention: for a 4-nt PAM they are `[1,0]` and `[2,3]`, which prioritize the two inner bases, as the code comment says.

This affects only the order of rows and columns. Each heatmap cell is still labeled with the correct PAM and holds the correct rate.

**Correction:** documentation only. The code is unchanged, so existing heatmaps are reproduced exactly. The description now says how the list is interpreted, and the examples use the defaults for a 4-nt PAM (`[1,0]` for the y-axis, `[2,3]` for the x-axis), with a description of the order each produces. The same text is applied to `README.md` and `code/inputs.py`.

---

## Summary of changes on `proposed-fixes`

```
 README.md      |  6 +++---
 code/PAMDA.py  | 17 +++++++++--------
 code/inputs.py |  6 +++---
```

Using the bundled example data and default settings, the only change in output comes from Issue 2: different normalized counts and rate constants, and therefore different heatmaps. Issues 1, 3 and 4 affect only configurations not used by the example, and Issue 5 changes documentation only.
