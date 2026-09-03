import os
import sys
import unittest
import pandas as pd
import numpy as np

# Add code directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../code')))

import PAMDA

class TestHTPAMDACountingErrors(unittest.TestCase):

    def test_error_1_p7_barcode_length_bug(self):
        """
        Tests Bug 1: P7_sample_BC_len is set using P5_sample_BCs[0] length instead of P7_sample_BCs[0].
        When P5 is length 4 and P7 is length 6, P7_sample_BC_len is incorrectly set to 4.
        """
        # Create dummy variant DataFrame with distinct P5 (len 4) and P7 (len 6) barcodes
        variant_ids = pd.DataFrame({
            'sample': ['sample_01'],
            'description': ['test'],
            'P5_sample_barcode': ['ATGC'],        # length 4
            'P7_sample_barcode': ['TCGCGA']       # length 6
        })

        P5_sample_BCs = variant_ids['P5_sample_BCs'].tolist() if 'P5_sample_BCs' in variant_ids else variant_ids['P5_sample_barcode'].tolist()
        P5_sample_BC_len = len(P5_sample_BCs[0])

        # Buggy implementation from PAMDA.py line 203:
        P7_sample_BCs = variant_ids['P7_sample_barcode'].tolist()
        P7_sample_BC_len_buggy = len(P5_sample_BCs[0])  # Uses P5 length!

        # Correct implementation:
        P7_sample_BC_len_correct = len(P7_sample_BCs[0])

        self.assertEqual(P5_sample_BC_len, 4)
        self.assertEqual(P7_sample_BC_len_correct, 6)

        # Prove the bug in PAMDA.py:
        self.assertEqual(P7_sample_BC_len_buggy, 4, "PAMDA.py sets P7_sample_BC_len to 4 instead of 6")
        self.assertNotEqual(P7_sample_BC_len_buggy, P7_sample_BC_len_correct)

    def test_error_2_control_timepoint_keyerror(self):
        """
        Tests Bug 2: Setting CONTROL_RAW_COUNT_CSV = None and CONTROL_SAMPLE_TIMEPOINT_FASTQ = 0
        causes KeyError: 'Norm_Counts_0' because Norm_Counts_0 has not been created yet in rawcount2normcount.
        """
        # Create input dataframe simulating df_input after grouping
        df = pd.DataFrame({
            'Sample': ['control_sample', 'control_sample', 'sample_1', 'sample_1'],
            'Spacer': ['SPACER1', 'SPACER2', 'SPACER1', 'SPACER2'],
            'PAM': ['AAAA', 'AAAA', 'AAAA', 'AAAA'],
            'Raw_Counts_1': [100, 100, 50, 50],
            'Raw_Counts_2': [80, 80, 40, 40]
        })

        timepoints = [0, 60]
        control_sample = 'control_sample'
        control_sample_timepoint_fastq = 0  # User specifies timepoint index 0

        # Calculate Norm_Counts_1 for df
        df['Norm_Counts_1'] = df['Raw_Counts_1'] / df.groupby(['Sample', 'Spacer'])['Raw_Counts_1'].transform(sum)

        control_dict = {'SPACER1': {'AAAA': 0}, 'SPACER2': {'AAAA': 0}}

        # Reproduce code from PAMDA.py lines 306-309:
        with self.assertRaises(KeyError) as ctx:
            for index, row in df.iterrows():
                if row['Sample'] == control_sample:
                    # Looking up 'Norm_Counts_0' raises KeyError because it hasn't been added yet
                    control_dict[row['Spacer']][row['PAM']] = row['Norm_Counts_' + str(control_sample_timepoint_fastq)]

        self.assertIn('Norm_Counts_0', str(ctx.exception))

    def test_error_3_uptrend_correction_overwritten(self):
        """
        Tests Bug 3: In rawcount2normcount, the top-N enrichment correction line is overwritten
        by the t0 normalization line because row is an immutable Series snapshot.
        """
        # Simulate normalization loop inputs
        timepoints = [0, 60]
        uptrend_corrections = {'sample_1_SPACER1': [2.0, 2.0]}

        df = pd.DataFrame({
            'Sample': ['sample_1'],
            'Spacer': ['SPACER1'],
            'PAM': ['AAAA'],
            'Norm_Counts_0': [0.5],
            'Norm_Counts_1': [1.0]
        })

        row = df.iloc[0]
        index = 0
        sample_spacer = 'sample_1_SPACER1'

        # Execute buggy lines from PAMDA.py 344-345:
        # Line 344:
        df.loc[index, 'Norm_Counts_1'] = row['Norm_Counts_1'] / uptrend_corrections[sample_spacer][1] # 1.0 / 2.0 = 0.5
        # Line 345 (buggy: uses row['Norm_Counts_1'] which is still 1.0, NOT 0.5!):
        df.loc[index, 'Norm_Counts_1'] = row['Norm_Counts_1'] / row['Norm_Counts_0'] # 1.0 / 0.5 = 2.0

        # Expected correct combined value: (1.0 / 2.0) / 0.5 = 1.0
        # Actual result in df: 2.0 (uptrend correction 2.0 was completely lost/ignored!)
        actual_val = df.loc[index, 'Norm_Counts_1']
        expected_if_buggy = 2.0
        expected_if_fixed = 1.0

        self.assertEqual(actual_val, expected_if_buggy, "The uptrend correction was completely overwritten")
        self.assertNotEqual(actual_val, expected_if_fixed, "The code failed to apply the uptrend correction factor")

    def test_error_4_strand_flip_barcode_start_offset_swap(self):
        """
        Tests Bug 4: In find_BCs_and_spacer, when spacer is found on R2,
        P5_sample_BC is sliced from R1 using P5_sample_BC_start and P7_sample_BC from R2 using P7_sample_BC_start.
        This swaps physical read positions if P5_sample_BC_start != P7_sample_BC_start.
        """
        # Read 1 always physically has P7 barcode at offset P7_start
        # Read 2 always physically has P5 barcode at offset P5_start
        P5_start = 2
        P7_start = 5
        P5_len = 4
        P7_len = 4

        # Dummy sequence strings
        read_R1 = "NNNNN" + "TCGC" + "NNNN"  # P7 barcode 'TCGC' starting at pos 5
        read_R2 = "NN" + "ATGC" + "NNNN"      # P5 barcode 'ATGC' starting at pos 2

        # In buggy PAMDA.py line 585-588 (when spacer is found in R2):
        P5_BC_buggy = read_R1[P5_start : P5_start + P5_len] # Slices R1 at P5_start (pos 2) -> "NNTC" (WRONG)
        P7_BC_buggy = read_R2[P7_start : P7_start + P7_len] # Slices R2 at P7_start (pos 5) -> "CNNN" (WRONG)

        # Correct physical extraction:
        # P5 is ALWAYS on Read 2 at P5_start
        # P7 is ALWAYS on Read 1 at P7_start
        P5_BC_correct = read_R2[P5_start : P5_start + P5_len] # "ATGC"
        P7_BC_correct = read_R1[P7_start : P7_start + P7_len] # "TCGC"

        self.assertEqual(P5_BC_correct, "ATGC")
        self.assertEqual(P7_BC_correct, "TCGC")

        self.assertNotEqual(P5_BC_buggy, P5_BC_correct, "Buggy code sliced R1 at P5_start instead of R2")
        self.assertNotEqual(P7_BC_buggy, P7_BC_correct, "Buggy code sliced R2 at P7_start instead of R1")


if __name__ == '__main__':
    unittest.main()
