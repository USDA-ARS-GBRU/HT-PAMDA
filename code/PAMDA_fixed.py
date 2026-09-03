import glob
import gzip
import itertools
import os
import sys
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.stats import skew
from tqdm.autonotebook import tqdm

# PAMDA_fixed.py - Contains proposed fixes for bugs identified in original PAMDA.py

def check_inputs(RUN_NAME,
                 BARCODE_CSV,
                 FASTQ_DIR,
                 TIMEPOINT_FASTQ,
                 PAM_ORIENTATION,
                 PAM_LENGTH,
                 PAM_START,
                 CONTROL_RAW_COUNT_CSV,
                 CONTROL_SAMPLE,
                 CONTROL_SAMPLE_TIMEPOINT_FASTQ,
                 TIMEPOINTS,
                 MAX_PAM_LENGTH,
                 SPACERS,
                 P5_SAMPLE_BARCODE_START,
                 P7_SAMPLE_BARCODE_START,
                 USE_TIMEPOINTS,
                 TOP_N_NORMALIZE,
                 INIT_RATE_EST,
                 READ_SUM_MIN,
                 TPS_SUM_MIN,
                 PAM1_NT_RANK,
                 PAM2_NT_RANK,
                 PAM1_INDEX_RANK,
                 PAM2_INDEX_RANK,
                 AVERAGE_SPACER,
                 HEATMAP_FIXED_MIN,
                 HEATMAP_FIXED_MAX,
                 LOG_SCALE_HEATMAP):
    """
    perform some checks for input parameters
    """
    if not os.path.exists(BARCODE_CSV):
        raise Exception('BARCODE_CSV "%s" not found' % BARCODE_CSV)
    if not os.path.exists(FASTQ_DIR):
        raise Exception('fastq directory "%s" not found' % FASTQ_DIR)
    if CONTROL_RAW_COUNT_CSV != None:
        if not os.path.exists(CONTROL_RAW_COUNT_CSV):
            raise Exception('CONTROL_RAW_COUNT_CSV "%s" not found' % CONTROL_RAW_COUNT_CSV)

    fastqs = glob.glob(FASTQ_DIR + '/**/*R1*.fastq.gz', recursive=True)
    if len(fastqs) == 0:
        raise Exception('no fastq files found')
    for fastq in fastqs:
        fastqR1 = fastq
        fastqR2 = fastq.replace('R1', 'R2')
        fastq_name = fastqR1.split('/')[-1]
        fastq_name = fastq_name.split('_L00')[0]
        try:
            TIMEPOINT_FASTQ[fastq_name]
        except:
            warnings.warn('%s not found in TIMEPOINT_FASTQ. This fastq will be ignored.' % fastq_name)

    if not isinstance(MAX_PAM_LENGTH, int):
        raise Exception('MAX_PAM_LENGTH should be an integer value, you entered: %s' % MAX_PAM_LENGTH)
    if not isinstance(P5_SAMPLE_BARCODE_START, int):
        raise Exception('P5_SAMPLE_BARCODE_START should be an integer value, you entered: %s' % P5_SAMPLE_BARCODE_START)
    if not isinstance(P7_SAMPLE_BARCODE_START, int):
        raise Exception('P7_SAMPLE_BARCODE_START should be an integer value, you entered: %s' % P7_SAMPLE_BARCODE_START)
    if not isinstance(PAM_LENGTH, int):
        raise Exception('PAM_LENGTH should be an integer value, you entered: %s' % PAM_LENGTH)
    if not isinstance(PAM_START, int):
        raise Exception('PAM_START should be an integer value, you entered: %s' % PAM_START)
    if (CONTROL_RAW_COUNT_CSV == None) and (CONTROL_SAMPLE_TIMEPOINT_FASTQ == None):
        raise Exception('Either CONTROL_RAW_COUNT_CSV or CONTROL_SAMPLE_TIMEPOINT_FASTQ must be specified')
    if not (isinstance(CONTROL_SAMPLE_TIMEPOINT_FASTQ, int) or (CONTROL_SAMPLE_TIMEPOINT_FASTQ == None)):
        raise Exception(
            'CONTROL_SAMPLE_TIMEPOINT_FASTQ should be "None" or an integer value, you entered: %s' % CONTROL_SAMPLE_TIMEPOINT_FASTQ)
    if not isinstance(TOP_N_NORMALIZE, int):
        raise Exception('TOP_N_NORMALIZE should be an integer value, you entered: %s' % TOP_N_NORMALIZE)
    if not isinstance(READ_SUM_MIN, int):
        raise Exception('READ_SUM_MIN should be an integer value, you entered: %s' % READ_SUM_MIN)
    if not isinstance(TPS_SUM_MIN, int):
        raise Exception('TPS_SUM_MIN should be an integer value, you entered: %s' % TPS_SUM_MIN)
    if not isinstance(TIMEPOINTS, list):
        raise Exception('TIMEPOINTS should be a list, you entered: %s' % TIMEPOINTS)
    if not isinstance(INIT_RATE_EST, list):
        raise Exception('INIT_RATE_EST should be a list, you entered: %s' % INIT_RATE_EST)
    if not isinstance(AVERAGE_SPACER, bool):
        raise Exception("AVERAGE_SPACER should be 'True' or 'False', you entered: %s" % AVERAGE_SPACER)
    if not isinstance(LOG_SCALE_HEATMAP, bool):
        raise Exception("LOG_SCALE_HEATMAP should be 'True' or 'False', you entered: %s" % LOG_SCALE_HEATMAP)
    if (HEATMAP_FIXED_MIN != False and not (
            isinstance(HEATMAP_FIXED_MIN, int) or isinstance(HEATMAP_FIXED_MIN, float))):
        raise Exception(
            "HEATMAP_FIXED_MIN should be 'False' or a float or an integer, you entered: %s" % HEATMAP_FIXED_MIN)
    if (HEATMAP_FIXED_MAX != False and not (
            isinstance(HEATMAP_FIXED_MAX, int) or isinstance(HEATMAP_FIXED_MAX, float))):
        raise Exception(
            "HEATMAP_FIXED_MAX should be 'False' or a float or an integer, you entered: %s" % HEATMAP_FIXED_MAX)
    if PAM_ORIENTATION not in ['three_prime', 'five_prime']:
        raise Exception("please enter 'three_prime' or 'five_prime' for PAM_ORIENTATION")
    if PAM_LENGTH > 8:
        warnings.warn('PAM lengths longer than 8 are not recommended')
    if PAM_LENGTH < 1:
        raise Exception('please choose a PAM length >0')
    if PAM_START + PAM_LENGTH > MAX_PAM_LENGTH:
        raise Exception('PAM_START (%s) + PAM_LENGTH (%s) is greater than MAX_PAM_LENGTH (%s)'
                        % (PAM_START, PAM_LENGTH, MAX_PAM_LENGTH))
    if (PAM1_INDEX_RANK != None and PAM2_INDEX_RANK != None):
        if not isinstance(PAM1_INDEX_RANK, list):
            raise Exception('PAM1_INDEX_RANK should be a list, you entered: %s' % PAM1_INDEX_RANK)
        if not isinstance(PAM2_INDEX_RANK, list):
            raise Exception('PAM2_INDEX_RANK should be a list, you entered: %s' % PAM2_INDEX_RANK)
        if len(PAM1_INDEX_RANK) + len(PAM2_INDEX_RANK) != PAM_LENGTH:
            raise Exception(
                'The number of ranked PAM positions in PAM1_INDEX_RANK and PAM2_INDEX_RANK is not equal to PAM_LENGTH')
    elif ((PAM1_INDEX_RANK == None and PAM2_INDEX_RANK != None) or
          (PAM1_INDEX_RANK != None and PAM2_INDEX_RANK == None)):
        raise Exception('Please specify both PAM1_INDEX_RANK and PAM2_INDEX_RANK or leave both as "None".')

#-----------------------------------------------------------------------------------------------------------------------------#

def fastq2count(run_name,
                barcode_csv,
                fastq_dir,
                timepoint_fastq,
                pam_orientation,
                timepoints=[0, 60, 480, 1920],
                max_pam_len=8,
                spacers={'SPACER1': 'GGGCACGGGCAGCTTGCCGG', 'SPACER2': 'GTCGCCCTCGAACTTCACCT'},
                P5_sample_BC_start=2,
                P7_sample_BC_start=2):
    """
    generate raw PAM read counts from fastq files
    """
    try:
        variant_ids = pd.read_csv(barcode_csv)
    except:
        raise Exception('BARCODE_CSV "%s" not found' % barcode_csv)

    fastqs = glob.glob(fastq_dir + '/**/*R1*.fastq.gz', recursive=True)
    if len(fastqs) == 0:
        raise Exception('no fastq files found')

    if pam_orientation not in ['three_prime', 'five_prime']:
        raise Exception("please enter 'three_prime' or 'five_prime' for PAM_ORIENTATION")

    P5_sample_BCs = variant_ids['P5_sample_barcode'].tolist()
    P5_sample_BC_len = len(P5_sample_BCs[0])
    P7_sample_BCs = variant_ids['P7_sample_barcode'].tolist()
    # FIX BUG 1: Set P7 barcode length from P7_sample_BCs instead of P5_sample_BCs
    P7_sample_BC_len = len(P7_sample_BCs[0])

    nt_complement = dict({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '_': '_', '-': '-'})

    nucleotides = ['A', 'T', 'C', 'G']
    total_pam_space = [''.join(p) for p in itertools.product(nucleotides, repeat=max_pam_len)]

    variant_dict = {}
    for index, row in variant_ids.iterrows():
        variant_dict[str(row['P5_sample_barcode']) + '_' + str(row['P7_sample_barcode'])] = row['sample']

    store_all_data = {}

    for sample in variant_ids['sample']:
        store_all_data[sample] = {spacer: {x: [0] * (len(timepoints) - 1) for x in total_pam_space}
                                  for spacer in spacers}

    pbar1 = tqdm(desc='fastq files: ', total=len(fastqs))
    pbar2 = tqdm(desc='reads: ')

    for fastq in fastqs:

        fastqR1 = fastq
        fastqR2 = fastq.replace('R1', 'R2')
        fastq_name = fastqR1.split('/')[-1]
        fastq_name = fastq_name.split('_L00')[0]

        try:
            timepoint = timepoint_fastq[fastq_name]
        except:
            continue

        if fastqR1.endswith('.gz'):
            infileR1 = gzip.open(fastqR1, 'rt')
            infileR2 = gzip.open(fastqR2, 'rt')
        else:
            infileR1 = open(fastqR1, 'r')
            infileR2 = open(fastqR2, 'r')

        wrong_barcode = 0
        wrong_spacer = 0
        total_reads = 0
        counted_reads = 0

        while infileR1.readline() and infileR2.readline():
            read_sequenceR1 = infileR1.readline().strip()
            infileR1.readline()
            infileR1.readline()
            read_sequenceR2 = infileR2.readline().strip()
            infileR2.readline()
            infileR2.readline()

            total_reads += 1

            top_read, bot_read, spacer, spacer_loc, P5_sample_BC, P7_sample_BC = \
                find_BCs_and_spacer(spacers, read_sequenceR1, read_sequenceR2,
                                    P5_sample_BC_start, P5_sample_BC_len,
                                    P7_sample_BC_start, P7_sample_BC_len)

            if spacer_loc == -1:
                wrong_spacer += 1
                continue

            if P5_sample_BC in P5_sample_BCs and P7_sample_BC in P7_sample_BCs:
                barcode_pair = P5_sample_BC + '_' + P7_sample_BC
                if barcode_pair in variant_dict.keys():
                    if pam_orientation == 'three_prime':
                        spacer3p = spacer_loc + len(spacers[spacer])
                        PAM = top_read[spacer3p: spacer3p + max_pam_len]
                        try:
                            store_all_data[variant_dict[barcode_pair]][spacer][PAM][timepoint] += 1
                            counted_reads += 1
                        except:
                            pass
                    elif pam_orientation == 'five_prime':
                        PAM = top_read[spacer_loc - max_pam_len: spacer_loc]
                        try:
                            store_all_data[variant_dict[barcode_pair]][spacer][PAM][timepoint] += 1
                            counted_reads += 1
                        except:
                            pass
            else:
                wrong_barcode += 1

            pbar2.update()

        pbar2.reset()
        pbar1.update()

        write_out = str(round(float(counted_reads) / float(total_reads) * 100, 2)) \
                    + '% of reads mapped from ' + str(fastq_name) + ' (' + str(counted_reads) + ' reads)'
        tqdm.write(write_out, file=sys.stdout)

    pbar1.close()
    pbar2.close()

    if not os.path.exists('output/%s' % run_name):
        os.makedirs('output/%s' % run_name)

    with gzip.open('output/%s/PAMDA_1_raw_counts.csv.gz' % (run_name), mode='wb') as f_out:
        f_out.write((','.join(map(str, ['Sample', 'Spacer', 'PAM'] +
                                  ['Raw_Counts_' + str(x)
                                   for x in range(1, len(timepoints))])) + '\n').encode('utf-8'))
        for fastq in store_all_data:
            for spacer in store_all_data[fastq]:
                for pam in store_all_data[fastq][spacer]:
                    total_info = [fastq, spacer, pam] + store_all_data[fastq][spacer][pam]
                    f_out.write((','.join(map(str, total_info)) + '\n').encode('utf-8'))

    raw_count_summary(run_name)

#-----------------------------------------------------------------------------------------------------------------------------#

def rawcount2normcount(run_name,
                       control_rawcount_csv,
                       control_sample,
                       control_sample_timepoint_fastq,
                       pam_orientation,
                       pam_length,
                       pam_start,
                       spacers={'SPACER1': 'GGGCACGGGCAGCTTGCCGG', 'SPACER2': 'GTCGCCCTCGAACTTCACCT'},
                       timepoints=[0, 60, 480, 1920],
                       max_pam_length=8,
                       top_n=5,
                       input_csv=None):
    """
    generate normalized PAM read counts from raw counts (FIXED VERSION)
    """
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    nucleotides = ['A', 'T', 'C', 'G']
    total_pam_space = [''.join(p) for p in itertools.product(nucleotides, repeat=pam_length)]

    if input_csv is None:
        df_input = pd.read_csv('output/%s/PAMDA_1_raw_counts.csv.gz' % (run_name))
    else:
        df_input = pd.read_csv(input_csv)

    control_timepoint_col = 'Norm_Counts_1'
    if control_rawcount_csv is not None:
        df_control = pd.read_csv(control_rawcount_csv)
        df_input = pd.concat([df_input, df_control], sort=False)
        control_sample_timepoint_fastq = 1
    else:
        # FIX BUG 2: Map 0-based control_sample_timepoint_fastq to 1-based Raw_Counts column index
        control_timepoint_col = 'Norm_Counts_' + str(control_sample_timepoint_fastq + 1)

    column_sort = {'Sample': 'first', 'Spacer': 'first'}
    count_columns = df_input.columns.values[3:]
    new_columns = ['Sample', 'Spacer', 'PAM']
    for count_column in count_columns:
        column_sort[count_column] = 'sum'
        new_columns.append(count_column)
    if pam_orientation == 'three_prime':
        df_input['selected_PAM'] = df_input['PAM'].str[pam_start:pam_start + pam_length]
    else:
        df_input['selected_PAM'] = df_input['PAM'].str[max_pam_length - pam_length - pam_start:
                                                       max_pam_length - pam_start]
    pbar = tqdm(desc='sample:', total=df_input['Sample'].nunique())
    df_list = []
    for sample in df_input['Sample'].unique().tolist():
        for spacer in df_input['Spacer'].unique().tolist():
            temp_df = df_input[(df_input['Sample'] == sample) & (df_input['Spacer'] == spacer)] \
                .groupby(['selected_PAM'], as_index=False).agg(column_sort)
            df_list.append(temp_df)
        pbar.update()
    pbar.close()
    df = pd.concat(df_list)
    df = df.rename(columns={'selected_PAM': 'PAM'})
    df = df.loc[:, new_columns]
    df = df.reset_index(drop=True)

    for i in range(1, len(timepoints)):
        df['Norm_Counts_' + str(i)] = df['Raw_Counts_' + str(i)] / \
                                      df.groupby(['Sample', 'Spacer'])['Raw_Counts_' + str(i)].transform(sum)

    # determine t0 read counts for control sample
    control_dict = {}
    for spacer in spacers:
        control_dict[spacer] = {}
        for PAM in total_pam_space:
            control_dict[spacer][PAM] = 0

    for index, row in df.iterrows():
        if row['Sample'] == control_sample:
            # FIX BUG 2: Use correct control timepoint column
            control_dict[row['Spacer']][row['PAM']] = row[control_timepoint_col]

    norm_counts_0 = []
    for index, row in df.iterrows():
        norm_counts_0.append(control_dict[row['Spacer']][row['PAM']])
    df['Norm_Counts_0'] = norm_counts_0
    df = df[df['Sample'] != control_sample]

    print('determining most enriched PAMs per sample')
    pbar = tqdm(desc='samples: ', total=df['Sample'].nunique(), file=sys.stdout)
    sample_last = None
    uptrends = {}
    x = range(len(timepoints))
    for index, row in df.iterrows():
        sample_current = row['Sample']
        if sample_current != sample_last:
            pbar.update()
        sample_spacer = str(row['Sample']) + '_' + str(row['Spacer'])
        y = [row['Norm_Counts_' + str(i)] for i in range(len(timepoints))]
        slope = linregress(x, y)
        if sample_spacer in uptrends:
            uptrends[sample_spacer].append([slope[0], y])
        else:
            uptrends[sample_spacer] = [[slope[0], y]]
        sample_last = sample_current
    pbar.close()

    uptrend_corrections = {}
    for u in uptrends:
        uptrends[u] = sorted(uptrends[u])
        top_n_entries = [x[1] for x in uptrends[u][-top_n:]]
        top_n_entries_reformat = map(list, zip(*top_n_entries))
        top_n_entries_median = [np.median(x) for x in top_n_entries_reformat]
        uptrend_corrections[u] = top_n_entries_median

    print('normalizing each sample:')
    pbar = tqdm(desc='samples: ', total=df['Sample'].nunique(), file=sys.stdout)
    sample_last = None
    for index, row in df.iterrows():
        sample_current = row['Sample']
        if sample_current != sample_last:
            pbar.update()

        sample_spacer = str(row['Sample']) + '_' + str(row['Spacer'])
        for i in range(len(timepoints)):
            # FIX BUG 3: Calculate intermediate normalized value without overwriting in single step
            norm_val = row['Norm_Counts_' + str(i)] / uptrend_corrections[sample_spacer][i]
            df.loc[index, 'Norm_Counts_' + str(i)] = norm_val / row['Norm_Counts_0']
        sample_last = sample_current
    pbar.close()

    if not os.path.exists('output/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length)):
        os.makedirs('output/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length))

    df.to_csv('output/%s/PAM_start_%s_length_%s/PAMDA_2_norm_counts.csv' %
              (run_name, pam_start, pam_length), index=False)

#-----------------------------------------------------------------------------------------------------------------------------#

def reverse_complement(seq):
    nt_complement = dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_','-':'-'})
    return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

def find_BCs_and_spacer(spacers, read_sequenceR1, read_sequenceR2,
                        P5_sample_BC_start, P5_sample_BC_len,
                        P7_sample_BC_start, P7_sample_BC_len):
    """
    find the sample barcodes and the spacer, orientation, and location (FIXED VERSION)
    """
    spacer_loc = -1
    spacer_loc_rc = -1
    spacer = None
    top_read = None
    bot_read = None
    P5_sample_BC = None
    P7_sample_BC = None

    for sp in spacers:
        spacer_loc = read_sequenceR1.find(spacers[sp])
        spacer_loc_rc = reverse_complement(read_sequenceR2).find(spacers[sp])
        if (spacer_loc != -1 or spacer_loc_rc != -1):
            spacer = sp
            P5_sample_BC = read_sequenceR2[P5_sample_BC_start:
                                           P5_sample_BC_start + P5_sample_BC_len]
            P7_sample_BC = read_sequenceR1[P7_sample_BC_start:
                                           P7_sample_BC_start + P7_sample_BC_len]
            if spacer_loc == -1:
                spacer_loc = spacer_loc_rc
                top_read = reverse_complement(read_sequenceR2)
                bot_read = read_sequenceR1
                break
            elif spacer_loc_rc == -1:
                top_read = read_sequenceR1
                bot_read = reverse_complement(read_sequenceR2)
                break
            elif spacer_loc > spacer_loc_rc:
                spacer_loc = spacer_loc_rc
                top_read = reverse_complement(read_sequenceR2)
                bot_read = read_sequenceR1
                break
            else:
                top_read = read_sequenceR1
                bot_read = reverse_complement(read_sequenceR2)
                break
        else:
            spacer_loc = read_sequenceR2.find(spacers[sp])
            spacer_loc_rc = reverse_complement(read_sequenceR1).find(spacers[sp])
            if (spacer_loc != -1 or spacer_loc_rc != -1):
                spacer = sp
                # FIX BUG 4: P5 barcode is always on Read 2, P7 barcode is always on Read 1
                P5_sample_BC = read_sequenceR2[P5_sample_BC_start:
                                               P5_sample_BC_start + P5_sample_BC_len]
                P7_sample_BC = read_sequenceR1[P7_sample_BC_start:
                                               P7_sample_BC_start + P7_sample_BC_len]
                if spacer_loc == -1:
                    spacer_loc = spacer_loc_rc
                    top_read = reverse_complement(read_sequenceR1)
                    bot_read = read_sequenceR2
                    break
                elif spacer_loc_rc == -1:
                    top_read = read_sequenceR2
                    bot_read = reverse_complement(read_sequenceR1)
                    break
                elif spacer_loc > spacer_loc_rc:
                    spacer_loc = spacer_loc_rc
                    top_read = reverse_complement(read_sequenceR1)
                    bot_read = read_sequenceR2
                    break
                else:
                    top_read = read_sequenceR2
                    bot_read = reverse_complement(read_sequenceR1)
                    break

    return top_read, bot_read, spacer, spacer_loc, P5_sample_BC, P7_sample_BC

def raw_count_summary(run_name, input_csv=None):
    if input_csv is None:
        df_input = pd.read_csv('output/%s/PAMDA_1_raw_counts.csv.gz' % run_name)
    else:
        df_input = pd.read_csv(input_csv)

    df_output = df_input.groupby(['Sample', 'Spacer']).sum()

    if not os.path.exists('output/%s' % run_name):
        os.makedirs('output/%s/' % run_name)

    df_output.to_csv('output/%s/PAMDA_1_raw_counts_summary.csv' % run_name)
