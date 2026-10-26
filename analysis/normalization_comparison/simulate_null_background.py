"""
Estimate the background rate that the HT-PAMDA normalization produces from
sampling noise alone, and compare it with the observed background.

Simulated data use the real control-library PAM composition and the real
per-timepoint read depths of the WT SpCas9 sample. Only NGG PAMs are cleaved
(at the rates estimated by the corrected code); every other PAM has k = 0.
The simulated counts are run through both the original and the corrected
rawcount2normcount / normcount2rate.

Usage (from the repository root, inside the htpamdaenv environment):
    python analysis/normalization_comparison/simulate_null_background.py \
        <original_worktree> <corrected_worktree> <work_dir>

Each worktree must already contain the outputs of library_QC.py and
PAMDA_complete.py. The summary is written to <work_dir>/summary.txt and the
95th percentile of the null background (log10 k) to <work_dir>/detection_limit.txt.
"""
import importlib.util
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')

TIMEPOINTS = [0, 60, 480, 1920]
SPACERS = {'SPACER1': 'GGGCACGGGCAGCTTGCCGG', 'SPACER2': 'GTCGCCCTCGAACTTCACCT'}
SAMPLE = 'sample_01'  # WT SpCas9
N_REPLICATES = 5
SEED = 1
RESULT_DIR = 'output/example_PAMDA_data/PAM_start_0_length_4'


def load_pamda(worktree, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(worktree, 'code', 'PAMDA.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def is_ngg(pam):
    return pam[1:3] == 'GG'


def background(df):
    return df[(df['Sample'] == SAMPLE) & ~df['PAM'].map(is_ngg)]


def spacer_spearman(df):
    wide = background(df).pivot(index='PAM', columns='Spacer', values='Rate_Constant_k')
    return spearmanr(wide['SPACER1'], wide['SPACER2']).correlation


def simulate_counts(control, observed_counts, observed_rates, rng):
    sample_rows, control_rows = [], []
    for spacer in SPACERS:
        ctl = control[control['Spacer'] == spacer].set_index('PAM')['Raw_Counts_1']
        p0 = ctl / ctl.sum()
        pams = p0.index
        control_rows += [dict(Sample='control_sample', Spacer=spacer, PAM=pam, Raw_Counts_1=n)
                         for pam, n in zip(pams, rng.multinomial(ctl.sum(), p0.values))]
        obs = observed_counts[(observed_counts['Sample'] == SAMPLE) &
                              (observed_counts['Spacer'] == spacer)]
        k_obs = observed_rates[(observed_rates['Sample'] == SAMPLE) &
                               (observed_rates['Spacer'] == spacer)].set_index('PAM')['Rate_Constant_k']
        k_true = np.where([is_ngg(p) for p in pams], k_obs.reindex(pams).values, 0.0)
        counts = {}
        for i, t in enumerate(TIMEPOINTS[1:], 1):
            abundance = p0.values * np.exp(-k_true * t)
            counts[i] = rng.multinomial(int(obs['Raw_Counts_%d' % i].sum()), abundance / abundance.sum())
        for j, pam in enumerate(pams):
            sample_rows.append(dict(Sample=SAMPLE, Spacer=spacer, PAM=pam,
                                    **{'Raw_Counts_%d' % i: counts[i][j] for i in (1, 2, 3)}))
    return pd.DataFrame(sample_rows), pd.DataFrame(control_rows)


def main(original_worktree, corrected_worktree, work_dir):
    os.makedirs(work_dir, exist_ok=True)
    versions = {'original': load_pamda(original_worktree, 'pamda_original'),
                'corrected': load_pamda(corrected_worktree, 'pamda_corrected')}

    observed_counts = pd.read_csv(os.path.join(corrected_worktree, RESULT_DIR, 'PAMDA_2_norm_counts.csv'))
    observed_rates = pd.read_csv(os.path.join(corrected_worktree, RESULT_DIR, 'PAMDA_3_rates.csv'))
    control = pd.read_csv(os.path.join(original_worktree, 'output/example_PAM_library/PAMDA_1_raw_counts.csv.gz'))
    control['PAM'] = control['PAM'].str[:4]
    control = control.groupby(['Sample', 'Spacer', 'PAM'], as_index=False)['Raw_Counts_1'].sum()

    os.chdir(work_dir)
    rng = np.random.default_rng(SEED)
    summary = {name: [] for name in versions}
    null_log_k = {name: [] for name in versions}
    for _ in range(N_REPLICATES):
        sample_df, control_df = simulate_counts(control, observed_counts, observed_rates, rng)
        sample_df.to_csv('sim_raw_counts.csv', index=False)
        control_df.to_csv('sim_control_counts.csv', index=False)
        for name, pamda in versions.items():
            run = 'null_' + name
            pamda.rawcount2normcount(run, 'sim_control_counts.csv', 'control_sample', None, 'three_prime',
                                     4, 0, SPACERS, TIMEPOINTS, 4, 5, input_csv='sim_raw_counts.csv')
            pamda.normcount2rate(run, 4, 0, TIMEPOINTS, [0.0001, 0.001, 0.01], 4, 1, None)
            rates = pd.read_csv('output/%s/PAM_start_0_length_4/PAMDA_3_rates.csv' % run)
            bg = background(rates)['Rate_Constant_k']
            log_k = np.log10(bg.dropna())
            null_log_k[name].extend(log_k.tolist())
            summary[name].append((log_k.median(),
                                  (bg == rates['Rate_Constant_k'].min()).mean(),
                                  spacer_spearman(rates)))

    obs_bg = np.log10(background(observed_rates)['Rate_Constant_k'])
    detection_limit = np.percentile(null_log_k['corrected'], 95)
    above = background(observed_rates)[(obs_bg > detection_limit).values]

    lines = ['WT SpCas9, non-NGG PAMs (n = %d per run)' % len(obs_bg),
             'Observed, corrected code: median log10 k %.2f; spacer Spearman %.2f'
             % (obs_bg.median(), spacer_spearman(observed_rates))]
    for name in versions:
        s = np.array(summary[name])
        lines.append('Null simulation, %s code (%d replicates): median log10 k %.2f; fraction at floor %.2f; '
                     'spacer Spearman %.2f' % (name, N_REPLICATES, s[:, 0].mean(), s[:, 1].mean(), s[:, 2].mean()))
    lines.append('Detection limit (95th percentile of corrected-code null background): log10 k = %.2f'
                 % detection_limit)
    lines.append('Observed PAMs above the detection limit: %d of %d (about %.0f expected by chance)'
                 % (len(above), len(obs_bg), 0.05 * len(obs_bg)))
    lines.append('PAMs above the detection limit: ' + ', '.join(sorted(set(above['PAM']))))

    print('\n'.join(lines))
    with open('summary.txt', 'w') as f:
        f.write('\n'.join(lines) + '\n')
    with open('detection_limit.txt', 'w') as f:
        f.write('%.1f\n' % detection_limit)


if __name__ == '__main__':
    main(*[os.path.abspath(p) for p in sys.argv[1:4]])
