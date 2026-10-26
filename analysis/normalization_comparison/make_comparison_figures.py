"""
Assemble one four-panel comparison image per sample:
  A. published figure (figures/github_example_figures)
  B. original code, regenerated
  C. corrected code, colour floor -5.0 (as in code/inputs.py)
  D. corrected code, colour floor at the detection limit

Usage:
    python analysis/normalization_comparison/make_comparison_figures.py \
        <repo_root> <png_dir> <output_dir> <detection_limit>

<png_dir> must contain original/, corrected/ and corrected_lod/ folders with
PAMDA_HEATMAP_<sample>.png files rendered from the pipeline's PDF output.
"""
import os
import sys

import matplotlib
matplotlib.use('agg')
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

SAMPLES = {'sample_01_WT_SpCas9': 'WT SpCas9',
           'sample_02_SpG': 'SpG',
           'sample_03_SpRY': 'SpRY'}


def main(repo, png_dir, out_dir, detection_limit):
    os.makedirs(out_dir, exist_ok=True)
    for sample, label in SAMPLES.items():
        name = 'PAMDA_HEATMAP_%s.png' % sample
        panels = [('A. Published figure (figures/github_example_figures)',
                   os.path.join(repo, 'figures', 'github_example_figures', name)),
                  ('B. Original code, regenerated (master)', os.path.join(png_dir, 'original', name)),
                  ('C. Corrected code (proposed-fixes), colour floor -5.0', os.path.join(png_dir, 'corrected', name)),
                  ('D. Corrected code, colour floor at detection limit (%s)' % detection_limit,
                   os.path.join(png_dir, 'corrected_lod', name))]
        fig, axes = plt.subplots(2, 2, figsize=(16, 13))
        for ax, (title, path) in zip(axes.ravel(), panels):
            image = mpimg.imread(path)
            h, w = image.shape[:2]
            ax.imshow(image[int(h * 0.08):int(h * 0.85), :int(w * 0.92)])
            ax.set_title(title, fontsize=13, loc='left')
            ax.axis('off')
        fig.suptitle('%s: effect of the corrections on the PAM heatmap (log10 rate)' % label, fontsize=16)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(os.path.join(out_dir, 'comparison_%s.png' % sample), dpi=100)
        plt.close(fig)


if __name__ == '__main__':
    main(*sys.argv[1:5])
