# draw_heatmap_general.py
# Draw heatmaps from tabular formats using Seaborn.
# Originally developed for CF Pipeline data.
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 28 Mar 2017
# Updated on 28 Mar 2017


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import seaborn as sns
import pandas as pd
import argparse
import os
import sys


###############################################################################
## ARGUMENT PARSING
###############################################################################
parser = argparse.ArgumentParser(description="Draw heatmap from a "
                                 "tabular file.",
                                 formatter_class=
                                 argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("file", help="Tabular input file")
parser.add_argument("outdir", help="Directory for output files")
parser.add_argument("-f", "--fmt",
                    choices=["png", "svg", "eps"], default="png",
                    help="Image format")
parser.add_argument("-d", "--dpi", type=int, default=100,
                    help="Set image DPI")
parser.add_argument("--delim", default="\t",
                    help="Set input file delimiter")
parser.add_argument("-g", "--numgenes", type=int, default=50,
                    help="Number of genes to plot")
parser.add_argument("--outname", help="Custom ouput file name")
parser.add_argument("-t", "--title", default="", help="Plot title")

args = parser.parse_args()

input_file = args.file
outdir = args.outdir if args.outdir.endswith("/") else args.outdir + "/"

# Check file exists
if not os.path.isfile(input_file):
    util.printStatus("Intput file '" + input_file + "' does not exist.")
    util.exitScript()

# Check output directory exists
if not os.path.isdir(outdir):
    util.printStatus("Output directory '" + input_file + "' does not exist.")
    util.exitScript()

outname = args.outname if args.outname else \
    os.path.splitext(os.path.basename(input_file))[0]

###############################################################################
## BEGIN PROCESSING
###############################################################################
# Read in data file as a Panda DataFrame
# Data is assumed to be in a pivot tabular format where
# - First line holds column headers
# - Column headers indicate: gene | <columns of samples> | annotation
# - Each row corresponds to a gene
# - Cells indicate: gene ID | <expression levels> | annotation label
data = pd.read_csv(input_file, sep=args.delim, header=0)

# Rename "-" column to "gene"
data.rename(columns={"-": "gene"}, inplace=True)

# Set gene column as the index -- required for sum across rows
data.set_index(["gene"], inplace=True)

# Calculate sums of each gene to use for sort
# numeric_only used since annotation column needs to be ignored
sum_order = data.sum(axis=1, numeric_only=True)
sum_order.sort_values(ascending=False, inplace=True)

# Reindex genes based on new sorted index
data = data.reindex_axis(sum_order.index)

# Set new index to annotation name instead of genes
# This is for heatmap left axis labels
# Do not do this step if you still want to use gene IDs
data.set_index(["annotation"], inplace=True)

###############################################################################
## DRAW HEATMAP
###############################################################################
sns.set()
# Create colormap for heatmap
# Blue to black palette
cmap = sns.cubehelix_palette(n_colors=50, start=2.8, rot=0,
                             light=0.9, dark=0, as_cmap=True)

# Create plot
# Calculate size based on number of genes and number of samples
# [3 samples = 1 unit, 50 genes = 15 units]
nhoriz = int(np.ceil(len(data.columns) / 3))
nverti = int(np.ceil(15 / 50 * args.numgenes))
fig = plt.figure(figsize=[nhoriz, nverti])

# Use GridSpec for sizing of the color bar so it does not span the height of
# the heatmap.
# This will create a 3 x 2 grid. Blank elements will be made for
# [row 0, col 1] and [row 2, col 1].
# GridSpec arguments:
# width_ratios: heatmap column will be 10x size of color bar column
gs = gridspec.GridSpec(nrows=3, ncols=2, width_ratios=[10, 1], wspace=0.15)
hm_ax = plt.subplot(gs[:, 0])
cbar_ax = plt.subplot(gs[1, 1])

# Blank axes
blank_ax1 = plt.subplot(gs[0, 1])
blank_ax2 = plt.subplot(gs[2, 1])
blank_ax1.axis("off")
blank_ax2.axis("off")

# Plot heatmap
sns.heatmap(data.head(n=args.numgenes),
            linewidths=0.5, cmap=cmap, ax=hm_ax, cbar_ax=cbar_ax)

# Turn off left axis label of heatmap
hm_ax.set_ylabel("")

# Rotate axis tick labels
hm_ax.set_xticklabels(hm_ax.get_xticklabels(), rotation=90)
hm_ax.set_yticklabels(hm_ax.get_yticklabels(), rotation=0)

# Titles
hm_ax.set_title(args.title)
cbar_ax.set_title("Normalized\nexpression")

plt.savefig(os.path.join(outdir, outname + "." + args.fmt),
            dpi=args.dpi,
            bbox_inches="tight")
