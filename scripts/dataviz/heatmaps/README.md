# Python Seaborn heatmaps

_Python Seaborn_ script to create heatmaps for viewing expression levels.

### DEPENDENCIES
Required _Python3.+_ packages.

1. Matplotlib
2. Numpy
3. Pandas
4. Seaborn


### SAMPLE PLOTS
![AR](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/dataviz/heatmaps/sample_output/AR.png "Antibiotic Resistance Plot")

![VF](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/dataviz/heatmaps/sample_output/VF.png "Virulence Factors Plot")

![VH](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/dataviz/heatmaps/sample_output/VH.png "Viral Hits Plot")


### USAGE
```
usage: draw_heatmap.py [-h] [-f {png,svg,eps}] [-d DPI] [--delim DELIM]
                       [-g NUMGENES] [--outname OUTNAME] [-t TITLE]
                       file outdir

Draw heatmap from a tabular file.

positional arguments:
  file                  Tabular input file
  outdir                Directory for output files

optional arguments:
  -h, --help            show this help message and exit
  -f {png,svg,eps}, --fmt {png,svg,eps}
                        Image format (default: png)
  -d DPI, --dpi DPI     Set image DPI (default: 100)
  --delim DELIM         Set input file delimiter (default: )
  -g NUMGENES, --numgenes NUMGENES
                        Number of genes to plot (default: 50)
  --outname OUTNAME     Custom ouput file name (default: None)
  -t TITLE, --title TITLE
                        Plot title (default: )
```
