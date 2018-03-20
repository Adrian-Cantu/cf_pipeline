# Sequence statistics plots

_Python_ and _R_ scripts to create bar charts and tables of sequencing statistics output.

### DEPENDENCIES
Required _Python2.7_ packages.

1. Numpy

Required _R_ packages.

1. getopt
2. ggplot2

Each _R_ package can be installed from CRAN by running the command:
```
install.packages(c("getopt", "ggplot2"))
```

### SAMPLE PLOTS
![qualbox](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/seqstats/sample/Sample_qualities_boxplots.png "Quality: Box Plot")

![qualdensity](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/seqstats/sample/Sample_qualities_density.png "Quality: Density Plot")

![gcbox](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/seqstats/sample/Sample_gcratios_boxplots.png "GC: Box Plot")

![gcdensity](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/seqstats/sample/Sample_gcratios_density.png "GC: Density Plot")

### USAGE
```
seqstats.sh version 0.1

usage: seqstats.sh -f fastq -d output_dir [Options]

Required
   -f [fastq]              : FASTQ file
   -o [output_dir]         : Directory for output files

Optional
   --gc                    : Flag to calculate GC
   --gz                    : Flag for gzipped compressed files
   --fasta                 : Flag for FASTA file instead of FASTQ
   -t                      : Title for plots
   -v                      : Verbose output
   -h, -?, --help          : This help message

Notes
   - None
```
