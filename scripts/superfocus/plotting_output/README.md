# SUPER-FOCUS plots

_R_ scripts to create bar charts and tables of SUPER-FOCUS output.

### DEPENDENCIES
Required _R_ packages.

1. getopt
2. ggplot2
3. reshape2
4. plyr
5. gridExtra

Each package can be installed from CRAN by running the command:
```
install.packages(c("getopt", "ggplot2", "reshape2", "plyr", "gridExtra"))
```

### SAMPLE PLOTS
![alldata](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/superfocus/sample/all_top_functions.png "All Top Functions")

![samplebar](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/superfocus/sample/sample1.fasta_top_functions.png "Sample Top Functions Barchart")

![sampletab](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/superfocus/sample/sample1.fasta_top_functions_table.png "Sample Top Functions Table")

### USAGE
```
usage: superfocus_plots.sh -d SF_dir [Options]

Required
   -d [SF_dir]             : SUPER-FOCUS directory of files

Optional
   --skip [NUM]            : Number of non-blank lines to skip before columns [Default: 2]
   --vir                   : Create virulence-specific plots
   -h, -?, --help          : This help message
   -v                      : Verbose output

Notes
   - This program specifically uses the SUPER-FOCUS output file:
      *_____results__all_levels_and_function.xls

   - R scripts specifically look for the columns:
      1) Subsystem Level 3 (case specific)
      2) SEED Function (case specific)
      3) Columns with the word 'abundance' (not case specific)
```
