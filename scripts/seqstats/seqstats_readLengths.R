#!/usr/bin/Rscript
# seqstats_readLengths.R
# Read length plots for sequencing stats.
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 28 Nov 2016
# Updated on 28 Nov 2016

# Import necessary packages
# These may need to be installed first
if ("getopt" %in% rownames(installed.packages()) == F) {
    install.packages("getopt")
}
if ("ggplot2" %in% rownames(installed.packages()) == F) {
    install.packages("ggplot2")
}

suppressMessages(require("getopt"))
suppressMessages(require("ggplot2"))

# Suppress warning messages
options(warn=-1)

#################################################################
# UTILITY FUNCTIONS
#################################################################
# Set theme
my.theme <-
    theme(axis.text=element_text(colour="black", size=12),
          axis.title=element_text(face="bold", size=15),
          axis.title.x=element_text(margin=margin(t=10, b=5)),
          axis.ticks=element_blank(),
          axis.line.x=element_line(colour="black"),
          axis.line.y=element_line(colour="black"),
          legend.key=element_rect(fill=NA),
          legend.text=element_text(size=12),
          plot.title=element_text(size=16, face="bold"),
          panel.background=element_blank(),
          panel.margin=unit(3, "mm"),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank())


#################################################################
# ARGUMENT PARSING
#################################################################
spec <- matrix(c(
    "datafile",   "i", 1, "character",    "Data file (required)",
    "out_dir",    "d", 1, "character",    "Output directory (required)",
    "header",     "e", 0, "logical",      "Data file contains header info",
    "title",      "t", 1, "character",    "Title for plot",
    "help",       "h", 0, "logical",      "This help message"
    ), ncol=5, byrow=T)

opt <- getopt(spec)

# Check if help flag was given
if (!is.null(opt$help)) {
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
}
# Check for data file
if (is.null(opt$datafile)) {
    cat("\nDatafile not supplied. Use the '-i' option.\n\n")
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
} else {
    fp <- opt$datafile
}
# Check for output directory
if (is.null(opt$out_dir)) {
    cat("\nOutput directory not supplied. Use the '-d' option.\n\n")
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
} else {
    outdir <- opt$out_dir
}
# Check if header flag was given
if (is.null(opt$header)) {
    headerFlag <- F
} else {
    headerFlag <- opt$header
}
# Check title
if (is.null(opt$title)) {
    title <- ""
    file_title <- ""
} else {
    title <- opt$title
    file_title <- paste(title, "_", sep="")
}

#################################################################
# DATA PROCESSING
#################################################################
# Load data
data <- read.delim(fp, header=headerFlag)

# Replace column header if given
if (headerFlag) {
    colnames(data) <- c("V2")
}
data$V1 <- "1"
metricName <- "Read length (bp)"

# Plot density
mean.x <- mean(data$V2)
plot.breaks <- seq(min(data$V2), max(data$V2), length=5)
pl <- ggplot(data, aes(V2)) +
    geom_histogram(aes(y=..density..), binwidth=0.5, fill=NA, colour="#a7a7a7") +
    geom_density(adjust=0.25, alpha=0.3) +
    geom_vline(aes(xintercept=mean.x), colour="#d62728") +
    my.theme +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(breaks=c(plot.breaks, mean.x), labels=c(plot.breaks, round(mean.x, digits=1))) +
    xlab(metricName) + ylab("Density") +
    ggtitle(title)

ggsave(paste(outdir, "/", file_title ,"readlength_density.png", sep=""),
       plot=pl,
       width=20,
       height=15,
       units="cm",
       dpi=300)

# Box plots
pl <- ggplot(data, aes(x=V1, y=V2)) +
    geom_boxplot() +
    my.theme +
    scale_y_continuous(breaks=c(10, 20, 30, 40)) +
    xlab("") + ylab(metricName) +
    ggtitle(title)

ggsave(paste(outdir, "/", file_title ,"readlength_boxplot.png", sep=""),
       plot=pl,
       width=30,
       height=15,
       units="cm",
       dpi=300)
