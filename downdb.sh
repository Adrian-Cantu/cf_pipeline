#!/bin/bash
#Get databases for metagenomics pipeline

#UniVec:
wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/README.uv

#Human Genome 38: GRCh38.p7
mkdir -p HG38 
cd HG38
wget ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq/*.fa.gz
gunzip *
wget ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/README_CURRENT_RELEASE
cat *.fa > HG38.fasta 
rm hs_alt*
cd ..

#Greengenes
mkdir -p greengenes
cd greengenes
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5.fasta.gz
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5_with_header.fasta.gz
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/00README
cd ..

#RefSeq Viral
mdir -p viral
cd viral
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz
gunzip viral.1.1.genomic.fna.gz
gunzip viral.2.1.genomic.fna.gz
cat viral.1.1.genomic.fna viral.2.1.genomic.fna > viral_refseq.fasta
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/README
cd ..

#RefSeq Bacteria
mkdir -p refseq_bacteria
cd refseq_bacteria
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/*.fna.gz
gunzip *
cat *.fna > bacteria_refseq.fasta
cd ..

