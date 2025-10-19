#!/bin/bash

# Workflow for processing Illumina ITS2 data using AMPtk pipeline.
# Raw sequences are available from NCBI’s SRA under BioProject ID PRJNA1344880
# Install AMPtk: https://amptk.readthedocs.io/en/latest/



################# Set directory #################


# Create new directory 
mkdir -p oil_project

# Move your downloaded raw Illumina sequence archive into this directory
# Example: mv downloaded_sequences.zip oil_project/oil_raw_sequences.zip

cd oil_project

# Unzip the raw sequences
unzip oil_raw_sequences.zip
mv */ oil_raw_sequences



################# AMPtk PIPELINE #################

#### 1. Import sequences ####

amptk illumina -i oil_raw_sequences -o  oil -f AACTTTYRRCAAYGGATCWCT -r AGCCTCCGCTTATTGATATGCTTAART --require_primer off --rescue_forward off --primer_mismatch 6


#### 2. Clustering data ####

amptk dada2 -i oil.demux.fq.gz -o oil 


#### 3. Index bleeding #### 

amptk filter -i oil.cluster.otu_table.txt -f oil.cluster.otus.fa -p 0.005


#### 4. Assign taxonomy (using ITS database) ####

amptk install -i ITS # Installing ITS database
amptk taxonomy -f oil.cluster.filtered.otus.fa -i oil.cluster.final.txt -m oil.mapping_file.txt -d ITS2 -o oil


## 5. FUNguild: assign functional groups

amptk funguild -i oil.otu_table.taxonomy.txt -d fungi -o funguild_oil



################# OUTPUTS ################# 

# Taxonomy finished: oil.taxonomy.txt
# Classic OTU table with taxonomy: oil.otu_table.taxonomy.txt
# BIOM OTU table created: oil.biom
# OTUs with taxonomy: oil.otus.taxonomy.fa
# OTU phylogeny: oil.tree.phy
# FUNguild annotation: funguild_oil.txt
