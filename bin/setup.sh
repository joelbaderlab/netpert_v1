#!/usr/bin/env bash

# Create conda environment
conda env create -f bin/netpert_env.yml

# Make database directory
mkdir databases

# Download and unzip large files
cd ./databases
curl -O http://iid.ophid.utoronto.ca/static/download/mouse_annotated_PPIs.txt.gz ; gunzip mouse_annotated_PPIs.txt.gz
curl -O http://iid.ophid.utoronto.ca/static/download/human_annotated_PPIs.txt.gz ; gunzip human_annotated_PPIs.txt.gz
curl -O https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv
curl -O https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv
curl -O https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz ; gunzip MOUSE_10090_idmapping.dat.gz
curl -O https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt
curl -O https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt
cd ..
