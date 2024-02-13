# Landscape and Selection of Vaccine Epitopes in SARS-CoV-2

* This repository was initially hosted at https://github.com/Benjamin-Vincent-Lab/Landscape-and-Selection-of-Vaccine-Epitopes-in-SARS-CoV-2 and was moved here on 2/13/2024

Preprint submitted to biorxiv: https://www.biorxiv.org/content/10.1101/2020.06.04.135004v1

Several larger files which were used part of these analyses are included alongside supplemental data on Mendeley due to GitHub file size limits: https://data.mendeley.com/datasets/c6pdfrwxgj/3.

## Dependent packages
R packages:
   caret 6.0-84, cowplot 0.9.4, data.table 1.12.8, DESeq2 1.22.2, doMC 1.3.6, dplyr 0.8.4, forcats 0.4.0, GenomicRanges 1.34.0, ggallin 0.1.1, ggbeeswarm 0.6.0, ggnewscale 0.4.1, ggplot2 3.3.0, ggpubr 0.2, ggrepel 0.8.1, gplots 3.0.3, gridExtra 2.3, huxtable 4.7.1, magrittr 1.5, officer 0.3.10, pROC 1.16.2, RColorBrewer 1.1-2, readxl 1.3.1, scales 1.1.0, seqinr 3.6-1, stringr 1.4.0, venneuler 1.1-0, viridis 0.5.1 

## Instructions

All code used to analyze data and generate figures from this manuscript have been released within this github page, with some larger files (>100mb) stored on Mendeley data (https://data.mendeley.com/datasets/c6pdfrwxgj/3). In order to replicate analysis presented in this manuscript, follow the steps listed below.  R v3.5.2 and Python v3.7.5 were used for these analyses:

1. Clone this repository onto your local computer
2. Download the large data files (along with supplemental tables) from the above Mendeley data link.
3. Unpack the data files from above.  The subdirectories contained in the "Large_files" directory correspond to the same order as those within the github repository.  Merge these subdirectories together with your local repository files.
4. Run the following command in the repository root directory: 
   >tar -xvf ./SARS-CoV-2_epitope_landscape/Working/large_files_040521.tar.gz
5. For Figures 2-3, Fig. S1-S6, and Table S1-S8: Open the R file "Main_figures_resubmission.R" in the repository root directory, changing the "WORKING_ROOT" variable (line 2) to the path of the repository.  This file contains a step-by-step workflow for recreating these above listed figures.
6. Run `pip install -r requirements.txt`
7. To generate Figures 4-6, run `collect-and-generate-python-figures.sh`, which runs a series of IPython notebooks to select B-cell epitope regions, T-cell epitopes, vaccine peptides, and the generation of all related plots and figures. This script has a long time-out setting because it can take up to an hour.  
8. To copy all figure images to the `figure-images` directory, run `collect-all-figures.sh`. This requires [pdftk](https://stackoverflow.com/questions/20804441/how-to-install-pdftk-on-mac-os-x/22054037) to be installed. 
