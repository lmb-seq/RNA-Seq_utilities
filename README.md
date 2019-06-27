# RNA-Seq utilities
Tools and utilities to supplement [PRAGUI](https://github.com/lmb-seq/PRAGUI).
Mainly pre-processing data preparation.

## Getting started

### Prerequisites
The [cell_bio_util](https://github.com/lmb-seq/cell_bio_util) repository will need 
to be installed next to the directory holding this project. In the case of MRC-LMB's 
Mario-Xeon machine, this has already been done.

For Mario-Xeon, RNA-Seq utilities is located in:  
> /data2/utilities/RNA-Seq_utilities/  
  
## Running the scripts  

-----------------------------------------------
## Removing ribosomal RNA from .fastaq files  
The **rRNA_remover.py** script will achieve this. In the terminal, simply run:  
> python3 /data2/utilities/RNA-Seq_utilities/rRNA_remover.py  

with the following arguments:

| Flag                                  	| Description                                                                                                                                                                                                           	|
|---------------------------------------	|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|
| `-h`, `--help`                        	| Show this help message and exit                                                                                                                                                                                       	|
| `-d`, `--directory` &lt;directory&gt; 	| Specify the location of the RNA-Seq data                                                                                                                                                                              	|
| `-l`, `--rRNA_library` &lt;file&gt;   	| Specify location of the rRNA genome library. i.e. path to the **.fa** file. Default is the _C. elegans_ library.                                                                                                      	|
| `-s`, `--single_end`                  	| Flag if RNA-Seq data are single end reads. Mutually exclusive with the `-p` / `--paired_end` argument.                                                                                                                	|
| `-p`, `--paired_end` &lt;pair_tag&gt; 	| Flag if RNA-Seq data are paired end reads. Mutually exclusive with the `-s` / `--single_end` argument. Provide  space separated pair tags, this will be the same as PRAGUI's "`pair_tags`" argument (e.g. `r_1 r_2`). 	|

#### Example  
> **python3** /data2/utilities/RNA-Seq_utilities/rRNA_remover.py **-d** /scratch/gurpreet/data/ **-l** /scratch/ribosomal_rna/worm/c_elegans_concat_rDNA.fa **-p** r_1 r_2  

-----------------------------------------------
## Merging RNA-Seq data into one file per sample, per lane
The **rna_seq_lane_merger.py** script will achieve this. In the terminal, simply run:
> python3 /data2/utilities/RNA-Seq_utilities/rna_seq_lane_merger.py

with the following arguments:

| Flag                                   	| Description                                                                                                                                                                                       	|
|----------------------------------------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|
| `-h`, `--help`                         	| Show this help message and exit                                                                                                                                                                   	|
| `-f`, `--submission_form` &lt;file&gt; 	| Path to the submission form provided (e.g. CRUKCI_SLX_Submission.xlsx) - Please ensure this file is in same folder as the RNA-Seq files.                                                          	|
| `-l`, `--lane_tags` &lt;lane_tag&gt;   	| Tags (space separated) that identify samples' RNA-Seq lanes e.g. `s_1 s_2`.                                                                                                                       	|
| `-s`, `--single_end`                   	| Flag if RNA-Seq data are single end reads. Mutually exclusive with the `-p` / `--paired_end` argument.                                                                                            	|
| `-p`, `--paired_end` &lt;pair_tag&gt;  	| Flag if RNA-Seq data are paired end reads. Mutually exclusive with the `-s` / `--single_end` argument. Provide space separated pair tags, this will be the same as PRAGUI's "pair_tags" argument. 	|

#### Example  
> **python3** /data2/utilities/RNA-Seq_utilities/rna_seq_lane_merger.py **-f** /scratch/gurpreet/rna_seq_data/CRUKCI_SLX_Submission.xlsx **-l** s_1 s_2 **-p** r_1 r_2

-----------------------------------------------
## Calculating mean and standard deviation of the TPM values
The **tpm_standard_deviation_mean_calculator.py** script will achieve this. In the terminal, simply run:
> python3 /data2/utilities/RNA-Seq_utilities/tpm_standard_deviation_mean_calculator.py

with the following arguments:

| Flag                            	|  Description                                                                                                                                    	|
|---------------------------------	|-------------------------------------------------------------------------------------------------------------------------------------------------	|
| `-t`, `--tpm_file` &lt;file&gt; 	| Full path for the tpm.txt file of interest (often uses "samples.csv_tpm.txt" filename - where "samples.csv" refers to PRAGUI's input csv file). 	|

#### Example
> **python3** /data2/utilities/RNA-Seq_utilities/tpm_standard_deviation_mean_calculator.py **-t** /scratch/gurpreet/rna_seq_data/samples.csv_tpm.txt

-----------------------------------------------

## Downloading your files from the CRUK FTP server
The **cruk_downloader.py** script will achieve this. In the terminal, simply run:
> python3 /data2/utilities/RNA-Seq_utilities/cruk_downloader.py

with the following arguments:

| Flag                            	|  Description                                                                                                                                    	|
|---------------------------------	|-------------------------------------------------------------------------------------------------------------------------------------------------	|
| `-f`, `--submission_form` &lt;file&gt; 	| Full path to the submission form (e.g. CRUKCI_SLX_Submission.xlsx) - Please ensure this file is in same folder as where you wish to download the RNA-Seq files to. 	|

This script reads in the CRUKCI_SLX_Submission.xlsx form and automatically retrieves the SLX ID and list of your files with which it will download to a directory of your choosing.

#### Example
> **python3** /data2/utilities/RNA-Seq_utilities/cruk_downloader.py **-f** /scratch/gurpreet/rna_seq_data/CRUKCI_SLX_Submission.xlsx

-----------------------------------------------