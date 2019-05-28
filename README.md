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
### Removing ribosomal RNA from .fastaq files
The **rRNA_remover.py** script will achieve this. In the terminal, simply run:  
> python3 /data2/utilities/RNA-Seq_utilities/rRNA_remover.py  

with the following arguments:


`-h`, `--help`  
show this help message and exit
  
`-d` &lt;DIRECTORY&gt;, `--directory` &lt;DIRECTORY&gt;  
Specify the location of the RNA-Seq data

`-l` &lt;FILE&gt;, `--rRNA_library` &lt;FILE&gt;  
Specify location of the rRNA genome library. i.e. path to the **.fa file**. Default is the _C. elegans_ library.

`-s`, `--single_end`  
Flag if RNA-Seq data are single end reads. Mutually exclusive with the `-p` / `--paired_end` argument. 
 
`-p` &lt;PAIR_TAG&gt; &lt;PAIR_TAG&gt;, `--paired_end` &lt;PAIR_TAG&gt; &lt;PAIR_TAG&gt;   
Flag if RNA-Seq data are paired end reads. Mutually exclusive with the `-s` / `--single_end` argument. Provide 
pair tags, this will be the same as PRAGUI's "`pair_tags`" argument (e.g. `r_1 r_2`).

#### Example
> **python3** /data2/utilities/RNA-Seq_utilities/rRNA_remover.py **-d** /scratch/gurpreet/data/ **-l** /scratch/ribosomal_rna/worm/c_elegans_concat_rDNA.fa **-p** r_1 r_2

### Merging RNA-Seq data into one file per lane
THe **rna_seq_lane_merger.py** script will achieve this. In the terminal, simply run:
> python3 /data2/utilities/RNA-Seq_utilities/rna_seq_lane_merger.py

with the following arguments:

`-f`, `--submission_form` &lt;FILE&gt;
Path to the submission form provided (e.g. CRUKCI_SLX_Submission.xlsx) - Please provide full path and ensure this file is in same folder as the RNA-Seq files.

`-l`, `--lane_tags` &lt;lane_tag&gt; &lt;lane_tag&gt;
Tags that identify samples' RNA-Seq lanes e.g. `s_1 s_2`.

`-s`, `--single_end`
Flag if RNA-Seq data are single end reads. Mutually exclusive with the `-p` / `--paired_end` argument.

`-p` &lt;PAIR_TAG&gt; &lt;PAIR_TAG&gt;, `--paired_end` &lt;PAIR_TAG&gt; &lt;PAIR_TAG&gt;
Flag if RNA-Seq data are paired end reads. Mutually exclusive with the `-s` / `--single_end` argument. Provide pair tags, this will be the same as PRAGUI's "pair_tags" argument.

#### Example
> **python3** /data2/utilities/RNA-Seq_utilities/rna_seq_lane_merger.py **-f** /scratch/gurpreet/rna_seq_data/CRUKCI_SLX_Submission.xlsx **-l** s_1 s_2 **-p** r_1 r_2

