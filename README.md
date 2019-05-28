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
The `rRNA_remover.py` script will achieve this. In the terminal, simply run:  
> python3 /data2/utilities/RNA-Seq_utilities/rRNA_remover.py  

with the following arguments:


`-h`, `--help`  
show this help message and exit
  
`-d <DIRECTORY>`, `--directory <DIRECTORY>`  
Specify the location of the RNA-Seq data

`-l <FILE>`, `--rRNA_library <FILE>`  
Specify location of the rRNA genome library. i.e. path to the .fa file. Default is the C. elegans library.

`-s`, `--single_end`  
Flag if RNA-Seq data are single end reads. Mutually exclusive with the `-p`/`--paired_end` argument. 
 
`-p <PAIR_TAG> <PAIR_TAG>`, `--paired_end <PAIR_TAG> <PAIR_TAG>`  
Flag if RNA-Seq data are paired end reads. Mutually exclusive with the `-s`/`--single_end` argument. Provide 
pair tags, this will be the same as PRAGUI's "`pair_tags`" argument (e.g. `r_1 r_2`).
