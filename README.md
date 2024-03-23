# NanoSort-Bonavita-Wells-Anthony-2024
This repository contains code for NanoSort published in Bonavita, Wells, and Anthony 2024: _Frequency of recombination between coronaviruses in vitro is determined by rates of cellular co-infection_.

## Requirements
The following programs are required to run NanoSort:
- BLAST (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
- seqkit (https://bioinf.shenwei.me/seqkit/)
- minimap2 (https://github.com/lh3/minimap2)
- samtools (https://www.htslib.org/)
- GNU Parallel (https://www.gnu.org/software/parallel/)
  
In addition, the following python packages should be installed:
- numpy
- pandas

## Installation
To use NanoSort, download the repository and add the location to your PATH environment variable: 

```
export PATH='/path/to/NanoSort/:$PATH'
```
Within the NanoSort folder, you additionally must edit the file `nanosort.sh` with the same location for the PATH_TO_NANOSORT variable (line 6): 
```
export PATH_TO_NANOSORT="/path/to/NanoSort"
```

## Usage
The program is called within the shell script `nanosort.sh` and requires 3 input items: 
- the first parental virus
- the second parental virus
- the fastq file

These items must be entered in the exact order. NanoSort currently supports three pairs of viruses:
| Parent 1      | Parent 2 |
| ----------- | ----------- |
| CaCoV_1-71      | FeCoV_WSU-79-1683       |
| CaCoV_UCD-1   | FIPV_WSU-79-1146        |
| RCoV_Parker  | MHV_JHM  |

A test file containing recombinant, subgenomic, and structural variant reads is included in the `test` folder. You can run the test with the following command:
```
nanosort.sh CaCoV_1-71 FeCoV_WSU-79-1683 test_reads.fastq.gz
```

## Output
NanoSort provides several output files. A single `overall_summary.tsv` file is generated containing a table of each viral read and a 1/0 column for whether the read was identified as recombinant, subgenomic, or structural variant. For each category, a list of read names assigned to that category is also generated. Finally, an output reads folder is created that contains a `.json` file for each read with information about the mapping position, likelihood score, and classification for each sliding sub-read.

The total number of reads in each category can be obtained by counting the number of reads in each read names file or by summing the 1/0 columns in `overall_summary.tsv`.
