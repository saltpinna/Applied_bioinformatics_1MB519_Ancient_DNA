# Applied_bioinformatics_1MB519_Ancient_DNA
Repository for Ancient DNA project in Applied bioinformatics course (1MB519) at Uppsala university 2022.

## Required packages
- Numpy
- Pandas
- Biopython

## Scripts
Following are brief descriptions of all scripts in the repository

### Incrementalstatistics

### GC_analysis_Mammoth
Writes a file with the average comparative depth ratio for each GC% as well as the standard deviation and the number of sliding windows contributing to the GC%. Used on the Mammoth coverage data files.

### GC_analysis_Neanderthal
Writes a file with the average comparative depth ratio for each GC% as well as the standard deviation and the number of sliding windows contributing to the GC%. Used on the Neanderthal coverage data files.

### find_telomere
Finds the last chromosome position of the telomere in a FASTA sequence file, where telomeres are denoted as a sequence of 'N's in the beginning of the FASTA sequence

### gene_nongene_depth_comparison
Returns two files of sequence depth values for windows in genic regions and intergenic regions respectively, based on a reference genome FASTA file and gene position information in that reference

### incremental_statistics
Calculates the genome wide average coverage and the standard deviation, including or excluding outliers
