# scripts
Various scripts

## Content

### vcf_assign_maf_bins.py

Assigns MAF bin to each record from input VCF/BCF file. Requires `python 3` and `pysam` library. Type `python vcf_assign_maf_bins.py --help` to get the description of all parameters. You can also pipe input/ouput files from/into standard input/output by specifying `-` instead of file names.

### vcf_assign_af_bins.py

Same as **vcf_assign_maf_bins.py**, but assigns AF bins.

### vcf_extract_lof.py

Extracts all Loss-of-Function (LoF) variants from VCF/BCF file annotated using VEP. Requires `python 3` and `pysam` library. Type `python vcf_extract_lof.py --help` to get the description of all parameters.

### vcf_add_cadd_scores.py

Annotates VCF file with CADD scores from https://cadd.gs.washington.edu/download. Note: you must download also the corresponding index files. Type `python vcf_add_cadd_scores.py --help` to get the description of all parameters.

### vcf_assign_depth_bins.py

Assigns AD (allele depth) and DP (total depth) bin to each record from input **single sample** VCF/BCF file. Requires `python 3` and `pysam` library. Type `python vcf_assign_depth_bins.py --help` to get the description of all parameters. You can also pipe input/ouput files from/into standard input/output by specifying `-` instead of file names.
