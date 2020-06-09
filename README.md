# scripts
Various scripts

## Content

### vcf_assign_maf_bins.py

Assigns MAF bin to each record from input VCF/BCF file. Requires `python 3` and `pysam` library. Type `python vcf_assign_maf_bins.py --help` to get the description of all parameters. You can also pipe input/ouput files from/into standard input/output by specifying `-` instead of file names.

### vcf_extract_lof.py

Extracts all Loss-of-Function (LoF) variants from VCF/BCF file annotated using VEP. Requires `python 3` and `pysam` library. Type `python vcf_assign_maf_bins.py --help` to get the description of all parameters.
