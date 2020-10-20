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


### vcf_lof_downsample.py

Subsets random sample from the specified VCF/BCF file and counts number of putative Loss-of-Function variants in the sample. Requires `python 3` and `pysam` library. Requires VCF/BCF file annotated using VEP+LOFTEE. Type `python vcf_lof_downsample.py --help` to get the description of all parameters.

### gencode_cds.py

Extracts CDS (coding sequence regions) for each protein coding gene from GENCODE GTF file. For each gene, it uses a union of CDS regions from all protein coding transcripts. Requires `python 3` and `intervaltree` library. Type `python gencode_cds.py --help` to get the description of all parameters.

### gencode_utr.py

Extracts 5' UTR and 3' UTR regions for each protein coding gene from GENCODE GTF file. For each gene, it uses a union of UTR regions from all protein coding transcripts. Requires `python 3` and `intervaltree` library. Type `python gencode_utr.py --help` to get the description of all parameters.

### bed_compute_gc.py

Computes number of bases (A, C, G, T) and GC-content for regions in BED file. Requires `python 3` and `pysam` library. Type `python bed_compute_gc.py --help` to get the description of all parameters.

### mendelian_per_trio.py

For each variant within trio (mother, father, offspring) identifies possible Mendelian error. Requires `python 3` and `pysam` library. Trios (or nuclear families) must be provided in PED file format (no header). Genotypes for all samples from trios (or nuclear families) must be in the VCF/BCF files split by chromosome.  Type `mendelian_per_trio.py --help` to get the description of all parameters.
