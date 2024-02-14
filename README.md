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

For each variant within trio (mother, father, offspring) identifies possible Mendelian error. Requires `python 3` and `pysam` library. Trios (or nuclear families) must be provided in PED file format (no header). Genotypes for all samples from trios (or nuclear families) must be in the same VCF/BCF file (multiple VCFs/BCFs split by chromosome are allowed).  Type `mendelian_per_trio.py --help` to get the description of all parameters.

### mendelian_per_trio_split.py

Same as **mendelian_per_trio.py**, but requires three sigle-sample input VCF/BCF files: for father, for mother, and for child. Type `mendelian_per_trio_split.py --help` to get the description of all parameters.

### vcf_gene_set.py
Create gene sets file from the VEP annotated VCF/BCF for gene-burden association tests (e.g. using SAIGE or EPACTS). Requires `python 3` and `pysam` library. Type `vcf_gene_set.py --help` to get the description of all parameters.

### variants_per_individual.py
Counts number of variants per individual from VCF/BCF file. Variants are grouped by most severe consequence based on VEP annotations. Requires `python 3` and `pysam` library. Type `variants_per_individual.py --help` to get the description of all parameters.

### plink2reference.py
Helps to align A1/A2 alleles in PLINK to the human genome reference. Keeps only SNPs with A,C,G, and T alleles. Removes palindromic, A/T and G/C, SNPs. Flips everything on the `+` strand. Type `plink2reference.py --help` to get the description of all parameters.

### plink_dupvar_prioritize.py
Helps to remove duplicated variants with the highest missigness from PLINK files. Type `plink_dupvar_prioritize.py --help` to get the description of all parameters.

### plink_freq_concordance.py
Compares genotype (allele) frequency between two PLINK files using Fisher's exact test. Type `plink_freq_concordance.py --help` to get the description of all parameters. You may need to load R (use `module spider r` and `module load r`) before running the script.

### plink_freq_strat_concordance.py
Similar to `plink_freq_concordance.py`, but allows stratification by the groups (e.g. ancestry). Type `plink_freq_strat_concordance.py --help` to get the description of all parameters. You may need to load R (use `module spider r` and `module load r`) before running the script.

### Mahalanobis.py
Performs a Mahalanobis distance of projected genomes. Type `Mahalanobis.py --help` to get the description of all parameters. Produces for each study sample: distance and its derived p-value per reference categories, and a list of non-rejected reference categories (at input threshold).

### merge_minimac_dosages.py
Merges imputed dosages when the genotype imputation was performed in subsets of samples. Subsets can have shared samples. Type `merge_minimac_dosages.py --help` to get the description of all parameters.

### imputation_dosage_comparison.py 
Produce sum of absolute differences of Dosages in overlap samples of two impute VCF file (per samples and per variant). Requires `python 3`, `pysam`, 'numpy' library.

### plink_freq_vs_topmed.py
Compares allele frequencies in PLINK file to TOPMed reference panel. Type `plink_freq_vs_topmed.py --help` to get the description of all parameters, and look inside the script for the additional hints.

### {GRAPH}_Manhattan_Miami.py
Plots a Miami or Manahttan plot (depending on options selected) form Regenie Output.

### {GRAPH}_Plot_PCA_Ancestry.py
Plots PCA Projections.

