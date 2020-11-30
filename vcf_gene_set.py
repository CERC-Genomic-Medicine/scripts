import argparse
from collections import OrderedDict
import pysam

#
# Note: python3 + pysam library are required
#
# Example:
#   python vcf_gene_set.py -i input.vcf.gz -o output.txt
#

argparser = argparse.ArgumentParser(description = 'Constructs gene sets from VCF/BCF annotated with VEP+LoFtee for gene burden association tests.')
argparser.add_argument('-i', '--in-vcf', metavar = 'file', dest = 'in_VCF', required = True, help = 'VEP+LoFtee annotated input VCF/BCF file.')
argparser.add_argument('-c', '--max-ac', metavar = 'number', dest = 'max_ac', required = True, type = int, help = 'Maximal alternate allele count (AC).')
argparser.add_argument('-l', '--lof-only', dest = 'lof_only', action = 'store_true', help = 'Include only LoF variants.')
argparser.add_argument('-o', '--out', metavar = 'file', dest= 'out_filename', required = True, help = 'Output tab-delimited file.')

cds_variant_types = [
    'start_lost',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'missense_variant',
    'inframe_insertion',
    'inframe_deletion',
    'stop_gained',
    'stop_lost',
    'frameshift_variant',
    'splice_donor_variant',
    'splice_acceptor_variant'
]

if __name__ == '__main__':
    args = argparser.parse_args() 
    with pysam.VariantFile(args.in_VCF, 'r') as vcf_in:
        # check if CSQ INFO field is present
        csq_meta = vcf_in.header.info.get('CSQ', None)
        if csq_meta is None:
            raise Exception('No meta-information entry about CSQ INFO field found!')

        # get VEP header
        csq_header = csq_meta.description.split(':', 1)[1].strip().split('|')
        genes = OrderedDict()

        # read input VCF records
        for record_in in vcf_in:

            # get all VEP annotations for each alternate allele
            annotations = dict()
            for annotation in record_in.info['CSQ']:
                annotation = annotation.split('|')
                assert len(csq_header) == len(annotation)
                annotation = dict(zip(csq_header, annotation))
                annotations.setdefault(int(annotation['ALLELE_NUM']), []).append(annotation)

            # iterate over all alternate alleles (i.e. in case variants is multi-allelic)
            for i, alt_allele in enumerate(record_in.alts, 1):
                if not 'AC' in record_in.info:
                    print(f'Warning: variant {record_in.chrom}:{record_in.pos}_{record_in.ref}/{alt_allele} does\'t have AC INFO field.')
                    continue
                ac = record_in.info['AC'][i - 1]
                if ac > args.max_ac:
                    continue
                allele_annotations = annotations[i]
                # iterate over all affected transcripts
                for transcript in allele_annotations:
                    if transcript['BIOTYPE'] != 'protein_coding': # we want only variants in/near protein coding regions
                        continue
                    if transcript['Feature_type'] != 'Transcript': # we are interested in transcripts only
                        continue
                    consequences = transcript['Consequence'].split('&')
                    if not any(c in cds_variant_types for c in consequences):
                        continue
                    is_lof = transcript['LoF'] == 'HC'
                    if args.lof_only and not is_lof:
                        continue
                    variant_name = f'{record_in.chrom}:{record_in.pos}_{record_in.ref}/{alt_allele}'
                    gene_name = transcript['Gene']
                    gene_variants = genes.setdefault(gene_name, OrderedDict())
                    if variant_name not in gene_variants:
                        gene_variants[variant_name] = None
             
    with open(args.out_filename, 'w') as ofile:
        sep = '\t'
        for gene, variants in genes.items():
            ofile.write(f'{gene}{sep}{sep.join(variants.keys())}\n')
                    
