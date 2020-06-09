import argparse
import pysam

#
# Note: python3 + pysam library are required
#
# Example:
#   python vcf_extract_lof.py -i input.vcf.gz -o output.txt
#

argparser = argparse.ArgumentParser(description = 'Extracts LoF variants from VCF/BCF annotated with VEP+LoFtee.')
argparser.add_argument('-i', '--in-vcf', metavar = 'file', dest = 'in_VCF', required = True, help = 'VEP+LoFtee annotated input VCF/BCF file.')
argparser.add_argument('-o', '--out', metavar = 'file', dest= 'out_filename', required = True, help = 'Output tab-delimited file.')


output_vep_fields = ['SYMBOL', 'Gene', 'Feature', 'Consequence', 'CLIN_SIG', 'PUBMED', 'gnomAD_AF', 'gnomAD_NFE_AF', 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info'] 


if __name__ == '__main__':
    args = argparser.parse_args()
    

    with pysam.VariantFile(args.in_VCF, 'r') as vcf_in, open(args.out_filename, 'wt') as file_out:
        # check if CSQ INFO field is present
        csq_meta = vcf_in.header.info.get('CSQ', None)
        if csq_meta is None:
            raise Exception('No meta-information entry about CSQ INFO field found!')

        # get VEP header
        csq_header = csq_meta.description.split(':', 1)[1].strip().split('|')

        # write output header
        file_out.write('CHROM\tPOS\tREF\tALT\tAF')
        for field in output_vep_fields:
            if field in csq_header:
                file_out.write(f'\t{field}')
        file_out.write('\n')

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
                allele_annotations = annotations[i]
 
                # iterate over all affected transcripts
                for transcript in allele_annotations:
                    if transcript['BIOTYPE'] != 'protein_coding': # we want only variants in/near protein coding regions
                        continue
                    if transcript['Feature_type'] != 'Transcript': # we are interested in transcripts only
                        continue
                    if transcript['LoF'] != '': # if this flag is not empty, then we have LoF variant
                        file_out.write(f'{record_in.chrom}\t{record_in.pos}\t{record_in.ref}\t{alt_allele}\t{record_in.info["AF"][i - 1]}')
                        for field in output_vep_fields:
                            if field in csq_header:
                                file_out.write(f'\t{"NA" if transcript[field] == "" else transcript[field]}')
                        file_out.write('\n')



