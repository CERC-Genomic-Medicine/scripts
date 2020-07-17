import argparse
import pysam
import random

argparser = argparse.ArgumentParser(description = 'Downsample VCF and count LoF variants.')
argparser.add_argument('-l', '--lof', metavar = 'file', dest = 'in_LoF_VCFs', required = True, nargs='+', help = 'VEP+LoFtee annotated input VCF/BCF file (genotypes are not required).')
argparser.add_argument('-g', '--gt', metavar = 'file', dest = 'in_GT_VCFs', required = True, nargs='+', help = 'VCF/BCF file with genotypes.')
argparser.add_argument('-k', '--k-random-individuals', metavar = 'number', dest = 'k_random', required = True, type = int, help = 'Size of random sample.')
argparser.add_argument('-i', '--iterations', metavar = 'number', dest = 'n_iterations', required = True, type = int, help = 'Number of sampling iterations.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_filename', required = True, help = 'Output file name.')


def read_lof(in_VCFs):
    for in_VCF in in_VCFs:
        with pysam.VariantFile(in_VCF, 'r') as vcf_in:
            # check if CSQ INFO field is present
            csq_meta = vcf_in.header.info.get('CSQ', None)
            if csq_meta is None:
                raise Exception('No meta-information entry about CSQ INFO field found!')
            # get VEP header
            csq_header = csq_meta.description.split(':', 1)[1].strip().split('|')

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
                            yield f'{record_in.chrom}-{record_in.pos}-{record_in.ref}-{alt_allele}'


def get_sample_list(in_VCFs):
    samples = list()
    for in_VCF in in_VCFs:
        with pysam.VariantFile(in_VCF, 'r') as vcf_in:
            if not samples:
                samples = list(vcf_in.header.samples)
            else:
                if not all(s1 == s2 for s1, s2 in zip(samples, list(vcf_in.header.samples))):
                    raise Exception('Samples in VCFs don\'t match')
    return samples


def sample_vcf(in_VCFs, samples, variants):
    n_variants_sampled = 0
    for in_VCF in in_VCFs:
        with pysam.VariantFile(in_VCF, 'r') as vcf_in:
            vcf_in.subset_samples(samples)
            for record_in in vcf_in:
                for i, alt_allele in enumerate(record_in.alts, 1):
                    variant_name = f'{record_in.chrom}-{record_in.pos}-{record_in.ref}-{alt_allele}'
                    if not variant_name in variants:
                        continue
                    for sample_format in record_in.samples.itervalues():
                        if i in sample_format['GT']:
                            n_variants_sampled += 1
                            break
    return n_variants_sampled



if __name__ == '__main__':
    args = argparser.parse_args()

    lof_variants = set()

    samples = get_sample_list(args.in_GT_VCFs)
    print(f'Total {len(samples)} samples.')

    for variant_name in read_lof(args.in_LoF_VCFs):
        lof_variants.add(variant_name)
    print(f'Total {len(lof_variants)} LoF variants.')

    random.seed()

    with open(args.out_filename, 'wt') as out:
        out.write(f'N_INDIVIDUALS\tN_LOF\n')
        for i in range(0, args.n_iterations):
            print(f'Iteration {i + 1} with {args.k_random} sampled individuals')
            random_samples = random.sample(samples, args.k_random)
            n_sampled_variants = sample_vcf(args.in_GT_VCFs, random_samples, lof_variants)
            out.write(f'{args.k_random}\t{n_sampled_variants}\n')


