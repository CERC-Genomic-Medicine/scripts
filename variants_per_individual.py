import argparse
import pysam
from collections import Counter

argparser = argparse.ArgumentParser(description = 'Variant counts by individual from VCF/BCF files.')
argparser.add_argument('-g', '--genotypes', metavar = 'file', dest = 'in_VCF_GT', required = True, help = 'Input VCF/BCF file with genotypes in GT format field.')
argparser.add_argument('-a', '--annotations', metavar = 'file', dest = 'in_VCF_CSQ', required = True, help = 'Input VCF/BCF file with annotations in CSQ INFO field (e.g. from Variant Effect Predictor). Requires tabix index. Genotypes are not required.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_file', required = True, help = 'Output tab-separated file.')

# Consequences in order of severity from http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequence_type_table
CONSEQUENCES = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'incomplete_terminal_codon_variant',
    'stop_retained_variant',
    'start_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant',
]
CONSEQUENCES_SEVERITY = {consequence:i for i, consequence in enumerate(CONSEQUENCES)}


def read_vcf_csq(in_vcf, chrom, start_bp, stop_bp):
    annotated_variants = dict()
    if start_bp is not None and start_bp > 0: # 0-based start position
        start_bp -= 1
    with pysam.VariantFile(in_vcf, drop_samples = True) as ifile:
        csq_meta = ifile.header.info.get('CSQ', None)
        if csq_meta is None:
            raise Exception('No meta-information entry about CSQ INFO field found!')
        csq_header = csq_meta.description.split(':', 1)[1].strip().split('|')
        csq_allele_num_index = csq_header.index('ALLELE_NUM')
        csq_consequence_index = csq_header.index('Consequence')
        if ifile.get_tid(chrom) < 0: # append or remove 'chr' prefix if neccessary
            if chrom.startswith('chr'):
                fetched = ifile.fetch(chrom[3:], start_bp, stop_bp)
            else:
                fetched = ifile.fetch(f'chr{chrom}', start_bp, stop_bp)
        else:
            fetched = ifile.fetch(chrom, start_bp, stop_bp)
        for record in fetched:
            name_prefix = f'{record.chrom}_{record.pos}_{record.ref}'
            severity_by_allele = dict()
            for csq in record.info.get('CSQ'):
                csq = csq.split('|')
                severity = min(CONSEQUENCES_SEVERITY[x] for x in csq[csq_consequence_index].split('&'))
                allele_num = int(csq[csq_allele_num_index]) - 1
                allele_severity = severity_by_allele.get(allele_num, None)
                if allele_severity is None:
                    severity_by_allele[allele_num] = severity
                else:
                    severity_by_allele[allele_num] = min(severity, allele_severity)
            for alt_index, severity in severity_by_allele.items():
                annotated_variants[f'{name_prefix}_{record.alts[alt_index]}'] = CONSEQUENCES[severity]
    return annotated_variants


def read_gt_vcf(in_gt_vcf, in_csq_vcf, out_file):
    annotations_chrom = ''
    annotations_end_bp = -1
    annotations = None
    annotations_size = 500000
    individual_counts = dict()
    total_counts = Counter()
    with pysam.VariantFile(in_gt_vcf) as ifile:
        for individual in ifile.header.samples:
            individual_counts[individual] = Counter()
        for record in ifile:
            if 'PASS' not in record.filter: # do not use non-PASS
                continue
            # load a chunk of annotated variants
            if annotations_chrom != record.chrom or annotations_end_bp < record.pos:
                annotations_chrom = record.chrom
                annotations_end_bp = record.pos + annotations_size
                annotations = read_vcf_csq(in_csq_vcf, annotations_chrom, record.pos, annotations_end_bp)

            # count alt carriers
            alt_carriers = dict((i, Counter()) for i, alt in enumerate(record.alts, 1)) # we assume that all variants are normalized such that ALT is never equal to '.'
            for individual, genotype in record.samples.iteritems():
                for i in genotype.allele_indices:
                    if i > 0:
                        alt_carriers[i].update((individual,))
            for alt_index, individuals in alt_carriers.items():
                j += 1
                ac = sum(individuals.values())
                if ac == 0: # skip monomorphic
                    continue
                consequence = annotations.get(f'{record.chrom}_{record.pos}_{record.ref}_{record.alts[alt_index - 1]}', None)
                if consequence is None:
                    continue
                for individual in individuals:
                    individual_counts[individual][consequence] += 1
                total_counts[consequence] += 1
    with open(out_file, 'w') as ofile:
        ofile.write('ID\t{}\n'.format('\t'.join(CONSEQUENCES)))
        for individual, counts in individual_counts.items():
            ofile.write('{}\t{}\n'.format(individual, '\t'.join([str(counts[c]) for c in CONSEQUENCES])))
        ofile.write('TOTAL\t{}\n'.format('\t'.join(str(total_counts[c]) for c in CONSEQUENCES)))



if __name__ == '__main__':
    args = argparser.parse_args()
    read_gt_vcf(args.in_VCF_GT, args.in_VCF_CSQ, args.out_file)

