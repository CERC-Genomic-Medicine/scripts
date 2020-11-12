#!/usr/bin/env python3

import argparse
import pysam
import gzip
from collections import namedtuple, OrderedDict

argparser = argparse.ArgumentParser(description = 'Counts per-trio Mendelian errors (in autosomals only).')
argparser.add_argument('-f', '--father', metavar = 'file', dest = 'in_father_vcf', type = str, required = True, help = 'VCF/BCF file with genotypes for father.')
argparser.add_argument('-m', '--mother', metavar = 'file', dest = 'in_mother_vcf', type = str, required = True, help = 'VCF/BCF file with genotypes for mother.')
argparser.add_argument('-c', '--child', metavar = 'file', dest = 'in_child_vcf', type = str, required = True, help = 'VCF/BCF file with genotypes for child.')
argparser.add_argument('-a', '--all', action = 'store_true', help = 'Count all PASS and non-PASS variants. Otherwise count only PASS variants.')
argparser.add_argument('-i', '--indels', action = 'store_true', help = 'Count only Indels. Otherwise count only SNPs.')
argparser.add_argument('-o', '--output', metavar = 'name', dest = 'output_prefix', type = str, required = True, help = 'Output file prefix.')

Trio = namedtuple('Trio', ['iid', 'fid', 'mid', 'family', 'counts'])


def is_mendel_inconsistent(father, mother, child):
    # based on PLINK documentation
    if father == 0 and mother == 0 and child == 1:
        return True # all
    elif father == 0 and mother == 0 and child == 2:
        return True # child
    elif father == 0 and mother > 0 and child == 2:
        return True # father, child
    elif father == 2 and mother == 2 and child == 1:
        return True # all
    elif father == 2 and mother == 2 and child == 0:
        return True # child
    elif father == 2 and mother < 2 and child == 0:
        return True # father, child
    elif father < 2 and mother == 2 and child == 0:
        return True # mother, child
    elif father > 0 and mother == 0 and child == 2:
        return True # mother, child
    return False


def count_independent(ordered_variants, max_dist):
    prev_pos = -1 * max_dist - 1
    group = 0
    groups = {}
    for variant, info in ordered_variants.items():
        pos = int(variant.split('-', 2)[1])
        if pos - prev_pos <= max_dist:
            groups[group].append(variant)
        else:
            group = pos
            groups[group] = [variant]
        prev_pos = pos

    n_independent = 0
    n_independent_with_errors = 0
    n_nogroup = 0
    n_nogroup_with_errors = 0
    for group, variants in groups.items():
        if len(variants) > 1:
            for v in variants:
                ordered_variants[v]['group'] = group
        else:
            n_nogroup += 1
            if ordered_variants[variants[0]]['errors'] > 0:
                n_nogroup_with_errors += 1

        if all(ordered_variants[v]['errors'] == 0 for v in variants): # count as single independent variant without error
            n_independent += 1
        elif all(ordered_variants[v]['errors'] > 0 for v in variants): # count as single independent variant with error
            n_independent += 1
            n_independent_with_errors += 1
        else: # count all as independent
            n_independent += len(variants)
            n_independent_with_errors += sum(ordered_variants[v]['errors'] > 0 for v in variants)

    return (n_independent, n_independent_with_errors, n_nogroup, n_nogroup_with_errors)


def load_variants(in_vcf, only_indels, variants, trio_sample):
    auto_chroms = set(f'{c}' for c in range(0, 23))
    with pysam.VariantFile(in_vcf) as ifile:
        if len(ifile.header.samples) > 1:
            raise Exception(f'Multiple samples in {in_vcf} found. Only one sample is allowed.')
        for record in ifile.fetch():
            if len(record.alts) > 1: # skip multi-allelic
                continue
            alt = record.alts[0]
            if alt == '*': # skip spanning deletions
                continue
            if not only_indels:
                if len(record.ref) != 1 or len(alt) != 1: # keep only SNPs
                    continue
            else:
                if len(record.ref) == 1 and len(alt) == 1: # keep only Indels
                    continue
            chrom = record.chrom[3:] if record.chrom.startswith('chr') else record.chrom
            if chrom not in auto_chroms:
                continue
            name = f'{chrom}-{record.pos}-{record.ref}-{alt}'
            is_pass = 'PASS' in list(record.filter)
            gt = record.samples[0]['GT']
            if not None in gt:
                gt = sum(gt)
            variant = variants.setdefault(name, [None, None, None, None, None, None])
            if trio_sample == 'father':
                variant[0] = gt
                variant[3] = is_pass
            elif trio_sample == 'mother':
                variant[1] = gt
                variant[4] = is_pass
            elif trio_sample == 'child':
                variant[2] = gt
                variant[5] = is_pass
            else:
                raise Exception(f'Unsopported argument.')


if __name__ == '__main__':
    args = argparser.parse_args()
    variants = dict()

    load_variants(args.in_father_vcf, args.indels, variants, 'father')
    load_variants(args.in_mother_vcf, args.indels, variants, 'mother')
    load_variants(args.in_child_vcf, args.indels, variants, 'child')

    filename = f'{args.output_prefix}.txt.gz'
    with gzip.open(filename, 'wt') as ofile:
        ofile.write('VARIANT\tALL_PASS\tMENDELIAN_ERROR\n')
        for variant, data in variants.items():
            fgt, mgt, cgt = data[:3]
            if fgt is None or mgt is None or cgt is None:
                continue
            fpass, mpass, cpass = data[3:]
            is_pass = fpass and mpass and cpass
            if not args.all:
                if not is_pass:
                    continue
            merror = is_mendel_inconsistent(fgt, mgt, cgt)
            ofile.write(f'{variant}\t{1 if is_pass else 0}\t{1 if merror else 0}\n')

