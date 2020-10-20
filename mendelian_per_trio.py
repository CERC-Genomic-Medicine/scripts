#!/usr/bin/env python3

import argparse
import pysam
import gzip
from collections import namedtuple, OrderedDict

argparser = argparse.ArgumentParser(description = 'Counts per-trio Mendelian errors (in autosomals only).')
argparser.add_argument('-p', '--pedigree', metavar = 'file', dest = 'in_ped_file', type = str, required = True, help = 'PED file with nuclear families (no header).')
argparser.add_argument('-g', '--genotypes', metavar = 'file', dest = 'in_vcf_files', type = str, required = True, nargs = '+', help = 'VCF/BCF file with genotypes. One file per chromosome.')
argparser.add_argument('-n', '--min-an', metavar = 'number', dest = 'min_an', type = int, required = True, help = 'Minimal number of alleles per variant (AN INFO field).')
argparser.add_argument('-i', '--indels', action = 'store_true', help = 'Count only Indels. Otherwise count only SNPs.')
argparser.add_argument('-o', '--output', metavar = 'name', dest = 'output_prefix', type = str, required = True, help = 'Output file prefix.')

Trio = namedtuple('Trio', ['iid', 'fid', 'mid', 'family', 'counts'])

def read_ped(filename):
    with open(filename, 'r') as ifile:
        for line in ifile:
            fields = line.rstrip().split()
            if len(fields) < 4:
                raise Exception('PED file entry with <4 columns was detected! Required >=4 columns.')
            family_id = fields[0]
            individual_id = fields[1]
            father_id = fields[2]
            mother_id = fields[3]
            if father_id == '0' or mother_id == '0':
                continue
            yield Trio(individual_id, father_id, mother_id, family_id, { 'n_snps': 0, 'n_snps_errors': 0, 'variants': OrderedDict() } )


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


if __name__ == '__main__':
    args = argparser.parse_args()
    trios = [ t for t in read_ped(args.in_ped_file) ]
    trios_variants = dict()
    for trio in trios:
        trios_variants[trio.iid] = []
    for filename in args.in_vcf_files:
        with pysam.VariantFile(filename) as ifile:
            for info_field in ['AC', 'AN']:
                if info_field not in ifile.header.info:
                    raise Exception(f'Meta-information line for {info_field} INFO field is missing in VCF/BCF header.')
            for record in ifile.fetch():
                if len(record.alts) > 1: # skip multi-allelic
                    continue
                alt = record.alts[0]
                if alt == '*': # skip spanning deletions
                    continue
                if not args.indels:
                    if len(record.ref) != 1 or len(alt) != 1: # keep only SNPs
                        continue
                else:
                    if len(record.ref) == 1 and len(alt) == 1: # keep only Indels
                        continue
                an = record.info['AN']
                ac = record.info['AC'][0]
                if ac == 0 or ac == record.info['AN']: # skip monomorphic
                    continue
                chrom = record.chrom[3:] if record.chrom.startswith('chr') else record.chrom
                name = f'{chrom}-{record.pos}-{record.ref}-{alt}'
                is_pass = 'PASS' in list(record.filter) and an >= args.min_an
                af = ac / an
                if af <= 0.5:
                    mac = ac
                    maf = af
                else:
                    mac = an - ac
                    maf = 1.0 - af
                for trio in trios:
                    igt = record.samples[trio.iid]['GT']
                    if None in igt:
                        continue
                    igt = sum(igt)
                    fgt = record.samples[trio.fid]['GT']
                    if None in fgt:
                        continue
                    fgt = sum(fgt)
                    mgt = record.samples[trio.mid]['GT']
                    if None in mgt:
                        continue
                    mgt = sum(mgt)
                    if igt == 0 and fgt == 0 and mgt == 0:
                        continue
                    merror = is_mendel_inconsistent(fgt, mgt, igt)
                    trios_variants[trio.iid].append((name, an, ac, mac, maf, is_pass, merror))
    for iid, variants in trios_variants.items():
        filename = f'{args.output_prefix}.{iid}.txt.gz'
        with gzip.open(filename, 'wt') as ofile:
            ofile.write('VARIANT\tAN\tAC\tMAC\tMAF\tPASS\tMERROR\n')
            for v in variants:
                ofile.write(f'{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[4]}\t{1 if v[5] else 0}\t{1 if v[6] else 0}\n')

