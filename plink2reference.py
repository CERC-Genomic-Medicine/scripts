#!/usr/bin/env python3

import argparse
import pysam

# After running this script, you can align your PLINK files to reference using the following sequence of steps:
# 1) plink --bfile <original PLINK file> --exclude <prefix>.remove.txt --make-bed --out temp1
# 2) plink --bfile temp1 --flip <prefix>.strand_flip.txt --make-bed --out temp2
# 3) plink --bfile temp2 --a1-allele <prefix>.force_a1.txt --make-bed --out <new PLINK file>
# 4) rm temp1.* temp2.*

argparser = argparse.ArgumentParser(description = 'Creates lists of variants for aligning A1/A2 PLINK alleles to reference.')
argparser.add_argument('-b', '--plink-bim', metavar = 'file', dest = 'in_plink_bim', type = str, required = True, help = 'Input PLINK bim file.')
argparser.add_argument('-f', '--fasta', metavar = 'file', dest = 'in_ref_fasta', type = str, required = True, help = 'Input reference FASTA file (indexed).')
argparser.add_argument('-o', '--output-prefix', metavar = 'file', dest = 'out_prefix', type = str, required = True, help = 'Output file prefix.')


def strand_flip(a):
    return { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A' }[a] 


if __name__ == '__main__':
    args = argparser.parse_args()

    n_total_variants = 0
    n_non_snps = 0
    n_palindromic = 0
    n_flip_strand = 0
    n_force_ref_allele = 0
    n_no_ref_match = 0 

    with open(args.in_plink_bim, 'rt') as ibim, pysam.FastaFile(args.in_ref_fasta) as ifasta, \
        open(f'{args.out_prefix}.remove.txt', 'wt') as oremove, open(f'{args.out_prefix}.strand_flip.txt', 'wt') as oflip, open(f'{args.out_prefix}.force_a1.txt', 'wt') as oforce:

        fasta_chroms = set(list(ifasta.references))
        for line in ibim:
            fields = line.rstrip().split()
            chrom, varid, pos, a1, a2 = fields[0], fields[1], int(fields[3]), fields[4], fields[5]
            
            n_total_variants += 1

            if chrom not in fasta_chroms:
                chrom = chrom[3:] if chrom.startswith('chr') else f'chr{chrom}'
                if chrom not in fasta_chroms:
                    print(f'Warning: skipping chromosome {fields[0]} because it is not in FASTA file.')
                    continue

            if a1 not in {'A', 'C', 'G', 'T'} or a2 not in {'A', 'C', 'G', 'T'}:
                oremove.write(f'{varid}\n')
                n_non_snps += 1
                continue

            if (a1 in {'A', 'T'} and a2 in {'A', 'T'}) or (a1 in {'C', 'G'} and a2 in {'C', 'G'}):
                oremove.write(f'{varid}\n')
                n_palindromic += 1
                continue


            ref_base = None
            for base in ifasta.fetch(chrom, pos - 1, pos):
                ref_base = base

            if ref_base == a2:
                n_force_ref_allele += 1
                oforce.write(f'{varid}\t{ref_base}\n')
            elif ref_base != a1:
                flipped_a1 = strand_flip(a1)
                flipped_a2 = strand_flip(a2)
                if ref_base == flipped_a2:
                    n_force_ref_allele += 1
                    oforce.write(f'{varid}\t{ref_base}\n')
                    n_flip_strand += 1
                    oflip.write(f'{varid}\n')
                elif ref_base == flipped_a1:
                    n_flip_strand += 1
                    oflip.write(f'{varid}\n')
                else:
                    n_no_ref_match += 1
                    oremove.write(f'{varid}\n')

    print(f'Total variants: {n_total_variants:,}')
    print(f'Not valid SNPs: {n_non_snps:,}')
    print(f'Palindromic SNPs: {n_palindromic:,}')
    print(f'Strand flips: {n_flip_strand:,}')
    print(f'Force reference allele: {n_force_ref_allele:,}')
    print(f'A1/A2 didn\'t match reference allele: {n_no_ref_match:,}')
