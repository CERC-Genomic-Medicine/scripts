import pandas as pd
import argparse
from rpy2 import robjects
from rpy2.robjects.packages import importr

# Before running this script, generate *.freqx files:
# plink --bfile <PLINK files> --freqx --output-chr chrMT --out <output file>.freqx


argparser = argparse.ArgumentParser(description = 'Performs Fisher\'s exact test on genotype counts (i.e. frequencies) between two studies.')
argparser.add_argument('-f1', '--freqx-1', metavar = 'file', dest = 'in_freqx1', type = str, required = True, help = '*.freqx file 1.')
argparser.add_argument('-f2', '--freqx-2', metavar = 'file', dest = 'in_freqx2', type = str, required = True, help = '*.freqx file 2.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_prefix', type = str, required = True, help = 'Output files prefix.')

AUTOSOMALS = set(map(lambda x: f'chr{x}', range(1, 23)))
CHR_X = 'chrX'
CHR_Y = 'chrY'
PAR_XY = 'chrXY'
CHR_MT = 'chrMT'

stats = importr('stats')
rstats_fisher = robjects.r['fisher.test']


def load_freq(filename):
    df = pd.read_csv(filename, header = 0, sep = '\t', low_memory = False)
    print(f'File: {filename}')
    print(f'\tAutosomal SNPs: {len(df[df.CHR.isin(AUTOSOMALS)]):,}')
    print(f'\tAll X SNPs: {len(df[df.CHR == CHR_X]):,}')
    print(f'\tAll Y SNPs: {len(df[df.CHR == CHR_Y]):,}')
    print(f'\tAll MT SNPs: {len(df[df.CHR == CHR_MT]):,}')
    print(f'\tPAR X/Y SNPs: {len(df[df.CHR == PAR_XY]):,}\n')
    return df


def merge_freq(df1, df2):
    df = df1.merge(df2, on = ['CHR', 'SNP'], suffixes = ['_1', '_2'])
    print(f'Shared:')
    print(f'\tAutosomal SNPs: {len(df[df.CHR.isin(AUTOSOMALS)]):,}')
    print(f'\tAll X SNPs: {len(df[df.CHR == CHR_X]):,}')
    print(f'\tAll Y SNPs: {len(df[df.CHR == CHR_Y]):,}')
    print(f'\tAll MT SNPs: {len(df[df.CHR == CHR_MT]):,}')
    print(f'\tPAR X/Y SNPs: {len(df[df.CHR == PAR_XY]):,}\n')
    return df


def get_freq_diploid(hethom_counts):
    n = sum(hethom_counts) * 2
    freq = (hethom_counts[0] * 2 + hethom_counts[1]) / n
    return freq


def get_freq_haploid(hom_counts):
    n = sum(hom_counts)
    freq = hom_counts[0] / n
    return freq


def test_diploid(df):
    for index, row in df[df.CHR.isin(AUTOSOMALS) | (df.CHR == CHR_X) | (df.CHR == CHR_MT) | (df.CHR == PAR_XY)].iterrows():
        table2x3 = [row['C(HOM A1)_1'], row['C(HET)_1'], row['C(HOM A2)_1']]
        if row.A1_1 == row.A1_2:
            table2x3 += [row['C(HOM A1)_2'], row['C(HET)_2'], row['C(HOM A2)_2']]
        elif row.A1_1 ==  row.A2_2:
            table2x3 += [row['C(HOM A2)_2'], row['C(HET)_2'], row['C(HOM A1)_2']]
        else:
            raise Exception(f'Alleles don\'t match for {row.SNP}.')
        table2x3_r = robjects.r['matrix'](robjects.FloatVector(table2x3), nrow = 2, byrow = 'TRUE')
        res = rstats_fisher(table2x3_r, hybrid = 'TRUE', workspace = 20000000)
        p = res[0][0]
        freq1 = get_freq_diploid(table2x3[:3])
        freq2 = get_freq_diploid(table2x3[3:])
        yield {'SNP': row.SNP, 'FREQ1': freq1, 'FREQ2': freq2, 'PVALUE': p}


def test_haploid(df):
    for index, row in df[(df.CHR == CHR_X) | (df.CHR == CHR_Y)].iterrows():
        table2x2 = [row['C(HAP A1)_1'], row['C(HAP A2)_1']]
        if row.A1_1 == row.A1_2:
            table2x2 += [row['C(HAP A1)_2'], row['C(HAP A2)_2']]
        elif row.A1_1 == row.A2_2:
            table2x2 += [row['C(HAP A2)_2'], row['C(HAP A1)_2']]
        else:
            raise Exception(f'Alleles don\'t match for {row.SNP}.')
        table2x2_r = robjects.r['matrix'](robjects.FloatVector(table2x2), nrow = 2, byrow = 'TRUE')
        res = rstats_fisher(table2x2_r, hybrid = 'TRUE', workspace = 20000000)
        p = res[0][0]
        freq1 = get_freq_haploid(table2x2[:2])
        freq2 = get_freq_haploid(table2x2[2:])
        yield {'SNP': row.SNP, 'FREQ1': freq1, 'FREQ2': freq2, 'PVALUE': p}


if __name__ == '__main__':
    args = argparser.parse_args()

    df1 = load_freq(args.in_freqx1)
    df2 = load_freq(args.in_freqx2)
    merged = merge_freq(df1, df2)

    diploid_res = pd.DataFrame(test_diploid(merged))
    diploid_res.to_csv(f'{args.out_prefix}.diploid.txt', sep ='\t', index = False) 

    haploid_res = pd.DataFrame(test_haploid(merged))
    haploid_res.to_csv(f'{args.out_prefix}.haploid.txt', sep = '\t', index = False)

    print('Done.')

