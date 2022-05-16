import pandas as pd
import argparse
import math
import sys
from rpy2 import robjects
from rpy2.robjects.packages import importr


argparser = argparse.ArgumentParser(description = 'Performs Fisher\'s exact test on genotype counts (i.e. frequencies) between two studies stratified by the groups. Groups for each variant must be specified in the \"GROUP\" column preceding all outher columns in the PLINK *.freqx file. To generate *.freqx files use PLINK with options --freqx and --output-chr chrMT.')
argparser.add_argument('-f1', '--freqx-1', metavar = 'file', dest = 'in_freqx1', type = str, required = True, help = '*.freqx file 1.')
argparser.add_argument('-f2', '--freqx-2', metavar = 'file', dest = 'in_freqx2', type = str, required = True, help = '*.freqx file 2.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_prefix', type = str, required = True, help = 'Output files prefix.')


AUTOSOMALS = set(map(lambda x: f'chr{x}', range(1, 23)))
CHR_X = 'chrX'
CHR_Y = 'chrY'
PAR_XY = 'chrXY'


stats = importr('stats')
rstats_fisher = robjects.r['fisher.test']
rstats_pchisq = robjects.r['pchisq']
rstats_qnorm = robjects.r['qnorm']
rstats_pnorm = robjects.r['pnorm']


def load_freq(filename):
    df = pd.read_csv(filename, header = 0, sep = '\t', low_memory = False)
    print(f'File: {filename}')
    print(f'\tAutosomal SNPs: {len(df[df.CHR.isin(AUTOSOMALS)].SNP.unique()):,}')
    print(f'\tAll X SNPs: {len(df[df.CHR == CHR_X].SNP.unique()):,}')
    print(f'\tAll Y SNPs: {len(df[df.CHR == CHR_Y].SNP.unique()):,}')
    print(f'\tPAR X/Y SNPs: {len(df[df.CHR == PAR_XY].SNP.unique()):,}\n')
    return df


def merge_freq(df1, df2):
    df = df1.merge(df2, on = ['GROUP', 'CHR', 'SNP'], suffixes = ['_1', '_2'])
    print(f'Shared:')
    print(f'\tAutosomal SNPs: {len(df[df.CHR.isin(AUTOSOMALS)].SNP.unique()):,}')
    print(f'\tAll X SNPs: {len(df[df.CHR == CHR_X].SNP.unique()):,}')
    print(f'\tAll Y SNPs: {len(df[df.CHR == CHR_Y].SNP.unique()):,}')
    print(f'\tPAR X/Y SNPs: {len(df[df.CHR == PAR_XY].SNP.unique()):,}\n')
    return df


def get_freq_diploid(hethom_counts):
    n = sum(hethom_counts) * 2
    if n == 0:
        return 0.0, 0, 0
    n_alt = hethom_counts[0] * 2 + hethom_counts[1]
    freq = n_alt / n
    return freq, n_alt, n


def get_freq_haploid(hom_counts):
    n = sum(hom_counts)
    if n == 0:
        return 0.0, 0, 0
    freq = hom_counts[0] / n
    return freq, hom_counts[0] * 2, n * 2


def test_diploid(df):
    for index, row in df[df.CHR.isin(AUTOSOMALS) | (df.CHR == CHR_X) | (df.CHR == PAR_XY)].iterrows():
        table2x3 = [row['C(HOM A1)_1'], row['C(HET)_1'], row['C(HOM A2)_1']]
        if row.A1_1 == row.A1_2:
            table2x3 += [row['C(HOM A1)_2'], row['C(HET)_2'], row['C(HOM A2)_2']]
        elif row.A1_1 ==  row.A2_2:
            table2x3 += [row['C(HOM A2)_2'], row['C(HET)_2'], row['C(HOM A1)_2']]
        else:
            raise Exception(f'Alleles don\'t match for {row.SNP}.')
        table2x3_r = robjects.r['matrix'](robjects.FloatVector(table2x3), nrow = 2, byrow = 'TRUE')
        res = rstats_fisher(table2x3_r, hybrid = 'TRUE', workspace = 20000000)
        p = max(res[0][0], sys.float_info.min)
        freq1, n_alt1, n1 = get_freq_diploid(table2x3[:3])
        freq2, n_alt2, n2 = get_freq_diploid(table2x3[3:])
        yield {'GROUP': row.GROUP, 'SNP': row.SNP, 'N_ALT1': n_alt1, 'N1': n1, 'N_ALT2': n_alt2, 'N2': n2, 'PVALUE': p}


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
        p = max(res[0][0], sys.float_info.min)
        freq1, n_alt1, n1 = get_freq_haploid(table2x2[:2])
        freq2, n_alt2, n2 = get_freq_haploid(table2x2[2:])
        yield {'GROUP': row.GROUP, 'SNP': row.SNP, 'N_ALT1': n_alt1, 'N1': n1, 'N_ALT2': n_alt2, 'N2': n2, 'PVALUE': p}


def meta_analysis(df):
    for name, group in df.groupby('SNP'):
        k = len(group)
        if sum(group.N1) !=0 :
            freq1 = sum(group.N_ALT1) / sum(group.N1)
        else :
            freq1 = 0.0
        if sum(group.N2) !=0 :
            freq2 = sum(group.N_ALT2) / sum(group.N2)
        else :
            freq2 = 0.0

        # Fisher's meta-analysis
        sumlog = sum(group['PVALUE'].apply(lambda x: -2 * math.log10(x))) 
        sumlog_p = rstats_pchisq(sumlog, 2 * k, lower_tail = False)[0]
    
        # Stouffer's meta-analysis (p-values weighted by sample size)
        weights = list(group.apply(lambda r: r.N1 + r.N2, axis = 1))
        z_scores = list(group['PVALUE'].apply(lambda x: rstats_qnorm(x, lower_tail = False)[0]))
        sumz = sum( math.sqrt(w) * z for w, z in zip(weights, z_scores)) / math.sqrt(sum(weights))
        sumz_p = rstats_pnorm(sumz, lower_tail = False)[0]

        yield {'SNP': name, 'K': k, 'FREQ1': freq1, 'FREQ2': freq2, 'SUMLOG': sumlog, 'SUMLOG_PVALUE': sumlog_p, 'SUMZ': sumz, 'SUMZ_PVALUE': sumz_p}



if __name__ == '__main__':
    args = argparser.parse_args()

    df1 = load_freq(args.in_freqx1) 
    df2 = load_freq(args.in_freqx2)

    merged = merge_freq(df1, df2)

    diploid_res = pd.DataFrame(test_diploid(merged))
    diploid_res.to_csv(f'{args.out_prefix}.diploid.txt', sep ='\t', index = False) 

    diploid_meta = pd.DataFrame(meta_analysis(diploid_res))
    diploid_meta.to_csv(f'{args.out_prefix}.diploid.meta.txt', sep = '\t', index = False)

    haploid_res = pd.DataFrame(test_haploid(merged))
    haploid_res.to_csv(f'{args.out_prefix}.haploid.txt', sep = '\t', index = False)

    haploid_meta = pd.DataFrame(meta_analysis(haploid_res))
    haploid_meta.to_csv(f'{args.out_prefix}.haploid.meta.txt', sep = '\t', index = False)

    print('Done.')

