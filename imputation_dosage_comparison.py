#!/usr/bin/env python3
import argparse
import pysam
import numpy as np
import gzip
import pandas as pd

argparser = argparse.ArgumentParser(description = 'Compares imputation dosage of overlap samples per samples and per records')
argparser.add_argument('-i1', '--in-vcf1', metavar = 'file', dest = 'in_VCF1', required = True, help = 'First imputed VCF file.')
argparser.add_argument('-i2', '--in-vcf2', metavar = 'file', dest = 'in_VCF2', required = True, help = 'Second imputed VCF file.')
argparser.add_argument('-s', '--sample', metavar = 'file',  dest = 'sample', required = False, help = 'Sample file. (Optional) (default all shared samples are considered)')
argparser.add_argument('-c', '--chromosomes', metavar = 'name', dest = 'chrom', required = True, type = str, help = 'Chromosome name.')
argparser.add_argument('-b', '--begin', metavar = 'position', dest = 'begin', required = True, type = int, help = 'Begin position.')
argparser.add_argument('-e', '--end', metavar = 'position', dest = 'end', required = True, type = int , help = 'End position.')
argparser.add_argument('-o', '--out', metavar = 'file', dest = 'out', required = True, help = 'Output perfix.')

### Output
# o prefix_per_rec produces list of SNP and the sum of the absolute dosage difference per samples (\tab delimited)
# o prefix_per_samples produces 1 column per sample with a single row of the associated sum of absolute difference for the region

if __name__ == '__main__':
    args = argparser.parse_args()
    with pysam.VariantFile(args.in_VCF1, 'r') as vcf1, pysam.VariantFile(args.in_VCF2, 'r') as vcf2, open(args.out+"_per_rec", 'w') as R , open(args.out+"_per_sample", 'w') as f:
        ### create overlap list ###
        samples1 = list(vcf1.header.samples)
        samples2 = list(vcf2.header.samples)
        samples_shared = [s1 for s1 in samples1 if s1 in samples2]
        ### take the potential list of wanted samples ###
        if args.sample :
            df_samples=pd.read_csv(args.sample, header=None, names=['samples'])
            samples_shared = [str(s) for s in df_samples['samples'] if str(s) in samples_shared]
        f.write('\t' + '\t'.join(samples_shared) + '\n' ) ### header of per samples output
        sampleDosages = [0 for number in range(len(samples_shared))] ### baseline of samples
        ### iter per record
        for record1, record2 in zip(vcf1.fetch(contig = args.chrom, start = args.begin, stop = args.end), vcf2.fetch(contig = args.chrom, start = args.begin, stop = args.end)):
            ### test for possible incorrect input file
            assert record1.chrom == record2.chrom and record1.pos == record2.pos and record1.ref == record2.ref and record1.alts[0] == record2.alts[0], f'{record1.chrom}:{record1.pos}:{record1.ref}:{record1.alts[0]} vs {record2.chrom}:{record2.pos}:{record2.ref}:{record2.alts[0]}'
            dosages = [] # initiallise
            if record1.pos < args.begin.: # correcting behaviour for border overlapping variants 
                continue
            ## iter per samples
            for i, sample in enumerate(samples_shared):
                frmt1 = record1.samples.get(sample, None)
                frmt2 = record2.samples.get(sample, None)
                dosages.append(abs(frmt2['DS']-frmt1['DS'])) ## list of dosage differences per samples
                #write each variants
            R.write(record1.id + '\t' +str(sum(dosages))  + '\t' + str(sum(dosages)/len(dosages))+'\t'+str(len(samples_shared)) +'\n')
            ### add variant result to samples dosage differences.
            sampleDosages=np.add(sampleDosages,dosages).tolist()
        f.write('\t'.join([str(q) for q in sampleDosages]) + '\n')
