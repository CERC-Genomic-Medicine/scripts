import argparse
import pysam
import numpy as np
import gzip

argparser = argparse.ArgumentParser(description = 'Merges minimac4 imputation results and recomputes Rsq imputation quality score.')
argparser.add_argument('-i1', '--in-vcf1', metavar = 'file', dest = 'in_VCF1', required = True, help = 'First imputed VCF file.')
argparser.add_argument('-i2', '--in-vcf2', metavar = 'file', dest = 'in_VCF2', required = True, help = 'Second imputed VCF file.')
argparser.add_argument('-c', '--chromosomes', metavar = 'name', dest = 'chrom', required = True, type = str, help = 'Chromosome name.')
argparser.add_argument('-b', '--begin', metavar = 'position', dest = 'begin', required = True, type = int, help = 'Begin position.')
argparser.add_argument('-e', '--end', metavar = 'position', dest = 'end', required = True, type = int , help = 'End position.')
argparser.add_argument('-m', '--min-rsq', metavar = 'float', dest = 'min_rsq', required = True, type = float, help = 'Minimal combined R2 to write to final VCF.')
argparser.add_argument('-o', '--out', metavar = 'file', dest = 'out_VCF', required = True, help = 'Output merged VCF file.')


if __name__ == '__main__':
   args = argparser.parse_args()

   with pysam.VariantFile(args.in_VCF1, 'r') as vcf1, pysam.VariantFile(args.in_VCF2, 'r') as vcf2, gzip.open(args.out_VCF, "wt") as vcf_merged:
      if not 'R2' in vcf1.header.info:
         raise Exception(f'No meta-information entry about R2 INFO field found in {args.in_VCF1}!')
      if not 'R2' in vcf2.header.info:
         raise Exception(f'No meta-information entry about R2 INFO field found in {args.in_VCF2}!')
      if not 'TYPED_ONLY' in vcf1.header.info:
         raise Exception(f'No meta-information entry about TYPED_ONLY INFO field found in {args.in_VCF1}!')
      if not 'TYPED_ONLY' in vcf2.header.info:
         raise Exception(f'No meta-information entry about TYPED_ONLY INFO field found in {args.in_VCF2}!')

      samples1 = list(vcf1.header.samples)
      samples2 = list(vcf2.header.samples)
      samples_shared = set([s1 for s1 in samples1 if s1 in samples2])
      samples_combined = samples1 + [s2 for s2 in samples2 if s2 not in samples1]

      n = len(samples_combined)

      print(f'Samples in {args.in_VCF1}: {len(samples1)}')
      print(f'Samples in {args.in_VCF2}: {len(samples2)}')
      print(f'Samples shared: {len(samples_shared)}')
      print(f'Samples combined: {n}')

      hds_combined = np.zeros(n * 2, dtype = np.float64)

      vcf_merged.write('##fileformat=VCFv4.1\n')
      for c in vcf1.header.contigs:
         vcf_merged.write(f'##contig=<ID={c}>\n')
      vcf_merged.write('##INFO=<ID=VCF1_AF,Number=1,Type=Float,Description="Estimated Alternate Allele Frequency in VCF 1">\n')
      vcf_merged.write('##INFO=<ID=VCF1_R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square) in VCF 1">\n')
      vcf_merged.write('##INFO=<ID=VCF2_AF,Number=1,Type=Float,Description="Estimated Alternate Allele Frequency in VCF 2">\n')
      vcf_merged.write('##INFO=<ID=VCF2_R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square) in VCF 2">\n')
      vcf_merged.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Combined Estimated Alternate Allele Frequency">\n')
      vcf_merged.write('##INFO=<ID=MAF,Number=1,Type=Float,Description="Combined Estimated Minor Allele Frequency">\n')
      vcf_merged.write('##INFO=<ID=R2,Number=1,Type=Float,Description="Combined Estimated Imputation Accuracy (R-square)">\n')
      vcf_merged.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
      vcf_merged.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]">\n')
      vcf_merged.write(f'##vcf1={args.in_VCF1}\n')
      vcf_merged.write(f'##vcf2={args.in_VCF2}\n')
      vcf_merged.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format('\t'.join(samples_combined)))

      for record1, record2 in zip(vcf1.fetch(contig = args.chrom, start = args.begin, stop = args.end), vcf2.fetch(contig = args.chrom, start = args.begin, stop = args.end)):
         assert record1.chrom == record2.chrom and record1.pos == record2.pos and record1.ref == record2.ref and record1.alts[0] == record2.alts[0], f'{record1.chrom}:{record1.pos}:{record1.ref}:{record1.alts[0]} vs {record2.chrom}:{record2.pos}:{record2.ref}:{record2.alts[0]}'         

         af1 = record1.info['AF']
         af2 = record2.info['AF']
         rsq1 = record1.info['R2']
         rsq2 = record2.info['R2']

         genotypes = []
         dosages = []
 
         for i, sample in enumerate(samples_combined):
            frmt1 = record1.samples.get(sample, None)
            frmt2 = record2.samples.get(sample, None)
 
            if frmt1 is None:
               phased = frmt2.phased
               gt = frmt2['GT']
               ds = frmt2['DS']
               hds = frmt2['HDS']
            elif frmt2 is None:
               phased = frmt1.phased
               gt = frmt1['GT']
               ds = frmt1['DS']
               hds = frmt1['HDS']
            else:
               if rsq1 >= rsq2:
                  phased = frmt1.phased
                  gt = frmt1['GT']
                  ds = frmt1['DS']
                  hds = frmt1['HDS']
               else:
                  phased = frmt2.phased
                  gt = frmt2['GT']
                  ds = frmt2['DS']
                  hds = frmt2['HDS']
             
            genotypes.append(f'{gt[0]}|{gt[1]}' if phased else f'{gt[0]}/{gt[1]}')
            dosages.append(ds)
            hds_combined[i * 2] = hds[0]
            hds_combined[i * 2 + 1] = hds[1]
         
         p = np.mean(hds_combined)
         if p > 0:
            rsq_combined = np.mean(np.square(hds_combined - p)) / (p * (1.0 - p))
         else:
            rsq_combined = 0.0

         if rsq_combined < args.min_rsq:
            continue

         vcf_merged.write('{chrom}\t{pos}\t{vid}\t{ref}\t{alts}\t{qual}\t{filters}\tVCF1_AF={af1:.5f};VCF2_AF={af2:.5f};VCF1_R2={rsq1:.5f};VCF2_R2={rsq2:.5f};AF={af_combined:.5f};MAF={maf_combined:.5f};R2={rsq_combined:.5f}\t{formats}\t{gt_ds}\n'.format(
               chrom = record1.chrom,
               pos = record1.pos,
               vid = record1.id,
               ref = record1.ref,
               alts = ','.join(record1.alts),
               qual = '.',
               filters = 'PASS',
               af1 = af1,
               af2 = af2,
               rsq1 = rsq1,
               rsq2 = rsq2,
               af_combined = p,
               maf_combined = p if p < 0.5 else 1.0 - p,
               rsq_combined = rsq_combined,
               formats = 'GT:DS',
               gt_ds = '\t'.join(f'{gt}:{ds:.3f}' for gt, ds in zip(genotypes, dosages)  )
               ))
