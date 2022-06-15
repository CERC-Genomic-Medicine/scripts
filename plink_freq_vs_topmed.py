import argparse
import pandas as pd
import pysam
import matplotlib.pyplot as plt

# 1) To create *.freq files, run `plink --bfile <PLINK file> --keep-allele-order -freq --out <output file>`
# 2) TOPMed AF are located in VCFs at https://bravo.sph.umich.edu/freeze8/hg38/downloads. 
#    Once downloaded, merge single chromosome VCFs into one VCF and index using `tabix` or `bcftools index`.

argparser = argparse.ArgumentParser(description = 'Compares allele frequencies to TOPMed reference panel.')
argparser.add_argument('-s', '--study', metavar = 'file', dest = 'in_freq', type = str, required = True, help = 'PLINK *.freq file.')
argparser.add_argument('-t', '--topmed', metavar = 'file', dest = 'in_topmed', type = str, required = True, help = 'TOPMed *.vcf.gz file.')
argparser.add_argument('-m', '--max-af-diff', metavar = 'number', dest = 'max_af_diff', type = float, required = True, help = 'Max AF difference e.g. 0.2')
argparser.add_argument('-o', '--out', metavar = 'file', dest = 'out_prefix', type = str, required = True, help = 'Output prefix.')


CACHE = dict()
CACHE_CHROM = ''
CACHE_START_BP = 0
CACHE_STOP_BP = 0


def get_topmed_freq(filename, variant_name):
   global CACHE_CHROM
   global CACHE_START_BP
   global CACHE_STOP_BP

   chrom, position, ref, alt = variant_name.split(':')
   position = int(position)
  
   if CACHE_CHROM != chrom or CACHE_START_BP > position or CACHE_STOP_BP < position:
      CACHE.clear()
      CACHE_CHROM = chrom
      CACHE_START_BP = max(position - 100000, 0)
      CACHE_STOP_BP = position + 5000000
      print(f'LOAD: {CACHE_CHROM}:{CACHE_START_BP}-{CACHE_STOP_BP}')
      with pysam.VariantFile(filename, 'r') as vcf:
         if vcf.get_tid(CACHE_CHROM) < 0: # append or remove 'chr' prefix if neccessary
            if CACHE_CHROM.startswith('chr'):
               fetched = vcf.fetch(CACHE_CHROM[3:], CACHE_START_BP, CACHE_STOP_BP)
            else:
               fetched = vcf.fetch(f'chr{CACHE_CHROM}', CACHE_START_BP, CACHE_STOP_BP)
         else:
            fetched = vcf.fetch(CACHE_CHROM, CACHE_START_BP, CACHE_STOP_BP)
         for record in fetched:
            if len(record.alts) > 1:
               continue
            if 'PASS' not in record.filter:
               continue
            CACHE[f'{CACHE_CHROM}:{record.pos}:{record.ref}:{record.alts[0]}'] = record.samples['TOPMED']['FRQ'][0]
   return CACHE.get(variant_name, None)


def plot(df, max_af_diff, output_filename):
   fig = plt.figure(figsize=(5, 5), dpi = 100)
   ax = fig.add_subplot(1, 1, 1)

   plt.grid(linestyle = '--', linewidth = 0.5)

   freqs_ok = df[df['AF_DIFF'] < max_af_diff]
   freqs_poor = df[df['AF_DIFF'] >= max_af_diff]

   ax.set_title('Alternate Allele Frequencies (AF)')
   ax.scatter(freqs_ok.AF, freqs_ok.TOPMED_AF, s = 5, label = r'$\Delta$ AF <' + f'{max_af_diff} (N={len(freqs_ok):,})')
   ax.scatter(freqs_poor.AF, freqs_poor.TOPMED_AF, s = 5, label = r'$\Delta$ AF $\geq$' + f'{max_af_diff} (N={len(freqs_poor):,})')

   ax.set_xlabel('Study AF')
   ax.set_ylabel('TOPMed AF')
   ax.legend()

   plt.savefig(output_filename)


if __name__ == '__main__':
   args = argparser.parse_args()

   df = pd.read_csv(args.in_freq, header = 0, sep = '\s+', low_memory = False, dtype = {'CHR': str})
   df['AF'] = 1.0 - df['MAF']

   df['TOPMED_AF'] = df['SNP'].apply(lambda x: get_topmed_freq(args.in_topmed, x))
   df = df[~df.TOPMED_AF.isnull()]
   df['AF_DIFF'] = abs(df['AF'] - df['TOPMED_AF'])

   plot(df, args.max_af_diff, f'{args.out_prefix}.png')

   df[df['AF_DIFF'] >= args.max_af_diff][['SNP']].to_csv(f'{args.out_prefix}.af_diff.txt', index = False, header = False, sep = '\t')

