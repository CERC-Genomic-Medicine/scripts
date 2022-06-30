import argparse
import pysam
import numpy as np
import scipy
import statsmodels.api as sm

argparser = argparse.ArgumentParser(description = 'Test for the significant differences in alternate allele frequencies between samples adjusting for PCs.')
argparser.add_argument('-v', '--vcf', metavar = 'file', dest = 'in_vcf', type = str, required = True, help = 'VCF with genotypes.')
argparser.add_argument('-l', '--labels', metavar = 'name', dest = 'in_labels', type = str, required = True, help = 'Tab-delimited text file with sample labels. Two columns (no header): <sample id>, <label name>.')
argparser.add_argument('-r', '--region', metavar = 'chr:start-stop', dest = 'in_region', type = str, required = True, help = '1-based region coordinates in the format: CHROM:START-STOP. If entire chromosome is needed, then specify just the chromosome name.')
argparser.add_argument('-p', '--pca', metavar = 'file', dest = 'in_pca', type = str, required = True, help = 'PCA file.')
argparser.add_argument('-k', '--k-pcs', metavar = 'number', dest = 'k_pcs', type = int, required = True, help = 'Number of PCs to use. E.g. 4')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_filename', type = str, required = True, help = 'Output filename.')


def load_pca(filename, k_pcs):
   sample2pca = dict()
   with open(filename, 'rt') as pca:
      line = pca.readline()
      header = line.strip().split()
      for line in pca:
         fields = dict(zip(header, line.strip().split()))
         sample_name = fields['IID']
         for k in range(1, k_pcs + 1):
            sample2pca.setdefault(sample_name, []).append(float(fields[f'PC{k}']))
   return sample2pca


def load_variables(filename):
   sample2variables = dict()
   with open(filename, 'rt') as variables:
      line = variables.readline()
      header = line.strip().split()
      for line in variables:
          fields = dict(zip(header, line.strip().split()))
          sample_name = fields['IID']
          sample_variables = [float(field_value) for field_name, field_value in fields.items() if field_name not in {'IID', 'FID'}]
          sample2variables[sample_name] = sample_variables
   return sample2variables


def test_variant(sample2gt, sample2pca, sample2label, k_pcs):
   result = {
      'AF_MODEL': None,
      'LRT_PVALUE': None,
      'FTEST_PVALUE': None
   }

   n_labels = len(set(sample2label.values()))
   n = [0] * n_labels
   an = [0] * n_labels

   y = []
   x = []
   x_null = []

   for sample_name, sample_gt in sample2gt.items():
      if sample_name not in sample2pca or sample_name not in sample2label:
         continue
      pcs = sample2pca[sample_name]
      label_code = sample2label[sample_name]
      n[label_code] += 2
      an[label_code] += sample_gt
      y.append([sample_gt, 2])
      dummy_coding = [0] * (n_labels - 1)
      if label_code > 0:
         dummy_coding[label_code - 1] = 1
      x.append(pcs + dummy_coding)
      x_null.append(pcs)

   af = [ a / b for a, b in zip(an, n) ]
   for i in range(0, n_labels):
      result[f'AF_{i}'] = af[i]
      result[f'N_{i}'] = n[i] / 2

   if sum(an) <= 1: # skip monomorphic and singletons
      return result

   x = sm.add_constant(x)
   x_null = sm.add_constant(x_null)

   model = sm.GLM(y, x, family=sm.families.Binomial())
   model_null = sm.GLM(y, x_null, family=sm.families.Binomial())
  
   try:
      results = model.fit()
      results_null = model_null.fit()
   except sm.tools.sm_exceptions.PerfectSeparationError as err:
      print("Warning:", err)
      return result

   #print(results.summary())
   #print(results_null.summary())

   h1 = ' = '.join( [f'x{i}' for i in range(k_pcs + 1, k_pcs + n_labels)] + ['0']  )
   f_test = results.f_test(h1)
   #print(h1, f_test)
   result['FTEST_PVALUE'] =  f_test.pvalue
    
   lrt = -2 * (results_null.llf - results.llf)
   lrt_pvalue = scipy.stats.chi2.sf(lrt, n_labels - 1) 
   b0 = results.params[0]
   p = np.exp(b0) / (1 + np.exp(b0))

   result['LRT_PVALUE'] =  lrt_pvalue
   result['AF_MODEL'] = p

   return result
   


if __name__ == '__main__':
   args = argparser.parse_args()

   print(f'Region: {args.in_region}')

   # Load sample names and labels
   label2code = dict()
   code2label = dict()
   sample2label = dict()
   label2samples = dict()
   with open(args.in_labels, 'rt') as ifile:
      for line in ifile:
         sample_name, label = line.rstrip().split('\t')
         label_code = label2code.get(label, None)
         if label_code is None:
            label_code = len(label2code)
            label2code[label] = label_code
            code2label[label_code] = label
         sample2label[sample_name] = label_code
         
   print('\nCoding:')
   for label, label_code in label2code.items():
      print(f'{label_code}: {label}')
    
   with pysam.VariantFile(args.in_vcf, 'r') as ivcf:
      vcf_samples = set(ivcf.header.samples)
      samples_without_gt = []
      for sample_name, label in sample2label.items():
          if sample_name not in vcf_samples:
              samples_without_gt.append(sample_name)
      for sample_name in samples_without_gt:
         sample2label.pop(sample_name)

   for sample_name, label_code in sample2label.items():
      label2samples.setdefault(label_code, []).append(sample_name)

   print('\nSamples with GT:')
   for label_code, samples in label2samples.items():
      print(f'{label_code}: N={len(samples):,}')

           
   # Load PCs
   sample2pca = load_pca(args.in_pca, args.k_pcs)

   print('\nAnalyzing... ')
   header = ['VARIANT']
   for i in range(0, len(code2label)):
      header.append(f'AF_{code2label[i]}')
      header.append(f'N_{code2label[i]}')
   header += ['AF_MODEL', 'LRT_PVALUE', 'FTEST_PVALUE']

   with open(args.out_filename, 'wt', buffering=1) as ofile, pysam.VariantFile(args.in_vcf, 'r') as ivcf:
      ofile.write('{}\n'.format('\t'.join(header)))
      for record in ivcf.fetch(region = args.in_region):
         if len(record.alts) > 1:
            continue
         variant_name = f'{record.chrom}:{record.pos}:{record.ref}:{record.alts[0]}'
         sample2gt = dict()
         for sample, sample_frmt in record.samples.iteritems():
            if None in sample_frmt['GT']:
               continue
            sample2gt[sample] = sum(sample_frmt['GT'])
         r = test_variant(sample2gt, sample2pca, sample2label, args.k_pcs)
         ofile.write(f'{variant_name}\t')
         for i in range(0, len(code2label)):
            ofile.write('{}\t{:.0f}\t'.format(r[f'AF_{i}'], r[f'N_{i}']))
         ofile.write(f'{r["AF_MODEL"]}\t{r["LRT_PVALUE"]}\t{r["FTEST_PVALUE"]}\n')
   print('Done.')
