#!/usr/bin/python3


'''
AUTHOR: Daniel Taliun
YEAR: 2022 ?

Modified by Vincent Chapdelaine (2023)

-> introduction of outlier filters (sd and iq)

'''

import pandas as pd
import math
import warnings
import argparse
import numpy as np
import sys

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


argparser = argparse.ArgumentParser(description = 'Prepares summary stats about continuous phenotypes in the CARTaGENE study.')
argparser.add_argument('-c', '--catalog', metavar = 'name', dest = 'in_catalog_spreadsheet', type = str, required = True, help = 'Excel spreadsheet "COMBINED_CATALOG_*.xlsx". It contains meta-information about all variables and the following sheets: "COMBINED_CATALOG", "Categories", "Linear data missing codes".')
argparser.add_argument('-p', '--phenotypes', metavar = 'name', dest = 'in_phenotypes', type = str, required = True, help = 'CSV file with phenotype values for each sample in the CARTaGENE study.')
argparser.add_argument('-s', '--samples', metavar = 'name', dest = 'in_samples', type = str, required = True, help = 'List of individual IDs to consider (e.g. PLINK\'s *.psam file). The individual IDs must be in the column named "IID".')
argparser.add_argument('-o', '--output', metavar = 'name', dest = 'out_prefix', type = str, required = True, help = 'Prefix for output files.') 
argparser.add_argument('-t', '--fold_outlier_threshold', metavar = 'value', dest = 'outlier_value', type = float, required = '--sd' in sys.argv or '--iq' in sys.argv, help = 'Fold value of the given outlier detection threshold (default 1.5)') 
argparser.add_argument('--sd', dest = 'sd', action='store_true', required = False, help = 'standard deviation outlier detection method, i.e. threshold value is x fold standard deviation above and below mean')
argparser.add_argument('--iq', dest = 'iq', action='store_true', required = False, help = 'interquartile outlier detection method, i.e. threshold value is x fold interquartile range below and above  quartile 1 and quartile 3 respectively')

catalog_sheet_categories = 'Categories'
catalog_sheet_missing_codes = 'Linear data missing codes'
catalog_sheet_main = 'COMBINED_CATALOG'

missing_categories = {'no answer', 'missing', 'missing answer', 'not available', 'not applicable', '-9', '-7', '77', '88', '99'}



def filter_binary_variables(df):
   for name, df_group in df.groupby(["SURVEY", "DOMAIN", "VARIABLE"]):
      df_only_codes = df_group[~df_group.CATEGORY.str.lower().isin(missing_categories)]
      if len(df_only_codes) != 2:
         continue
      yield name[2]


def filter_categorical_variables(df):
   for name, df_group in df.groupby(["SURVEY", "DOMAIN", "VARIABLE"]):
      df_only_codes = df_group[~df_group.CATEGORY.str.lower().isin(missing_categories)]
      if len(df_only_codes) < 2:
         continue
      yield name[2]


def recursive_detect_missing_codes_9(variable, missing_codes, df_pheno):
   has_missing_code = False
   max_value = df_pheno[variable].max()
   min_value = df_pheno[variable].min()
   if max_value in [9, 99, 999, 9999, 99999, 999999]:
      has_missing_code = True
      missing_codes.add(max_value)
   if min_value in [-9]:
      has_missing_code = True
      missing_codes.add(min_value)
   return has_missing_code


def recursive_detect_missing_codes_7(variable, missing_codes, df_pheno):
   has_missing_code = False
   max_value = df_pheno[variable].max()
   min_value = df_pheno[variable].min()
   if max_value in [7, 77, 777, 7777, 77777, 777777]:
      has_missing_code = True
      missing_codes.add(max_value)
   if min_value in [-7]:
      has_missing_code = True
      missing_codes.add(min_value)
   return has_missing_code



def filter_continuous_variables(df_variables, df_missing_codes, df_pheno, df_pheno_final):
   for index, row in df_variables.iterrows():
      variable = row['Varname']

      if variable not in df_pheno.columns:
         continue
      if variable in df_pheno_final.columns:
         continue
      if df_pheno.dtypes[variable] != 'float64':
         continue
   
      unit = row['UNIT_EN']
      missing_codes = variable_missing_codes.get(variable, set())
   
      df_pheno_cut = df_pheno[['SEX_BIRTH', variable]].copy()
      df_pheno_cut['RECODED'] = df_pheno_cut[variable].apply(lambda x: None if x in missing_codes else x )

      while recursive_detect_missing_codes_9('RECODED', missing_codes, df_pheno_cut):
         df_pheno_cut['RECODED'] = df_pheno_cut[variable].apply(lambda x: None if x in missing_codes else x )
      while recursive_detect_missing_codes_7('RECODED', missing_codes, df_pheno_cut):
         df_pheno_cut['RECODED'] = df_pheno_cut[variable].apply(lambda x: None if x in missing_codes else x )

      df_all = df_pheno_cut[~df_pheno_cut.RECODED.isna()]
      if df_all.empty:
         continue

      if args.iq:
         #interquartile method of filtering
         Q1,Q3 = np.percentile(df_all.RECODED , [25,75])
         IQR = Q3 - Q1
         ul = Q3+args.outlier_value*IQR
         ll = Q1-args.outlier_value*IQR
         df_pheno_cut['RECODED'] = df_pheno_cut[variable].apply(lambda x: None if (x > ul or x < ll) else x )
         n_outliers=sum([i > ul or i < ll for i in df_all.RECODED])

      if args.sd:
         #Standard deviation  method of filering
         mean = df_all.RECODED.mean()
         sd=np.std(df_all.RECODED)
         ul = mean+args.outlier_value*sd
         ll = mean-args.outlier_value*sd
         df_pheno_cut['RECODED'] = df_pheno_cut[variable].apply(lambda x: None if (x > ul or x < ll) else x )
         n_outliers=sum([i > ul or i < ll for i in df_all.RECODED])

      df_all = df_pheno_cut[~df_pheno_cut.RECODED.isna()]
      min_value = df_all.RECODED.min()
      max_value = df_all.RECODED.max()
      unique_values = len(df_all.RECODED.unique())

      if min_value == max_value:
         continue
      if unique_values <= 10:
         continue  
  
      df_males = df_all[df_all.SEX_BIRTH == 0]
      df_females = df_all[df_all.SEX_BIRTH == 1]
      mean_value = df_all.RECODED.mean()
      median_value = df_all.RECODED.median()
      value_counts = df_all.RECODED.value_counts()
      mode_freq = value_counts.values[0]
      mode_value = value_counts.index[0]
      n_total = len(df_all)
      n_males = len(df_males)
      n_females = len(df_females)

      df_pheno_final[variable] = df_pheno_cut.RECODED
      if args.sd or args.iq :
         yield {'DOMAIN': row['database'], 'VARIABLE': variable, 'UNIT': unit, 'N': n_total, 'MALES': n_males, 'FEMALES': n_females, 'MIN_VALUE': min_value, 'MAX_VALUE': max_value, 'UNIQUE_VALUES': unique_values, 'MEAN_VALUE': mean_value, 'MEDIAN_VALUE': median_value, 'MODE_VALUE': mode_value, 'MODE_FREQ': mode_freq, 'outliers' : n_outliers}
      else :
         yield {'DOMAIN': row['database'], 'VARIABLE': variable, 'UNIT': unit, 'N': n_total, 'MALES': n_males, 'FEMALES': n_females, 'MIN_VALUE': min_value, 'MAX_VALUE': max_value, 'UNIQUE_VALUES': unique_values, 'MEAN_VALUE': mean_value, 'MEDIAN_VALUE': median_value, 'MODE_VALUE': mode_value, 'MODE_FREQ': mode_freq}


if __name__ == '__main__':
   args = argparser.parse_args()
   df = pd.read_excel(args.in_catalog_spreadsheet, sheet_name = catalog_sheet_categories)
   df = df[["SURVEY", "DOMAIN", "VARIABLE", "CODE", "CATEGORY"]]
   df['CATEGORY'] = df['CATEGORY'].astype("string")
   
   binary_variables = set(filter_binary_variables(df))
   categorical_variables = set(filter_categorical_variables(df))

   df_missing_codes = pd.read_excel(args.in_catalog_spreadsheet, sheet_name = catalog_sheet_missing_codes)
   variable_missing_codes = {row['varname']: set(int(x) for x in row['Missing codes'].split(',')) for index, row in df_missing_codes.iterrows()}

   df_variables = pd.read_excel(args.in_catalog_spreadsheet, sheet_name = catalog_sheet_main, header = 0, skiprows = lambda x: x == 1, dtype = {'UNIT_EN': str})
   df_variables['UNIT_EN'] = df_variables['UNIT_EN'].fillna('') 

   df_variables = df_variables[(df_variables['Type '] == 1) & \
      (df_variables['Survey'].isin({'Phase A and B', 'Phase A', 'Phase B'})) & \
      (~df_variables['database'].isin({'ETHNIC', 'RESIDENTIAL_HISTORY', 'IDENTITY'})) & \
      (~df_variables['Format'].isin({'DATETIME', 'YYMMDD'})) & \
      (~df_variables['Varname'].isin(binary_variables)) & \
      (~df_variables['Varname'].isin(categorical_variables))]
   df_variables =  df_variables[~df_variables['Varname'].str.lower().str.endswith('_onset_year')]

   df_samples = pd.read_csv(args.in_samples, sep = '\t', header = 0, low_memory = False, usecols = ['IID'])
   df_pheno = pd.read_csv(args.in_phenotypes, sep = ',', header = 0, low_memory = False)
   df_pheno = df_pheno[df_pheno.PROJECT_CODE.isin(df_samples['IID'])].copy()
   assert len(df_pheno) == len(df_samples)

   df_pheno_final = pd.DataFrame({'FID': df_pheno.PROJECT_CODE, 'IID': df_pheno.PROJECT_CODE})
   df_variables_final = pd.DataFrame(filter_continuous_variables(df_variables, variable_missing_codes, df_pheno, df_pheno_final))
 
   df_pheno_final.to_csv(f'{args.out_prefix}.pheno', sep = '\t', header = True, index = False, na_rep = 'NA')
   df_variables_final.to_csv(f'{args.out_prefix}.summary.tsv', sep = '\t', header = True, index = False)
