import pandas as pd
import argparse

argparser = argparse.ArgumentParser(description = 'Prepares summary stats about binary phenotypes in the CARTaGENE study.')
argparser.add_argument('-c', '--catalog', metavar = 'name', dest = 'in_catalog_spreadsheet', type = str, required = True, help = 'Excel spreadsheet "COMBINED_CATALOG_*.xlsx". It contains meta-information about all variables and the following sheets: "COMBINED_CATALOG", "Categories", "Linear data missing codes".')
argparser.add_argument('-p', '--phenotypes', metavar = 'name', dest = 'in_phenotypes', type = str, required = True, help = 'CSV file with phenotype values for each sample in the CARTaGENE study.')
argparser.add_argument('-s', '--samples', metavar = 'name', dest = 'in_samples', type = str, required = True, help = 'List of individual IDs to consider (e.g. PLINK\'s *.psam file). The individual IDs must be in the column named "IID".')
argparser.add_argument('-o', '--output', metavar = 'name', dest = 'out_prefix', type = str, required = True, help = 'Prefix for output files.')


catalog_sheet = 'Categories'
categories = {'yes', 'no', 'missing', 'no answer'}


def filter_binary_variables(df):
   for name, group in df.groupby(["SURVEY", "DOMAIN", "VARIABLE"]):
      if len(group) != 4: # Skip if not 4 values
         continue
      if not all(str(c).lower() in categories for c in group.CATEGORY):
         continue
      yes_code = group[group.CATEGORY.str.lower() == 'yes'].CODE.values[0]
      no_code = group[group.CATEGORY.str.lower() == 'no'].CODE.values[0]
      noanswer_code = group[group.CATEGORY.str.lower() == 'no answer'].CODE.values[0]
      missing_code = group[group.CATEGORY.str.lower() == 'missing'].CODE.values[0]
      yield {'SURVEY': name[0], 'DOMAIN': name[1], 'VARIABLE': name[2], 'CODE_YES': yes_code, 'CODE_NO': no_code, 'CODE_NOANSWER': noanswer_code, 'CODE_MISSING': missing_code}


def lookup_cases_controls(df_variables, df_pheno, df_pheno_final):
   for index, row in df_variables.iterrows():
      variable = row.VARIABLE
      if variable not in df_pheno.columns:
         continue
      if variable in df_pheno_final.columns:
         continue

      codes_map = {
         row.CODE_YES: 1,
         row.CODE_NO: 0
      }

      df_pheno_cut = df_pheno[['SEX_BIRTH', variable]].copy()
      df_pheno_cut['RECODED'] = df_pheno_cut[variable].apply(lambda x: codes_map.get(x, None))
      df_pheno_final[variable] = df_pheno_cut.RECODED
     
      df_all = df_pheno_cut[df_pheno_cut.RECODED.isin({0, 1})]
      df_males = df_all[df_all.SEX_BIRTH == 0]
      df_females = df_all[df_all.SEX_BIRTH == 1]

      n_total = len(df_all)
      n_males = len(df_males)
      n_females = len(df_females)
      
      n_cases = len(df_all[df_all.RECODED == 1])
      n_cases_males = len(df_males[df_males.RECODED == 1])
      n_cases_females = len(df_females[df_females.RECODED == 1])

      n_controls = len(df_all[df_all.RECODED == 0])

      yield({'DOMAIN': row.DOMAIN, 'VARIABLE': variable, 'N': n_total, 'MALES': n_males, 'FEMALES': n_females, 'CASES': n_cases, 'CONTROLS': n_controls, 'CASES_MALES': n_cases_males, 'CASES_FEMALES': n_cases_females})



if __name__ == '__main__':
   args = argparser.parse_args()
   df = pd.read_excel(args.in_catalog_spreadsheet, sheet_name = catalog_sheet)

   df_codes_only = df[["SURVEY", "DOMAIN", "VARIABLE", "CODE", "CATEGORY"]]
   df_variables = pd.DataFrame(filter_binary_variables(df_codes_only))

   df_samples = pd.read_csv(args.in_samples, sep = '\t', header = 0, low_memory = False, usecols = ['IID'])
   df_pheno = pd.read_csv(args.in_phenotypes, sep = ',', header = 0, low_memory = False)
   df_pheno = df_pheno[df_pheno.PROJECT_CODE.isin(df_samples['IID'])].copy()

   assert len(df_pheno) == len(df_samples)

   df_pheno_final = pd.DataFrame({'FID': df_pheno.PROJECT_CODE, 'IID': df_pheno.PROJECT_CODE})
   df_variables_final = pd.DataFrame(lookup_cases_controls(df_variables, df_pheno, df_pheno_final))

   df_pheno_final.to_csv(f'{args.out_prefix}.pheno', sep = '\t', header = True, index = False, na_rep = 'NA', float_format = '%.0f')
   df_variables_final.to_csv(f'{args.out_prefix}.summary.tsv', sep = '\t', header = True, index = False)
