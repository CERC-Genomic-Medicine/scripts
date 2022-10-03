import argparse
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import chi2

argparser = argparse.ArgumentParser(description = 'Calculate Mahalanobis distance and associated outlier probability')
argparser.add_argument('-R', '--Ref', metavar = 'file', dest = 'RefPC', required = True, help = 'Reference Sample PCA : minimally an ID and a number of PC concordant with PC argument ( considered fields indivID and PCn tab or space delimited)')
argparser.add_argument('-S', '--Study', metavar = 'file', dest = 'SPC', required = True, help = 'Projected Study Samples : minimally an ID and a number of PC concordant with PC argument ( considered fields indivID and PCn tab or space delimited) ')
argparser.add_argument('-E', '--Ethnicity', metavar = 'file', dest = 'Ethn', required = True, help = 'Ethnicity : minimally an ID and known ancestry (considered fields indivID and Ethnicity tab or space delimited)')
argparser.add_argument('-PC', metavar = 'int', dest = 'pc', type = int, required = True, help = 'number of PC to consider')
argparser.add_argument('-T','--Threshold', metavar = 'float', dest = 'thresh', type = float, required = False, help = 'threshold of the pvalue derived of mahalanobis distance for rejection from an ancestry')
argparser.add_argument('-O', '--Out', metavar = 'str', dest = 'Out', required = True, help = 'Output file name')

### Output file Description
# o) tab delimited dataframe with 2 field per ancestry plus index
# o) index correspond to one Study individuals in -S argument
# o) Per ancestry : Ancestry_Maha -> Distance from ancestry ||  distribution Ancestry_Pval -> Outlier p-value

def Maha(pca_ref,pca_study):
	mean_ref = pca_ref.mean()
	pca_study_minus_mean = pca_study - mean_ref # Distance verctor from the mean
	cov_ref = np.cov(pca_ref.values.T) # covariance matrix
	inv_cov_ref = np.linalg.inv(cov_ref) #inverse covariance matrix
	left_term = np.dot(pca_study_minus_mean, inv_cov_ref) #inverse covariance matrix * Distance verctor from the mean
	mahalanobis = np.dot(left_term, pca_study_minus_mean.T) # multiplied again by Distance verctor from the mean
	MAHA = mahalanobis.diagonal()
	## Calculate Oulier Probability
	PVALUE = 1 - chi2.cdf(MAHA, args.pc-1)
	df=pd.DataFrame(list(zip(MAHA, PVALUE)),columns=['Maha','Pval'],index=pca_study.index)
	return df

if __name__ == '__main__':
	args = argparser.parse_args()
	#input only relevant Data
	ref= pd.read_table(args.RefPC,header=0, sep='\s+',low_memory=False,usecols = ['indivID'] + [f'PC{i}' for i in range(1, args.pc+1)],index_col='indivID')
	Eth= pd.read_table(args.Ethn,header=0, sep='\s+',low_memory=False,usecols = ['indivID', 'Ethnicity'],index_col='indivID')
	study= pd.read_table(args.SPC,header=0, sep='\s+',low_memory=False,usecols = ['indivID'] + [f'PC{i}' for i in range(1, args.pc+1)],index_col='indivID')
	res=pd.DataFrame(index=study.index)
	#Calulate Mahalanobis for each ancestry seperatly
	for i in Eth.Ethnicity.unique():
		df=Maha(ref.loc[Eth[Eth.Ethnicity==i].index],study)
		df.columns=[i+'_'+j for j in df.columns]
		res=pd.concat([res,df],axis=1,copy=False)
	#Determine list of non rejected Ancestry
	pvalues=res[[i for i in res.columns if 'Pval' in i]] 
	if args.thresh :
		threshold=args.thresh
	else : 
		threshold=0.001 
    	Eth_list=[list(pvalues.loc[i][pvalues.loc[i]>threshold].index) for i in pvalues.index]
	Eth_list=[['None'] if not i else i for i in Eth_list]
	res['Ethnicity']=[','.join(q).replace('_Pval', '') for q in Eth_list]
	res.to_csv(args.Out, sep= '\t', index=True)
