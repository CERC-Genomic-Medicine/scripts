import argparse
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import chi2


### Inputs files  Trace Output compatible
argparser = argparse.ArgumentParser(description = 'Calculate Mahalanobis distance')
argparser.add_argument('-R', '--Ref', metavar = 'file', dest = 'RefPC', required = True, help = 'Reference PCA')
argparser.add_argument('-S', '--Study', metavar = 'file', dest = 'SPC', required = True, help = 'Study projected PCA')
argparser.add_argument('-PC', metavar = 'int', dest = 'pc', type = int, required = True, help = 'number of PC to consider')
argparser.add_argument('-O', '--Out', metavar = 'str', dest = 'Out', required = True, help = 'Name output')

if __name__ == '__main__':
	args = argparser.parse_args()
	#input only relevant Data
	pca_ref= pd.read_table(args.RefPC,header=0, sep='\s+',low_memory=False,usecols = ['indivID'] + [f'PC{i}' for i in range(1, args.pc+1)],index_col='indivID')
	pca_study= pd.read_table(args.SPC,header=0, sep='\s+',low_memory=False,usecols = ['indivID'] + [f'PC{i}' for i in range(1, args.pc+1)],index_col='indivID')
	## Calculate Mahalanobis Distance
	mean_ref = pca_ref.mean()
	pca_study_minus_mean = pca_study - mean_ref
	cov_ref = np.cov(pca_ref.values.T)
	inv_cov_ref = np.linalg.inv(cov_ref)
	left_term = np.dot(pca_study_minus_mean, inv_cov_ref)
	mahalanobis = np.dot(left_term, pca_study_minus_mean.T)
	dist_mahalanobis = mahalanobis.diagonal()
	## Calculate Oulier Probability
	pval = 1 - chi2.cdf(dist_mahalanobis, args.pc-1)
	df=pd.DataFrame(list(zip(dist_mahalanobis, pval)),columns=['Diag','Pval'],index=pca_study.index)
	df.to_csv(args.Out, sep= '\t', index=True)

