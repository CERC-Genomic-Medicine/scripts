#!/usr/bin/env python3
/*
* Author: Vincent Chapdelaine <vincent.chapdelaine@mail.mcgill.ca>
* Version: 1.0
* Year: 2023
*/

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
from scipy import stats
from ncls import NCLS

argparser = argparse.ArgumentParser(description = 'Script to plot either a Manhattan or Miami plot')
argparser.add_argument('-i1', metavar = 'name', dest = 'in_file1', type = str, required = True, help = 'Regenie Input to Mahattan/Miami plot')
argparser.add_argument('-T1', metavar = 'name', dest = 'in_file_Title1', type = str, required = True, help = 'Title of Manhattan plot/Subtitle of Miami plot')
argparser.add_argument('-T2', metavar = 'name', dest = 'in_file_Title2', type = str, required = False, help = 'Subtitle of Miami plot')
argparser.add_argument('-i2', metavar = 'name', dest = 'in_file2', type = str, required = False, help = 'Regenie Second Input to Miami plot.')
argparser.add_argument('-T', metavar = 'name', dest = 'trait', type = str, required = True, help = 'Name of the trait space = \'_\' | base Name of output ')


plt.rcParams.update({'font.size': 18})

def manhattan(Fig, file,inverted,title):
	df= pd.read_table(file, header=0,sep="\s+",low_memory=False,usecols=['CHROM','LOG10P','GENPOS']) ## necessary columns only
	df=df[df['LOG10P']!='TEST_FAIL'] ## removes any failed test
	df['minuslog10pvalue'] = df['LOG10P'].astype(float) ## type conversion
	df.loc[df.minuslog10pvalue>11,'minuslog10pvalue'] = 11 ## establish a limit for representation
	df["chromosome"]=df["CHROM"]
	df=df[~df.chromosome.isin(["Y"])] ## May be subject to change, Y chromosome was not considered
  ## Order and group chromosomes
	df.chromosome = df.chromosome.astype('category')
	df.chromosome = df.chromosome.cat.set_categories([ i for i in set(df.chromosome)], ordered=True)
	df = df.sort_values('chromosome')
  ## converts position per chrom to position total
	maxi=0
	for  i in df.chromosome.unique() :
		df.loc[df.chromosome==i,'GENPOS']=df.loc[df.chromosome==i,'GENPOS'] + maxi
		maxi=df.loc[df.chromosome==i,'GENPOS'].max()
	df['ind'] = df['GENPOS']
	df_grouped = df.groupby('chromosome')
	ax = Fig
	colors = ['black','gray']
	x_labels = []
	x_labels_pos = []
  #plot each chromosome
	for num, (name, group) in enumerate(df_grouped):
		group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=colors[name % len(colors)], s=9,ax=ax)
		x_labels.append(name)
		x_labels_pos.append((group['ind'].max() - (group['ind'].max() - group['ind'].min())/2))
	ax.set_xticks(x_labels_pos)
	ax.set_xticklabels(x_labels)
	# set axis limits
	ax.set_xlabel('Chromosome')
	ax.set_ylabel('-'+r'$\log_{10}$' +'(P-value)')
	ax.set_ylim([0, 12])
	if inverted :
		ax.invert_yaxis()
		ax.xaxis.tick_top()
		ax.xaxis.set_label_position('top')
		ax.title.set_text(title)
	else :
		ax.title.set_text(title)
	ax.axhline(y=7.30102999566 , color='r', linestyle='-') ##  GWAS significance in -log10 scale
	ax.title.set_text(title)


if __name__ == '__main__':
	args = argparser.parse_args()
	trait=args.trait
	if args.in_file2 is not None:
		fig, axs = plt.subplots(2, 1, constrained_layout=True,figsize=(32,15))
		manhattan(axs[0],args.in_file1,False, args.in_file_Title1)
		manhattan(axs[1],args.in_file2,True, args.in_file_Title2)
		plt.suptitle(trait.replace('_',' '),fontsize=22)
		fig.savefig(trait + '.png', dpi=300)
	if args.in_file2 is None:
		fig, axs = plt.subplots(1, 1, constrained_layout=True,figsize=(32,8))
		manhattan(axs,args.in_file1,False, args.in_file_Title1)
		fig.savefig(trait + '.png', dpi=300)
