#!/usr/bin/env python3


import matplotlib.colors as mcolors
import pandas as pd
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Ancestry_projection')

parser.add_argument('-P', '--projected', metavar='FILE', dest='Proj', required=True, type=str, help='Projected study data PCA coordinates in ProPC.coord file from TRACE output.')
parser.add_argument('-R','--reference', metavar='FILE', dest='Ref', required=True, type=str, help='Reference data PCA coordinates in RefPC.coord file from TRACE output.')
parser.add_argument('-S', '--study-ancestry' , metavar='FILE', dest='Study', required=True, type=str, help='Predicted ancestry labels for study samples in tab-separated file (predicted_ancestry.txt).')
parser.add_argument('-A', '--ref-ancestry', metavar='FILE', dest='Ancestry', type=str, required=True, help='Population labels for reference data.')
parser.add_argument('-c', '--selected', dest='Selected', type=str, required=True, help='Predicted ancestry group in study data to subset and plot.')
parser.add_argument('-T', '--threshold', dest='Threshold', type=float,required=True, help='Threshold of probability to use when subsetting by predicted ancestry.')
parser.add_argument('-n', '--nPC', dest='n', type=int, required=True, help='Number of PCs to plot.')
parser.add_argument('-l', '--label', metavar='STUDY', dest='study_name', type=str, required=False, default='Study', help='Study label to use when titling plot and legend.')
parser.add_argument('--out', dest='output', default='output', type=str, help='Prefix of output image file name (output.png).')

plt.rcParams.update({'font.size': 14})

def Name_dictionary(Ancestries):
  print(Ancestries)
  #NameDict_default={'CSA': 'Central & South Asia','SAS': 'South Asia', 'EUR': 'Europe', 'EAS': 'East Asia', 'AMR': 'Americas', 'AFR': 'Africa', 'MID' : 'Middle East', 'OCE' : 'Oceania'}
  NameDict_default={'CSA': 'Reference: Central & South Asia (CSA)','SAS': 'Reference: South Asia (SAS)', 'EUR': 'Reference: Europe (EUR)', 'EAS': 'Reference: East Asia (EAS)', 'AMR': 'Reference: Americas (AMR)', 'AFR': 'Reference: Africa (AFR)', 'MID' : 'Reference: Middle East (MID)', 'OCE' : 'Reference: Oceania (OCE)'}
  if Ancestries.issubset(set(NameDict_default.keys())) :
    return {key: NameDict_default[key] for key in Ancestries if key in NameDict_default}
  else :
    NameDict=dict(zip(Ancestries,Ancestries))
    return NameDict

def Color_dict(NameDict):
  colors=mcolors.TABLEAU_COLORS
  colors=list(colors.keys())[slice(0,len(NameDict.keys()),1)]
  color_dict={list(set(NameDict.keys()))[i]: colors[i] for i in range(0,len(set(NameDict.keys())))}
  return color_dict


def plot_Projection(Ax,study,reference,PCa,PCb) :
  Ax.scatter(reference['PC'+str(PCa)],reference["PC"+str(PCb)],alpha=0.4,s=5,facecolors='none', edgecolors=[color_dict[dic[i]] for i in list(reference.index)])
  Ax.scatter(study["PC"+str(PCa)],study["PC"+str(PCb)],s=5,facecolor='none',edgecolor='black')
  Ax.set(adjustable='box')
  Ax.set_xlabel('PC'+str(PCa))
  Ax.set_ylabel('PC'+str(PCb))

def suite(n):
    summ=0
    for i in range(0,n):
        summ=summ+i
    return summ

def is_prime(n):
    """Check if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def find_factorial_grid_approx_prime(product):
    # If the product is a prime number greater than 3, approximate to the next value
    if is_prime(product) and product > 3:
        product += 1  # Approximate to the next number        
    
    a, b = 1, product
    
    # Loop to find the factors of product that are closest to each other
    for i in range(1, int(product ** 0.5) + 1):
        if product % i == 0:
            a, b = i, product // i
            
    return a, b

if __name__ == '__main__':
  args = parser.parse_args()
  Ref_proj=pd.read_csv(args.Ref, header=0,sep="\s+",low_memory=False,usecols=['indivID']+[f"PC{i}" for i in range(1, args.n + 1)],index_col='indivID')
  Projected_whole = pd.read_csv(args.Proj, header=0,sep="\s+",low_memory=False,usecols=['indivID']+[f"PC{i}" for i in range(1, args.n + 1)],index_col='indivID')
  Inferred = pd.read_csv(args.Study, header=0,sep="\s+",low_memory=False,usecols=['indivID',"predicted_ancestry",'max_prob'],index_col='indivID')
  ID_plotted = Inferred.loc[(Inferred['max_prob']>=args.Threshold) & (Inferred['predicted_ancestry']==args.Selected)].index
  Projected = Projected_whole.loc[ID_plotted,:]
  Ancestry_Ref = pd.read_csv(args.Ancestry,header=0, index_col=0)
  NameDict=Name_dictionary(set(Ancestry_Ref.genetic_region))
  color_dict=Color_dict(NameDict)
  dic=dict(zip(Ancestry_Ref.index,Ancestry_Ref.iloc[:, 0]))
  grid=find_factorial_grid_approx_prime(suite(args.n))
  print(grid)
  if grid[0] ==1 :
    if grid[1] !=1:
      fig, axs = plt.subplots(grid[0],grid[1], figsize=(grid[1]*8,grid[0]*8), dpi = 300)
      i=0
      for x in range(1, args.n):
        for y in range(x+1, args.n+1):
          plot_Projection(axs[i],Projected,Ref_proj,x,y)
          i=i+1
    else :
      fig, axs = plt.subplots(grid[0],grid[1], figsize=(grid[1]*8,grid[0]*8), dpi = 300)
      plot_Projection(axs,Projected,Ref_proj,1,2)
  else:
    fig, axs = plt.subplots(grid[0], grid[1], figsize=(grid[1]*8,grid[0]*8), dpi = 300)
    i=0
    for x in range(1, args.n):
      for y in range(x+1, args.n+1):
        plot_Projection(axs[i//grid[1],i%grid[1]],Projected,Ref_proj,x,y)
        i=i+1

patches=[]
for q in color_dict.keys() :
  patches.append(mpatches.Patch(color=color_dict[q], label=NameDict[q]))
patches.append(mpatches.Patch(color='black', label=f"{args.study_name} predicted {args.Selected} samples \n (Probability ≥ {args.Threshold}, N = {len(Projected.index)})"))
plt.legend(handles=patches,title= "Genetic Ancestry",loc='center left', bbox_to_anchor=(1.0, 0.5),fontsize="18")
fig.suptitle(f"Projection of {args.study_name} {args.Selected} samples (N={len(Projected.index)}) against reference PCA", fontsize=20)
fig.tight_layout()
plt.savefig(args.output + '.png')
