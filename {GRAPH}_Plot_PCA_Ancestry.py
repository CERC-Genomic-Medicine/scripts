#!/usr/bin/env python3


import matplotlib.colors as mcolors
import pandas as pd
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Ancestry_projection')

parser.add_argument('-P', '--Projected', metavar = 'File',dest='Proj',required = True, type=str, help='Projected samples position')
parser.add_argument('-R','--Reference', metavar = 'File',dest='Ref',required = True, type=str, help='Reference samples position')
parser.add_argument('-S', '--Study' , metavar = 'File',dest='Study',required = True, type=str, help='Infered Ethnicity projected samples')
parser.add_argument('-A', '--Ancestry', metavar = 'File',dest='Ancestry', type=str,required = True, help='Reference Samples ethnicity')
parser.add_argument('-c', '--selected', metavar = 'str',dest='Seleted', type=str,required = True, help='Ethnicty to plot')
parser.add_argument('-T', '--Threshold', metavar = 'number',dest='Threshold', type=float,required = True, help='Reference Samples ethnicity')
parser.add_argument('-n', '--PC', metavar = 'number',dest='n', type=int,required = True, help='number of PC to plot')
parser.add_argument('-l', '--label', metavar = 'string',dest='study_name', type=str,required = False, default='Study', help='label of the study')
parser.add_argument('--out', dest='output', default='output', type=str, help='output')

plt.rcParams.update({'font.size': 12})

def Name_dictionnary(Ancestries):
  print(Ancestries)
  NameDict_default={'CSA': 'Central South Asian','SAS': 'South Asian', 'EUR': 'European', 'EAS': 'East Asian', 'AMR': 'Ad Mixed American', 'AFR': 'African', 'MID' : 'Middle eastern', 'OCE' : 'Oceanian'}
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
  Ax.scatter(reference['PC'+str(PCa)],reference["PC"+str(PCb)],alpha=0.5,facecolors='none', edgecolors=[color_dict[dic[i]] for i in list(reference.index)])
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
  Infered = pd.read_csv(args.Study, header=0,sep="\s+",low_memory=False,usecols=['indivID',"predicted_ancestry",'max_prob'],index_col='indivID')
  ID_ploted = Infered.loc[(Infered['max_prob']>=args.Threshold) & (Infered['predicted_ancestry']==args.Seleted)].index
  Projected = Projected_whole.loc[ID_ploted,:]
  Ancestry_Ref = pd.read_csv(args.Ancestry,header=0, index_col=0)
  NameDict=Name_dictionnary(set(Ancestry_Ref.genetic_region))
  color_dict=Color_dict(NameDict)
  dic=dict(zip(Ancestry_Ref.index,Ancestry_Ref.iloc[:, 0]))
  grid=find_factorial_grid_approx_prime(suite(args.n))
  fig, axs = plt.subplots(grid[0], grid[1], figsize=(grid[0]*8+5,grid[0]*8), dpi = 100*args.n*2)
  i=0
  for x in range(1, args.n):
    for y in range(x+1, args.n+1):
        plot_Projection(axs[i//grid[1],i%grid[1]],Projected,Ref_proj,x,y)
        i=i+1
patches=[]
for q in color_dict.keys() :
  patches.append(mpatches.Patch(color=color_dict[q], label=NameDict[q]))
patches.append(mpatches.Patch(color='black', label=f"{args.study_name} of {NameDict[args.Seleted]} predicted Ancestry \n (Probability ≥ {args.Threshold}, N = {len(Projected.index)})"))
plt.legend(handles=patches,title= "Genetic Ancestry",loc='center left', bbox_to_anchor=(1.0, 0.5))
fig.suptitle("Genetic ancestry projection", fontsize=16)
fig.tight_layout()
plt.savefig(args.output + '.png')
