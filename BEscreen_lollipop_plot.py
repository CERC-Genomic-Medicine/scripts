#!/usr/bin/env python3


#
# Author : Vincent Chapdelaine (vincent.chapdelaine@mcgill.ca)
#
# This scripts includes functions to :
#     - Create a representation of a protein or gene with feature/domains
#     - A methodologie to collapse intron (optional)
#     - Creat a lolliplot
#
# The main goal is to represent two results from Base Editing Crispr screening. (as pos/neg lolliplot separated by a diagram of the protein)
#
# Example:
#   python3 .py -b path/to/bedfile.bed -r path/to/crispr_screen.tsv --highlight_region Tetramerization_Domain -m Full --out Lolliplot -c 100
#

import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.lines as mlines
from matplotlib.patches import Patch
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import argparse
from scipy.stats import permutation_test
from scipy.stats import wilcoxon
from scipy.stats import rankdata
from scipy.stats import binom_test


parser = argparse.ArgumentParser(description='A plotting script for loolipop plot')
parser.add_argument('-b', '--bed', metavar='FILE', dest='bed_file', required=True, type=str, help='space/tab-delimited bed file of features, must contains columns : start, stop and name')
parser.add_argument('-t', '--target', metavar='FILE', dest='result_file', required=True, type=str, help='space/tab-delimited result file TBD')
parser.add_argument('-v', '--vep', metavar='FILE', dest='vep_file', required=True, type=str, help='tab delimited variant effect prediction file (VEP)')
parser.add_argument('-n', '--negative_controls', metavar='FILE', dest='negative_control_file', required=False, type=str, help='space/tab-delimited result file TBD')
parser.add_argument('-c', '--intron_collapse', metavar='number', dest='collapse_factor', required=False,default=1 , type=int, help='factor by which to collapse introns (default uncollapsed), intron are defined by name containing "intron"')
parser.add_argument('--stat_method',metavar='str', dest='stat_method', required=False,default='Rank' , choices={"binom_sign", "Rank",'sign_test'}, type=str, help='Statistical method to calculate p-value of Biological significance')
parser.add_argument('--Bio_threshold', dest='Biological_threshold', metavar='float',type=float, required=False, default='0.05', help='Biological threshold (two-sided)')
parser.add_argument('--scheme_location', metavar='str', dest='scheme_loc', required=False,default='middle' , type=str, choices={"top",'bottom', "middle"}, help='Location of the protein/gene scheme ("top"/"middle"/"bottom") (default : middle)')
parser.add_argument('--histogram', dest='histogram', required=False,action='store_true', help='Flag for wheather a histogram should represent the coverage')
parser.add_argument('-p', '--Prob_Threshold', metavar='float',type=float, dest='p_thresh', required=False, default='1', help='P-Value Threshold for loolipop representation')
parser.add_argument('-B', '--Biological_Sig', dest='Biological_REP', action='store_true', help='Flag for Biological significance requirerement in lolipop plot')
parser.add_argument('--no_stem', dest='no_stem', action='store_true', help='Flag to remove the stemlines')
parser.add_argument('-m', '--legend_mode', metavar='str', dest='legend_mode', required=False,default='Mix' , choices={"Full", "Mix"}, type=str, help='legend mode (possible values Mix (default) or Full), full all features are in legend, mix of on-graph and on-legend')
parser.add_argument('--highlight_region', metavar='str', dest='highlight', required=False,default=None , type=str, help='Feature to be highlighted')
parser.add_argument('--out', dest='output', default='output', type=str, help='Prefix of output image file name (output.png).')

# Constante

#def wilcoxon_test(sample, hypothesized_median=0):
#    # Calculate the differences from the hypothesized median
#    differences = [x - hypothesized_median for x in sample]
#    
#    # Perform the Wilcoxon signed-rank test
#    stat, p_value = wilcoxon(differences,alternative='two-sided')
#    
#    return p_value


# Mapping for markers and colors based on Consequence
consequence_mapping = {
    'synonymous': ('g', 'D'),  # green, diamond
    'missense': ('purple', 'o'),  # purple, circle
    'non-sense': ('red', 's'),  # red, square
    'splice': ('gold', '^')  # gold, triangle
}
variant_consequences_mapping = {
'missense_variant': 'missense',
'intron_variant': 'none',
'downstream_gene_variant': 'none',
'NMD_transcript_variant': 'nonsense',
'upstream_gene_variant': 'none',
'3_prime_UTR_variant': 'none',
'synonymous_variant': 'synonymous',
'non_coding_transcript_exon_variant': 'none',
'splice_region_variant': 'splice',
'splice_polypyrimidine_tract_variant': 'splice',
'stop_gained': 'nonsense',
'coding_sequence_variant': 'none',
'5_prime_UTR_variant': 'none',
'regulatory_region_variant': 'none',
'splice_donor_variant': 'splice',
'splice_acceptor_variant': 'splice',
'non_coding_transcript_variant': 'none',
'splice_donor_region_variant': 'splice',
'splice_donor_5th_base_variant': 'splice',
'TF_binding_site_variant': 'none',
'start_lost': 'nonsense',
'incomplete_terminal_codon_variant': 'nonsense'
}

def convert_to_int(s):
    try:
        return int(s)
    except ValueError:
        return None

def Transform_MaGeCK(filename):
    Mage = pd.read_csv(filename,sep='\t',header=0)
    returned = pd.DataFrame({
    'id': Mage.id,
    'LFC': [None] * len(Mage.id),
    'p_value': [None] * len(Mage.id)
    })
    returned['LFC'] = Mage['pos|lfc']
    returned['p_value'] = [row['neg|p-value'] if row['pos|lfc']<0 else row['pos|p-value'] for _, row in Mage.iterrows()]
    return returned


def parse_VEP(filename, variant_consequences_mapping):
    # From Oligomer.py https://github.com/CERC-Genomic-Medicine/CRISPR_Library_prep.git
    with open(filename) as f:
        for line in f:
            if line[0:2]=='##':
                continue
            elif line[0:1]=='#':
                header = line.strip('\n').split("\t")
                colConsequence = header.index('Consequence')
                colAA = header.index('Protein_position')
            else :
                if 'header' not in locals():
                    raise Exception("VEP format not as expected") 
                else :
                    record = line.strip('\n').split("\t")
                    ID = record[0]
                    POSi = record[colAA].split('-')[0]
                    POS = convert_to_int(POSi)
                    if 'non-sense' in [variant_consequences_mapping[i] for i in record[colConsequence].split(',')]:
                        consequence = 'non-sense'
                    elif 'splice' in [variant_consequences_mapping[i] for i in record[colConsequence].split(',')]:
                        consequence = 'splice'
                    elif 'missense' in [variant_consequences_mapping[i] for i in record[colConsequence].split(',')]:
                        consequence = 'missense'
                    elif 'synonymous' in [variant_consequences_mapping[i] for i in record[colConsequence].split(',')]:
                        consequence = 'synonymous'
                    else :
                        consequence = None
                    yield ID, POS, consequence

# Function to generate a random color
def is_valid_color(color):
    """
    Determine if color is a valid value.

    Parameters:
    color (str) : color to be validated

    Returns:
    bool : validity of the color (in matplotlib)
    """
    try:
        mcolors.to_rgba(color)
        return True
    except ValueError:
        return False


def calc_biological_significance(value, distribution,method,thresh):
    if method == 'Rank':
        MAX = distribution.quantile(q=1-thresh)
        MIN = distribution.quantile(q=thresh)
        sign=value<=MIN or value>=MAX
    if method == 'binom_test':
        # Count how many values are greater than or equal to the test value
        greater_equal_count = sum(v >= value for v in distribution)
        n = len(distribution)
        
        # Perform the binom test
        p_value = binom_test(greater_equal_count, n, 0.5, alternative='two-sided')
        sign= thresh>p_value
    
    elif method == 'sign_test':
        # Count positive and negative differences
        positive_differences = sum((v - value) > 0 for v in distribution)
        negative_differences = sum((v - value) < 0 for v in distribution)
        
        # Perform a binomial test for the sign test
        n = positive_differences + negative_differences
        p_value = binom_test(positive_differences, n, 0.5, alternative='two-sided')
        sign= thresh>p_value
    return sign


# Function to get color palettes
def get_color_palette(num_entries):
    """
    Determine the appropriate palette for the number of entries.

    Parameters:
    num_entries (int) : number of entries
    
    Returns:
    str: color palette
    """
    pastel1 = plt.colormaps.get_cmap('Pastel1').colors
    tab20 = plt.colormaps.get_cmap('tab20').colors
    if num_entries <= len(pastel1):
        return pastel1
    elif num_entries <= len(tab20):
        return tab20
    else:
        print("Warning: Number of entries exceeds the available colors in tab20 palette.")
        return tab20

def get_color(index, palette, existing_colors):
    """
    Generate a color from the palette that does not reoccur.
    
    Parameters:
    index (int): Index of the current entry.
    palette (list): List of colors in the palette.
    existing_colors (set): Set of existing colors to avoid duplicates.
    
    Returns:
    str: Hex color string.
    """
    color = mcolors.to_hex(palette[index % len(palette)])
    while color in existing_colors:
        index += 1
        color = mcolors.to_hex(palette[index % len(palette)])
    existing_colors.add(color)
    return color


def create_color_dict(df):
    """
    Create a dictionary of {name: color} for all names in the BED file
    that do not already have a valid color and skip the row if 'intron' is within the name of the row.
    
    Parameters:
    df (pd.DataFrame): BED DataFrame with at least 'name' and optionally 'color' columns.
    
    Returns:
    dict: Dictionary mapping region names to their colors.
    """
    if 'color' in df.columns:
        existing_colors = {color for color in df.get('color', []) if pd.notnull(color) and is_valid_color(color)}
        filtered_df = df[~df['name'].str.contains('intron') & (~df['color'].apply(is_valid_color) | df['color'].isnull())]
    else :
        existing_colors=set()
        filtered_df = df[~df['name'].str.contains('intron')]

    
    # Filter out rows with 'intron' in the name and with valid colors
    num_entries = len(filtered_df)

    color_palette = get_color_palette(num_entries)
    color_dict = {}

    index = 0
    for _, row in df.iterrows():
        if 'intron' in row['name']:
            continue
        
        if 'color' in df.columns and pd.notnull(row['color']) and is_valid_color(row['color']):
            color_dict[row['name']] = row['color']
        else:
            color_dict[row['name']] = get_color(index, color_palette, existing_colors)
            index += 1

    return color_dict

def adjust_position(pos, df, b):
    """
    Adjust the position based on the BED file regions and a condensation factor.
    
    Parameters:
    pos (int): Original position.
    bed (pd.DataFrame): BED DataFrame.
    b (int): condensation factor.
    
    Returns:
    int: Adjusted position.
    """
    c = 1 - (1 / b)
    cumulative_intron_length = 0
    for _, row in df.iterrows():
        start, end = row['start'], row['end']
        if 'intron' in row['name']:
            if start < pos:
                intron_length = min(end, pos) - start
                cumulative_intron_length += intron_length
            else:
                break
    adjusted_position = pos - c * cumulative_intron_length
    return adjusted_position

def reverse_adjust_position(pos, df, b):
    """
    Reverse the adjustement of a position based on the BED file regions and a condensation factor.
    
    Parameters:
    pos (int): Original position.
    bed (pd.DataFrame): BED DataFrame.
    b (int): condensation factor.
    
    Returns:
    int: reverse adjusted position.
    """
    cumulative_intron_length = 0
    for _, row in df.iterrows():
        start, end = row['start'], row['end']
        adjusted_start, adjusted_end = adjust_position(row['start'],df,b), adjust_position(row['end'],df,b)
        if 'intron' in row['name']:
            if adjusted_start < pos :
                intron_length = min(adjusted_end,pos) - adjusted_start
                cumulative_intron_length += intron_length
            else:
                break
    original_position = pos + (b-1) * cumulative_intron_length
    return original_position

def add_text_or_legend(ax, start, length,size, text, color, legend_dict, legend_mode):
    """
    Add legend hande to dictionnary for plotting
    
    Parameters:
    ax (matplotlib.axes._subplots.AxesSubplot): Axis to add the legend to.
    start (int) : start position of feature
    length (int) : length of the feature
    size (float) : Width of the rendered feature
    text (str) : Feature name
    color (str) : Color of the feature
    legend_dict (dict) : dictionnary to be updated
    legend_mode (str) : full -> all feature in legend, mix -> Feature with short name directly on figure
    
    Returns:
    None
    """
    renderer = ax.figure.canvas.get_renderer()
    
    # Create the temporary text object for width measurement
    temp_text = ax.text(0, 0, text,fontsize=8)
    
    # Get the width of the text
    text_width = temp_text.get_window_extent(renderer=renderer).width
    
    # Remove the temporary text object
    temp_text.remove()
    
    # Check if there is enough space to display the text
    if size > text_width and legend_mode == 'Mix':
        ax.text(start + length / 2, 0.5, text, ha='center', va='bottom', fontsize=8, color='black')
    else:
        if text not in legend_dict:
            legend_dict[text] = {'color': color}


# Function to adjust tick positions and labels
def adjusted_ticks_and_labels(minimum, maximum, bed, b,num_ticks):
    """
    Adjust ticks and labels for the X-axis based on BED file regions.
    
    Parameters:
    bed (pd.DataFrame): BED DataFrame.
    
    Returns:
    list, list: Adjusted ticks and labels.
    """
    tick_step = (maximum - minimum) / num_ticks
    equidistant_ticks = np.linspace(minimum, maximum, num_ticks + 1)
    adj_labels = []
    for tick in equidistant_ticks:
        adj_labels.append(round(reverse_adjust_position(tick, bed, b)))
    return equidistant_ticks, adj_labels


def plot_genomic_regions(df, ax, legend_loc='upper left', xlabel='', title='', legend_title='Legend',legend_mode='Mix', color_dict=None, Maximum=None, minumum=None):
    """
    Plot genomic regions on the given axis.
    
    Parameters:
    bed (pd.DataFrame): BED DataFrame.
    ax (matplotlib.axes._subplots.AxesSubplot): Axis to plot on.
    b (int): Baseline adjustment.
    legend_loc (str): Legend location.
    num_ticks (int): Number of ticks on the X-axis.
    xlabel (str): Label for the X-axis.
    title (str): Title of the plot.
    
    Returns:
    None
    """    
    if color_dict ==None :
        color_dict=create_color_dict(bed)
    # Dictionary to store legend entries
    legend_dict = {}

    if legend_mode not in ['Mix','Full']:
        raise ValueError(f"option '{legend_mode}' is not valid (Valid values : ['Mix','Full'])")

    # Plot each region
    previous_end = None
    #ax.plot([0, 0], [0.75, 0.25], color='grey', linestyle='-')
    #ax.plot([Maximum, Maximum], [0.75, 0.25], color='grey', linestyle='-')
    for index, row in df.iterrows():
        start, end = row['start'], row['end']
        length = end - start
        if previous_end is None and start > 0:
            ax.plot([0, start], [0.5, 0.5], color='grey', linestyle='-')
        # Plot regions not described in the BED file as straight lines
        if previous_end is not None and start > previous_end:
            ax.plot([previous_end, start], [0.5, 0.5], color='grey', linestyle='-')

        # Plot described regions
        if 'intron' in row['name']:
            # Plot as a broken line joined upward
            mid_point = start + (end - start) / 2
            ax.plot([start, mid_point], [0.5, 0.5 + 0.5], color='black', linestyle='--')
            ax.plot([mid_point, end], [0.5 + 0.5, 0.5], color='black', linestyle='--')
        elif row['name'] != None:
            color = color_dict[row['name']]
            renderer = ax.figure.canvas.get_renderer()
            #rect = plt.Rectangle((start, 0), end - start, 1, color=color, alpha=0.5)
            rect = patches.FancyBboxPatch((start, 0), end - start, 1, color=color, alpha=0.5,boxstyle="round4")
            box=ax.add_patch(rect)
            rect_extent = rect.get_window_extent(renderer=renderer)
            rect_width = rect_extent.width
            # Add text inside the box or to the legend
            add_text_or_legend(ax, start, end - start,rect_width , row['name'], color, legend_dict, legend_mode)
        else:
            ax.plot([start, end], [0.5, 0.5], color=color, linewidth=2)

        previous_end = end
    if previous_end < Maximum :
        ax.plot([previous_end, Maximum], [0.5, 0.5], color='grey', linestyle='-')

    # Customize plot
    ax.set_ylim([0, 1])
    ax.axis('off')
    ax.set_yticks([])
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.xaxis.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_frame_on(False)
    if legend_dict:
        legend_handles = [Patch(color=info['color'], label=f"{name}") for name, info in legend_dict.items()]
        return legend_handles


def create_lollipop_plot(ax, x, y, color='b', marker='o', line_style='-', line_width=2, alpha=1.0, size=6, Custom_Xaxis=False,stemline_remove=False):
    """
    Create a lollipop plot on the specified axes.
    
    Parameters:
    ax (matplotlib.axes.Axes): The axes to plot on.
    x (list or array): The x values of the points.
    y (list or array): The y values of the points.
    color (str or list): Color of the markers. Can be a single color or a list of colors. Default is 'blue'.
    marker (str or list): Marker style. Can be a single marker style or a list of marker styles. Default is 'o'.
    line_style (str): Line style. Default is '-' (solid line).
    line_width (float): Line width. Default is 2.
    alpha (float or list): Transparency of the markers. Can be a single value or a list of values. Default is 1.0.
    size (int or list): Size of the markers. Can be a single value or a list of values. Default is 6.
    Custom_Xaxis ([[int],[str]]) : list of xticks and x labels 
    
    Returns:
    None
    """
    if isinstance(color, str):
        color = [color] * len(x)
    
    if isinstance(marker, str):
        marker = [marker] * len(x)
    
    if isinstance(alpha, (int, float)):
        alpha = [alpha] * len(x)
    
    if isinstance(size, (int, float)):
        size = [size] * len(x)
    
    (markers, stemlines, baseline) = ax.stem(x, y, linefmt='gray', markerfmt=" ", basefmt=" ")
    if stemline_remove :
        stemlines.remove()
        baseline.remove()
    else :
        plt.setp(stemlines, alpha=0.15)
        plt.setp(baseline, alpha=0.15)
    for xi, yi, ci, mi, ai, si in zip(x, y, color, marker, alpha, size):
        ax.plot([xi], [yi], marker=mi, color=ci, linestyle='None', alpha=ai, markersize=si)
    
    # Customize the plot appearance
    ax.set_ylabel('Log Fold Change')
    
    # Remove grid
    ax.grid(False)
    if Custom_Xaxis :
        ax.set_xticks(Custom_Xaxis[0])
        ax.set_xticklabels(Custom_Xaxis[1])
        ax.set_xlabel('Amino Acid Position')
    else :
        ax.xaxis.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_frame_on(False)

def add_legend(ax, consequence_mapping, pvalue_mapping,transparency=False, add_legend=False ,loc='upper left'):
    """
    Add a legend to the plot for the given consequence mapping and P-value sizes.
    
    Parameters:
    ax (matplotlib.axes.Axes): The axes to add the legend to.
    consequence_mapping (dict): A dictionary mapping consequence types to their colors and markers.
    pvalue_mapping (dict): A dictionary mapping P-value levels to their sizes.
    
    Returns:
    None
    """
    # Create custom legend handles
    consequence_legend_elements = [
        mlines.Line2D([], [], color=value[0], marker=value[1], linestyle='None', markersize=10, label=key)
        for key, value in consequence_mapping.items()
    ]
    
    size_legend_elements = [
        mlines.Line2D([], [], color='gray', marker='o', linestyle='None', markersize=size, label=f'{pvalue}')
        for pvalue, size in pvalue_mapping.items()
    ]
    
    transparency_legend_elements = [
        mlines.Line2D([], [], color='purple', marker='o', linestyle='None', markersize=10, alpha=alpha, label=label)
        for alpha, label in zip([0.2, 1.0], ['p-value<0.05', 'p-value>0.05'])
    ]

    # Add subtitles
    subtitle_fontsize = 'medium'
    leg_a=ax.legend(handles=consequence_legend_elements,loc='upper left', title='Consequences', handlelength=1,bbox_to_anchor=(1, 1), fontsize=subtitle_fontsize, labelspacing=1.25,frameon=False, alignment='left') #,bbox_transform=fig.transFigure
    if transparency :
        leg_c=ax.legend(handles=transparency_legend_elements,loc='center left', title='Biologically Significance', handlelength=1,bbox_to_anchor=(1, 0.80), fontsize=subtitle_fontsize, labelspacing=1.25, frameon=False, alignment='left') #,bbox_transform=fig.transFigure
    leg_b=ax.legend(handles=size_legend_elements, loc='lower left',title='P-values', handlelength=1,bbox_to_anchor=(1, 0.57), fontsize=subtitle_fontsize, labelspacing=1.25,frameon=False, alignment='left') #,bbox_transform=fig.transFigure
    ax.add_artist(leg_a)
    ax.add_artist(leg_b)
    if transparency :
        ax.add_artist(leg_c)
    if add_legend :
        leg_d=ax.legend(handles=add_legend, handlelength=1, loc='upper left',bbox_to_anchor=(1, 0.5), labelspacing=1.25,frameon=False, alignment='left',title='Domains')
        ax.add_artist(leg_b)


def highlight_region(df, ax, name,color_dict=None, negative=False,):
    """
    Add a rectangle to highlight a region in the BED DataFrame by name.
    
    Parameters:
    df (pd.DataFrame): BED DataFrame with at least 'name', 'start', and 'end' columns.
    ax (matplotlib.axes.Axes): Axes object to draw the rectangle on.
    name (str): Name of the region to highlight.
    """
    # Check if the name exists in the DataFrame
    if name not in df['name'].values:
        raise ValueError(f"Name '{name}' not found in the DataFrame")
    if color_dict ==None :
        color_dict=create_color_dict(bed)
    # Get the start and end positions
    region = df[df['name'] == name].iloc[0]
    start = region['start']
    end = region['end']
    bottom, top = ax.get_ylim()
    origin= min(bottom, top) if negative else 0
    length= abs(origin) if negative else top
    # Add a rectangle covering the whole y-axis
    rect = plt.Rectangle((start, origin),
                             end - start, 
                             length,
                         color=color_dict[name], alpha=0.3)
    ax.add_patch(rect)

############################################################# Read BED file



####################################
# Plot each region
if __name__ == '__main__':
    args = parser.parse_args()
    bed = pd.read_csv(args.bed_file, sep='\t', header=0)
    data = Transform_MaGeCK(args.result_file)
    bed_ajusted=bed.copy()
    plt.rcParams.update({'font.size': 18})
    if args.negative_control_file :
        negative_controls = Transform_MaGeCK(args.negative_control_file)
        distribution = negative_controls['LFC']
        data['Bio_Sig']=[calc_biological_significance(i, distribution,args.stat_method,args.Biological_threshold) for i in data['LFC']]
        figB, axB= plt.subplots(1, 1,figsize=(15,15))
        axB.hist(distribution, color=plt.cm.Paired(0))
        figB.savefig(args.output + 'negative_controls_distribution.pdf',format="pdf",bbox_inches="tight")
        MAX = distribution.quantile(q=0.95)
        MIN = distribution.quantile(q=0.05)
        print((MAX,MIN))
    b = args.collapse_factor  # Example value for b
    if args.collapse_factor !=1 and args.histogram :
        raise ValueError('Collapse_factor cannot be used in with histogram flag')
    VEP= pd.DataFrame(parse_VEP(args.vep_file, variant_consequences_mapping))
    Variant_effect=dict(zip(VEP[0],VEP[2]))
    Position_dic=dict(zip(VEP[0],VEP[1]))
    data['Consequence'] = [Variant_effect[i] for i in data['id']]
    data['Position']=[Position_dic[i] for i in data['id']]
    data=data.drop(data[data['Position'].isna()].index)
    data=data.loc[[i!=None for i in data['Consequence']],:]
    data['Position']=[adjust_position(i,bed,b) for i in data['Position']]
    bed_ajusted['start']=[adjust_position(i,bed,b) for i in bed['start']]
    bed_ajusted['end']=[adjust_position(i,bed,b) for i in bed['end']]
    maximum=max(pd.concat([bed_ajusted['end'],data['Position']]))
    minimum=min(pd.concat([bed_ajusted['start'],data['Position']]))

    # Constants

    if args.scheme_loc == 'top' :
        ratios= [1, 20]
        nfigure = 2
    elif args.scheme_loc == 'middle' :
        ratios= [10, 1, 10]
        nfigure = 3
    else :
        ratios= [20, 1]
        nfigure = 2       

    fig = plt.figure(figsize=(25, 25 if args.histogram else 24,))
    if args.histogram :
        gs = fig.add_gridspec(nfigure + args.histogram, hspace=0.05, wspace=0, height_ratios=[2] + ratios, left=0.05, right=1, top=0.95, bottom=0.05)
        ax = gs.subplots(sharex=True)
        ax[0].hist(data['Position'], bins=300, color=plt.cm.Paired(0), edgecolor='none')
        ax[0].xaxis.set_visible(False)
        for spine in ax[0].spines.values():
            spine.set_visible(False)
    else :  
        gs = fig.add_gridspec(nfigure + args.histogram, hspace=0,wspace=0, height_ratios=ratios)
        ax = gs.subplots(sharex=True)

    if args.Biological_REP:
        data=data.drop(data[~data['Bio_Sig']].index)
    if args.p_thresh:
        data=data.drop(data[data['p_value']>args.p_thresh].index)
    # Apply mappings and conditions
    colors = data['Consequence'].map(lambda x: consequence_mapping[x][0])
    markers = data['Consequence'].map(lambda x: consequence_mapping[x][1])
    if args.negative_control_file :
        alphas = data['Bio_Sig'].map(lambda x: 0.2 if not x else 1.0)
    data['logpvalue']=data['p_value'].apply(lambda x:-np.log10(x))
    sizes = data['logpvalue'].apply(lambda x: 2 + 10 * (x / max(data['logpvalue'])))

    pvalue_levels = [0.05, 0.01, 0.001]
    max_pvalue = np.max(data['logpvalue'])
    pvalue_mapping = {p: 2 + 10 * (-np.log10(p) / max_pvalue) for p in pvalue_levels}

    ticks, labels = adjusted_ticks_and_labels(minimum, maximum, bed, b, 10)

    # Plot
    if args.scheme_loc == 'middle' :
        create_lollipop_plot(ax[0 + args.histogram], data.loc[data['LFC']>0,'Position'], data.loc[data['LFC']>0,'LFC'], color=colors[data['LFC']>0], marker=markers[data['LFC']>0], alpha=alphas[data['LFC']>0], size=sizes[data['LFC']>0],stemline_remove=args.no_stem)
        if args.highlight:
            highlight_region(bed_ajusted,ax[0 + args.histogram],args.highlight)
        create_lollipop_plot(ax[2 + args.histogram], data.loc[data['LFC']<0,'Position'], data.loc[data['LFC']<0,'LFC'], color=colors[data['LFC']<0], marker=markers[data['LFC']<0], alpha=alphas[data['LFC']<0], size=sizes[data['LFC']<0],Custom_Xaxis=[ticks, labels],stemline_remove=args.no_stem)
        if args.highlight:
            highlight_region(bed_ajusted,ax[2 + args.histogram],args.highlight ,None,True)
        leg=plot_genomic_regions(bed_ajusted,ax[1 + args.histogram], legend_loc='upper left', xlabel='', title='',legend_title='Features',legend_mode=args.legend_mode,Maximum=maximum)
        add_legend(fig, consequence_mapping, pvalue_mapping,transparency=(not args.Biological_REP),add_legend=leg)
    else :
        create_lollipop_plot(ax[(args.scheme_loc == 'bottom') + 1 + args.histogram], data['Position'], data['LFC'], color=colors, marker=markers, alpha=alphas, size=sizes,stemline_remove=args.no_stem)
        leg=plot_genomic_regions(bed_ajusted,ax[(args.scheme_loc == 'bottom') + 0 + args.histogram], legend_loc='upper left', xlabel='', title='',legend_title='Features',legend_mode=args.legend_mode,Maximum=maximum)
        add_legend(fig, consequence_mapping, pvalue_mapping,transparency=(not args.Biological_REP),add_legend=leg)
        if args.highlight:
            highlight_region(bed_ajusted,ax[(args.scheme_loc == 'bottom') + 1 + args.histogram],args.highlight)
            highlight_region(bed_ajusted,ax[(args.scheme_loc == 'bottom') + 1 + args.histogram],args.highlight,None,True)
    plt.tight_layout()
    fig.savefig(args.output + '.pdf',format="pdf",bbox_inches="tight")
