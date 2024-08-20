#!/usr/bin/env python3


'''
AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>
VERSION: 1.1
YEAR: 2024

Goal :

To produce lollipop plot with schema showing features at location. Lollipop graphic parameters represented p-value/fdr (size), lfc (y axis), position (x axis), Biological significance (optional)

table of Content : 

    Parsers:
      parse_VEP(filename, variant_consequences_mapping, position): Parses a VEP file and yields specific variant information.
      Transform_MaGeCK(filename): Transforms a MaGeCK output file into a specific DataFrame format.
    Color Pack (High reusability):
      get_color_palette(num_entries): Returns a color palette suitable for a given number of entries.
      get_color(index, palette, existing_colors): Generates a unique color from a palette, avoiding duplicates.
      create_color_dict(df): Creates a dictionary mapping names to colors, excluding rows with "intron" in their name.
      is_valid_color(color): Checks if a color string is valid in matplotlib.
    Collapse intron :
      adjust_position(pos, df, b): Adjusts a position based on BED file regions and a condensation factor.
      reverse_adjust_position(pos, df, b): Reverses the adjustment of a position based on BED file regions and a condensation factor.
      adjusted_ticks_and_labels(minimum, maximum, bed, b, num_ticks): Adjusts tick positions and labels on a plot based on BED file regions.
    Graph Pack :
      plot_genomic_regions(df, ax, legend_loc, title, legend_title, legend_mode, color_dict, Maximum, minumum, Custom_Xaxis, xlabel): Plots genomic regions on a given axis, with custom legend and axis options.
      add_text_or_legend(ax, start, length, size, text, color, legend_dict, legend_mode): Adds text or legend entry to a plot for a genomic feature.
      create_lollipop_plot(ax, x, y, color, marker, line_style, line_width, alpha, size, Custom_Xaxis, stemline_remove, fdr, xlabel): Creates a lollipop plot on the specified axes with various customization options.
      add_legend(ax, consequence_mapping, pvalue_mapping, transparency, add_legend, loc, fdr): Adds a legend to a plot for consequence types and p-value levels.
      highlight_region(df, ax, name, color_dict, negative): Highlights a region in a plot based on a specified name from a BED DataFrame.
    Auxiliary: 
      convert_to_int(s): Converts a string to an integer, returns None if conversion fails.
      calc_biological_significance(value, distribution, method, thresh): Determines the biological significance of a value against a distribution using a specified method.

 Example:
   python3 BEscreen_lollipop_plot.py -b path/to/features.bed -i gene_summary.txt -v path/to/VEP.tsv --histogram -m Full --no_stem --scheme_location top -n 'custom:chr19:55113873:55117983' --stat_method sign_test


'''


import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.lines as mlines
from matplotlib.patches import Patch
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import argparse
from scipy.stats import rankdata
from scipy.stats import binomtest
import sys
import re
import warnings




parser = argparse.ArgumentParser(description='This script produce a lolipop plot (with domain/feature annotations) base on a MaGeCK, a Variant Effect Prediction file (VEP), and bed style file')
parser.add_argument('-b', '--bed', metavar='FILE', dest='bed_file', required=True, type=str, help='space/tab-delimited bed file of features, must contains columns : protein, start, stop and name (optional color)')
parser.add_argument('-i', '--input', metavar='FILE', dest='result_file', required=True, type=str, help='MaGeCK output or file with tab/space delimited data containing columns : (id,lfc, p-value,fdr) id requiere format Protein_number')
parser.add_argument('-v', '--vep', metavar='FILE', dest='vep_file', required=True, type=str, help='tab delimited variant effect prediction file (VEP) with id corresponding to input')
parser.add_argument('-c', '--intron_collapse', metavar='number', dest='collapse_factor', required=False,default=1 , type=int, help='factor by which to collapse introns (default uncollapsed), intron are defined by name containing "intron"')
parser.add_argument('-n', '--negative_Control', metavar='str', dest='Control', required=False, nargs='+', type=str, help='list of negative control protein/region (if no negative controls (empty widows/regions) Biological significance ignored)')
parser.add_argument('--empty_controls',  metavar='str', dest='Control_empty', required=False, help='Use the list of empty windows sgrna (default empty widows are ignored)')
parser.add_argument('--stat_method',metavar='str', dest='stat_method', required=False,default='quantile' , choices={"binom_sign", "quantile",'sign_test'}, type=str, help='Statistical method to calculate p-value of Biological significance, quantile uses p-value as q-quantile limit')
parser.add_argument('--Bio_threshold', dest='Biological_threshold', metavar='float',type=float, required=False, default='0.05', help='Biological threshold (two-sided)')
parser.add_argument('--scheme_location', metavar='str', dest='scheme_loc', required=False,default='middle' , type=str, choices={"top",'bottom', "middle"}, help='Location of the protein/gene scheme ("top"/"middle"/"bottom") (default : middle)')
parser.add_argument('--histogram', dest='histogram', required=False,action='store_true', help='Flag for wheather a histogram should represent the coverage')
parser.add_argument('-p', '--Prob_Threshold', metavar='float',type=float, dest='p_thresh', required=False, default='1', help='P-Value Threshold for loolipop representation')
parser.add_argument('-B', '--Biological_Sig', dest='Biological_REQ', action='store_true', help='Flag for Biological significance requirerement in lolipop plot (not present not represented)')
parser.add_argument('--no_stem', dest='no_stem', action='store_true', help='Flag to remove the stemlines')
parser.add_argument('-F', dest='fdr', action='store_true', help='Flag to use FDR (default p-value)')
parser.add_argument('--X', metavar='str', dest='X_axis', required=False,default='AA' , type=str, choices={"AA",'Nuc'}, help='X axis label')
parser.add_argument('-m', '--legend_mode', metavar='str', dest='legend_mode', required=False,default='Mix' , choices={"Full", "Mix"}, type=str, help='legend mode (possible values Mix (default) or Full), full all features are in legend, mix of on-graph and on-legend')
parser.add_argument('--highlight_region', metavar='str', dest='highlight', required=False,default=None ,nargs='+', type=str, help='Feature to be highlighted expect format protein-feature (can be multiple)')
parser.add_argument('--out', dest='output', default='output', type=str, help='Prefix of output image file name (output.png).')

plt.rcParams.update({'font.size': 18})


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
'NMD_transcript_variant': 'non-sense',
'upstream_gene_variant': 'none',
'3_prime_UTR_variant': 'none',
'synonymous_variant': 'synonymous',
'non_coding_transcript_exon_variant': 'none',
'splice_region_variant': 'splice',
'splice_polypyrimidine_tract_variant': 'splice',
'stop_gained': 'non-sense',
'coding_sequence_variant': 'none',
'5_prime_UTR_variant': 'none',
'regulatory_region_variant': 'none',
'splice_donor_variant': 'splice',
'splice_acceptor_variant': 'splice',
'non_coding_transcript_variant': 'none',
'splice_donor_region_variant': 'splice',
'splice_donor_5th_base_variant': 'splice',
'TF_binding_site_variant': 'none',
'start_lost': 'non-sense',
'incomplete_terminal_codon_variant': 'non-sense'
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
    'lfc': [None] * len(Mage.id),
    'p_value': [None] * len(Mage.id),
    'fdr': [None] * len(Mage.id)
    })
    returned['lfc'] = Mage['pos|lfc']
    returned['p_value'] = [row['neg|p-value'] if row['pos|lfc']<0 else row['pos|p-value'] for _, row in Mage.iterrows()]
    returned['fdr'] = [row['neg|fdr'] if row['pos|lfc']<0 else row['pos|fdr'] for _, row in Mage.iterrows()]
    return returned

def parse_VEP(filename, variant_consequences_mapping, position):
    # From Oligomer.py https://github.com/CERC-Genomic-Medicine/CRISPR_Library_prep.git
    with open(filename) as f:
        for line in f:
            if line[0:2]=='##':
                continue
            elif line[0:1]=='#':
                header = line.strip('\n').split("\t")
                colConsequence = header.index('Consequence')
                colAA = header.index('Protein_position')
                colIMPACT= header.index('IMPACT')
                colAmino = header.index('Amino_acids')
                colLoc=header.index('Location')
            else :
                if 'header' not in locals():
                    raise Exception("VEP format not as expected") 
                else :
                    record = line.strip('\n').split("\t")
                    ID = record[0]
                    if position == 'AA':
                        POSi = record[colAA].split('-')[0] if record[colAA].split('-')[0] != '?' else record[colAA].split('-')[1]
                        POS = convert_to_int(POSi)
                    elif position == 'Nuc':
                        Pos= re.search(r'chr\d+:(\d+)-\d+', record[colLoc]).group(1)
                    if 'splice' in [variant_consequences_mapping[i] for i in record[colConsequence].split(',')]:
                        consequence = 'splice'
                    elif 'non-sense' in [variant_consequences_mapping[i] for i in record[colConsequence].split(',')]:
                        consequence = 'non-sense'
                    elif 'missense' in [variant_consequences_mapping[i] for i in record[colConsequence].split(',')]:
                        consequence = 'missense'
                    elif 'synonymous' in [variant_consequences_mapping[i] for i in record[colConsequence].split(',')]:
                        consequence = 'synonymous'
                    else :
                        consequence = None
                    yield ID, POS, consequence, record[colConsequence], record[colIMPACT],record[colAmino]

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
    if method == 'quantile':
        MAX = distribution.quantile(q=1-thresh)
        MIN = distribution.quantile(q=thresh)
        sign=value<=MIN or value>=MAX
    if method == 'binom_test':
        # Count how many values are greater than or equal to the test value
        greater_equal_count = sum(v >= value for v in distribution)
        n = len(distribution)
        
        # Perform the binom test
        p_value = binomtest(greater_equal_count, n, 0.5, alternative='two-sided').pvalue
        sign= thresh>=p_value
    
    elif method == 'sign_test':
        # Count positive and negative differences
        positive_differences = sum((v - value) > 0 for v in distribution)
        negative_differences = sum((v - value) < 0 for v in distribution)
        
        # Perform a binomial test for the sign test
        n = positive_differences + negative_differences
        p_value = binomtest(positive_differences, n, 0.5, alternative='two-sided').pvalue
        sign= thresh>=p_value
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
    pastel1 = plt.get_cmap('Pastel1').colors
    tab20 = plt.get_cmap('tab20').colors
    if num_entries <= len(pastel1):
        return pastel1
    elif num_entries <= len(tab20):
        return tab20
    else:
        warnings.warn("Warning: Number of entries exceeds the available colors in tab20 palette.")
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
    df=df[~df['name'].isna()].copy()
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
        if not row['name'] and 'intron' in row['name']:
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
        if not row['name'] and 'intron' in row['name']:
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


def plot_genomic_regions(df, ax, legend_loc='upper left', title='', legend_title='Legend',legend_mode='Mix', color_dict=None, Maximum=None, minumum=None, Custom_Xaxis=False, xlabel='Amino acid position'):
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
        if not row['name'] and 'intron' in row['name']:
            # Plot as a broken line joined upward
            mid_point = start + (end - start) / 2
            ax.plot([start, mid_point], [0.5, 0.5 + 0.5], color='black', linestyle='--')
            ax.plot([mid_point, end], [0.5 + 0.5, 0.5], color='black', linestyle='--')
        elif not pd.isna(row['name']):
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
            ax.plot([start, end], [0.5, 0.5], color='grey', linestyle='-')

        previous_end = end
    if previous_end < Maximum :
        ax.plot([previous_end, Maximum], [0.5, 0.5], color='grey', linestyle='-')

    # Customize plot
    ax.set_ylim([0, 1])
    ax.grid(False)
    ax.yaxis.set_visible(False)
    if Custom_Xaxis :
        ax.set_xticks(Custom_Xaxis[0])
        ax.set_xticklabels(Custom_Xaxis[1])
        ax.set_xlabel(xlabel)
    else :
        ax.xaxis.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_frame_on(False)
    if legend_dict:
        legend_handles = [Patch(color=info['color'], label=f"{name}") for name, info in legend_dict.items()]
        return legend_handles


def create_lollipop_plot(ax, x, y, color='b', marker='o', line_style='-', line_width=2, alpha=1.0, size=6, Custom_Xaxis=False,stemline_remove=False, fdr=False, xlabel='Amino acid position'):
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

    (markers, stemlines, baseline) = ax.stem(x.to_numpy(), y.to_numpy(), linefmt='gray', markerfmt=" ", basefmt=" ")
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
        ax.set_xlabel(xlabel)
    else :
        ax.xaxis.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_frame_on(False)

def add_legend(ax, consequence_mapping, pvalue_mapping,transparency=None, add_legend=False ,loc='upper left', fdr=False):
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
    if transparency:
        transparency_legend_elements = [
            mlines.Line2D([], [], color='purple', marker='o', linestyle='None', markersize=10, alpha=alpha, label=label)
            for alpha, label in zip([0.4, 1.0], transparency)
        ]

    # Add subtitles
    subtitle_fontsize = 'medium'
    leg_a=ax.legend(handles=consequence_legend_elements,loc='upper left', title='Consequences', handlelength=1,bbox_to_anchor=(1, 1), fontsize=subtitle_fontsize, labelspacing=1.25,frameon=False) #,bbox_transform=fig.transFigure
    if transparency :
        leg_c=ax.legend(handles=transparency_legend_elements,loc='center left', title='Biologically Significance', handlelength=1,bbox_to_anchor=(1, 0.80), fontsize=subtitle_fontsize, labelspacing=1.25, frameon=False) #,bbox_transform=fig.transFigure
    leg_b=ax.legend(handles=size_legend_elements, loc='lower left',title='FDR' if fdr else 'P-value', handlelength=1,bbox_to_anchor=(1, 0.57), fontsize=subtitle_fontsize, labelspacing=1.25,frameon=False) #,bbox_transform=fig.transFigure
    ax.add_artist(leg_a)
    ax.add_artist(leg_b)
    if transparency :
        ax.add_artist(leg_c)
    if add_legend :
        leg_d=ax.legend(handles=add_legend, handlelength=1, loc='upper left',bbox_to_anchor=(1, 0.5), labelspacing=1.25,frameon=False,title='Domains')
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
    #### Input Parsing
    bed_full = pd.read_csv(args.bed_file, sep='\t', header=0)
    if not set(['start','end','name','proteins']).issubset(bed_full.columns):
        raise ValueError('Bed-like file does not contain the relevant columns')
    try :
        data_full = pd.read_csv(args.result_file, sep='\s+',header=0, usecols=['id','lfc', 'p-value','fdr'])
    except:
        try :
            data_full = Transform_MaGeCK(args.result_file)
        except: 
            raise ValueError('Main input data (-i) does not have a supported format')
    data_full['proteins']=[i.split('_')[0] for i in data_full['id']]
    b = args.collapse_factor  # Example value for b
    Biological_REP = args.Control or args.Control_empty # Whether or not Biological significance is to be analysed

    if args.collapse_factor !=1 and args.histogram :
        raise ValueError('Collapse_factor cannot be used in with histogram flag')

    if args.Control :
        ### Treating the controls
        protein=set(data_full['proteins'].values)
        if not set(args.Control).issubset(protein) :
            raise ValueError(f'specified negative control(s) not found in main input data {set(args.Control)-protein}')
        neg_controls_true = [i in args.Control for i in data_full['proteins']]
        negative_controls = data_full.loc[neg_controls_true,:]
        data_full=data_full.loc[~np.array(neg_controls_true),:] 

    ### Analysing the Vep
    VEP= pd.DataFrame(parse_VEP(args.vep_file, variant_consequences_mapping,position=args.X_axis))
    Variant_effect=dict(zip(VEP[0],VEP[2]))
    Position_dic=dict(zip(VEP[0],VEP[1]))
    Variant_effect_full=dict(zip(VEP[0],VEP[3]))
    Variant_impact=dict(zip(VEP[0],VEP[4]))
    Variant_Amino=dict(zip(VEP[0],VEP[5]))

    ### treating Empty controls
    if args.Control_empty :
        empties = pd.read_csv(args.Control_empty, sep='\t', header=None)[0]
        empties_negative = data_full[data_full['id'].isin(empties)]
        negative_controls=pd.concat([negative_controls,empties_negative])
        data_full = data_full[~data_full['id'].isin(empties)]
    
    ### Analysing Biological Significance
    if Biological_REP:
        distribution = negative_controls['lfc']
        data_full['Bio_Sig']=[calc_biological_significance(i, distribution,args.stat_method,args.Biological_threshold) for i in data_full['lfc']]
        figB, axB= plt.subplots(1, 1,figsize=(15,15))
        axB.hist(distribution, color=plt.cm.Paired(0))
        figB.savefig(args.output + '_negative_controls_distribution.pdf',format="pdf",bbox_inches="tight")
    ### If we don't analyse Biological significance all true
    else :
        data_full['Bio_Sig']=True

    ## trouble shoot warning if sgRNA absent from VEP
    if (~data_full['id'].isin(Variant_effect.keys())).any():
        warnings.warn(f'warning : {sum(~data_full["id"].isin(Variant_effect.keys()))} were filtered out because they were not in found in VEP file')
    data_full=data_full[data_full['id'].isin(Variant_effect.keys())]

    ### Adding VEP data to dataset
    proteins = set([i.split('_')[0] for i in data_full['id']])
    data_full['Consequence'] = [Variant_effect[i] for i in data_full['id']]
    data_full['Position']=[Position_dic[i] for i in data_full['id']]    

     ## trouble shoot warning if sgRNA has no clear position in VEP
    if data_full['Position'].isna().any() and  args.X_axis == 'AA':
        warnings.warn(f'warning : {sum(data_full["Position"].isna())} were filtered out due unclear position \n this only a warning du to X axis is AA (sgRNA can ve outside coding regions) ')
    else :
        raise ValueError(f'warning : {sum(data_full["Position"].isna())} were filtered out due unclear position \n this is an error because X axis is Nucleotide (and thus sgRNA should have a position) \n list {data_full[data_full["Position"].isna(),"id"]}')
    data_full=data_full.drop(data_full[data_full["Position"].isna()].index)

    ### trouble shoot warning if sgRNA has no clear consequence in VEP
    if sum([i==None for i in data_full['Consequence']]) != 0 :
        warnings.warn(f'warning : {sum([i==None for i in data_full["Consequence"]])} were filtered out due to unclear consequence (i.e. no base category could be found)\n')

    data_full=data_full.loc[[i!=None for i in data_full['Consequence']],:].copy()

    ### trouble shoot user error proteins bed and protein input file
    if not proteins.issubset(set(bed_full.proteins)): 
        raise ValueError(f'Not all proteins/regions in the main input file (-i) is present in the Bed-like file:{set(proteins)-set(bed_full.proteins)} \n')

    ### Test highlight list
    if args.highlight:
        if not set(args.highlight).issubset(set(bedfull['proteins'] + '-' + bedfull['name'])):
            raise ValueError(f'Not all proteins/regions to be highlighted are in the bed-like file under the protein-feature format :{set(args.highlight)-set(bedfull["proteins"] + "-" + bedfull["name"])} \n')
    
    for protein in proteins:
        ### Produce a graph per protein for all protein
        ### Produce 
        data=data_full.loc[data_full['proteins']==protein,:]
        bed_adjusted=bed_full.loc[bed_full['proteins']==protein,:]
        bed=bed_adjusted.copy()
        data.loc[:,'Position']=[adjust_position(i,bed,b) for i in data['Position']].copy() #adjust if intron collapse
        bed_adjusted.loc[:,'start']=[adjust_position(i,bed,b) for i in bed['start']].copy()  #adjust if intron collapse
        bed_adjusted.loc[:,'end']=[adjust_position(i,bed,b) for i in bed['end']].copy()  #adjust if intron collapse
        maximum=max(pd.concat([bed_adjusted['end'],data['Position']])) # max in graph
        minimum=min(pd.concat([bed_adjusted['start'],data['Position']])) # min in graph
        
        # Figure Configuration
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
        ### Histogram plot and figure grid config (according to histogram yes/no)
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
        
        if args.fdr :
            data['p_value']=data['fdr']

        ## Represent only biological significant
        if args.Biological_REQ:
            data=data.drop(data[~data['Bio_Sig']].index)
        ## Represent only p-val significant (default 1 so no filter)
        if args.p_thresh:
            data=data.drop(data[data['p_value']>args.p_thresh].index)

        ### parameter data representation
        colors = data['Consequence'].map(lambda x: consequence_mapping[x][0])
        markers = data['Consequence'].map(lambda x: consequence_mapping[x][1])
        alphas = data['Bio_Sig'].map(lambda x: 0.4 if not x else 1.0)
        if (not args.Biological_REQ and Biological_REP):
            if args.stat_method=='quantile':
                Biosig_labels=[fr"$Q_{{{args.Biological_threshold:.2f}}} \leq x \leq Q_{{{1 - args.Biological_threshold:.2f}}}$", fr"$x < Q_{{{args.Biological_threshold:.2f}}}$ or $x > Q_{{{1 - args.Biological_threshold:.2f}}}$"]
            else :
                Biosig_labels = [fr"$p$-value $> {args.Biological_threshold:.2f}$", fr"$p$-value $\leq {args.Biological_threshold:.2f}$"]
        data['logpvalue']=data['p_value'].apply(lambda x:-np.log10(x))
        sizes = data['logpvalue'].apply(lambda x: 2 + 20 * (x / max(data['logpvalue'])))
        pvalue_levels = [0.05, 0.01, 0.001]
        max_pvalue = np.max(data['logpvalue'])
        pvalue_mapping = {p: 2 + 20 * (-np.log10(p) / max_pvalue) for p in pvalue_levels}
        ticks, labels = adjusted_ticks_and_labels(minimum, maximum, bed, b, 10) ## change only if intron collapse
        if args.scheme_loc == 'middle' :
            create_lollipop_plot(ax[0 + args.histogram], data.loc[data['lfc']>=0,'Position'], data.loc[data['lfc']>=0,'lfc'], color=colors[data['lfc']>=0], marker=markers[data['lfc']>=0], alpha=alphas[data['lfc']>=0], size=sizes[data['lfc']>=0],stemline_remove=args.no_stem)
            if args.highlight:
                for hl in args.highlight : 
                    protein_hl, feature=hl.split('-')
                    if protein_hl == protein :
                        highlight_region(bed_adjusted,ax[0 + args.histogram],feature)
            create_lollipop_plot(ax[2 + args.histogram], data.loc[data['lfc']<0,'Position'], data.loc[data['lfc']<0,'lfc'], color=colors[data['lfc']<0], marker=markers[data['lfc']<0], alpha=alphas[data['lfc']<0], size=sizes[data['lfc']<0],Custom_Xaxis=[ticks, labels],stemline_remove=args.no_stem, xlabel = 'Nucleotide position' if args.X_axis == 'nuc' else 'Amino acid position' if args.X_axis == 'AA' else None)
            if args.highlight:
                for hl in args.highlight : 
                    protein_hl, feature=hl.split('-')
                    if protein_hl == protein :
                        highlight_region(bed_adjusted,ax[2 + args.histogram],feature ,None,True)
            leg=plot_genomic_regions(bed_adjusted,ax[1 + args.histogram], legend_loc='upper left', title='',legend_title='Features',legend_mode=args.legend_mode,Maximum=maximum)
            add_legend(fig, consequence_mapping, pvalue_mapping,transparency=Biosig_labels,add_legend=leg,fdr=args.fdr)
        else :
            create_lollipop_plot(ax[(args.scheme_loc != 'bottom') + args.histogram], data['Position'], data['lfc'], color=colors, marker=markers, alpha=alphas, size=sizes,stemline_remove=args.no_stem,Custom_Xaxis=[ticks, labels] if args.scheme_loc == 'top' else None, xlabel = 'Nucleotide position' if args.X_axis == 'nuc' else 'Amino acid position' if args.X_axis == 'AA' else None)
            leg=plot_genomic_regions(bed_adjusted,ax[(args.scheme_loc == 'bottom') + args.histogram], legend_loc='upper left', title='',legend_title='Features',legend_mode=args.legend_mode,Maximum=maximum,Custom_Xaxis=[ticks, labels] if args.scheme_loc == 'bottom' else None ,xlabel = 'Nucleotide position' if args.X_axis == 'nuc' else 'Amino acid position' if args.X_axis == 'AA' else None)
            add_legend(fig, consequence_mapping, pvalue_mapping,transparency=Biosig_labels,add_legend=leg,fdr=args.fdr)
            if args.highlight:
                for hl in args.highlight : 
                    protein_hl, feature=hl.split('-')
                    if protein_hl == protein :
                        highlight_region(bed_adjusted,ax[(args.scheme_loc != 'bottom') + args.histogram],feature)
        fig.savefig(args.output +'_'+protein + '.pdf',format="pdf",bbox_inches="tight")
