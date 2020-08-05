import argparse
import gzip
import pysam
import sys
from collections import Counter
from intervaltree import IntervalTree

argparser = argparse.ArgumentParser(description = 'Extracts CDS from GENCODE GTF.')
argparser.add_argument('-g', '--gencode', metavar = 'file', dest = 'in_gencode_file', required = True, help = 'Input GENCODE file in GTF format.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_BED', required = True, help = 'Output file in BED format.')


def get_genes_cds(in_gencode_file):
    genes = dict()
    with gzip.open(in_gencode_file, 'rt') as ifile:
        for line in ifile:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            feature_type = fields[2]
            # start codon is already included to CDS
            # stop codon is not included to CDS by GENCODE, but lets not consider it at this point
            if feature_type != 'CDS':
                continue
            chrom = fields[0]
            start_bp, stop_bp = map(int, fields[3:5])
            strand = fields[6]
            attributes = dict(map(lambda x: x.strip('"'), x.strip().split()) for x in fields[8].split(';') if x != '')
            if attributes['gene_type'] != 'protein_coding' or attributes['transcript_type'] != 'protein_coding':
                continue
            # skip if the coding region end or start could not be confirmed
            if any(x in ['cds_end_NF', 'cds_start_NF', 'mRNA_end_NF', 'mRNA_start_NF'] for x in attributes.get('tag', [])):
                continue
            gene_id = attributes['gene_id']
            gene_name = attributes['gene_name']
            gene_data = genes.setdefault(gene_id, { 'name': gene_name, 'chrom': chrom, 'cds': IntervalTree() })
            gene_data['cds'].addi(start_bp, stop_bp + 1)
    return genes


if __name__ == '__main__':
    args = argparser.parse_args()
    genes = get_genes_cds(args.in_gencode_file)
    with open(args.out_BED, 'w') as ofile:
        for gene_id, gene_data in genes.items():
            gene_data['cds'].merge_overlaps() # merge overlapping CDS from multiple transcripts
            for cds in sorted(gene_data['cds']):
                ofile.write('{}\t{}\t{}\t{}\t{}\n'.format(gene_data['chrom'], cds.begin, cds.end - 1, gene_id, gene_data['name']))

