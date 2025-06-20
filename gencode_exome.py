import argparse
import gzip
import pysam
import sys
from collections import Counter
from intervaltree import IntervalTree

argparser = argparse.ArgumentParser(description = 'Extracts exome from GENCODE GTF.')
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
            if feature_type != 'exon':
                continue
            chrom = fields[0]
            start_bp, stop_bp = map(int, fields[3:5])
            strand = fields[6]
            attributes = dict()
            for x in fields[8].split(';'):
                if x == '':
                    continue
                key, value = map(lambda x: x.strip('"'), x.strip().split())
                if key in {'tag', 'ont'}:
                    attributes.setdefault(key, set()).add(value)
                else:
                    assert key not in attributes, key
                    attributes[key] = value
            if attributes['gene_type'] != 'protein_coding' or attributes['transcript_type'] != 'protein_coding':
                continue
            # skip if automatically annotated locus
            if attributes['level'] not in ['1', '2']:
                continue
            tags = attributes.get('tag', {})
            if 'CCDS' not in tags:
                continue
            if 'Ensembl_canonical' not in tags:
                continue
            if 'basic' not in tags:
                continue
            gene_id = attributes['gene_id']
            gene_name = attributes['gene_name']
            gene_data = genes.setdefault(gene_id, { 'name': gene_name, 'chrom': chrom, 'exons': IntervalTree() })
            gene_data['exons'].addi(start_bp, stop_bp + 1)
    return genes


if __name__ == '__main__':
    args = argparser.parse_args()
    genes = get_genes_cds(args.in_gencode_file)
    print(f'Selected {len(genes)} gene(s)')
    chromosomes = dict()
    for gene_id, gene_data in genes.items():
        gene_data['exons'].merge_overlaps() # merge overlapping exons from multiple transcripts (this is needed if also non-canonical trancstipts were used)
        chromosome_exons = chromosomes.setdefault(gene_data['chrom'], IntervalTree())
        for exon in gene_data['exons']:
            chromosome_exons.addi(exon.begin, exon.end)
    with open(args.out_BED, 'w') as ofile:
        for chromosome_name, chromosome_exons in chromosomes.items():
            chromosome_exons.merge_overlaps() # merge overlapping exons from multiple genes
            for exon in sorted(chromosome_exons):
                ofile.write('{}\t{}\t{}\n'.format(chromosome_name, exon.begin, exon.end - 1))

