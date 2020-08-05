import argparse
import gzip
import pysam
import sys
from collections import Counter
from intervaltree import IntervalTree

argparser = argparse.ArgumentParser(description = 'Extracts 5\'UTR and 3\'UTR regions from GENCODE GTF.')
argparser.add_argument('-g', '--gencode', metavar = 'file', dest = 'in_gencode_file', required = True, help = 'Input GENCODE file in GTF format.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_BED', required = True, help = 'Output file in BED format.')


def get_transcripts_coords(in_gencode_file):
    transcripts = dict()
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
            attributes = dict()
            for key, value in (map(lambda x: x.strip('"'), x.strip().split()) for x in fields[8].split(';') if x != ''):
                if key == 'tag':
                    attributes.setdefault(key, []).append(value)
                else:
                    attributes[key] = value
            if attributes['gene_type'] != 'protein_coding' or attributes['transcript_type'] != 'protein_coding':
                continue
            transcript = transcripts.setdefault(attributes['transcript_id'], IntervalTree())
            transcript.addi(start_bp, stop_bp + 1)
    return dict((transcript_id, [cds.begin(), cds.end() - 1]) for transcript_id, cds in transcripts.items())


def get_genes_utr(in_gencode_file, transcripts_coords):
    genes = dict()
    with gzip.open(in_gencode_file, 'rt') as ifile:
        for line in ifile:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            feature_type = fields[2]
            # start codon is already included to CDS
            # stop codon is not included to CDS by GENCODE, but lets not consider it at this point
            if feature_type != 'UTR':
                continue
            chrom = fields[0]
            start_bp, stop_bp = map(int, fields[3:5])
            strand = fields[6]
            attributes = dict(map(lambda x: x.strip('"'), x.strip().split()) for x in fields[8].split(';') if x != '')
            if attributes['gene_type'] != 'protein_coding' or attributes['transcript_type'] != 'protein_coding':
                continue
            cds_start_bp, cds_stop_bp = transcripts_coords[attributes['transcript_id']]
            if strand == '+':
                if stop_bp <= cds_start_bp:
                    utr_type = '5\'UTR'
                else:
                    utr_type = '3\'UTR'
            else:
                if start_bp >= cds_stop_bp:
                    utr_type = '5\'UTR'
                else:
                    utr_type = '3\'UTR'
            gene_id = attributes['gene_id']
            gene_data = genes.setdefault(gene_id, { 'strand': strand, 'chrom': chrom, 'name': attributes['gene_name'], '5\'UTR': IntervalTree(), '3\'UTR': IntervalTree() })
            gene_data[utr_type].addi(start_bp, stop_bp + 1)
    return genes


if __name__ == '__main__':
    args = argparser.parse_args()
    genes_utrs = get_genes_utr(args.in_gencode_file, get_transcripts_coords(args.in_gencode_file))
    with open(args.out_BED, 'w') as ofile:
        for gene_id, gene_data in genes_utrs.items():
            for utr_type in ['5\'UTR', '3\'UTR']:
                gene_data[utr_type].merge_overlaps()
                for utr in sorted(gene_data[utr_type]):
                    ofile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_data['chrom'], utr.begin, utr.end - 1, gene_data['strand'], utr_type, gene_id, gene_data['name']))

