import argparse
import pysam
from collections import Counter

argparser = argparse.ArgumentParser(description = '')
argparser.add_argument('-b', '--bed', metavar = 'name', dest = 'in_bed', type = str, required = True, help = 'Input BED file.')
argparser.add_argument('-f', '--fasta', metavar = 'file', dest = 'in_ref_fasta', type = str, required = True, help = 'Input reference FASTA file (indexed).')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_bed', required = True, help = 'Output BED file.')


if __name__ == '__main__':
    args = argparser.parse_args()
    with open(args.in_bed, 'rt') as ibed, pysam.FastaFile(args.in_ref_fasta) as ifasta, open(args.out_bed, 'wt') as obed:
        fasta_chroms = set(list(ifasta.references))
        obed.write('# columns appended: no. of A bases, no. of C bases, no. of G bases, no. of T bases, %GC\n')
        for line in ibed:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split()
            chrom = fields[0]
            if chrom not in fasta_chroms:
                chrom = chrom[3:] if chrom.startswith('chr') else f'chr{chrom}'
                if chrom not in fasta_chroms:
                    print(f'Warning: skipping chromosome {fields[0]} because it is not in FASTA file.')
                    continue
            start_bp = int(fields[1])
            stop_bp = int(fields[2])
            base_counts = Counter()
            for base in ifasta.fetch(chrom, start_bp - 1, stop_bp):
                base_counts[base] += 1
            obed.write('{}'.format('\t'.join(fields)))
            n_all = 0
            for base in ['A', 'C', 'G', 'T']:
                n_all += base_counts[base]
                obed.write(f'\t{base_counts[base]}')
            if n_all > 0:
                gc_content = (base_counts['C'] + base_counts['G']) / float(n_all)
            else:
                gc_content = 'NA'
            obed.write(f'\t{gc_content}\n')

