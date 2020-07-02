import argparse
import pysam

#
# Note: python3 + pysam library are required
#
# Example:
#   python vcf_assign_af_bins.py -i input.vcf.gz -o output.vcf.gz -G -- bins variants by specified AF INFO field into "<=0", "(0, 0.001]", "(0.001, 0.005]", "(0.005, 0.01]", "(0.01, 0.05]", "(0.05, 0.1]", "(0.1, 0.2]", ...,  "(0.99, 0.995]", "(0.995, 0.999]", "(0.999, 1.0]" and drops genotype information
#

argparser = argparse.ArgumentParser(description = 'Assigns MAF bin to each alternate allele in VCF/BCF.')
argparser.add_argument('-i', '--in-vcf', metavar = 'file', dest = 'in_VCF', required = True, help = 'Input VCF/BCF file.')
argparser.add_argument('-f', '--info-field', metavar = 'name', dest = 'af_key', type = str, default = 'AF', help = 'Name of the INFO field with alternate allele frequency. Default is "AF"')
argparser.add_argument('-G', '--drop-genotypes', dest = 'drop_gt', action='store_true', help = 'Don\'t include genotype information in the output VCF/BCF file.')
argparser.add_argument('-b', '--af-bins', metavar = 'number', dest = 'af_bins', type = float, nargs = '+', default = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.995, 0.999, 1.0], help = 'List of MAF bin boundaries e.g. 0, 0.1, 0.5, 1.0 => "<=0", "(0, 0.1]", "(0.1, 0.5]", "(0.5, 1.0]". Default: 0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.995, 0.999, 1.0.')
argparser.add_argument('-o', '--out-vcf', metavar = 'file', dest= 'out_VCF', required = True, help = 'Output VCF/BCF file.')


def assign_bin(bins, value):
    if value is not None:
        for min_value, max_value, bin_label in bins:
            if value > min_value and value <= max_value:
                return bin_label
    return None


if __name__ == '__main__':
    args = argparser.parse_args()
    
    assert len(args.af_bins) > 1

    args.af_bins.sort()
    af_bins = [ (float('-inf'), args.af_bins[0], f'<={args.af_bins[0]}')  ]
    for min_af, max_af in zip(args.af_bins, args.af_bins[1:]):
        af_bins.append((min_af, max_af, f'({min_af}, {max_af}]'))


    with pysam.VariantFile(args.in_VCF, 'r', drop_samples = args.drop_gt) as vcf_in, pysam.VariantFile(args.out_VCF, 'w') as vcf_out:
        # check if requested INFO field is present
        if not args.af_key in vcf_in.header.info:
            raise Exception(f'No meta-information entry about {args.af_key} INFO field found!')
        
        # add new meta-information line to the VCF, describing new info field we will be computing
        vcf_in.header.info.add('AF_BIN', number = 'A', type = 'String', description = 'List of AF bins: one per alternate allele.')
       
        # copy all meta-information lines from input VCF header to output VCF header
        for record in vcf_in.header.records:
            vcf_out.header.add_record(record)

        # if output VCF will include sample-level GT information, then copy all samples from input VCF file
        if not args.drop_gt:
            for sample in vcf_in.header.samples:
                vcf_out.header.add_sample(sample)
        
        # read input VCF records
        for record_in in vcf_in:
            record_out = record_in.copy()

            if args.af_key not in record_out.info:
                continue

            # Assign MAF bin and save it to the INFO field
            assigned_af_bins = []
            for af in record_out.info[args.af_key]:
                assigned_af_bins.append(assign_bin(af_bins, af))
            record_out.info['AF_BIN'] = assigned_af_bins if assigned_af_bins else None

            # Write new VCF record
            vcf_out.write(record_out)
