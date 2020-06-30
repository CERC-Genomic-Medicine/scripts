import argparse
import pysam

#
# Note: python3 + pysam library are required
#
# Example:
#   python vcf_assign_depth_bins.py -i input.vcf.gz -o output.vcf.gz -b 0 1 5 10 20 30 -- bins into "<=0", "(0, 1]", "(1, 5]", "(5, 10]", "(10, 20]", "(20, 30]", ">30".
#

argparser = argparse.ArgumentParser(description = 'Assigns DP (depth) bin to each variant in VCF/BCF.')
argparser.add_argument('-i', '--in-vcf', metavar = 'file', dest = 'in_VCF', required = True, help = 'Input single sample VCF/BCF file.')
argparser.add_argument('-b', '--dp-bins', metavar = 'number', dest = 'dp_bins', type = float, nargs = '+', default = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], help = 'List of DP bin boundaries e.g. 0, 1, 2 => "<=0", "(0, 1]", "(1, 2]", ">2". Default: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10.')
argparser.add_argument('-o', '--out-vcf', metavar = 'file', dest= 'out_VCF', required = True, help = 'Output VCF/BCF file.')


def assign_bin(bins, value):
    if value is not None:
        for min_value, max_value, label in bins:
            if value > min_value and value <= max_value:
                return label
    return None


if __name__ == '__main__':
    args = argparser.parse_args()
    
    assert len(args.dp_bins) > 1

    args.dp_bins.sort()
    dp_bins = [ (float('-inf'), args.dp_bins[0], f'<={args.dp_bins[0]}') ]
    for min_dp, max_dp in zip(args.dp_bins, args.dp_bins[1:]):
        dp_bins.append((min_dp, max_dp, f'({min_dp}, {max_dp}]'))
    dp_bins.append((args.dp_bins[-1], float('inf'), f'>{args.dp_bins[-1]}'))

    with pysam.VariantFile(args.in_VCF, 'r') as vcf_in, pysam.VariantFile(args.out_VCF, 'w') as vcf_out:
        # check if requested DP FORMAT field is present
        if not 'DP' in vcf_in.header.formats:
            raise Exception(f'No meta-information entry about DP FORMAT field found!')

        # check if requested AD FORMAT field is present
        if not 'AD' in vcf_in.header.formats:
            raise Exception(f'No meta-information entry about AD FORMAT field found!')

        # add new meta-information line to the VCF, describing new info field we will be computing
        vcf_in.header.info.add('DP_BIN', number = '1', type = 'String', description = 'DP bin: one per variant.')
        vcf_in.header.info.add('AD_BIN', number = 'R', type = 'String', description = 'AD bin: one per allele (including reference).')
      
        # copy all meta-information lines from input VCF header to output VCF header
        for record in vcf_in.header.records:
            vcf_out.header.add_record(record)

        # copy sample from input VCF file
        for sample in vcf_in.header.samples:
            vcf_out.header.add_sample(sample)
        
        # read input VCF records
        for record_in in vcf_in:
            record_out = record_in.copy()

            dp = record_out.samples[0]['DP']
            record_out.info['DP_BIN'] = assign_bin(dp_bins, dp)

            ad = record_out.samples[0]['AD']
            assert len(record_out.alleles) == len(ad)
            ad_bins = []
            for allele_dp in ad:
                ad_bins.append(assign_bin(dp_bins, allele_dp))
            record_out.info['AD_BIN'] = ad_bins if ad_bins else None

            # Write new VCF record
            vcf_out.write(record_out)
