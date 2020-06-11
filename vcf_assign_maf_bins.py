import argparse
import pysam

#
# Note: python3 + pysam library are required
#
# Example:
#   python vcf_assign_maf_bins.py -i input.vcf.gz -o output.vcf.gz -b 0 0.001 0.005 0.01 0.05 0.1 0.5 -G -- bins into (0,0.001], (0.001,0.005], (0.005,0.01], (0.01,0.05], (0.05,0.1], (0.1,0.5] and drops genotype information
#

argparser = argparse.ArgumentParser(description = 'Assigns MAF bin to each alternate allele in VCF/BCF.')
argparser.add_argument('-i', '--in-vcf', metavar = 'file', dest = 'in_VCF', required = True, help = 'Input VCF/BCF file.')
argparser.add_argument('-f', '--info-field', metavar = 'name', dest = 'af_key', type = str, default = 'AF', help = 'Name of the INFO field with alternate allele frequency. Default is "AF"')
argparser.add_argument('-G', '--drop-genotypes', dest = 'drop_gt', action='store_true', help = 'Don\'t include genotype information in the output VCF/BCF file.')
argparser.add_argument('-b', '--maf-bins', metavar = 'number', dest = 'maf_bins', type = float, nargs = '+', default = [0, 0.01, 0.05, 0.1, 0.5], help = 'List of MAF bin boundaries e.g. 0, 0.1, 0.5 => (0, 0.1], (0.1, 0.5]. Default: 0, 0.01, 0.05, 0.1, 0.5.')
argparser.add_argument('-o', '--out-vcf', metavar = 'file', dest= 'out_VCF', required = True, help = 'Output VCF/BCF file.')

if __name__ == '__main__':
    args = argparser.parse_args()
    
    assert len(args.maf_bins) > 1

    args.maf_bins.sort()
    maf_bins = []
    for min_maf, max_maf in zip(args.maf_bins, args.maf_bins[1:]):
        maf_bins.append((min_maf, max_maf, f'({min_maf}, {max_maf}]'))

    with pysam.VariantFile(args.in_VCF, 'r', drop_samples = args.drop_gt) as vcf_in, pysam.VariantFile(args.out_VCF, 'w') as vcf_out:
        # check if requested INFO field is present
        if not args.af_key in vcf_in.header.info:
            raise Exception(f'No meta-information entry about {args.af_key} INFO field found!')
        
        # add new meta-information line to the VCF, describing new info field we will be computing
        vcf_in.header.info.add('MAF_BIN', number = 'A', type = 'String', description = 'List of AF bins: one per alternate allele.')
       
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
            assigned_maf_bins = []
            for af in record_out.info[args.af_key]:
                bin_assigned = False
                if af is not None:
                    maf = 1.0 - af if af > 0.5 else af
                    for min_maf, max_maf, bin_label in maf_bins:
                        if maf > min_maf and maf <= max_maf:
                            assigned_maf_bins.append(bin_label)
                            bin_assigned = True
                            break
                if not bin_assigned:
                    assigned_maf_bins.append(None)
            record_out.info['MAF_BIN'] = assigned_maf_bins if assigned_maf_bins else None

            # Write new VCF record
            vcf_out.write(record_out)
