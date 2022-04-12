#!/usr/bin/env python3
### Objective : Comparing VCFs file (Target to Reference)
### All target variant are compared to the Reference (values, samples)
# Necessity : Equal versions (issue tracted, VCT)
# Efficiency improvement welcomed


import pysam
import argparse
import warnings


argparser = argparse.ArgumentParser(description = 'Comparison of 2 VCFs')
argparser.add_argument('-R', '--Ref-vcf', metavar = 'file', dest = 'Rvcf', required = True, help = 'Reference Variant File, Origin VCFs')
argparser.add_argument('-T', '--Test-vcf', metavar = 'file', dest = 'Cvcf', required = True, help = 'File to compare')

def test_samples(R_VCF,C_VCF):
	if set(R_VCF.header.samples) == set(C_VCF.header.samples):
		print("OK .................... Sample are Equivalent")
		return set(C_VCF.header.formats.keys()).intersection(R_VCF.header.formats.keys())
	elif '_' in R_VCF.header.samples[1] and '_' not in C_VCF.header.samples[1]:
		warnings.warn('Delimiter might cause problem')
		return True
	elif '_' in C_VCF.header.samples[1] and '_' not in R_VCF.header.samples[1]:
		warnings.warn("Delimiter might cause problem")
		return True
	else :
		print('failed ........... Samples are not Equivalent')
		return False

def test_format(record_C,record_R):
	formats = set(record_R.format.keys()).intersection(record_C.format.keys())
	if not formats:
		print('failed ........... Reference format does not include test format')
	else :
		return formats


def test_Vequal(C_VCF,R_VCF):
	for rec_C in C_VCF.fetch():
		equality=False
    REF=R_VCF.fetch("chr" + rec_C.contig,rec_C.pos-1,rec_C.pos)       
		for rec_REF in REF:
			if rec_REF.alts == rec_C.alts and rec_REF.ref ==rec_C.ref :
				equality = equality or test_values_Vequal(rec_C,rec_REF)
			if equality:
				break
		if not equality:
			print("problem Position ", rec_C.pos)
			break
	print("Test Terminated")		

def test_values_Vequal(record_C,record_R):
	returned=True
	formats=test_format(record_C,record_R)
	for k in range(len(record_C.samples)):
		for q in formats:
			returned= returned and set(record_C.samples[k][q]) == set(record_R.samples[k][q])  # none is a possible value thus sum is non-valid 
			if returned == False:
				return returned
	return returned



if __name__ == "__main__":
	args = argparser.parse_args()
	C = pysam.VariantFile(args.Cvcf,'r')
	R = pysam.VariantFile(args.Rvcf,'r')
	if test_samples(C,R):
		test_Vequal(C,R)


