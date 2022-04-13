#!/usr/bin/env python3
### Objective : Comparing VCFs file (Target to Reference)
### All target variant are compared to the Reference (values, samples)
# Necessity : Equal versions (issue tracted, VCT)
# Requieres sorted VCF

import pysam
import argparse
import warnings
from datetime import datetime

argparser = argparse.ArgumentParser(description = 'Comparison of 2 VCFs')
argparser.add_argument('-R', '--Ref-vcf', metavar = 'file', dest = 'Rvcf', required = True, help = 'Reference Variant File, Origin VCFs')
argparser.add_argument('-T', '--Test-vcf', metavar = 'file', dest = 'Cvcf', required = True, help = 'File to compare')

def test_samples(R_VCF,C_VCF):
	if list(R_VCF.header.samples) == list(C_VCF.header.samples):
		print("OK .................... Sample are Equivalent")
		return list(C_VCF.header.samples)
	elif '_' in R_VCF.header.samples[1] and '_' not in C_VCF.header.samples[1]:
		warnings.warn('Delimiter might cause problem')
		return False
	elif '_' in C_VCF.header.samples[1] and '_' not in R_VCF.header.samples[1]:
		warnings.warn("Delimiter might cause problem")
		return False
	else :
		print('failed ........... Samples are not Equivalent')
		return False


def test_format(record_C,record_R):
	formats = set(record_R.format.keys()).intersection(record_C.format.keys())
	if not formats:
		print('failed ........... Reference format does not include test format')
	else :
		return formats


def test_values_Vequal(record_C,record_R):
	returned=True
	formats=test_format(record_C,record_R)
	for k in Samples:
		for q in formats:
			returned= returned and set(record_C.samples[k][q])==set(record_R.samples[k][q]) # none is a possible value thus sum is non-valid 
			if returned == False:
				return [k,q]
				break
	return returned


def test_Vequal(C_VCF,R_VCF):
###	timer=-1 #tracking information
###	counter=0
	REF=R_VCF.fetch()
	current=next(REF)
	for rec_C in C_VCF.fetch():
		equality=False
		while(current.pos<=rec_C.pos): # if Reference position is supperior to target we assume it is absent or files are not in order
			if current.alts == rec_C.alts and current.ref==rec_C.ref and current.pos==rec_C.pos:
				test=test_values_Vequal(rec_C,current)
				if test==True :
					equality=True
					break
				else :
					break      
			else:
				current=next(REF)
		if not equality and current.pos<=rec_C.pos: #broken out of loop without surpassing the target position
			print("problem Position ", rec_C.pos, current.pos,"at the", test[0],"samples, in the format", test[1])
			print("Values are",rec_C.samples[test[0]][test[1]], ' (target) and ', current.samples[test[0]][test[1]], "(Reference)")
			break
		elif not equality and current.pos>rec_C.pos: 
			print("problem Position ", rec_C.pos, "is absent in the Reference VCFs")
			break
##		if datetime.now().second!=timer: ##tracking information
#			timer=datetime.now().second
#			print('current position :', rec_C.pos, " processed per second: ", counter, end = "\r")
#			counter=0
#		else :
#			counter=counter+1
	print("Test Terminated")		




if __name__ == "__main__":
	args = argparser.parse_args()
	C = pysam.VariantFile(args.Cvcf,'r')
	R = pysam.VariantFile(args.Rvcf,'r')
	Samples=test_samples(C,R)
	if Samples:
		test_Vequal(C,R)
