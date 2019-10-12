from pysam import VariantFile
import re
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True, help="input vcf")
parser.add_argument("-o", "--output", required=True, type=str, help="output bedpe")

def get_bedpe(sv):
    chr1 = sv.chrom
    pos1 = int(sv.pos)
    pos2= int(sv.stop)
    strand1 = strand2 = '+'
    sv_type = sv.info['SVTYPE']
    if sv_type == "TRA":
        chr2=sv.info["CHR2"]
    else:
        chr2=chr1
    score = 100
    name = '%s(%s:%s-%s:%s)' % (sv_type, chr1, pos1, chr2, pos2)
        
    # substruct 1 from first coord since VCF is 1 based while bedpe first position is 0 based

    return (chr1, pos1-1, pos1,
            chr2, pos2-1, pos2,
            name, score, strand1, strand2, sv_type)

def convert(vcf_file):
    vcf = VariantFile(vcf_file)
    bedpe = pd.DataFrame([get_bedpe(sv) for sv in vcf.fetch()],
                         columns = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2', 'type'])
    return bedpe

if __name__ == "__main__":
    args = parser.parse_args()
    bedpe = convert(args.input)
    bedpe.to_csv(args.output, sep = '\t', header = None, index = None)