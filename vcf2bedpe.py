#!/ifs/work/leukgen/bin/python/.virtualenvs/isabel3.6/bin/python
import pysam
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
    strand1 = strand2 = '+'
    alt = sv.alts[0]
    bracket, chr2, pos2 = parse_alt(alt)
    pos2 = int(pos2)
    
    # zero based and one based??
    
    if alt.startswith(bracket):
        strand1 = '-'

    if bracket == '[':
        strand2 = '-'
        
    score = 100
    
    sv_type = sv.info['SVTYPE']
        
    name = '%s(%s:%s-%s:%s)' % (sv_type, chr1, pos1, chr2, pos2)
        
    # substruct 1 from first coord since VCF is 1 based while bedpe first position is 0 based

    return (chr1, pos1-1, pos1,
            chr2, pos2-1, pos2,
            name, score, strand1, strand2, sv_type)

def parse_alt(alt):
    result = re.findall(r'([][])(.+?)([][])', alt)
    bracket1, region, bracket2 = result[0]
    chr2, pos2 = region.rsplit(':', 1)
    pos2 = pos2
    return (bracket1, chr2, pos2)

def convert(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    bedpe = pd.DataFrame([get_bedpe(sv) for sv in vcf.fetch()],
                         columns = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2', 'type'])
    return bedpe

if __name__ == "__main__":
    args = parser.parse_args()
    bedpe = convert(args.input)
    bedpe.to_csv(args.output, sep = '\t', header = None, index = None)