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
    
    return (chr1, pos1-1, pos1,
            chr2, pos2-1, pos2,
            name, score, strand1, strand2, sv_type)

def chromosome(chrom):
    chrs={"chr1":1,"chr2":2,"chr3":3,"chr4":4,"chr5":5,"chr6":6,"chr7":7,"chr8":8,"chr9":9,"chr10":10,"chr11":11,"chr12":12,"chr13":13,"chr14":14,"chr15":15,"chr16":16,"chr17":17,"chr18":18,"chr19":19,"chr20":20,"chr21":21,"chr23":22,"chrX":"X","chrY":"Y","NC_012920.1":"MT"}
    ch=chrs[chrom]
    return ch


def convert(vcf_file):
    vcf = VariantFile(vcf_file)
    bed = pd.DataFrame([get_bedpe(sv) for sv in vcf.fetch()],
                         columns = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2', 'type'])
    bed['chr1'] = bed['chr1'].apply(chromosome)
    bed['chr2'] = bed['chr2'].apply(chromosome)
    return bed

if __name__ == "__main__":
    args = parser.parse_args()
    bedpe = convert(args.input)
    bedpe.to_csv(args.output, sep = '\t', header = None, index = None)