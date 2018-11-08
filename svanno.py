import pysam
import importlib
import sys
# sys.path.append('./')
import vcf2bedpe
importlib.reload(vcf2bedpe)
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True, help="input vcf")
parser.add_argument("-o", "--output", required=True, type=str, help="output bedpe")

def query_tab(tab, chr, start, end, col = 3):
    hit = next(tab.fetch(chr, start, end), '')
    if hit:
        return hit.split('\t')[col]
    else:
        return ''
    
def annotate(bedpe):
    
    bedpe = bedpe.copy()
    
    script_path = os.path.dirname(os.path.realpath(__file__))
    
    tab_files = [os.path.join(script_path, 'oncokb_genes.bed.gz'),
                 os.path.join(script_path, 'oncokb_exons.bed.gz'),
                 os.path.join(script_path, 'oncokb_utrs.bed.gz'),
                 os.path.join(script_path, 'oncokb_introns.bed.gz')]

    for feature, tab_file in zip(['gene', 'exon', 'utr', 'intron'], tab_files):
        tab = pysam.TabixFile(tab_file)
        bedpe['%s1' % feature] = bedpe.apply(lambda row: query_tab(tab, row['chr1'], row['start1'], row['end1']), axis = 1)
        bedpe['%s2' % feature] = bedpe.apply(lambda row: query_tab(tab, row['chr2'], row['start2'], row['end2']), axis = 1)
        tab.close()

    bedpe['gene1'] = bedpe['gene1'].apply(lambda x: x if x else 'Intergenic')
    bedpe['gene2'] = bedpe['gene2'].apply(lambda x: x if x else 'Intergenic')

    bedpe['exon1'] = bedpe['exon1'].apply(lambda x: 'exon' + x if x else '')
    bedpe['exon2'] = bedpe['exon2'].apply(lambda x: 'exon' + x if x else '')

    bedpe['intron1'] = bedpe['intron1'].apply(lambda x: 'intron' + x if x else '')
    bedpe['intron2'] = bedpe['intron2'].apply(lambda x: 'intron' + x if x else '')

    bedpe['utr1'] = bedpe['utr1'].str.replace('five_prime_', '5p').str.replace('three_prime_', '3p')
    bedpe['utr2'] = bedpe['utr2'].str.replace('five_prime_', '5p').str.replace('three_prime_', '3p')

    bedpe['annot'] = bedpe.apply(lambda row: '%s%s>>%s%s' % (row['gene1'] + ':' if row['gene1'] != 'Intergenic' and (row['exon1'] or row['intron1'] or row['utr1']) else row['gene1'],
                                                            row['utr1'] if row['utr1'] else row['exon1'] if row['exon1'] else row['intron1'],
                                                            row['gene2'] + ':' if row['gene2'] != 'Intergenic' and (row['exon2'] or row['intron2'] or row['utr2']) else row['gene2'],
                                                            row['utr2'] if row['utr2'] else row['exon2'] if row['exon2'] else row['intron2']), axis = 1)

    bedpe = bedpe[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2', 'type', 'annot']]
    
    return bedpe

if __name__ == "__main__":
    args = parser.parse_args()
    bedpe = annotate(vcf2bedpe.convert(args.input))
    bedpe.to_csv(args.output, index = False, header = True, sep = '\t')