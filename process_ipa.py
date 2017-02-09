#!usr/bin/python
import csv
import os
import argparse

"""
To process eQTL SNPs, genes and disease of biological functions
from IPA.
"""

def get_genes(ipa_file):
    gene_list = []
    with open(ipa_file) as ifile:
        reader = csv.reader(ifile, delimiter = '\t')
        next(reader, None)
        next(reader, None)
        next(reader, None)
        for line in reader:
            genes = line[3]
            genes = genes.replace(',', ' ')
            genes = genes.split(' ')
            for gene in genes:
                if gene not in gene_list:
                    gene_list.append(gene)
    ifile.close()
    return gene_list

def match_snp(genes, eqtl_file):
    snp_genes = []
    done = []
    efile = open(eqtl_file, 'rb')
    reader = csv.reader(efile, delimiter = '\t')
    next(reader, None)
    for line in reader:
        for gene in genes:
            snp = line[0]
            tester = snp+gene
            if gene == line[3]:
                if tester not in done:
                    done.append(tester)
                    snp_genes.append([snp, gene])

    efile.close()

    return snp_genes

def process(gene_list, diabetes_eqtls, obesity_eqtls):
    snp_genes = []
    diabetes = match_snp(gene_list, diabetes_eqtls)
    for row in diabetes:
        snp_genes.append([row[0], row[1], 'diabetes'])
    obesity = match_snp(gene_list, obesity_eqtls)
    for row in obesity:
        snp_genes.append([row[0], row[1], 'obesity'])
    
    print len(snp_genes), len(gene_list)
    done_file = args.input[:len(args.input)-4] + '_genes.txt'
    dfile = open(done_file, 'wb')
    writer = csv.writer(dfile, delimiter = '\t')
    writer.writerows(snp_genes)

if __name__== '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, \
                            help = "Exported file from IPA")

    diabetes_eqtls = '/mnt/3dgenome/projects/tfad334/hvn/diabetes/results/' + \
        'dhs_results/sig_SNP-gene_eqtls.txt'
    obesity_eqtls = '/mnt/3dgenome/projects/tfad334/hvn/obesity/results/' + \
        'codes3d_results/dhs_results/sig_SNP-gene_eqtls.txt'

    args = parser.parse_args()
    gene_list = get_genes(args.input)
    process(gene_list, diabetes_eqtls, obesity_eqtls)

