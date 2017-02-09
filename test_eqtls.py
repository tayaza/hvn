#!usr/bin/python
import csv
import os


def get_genes(genes_file):
    genes = []
    with open(genes_file) as gfile:
        reader = csv.reader(gfile, delimiter = '\t')
        for line in reader:
            genes.append(line[0])
    gfile.close()
    return genes

def get_eqtls(genes, obesity_eqtls, diabetes_eqtls):
    to_gtex_obesity = []
    to_gtex_diabetes = []

    with open(obesity_eqtls) as ofile:
        reader = csv.reader(ofile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            snp = line[0]
            e_gene = line[3]
            tissue = line[7]
            for gene in genes:
                if e_gene == gene:
                    to_gtex_obesity.append([snp,gene,tissue])
    ofile.close()
    with open(diabetes_eqtls) as dfile:
        reader = csv.reader(dfile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            snp = line[0]
            e_gene = line[3]
            tissue = line[7]
            for gene in genes:
                if e_gene == gene:
                    to_gtex_diabetes.append([snp,gene,tissue])
    dfile.close()

    gtex_ofile = open('../analysis/insulin_eqtls_obesity.txt', 'wb')
    writer = csv.writer(gtex_ofile, delimiter = ',')
    writer.writerows(to_gtex_obesity)
    gtex_ofile.close()

    gtex_dfile = open('../analysis/insulin_eqtls_diabetes.txt', 'wb')
    writer = csv.writer(gtex_dfile, delimiter = ',')
    writer.writerows(to_gtex_diabetes)
    gtex_dfile.close()


if __name__== '__main__':
    obesity_eqtls = '../obesity/results/codes3d_results/dhs_results/sig_SNP-gene_eqtls.txt'
    diabetes_eqtls = '../diabetes/results/dhs_results/sig_SNP-gene_eqtls.txt'
    genes_file = '../analysis/leptin_genes.txt'
    genes = get_genes(genes_file)
    get_eqtls(genes, obesity_eqtls, diabetes_eqtls)
