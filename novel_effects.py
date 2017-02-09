#!usr/bin/python
"""
Extracts eGenes with RPKM > 1.0 across liver, pancreas, subcutaneous and 
visceral adipose, and skeletal muscle
"""
import csv

def openfile(sfile):
    sfile = open(sfile, 'rb')
    s_reader = csv.reader(sfile, delimiter = '\t')
    next(s_reader, None)
    snp_list = []

    for row in s_reader:
        snp_list.append(row[0])
    return snp_list

def find_records(snp_list, gwas_file):
    gwas = open(gwas_file, 'rb')
    greader = csv.reader(gwas, delimiter = '\t')
    next(greader, None)
    snps = {}
    for row in greader:
        snp = row[21]
        genes = row[14]
        trait = row[7]
        author = row[2]
        pubdate = row[3]
        trait_desc = row[29].replace(')', '')
        trait_desc = trait_desc.replace('(', '')
        or_beta = row[30]
        ci = row[31]
        direction = ''
        disease = ''
        d = gwas_file.split('-')
        if 'obesity.tsv' in d:
            disease = 'obesity'
        if 'type_2_diabetes.tsv' in d:
            disease = 'diabetes'
        if 'increase' in ci:
            direction = 'increase'
        if 'decrease' in ci:
            direction = 'decrease'

        if snp not in snps.keys():
            to_snps = []
            to_snps.append([snp, genes, trait, author, pubdate, trait_desc, \
                                or_beta, ci, direction, disease])
            snps[snp] = to_snps
        else:
           to_snps = snps[snp]
           to_snps.append([snp, genes, trait, author, pubdate, trait_desc, \
                                or_beta, ci, direction, disease])
           snps[snp] = to_snps

    return snps

def writefile(snps, outfile):
    writer = csv.writer(outfile, delimiter = '\t')
    writer.writerow(['SNP', 'Genes', 'Trait', 'Author', 'Pub_Date', \
                         'Trait_desc', 'OR/Beta', 'Effect_direction', \
                         'Disease'])
    for snp in snps:
        for row in snps[snp]:
            if row[9] == 'diabetes':
                writer.writerow(row)
            if row[9] == 'obesity':
                writer.writerow(row)
            print row[9]
    
def extract_eqtls(genes, t2d_gwas, obesity_gwas):
    #t2d_snps = openfile(t2d_gwas)
    obesity_snps = openfile(obesity_gwas)
    #t2d_rec = find_records(t2d_snps, t2d_gwas)
    obesity_rec = find_records(obesity_snps, obesity_gwas)
    #t2d_genes = sort_genes(t2d_eqtls)
    #obesity_genes = sort_genes(obesity_eqtls)
    #t2d_records = sort_tissues(t2d_genes)
    #obesity_records = sort_tissues(obesity_genes)
    
    #t2d_outfile = open('../analysis/novel5_effect_t2d.txt', 'wb')
    obesity_outfile = open('../analysis/novel5_effect_obesity.txt', 'wb')
    #writefile(t2d_rec, t2d_outfile)
    writefile(obesity_rec, obesity_outfile)


if __name__== '__main__':
    t2d_gwas = '../diabetes/data/gwas-association-downloaded_2016-08-26-type_2_diabetes.tsv'
    obesity_gwas = '../obesity/data/gwas-association-downloaded_2016-07-13-obesity.tsv'
    novel_genes = '../analysis/novel_genes.txt'

    extract_eqtls(novel_genes, t2d_gwas, obesity_gwas)
    #openfile(t2d_eqtl_file)
