#!usr/bin/python
"""
Extracts eGenes with RPKM > 1.0 across liver, pancreas, subcutaneous and 
visceral adipose, and skeletal muscle
"""
import csv

def openfile(efile):
    efile = open(efile, 'rb')
    e_reader = csv.reader(efile, delimiter = '\t')
    next(e_reader, None)
    eqtls = {}
    tissues = ['Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 
               'Muscle_Skeletal', 'Liver', 'Pancreas']
    for tissue in tissues:
        eqtls[tissue] = {}
    ex_gene_file = open(genes_exclude, 'rb')
    ex_reader = csv.reader(ex_gene_file)
    exclude_genes = []
    for row in ex_reader:
        if row[0] not in exclude_genes:
            exclude_genes.append(row[0])


    for row in e_reader:
        if row[3] not in exclude_genes:
            if row[7] in tissues and float(row[14]) >= 1.0:
                if len(eqtls[row[7]]) == 0:
                    to_eqtls = []
                    to_eqtls.append(row)
                    eqtls[row[7]] = to_eqtls
                else:
                    to_eqtls = eqtls[row[7]]
                    to_eqtls.append(row)
                    eqtls[row[7]] = to_eqtls
    return eqtls

def sort_genes(dic):
    genes = {}
    for tissue in dic:
        for row in dic[tissue]:
            gene = row[3]
            snp = row[0]
            rpkm = row[14]
            if gene not in genes.keys():
                genes[gene] = {}
                genes[gene][snp] = {'Adipose_Subcutaneous':0.0, 
                               'Adipose_Visceral_Omentum':0.0, 
                               'Muscle_Skeletal':0.0, 
                               'Liver':0.0, 
                               'Pancreas':0.0}
                genes[gene][snp][tissue] = rpkm
            else:
                if snp not in genes[gene].keys():
                    genes[gene][snp] = {'Adipose_Subcutaneous':0.0, 
                                   'Adipose_Visceral_Omentum':0.0, 
                                   'Muscle_Skeletal':0.0, 
                                   'Liver':0.0, 
                                   'Pancreas':0.0}
                    genes[gene][snp][tissue] = rpkm
                else:
                    genes[gene][snp][tissue] = rpkm
    return genes
                
def sort_tissues(genes_dict):
    records = []
    for gene in genes_dict:
        for snp in genes_dict[gene]:
            #print gene, snp, genes_dict[gene][snp]
            row = [snp, \
                gene, \
                genes_dict[gene][snp]['Adipose_Subcutaneous'], \
                genes_dict[gene][snp]['Adipose_Visceral_Omentum'], \
                genes_dict[gene][snp]['Liver'], \
                genes_dict[gene][snp]['Muscle_Skeletal'], \
                genes_dict[gene][snp]['Pancreas']] 
            for i in xrange (0, len(row)):
                if row[i] == 0.0:
                    row[i] = ' '
            records.append(row)
        
    return records

def writefile(records, outfile):
    outwriter = csv.writer(outfile, delimiter = '\t')
    outwriter.writerow(('SNP', 'Gene', 'Adipose_Subcutaneous', \
                           'Adipose_Visceral_Omentum', 'Liver', \
                           'Muscle_Skeletal', 'Pancreas'))
    outwriter.writerows(records)
    
def extract_eqtls(t2d_file, obesity_file):
    t2d_eqtls = openfile(t2d_file)
    obesity_eqtls = openfile(obesity_file)
    t2d_genes = sort_genes(t2d_eqtls)
    obesity_genes = sort_genes(obesity_eqtls)
    t2d_records = sort_tissues(t2d_genes)
    obesity_records = sort_tissues(obesity_genes)
    
    t2d_outfile = open('../analysis/novel5_eqtls_t2d.txt', 'wb')
    obesity_outfile = open('../analysis/novel5_eqtls_obesity.txt', 'wb')
    writefile(t2d_records, t2d_outfile)
    writefile(obesity_records, obesity_outfile)


if __name__== '__main__':
    t2d_eqtl_file = '../diabetes/results/dhs_results/sig_SNP-gene_eqtls.txt'
    obesity_eqtl_file = '../obesity/results/codes3d_results/dhs_results/'
    obesity_eqtl_file += 'sig_SNP-gene_eqtls.txt'
    genes_exclude = '../analysis/lipid_genes_exclude.txt'

    extract_eqtls(t2d_eqtl_file, obesity_eqtl_file)
    #openfile(t2d_eqtl_file)
