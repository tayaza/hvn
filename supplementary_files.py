#!usr/bin/python
import csv
import os


def process():
    eqtls = []
    gwas_pairs = []
    gwas_record = []
    novel_record = []
    cis_record = []
    trans_record = []

    with open(gwas_lined) as gfile:
        reader = csv.reader(gfile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            gwas_pairs.append(line)
    gfile.close()

    with open(eqtl_file) as efile:
        reader = csv.reader(efile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            snp = line[0]
            gene = line[3]
            cis = line[12]
            pair = gene,  snp
            for row in gwas_pairs:
                if row[0] == gene and row[1]==snp:
                    if line not in gwas_record:
                        gwas_record.append(line)
            if cis == 'True':
                cis_record.append(line)
            else:
                trans_record.append(line)
            eqtls.append(line)
    for row in eqtls:
        if row not in gwas_record:
            novel_record.append(row)

    header = 'SNP', 'SNP_Chromosome', 'SNP-Locus', 'Gene_Name', \
        'Gene_Chromosome', 'Gene_Start', 'Gene_End', 'Tissue', 'p-value', \
        'q-value', 'Cell_Lines', 'GTEx_cis_p_Threshold', \
        'cis_SNP_Gene_Interaction', 'SNP-gene_Distance', \
        'Expression_Level_in_eQTL_Tissue', 'Max_Expressed_Tissue', \
        'Maximum_Expression_Level', 'Min_Expressed_Tissue', \
        'Min_Expressed_Level'

    gwas_file = open(output_dir + '/gwas_predicted_interactions.txt', 'wb')
    gwas_filer = csv.writer(gwas_file, delimiter = '\t')
    gwas_filer.writerow(header)
    gwas_filer.writerows(gwas_record)
    gwas_file.close()
    
    novel_file = open(output_dir + '/novel_interactions.txt', 'wb')
    novel_filer = csv.writer(novel_file, delimiter = '\t')
    novel_filer.writerow(header)
    novel_filer.writerows(novel_record)
    novel_file.close()

    trans_file = open(output_dir + '/trans_interactions.txt', 'wb')
    trans_filer = csv.writer(trans_file, delimiter = '\t')
    trans_filer.writerow(header)
    trans_filer.writerows(trans_record)
    trans_file.close()

    cis_file = open(output_dir + '/cis_interactions.txt', 'wb')
    cis_filer = csv.writer(cis_file, delimiter = '\t')
    cis_filer.writerow(header)
    cis_filer.writerows(cis_record)
    cis_file.close()


if __name__== '__main__':
    results_dir = '../obesity/results/codes3d_results/dhs_results'
    eqtl_file = results_dir + '/sig_SNP-gene_eqtls.txt'
    gwas_lined = results_dir + '/gwas_in_line.txt'
    novel_genes_file = '../analysis/codes3d_novel_associations.txt'
    output_dir = results_dir + '/supplementary_files'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    process()
