#!usr/bin/python
import csv

def get_associations(novel):
    associations = []
    with open(novel, 'rb') as nfile:
        reader = csv.reader(nfile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            gene = line[0]
            snps = line[1]
            disease = line[2]
            cis = line[3]
            snps = snps.replace('[', '')
            snps = snps.replace(']', '')
            snps = snps.replace("'", "")
            snps = snps.split(',')
            for snp in snps:
                snp = snp.strip()
                associations.append((gene, snp, disease, cis))
    nfile.close()
    with open('../analysis/new_novel_associations.txt', 'wb') as wfile:
        writer = csv.writer(wfile, delimiter = '\t')
        writer.writerow(['Gene', 'SNP', 'Disease', 'Interaction'])
        writer.writerows(associations)
    wfile.close()
    return associations

def get_novel_eqtls(associations, diabetes_eqtls, obesity_eqtls):
    diabetes_tissues = {}
    obesity_tissues = {}
    diabetes_tissue_count = 0
    obesity_tissue_count = 0
    diabetes_proportions = []
    obesity_proportions = []
    exclude = ['Cells_Transformed_fibroblasts', \
                   'Cells_EBV-transformed_lymphocytes']
    for row in associations:
        disease = row[2]
        tissue = ''
        if disease == 't2d':
            tissue = get_tissue(row, diabetes_eqtls)
            diabetes_tissue_count += 1 
            if tissue not in exclude:
                if tissue not in diabetes_tissues.keys():
                    diabetes_tissues[tissue] = 0
                else:
                    diabetes_tissues[tissue] += 1
        elif disease == 'obesity':
            tissue = get_tissue(row, obesity_eqtls)
            obesity_tissue_count += 1
            if tissue not in exclude:
                if tissue not in obesity_tissues.keys():
                    obesity_tissues[tissue] = 0
                else:
                    obesity_tissues[tissue] += 1
    for tissue in diabetes_tissues.keys():
        proportion = diabetes_tissues[tissue] / float(diabetes_tissue_count)
        diabetes_proportions.append([tissue, diabetes_tissues[tissue], \
                                         proportion, 'diabetes'])
    for tissue in obesity_tissues.keys():
        proportion = obesity_tissues[tissue] / float(obesity_tissue_count)
        obesity_proportions.append([tissue, obesity_tissues[tissue], \
                                         proportion, 'obesity'])

    with open('../analysis/novel_eqtl_tissues.txt', 'wb') as tfile:
        writer = csv.writer(tfile, delimiter = '\t')
        writer.writerow(['tissue', '#eQTLs', 'percentage', 'disease'])
        writer.writerows(obesity_proportions)
        writer.writerows(diabetes_proportions)

def get_tissue(interaction, eqtl_file):
    gene = interaction[0]
    snp  = interaction[1]

    with open(eqtl_file, 'rb') as dfile:
        reader = csv.reader(dfile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            tissue = line[7]
            if snp == line[0] and gene == line[3]:
                return tissue
    dfile.close()

if __name__ == '__main__':
    novel = '../analysis/codes3d_novel_associations.txt'
    diabetes_eqtls = '../diabetes/results/dhs_results/sig_SNP-gene_eqtls.txt'
    obesity_eqtls = '../obesity/results/codes3d_results/sig_SNP-gene_eqtls.txt'
    associations = get_associations(novel)
    get_novel_eqtls(associations, diabetes_eqtls, obesity_eqtls)

