#!usr/bin/python

import sys
import os
import csv
import argparse

def get_common_snps(bedfile1, bedfile2):
    trait1_snps = []
    trait2_snps = []
    common_snps = []
    if os.path.isfile(bedfile1):
        with open(bedfile1, 'rb') as bed1:
            reader = csv.reader(bed1, delimiter = '\t')
            next(reader, None)
            for line  in reader:
                snp = line[3]
                if snp not in trait1_snps:
                    trait1_snps.append(snp)
        bed1.close()
    if os.path.isfile(bedfile2):
        with open(bedfile2, 'rb') as bed2:
            reader = csv.reader(bed2, delimiter = '\t')
            next(reader, None)
            for line  in reader:
                snp = line[3]
                if snp not in trait2_snps:
                    trait2_snps.append(snp)
        bed2.close()
    for snp1 in trait1_snps:
        for snp2 in trait2_snps:
            if snp1 == snp2:
                common_snps.append(snp1)
    print len(trait1_snps), len(trait2_snps)
    with open(output_fp + '/t2d_snps.txt', 'wb') as afile:
        writer = csv.writer(afile)
        for snp in trait1_snps:
            afile.write(snp + '\n')
    return common_snps


def get_gwas_details(snps, gwas1, gwas2):
    gwas_details = {}
    for snp in snps:
        if os.path.isfile(gwas1):
            with open(gwas1, 'rb') as file1:
                reader = csv.reader(file1, delimiter = '\t')
                next(reader, None)
                for line in reader:
                    file_snp = line[21]
                    if snp == file_snp:
                        to_dict = []
                        pubmed_id = line[1]
                        author = line[2]
                        year = line[3][:4]
                        author = author + ', ' + year
                        trait = line[7]
                        reported_genes = line[13]
                        mapped_genes = line[14]
                        context = line[24]
                        details = (snp, context, trait, mapped_genes, author, \
                                       pubmed_id)
                        if snp not in gwas_details.keys():
                            to_dict.append(details)
                        else:
                            same = False
                            to_dict = gwas_details[snp]
                            for row in to_dict:
                                if pubmed_id == row[5]:
                                    same = True
                            if same == False:
                                to_dict.append(details)
                        gwas_details[snp] = to_dict
            file1.close()
        if os.path.isfile(gwas2):
            with open(gwas2, 'rb') as file2:
                reader = csv.reader(file2, delimiter = '\t')
                next(reader, None)
                for line in reader:
                    file_snp = line[21]
                    if snp == file_snp:
                        to_dict = []
                        pubmed_id = line[1]
                        author = line[2]
                        year = line[3][:4]
                        trait = line[7]
                        reported_genes = line[13]
                        mapped_genes = line[14]
                        context = line[24]
                        details = (snp, context, trait, mapped_genes, author, \
                                       pubmed_id)
                        if snp not in gwas_details.keys():
                            to_dict.append(details)
                        else:
                            same = False
                            to_dict = gwas_details[snp]
                            for row in to_dict:
                                if pubmed_id == row[5]:
                                    same = True
                            if same == False:
                                to_dict.append(details)
                        gwas_details[snp] = to_dict
            file2.close()
    with open(output_fp + '/common_gwas_snps.txt', 'wb') as outfile:
        writer = csv.writer(outfile, delimiter = '\t')
        file_head = 'SNP', 'CONTEXT', 'TRAIT', 'MAPPED_GENES', 'AUTHOR', \
            'PUBMED_ID'
        writer.writerow(file_head)
        for snp in gwas_details:
            for row in gwas_details[snp]:
                to_file = row[0], row[1], row[2], row[3], row[4], row[5]
                writer.writerow(to_file)
    print len(gwas_details.keys())
                            
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--gwas1', required = True,
                        help = 'Filepath of GWAS associations of first disease')
    parser.add_argument('-b', '--gwas2', required = True,
                        help = 'Filepath of GWAS associations of second disease')
    parser.add_argument('-i', '--bed1', required = True,
                        help = 'Filepath of SNP bed file of first disease')
    parser.add_argument('-j', '--bed2', required = True,
                        help = 'Filepath of SNP bed file of second disease')

    args = parser.parse_args()
    output_fp = '../analysis'
    common_snps = get_common_snps(args.bed1, args.bed2)
    #get_gwas_details(common_snps, args.gwas1, args.gwas2)
    
