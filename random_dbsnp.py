#!usr/bin/python

import csv
import random
import sqlite3

# Generate a list of random genomic positions
chr_file = './chromosomes.txt'
chr_list = {}
"""
with open(chr_file, 'rb') as cfile:
    reader = csv.reader(cfile, delimiter='\t')
    next(reader, None)
    for line in reader:
        chrom = line[0]
        if chrom != 'Total':
            chr_list[chrom] = line[4]
cfile.close()
"""
snp_list = {}
conn = sqlite3.connect('/mnt/3dgenome/projects/cekb635/codes3d/lib/snp_index_dbSNP_b147.db')
cur = conn.cursor()
conn.text_factory = str
cur.execute("SELECT COUNT(*) FROM snps")
snp_total = cur.fetchone()
print snp_total
rand_snps = []
for i in xrange(0,483):
    rand_id = random.randrange(1, int(snp_total[0]))
    cur.execute("SELECT * FROM snps WHERE rowid = ?;", (rand_id,))
    snp_data = cur.fetchone()
    print rand_id, snp_data
    rand_snps.append(snp_data[0])


with open('random_dbsnps.txt', 'wb') as rfile:
    #writer = csv.writer(rfile, delimiter = '\t')
    for snp in rand_snps:
        rfile.write(snp + '\n')
rfile.close
