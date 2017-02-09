#!usr/bin/python

import csv
import random


# Generate a list of random genomic positions
chr_file = './chromosomes.txt'
chr_list = {}
with open(chr_file, 'rb') as cfile:
    reader = csv.reader(cfile, delimiter='\t')
    next(reader, None)
    for line in reader:
        chrom = line[0]
        if chrom != 'Total':
            chr_list[chrom] = line[4]
cfile.close()

rand_bases = []
for i in xrange(0,483):
    rand_chrom = random.choice(chr_list.keys())
    rand_base =  random.randrange(1, int(chr_list[rand_chrom]))
    rand_bases.append(['chr' + rand_chrom + ':' + str(rand_base)])

with open('random_bases.txt', 'wb') as rfile:
    writer = csv.writer(rfile, delimiter = '\t')
    #writer.writerow(['Chr', 'Base'])
    writer.writerows(rand_bases)

