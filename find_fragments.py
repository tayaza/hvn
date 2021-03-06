#!/usr/bin/python
from configobj import ConfigObj
from itertools import cycle
from sets import Set
from wikipathways_api_client import WikipathwaysApiClient
import argparse,ast,bisect,json,multiprocessing,os,pandas,pybedtools,re,requests,sqlite3,time

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import style
from matplotlib.ticker import FuncFormatter
import csv

def process_inputs(snp_genes,snp_database_fp,fragment_database_fp):
    print "Processing input..."
    snp_db = sqlite3.connect(snp_database_fp)
    snp_db.text_factory = str
    snp_index = snp_db.cursor()
    fragment_index_db = sqlite3.connect(fragment_database_fp)
    fragment_index_db.text_factory = str
    fragment_index = fragment_index_db.cursor()
    snps = {}
    for line in snp_genes:
        id = line[0]
	snp = None
	snp_index.execute("SELECT * FROM snps WHERE rsID=?",(id,))
	snp = snp_index.fetchone()
	if snp == None:
            print "Warning: %s does not exist in SNP database." % id
	else:
            #Query fragmentIndex to find out to which fragment the SNP belongs
            fragment_index.execute("SELECT fragment FROM fragments WHERE \
                                       chr=? AND start<=? AND end>=?",\
					   ["chr" + snp[1],snp[2],snp[2]])
	    snp_fragment_result = fragment_index.fetchone()
	    if snp_fragment_result == None:
                print "Warning: error retrieving SNP fragment for SNP " + snp
	    else:
                snps[snp[0]]={ "chr": snp[1], "locus": snp[2], \
				       "frag": snp_fragment_result[0] }
    #for snp in snps.keys():
    #    print snp, snps[snp]
    return snps

def find_interactions(snp_genes, snps, gene_fragments, hic_data_dir):
    print "Finding fragment interactions for ..."
    #Look for all interactions involving SNP fragments in the HiC databases
    interactions = {} #A mapping of each SNP to the fragments with which \
                      #     the fragment it is on interacts
    for snp in snps.keys():
        interactions[snp] = {}

    detail_file = open(output_fp + '/interactions_fragments.txt', 'wb')
    detail = csv.writer(detail_file, delimiter = '\t')
    detail.writerow(['snp', 'snp_chromosome', 'snp_locus', 'snp_fragment', \
                        'gene', 'gene_chromosome', 'gene_start', 'gene_stop',\
                        'cis', 'fragment', 'fragment_start', 'fragment_stop', \
                         'cell_line', '#replicates', '#interactions'])

    summary_file = open(output_fp + '/interactions_summary.txt', 'wb')
    summary = csv.writer(summary_file, delimiter = '\t')
    summary.writerow(['snp', 'gene', 'cis', 'cell_lines', '#fragments', \
                         'total_interactions'])

    genes = find_gene_loci(snp_genes, gene_database_fp)
    for pair in snp_genes:
        snp = pair[0]
        for gene in gene_fragments.keys():
            if gene == pair[1]:
                print '\t', snp, '\t', gene 
                interactions[snp] = {}
                interactions[snp][gene] = {}
                fragment_count = 0
                total_inter = 0
                cell_line_list = ''
                for fragment in gene_fragments[gene]:
                    gene_chr = fragment[0][3:]
                    frag_start = int(fragment[1])
                    frag_stop = int(fragment[2])
                    gene_frag = fragment[3]
                    cell_line_count = 0
                    fragment_inter = 0
                    for cell_line in os.listdir(hic_data_dir):
                        if cell_line != 'GM12878_TEST' and \
                                os.path.isdir(hic_data_dir + '/' + cell_line):
                            #print "\tSearching cell line " + cell_line
                            count_replicates = 0
                            cell_line_inter = 0
                            for replicate in os.listdir(hic_data_dir + '/' + cell_line):
                                if replicate.endswith(".db"):
                                    rep_db = sqlite3.connect(hic_data_dir + '/' + \
                                                                 cell_line + '/' + replicate)
                                    rep_db.text_factory = str
                                    rep_ints = rep_db.cursor()
                                    #print "\t\tSearching replicate " + replicate
                                    interactions[snp][cell_line] = Set([])
                                    #print "\t\t\tFinding interactions for " + snp
                                    count_inter = 0
                                    for interaction in rep_ints.execute("SELECT * FROM \
                                                 interactions WHERE chr2 =? AND fragment2 =? \
                                                 AND chr1=? AND fragment1=?", \
                                                 [gene_chr, gene_frag,snps[snp]["chr"],\
                                                 snps[snp]["frag"]]):
                                        interactions[snp][cell_line].add(interaction)
                                        count_inter += 1
                                    if count_inter > 0:
                                        count_replicates += 1
                                        cell_line_inter += count_inter
                            if count_replicates > 0:
                                interactions[snp][gene][gene_frag] = {'replicates':count_replicates, \
                                                                      'interactions':cell_line_inter}
                                to_file = snp, snps[snp]['chr'], snps[snp]['locus'], \
                                    snps[snp]['frag'], gene, genes[gene]['chr'], \
                                    genes[gene]['start'], genes[gene]['stop'], pair[2], \
                                    gene_frag, frag_start, frag_stop, cell_line, \
                                    count_replicates, cell_line_inter
                                detail.writerow(to_file)
                                cell_line_count += 1
                                fragment_inter += cell_line_inter
                                if cell_line_list == '':
                                    cell_line_list = cell_line
                                else:
                                    if cell_line not in cell_line_list:
                                        cell_line_list = cell_line_list + ', ' + cell_line

                    if fragment_inter > 0:
                        total_inter += fragment_inter
                        fragment_count += 1
                summary.writerow([snp, gene, pair[2], cell_line_list, fragment_count, \
                                      total_inter])
    detail_file.close()
    summary_file.close()

def find_gene_loci(snp_genes, gene_database_fp):
    gene_db = sqlite3.connect(gene_database_fp)
    gene_db.text_factory = str
    gene_index = gene_db.cursor()
    gene_loci = {}
    for line in snp_genes:
        gene = line[1]
	#gene_index.execute('PRAGMA table_info(genes)')
	
	gene_index.execute('SELECT * FROM genes WHERE symbol = ?', \
				   (gene,))
	gene_data = gene_index.fetchone()
        gene_loci[gene] = {'chr':gene_data[1], 'start':gene_data[2], 'stop':gene_data[3]}

    return gene_loci


def find_gene_fragments(snp_genes, gene_database_fp, fragment_database_fp):
    fragment_db = sqlite3.connect(fragment_database_fp)
    fragment_db.text_factory = str
    fragment_index = fragment_db.cursor()
    genes = find_gene_loci(snp_genes, gene_database_fp)
    gene_fragments = {}
    for gene in genes.keys():
        gene_chr = genes[gene]['chr']
        gene_start = genes[gene]['start']
        gene_end = genes[gene]['stop']
        fragments = []
        """
        fragment_index.execute('SELECT * FROM fragments WHERE chr=? AND ' + \
                                   '((start <=? AND end <?) ' + \
                                   'OR (start <=? AND end >?)' + \
                                   'OR (start >? AND end <?)' + \
                                   'OR (start >? AND end <=?)' + \
                                   'OR (start >? AND end >?))', \
                                   ('chr'+gene_chr, \
                                        gene_start, gene_end, \
                                        gene_start, gene_end, \
                                        gene_start, gene_end, \
                                        gene_start, gene_end, \
                                        gene_start, gene_end))
        """
        fragment_index.execute('SELECT * FROM fragments WHERE chr=? AND ' + \
                                   'start >=? AND end <=?',  \
                                   ('chr'+gene_chr, \
                                        gene_start, gene_end)) 
        fragments = fragment_index.fetchall()
        fragment_index.execute('SELECT * FROM fragments WHERE chr=? AND ' + \
                                   'start <=? AND end >=?',  \
                                   ('chr'+gene_chr, \
                                        gene_start, gene_end)) 
        more_fragments = fragment_index.fetchall()
        for fragment in more_fragments:
            if fragment not in fragments:
                fragments.append(fragment)
        fragment_index.execute('SELECT * FROM fragments WHERE chr=? AND ' + \
                                   'start <? AND end >?',  \
                                   ('chr'+gene_chr, gene_start, gene_start)) 
        start_fragment = fragment_index.fetchone()
        if start_fragment is not None and start_fragment not in fragments:
                fragments.append(start_fragment)
        fragment_index.execute('SELECT * FROM fragments WHERE chr=? AND ' + \
                                   'start <? AND end >?',  \
                                   ('chr'+gene_chr, gene_end, gene_end)) 
        end_fragment = fragment_index.fetchone()
        if end_fragment is not None and end_fragment not in fragments:
                fragments.append(end_fragment)
        gene_fragments[gene] = fragments

    return gene_fragments
	
def parse_matchfile(matchfile):
    snp_genes =[]
    with open(matchfile, 'rb') as mfile:
        reader = csv.reader(mfile, delimiter = '\t')
	next(reader, None)
	for line in reader:
           snp = line[0]
	   gene = line[4]
           cis = line[6]
	   snp_genes.append([snp, gene, cis])
    mfile.close()

    return snp_genes


def set_filepath(input_fp, output):
    """ Set output directory filepath."""
    filepath = ''
    if output == 'default':
        filepath = input_fp
        while not (filepath.endswith('/')):
            filepath = filepath[:len(filepath)-1]

    else:
        filepath = output
        if not filepath.endswith('/'):
            filepath = filepath + '/'
    output_fp = filepath + 'interactions'
    if not os.path.isdir(output_fp):
        os.mkdir(output_fp)
    print output_fp
    return output_fp
                                                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', '--input', required = True, \
                            help = 'The \'match.txt\' file in \'dhs_results\'')
    parser.add_argument('-o', '--output', default = 'default', \
                            help = 'Directory to save result files')
    args = parser.parse_args()
    config_file = '/mnt/3dgenome/projects/tfad334/codes3d/docs/conf.py'
    config = ConfigObj(config_file)
    snp_database_fp = config["SNP_DATABASE_FP"]
    hic_data_dir = config["HIC_DATA_DIR"]
    fragment_bed_fp = config["FRAGMENT_BED_FP"]
    fragment_database_fp = config["FRAGMENT_DATABASE_FP"]
    gene_bed_fp = config["GENE_BED_FP"]
    gene_database_fp = config["GENE_DATABASE_FP"]
    eqtl_data_dir = config["EQTL_DATA_DIR"]
    expression_table_fp = config["EXPRESSION_TABLE_FP"]
    
    matchfile = args.input
    output_fp = set_filepath(args.input, args.output)
    snp_genes = parse_matchfile(matchfile)
    snp_fragments = process_inputs(snp_genes,snp_database_fp,fragment_database_fp)
    gene_fragments = find_gene_fragments(snp_genes, gene_database_fp, fragment_database_fp)
    interactions = find_interactions(snp_genes, snp_fragments,gene_fragments,hic_data_dir)
