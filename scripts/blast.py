""" Working with BLAST files from searching against the KEGG database.

This module contains functions for dealing with tab-seperated BLAST
files generated from searching against the KEGG database. Specifically,
the goal of these functions is to go from ORFs that
were blasted against KEGG to the abundance of KOs. Abundance information
is taken from bowtie files that were generated from mapping sample reads
against ORFs.
"""

import os
import re
import sqlite3

def get_orf_gene_dict(blast_file):
    """Docstring.
    

    """


    #Open the provided blast file and catch the potential IOError.
    try:
        blast = open(blast_file, 'rU')
    except IOError:
        print('Input Blast file could not be opened!')
    
    #The current variable is to track the current ORF id. The input
    #blast file could have multiple entries for a given ORF and we want
    #to only use the top hit.
    current = None
    count_top_hits = 0
    orf_gene_dict = {}
    for line in blast:
        strip_line = line.rstrip()
        blast_list = strip_line.split('\t')
        
        #Assert that the number of fields in the input BLAST file is correct. 
        assert len(blast_list) == 12, (
                '%s does not contain 12 tab seperated fields as'
                'expected' % blast_file)

        #The first entry in the blast_list is the ORF ID seperated by a
        #space from another identifier. The second entry is the gene ID.
        #Building a dictionary with ORF ID as the key and gene ID as
        #the value.
        if  blast_list[0] != current:
            count_top_hits += 1
            orf_list = blast_list[0].split()
            orf_gene_dict[orf_list[0]] = blast_list[1]
        
        current = blast_list[0]
    
    #This will print to screen the number of ORFs with hits.
    output = '\nTop hits in BLAST file: %d\n\n' % (count_top_hits)
    print(output)  
    return orf_gene_dict


def get_gene_ko_dict(gene_ko_mapping_file):
    """Docstring.


    """

    #Open the provided file that contains information regarding KO
    #annotations for each gene ID and catch potential IOError.
    try:
        gene_ko_map = open(gene_ko_mapping_file, 'rU')
    except IOError:
        print('Input gene-KO maping file could not be opened!')
    
    #The gene mapping file is a tab-sepearted file. The first field is
    #the KO ID and the second field in the gene ID. Going to build a 
    #dictionary in which the key in the gene ID and value is KO ID.
    gene_ko_dict = {}
    for line in gene_ko_map:
        strip_line = line.rstrip()
        gene_ko_list = strip_line.split('\t')

        #Assert that the number of fields in file is correct
        assert len(gene_ko_list) == 2 ,(
               '%s does not contain the expected number of fields,'
                'likely wrong input file' % gene_ko_mapping_file)

        #The KO ID has unnecessary info prior to the colon.
        ko_id_list = gene_ko_list[0].split(':')
        gene_ko_dict[gene_ko_list[1]] = ko_id_list[1]
    return gene_ko_dict

def count_sample_reads_bowtie(bowtie_file):
    """


    """

    #Open the provided bowtie file which should contain information
    #regarding how reads map to ORFs
    try:
        bowtie = open(bowtie_file, 'rU')
    except IOError:
        print('Input Bowtie file could not be opened!')
    read_count = 0
    map_count = 0
    for line in bowtie:
        strip_line = line.rstrip()
        #Lines in bowtie files beginning with '@' are not needed.
        if strip_line.startswith('@'):
            continue
        read_count += 1
        
        split_line = strip_line.split('\t')
        if (split_line[2] != "*"):
            map_count += 1

    return(read_count, map_count)

def get_orf_abundance_dict(bowtie_file):
    """

    """

    #Open the provided bowtie file which should contain information
    #regarding how reads map to ORFs
    try:
        bowtie = open(bowtie_file, 'rU')
    except IOError:
        print('Input Bowtie file could not be opened!')
    orf_abundance_dict = {}
    for line in bowtie:
        strip_line = line.rstrip()
        #Lines in bowtie files beginning with '@' are not needed.        
        if strip_line.startswith('@'):
            continue
        bowtie_list = strip_line.split('\t')
        if bowtie_list[2] != '*':
            if bowtie_list[2] in orf_abundance_dict:
                orf_abundance_dict[bowtie_list[2]] += 1
            else:
                orf_abundance_dict[bowtie_list[2]] = 1
    return orf_abundance_dict


def get_ko_abundance_dict(orf_gene_dict, gene_ko_dict,
                            orf_abundance_dict): 
    """

    """ 
    ko_abundance_dict = {}
    for orf, orf_ab in orf_abundance_dict.items():
        #Not all orfs will be in the orf_gene dict because not all orfs
        #had hits to the kegg db
        if orf in orf_gene_dict:
            gene = orf_gene_dict[orf] 
            #Not all genes will be in the dictionary because not all
            #genes are assigned KOs 
            if gene in gene_ko_dict:
                ko = gene_ko_dict[gene]
            
                if ko in ko_abundance_dict:
                    ko_abundance_dict[ko] = ko_abundance_dict[ko] + orf_ab
                else:
                    ko_abundance_dict[ko] = orf_ab
                
    
    return ko_abundance_dict

def get_viral_orf_taxa_dicts(orf_gene_dict, viral_db):
    """

    """ 

    try:
        db = open(viral_db, 'rU')
    except IOError:
        print('Input viral database could not be opened!')

    db_dict = {}
    for line in db:
        strip_line = line.rstrip()
        if strip_line.startswith('>'):
            space_split = strip_line.split()[0]
            blast_result = space_split[1:]
        
            taxa = strip_line[strip_line.rfind("[")+1:strip_line.rfind("]")]
            taxa_sub1 = re.sub(r'_+', ' ', taxa)
            taxa_sub2 = re.sub(r'-+', ' ', taxa_sub1)
            taxa_sub3 = re.sub(r'([.])(\d)', r' \2', taxa_sub2)
            taxa_sub4 = re.sub(r' MM 2014', r'', taxa_sub3)
            taxa_sub5 = re.sub(r'\/', r' ', taxa_sub4)
            taxa_sub6 = re.sub(r'\'', r'', taxa_sub5)
            taxa_sub7 = re.sub(r'Baboon endogenous virus M7', 'Baboon endogenous virus strain M7', taxa_sub6)

        
            final_taxa = taxa_sub7.strip()
            db_dict[blast_result] = final_taxa 
    
    family_dict = {}
    genus_dict = {}

    for orf, gene in orf_gene_dict.items():
        #make this an assertion in the future?
        taxa = db_dict[gene]
        conn = sqlite3.connect('biosqldb-sqlite.db')
        result = conn.execute('select * from taxon_name inner join taxon on taxon_name.taxon_id == taxon.taxon_id where taxon_name.name=?', (taxa,)).fetchall()
        if len(result) == 0:
            result = conn.execute('select * from taxon_name inner join taxon on taxon_name.taxon_id == taxon.taxon_id where upper(taxon_name.name)=upper(?)', (taxa,)).fetchall()
            name = result[0][1]
            parent_taxon = result[0][5]
            rank = result[0][6]
        else:
            name = result[0][1]
            parent_taxon = result[0][5]
            rank = result[0][6]
 
        if rank == 'family':
            family_dict[orf] = name            
                
        if rank == 'genus':
            genus_dict[orf] = name  

        while parent_taxon != 1:       
            result = conn.execute('select * from taxon_name inner join taxon on taxon_name.taxon_id == taxon.taxon_id where taxon.taxon_id=?', (parent_taxon,)).fetchone()
            rank = result[6]
            name = result[1]

            if rank == 'family':
                family_dict[orf] = name

            if rank == 'genus':
                genus_dict[orf] = name 
            
            parent_taxon = result[5]

    return(family_dict, genus_dict)


def setup_viral_taxa_sql_db():
    """


    """
    conn = sqlite3.connect('biosqldb-sqlite.db')
    sqlite_file = 'biosqldb-sqlite.sql'
    f = open(sqlite_file,'r')
    sql_contents = f.read()
    conn.executescript(sql_contents)
    os.system('perl load_ncbi_taxonomy.pl --dbname=biosqldb-sqlite.db --driver=SQLite --directory=taxdata/')

def get_sample_viral_taxa_abundances(
        orf_abundance_dict,viral_orf_fam_dict,viral_orf_genus_dict):
    """

    
    """
    
    sample_fam_dict = {}
    sample_genus_dict = {}
    for orf, orf_ab in orf_abundance_dict.items():
        #make this an assetion in the future? maybe not b/c not
        #every orf will have a hit
        if orf in viral_orf_fam_dict:
            fam = viral_orf_fam_dict[orf]
            if fam in sample_fam_dict:
                sample_fam_dict[fam] = sample_fam_dict[fam] + orf_ab
            else:
                sample_fam_dict[fam] = orf_ab

        if orf in viral_orf_genus_dict:
            genus = viral_orf_genus_dict[orf]
            if genus in sample_genus_dict:
                sample_genus_dict[genus] = sample_genus_dict[genus] + orf_ab
            else:
                sample_genus_dict[genus] = orf_ab

    return(sample_fam_dict, sample_genus_dict)
