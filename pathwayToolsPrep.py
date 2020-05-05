# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# 2020

import re

def read_file(path):
    f = open(path, "r")
    res = f.readlines()
    return res


def write_file(WD, filename, data):
    f = open(WD + filename, "w")
    for i in data:
        f.write(i)
    f.close()

def get_sequence_region(data):
    res = []
    for i in data:
        try:
            test = re.search('(?<=##sequence-region)[ \t]*\w+(\.\w+)*',i).group(0).strip()
            res.append(test)
        except AttributeError:
            pass
    return res


def parse_eggNog():
    print("TODO")


def make_pf():
    print("TODO")


def make_dat(WD, name, seq_regions):
    ###Type ??
    TYPE = ':CHRSM'
    ###Circular ??
    CIRC = 'N'

    datFile = []
    if TYPE != ":CONTIG":
        for i in seq_regions:
            datFile.append('ID\t%s\nNAME\t%s\nTYPE\t%s\nCIRCULAR?\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n'
                           %(i, i.upper(), TYPE, CIRC, WD + i + '.pf', WD + i + '.fsa'))
    write_file(WD, name, datFile)


# :CHRSM, :PLASMID, :MT (mitochondrial chromosome),:PT (chloroplast chromosome), or :CONTIG

# ID	TEST-CHROM-1
# NAME	Chromosome 1
# TYPE	:CHRSM
# CIRCULAR?	N
# ANNOT-FILE	chrom1.pf
# SEQ-FILE	chrom1.fsa
# //


if __name__=="__main__":
    ###Files and working directory###
    WDtom = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Tomato/'
    WDkiw = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Kiwi/'
    WDcuc = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Cucumber/'
    WDche = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Cherry/'
    WDcam = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Camelina/'
    
    ###Fasta files
    tomatoFasta = 'ITAG4.0_proteins.fasta'
    kiwiFasta = 'Actinidia_chinensis.Red5_PS1_1.69.0.pep.all.fa'
    cucumberFasta = 'Gy14_pep_v2.fa'
    cherryFasta = 'PRUAV_Regina_pep.fa'
    camelinaFasta = 'GCF_000633955.1_Cs_protein.faa'
    
    ###GFF files
    tomatoGFF = "ITAG4.0_gene_models.gff"
    kiwiGFF = "Actinidia_chinensis.Red5_PS1_1.69.0.44.gff3"
    cucumberGFF = "Gy14_gene_gff_v2.gff"
    cherryGFF = "PRUAV_Regina.gff3"
    camelinaGFF = "GCF_000633955.1_Cs_genomic.gff"
    
    ###eggNOG files
    tomatoEgg = "eggNOG_annotations.tsv"
    kiwiEgg = "eggNOG_annotations.tsv"
    cucumberEgg = "eggNOG_annotations.tsv"
    cherryEgg = "eggNOG_annotations.tsv"
    camelinaEgg = "eggNOG_annotations.tsv"
    
    ###Main###
    tomatoGFFfile = read_file(WDtom + tomatoGFF)
    tomatoSeqRegions = get_sequence_region(tomatoGFFfile)
    res = tomatoSeqRegions
    make_dat(WDtom, "Tomato.dat", res)
    # kiwiGFFfile = read_file(WDkiw + kiwiGFF)
    # kiwiSeqRegions = get_sequence_region(kiwiGFFfile)
    # print(kiwiSeqRegions)
    # get_sequence_region(read_file(WDche + cherryGFF))
    # get_sequence_region(read_file(WDcuc + cucumberGFF))
    # get_sequence_region(read_file(WDcam + camelinaGFF))