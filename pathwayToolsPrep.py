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
    regionsAndGenes = {}
    for i in data:
        try:
            region = re.search('(?<=##sequence-region)[ \t]*\w+(\.\w+)*',i).group(0).strip()
            regionsAndGenes[region] = []
        except AttributeError:
            pass
    get_regions_genes(data, regionsAndGenes)
    return regionsAndGenes


def get_regions_genes(data, dicoRegions):
    for i in data:
        for region in dicoRegions.keys():
            if region in i:
                try:
                    gene = re.search('(?<=ID=gene:)\w+(\.\w+)*', i).group(0)
                    dicoRegions[region].append(gene)
                except AttributeError:
                    pass
    print(dicoRegions.keys())


def parse_eggNog():
    print("TODO")


def make_pf():
    print("TODO")


def make_dat(WD, name, seq_regions, TYPE):
    """Function to create the .dat file.
    
    ARGS:
        WD --
        name --
        seq_regions --
        TYPE --
    RETURN:
    
    """
    
    print("\nWARNING ! :\n - If there are circular chromosomes in your data, you have to manually \
correct the field 'CIRCULAR?' in the .dat file by changing 'N' (no) with 'Y' (yes).\n")
    CIRC = 'N'
    datFile = []
    if TYPE == "NONE":
        for i in seq_regions:
            datFile.append('ID\t%s\nCIRCULAR?\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n'
                           %(i, CIRC, WD + i + '.pf', WD + i + '.fsa'))
    elif TYPE == ":CONTIG":
        for i in seq_regions:
            datFile.append('ID\t%s\nTYPE\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n'
                           %(i, TYPE, WD + i + '.pf', WD + i + '.fsa'))
    else:
        for i in seq_regions:
            datFile.append('ID\t%s\nTYPE\t%s\nCIRCULAR?\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n'
                           %(i, TYPE, CIRC, WD + i + '.pf', WD + i + '.fsa'))
    write_file(WD, name + ".dat", datFile)


def pipelinePT(WD, fileGFF, name, TYPE="NONE"):
    gffFile = read_file(WD + fileGFF)
    seqRegions = get_sequence_region(gffFile)
    make_dat(WD, name, seqRegions, TYPE)



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
    pipelinePT(WDtom, tomatoGFF, "Tomato", TYPE=":CHRSM")
    
    # kiwiGFFfile = read_file(WDkiw + kiwiGFF)
    # kiwiSeqRegions = get_sequence_region(kiwiGFFfile)
    # print(kiwiSeqRegions)
    
    # cherryGFFfile = read_file(WDche + cherryGFF)
    # cherrySeqRegions = get_sequence_region(cherryGFFfile)
    # make_dat(WDche, "Cherry.dat", cherrySeqRegions, TYPE=":CONTIG")
    
    # get_sequence_region(read_file(WDcuc + cucumberGFF))
    # get_sequence_region(read_file(WDcam + camelinaGFF))