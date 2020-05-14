# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# 2020
"""This file is used for the preparation of the required files for the 
Pathway Tools software reconstruction."""

import re
import numpy as np

def read_file(path):
    f = open(path, "r")
    res = f.readlines()
    f.close()
    return res


def write_file(WD, filename, data):
    f = open(WD + filename, "w")
    for i in data:
        f.write(i)
    f.close()


def get_sequence_region(data):
    dicoRegions = {}
    gene = ""
    for i in data:
        if "\tgene\t" in i:
            spl = i.split("\t")
            region = spl[0]
            gene = re.search('(?<=Name=)\w+(\.\w+)*', i).group(0)
            if region not in dicoRegions.keys():
                dicoRegions[region] = {}
            dicoRegions[region][gene] = {"Begin": spl[3], "End": spl[4], "Transcribes" : []}
        # if "\tmRNA\t" in i: ###TODO : on voulait prendre le nom des mRNA, mais marche pas pour Cameline
            # 
    return dicoRegions


def make_dat(WD, name, dicoRegions, TYPE):
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
        for i in dicoRegions.keys():
            datFile.append('ID\t%s\nCIRCULAR?\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n'
                           %(i, CIRC, WD + i + '.pf', WD + i + '.fsa'))
    elif TYPE == ":CONTIG":
        for i in dicoRegions.keys():
            datFile.append('ID\t%s\nTYPE\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n'
                           %(i, TYPE, WD + i + '.pf', WD + i + '.fsa'))
    else:
        for i in dicoRegions.keys():
            datFile.append('ID\t%s\nTYPE\t%s\nCIRCULAR?\t%s\nANNOT-FILE\t%s\nSEQ-FILE\t%s\n//\n'
                           %(i, TYPE, CIRC, WD + i + '.pf', WD + i + '.fsa'))
    write_file(WD, name + ".dat", datFile)


def make_fsa(WD, fileFASTA, dicoRegions):
    with open(WD + fileFASTA, "r") as file:
        fasta = file.read()
    fasta = fasta.split(">")
    fasta = list(filter(None, fasta))
    for i in fasta:
        region = re.search("\w+(\.\w+)*", i).group(0)
        listRegions = list(dicoRegions.keys())
        if region in listRegions:
            listRegions.remove(region)
            write_file(WD, region + ".fsa", i)


def make_pf(WD, fileEggNOG, dicoRegions):
    tsv = read_file(WD + fileEggNOG)
    list_index = list(np.arange(0, len(tsv)))
    for region in dicoRegions.keys():
        subPf = []
        for gene in dicoRegions[region].keys():
            for transcribe in dicoRegions[region][gene]["Transcribes"]:
                found = False
                for i in list_index:
                    if found:
                        break
                    if transcribe in tsv[i]:
                        list_index.remove(i)
                        found = True
                        ###TODO : parse_eggNog()
                        # subPf.append(parse_eggNog(i)) ## TODO


def parse_eggNog():
    print("TODO")


def pipelinePT(WD, fileGFF, fileFASTA, fileEggNOG, name, TYPE="NONE"):
    gffFile = read_file(WD + fileGFF)
    dicoRegions = get_sequence_region(gffFile)
    make_dat(WD, name, dicoRegions, TYPE)
    make_fsa(WD, fileFASTA, dicoRegions)
    # make_pf(WD, fileEggNOG, dicoRegions)


if __name__=="__main__":
    ###Files and working directory###
    WDtom = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Tomato/'
    WDkiw = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Kiwi/'
    WDcuc = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Cucumber/'
    WDche = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Cherry/'
    WDcam = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Camelina/'
    
    ###Fasta files
    tomatoFasta = 'S_lycopersicum_chromosomes.4.00.faa'
    kiwiFasta = 'Actinidia_chinensis.Red5_PS1_1.69.0.dna.toplevel.fa'
    cucumberFasta = 'Gy14_genome_v2.fa'
    cherryFasta = 'PRUAV_Regina.fa'
    camelinaFasta = 'GCF_000633955.1_Cs_genomic.fna'
    
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
    # pipelinePT(WDtom, tomatoGFF, tomatoFasta, tomatoEgg, "Tomato", TYPE=":CHRSM")
    # pipelinePT(WDkiw, kiwiGFF, kiwiFasta, kiwiEgg, "Kiwi", TYPE=":CONTIG")
    # pipelinePT(WDcuc, cucumberGFF, cucumberFasta, cucumberEgg, "Cucumber", TYPE=":CHRSM")
    # pipelinePT(WDche, cherryGFF, cherryFasta, cherryEgg, "Cherry", TYPE=":CONTIG")
    # pipelinePT(WDcam, camelinaGFF, camelinaFasta, camelinaEgg, "Camelina", TYPE=":CONTIG")