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
    return res


def write_file(WD, filename, data):
    f = open(WD + filename, "w")
    for i in data:
        f.write(i)
    f.close()


def get_sequence_region(data):
    dicoRegions = {}
    for i in data:
        try:
            region = re.search('(?<=##sequence-region)[ \t]*\w+(\.\w+)*',i).group(0).strip()
            dicoRegions[region] = {}
        except AttributeError:
            pass
    get_regions_genes(data, dicoRegions)
    return dicoRegions


def get_regions_genes(data, dicoRegions):
    ###TODO : optimize the reading of the file when many regions => contigs
    for i in data:
        for region in dicoRegions.keys():
            if region in i and "\tgene\t" in i:
                try:
                    gene = re.search('(?<=ID=)(gene:)*(gene-)*\w+(\.\w+)*', i).group(0)
                    if "gene" in gene:
                        gene = gene[5:]
                    spl = i.split("\t")
                    dicoRegions[region][gene] = {"Begin": spl[3], "End": spl[4], "Transcribes" : []}
                except AttributeError:
                    pass


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
    fasta = read_file(WD + fileFASTA)
    list_index = list(np.arange(0,len(fasta),2))
    for region in dicoRegions.keys():
        subFasta = []
        for gene in dicoRegions[region].keys():
            found = False
            for i in list_index:
                if found:
                    break
                if gene in fasta[i]:
                    subFasta.append(fasta[i])
                    subFasta.append(fasta[i + 1])
                    dicoRegions[region][gene]["Transcribes"].append(re.search('(?<=>)\w+(\.\w+)*', fasta[i]).group(0))
                    list_index.remove(i)
                    try:
                        if gene not in fasta[i + 2]:
                            found = True
                    except IndexError:
                        pass
        write_file(WD, region + ".fsa", subFasta)


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
                        ###TODO : parse_eggNog()
                        # subPf.append(parse_eggNog(i)) ## TODO
                        found = True


def parse_eggNog():
    print("TODO")


def pipelinePT(WD, fileGFF, fileFASTA, fileEggNOG, name, TYPE="NONE"):
    gffFile = read_file(WD + fileGFF)
    dicoRegions = get_sequence_region(gffFile)
    make_dat(WD, name, dicoRegions, TYPE)
    make_fsa(WD, fileFASTA, dicoRegions)
    make_pf(WD, fileEggNOG, dicoRegions)


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
    pipelinePT(WDtom, tomatoGFF, tomatoFasta, tomatoEgg, "Tomato", TYPE=":CHRSM")
    # pipelinePT(WDkiw, kiwiGFF, kiwiFasta, kiwiEgg, "Kiwi", TYPE=":CONTIG")
    # pipelinePT(WDcuc, cucumberGFF, cucumberFasta, cucumberEgg, "Cucumber", TYPE=":CHRSM")
    # pipelinePT(WDche, cherryGFF, cherryFasta, cherryEgg, "Cherry", TYPE=":CONTIG")
    ##Too much regions without any genes
    # pipelinePT(WDcam, camelinaGFF, camelinaFasta, camelinaEgg, "Camelina", TYPE=":CONTIG")