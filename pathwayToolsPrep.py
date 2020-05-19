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
    for line in data:
        if "\tgene\t" in line:
            spl = line.split("\t")
            region = spl[0]
            try:
                gene = re.search('(?<=Name=)\w+(\.\w+)*(\-\w+)*', line).group(0)
            except AttributeError:
                print("The gene has non attribute 'Name='...")
                pass
            if region not in dicoRegions.keys():
                dicoRegions[region] = {}
            dicoRegions[region][gene] = {"Start": spl[3], "End": spl[4], "Transcribes" : []}
        if "RNA\t" in line:
            try:
                transcribe = re.search('(?<=Name=)\w+(\.\w+)*(\-\w+)*', line).group(0)
                dicoRegions[region][gene]["Transcribes"].append(transcribe)
            except AttributeError:
                print("The mRNA has no attribute 'Name='...")
                dicoRegions[region][gene]["Transcribes"].append("None")
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
        region = re.search("\w+(\.\w+)*(\-\w+)*", i).group(0)
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
                        subPf.append(parse_eggNog(gene,
                                                  dicoRegions[region][gene]["Start"],
                                                  dicoRegions[region][gene]["End"], 
                                                  tsv[i]))
        if subPf:
            f = open(WD + region + ".pf", "w")
            for i in subPf:
                for j in i:
                    f.write(j)
            f.close()


def parse_eggNog(id, start, end, line):
    info = []
    spl = line.split("\t")
    info.append("ID\t" + id + "\n")
    if spl[5]:
        info.append("NAME\t" + spl[5] + "\n")
    else:
        info.append("NAME\tORF\n")
    info.append("STARTBASE\t" + start + "\n")
    info.append("ENDBASE\t" + end + "\n")
    spl[21] = spl[21].replace("\n", "")
    if spl[21]:
        info.append("FUNCTION\t" + spl[21] + "\n")
    else:
        info.append("FUNCTION\tORF\n")
    info.append("PRODUCT-TYPE\tP\n")
    if spl[7]:
        info.append("EC\t" + spl[7] + "\n")
    if spl[6]:
        go = spl[6].split(",")
        for i in go:
            info.append("DBLINK\t" + i + "\n")
    info.append("//\n")
    return info


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
    # pipelinePT(WDcuc, cucumberGFF, cucumberFasta, cucumberEgg, "Cucumber", TYPE=":CHRSM")
    pipelinePT(WDche, cherryGFF, cherryFasta, cherryEgg, "Cherry", TYPE=":CONTIG")
    ##Don't work because don't have name or ID not corresponding to transcript
    # pipelinePT(WDkiw, kiwiGFF, kiwiFasta, kiwiEgg, "Kiwi", TYPE=":CONTIG")
    # pipelinePT(WDcam, camelinaGFF, camelinaFasta, camelinaEgg, "Camelina", TYPE=":CONTIG")