# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# 2020
"""This file is used for the preparation of the required files for the 
Pathway Tools software reconstruction."""

import re
import numpy as np
import string
import random
import configparser


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


def read_config(ini):
    """Runs the config file containing all the information to make a new model.
    ARGS :
        ini -- the path to the .ini file.
    RETURN :
        config -- the configuration in a python dictionary.
    """
    config = configparser.ConfigParser()
    config.read(ini)
    return config


def get_sequence_region(data, mRNA):
    """Function which browse the .gff file and get the name of each gene, 
    the position in the genome and each corresponding region and transcribe(s).
    
    ARGS:
        data -- the gff file to browse (already put in memory by the function "read_file").
        mRNA -- a boolean to know if the function must search the mRNA line or CDS.
    RETURN:
        dicoRegions -- a dictionary containing all the gathered 
        informations with regions as principale key.
    """
    
    dicoRegions = {}
    CDSfound = False
    for line in data:
        if "\tgene\t" in line:
            CDSfound = False
            spl = line.split("\t")
            region = spl[0]
            try:
                gene = re.search('(?<=Name=)\w+(\.\w+)*(\-\w+)*', line).group(0)
            except AttributeError:
                try:
                    gene = re.search('(?<=ID=)(gene:)*\w+(\.\w+)*(\-\w+)*', line).group(0)
                    if "gene:" in gene:
                        gene = gene[5:]
                except AttributeError:
                    print("The gene name hasn't been found...")
                    gene = ""
                    pass
            if region not in dicoRegions.keys():
                dicoRegions[region] = {}
            dicoRegions[region][gene] = {"Start": spl[3], "End": spl[4], "Transcribes" : []}
        if mRNA:
            if "RNA\t" in line:
                try:
                    transcribe = re.search('(?<=Name=)\w+(\.\w+)*(\-\w+)*', line).group(0)
                    dicoRegions[region][gene]["Transcribes"].append(transcribe)
                except AttributeError:
                    print("The mRNA has no attribute 'Name='...")
                    dicoRegions[region][gene]["Transcribes"].append("None")
        else:
        #In case the gff file needs to be looked at on the CDS 
        #and not the mRNA to corresponds to the TSV file
            if CDSfound == False and "CDS\t" in line:
                try: #Searching for CDS ID instead of mRNA.
                    transcribe = re.search('(?<=ID=)[CcDdSs]*[:-]*\w+(\.\w+)*', line).group(0)[4:]
                    dicoRegions[region][gene]["Transcribes"].append(transcribe)
                    CDSfound = True
                except AttributeError:
                    print("The CDS has no attribute 'ID='...")
                    dicoRegions[region][gene]["Transcribes"].append("None")
            
    return dicoRegions


def make_dat(WD, dicoRegions, TYPE):
    """Function to create the .dat file.
    
    ARGS:
        WD -- the working directory in which the file will be saved.
        dicoRegions -- the dictionary containing the data to create the .dat file.
        TYPE -- indication if the sequence of the organism are assembled 
        as chromosomes or contigs (or else, see Pathway Tools guide).
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
    write_file(WD, "genetic-elements" + ".dat", datFile)


def make_fsa(WD, fileFASTA, dicoRegions):
    """Function to make the .fsa files.
    
    ARGS:
        WD -- the working directory to save those files.
        fileFASTA -- the fasta file of the organism where to get the sequences.
        dicoRegions -- the dictionary containing the data to create 
        those files, the region and the corresponding genes in it. 
    """
    
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
    """Function to make the .pf files.
    
    ARGS:
        WD -- the working directory to save those files.
        fileEggNOG -- the .tsv file from EggNOG with most information for the .pf file.
        dicoRegions -- the dictionary containing the data to create 
        those files, the region and the corresponding transcribe(s) in it.
    """
    
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
    
    
def make_organism_params(WD, species, abbrev, rank, storage = "file", private = "NIL", tax = 2, codon = 1, mito_codon = 1):
    #Choose tax = 1(4) for Bacteria, 2(5) for Eukaryota and 3(6) for Archae (2 is default).
    dico_tax = {1 : "TAX-2", 2 : "TAX-2759", 3 : "TAX-2157",
                4 : "2", 5 : "2759", 6 : "2157"}
    info = []
    #Making the random ID
    ID = random.choice(string.ascii_lowercase)
    string_choice = string.ascii_lowercase + "0123456789"
    for loop in range(random.randint(1,10)):
        ID += random.choice(string_choice)
    info.append("ID\t" + ID + "\n")
    info.append("STORAGE\t" + storage + "\n")
    info.append("NAME\t" + species + "\n")
    info.append("ABBREV-NAME\t" + abbrev + "\n")
    info.append("PRIVATE?\t" + private + "\n")
    info.append("RANK\t" + str(rank) + "\n")
    info.append("ORG-COUNTER\t\n") ##Test with or without
    info.append("DOMAIN\t" + dico_tax[tax] + "\n")
    info.append("CODON-TABLE\t" + str(codon) + "\n")
    info.append("MITO-CODON-TABLE\t" + str(mito_codon) + "\n")
    info.append("DBNAME\t" + abbrev + "DBcyc\n")
    info.append("NCBI-TAXON-ID\t" + dico_tax[tax+3] + "\n")
    write_file(WD, "organism-params.dat", info)


def parse_eggNog(id, start, end, line):
    """Sub-function of make_pf() to write the info in the correct order for each transcribe.
    
    ARGS:
        id -- the gene name for the transcribe.
        start -- the start position of the sequence.
        end -- the end position of the sequence.
        line -- the line corresponding to the transcribe of the .tsv file.
    RETURN:
        info -- a string with all the information and with the correct page settings.
    """
    
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
        for res in spl[7].split(","):
            info.append("EC\t" + res + "\n")
    if spl[6]:
        go = spl[6].split(",")
        for i in go:
            info.append("DBLINK\t" + i + "\n")
    info.append("//\n")
    return info


def pipelinePT(ini, TYPE="NONE", mRNA = True):
    """Function to run all the process.
    
    ARGS:
        ini -- the initialistion file containing all the following parameters:
            WD -- the working directory where to find/save the files.
            fileGFF -- the .gff file for the organism.
            fileFASTA -- the .fasta file for the organism.
            fileEggNOG -- the .tsv file for the organism from EggNOG.
        TYPE -- indication if the sequence of the organism are assembled 
        as chromosomes or contigs (or else, see Pathway Tools guide).
        mRNA -- True by default but false means the .gff file must be read by 
        the CDS and not the mRNA to get the transcribes name.
    """
    #Reading the parameters from the ini file
    param = read_config(ini)
    WD = param["PATH"]["DIRECTORY"]
    fileGFF = param["FILES"]["GFF"]
    fileFASTA = param["FILES"]["FASTA"]
    fileEggNOG = param["FILES"]["EGGNOG"]
    
    gffFile = read_file(WD + fileGFF)
    dicoRegions = get_sequence_region(gffFile, mRNA)
    make_dat(WD, dicoRegions, TYPE)
    make_fsa(WD, fileFASTA, dicoRegions)
    make_pf(WD, fileEggNOG, dicoRegions)


if __name__=="__main__":
    #Lauching the program for the 5 organism on which I'm working
    pipelinePT("/home/asa/INRAE/Work/PathwayToolsData/Tomato/TomatoAracycPT.ini", TYPE=":CHRSM")
    # pipelinePT("/home/asa/INRAE/Work/PathwayToolsData/Kiwi/KiwiAracycPT.ini", TYPE=":CHRSM")
    # pipelinePT("/home/asa/INRAE/Work/PathwayToolsData/Cucumber/CucumberAracycPT.ini", TYPE=":CONTIG")
    # pipelinePT("/home/asa/INRAE/Work/PathwayToolsData/Cherry/CherryAracycPT.ini", TYPE=":CONTIG", mRNA = False)
    # pipelinePT("/home/asa/INRAE/Work/PathwayToolsData/Camelina/CamelinaAracycPT.ini", TYPE=":CONTIG", mRNA = False)
    
    # make_organism_params("/home/asa/INRAE/Work/PathwayToolsData/Tomato/", "Solanum lycopersicum ITAG4.0", "Tomato", 195583)