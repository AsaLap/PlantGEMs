# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used for the preparation of the required files for the 
Pathway Tools software reconstruction and launching of reconstruction 
using mpwt package from AuReMe."""

import re
import numpy as np
import string
import random
import mpwt
import subprocess
import multiprocessing

import utils


def get_sequence_region(data, mRNA):
    """Function which browse the .gff file and get the name of each gene, 
    the position in the genome and each corresponding region and transcript(s).
    
    ARGS:
        data (file) -- the gff file to browse (already put in memory by the function "read_file").
        mRNA (bool) -- decides if the function must search the mRNA (True) line or CDS (False).
    RETURN:
        dicoRegions -- a dictionary containing all the gathered 
        informations (see pipelinePT() for the structure).
    """
    
    dicoRegions = {}
    transcript_found = False
    for line in data:
        ##Searching the gene's informations
        if "\tgene\t" in line:
            transcript_found = False
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
            if spl[6] == "+":
                dicoRegions[region][gene] = {"Start": spl[3], "End": spl[4], "Transcripts" : {}}
            else:
                dicoRegions[region][gene] = {"Start": spl[4], "End": spl[3], "Transcripts" : {}}
        ##Searching the transcript's information
        if mRNA:
            if not transcript_found and "RNA\t" in line:
                try:
                    transcript = re.search('(?<=Name=)\w+(\.\w+)*(\-\w+)*', line).group(0)
                    dicoRegions[region][gene]["Transcripts"][transcript] = []
                    transcript_found = True
                except AttributeError:
                    print("The mRNA has no attribute 'Name='...")
                    dicoRegions[region][gene]["Transcripts"]["None"] = []
        else:
        #In case the gff file needs to be looked at on the CDS 
        #and not the mRNA to corresponds to the TSV file
            if not transcript_found and "CDS\t" in line:
                try: #Searching for CDS ID instead of mRNA.
                    transcript = re.search('(?<=ID=)[CcDdSs]*[:-]*\w+(\.\w+)*', line).group(0)[4:]
                    dicoRegions[region][gene]["Transcripts"][transcript] = []
                    transcript_found = True
                except AttributeError:
                    print("The CDS has no attribute 'ID='...")
                    dicoRegions[region][gene]["Transcripts"]["None"] = []
        ##Searching the exon's information
        if "\tCDS\t" in line:
            spl = line.split("\t")
            dicoRegions[region][gene]["Transcripts"][transcript].append([int(spl[3]),int(spl[4])])
    return dicoRegions


def make_dat(WD, dicoRegions, TYPE):
    """Function to create the .dat file.
    
    ARGS:
        WD (str) -- the path to the working directory in which the file (.dat) will be saved.
        dicoRegions -- the dictionary containing the data to create the .dat file 
        (see pipelinePT() for the structure).
        TYPE (str) -- indication if the sequence of the organism are assembled 
        as chromosomes or contigs (or else, see Pathway Tools guide).
    """
    
    print("\nWARNING ! :\n - If there are circular chromosomes in your data, you have to manually",\
        "correct the field 'CIRCULAR?' in the .dat file by changing 'N' (no) with 'Y' (yes).\n")
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
    utils.write_file(WD, "genetic-elements" + ".dat", datFile)


def make_fsa(WD, fileFASTA, dicoRegions):
    """Function to make the .fsa files.
    
    ARGS:
        WD (str) -- the path to the working directory to save the files.
        fileFASTA (str) -- the path to the fasta file of the organism.
        dicoRegions -- the dictionary containing the data to create 
        those files (see pipelinePT() for the structure). 
    """
    
    with open(fileFASTA, "r") as file:
        fasta = file.read()
    fasta = fasta.split(">")
    fasta = list(filter(None, fasta))
    listRegions = list(dicoRegions.keys())
    for i in fasta:
        region = re.search("\w+(\.\w+)*(\-\w+)*", i).group(0)
        if region in listRegions:
            listRegions.remove(region)
            utils.write_file(WD, region + ".fsa", i)


def make_pf(WD, fileEggNOG, dicoRegions):
    """Function to make the .pf files.
    
    ARGS:
        WD (str) -- the path of the working directory to save the files.
        fileEggNOG (str) -- the path to the .tsv file from EggNOG with most information 
        for the .pf file.
        dicoRegions -- the dictionary containing the data to create 
        the files (see pipelinePT() for the structure).
    """
    
    tsv = utils.read_file(fileEggNOG)
    list_index = list(np.arange(0, len(tsv)))
    for region in dicoRegions.keys():
        subPf = []
        for gene in dicoRegions[region].keys():
            for transcript in dicoRegions[region][gene]["Transcripts"].keys():
                found = False
                for i in list_index:
                    if found:
                        break
                    if transcript in tsv[i]:
                        list_index.remove(i)
                        found = True
                        subPf.append(parse_eggNog(gene,
                                                  dicoRegions[region][gene]["Start"],
                                                  dicoRegions[region][gene]["End"],
                                                  dicoRegions[region][gene]["Transcripts"][transcript], 
                                                  tsv[i]))
        if subPf:
            f = open(WD + region + ".pf", "w")
            for i in subPf:
                for j in i:
                    f.write(j)
            f.close()


def make_tsv(WD, taxon_name_list):
    """Function to make the taxon_id.tsv file.
    
    ARGS:
        WD (str) -- the path where to store this file.
        taxon_name_list (list) -- the list containing the name of the organism and its taxon id. 
    """
    
    res = "species\ttaxon_id\n"
    for i in taxon_name_list:
        res += i[0] + "\t" + str(i[1]) + "\n"
    utils.write_file(WD, "taxon_id.tsv", res)


def parse_eggNog(id, start, end, exon_pos, line):
    """Sub-function of make_pf() to write the info in the correct order for each transcript.
    
    ARGS:
        id (str) -- the gene name for the transcript.
        start (int) -- the start position of the sequence.
        end (int) -- the end position of the sequence.
        exon_pos (list of int) -- list of position of the exons.
        line (str) -- the line corresponding to the transcript in the .tsv file.
    RETURN:
        info (str) -- a string with all the information and with the correct 
        page settings for the .pf file.
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
    for exon in exon_pos:
        info.append("CODING-SEGMENT\t" + str(exon[0]) + "-" + str(exon[1]) + "\n")
    if spl[6]:
        go = spl[6].split(",")
        for i in go:
            info.append("DBLINK\t" + i + "\n")
    info.append("//\n")
    return info


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
    utils.write_file(WD, "organism-params.dat", info)


def main(data):
    """The function to make all the pipeline working, from creation of 
    the files to Pathway Tools via mpwt, and merging models.
    
    ARGS:
        data (str) -- the path to a index.txt file containing the path 
        to each ini file for each organism to run.
    -- Parameters taken from an ini file:
        fileGFF (str) -- the name of the .gff file for the organism.
        fileFASTA (str) -- the name of the .fasta file for the organism.
        fileEggNOG (str) -- the name of the .tsv file from EggNOG for the organism.
        TYPE (str) -- indication if the sequences of the organism are assembled 
        as chromosomes or contigs (or else, see Pathway Tools user guide).
        taxon_ID (int) -- the NCBI taxon ID of the species.
        mRNA (bool) -- decides if the function must search the mRNA line (True) 
        or CDS line (False) to get the name of the transcript.
        name (str) -- the name of the species database.
    
    -- Structure of dicoRegions (created in this function):
    {Region name (str):
        {Gene name (str):
            {"Start": int, "End": int, "Transcripts":
                {Transcript name (str):[[begin, end] the transcript(s)'s CDS positions (list of list of int]}
        }
    }
    """
    
    index = utils.read_file(data)
    WD = index.pop(0).rstrip()
    WDfiles = WD + "files/"
    WDinput = WD + "input/"
    WDoutput = WD + "output/"
    WDlog = WD + "log/"
    subprocess.run(["mkdir", WDinput, WDoutput, WDlog])
    WDpt = mpwt.find_ptools_path() + "/pgdbs/user/"
    taxon_name_list = []
    cpu = 0
    for ini in index:
        if ini:
            cpu += 1
            ###Reading of the parameters of the organism
            parameters = utils.read_config(ini.rstrip())
            fileGFF = parameters["FILES"]["GFF"]
            fileFASTA = parameters["FILES"]["FASTA"]
            fileEggNOG = parameters["FILES"]["EGGNOG"]
            TYPE = parameters["INFO"]["TYPE"]
            taxon_ID = int(parameters["INFO"]["NCBI_TAXON_ID"])
            mRNA = parameters.getboolean("INFO","mRNA")
            name = parameters["INFO"]["DATABASE_NAME"]
            
            ###Keeping some information for later
            taxon_name_list.append([name, taxon_ID])
            
            ###Preparing the files
            subprocess.run(["mkdir", WDinput + name])
            WDorg = WDinput + name + "/"
            print("------\n" + name + "\n------")
            gffFile = utils.read_file(WDfiles + fileGFF)
            dicoRegions = get_sequence_region(gffFile, mRNA) 
            make_dat(WDorg, dicoRegions, TYPE)
            make_fsa(WDorg, WDfiles + fileFASTA, dicoRegions)
            make_pf(WDorg, WDfiles + fileEggNOG, dicoRegions)

    ###Creating the tsv file for the taxon IDs
    make_tsv(WDinput, taxon_name_list) #Quite unspecific, could be improved (type, codon table)
    print("------\nCreation of the files finished\n------")
    
    ##Counting the number of cpu to use
    if cpu <= multiprocessing.cpu_count() - 2:
        nb_cpu = cpu
    else:
        nb_cpu = multiprocessing.cpu_count() - 2
    print("Number of CPU used : ", nb_cpu)
    
    ###Starting the mpwt script
    mpwt.multiprocess_pwt(input_folder = WDinput, output_folder = WDoutput,
                        patho_inference = True, patho_hole_filler = False,
                        patho_operon_predictor = False, pathway_score = 1,
                        dat_creation = True, dat_extraction = True,
                        number_cpu = nb_cpu, size_reduction = False,
                        patho_log = WDlog, ignore_error = False,
                        taxon_file = True, verbose = True)
    # fusion()


if __name__=="__main__":
    main("/home/asa/INRAE/Work/mpwt/index.txt")