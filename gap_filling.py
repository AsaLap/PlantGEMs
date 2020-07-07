# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used for the gap filling of the previously reconstructed models."""

import subprocess
import string
import os
import re
from pyasp.term import *
from urllib.request import urlopen
from meneco import query, utils, sbml, run_meneco

import utils as tools


def clean(path):
    """Function to get rid of the chevron character and uppercase fasta files.
    
    ARGS:
        path (str) -- the path to the fasta file.
    RETURN:
        res (str) -- the string uppercased and without the '>' char.
    """
    
    file = tools.read_file(path)
    res = ""
    for line in file:
        if ">" not in line:
            res += line.rstrip().upper()
    return res


def count(file, letters):
    """Function to count the occurence of the 'letters' given as argument.
    
    ARGS:
        file (str) -- the path to the file to count.
        letters (list of str) -- the char/str to count.
    RETURN:
        res (dic) -- a dictionary, each key is a str counted and the value 
        the number of occurence.
    """
    
    res = {}
    for letter in letters:
        res[letter] = 0
        for i in file:
            if letter == i:
                res[letter] += 1
    return res


def stats(res):
    """Function that transforms the dictionary of results from the 'count()' 
    function into percentages, print and return them.
    
    ARGS:
        res -- the dictionary from 'count()' (see count() for the structure).
    RETURN:
        str_res (str) -- a string containing the computed information line by line. 
    """
    
    total = 0
    str_res = ""
    for key in res.keys():
        total += res[key]
    for key in res.keys():
        str_res += str(key) + "\t" + str(res[key]) + "\t" + str(round(res[key] * 100 / total, 2)) + "\n"
        print("%s : %.0f = %.2f percent" %(key, res[key], res[key] * 100 / total))
    return str_res


def make_info(WD, DNA_file, RNA_file, prot_file, name):
    """Function to organize all the previous ones to write a .txt file for an organism.
    ARGS:
        WD (str) -- the working directory.
        DNA_file (str) -- the name of the fasta file of the total DNA.
        RNA_file (str) -- the name of the fasta file of only the transcripts.
        prot_file (str) -- the name of the fasta file of the proteins.
        name (str) -- the name for the output .txt file.
    """
    
    letterProt = "ACDEFGHIKLMNPQRSTVWY*"
    letterDNA = "CGTAN"
    str_res = "DNA stats :\n"
    DNA_file_clean = clean(WD + DNA_file)
    resDNA = count(DNA_file_clean, letterDNA)
    str_res += stats(resDNA)
    RNA_file_clean = clean(WD + RNA_file)
    resRNA = count(RNA_file_clean, letterDNA)
    str_res += "RNA stats :\n"
    str_res += stats(resRNA)
    prot_file_clean = clean(WD + prot_file)
    resProt = count(prot_file_clean, letterProt)
    str_res += "Prot stats :\n"
    str_res += stats(resProt)
    tools.write_file(WD, name + "_stats.txt", str_res)
    print(str_res)


def clean_sbml(WD, name):
    """Function to get rid of specific character COBRA puts into its sbml models.
    
    ARGS:
        WD (str) -- the working directory.
        name (str) -- the name of the sbml model.
    """
    
    file = tools.read_file(WD + name)
    new_file = []
    meta = re.compile('(<species id="M_)')
    sp_ref = re.compile('(<speciesReference species="M_)')
    reac = re.compile('(id="R_)')
    for i in file:
        i = meta.sub('<species id="', i)
        i = sp_ref.sub('<speciesReference species="', i)
        if "<reaction" in i:
            i = reac.sub('id="', i)
        new_file.append(i)
    tools.write_file(WD, "clean_" + name, new_file)



if __name__=="__main__":
    # make_info("/home/asa/INRAE/Work/Gap_filling/",
    # "testDNA.faa",
    # "testRNA.faa",
    # "testProt.faa",
    # "Tomato")
    # make_info("/home/asa/INRAE/Work/Gap_filling/",
    # "S_lycopersicum_chromosomes.4.00.faa",
    # "ITAG4.0_CDS.fasta",
    # "ITAG4.0_proteins.fasta",
    # "Tomato")
    
    ###Tomato
    clean_sbml("/home/antoine/INRAE/Work/Gap_filling/", "TomatoFusion.sbml")
    result = run_meneco(draftnet="/home/antoine/INRAE/Work/Gap_filling/clean_TomatoFusion.sbml",
                seeds="/home/antoine/INRAE/Work/Gap_filling/seedsPlants.sbml",
                targets="/home/antoine/INRAE/Work/Gap_filling/targetsTomato.sbml",
                repairnet="/home/antoine/INRAE/Work/Gap_filling/metacyc.sbml",
                enumeration=False,
                json=False)
    print(result)
    
    ###Kiwi
    # clean_sbml("/home/antoine/INRAE/Work/Gap_filling/", "KiwiFusion.sbml")
    # result = run_meneco(draftnet="/home/antoine/INRAE/Work/Gap_filling/clean_KiwiFusion.sbml",
    #             seeds="/home/antoine/INRAE/Work/Gap_filling/seedsPlants.sbml",
    #             targets="/home/antoine/INRAE/Work/Gap_filling/targetsKiwi.sbml",
    #             repairnet="/home/antoine/INRAE/Work/Gap_filling/metacyc.sbml",
    #             enumeration=False,
    #             json=False)
    # print(result)
    
    ###Cucumber
    # clean_sbml("/home/antoine/INRAE/Work/Gap_filling/", "CucumberFusion.sbml")
    # result = run_meneco(draftnet="/home/antoine/INRAE/Work/Gap_filling/clean_CucumberFusion.sbml",
    #             seeds="/home/antoine/INRAE/Work/Gap_filling/seedsPlants.sbml",
    #             targets="/home/antoine/INRAE/Work/Gap_filling/targetsCucumber.sbml",
    #             repairnet="/home/antoine/INRAE/Work/Gap_filling/metacyc.sbml",
    #             enumeration=False,
    #             json=False)
    # print(result)
    
    ###Cherry
    # clean_sbml("/home/antoine/INRAE/Work/Gap_filling/", "CherryFusion.sbml")
    # result = run_meneco(draftnet="/home/antoine/INRAE/Work/Gap_filling/clean_CherryFusion.sbml",
    #             seeds="/home/antoine/INRAE/Work/Gap_filling/seedsPlants.sbml",
    #             targets="/home/antoine/INRAE/Work/Gap_filling/targetsCherry.sbml",
    #             repairnet="/home/antoine/INRAE/Work/Gap_filling/metacyc.sbml",
    #             enumeration=False,
    #             json=False)
    # print(result)
    
    ###Camelina
    # clean_sbml("/home/antoine/INRAE/Work/Gap_filling/", "CamelinaFusion.sbml")
    # result = run_meneco(draftnet="/home/antoine/INRAE/Work/Gap_filling/clean_CamelinaFusion.sbml",
    #             seeds="/home/antoine/INRAE/Work/Gap_filling/seedsPlants.sbml",
    #             targets="/home/antoine/INRAE/Work/Gap_filling/targetsCamelina.sbml",
    #             repairnet="/home/antoine/INRAE/Work/Gap_filling/metacyc.sbml",
    #             enumeration=False,
    #             json=False)
    # print(result)