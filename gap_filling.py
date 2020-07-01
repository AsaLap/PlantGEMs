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
from pyasp.term import *
from urllib.request import urlopen
from meneco import query, utils, sbml, run_meneco

# import utils


def clean(path):
    file = utils.read_file(path)
    res = ""
    for line in file:
        if ">" not in line:
            res += line.rstrip().upper()
    return res


def count(file, letters):
    res = {}
    for letter in letters:
        res[letter] = 0
        for i in file:
            if letter == i:
                res[letter] += 1
    return res


def stats(res):
    total = 0
    str_res = ""
    for key in res.keys():
        total += res[key]
    for key in res.keys():
        str_res += str(key) + "\t" + str(res[key]) + "\t" + str(round(res[key] * 100 / total, 2)) + "\n"
        print("%s : %.0f = %.2f percent" %(key, res[key], res[key] * 100 / total))
    return str_res


def make_info(WD, DNA_file, RNA_file, prot_file, name):
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
    utils.write_file(WD, name + "_stats.txt", str_res)
    print(str_res)


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
    
    # draft_sbml = utils.open_file("/home/asa/INRAE/Work/Gap_filling/Gap_filling_test/TomatoFusion.sbml")
    # draftnet = sbml.readSBMLnetwork(draft_sbml, 'draft')
    # seeds_sbml = urlopen("https://raw.githubusercontent.com/bioasp/meneco/master/Ectodata/seeds.sbml")
    # seeds = sbml.readSBMLseeds(seeds_sbml)
    # targets_sbml = utils.open_file("/home/asa/INRAE/Work/Gap_filling/Gap_filling_test/targetsTomato.sbml")
    # targets = sbml.readSBMLtargets(targets_sbml)
    # unproducible = TermSet()
    # model = query.get_unproducible(draftnet, targets, seeds)
    # for a in model:
    #     target = str(a)[13:]
    #     t      = String2TermSet(target)
    #     unproducible = unproducible.union(t)
    # unproducible = TermSet(unproducible)
    # utils.print_met(unproducible.to_list())
    # repair_sbml = utils.open_file("/home/asa/INRAE/Work/Gap_filling/Gap_filling_test/Metacyc.sbml")
    # repairnet = sbml.readSBMLnetwork(repair_sbml, 'repairnet')
    
    # draft_sbml= urlopen('https://raw.githubusercontent.com/bioasp/meneco/master/Ectodata/ectocyc.sbml')
    # draftnet = sbml.readSBMLnetwork(draft_sbml, 'draft') 
    # seeds_sbml = urlopen('https://raw.githubusercontent.com/bioasp/meneco/master/Ectodata/seeds.sbml')
    # seeds = sbml.readSBMLseeds(seeds_sbml)
    # targets_sbml = urlopen('https://raw.githubusercontent.com/bioasp/meneco/master/Ectodata/targets.sbml')
    # targets = sbml.readSBMLtargets(targets_sbml)
    # unproducible = TermSet()
    # model = query.get_unproducible(draftnet, targets, seeds)
    # for a in model:
    #     target = str(a)[13:]
    #     t = String2TermSet(target)
    #     unproducible = unproducible.union(t)
    # unproducible = TermSet(unproducible)
    # utils.print_met(unproducible.to_list())
    
    result = run_meneco(draftnet="/home/asa/INRAE/Work/Gap_filling/Gap_filling_test/draft.sbml",
                seeds="/home/asa/INRAE/Work/Gap_filling/Gap_filling_test/seeds.sbml",
                targets="/home/asa/INRAE/Work/Gap_filling/Gap_filling_test/targets.sbml",
                repairnet="/home/asa/INRAE/Work/Gap_filling/Gap_filling_test/Metacyc.sbml",
                enumeration=False,
                json=True)