# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used for the gap filling of the previously reconstructed models."""

import cobra
import copy
import os
import re
import subprocess
import string
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


def add_filled_reactions(WD, reacs, repair, draft, json = False):
    """Function to add the retrieved reactions to the model.
    
    ARGS:
        WD (str) -- the path to the working directory where there are the repair
        model and the draft.
        reacs (list of str)-- the list of reaction names.
        repair -- the repair model (SBML format).
        draft -- the draft model (SBML format).
        json (boolean) -- True if you want a save of the model in the JSON
        format instead of SBML (default).
    """
    
    repair_model = cobra.io.read_sbml_model(WD + repair)
    draft_model = cobra.io.read_sbml_model(WD + draft)
    old_size = len(draft_model.reactions)
    for reac in reacs:
        reac = tools.cobra_compatibility(reac.rstrip())
        try:
            print(repair_model.reactions.get_by_id(reac))
            draft_model.add_reactions([repair_model.reactions.get_by_id(reac)])
        except KeyError:
            print("No match for this reaction : ", reac)
    print("Number of reactions of the unrepaired model : %i\nNumber of reactions of the repaired model : %i" %(old_size, len(draft_model.reactions)))
    if json:
        cobra.io.save_json_model(draft_model, WD + "filled_" + draft + ".json")
    else:
        cobra.io.write_sbml_model(draft_model, WD + "filled_" + draft + ".sbml")
    
    
def make_plantnetwork(WD, metacycPath, reactionsPath):
    """Function to create a subnetwork of Metacyc containing only the Viridiplantae taxon.
    
    ARGS:
        WD (str) -- the path to the working directory where to save the new model.
        metacycPath (str) -- the exact path and filename of the Metacyc model.
        reactionsPath (str) the exact path and filename of the reactions.dat file of the organism.  
    """
    reactionsFile = tools.read_file(reactionsPath)
    list_plant = []
    for line in reactionsFile:
        if "UNIQUE-ID" in line:
            try:
                unique_id = re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', line).group(0).rstrip()
            except AttributeError:
                print("No match for : ", line)
        if "TAXONOMIC-RANGE" in line and "33090" in line:
            list_plant.append(unique_id)
    
    ###Necessary step for the gathering of short/long IDs
    dico_matching, dico_matching_rev = tools.corres_dico("/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv")
    
    ###Transformation of all IDs in short + long to get them from Metacyc
    new_plant_list = []
    for reac in list_plant:
        try:
            for long_reac in dico_matching[reac]:
                new_plant_list.append(long_reac)
            new_plant_list.append(reac)
        except KeyError:
            try:
                new_plant_list.append(dico_matching_rev[reac])
                new_plant_list.append(reac)
            except KeyError:
                print("No match for reac : ", reac)
    
    ###Fetching the plant reactions
    metacyc = cobra.io.read_sbml_model(metacycPath)
    meta_plant = cobra.Model(name="Metacyc Plant")
    for reac in new_plant_list:
        try:
            plant_reac = copy.deepcopy(metacyc.reactions.get_by_id(reac))
            meta_plant.add_reactions([plant_reac])
        except KeyError:
            try:
                plant_reac = copy.deepcopy(metacyc.reactions.get_by_id(reac))
                meta_plant.add_reactions([plant_reac])
            except KeyError:
                print("Reaction not found in Metacyc : ", reac)
    cobra.io.write_sbml_model(meta_plant, WD + "meta_plant.sbml")


def pipeline_gap_filling(WD, draft, seeds, targets, repair, enumeration = False, json = False):
    """The main function to make all the process.
    
    ARGS:
        WD (str) -- the path to the directory where to find all the needed files.
        draft (str) -- the filename of the model to enhance.
        seeds (str) -- the filename of the seeds.
        targets (str) -- the filename of the targets.
        repair (str) -- the filename of the repair model.
        enumeration (boolean) -- Meneco choice to list all the reactions found or not.
        json (boolean) -- Meneco choice of getting the result as a JSON (not working).
    """
    
    clean_draft = tools.clean_sbml(WD, draft)
    clean_repair = tools.clean_sbml(WD, repair)
    result = run_meneco(draftnet = WD + clean_draft,
                        seeds = WD + seeds,
                        targets = WD + targets,
                        repairnet = WD + clean_repair,
                        enumeration = enumeration,
                        json = json)
    print(result)
    # add_filled_reactions(WD, result['Union of cardinality minimal completions'], repair, clean_draft)


if __name__=="__main__":  
    # make_info("/home/asa/INRAE/Work/Gap_filling/",
    # "S_lycopersicum_chromosomes.4.00.faa",
    # "ITAG4.0_CDS.fasta",
    # "ITAG4.0_proteins.fasta",
    # "Tomato")
    
    ##Tomato
    # pipeline_gap_filling("/home/asa/INRAE/Work/Gap_filling/",
    #                      "TomatoFusion.sbml",
    #                      "seedsPlants.sbml",
    #                      "targetsTomato.sbml",
    #                      "metacyc.sbml")
    
    ###Kiwi
    # pipeline_gap_filling("/home/asa/INRAE/Work/Gap_filling/",
    #                      "KiwiFusion.sbml",
    #                      "seedsPlants.sbml",
    #                      "targetsKiwi.sbml",
    #                      "metacyc.sbml")
    
    ###Cucumber
    # pipeline_gap_filling("/home/asa/INRAE/Work/Gap_filling/",
    #                      "CucumberFusion.sbml",
    #                      "seedsPlants.sbml",
    #                      "targetsCucumber.sbml",
    #                      "metacyc.sbml")
    
    ###Cherry
    # pipeline_gap_filling("/home/asa/INRAE/Work/Gap_filling/",
    #                      "CherryFusion.sbml",
    #                      "seedsPlants.sbml",
    #                      "targetsCherry.sbml",
    #                      "metacyc.sbml")
    
    ###Camelina
    # pipeline_gap_filling("/home/asa/INRAE/Work/Gap_filling/",
    #                      "CamelinaFusion.sbml",
    #                      "seedsPlants.sbml",
    #                      "targetsTomato.sbml",
    #                      "metacyc.sbml")
    
    # make_plantnetwork("/home/asa/INRAE/Work/Gap_filling/",
    #                 "/home/asa/INRAE/Logiciels/pathway-tools/metacyc.sbml",
    #                 "/home/asa/INRAE/Logiciels/pathway-tools/24.0/data/#reactions.dat#")
    
    ###Manual curation
    # pipeline_gap_filling("/home/asa/INRAE/Work/Gap_filling_tests/",
    #                      "clean_KiwiFusion.sbml",
    #                      "seedsPlantsCorrected.sbml",
    #                      "unrec_targetsKiwi.sbml",
    #                      "metacyc.sbml")
    
    ###Photosynthesis
    # pipeline_gap_filling("/home/asa/INRAE/Work/Gap_filling_tests/",
    #                     "clean_KiwiFusion.sbml",
    #                     "seeds_photo.sbml",
    #                     "targets_photo.sbml",
    #                     "metacyc.sbml")
    
    ###Adding manually curated reactions
    ###Tomato
    reacs = tools.read_file("/home/asa/INRAE/Work/Gap_filling/reactions_to_add_Tomato.txt")
    reacs = tools.trans_short_ID(reacs, "/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv")
    add_filled_reactions("/home/asa/INRAE/Work/Gap_filling/",
                         reacs,
                         "clean_metacyc.sbml",
                         "clean_FusionGenesTomato.sbml",
                         True)
    
    ###Kiwi
    reacs = tools.read_file("/home/asa/INRAE/Work/Gap_filling/reactions_to_add_Kiwi.txt")
    reacs = tools.trans_short_ID(reacs, "/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv")
    add_filled_reactions("/home/asa/INRAE/Work/Gap_filling/",
                         reacs,
                         "clean_metacyc.sbml",
                         "clean_FusionGenesKiwi.sbml",
                         True)
    
    ###Cucumber
    reacs = tools.read_file("/home/asa/INRAE/Work/Gap_filling/reactions_to_add_Cucumber.txt")
    reacs = tools.trans_short_ID(reacs, "/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv")
    add_filled_reactions("/home/asa/INRAE/Work/Gap_filling/",
                         reacs,
                         "clean_metacyc.sbml",
                         "clean_FusionGenesCucumber.sbml",
                         True)
    
    ###Cherry
    reacs = tools.read_file("/home/asa/INRAE/Work/Gap_filling/reactions_to_add_Cherry.txt")
    reacs = tools.trans_short_ID(reacs, "/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv")
    add_filled_reactions("/home/asa/INRAE/Work/Gap_filling/",
                         reacs,
                         "clean_metacyc.sbml",
                         "clean_FusionGenesCherry.sbml",
                         True)
    
    ###Camelina
    reacs = tools.read_file("/home/asa/INRAE/Work/Gap_filling/reactions_to_add_Camelina.txt")
    reacs = tools.trans_short_ID(reacs, "/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv")
    add_filled_reactions("/home/asa/INRAE/Work/Gap_filling/",
                         reacs,
                         "clean_metacyc.sbml",
                         "clean_FusionGenesCamelina.sbml",
                         True)