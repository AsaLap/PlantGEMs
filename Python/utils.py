# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file contains utility functions used in several scripts."""

import cobra
import configparser
import copy
import csv
import json
import re


def read_file(path):
    """Function to read and return a file line by line in a list."""
    
    f = open(path, "r")
    res = f.readlines()
    f.close()
    return res


def read_csv(path, delim):
    """Function to read and return a csv file in a list, choosing the delimiter."""
    
    f = open(path, "r")
    res = []
    for row in csv.reader(f, delimiter = delim):
        res.append(row)
    f.close()
    return res


def read_json(path):
    """Function to read a JSON file."""
    
    f = open(path, "r")
    res = f.read()
    data = json.loads(res)
    f.close()
    return data


def read_config(ini):
    """Runs the config file containing all the information to make a new model.
    
    ARGS :
        ini (str) -- the path to the .ini file.
    RETURN :
        config (dict of str) -- the configuration in a python dictionary object.
    """
    
    config = configparser.ConfigParser()
    config.read(ini)
    return config


def write_file(WD, filename, data):
    """Function to write a file from a list."""
    
    f = open(WD + filename, "w")
    for i in data:
        f.write(i)
    f.close()


def write_csv(WD, list_value, name, separator = ","):
    """Function to save a file as a CSV format, needs a list of lists, 
    first list as the column names."""
    
    with open(WD + name + '.csv', 'w', newline = '') as file:
        writer = csv.writer(file, delimiter = separator)
        for f in list_value:
            writer.writerow(f)


def get_reactions_PT(path):
    """Function to get the reactions in a reactions.dat file of Pathway Tools PGDB.
    
    ARGS:
        path (str) -- the path to the reactions.dat file.
    RETURN:
        liste_reac (list of str) -- the list containing all the reactions in this model.
    """
    
    liste_Reac = []
    PT_reac = open(path, "r")
    for line in PT_reac:
        if "UNIQUE-ID" in line and not "#" in line:
            try:
                liste_Reac.append(re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', line).group(0).rstrip())
            except AttributeError:
                print("No match for : ", line)
    return liste_Reac


def clean_sbml(WD, name):
    """Function to get rid of specific character COBRA puts into its 
    sbml models, and cannot read (you read it right...).
    
    ARGS:
        WD (str) -- the working directory.
        name (str) -- the name of the sbml model.
    RETURN:
        the name of the new file.
    """
    
    file = read_file(WD + name)
    new_file = []
    meta = re.compile('( id="M_)')
    sp_ref = re.compile('(<speciesReference species="M_)')
    reac = re.compile('(id="R_)')
    for i in file:
        i = meta.sub(' id="', i)
        i = sp_ref.sub('<speciesReference species="', i)
        if "<reaction" in i:
            i = reac.sub('id="', i)
        new_file.append(i)
    write_file(WD, "clean_" + name, new_file)
    return "clean_" + name


def cobra_compatibility(reac, side = True):
    """Function to transform a reaction ID into a cobra readable ID and vice versa.
    
    ARGS:
        reac (str) -- the reaction.
        side (boolean) -- True if you want to convert a COBRA ID into a readable ID, 
        False for the reverse.
    RETURN:
        reac (str) -- the transformed reaction.
    """
    
    if side:
        reac = reac.replace("__46__", ".").replace("__47__", "/").replace("__45__", "-").replace("__43__", "+").replace("__91__", "[").replace("__93__", "]")
        if re.search('(^_\d)', reac):
            reac = reac[1:]
    else:
        reac = reac.replace("/", "__47__").replace(".", "__46__").replace("-", "__45__").replace("+", "__43__").replace("[", "__91__").replace("]", "__93")
        if re.search('(\d)', reac[0]):
            reac = "_" + reac
    return reac


def metacyc_IDs(WD, path):
    """Function to make the correspondence file between short and long ID of Metacyc.
    
    ARGS:
        WD (str) -- the directory to save the correspondence file.
        path (str) -- the path to the metacyc model in JSON format.
    """
    
    data = read_json(path)
    res = []
    print(len(data["reactions"]))
    for reac in data["reactions"]:
        long_ID = reac["name"].split("/")[0]
        ###Getting rid of the brackets in the name sometimes!
        reac_pattern = re.compile('([[].*[]])')
        tmp_ID = reac_pattern.sub("", long_ID)
        short_ID = tmp_ID
        ###Regexp to "clean" the metabolite's names
        meta_pattern = re.compile('(_CC[OI]-.*)|(^[_])|([_]\D$)')
        meta_list = []
        # print(len(reac["metabolites"].keys()))
        if len(reac["metabolites"].keys()) != 0:
            for metabolite in reac["metabolites"].keys():
                ###The metabolite are "cleaned" here
                metabolite = meta_pattern.sub("", metabolite)
                len_ID, len_meta = len(tmp_ID), len(metabolite)
                diff = len_ID - len_meta
                ###Small trick to get only the end of the ID removed and not the beginning
                ###(metabolite's names can be in the reaction's name)
                test_ID = tmp_ID[:diff - 1] + tmp_ID[diff - 1:].replace("-" + metabolite, "")
                if len(test_ID) < len(short_ID):
                    short_ID = test_ID
            res.append([short_ID, reac["name"]])
    write_csv(WD, res, "MetacycCorresIDs", "\t")


def corres_dico(path, sep = "\t"):
    """Function to create a ditionary of correspondence between 
    short and long IDs from a correspondence file (Metacyc IDs).
    
    ARGS:
        path (str) -- the path to the file containing the correspondence information.
        sep (str) -- the separator of the correspondence file (default = tab).
    RETURN:
        dico_matching -- dictionary with short IDs as key and list of long IDs 
        as values (NB : one short ID can have several long IDs correspondence).
        dico_matching_rev -- dictionary with long IDs as key and the
        corresponding short ID as value (NB : one long ID as only one short ID correspondence).
    """
    
    matching = read_file(path)
    dico_matching = {}
    dico_matching_rev = {}
    for line in matching:
        if line:
            couple = line.rstrip().split(sep)
            if couple[0] in dico_matching.keys():
                dico_matching[couple[0]].append(couple[1])
            else:
                dico_matching[couple[0]] = [couple[1]]
            dico_matching_rev[couple[1]] = couple[0]
    return dico_matching, dico_matching_rev


def trans_short_ID(list_IDs, corres, short = True):
    """Function to transform short IDs that can be ambiguous
    into long ones thanks to the correspondence ID file.
    
    ARGS:
        list_IDs (list of str) --  the list of IDs to convert (must be Metacyc format IDs).
        corres (str) -- the path to the correspondence file of Metacyc IDs.
        short (boolean) -- True if you want to have short IDs becoming long,
        False if you want long IDs to become short (not implemented yet).
    RETURN:
        new_list (list of str) -- the list with the converted IDs.
    """
    
    dico_matching, dico_matching_rev = corres_dico(corres)
    new_list = []
    if short:
        for reac in list_IDs:
            reac = reac.rstrip()
            try:
                for long_reac in dico_matching[reac]:
                    new_list.append(long_reac)
            except KeyError:
                try:
                    dico_matching_rev[reac]
                    new_list.append(reac)
                except KeyError:
                    print("No match for reac : ", reac)
    else:
        for reac in list_IDs:
            reac = reac.rstrip()
            try:
                dico_matching[reac]
                new_list.append(reac)
            except KeyError:
                try:
                    new_list.append(dico_matching_rev[reac])
                except KeyError:
                    print("No match for reac : ", reac)
    return new_list


def transcript_to_gene(WD, model, transcript_corres, name):
    """Function to transform the transcripts in gene_reaction_rule into its corresponding genes.
    It creates a new model that will be rid of all the transcripts.
    
    ARGS:
        WD (str) -- the path where to find the model.
        model (str) -- the model file in json format.
        transcript_corres (str) -- the exact path of the csv file of the
        correspondance between a transcript and its gene.
        name (str) -- the new name for the new corrected model.
    """
    
    dico_corres, dico_corres_rev = corres_dico(transcript_corres) 
    transcript_model = cobra.io.load_json_model(WD + model)
    gene_model = cobra.Model(transcript_model.id)
    for reaction in transcript_model.reactions:
        genes = []
        list_transcripts = reaction.gene_reaction_rule.split(" or ")
        for transcript in list(filter(None, [trans.upper() for trans in list_transcripts])):
            try:
                genes.append(dico_corres_rev[transcript])
            except KeyError:
                try:
                    dico_corres[transcript]
                    genes.append(transcript)
                except KeyError:
                    print("No match for : ", transcript)
        new_reaction = copy.deepcopy(reaction)
        new_reaction.gene_reaction_rule = " or ".join(set(genes))
        gene_model.add_reactions([new_reaction])
    cobra.io.save_json_model(gene_model, WD + name + transcript_model.id + ".json")