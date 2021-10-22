# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file contains utility functions used in several PlantGEMs' scripts."""

import cobra
import configparser
import copy
import csv
import json
import pickle
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
    for row in csv.reader(f, delimiter=delim):
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


def write_file(wd, filename, data):
    """Function to write a file from a list."""

    f = open(wd + filename, "w")
    for i in data:
        f.write(i.rstrip() + "\n")
    f.close()


def write_csv(wd, list_value, name, separator=","):
    """Function to save a file as a CSV format, needs a list of lists, 
    first list as the column names."""

    with open(wd + name + '.csv', 'w', newline='') as file:
        writer = csv.writer(file, delimiter=separator)
        for f in list_value:
            writer.writerow(f)


def save_obj(obj, path):
    """Saves an object in a pickle file."""

    with open(path + '.pkl', 'wb+') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def load_obj(path):
    """Loads a pickle object."""

    with open(path, 'rb') as input_file:
        return pickle.load(input_file)


def get_pwt_reactions(path):
    """Function to get the reactions in a reactions.dat file of Pathway Tools PGDB.
    
    ARGS:
        path (str) -- the path to the reactions.dat file.
    RETURN:
        liste_reactions (list of str) -- the list containing all the reactions in this model.
    """

    list_reactions = []
    pwt_reactions = open(path, "r")
    for line in pwt_reactions:
        if "UNIQUE-ID" in line and "#" not in line:
            try:
                list_reactions.append(re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', line).group(0).rstrip())
            except AttributeError:
                print("No match for : ", line)
    return list_reactions


def clean_sbml(wd, name):
    """Function to get rid of specific character COBRA puts into its 
    sbml models, and cannot read (you read it right...).
    
    ARGS:
        wd (str) -- the working directory.
        name (str) -- the name of the sbml model.
    RETURN:
        the name of the new file.
    """

    file = read_file(wd + name)
    new_file = []
    meta = re.compile('( id="M_)')
    sp_ref = re.compile('(<speciesReference species="M_)')
    reaction = re.compile('(id="R_)')
    for i in file:
        i = meta.sub(' id="', i)
        i = sp_ref.sub('<speciesReference species="', i)
        if "<reaction" in i:
            i = reaction.sub('id="', i)
        new_file.append(i)
    write_file(wd, "clean_" + name, new_file)
    return "clean_" + name


def cobra_compatibility(reaction, side=True):
    """Function to transform a reaction ID into a cobra readable ID and vice versa.
    
    ARGS:
        reaction (str) -- the reaction.
        side (boolean) -- True if you want to convert a COBRA ID into a readable ID, 
        False for the reverse.
    RETURN:
        reaction (str) -- the transformed reaction.
    """

    if side:
        reaction = reaction.replace("__46__", ".").replace("__47__", "/").replace("__45__", "-") \
            .replace("__43__", "+").replace("__91__", "[").replace("__93__", "]")
        if re.search('(^_\d)', reaction):
            reaction = reaction[1:]
    else:
        reaction = reaction.replace("/", "__47__").replace(".", "__46__").replace("-", "__45__")\
            .replace("+", "__43__").replace("[", "__91__").replace("]", "__93")
        if re.search('(\d)', reaction[0]):
            reaction = "_" + reaction
    return reaction


def metacyc_ids(wd, path):
    """Function to make the correspondence file between short and long ID of Metacyc.
    
    ARGS:
        wd (str) -- the directory to save the correspondence file.
        path (str) -- the path to the metacyc model in JSON format.
    """

    data = read_json(path)
    res = []
    print(len(data["reactions"]))
    for reaction in data["reactions"]:
        long_id = reaction["name"].split("/")[0]
        # Getting rid of the brackets in the name sometimes!
        reaction_pattern = re.compile('([[].*[]])')
        tmp_id = reaction_pattern.sub("", long_id)
        short_id = tmp_id
        # Regexp to "clean" the metabolite's names
        meta_pattern = re.compile('(_CC[OI]-.*)|(^[_])|([_]\D$)')
        meta_list = []
        if len(reaction["metabolites"].keys()) != 0:
            for metabolite in reaction["metabolites"].keys():
                # The metabolite are "cleaned" here
                metabolite = meta_pattern.sub("", metabolite)
                len_id, len_meta = len(tmp_id), len(metabolite)
                diff = len_id - len_meta
                # Small trick to get only the end of the ID removed and not the beginning
                # (metabolite's names can be in the reaction's name)
                test_id = tmp_id[:diff - 1] + tmp_id[diff - 1:].replace("-" + metabolite, "")
                if len(test_id) < len(short_id):
                    short_id = test_id
            res.append([short_id, reaction["name"]])
    write_csv(wd, res, "MetacycCorresIDs", "\t")


def build_correspondence_dict(path, sep="\t"):
    """Function to create a ditionary of correspondence between 
    short and long IDs from a correspondence file (Metacyc IDs).
    
    ARGS:
        path (str) -- the path to the file containing the correspondence information.
        sep (str) -- the separator of the correspondence file (default = tab).
    RETURN:
        metacyc_matching_id_dict -- dictionary with short IDs as key and list of long IDs
        as values (NB : one short ID can have several long IDs correspondence).
        metacyc_reverse_id_dict -- dictionary with long IDs as key and the
        corresponding short ID as value (NB : one long ID as only one short ID correspondence).
    """

    matching = read_file(path)
    metacyc_matching_id_dict = {}
    metacyc_matching_id_dict_reversed = {}
    for line in matching:
        if line:
            couple = line.rstrip().split(sep)
            if couple[0] in metacyc_matching_id_dict.keys():
                metacyc_matching_id_dict[couple[0]].append(couple[1])
            else:
                metacyc_matching_id_dict[couple[0]] = [couple[1]]
            metacyc_matching_id_dict_reversed[couple[1]] = couple[0]
    return metacyc_matching_id_dict, metacyc_matching_id_dict_reversed


def trans_short_id(list_ids, correspondence, short=True, keep=False):
    """Function to transform short IDs that can be ambiguous
    into long ones thanks to the correspondence ID file.
    
    ARGS:
        list_ids (list of str) --  the list of IDs to convert (must be Metacyc format IDs).
        correspondence (str) -- the path to the correspondence file of Metacyc IDs.
        short (boolean) -- True if you want to have short IDs becoming long,
        False if you want long IDs to become short.
        keep (boolean) -- True if you want to keep the reactions even if they are not found.
    RETURN:
        new_list (list of str) -- the list with the converted IDs.
    """

    metacyc_matching_id_dict, metacyc_matching_id_dict_reversed = build_correspondence_dict(correspondence)
    new_list = []
    if short:
        for reaction in list_ids:
            reaction = reaction.rstrip()
            try:
                for long_reaction in metacyc_matching_id_dict[reaction]:
                    new_list.append(long_reaction)
            except KeyError:
                try:
                    if reaction in metacyc_matching_id_dict_reversed.keys():
                        new_list.append(reaction)
                except KeyError:
                    print("No match for reaction : ", reaction)
                    if keep:
                        new_list.append(reaction)
                        print("Keeping it anyway...")
    else:
        for reaction in list_ids:
            reaction = reaction.rstrip()
            try:
                if reaction in metacyc_matching_id_dict.keys():
                    new_list.append(reaction)
            except KeyError:
                try:
                    new_list.append(metacyc_matching_id_dict_reversed[reaction])
                except KeyError:
                    print("No match for reaction : ", reaction)
                    if keep:
                        new_list.append(reaction)
                        print("Keeping it anyway...")
    return new_list


def protein_to_gene(wd, model, protein_correspondence, name):
    """Function to transform the proteins in gene_reaction_rule into its corresponding genes.
    It creates a new model that will be rid of all the proteins.
    
    ARGS:
        wd (str) -- the path where to find the model.
        model (str) -- the model file in json format.
        protein_correspondence (str) -- the exact path of the csv file of the
        correspondence between a protein and its gene.
        name (str) -- the new name for the new corrected model.
    """

    metacyc_matching_id_dict, metacyc_matching_id_dict_reversed = build_correspondence_dict(protein_correspondence)
    protein_model = cobra.io.load_json_model(wd + model)
    gene_model = cobra.Model(protein_model.id)
    for reaction in protein_model.reactions:
        genes = []
        list_proteins = reaction.gene_reaction_rule.split(" or ")
        for protein in list(filter(None, [trans.upper() for trans in list_proteins])):
            try:
                genes.append(metacyc_matching_id_dict_reversed[protein])
            except KeyError:
                try:
                    if protein in metacyc_matching_id_dict.keys():
                        genes.append(protein)
                except KeyError:
                    print("No match for : ", protein)
        new_reaction = copy.deepcopy(reaction)
        new_reaction.gene_reaction_rule = " or ".join(set(genes))
        gene_model.add_reactions([new_reaction])
    cobra.io.save_json_model(gene_model, wd + name + protein_model.id + ".json")
