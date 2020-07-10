# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file contains utility functions used as is in several scripts."""

import configparser
import csv
import json
import re


def read_file(path):
    """Function to read and return a file line by line in a list."""
    
    f = open(path, "r")
    res = f.readlines()
    f.close()
    return res


def write_file(WD, filename, data):
    """Function to write a file from a list."""
    
    f = open(WD + filename, "w")
    for i in data:
        f.write(i)
    f.close()


def write_csv(WD, list_value, name):
    """Function to save a file as a CSV format, needs a list of lists, 
    first list as the column names."""
    
    with open(WD + name + '.csv', 'w', newline = '') as file:
        writer = csv.writer(file)
        for f in list_value:
            writer.writerow(f)


def write_tsv(WD, list_value, name):
    """Function to save a file as a TSV format, needs a list of lists, 
    first list as the column names."""
    
    with open(WD + name + '.tsv', 'w', newline = '') as file:
        writer = csv.writer(file, delimiter = "\t")
        for f in list_value:
            writer.writerow(f)


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
        if "UNIQUE-ID" in line:
            try:
                liste_Reac.append(re.search('(?<=UNIQUE-ID - )\w+(.*\w+)*(-*\w+)*', line).group(0).rstrip())
            except AttributeError:
                pass
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
        side (boolean) -- True if you want to convert an ID into a COBRA readable ID, 
        False for the reverse.
    RETURN:
        reac (str) -- the transformed reaction.
    """
    
    if side:
        reac = reac.replace("__47__", "/").replace("__46__", ".").replace("__45__", "-")
        if re.search('(_\d)', reac):
            reac = reac[1:]
    else:
        reac = reac.replace("/", "__47__").replace(".", "__46__").replace("-", "__45__")
        if re.search('(\d)', reac[0]):
            reac = "_" + reac
    return reac