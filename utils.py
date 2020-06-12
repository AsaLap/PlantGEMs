# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file contains utility functions"""

import configparser
import json
import re

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


def read_json(path):
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