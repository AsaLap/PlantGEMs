# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file contains utility functions used in several PlantGEMs' scripts."""

import sys
import configparser
import csv
import json
import os
import pickle
import subprocess
import re


def build_correspondence_dict(path, sep="\t"):
    """Function to create a dictionary of correspondence between
    elements from a correspondence file (Metacyc short and long IDs for example).

    PARAMS:
        path (str) -- the path to the csv/tsv file containing the corresponding information (only two elements).
        sep (str) -- the separator of the correspondence file (default = tab).
    RETURNS:
        matching_dict -- dictionary with an element as key that matches the value.
        matching_dict_reversed -- same thing as matching_dict but in reversed (value is key and vice-versa).
    """

    matching = read_file_listed(path)
    matching_dict = {}
    matching_dict_reversed = {}
    for line in matching:
        if line:
            couple = line.rstrip().split(sep)
            if couple[0] in matching_dict.keys():
                matching_dict[couple[0]].append(couple[1])
            else:
                matching_dict[couple[0]] = [couple[1]]
            matching_dict_reversed[couple[1]] = couple[0]
    return matching_dict, matching_dict_reversed


def check_path(path, sys_exit=False):
    """Function to check a file or folder existence.

    ARGS:
        path (str) -- the path to check.
        sys_exit (boolean) -- True if you want to stop the process and quit, False (default) if you juste want the
            return to be False if applicable (doesn't stop the process).
        """

    if os.path.exists(path):
        return True
    elif sys_exit:
        sys.exit("File or folder not found : " + path)
    else:
        return False


def clean_sbml(directory, name):
    """Function to get rid of specific character COBRA puts into its 
    sbml models, and cannot read (you read it right...).
    
    PARAMS:
        wd (str) -- the working directory.
        name (str) -- the name of the sbml model.
    RETURNS:
        the name of the new file.
    """

    file = read_file_listed(directory + name)
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
    write_file(directory + "clean_" + name, new_file)
    return "clean_" + name


def cobra_compatibility(reaction, side=True):
    """Function to transform a reaction ID into a cobra readable ID and vice versa.
    
    PARAMS:
        reaction (str) -- the reaction.
        side (boolean) -- True if you want to convert a COBRA ID into a readable ID, 
        False for the reverse.
    RETURNS:
        reaction (str) -- the transformed reaction.
    """

    if side:
        reaction = reaction.replace("__46__", ".").replace("__47__", "/").replace("__45__", "-") \
            .replace("__43__", "+").replace("__91__", "[").replace("__93__", "]")
        if re.search('(^_\d)', reaction):
            reaction = reaction[1:]
    else:
        reaction = reaction.replace("/", "__47__").replace(".", "__46__").replace("-", "__45__") \
            .replace("+", "__43__").replace("[", "__91__").replace("]", "__93")
        if re.search('(\d)', reaction[0]):
            reaction = "_" + reaction
    return reaction


def copy_file(start, end):
    if check_path(start):
        try:
            subprocess.run(["cp", start, end])
        except PermissionError:
            print("Permission to create this folder :\n" + end + "\nnot granted !")
        except FileNotFoundError:
            print("Path not found : ", end)
    else:
        print("File not found : " + start)


def find_file(directory, target, extension):  # TODO : take multiple file extensions in parameters
    """Search a file corresponding to the target in the 'files' directory. If no match is found,
    asks the user to input the exact path to the file he wants to use.

    PARAMS:
        directory (str) -- the directory where to find the file.
        target (str) -- the name of the file to search if correctly named.
        extension (str) -- the file's extension.

    RETURNS:
        file_path (str) -- the exact path to the file.
    """

    file_path = slash(directory) + target + dot(extension)
    if os.path.isfile(file_path):
        return file_path
    else:
        print("No corresponding file found here : " + directory)
        try:
            file_path = str(input("Path to " + target + dot(extension) + "'s file : "))
            if os.path.isfile(file_path):
                return file_path
            else:
                print("No file found with this path, make sure you entered it correctly... Restarting.")
                find_file(directory, target, extension)
        except ValueError:
            print("Please enter a string only... Restarting.")
            find_file(directory, target, extension)


def find_files(directory, extension):
    """Function to return the list of files in a given directory with a certain extension"""

    return [i for i in os.listdir(slash(directory)) if i.endswith(dot(extension))]


def get_list_directory(path):
    """Function to retrieve all the directories names at a specified location (path)

    ARGS:
        path (str) -- the folder path in which are searched the sub-folders.
    RETURNS:
        list_directory (list of str) -- the names of the found folders.
    """

    path = slash(path)
    list_directory = []
    for found in os.listdir(path):
        if os.path.isdir(path + found):
            list_directory.append(found)
    return list_directory


def get_metacyc_ids(metacyc_json_model_path):
    """Function to make the correspondence file between short and long ID of Metacyc."""

    data = read_json(metacyc_json_model_path)
    res = []
    print("Number of reactions found in the metacyc.json file : {}", format(len(data["reactions"])))
    for reaction in data["reactions"]:
        long_id = reaction["name"].split("/")[0]
        # Getting rid of the brackets in the name sometimes!
        reaction_pattern = re.compile('([[].*[]])')
        tmp_id = reaction_pattern.sub("", long_id)
        short_id = tmp_id
        # Regexp to "clean" the metabolite's names
        meta_pattern = re.compile('(_CC[OI]-.*)|(^[_])|([_]\D$)')
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
    write_csv(os.path.dirname(metacyc_json_model_path), "/metacyc_ids", res, "\t")


def get_sequence_region(gff_file_path):
    """Function that browses the .gff file and gets the name of each gene,
    the position in the genome and each corresponding region and protein(s).

    RETURNS:
        regions_dict -- a dictionary containing all the gathered information (see pipelinePT() for the structure).

    -- Structure of regions_dict:
    {Region name (str):
        {Gene name (str):
            {"Start": int, "End": int, "Proteins":
                {Protein name (str): [CDS's pos (tuple)]}}}}
    """

    regions_dict = {}
    gff_file = read_file_listed(gff_file_path)
    protein_found = False  # Boolean to avoid testing a protein on each line that has already been found.
    cds_break = False  # Boolean to avoid an error if the CDS's name hasn't been found.
    region = None  # Assignment before use
    gene = None  # Assignment before use
    protein = None  # Assignment before use
    rna = None  # Assignment before use
    for line in gff_file:
        if "RNA\t" in line:
            protein_found = False
            cds_break = False
            try:
                rna = re.search('(?<=ID=)[RrNnAaTtSsCcIiPpMm]*[:-]*\w+(\.\w+)*(\-\w+)*', line).group(0)
                garbage = re.search('([RrNnAaTtSsCcIiPpMm]*[:-])*', rna).group(0)
                rna = str.replace(rna, garbage, "")
            except AttributeError:
                print("RNA ID not found, might cause a problem if CDS name is different than RNA.\n", line)
                pass
        if "\tgene\t" in line:  # Searching the gene's information
            protein_found = False
            spl = line.split("\t")
            region = spl[0]
            try:
                gene = re.search('(?<=ID=)[GgEeNn]*[:-]*\w+(\.\w+)*(\-\w+)*', line).group(0)
                garbage = re.search('([GgEeNn]*[:-])*', gene).group(0)
                gene = str.replace(gene, garbage, "")
            except AttributeError:
                print("The gene name hasn't been found here : {}", format(line))
                break
            if region not in regions_dict.keys():
                regions_dict[region] = {}
            regions_dict[region][gene] = {"Start": spl[3], "End": spl[4], "Proteins": {}}
        if region and gene and not protein_found and "\tCDS\t" in line:  # Searching the protein's information
            try:
                protein = re.search('(?<=ID=)[CcDdSs]*[:-]*\w+(\.\w+)*', line).group(0)
                garbage = re.search('([CcDdSs]*[:-])*', protein).group(0)
                protein = str.replace(protein, garbage, "")
                if protein[:-2] == rna:  # I check if the GFF is compatible with the faa file for the tsv annotation
                    protein, rna = rna, protein
                regions_dict[region][gene]["Proteins"][protein] = []
                protein_found = True
                cds_break = False
            except AttributeError:
                print("The CDS has no attribute 'ID='...")
                cds_break = True
        if not cds_break and "\tCDS\t" in line:  # Searching the CDS' information
            spl = line.split("\t")
            regions_dict[region][gene]["Proteins"][protein].append([int(spl[3]), int(spl[4])])
    return regions_dict


def get_list_ids_reactions_cobra(model):
    """Function to gather all the reactions' id of a cobra model in a list.

    PARAMS:
        model -- a cobra model.
    RETURNS:
        res -- a list containing all the model's reactions' id.
    """

    res = []
    for reaction in model.reactions:
        res.append(reaction.id)
    return res


def get_list_reactions_cobra(model):
    """Function to gather all the reactions of a cobra model in a list.

    PARAMS:
        model -- a cobra model.
    RETURNS:
        res -- a list containing all the model's reactions.
    """

    res = []
    for reaction in model.reactions:
        res.append(reaction)
    return res


def load_obj(path):
    """Loads a pickle object."""

    if check_path(path):
        with open(path, 'rb') as input_file:
            return pickle.load(input_file)


def make_directory(directory):
    """Function to create a directory with a given path."""

    if not os.path.isdir(directory):
        new_dir = directory.strip(" /").split("/")[-1]
        base_dir = directory.strip(" /")[:-len(directory.strip(" /").split("/")[-1])]
        print("Creation of directory '" + new_dir + "' in '" + base_dir + "'")
        try:
            subprocess.run(["mkdir", directory])
        except PermissionError:
            print("Permission to create this folder :\n" + directory + "\nnot granted !")
        except FileNotFoundError:
            print("Path not found : ", directory)
    else:
        print("Directory already exists : ", directory)


def migrate(main_directory):
    blast_directory = main_directory + "blast/"
    merge_directory = main_directory + "merge/"
    mpwt_directory = main_directory + "mpwt/"
    make_directory(merge_directory)
    list_species = get_list_directory(blast_directory)
    for species in list_species:
        make_directory(merge_directory + species)
        copy_file(blast_directory + species + "/" + species + "_blast_draft.json",
                  merge_directory + species + "/" + species + "_blast_draft.json")
        copy_file(mpwt_directory + "/output/" + species + "/reactions.dat",
                  merge_directory + species + "/reactions.dat")
        copy_file(mpwt_directory + "/output/" + species + "/proteins.dat",
                  merge_directory + species + "/proteins.dat")
        copy_file(mpwt_directory + "/output/" + species + "/enzrxns.dat",
                  merge_directory + species + "/enzrxns.dat")


def dot(extension):
    """Utility function to put a dot before the extension if needed."""

    if extension[0] != ".":
        return "." + extension
    else:
        return extension


def read_config(ini):
    """Runs the config file containing all the information to make a new model.

    PARAMS :
        ini (str) -- the path to the .ini file.
    RETURNS :
        config (dict of str) -- the configuration in a python dictionary object.
    """

    if os.path.isfile(ini):
        config = configparser.ConfigParser()
        config.read(ini)
    else:
        print("File not found : " + ini + "\nEnding the process.")
        sys.exit()
    return config


def read_csv(path, delim):
    """Function to read and return a csv file in a list, choosing the delimiter."""

    f = open(path, "r")
    res = []
    for row in csv.reader(f, delimiter=delim):
        res.append(row)
    f.close()
    return res


def read_file_listed(path):
    """Function to read and return a file line by line in a list."""

    f = open(path, "r")
    res = f.readlines()
    f.close()
    return res


def read_file_stringed(path):
    """Function to read and return a file in a string."""

    f = open(path, "r")
    res = f.read()
    f.close()
    return res


def read_json(path):
    """Function to read a JSON file."""

    f = open(path, "r")
    res = f.read()
    data = json.loads(res)
    f.close()
    return data


def remove_directory(directory):
    try:
        subprocess.run(["rm", "-rf", directory])
    except FileNotFoundError:
        pass
    except PermissionError:
        print("Permission to erase this folder :\n" + directory + "\nnot granted !")


def save_obj(obj, path):
    """Saves an object in a pickle file."""

    with open(path + '.pkl', 'wb+') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def slash(directory):
    """Utility function to put a slash after a directory path if needed."""

    if directory[-1] != "/":
        return directory + "/"
    else:
        return directory


def trans_short_id(list_ids, correspondence, short=True, keep=False):
    """Function to transform short IDs that can be ambiguous
    into long ones thanks to the correspondence ID file.

    PARAMS:
        list_ids (list of str) --  the list of IDs to convert (must be Metacyc format IDs).
        correspondence (str) -- the path to the correspondence file of Metacyc IDs.
        short (boolean) -- True if you want to have short IDs becoming long,
        False if you want long IDs to become short.
        keep (boolean) -- True if you want to keep the reactions even if they are not found.
    RETURNS:
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


def write_csv(directory, name, list_value, separator=","):
    """Function to save a file as a CSV format, needs a list of lists,
    first list as the column names."""

    if separator == "\t":
        extension = ".tsv"
    else:
        extension = ".csv"
    with open(slash(directory) + name + extension, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=separator)
        for f in list_value:
            writer.writerow(f)


def write_file(path, data, strip=True):
    """Function to write a file from a list."""

    f = open(path, "w")
    if strip:
        for i in data:
            f.write(i.rstrip() + "\n")
    else:
        for i in data:
            f.write(i)
    f.close()


if __name__ == "__main__":
    sys.exit("Please launch the pipeline with PlantGEMs.py. Terminating the process.")
