# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used to merge two different GEMs, one coming from my own model based
reconstruction and another one from a Pathway Tools reconstruction with the AuReMe's 
mpwt package."""

import argparse
import cobra
import copy
import logging
import multiprocessing
import os
import re

from datetime import date

import graphs
import module
import utils


class Merging(module.Module):

    def __init__(self, _name, _main_directory):
        """Explanations here"""
        super().__init__(_name, _main_directory)
        self.directory = self.main_directory + "merge/" + self.name + "/"
        self.files_directory = self.main_directory + "files/"
        self.pwt_reactions_id_list = []
        self.pwt_metacyc_reactions_id_list = []
        self.json_reactions_list = []
        self.sbml_reactions_list = []
        self.dict_upsetplot_reactions = {}
        self.merged_model = cobra.Model(self.name, name=self.name + "_PlantGEMs_" + str(date.today()))

        # Metacyc files
        self.metacyc_file_path = utils.find_file(self.files_directory, "metacyc", "json")
        self.metacyc_model = cobra.io.load_json_model(self.metacyc_file_path)
        if not os.path.isfile(self.files_directory + "metacyc_ids.tsv"):
            utils.get_metacyc_ids(self.metacyc_file_path)
        self.metacyc_ids_file_path = utils.find_file(self.files_directory, "metacyc_ids", "tsv")
        self.metacyc_matching_id_dict, self.metacyc_matching_id_dict_reversed = \
            utils.build_correspondence_dict(self.metacyc_ids_file_path)

    def _search_metacyc_reactions_ids(self):
        """Function to search the reactions' ids for in the flat files from a Pathway Tools reconstruction. Then, adds
        them into the object's Pathway Tools' reactions list."""

        if self.pwt_reactions_id_list:
            pwt_metacyc_long_id_list = []
            pwt_metacyc_no_match_id_list = []
            self.dict_upsetplot_reactions["Pathway_Tools"] = []
            for reaction in self.pwt_reactions_id_list:
                try:
                    var = self.metacyc_matching_id_dict[reaction]
                    self.pwt_metacyc_reactions_id_list.append(reaction)
                    self.dict_upsetplot_reactions["Pathway_Tools"].append(reaction)
                except KeyError:
                    try:
                        self.pwt_metacyc_reactions_id_list.append(self.metacyc_matching_id_dict_reversed[
                                                                      reaction])
                        pwt_metacyc_long_id_list.append(reaction)
                        self.dict_upsetplot_reactions["Pathway_Tools"].append(reaction)
                    except KeyError:
                        pwt_metacyc_no_match_id_list.append(reaction)
            logging.info(
                "{} : List of despecialized reactions ({}): \n{}".format(self.name, len(pwt_metacyc_long_id_list), "\n".
                                                                         join([i for i in pwt_metacyc_long_id_list])))
            logging.info(
                "{} : List of unmatched reactions ({}): \n{}".format(self.name, len(pwt_metacyc_no_match_id_list), "\n".
                                                                     join([i for i in pwt_metacyc_no_match_id_list])))

    def _get_networks_reactions(self, extension):
        """
        Search all the reactions in a model and add them to the object's list of regarding the extension of the file.

        PARAMS:
            extension (str) -- Extension of the model. json/JSON or sbml/SBML only for the moment.
        """

        list_networks = utils.find_files(self.directory, extension)
        if list_networks:
            count = 0
            if extension in ["json", "JSON"]:
                count += 1
                for json_model_file in list_networks:
                    logging.info("{} : JSON model network found : {}".format(self.name, json_model_file))
                    json_model = cobra.io.load_json_model(self.directory + json_model_file)
                    self.json_reactions_list.extend(utils.get_list_reactions_cobra(json_model))
                    self.dict_upsetplot_reactions[json_model.id + "_j" + str(count)] = utils. \
                        get_list_ids_reactions_cobra(json_model)
            count = 0
            if extension in ["sbml", "SBML"]:
                count += 1
                for sbml_model_file in list_networks:
                    logging.info("{} : SBML model network found : {}".format(self.name, sbml_model_file))
                    sbml_model = cobra.io.read_sbml_model(self.directory + sbml_model_file)
                    self.sbml_reactions_list.extend(utils.get_list_reactions_cobra(sbml_model))
                    self.dict_upsetplot_reactions[sbml_model.id + "_s" + str(count)] = utils. \
                        get_list_ids_reactions_cobra(sbml_model)
        else:
            logging.info("{} : No {} file of draft network or model found".format(self.name, extension))

    def _get_pwt_reactions(self):
        """Function to get the reactions in a reactions.dat file of Pathway Tools PGDB.

        PARAMS:
            path (str) -- the path to the reactions.dat file.
        RETURNS:
            list_reactions (list of str) -- the list containing all the reactions in this model.
        """

        pwt_reactions = open(self.directory + "/reactions.dat", "r")
        count = 0
        for line in pwt_reactions:
            if "UNIQUE-ID" in line and "#" not in line:
                try:
                    self.pwt_reactions_id_list.append(
                        re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', line).group(0).rstrip())
                    count += 1
                except AttributeError:
                    logging.error("{} : No match for : {}".format(self.name, line))
        logging.info("{} : Number of reactions found in the Pathway Tools reconstruction files : {}".format(self.name,
                                                                                                            count))

    def _browse_pwt_dat_files(self, reaction):
        """Function to correct the gene reaction rule in each reaction taken from Metacyc/Pathway Tools
        to make it fit the organism for which the model is reconstructed.

        PARAMS:
            reaction: the reaction's gene reaction rule to change.
        Returns the corrected gene reaction rule.
        """

        reactions_file = utils.read_file_listed(self.directory + "/reactions.dat")
        enzrxns_file = utils.read_file_listed(self.directory + "/enzrxns.dat")
        proteins_file = utils.read_file_listed(self.directory + "/proteins.dat")

        # First step : gathering the ENZYME-REACTION fields in enzrxns (could be several or none for one ID).
        stop = False
        enzrxns = []
        unique_id = ""
        for reaction_line in reactions_file:
            if "UNIQUE-ID" in reaction_line and "#" not in reaction_line:
                if stop:
                    break
                try:
                    unique_id = re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', reaction_line).group(0).rstrip()
                    if unique_id == reaction.name or unique_id == self.metacyc_matching_id_dict_reversed[reaction.name]:
                        stop = True
                except AttributeError:
                    logging.error("{} : No UNIQUE-ID match for reactions.dat : {}".format(self.name, reaction_line))
            if stop and "ENZYMATIC-REACTION " in reaction_line and "#" not in reaction_line:
                try:
                    enzrxns.append(re.search('(?<=ENZYMATIC-REACTION - )[+-]*\w+(.*\w+)*(-*\w+)*',
                                             reaction_line).group(0).rstrip())
                except AttributeError:
                    logging.error("{} : No ENZYMATIC-REACTION match for reactions.dat : {}".format(self.name,
                                                                                                   reaction_line))

        # Second step : getting the corresponding ENZYME for each ENZYME-REACTION.
        stop = False
        gene_list = []
        enzyme = ""
        no_match_enzrxns = []
        if enzrxns:
            for enzrxn in enzrxns:
                for line_enzrxn in enzrxns_file:
                    if "UNIQUE-ID" in line_enzrxn and "#" not in line_enzrxn:
                        if stop:
                            stop = False
                            break
                        try:
                            unique_id_rxn = re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', line_enzrxn).group(
                                0).rstrip()
                            if unique_id_rxn == enzrxn:
                                stop = True
                        except AttributeError:
                            logging.error("{} : No UNIQUE-ID match for enzrxns.dat : {}".format(self.name, line_enzrxn))
                    if stop and "ENZYME " in line_enzrxn and "#" not in line_enzrxn:
                        try:
                            enzyme = re.search('(?<=ENZYME - )[+-]*\w+(.*\w+)*(-*\w+)*', line_enzrxn).group(0).rstrip()
                        except AttributeError:
                            logging.error("{} : No ENZYME match for enzrxns.dat : {}".format(self.name, line_enzrxn))
                # Third step into the second one : getting the corresponding GENE for each ENZYME and put it into
                # geneList (which contains all that we're looking for).
                for lineProt in proteins_file:
                    if "UNIQUE-ID " in lineProt and "#" not in lineProt:
                        if stop:
                            stop = False
                            break
                        try:
                            unique_id_prot = re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', lineProt).group(
                                0).rstrip()
                            if unique_id_prot == enzyme:
                                stop = True
                        except AttributeError:
                            logging.error("No UNIQUE-ID match for proteins.dat : {}".format(lineProt))
                    if stop and "GENE " in lineProt and "#" not in lineProt:
                        try:
                            gene_list.append(
                                re.search('(?<=GENE - )[+-]*\w+(.*\w+)*(-*\w+)*', lineProt).group(0).rstrip())
                        except AttributeError:
                            logging.error("No GENE match for proteins.dat : {}".format(lineProt))
        else:
            no_match_enzrxns.append(unique_id)
            pass
        reaction.gene_reaction_rule = " or ".join(set(gene_list))
        return reaction, [reaction.name, len(enzrxns)], no_match_enzrxns

    def _correct_pwt_reactions(self, verbose=True):
        """
        Function to correct the gene_reaction_rule of each reaction coming from the Pathway Tools' software as they are
        copied from Metacyc and therefore wrongfully linked to the proper genes of the species reconstructed.

        PARAMS:
            verbose (bool) -- toggles the printing of work progression.
        """

        count = 0
        list_no_match_correction = []
        list_no_match_enzrxns = []
        list_match_nb_enzymatic_reactions = []
        for reaction_id in self.pwt_reactions_id_list:
            if count % 100 == 0:
                if verbose:
                    print("%s : gene correction %s ou of %s" % (self.name, str(count),
                                                                str(len(self.pwt_reactions_id_list))))
            count += 1
            try:  # Changing pwt ids into long ids to match Metacyc ones
                long_reactions = self.metacyc_matching_id_dict[reaction_id]
                for reaction_id2 in long_reactions:
                    if reaction_id2[0] in ['0', '1', '2', '3', '4', '5', '6', '7', '8',
                                           '9']:  # Another Metacyc ids' specificity
                        reaction_id2 = "_" + reaction_id2
                    try:
                        added_reaction = copy.deepcopy(self.metacyc_model.reactions.get_by_id(reaction_id2))
                        added_reaction_corrected, tuple_nb_enzymatic_reactions_match, no_match_enzrxns = \
                            self._browse_pwt_dat_files(added_reaction)
                        list_match_nb_enzymatic_reactions.append(tuple_nb_enzymatic_reactions_match)
                        list_no_match_enzrxns.extend(no_match_enzrxns)
                        self.merged_model.add_reactions([added_reaction_corrected])
                    except KeyError:
                        list_no_match_correction.append(reaction_id2)
                        pass
            except KeyError:
                list_no_match_correction.append(reaction_id)
                pass
        logging.info("{} : No match in gene correction for reactions ({}) : \n{}".format(self.name, len(
            list_no_match_correction), "\n".join([i for i in list_no_match_correction])))
        logging.info("{} : No enzrxns entry for unique-ids ({}) :\n{}".format(self.name, len(list_no_match_enzrxns),
                                                                              "\n".
                                                                              join([i for i in list_no_match_enzrxns])))
        logging.info("{} : Number of enzymatic reaction(s) found associated to each reaction : \n{}".
                     format(self.name, "\n".
                            join([(str(i[0]) + " : " + str(i[1])) for i in list_match_nb_enzymatic_reactions])))

    def _conservative_merging(self, merging_reactions_list):
        """
        Function to merge reactions in the merging_reactions_list that are already in the merged model taking care of
        keeping every gene of the gene reaction rule, either from the merged model or the reactions' list. Just adds the
        reactions that were not already present.

        PARAMS:
            merging_reactions_list (cobra list reactions) -- list of cobra reactions you want to merge to the object's
            model.
        RETURNS:
            a list of reactions that are not already in the merged_model
        """

        temp_model = cobra.Model("temp_" + self.name)
        for reaction in merging_reactions_list:
            temp_model.add_reactions([reaction])
        merging_reactions_list_ids = utils.get_list_ids_reactions_cobra(temp_model)
        for reaction in self.merged_model.reactions:
            if reaction.id in merging_reactions_list_ids:
                list_genes = reaction.gene_reaction_rule.split(" or ")
                list_genes.extend(temp_model.reactions.get_by_id(reaction.id).gene_reaction_rule.split(" or "))
                list_genes = list(filter(None, set(list_genes)))
                reaction.gene_reaction_rule = " or ".join(list_genes)
                temp_model.remove_reactions([reaction.id])
        self.merged_model.add_reactions(utils.get_list_reactions_cobra(temp_model))

    def _merge(self):
        """Function to merge models, either from a model-based reconstruction (blasting module) or from
        Pathway Tools's Pathologic software."""

        self._correct_pwt_reactions()
        logging.info("{} : network size with only Pathway Tools' reactions : {}", format(self.name, str(len(
            self.merged_model.reactions))))
        self._conservative_merging(self.json_reactions_list)
        logging.info("{} : network size with the addition of JSON's reactions : {}", format(self.name, str(len(
            self.merged_model.reactions))))
        self._conservative_merging(self.sbml_reactions_list)
        logging.info("{} : network size with the addition of SBML's reactions : {}", format(self.name, str(len(
            self.merged_model.reactions))))

    def build(self):
        """Function to call the method in correct order for a complete merging."""

        if utils.check_path(self.directory + "reactions.dat"):
            self._get_pwt_reactions()
            self._search_metacyc_reactions_ids()
        else:
            logging.info("{} : No .dat files found, proceeding with drafts and models only.".format(self.name))
        if utils.check_path(self.directory, sys_exit=True):
            self._get_networks_reactions("json")
            self._get_networks_reactions("sbml")
        self._merge()
        cobra.io.save_json_model(self.merged_model, self.directory + self.name + "_merged.json")
        graphs.make_upsetplot(self.directory, self.name + "_merging_upsetplot", self.dict_upsetplot_reactions,
                              "Intersection of different sources' reactions")


def merging_multirun_first(main_directory):
    """
        Split of major function 'run', first part = creates a merge object for each individual found in the directory
        given and puts them in a list.
    """

    list_objects = []
    for species in utils.get_list_directory(main_directory + "merge"):
        list_objects.append(Merging(species, main_directory))
    return list_objects


def merging_multirun_last(list_objects):
    """
    Split of major function 'run', second part = launches the process on each given object with multiprocessing.
    """

    cpu = len(list_objects)
    p = multiprocessing.Pool(cpu)
    p.map(build_merge_objects, list_objects)


def build_merge_objects(organism):
    """Small function required for the multiprocessing reconstruction."""

    organism.build()


def run(main_directory):
    utils.check_path(main_directory)
    list_objects = merging_multirun_first(main_directory)
    merging_multirun_last(list_objects)


def merging_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("main_directory", help="The path to the main directory where the 'files/' directory is stored",
                        type=str)
    parser.add_argument("-m", "--migrate", help="Take the files previously created by blasting or mpwting modules and "
                                                "store them in a new merge folder, ready to be merged.",
                        action="store_true")
    parser.add_argument("-v", "--verbose", help="Toggle the printing of more information", action="store_true")
    parser.add_argument("-le", "--log_erase", help="Erase the existing log file to create a brand new one",
                        action="store_true")
    args = parser.parse_args()
    return args


def main():
    args = merging_arguments()
    if args.log_erase:
        logging.basicConfig(filename=args.main_directory + '/merging.log', filemode='w', level=logging.INFO,
                            format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(filename=args.main_directory + '/merging.log', level=logging.INFO,
                            format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    if args.verbose:
        logging.getLogger().addHandler(logging.StreamHandler())
    if args.migrate:
        utils.migrate(utils.slash(args.main_directory))
    logging.info("------ Merging module started ------")
    run(utils.slash(args.main_directory))


if __name__ == "__main__":
    main()
    # run("/home/asa/INRAE/These/Dev/tests/")