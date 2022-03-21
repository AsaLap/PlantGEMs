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
        self.pwt_metacyc_no_match_id_list = []
        self.pwt_metacyc_long_id_list = []
        self.json_reactions_list = []
        self.sbml_reactions_list = []
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
        for reaction in self.pwt_reactions_id_list:
            try:
                var = self.metacyc_matching_id_dict[reaction]
                self.pwt_metacyc_reactions_id_list.append(reaction)
            except KeyError:
                try:
                    self.pwt_metacyc_reactions_id_list.append(self.metacyc_matching_id_dict_reversed[
                                                                  reaction])  # despécialisation de la réaction (long à court)
                    self.pwt_metacyc_long_id_list.append(reaction)
                except KeyError:
                    self.pwt_metacyc_no_match_id_list.append(reaction)

    def _get_networks_reactions(self, extension):
        list_networks = utils.find_files(self.directory, extension)
        if list_networks:
            if extension == "json":
                for json_model in list_networks:
                    logging.info("JSON model network found : {}".format(json_model))
                    self.json_reactions_list += cobra.io.load_json_model(self.directory + json_model).reactions
            if extension == "sbml":
                for sbml_model in list_networks:
                    logging.info("SBML model network found : {}".format(sbml_model))
                    self.sbml_reactions_list += cobra.io.read_sbml_model(self.directory + sbml_model).reactions
        else:
            logging.info("------\nNo {} file of draft network or model found for {}\n------".format(extension,
                                                                                                    self.name))

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
                    count += 1
                    self.pwt_reactions_id_list.append(
                        re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', line).group(0).rstrip())
                except AttributeError:
                    count -= 1
                    logging.error("No match for : {}".format(line))
        logging.info("{} : Number of reactions found in the Pathway Tools reconstruction files : {}".format(self.name,
                                                                                                            count))

    def _correct_pwt_gene_reaction_rule(self, reaction, verbose=True):
        """Function to correct the gene reaction rule in each reaction taken from Metacyc/Pathway Tools
        to make it fit the organism for which the model is reconstructed.

        Args:
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
                    logging.error("No UNIQUE-ID match for reactions.dat : {}".format(reaction_line))
            if stop and "ENZYMATIC-REACTION " in reaction_line and "#" not in reaction_line:
                try:
                    enzrxns.append(re.search('(?<=ENZYMATIC-REACTION - )[+-]*\w+(.*\w+)*(-*\w+)*',
                                             reaction_line).group(0).rstrip())
                except AttributeError:
                    logging.error("No ENZYMATIC-REACTION match for reactions.dat : {}".format(reaction_line))
        logging.info("{} : {} enzymatic reaction(s) found associated to this reaction.".format(reaction.name,
                                                                                               len(enzrxns)))

        # Second step : getting the corresponding ENZYME for each ENZYME-REACTION.
        stop = False
        gene_list = []
        enzyme = ""
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
                            logging.error("No UNIQUE-ID match for enzrxns.Dat : {}".format(line_enzrxn))
                    if stop and "ENZYME " in line_enzrxn and "#" not in line_enzrxn:
                        try:
                            enzyme = re.search('(?<=ENZYME - )[+-]*\w+(.*\w+)*(-*\w+)*', line_enzrxn).group(0).rstrip()
                        except AttributeError:
                            logging.error("No ENZYME match for enzrxns.dat : {}".format(line_enzrxn))
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
                # logging.info(unique_id, "\n", " or ".join(set(gene_list))) Not needed, too much flooding
        else:
            logging.error("{} : No corresponding enzrxns entry for {}".format(self.name, unique_id))
            pass
        reaction.gene_reaction_rule = " or ".join(set(gene_list))
        return reaction

    def _merge(self, verbose=False):
        """Function to merge models, either from a model-based reconstruction (blasting module) or from
        Pathway Tools's Pathologic software.

        ARGS:
            trust_model (boolean) -- True if you trust the model(s) you gave in input (blasting/merging), taking all of
                their reactions even if they are not found in the Metacyc network.
                False if you only want the Metacyc reactions.
            verbose (boolean) -- print or not the protein matches in the terminal (won't change the logs).
        """

        count = 0
        for reaction in self.pwt_reactions_id_list:
            if count % 100 == 0:
                print("%s : gene correction %s ou of %s" % (self.name, str(count),
                                                            str(len(self.pwt_reactions_id_list))))
            count += 1
            if reaction[0] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                reaction = "_" + reaction
            try:
                added_reactions = copy.deepcopy(self.metacyc_model.reactions.get_by_id(reaction))
                added_reactions_corrected = self._correct_pwt_gene_reaction_rule(added_reactions, verbose)
                self.merged_model.add_reactions([added_reactions_corrected])
            except KeyError:
                logging.info("No match in gene correction for {} ({}'s correction)".format(reaction, self.name))
                pass
        for reaction in self.json_reactions_list:
            self.merged_model.add_reactions([reaction])
        for reaction in self.sbml_reactions_list:
            self.merged_model.add_reactions([reaction])

    def build(self):
        if utils.check_path(self.directory + "reactions.dat"):
            self._get_pwt_reactions()
            self._search_metacyc_reactions_ids()
        else:
            logging.info("------\nNo .dat files found, proceeding with drafts and models only.\n------")
        if utils.check_path(self.directory, sys_exit=True):
            self._get_networks_reactions("json")
            self._get_networks_reactions("sbml")
        self._merge()
        cobra.io.save_json_model(self.merged_model, self.directory + self.name + "_merged.json")
        # TODO : log the lists of reactions that are in both reconstruction method


def merging_multirun_first(main_directory):
    list_objects = []
    for species in utils.get_list_directory(main_directory + "merge"):
        list_objects.append(Merging(species, main_directory))
    return list_objects


def merging_multirun_last(list_objects):
    """
    Split of major function 'run', second part = launching the process on each given object with multiprocessing.
    """

    cpu = len(list_objects)
    p = multiprocessing.Pool(cpu)
    p.map(build_merge_objects, list_objects)


def build_merge_objects(organism):
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


# TODO : Check for PHOSCHOL-RXN if specificity is kept during merging process (cucumis_sativus)
if __name__ == "__main__":
    main()
