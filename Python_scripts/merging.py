# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used to merge two different GEMs, one coming from my own model based
reconstruction and another one from a Pathway Tools reconstruction with the AuReMe's 
mpwt package."""

import cobra
import copy
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
        self.pwt_reactions_id_list = []
        self.pwt_metacyc_reactions_id_list = []
        self.pwt_metacyc_no_match_id_list = []
        self.pwt_metacyc_long_id_list = []
        self.json_reactions_list = []
        self.sbml_reactions_list = []
        self.merged_model = cobra.Model(self.name, name=self.name + "_PlantGEMs_" + str(date.today()))
        self.metacyc_model = cobra.io.load_json_model(self.main_directory + "merge/metacyc.json")
        self.metacyc_ids_file_path = utils.find_file(self.main_directory + "merge/", "metacyc_ids", "tsv")
        if not os.path.isfile(self.metacyc_ids_file_path):
            utils.get_metacyc_ids(self._find_sbml_model(self.directory))
        self.metacyc_matching_id_dict, self.metacyc_matching_id_dict_reversed = \
            utils.build_correspondence_dict(self.metacyc_ids_file_path)

    def _get_pwt_reactions(self):
        self.pwt_reactions_id_list = utils.get_pwt_reactions(self.directory + "/reactions.dat")

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
                    self.json_reactions_list += cobra.io.load_json_model(self.directory + json_model).reactions
            if extension == "sbml":
                for sbml_model in list_networks:
                    self.sbml_reactions_list += cobra.io.read_sbml_model(self.directory + sbml_model).reactions
        else:
            print("------\nNo " + extension + " file of draft network or model found.\n------")

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
        for reaction_line in reactions_file:
            if "UNIQUE-ID" in reaction_line and "#" not in reaction_line:
                if stop:
                    break
                try:
                    unique_id = re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', reaction_line).group(0).rstrip()
                    if unique_id == reaction.name or unique_id == self.metacyc_matching_id_dict_reversed[reaction.name]:
                        stop = True
                except AttributeError:
                    print("No UNIQUE-ID match for reactions.dat : ", reaction_line)
            if stop and "ENZYMATIC-REACTION " in reaction_line and "#" not in reaction_line:
                try:
                    enzrxns.append(re.search('(?<=ENZYMATIC-REACTION - )[+-]*\w+(.*\w+)*(-*\w+)*',
                                             reaction_line).group(0).rstrip())
                except AttributeError:
                    print("No ENZYMATIC-REACTION match for reactions.dat : ", reaction_line)
        if verbose:
            print("%s : %i enzymatic reaction(s) found associated to this reaction." % (reaction.name, len(enzrxns)))

        # Second step : getting the corresponding ENZYME for each ENZYME-REACTION.
        stop = False
        gene_list = []
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
                            print("No UNIQUE-ID match for enzrxns.Dat : ", line_enzrxn)
                    if stop and "ENZYME " in line_enzrxn and "#" not in line_enzrxn:
                        try:
                            enzyme = re.search('(?<=ENZYME - )[+-]*\w+(.*\w+)*(-*\w+)*', line_enzrxn).group(0).rstrip()
                        except AttributeError:
                            print("No ENZYME match for enzrxns.dat : ", line_enzrxn)
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
                            print("No UNIQUE-ID match for proteins.dat : ", lineProt)
                    if stop and "GENE " in lineProt and "#" not in lineProt:
                        try:
                            gene_list.append(
                                re.search('(?<=GENE - )[+-]*\w+(.*\w+)*(-*\w+)*', lineProt).group(0).rstrip())
                        except AttributeError:
                            print("No GENE match for proteins.dat : ", lineProt)
            if verbose:
                print(unique_id, "\n", " or ".join(set(gene_list)))
        else:
            pass
        reaction.gene_reaction_rule = " or ".join(set(gene_list))
        return reaction

    def _merge(self, trust_model=True, verbose=False):
        """Function to merge models, either from a model-based reconstruction (blasting module) or from
        Pathway Tools's Pathologic software.

        ARGS:
            trust_model (boolean) -- True if you trust the model(s) you gave in input (blasting/merging), taking all of
                their reactions even if they are not found in the Metacyc network.
                False if you only want the Metacyc reactions.
            verbose (boolean) -- print or not the protein matches in the terminal (won't change the logs).
        """

        for reaction in self.pwt_reactions_id_list:
            if reaction[0] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                reaction = "_" + reaction
            try:
                added_reactions = copy.deepcopy(self.metacyc_model.reactions.get_by_id(reaction))
                added_reactions_corrected = self._correct_pwt_gene_reaction_rule(added_reactions, verbose)
                self.merged_model.add_reactions([added_reactions_corrected])
            except KeyError:
                pass
        if trust_model:
            for reaction in self.json_reactions_list:
                self.merged_model.add_reactions([reaction])
            for reaction in self.sbml_reactions_list:
                self.merged_model.add_reactions([reaction])
        # else:  # TODO

    def run(self):
        if utils.check_path(self.directory + "reactions.dat"):
            self._get_pwt_reactions()
            self._search_metacyc_reactions_ids()
        else:
            print("------\nNo .dat files found, proceeding with drafts and models only.\n------")
        if utils.check_path(self.directory, sys_exit=True):
            self._get_networks_reactions("json")
            self._get_networks_reactions("sbml")
        self._merge()
        cobra.io.save_json_model(self.merged_model, self.directory + self.name + "_merged.json")

    # TODO : Check for PHOSCHOL-RXN if specificity is kept during merging process (cucumis_sativus)


if __name__ == '__main__':
    main_directory = "/home/asa/INRAE/These/Tests/"
    test = Merging("chimera", main_directory)
    # test = Merging("actinidia_chinensis", main_directory)
    # test = Merging("cucumis_sativus", main_directory)
    test.run()
