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

import module
import utils


class Merging(module.Module):

    def __init__(self, _name, _main_directory):
        """Explanations here"""
        super().__init__(_name, _main_directory)
        self.directory = self.main_directory + "merge/" + self.name + "/"
        self.pwt_reactions_id_list = []
        self.json_reactions_list = []
        self.sbml_reactions_list = []

        self.metacyc_ids_file_path = utils.find_file(self.main_directory + "merge/", "metacyc_ids", "tsv")
        if not os.path.isfile(self.metacyc_ids_file_path):
            utils.get_metacyc_ids(self._find_sbml_model(self.directory))
        self.metacyc_matching_id_dict, self.metacyc_matching_id_dict_reversed = \
            utils.build_correspondence_dict(self.metacyc_ids_file_path)

    def _get_pwt_reactions(self):
        self.pwt_reactions_id_list = utils.get_pwt_reactions(self.directory + "/reactions.dat")

    def _get_json_models_reactions(self):
        list_json_model = utils.find_files(self.directory, "json")
        for json_model in list_json_model:
            self.json_reactions_list += cobra.io.load_json_model(self.directory + json_model).reactions

    def _get_sbml_models_reactions(self):
        list_sbml_model = utils.find_files(self.directory, "sbml")
        for sbml_model in list_sbml_model:
            self.sbml_reactions_list += cobra.io.read_sbml_model(self.directory + sbml_model).reactions

    def _search_metacyc_reactions_ids(self):
        reactions_ids = []
        long_reactions_ids = []
        no_match_ids = []
        for reaction in self.pwt_reactions_id_list:
            try:
                var = self.metacyc_matching_id_dict[reaction]
                reactions_ids.append(reaction)
            except KeyError:
                try:
                    reactions_ids.append(self.metacyc_matching_id_dict_reversed[
                                         reaction])  # despécialisation de la réaction (long à court)
                    long_reactions_ids.append(reaction)
                except KeyError:
                    no_match_ids.append(reaction)
        return reactions_ids, long_reactions_ids, no_match_ids

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

    # def merge(self, verbose=False):
    #     """Function to merge models, either from a model-based reconstruction (blasting module) or from
    #     Pathway Tools's Pathologic software.
    #
    #     ARGS:
    #         verbose (boolean) -- print or not protein matches.
    #     """
    #
    #     self._get_pwt_reactions()
    #     self._get_blast_model_reactions()
    #     # Here we take away the reactions of the blast reconstruction because we will deep-copy them as they are
    #     # already correctly linked to the good genes and we just need to correct the gene reaction rule from the
    #     # pathway tools reconstruction.
    #     metacyc_reactions_set = set(self.pwt_reactions_list) - set.intersection(set(self.pwt_reactions_list),
    #                                                                             set(self.blast_reactions_list))
    #     new_model = copy.deepcopy(self.blast_model)
    #     list_fail = []
    #     for reaction in metacyc_reactions_set:
    #         try:
    #             if reaction[0] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
    #                 reaction = "_" + reaction
    #             added_reactions = copy.deepcopy(self.metacyc_model.reactions.get_by_id(reaction))
    #             added_reactions_corrected = self._correct_gene_reaction_rule(added_reactions, verbose)
    #             new_model.add_reactions([added_reactions_corrected])
    #         except KeyError:
    #             list_fail.append(reaction)
    #     print("Nb of reactions not found in the metacyc model: ", len(list_fail))
    #     cobra.io.save_json_model(new_model, self.wd_log + "/../" + self.name + "_merged.json")
    #     print("Nb of reactions in the merged model : ", len(new_model.reactions))


if __name__ == '__main__':
    main_directory = "/home/asa/INRAE/These/Tests/"
    test = Merging("actinidia_chinensis", main_directory)
