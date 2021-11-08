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
import re

import utils


class Merging:

    def __init__(self, _name, _metacyc_correspondence_file_path, _wd_log, _wd_pgdb, _pwt_reactions_path,
                 _blast_sbml_path, _metacyc_sbml_path):
        self.name = _name
        self.metacyc_matching_id_dict, self.metacyc_matching_id_dict_reversed = \
            utils.build_correspondence_dict(_metacyc_correspondence_file_path)
        self.wd_log = _wd_log
        self.wd_pgdb = _wd_pgdb
        self.pwt_reactions = utils.get_pwt_reactions(_pwt_reactions_path)
        self.blast_model = cobra.io.load_json_model(_blast_sbml_path)
        self.metacyc_model = cobra.io.load_json_model(_metacyc_sbml_path)
        self.reactions_list = []
        self.pwt_reactions_list = []
        self.blast_reactions_list = []

    def _get_pwt_reactions(self):
        """Function to get the reaction's ID of Metacyc from a Pathway Tools reconstruction."""

        no_match_list = []

        for reaction in self.pwt_reactions:
            try:
                self.pwt_reactions_list += self.metacyc_matching_id_dict[reaction]
            except KeyError:
                if reaction in self.metacyc_matching_id_dict_reversed.keys():
                    self.pwt_reactions_list.append(reaction)
                else:
                    no_match_list.append(reaction + "\n")
                    print("No match for reaction :", reaction.id, " | ", reaction.name)
        print("Number of reactions from PT model : %i\nNumber of those reactions found in Metacyc : %i\nTotal of "
              "reactions not found : %i "
              % (len(self.pwt_reactions), len(set(self.pwt_reactions_list)), len(no_match_list)))
        no_match_list.append("------\nTotal no match : " + str(len(no_match_list)) + "\n------")
        utils.write_file(self.wd_log + self.name + "_error_reaction_pwt.log", no_match_list)

    def _get_aracyc_model_reactions(self):
        """Function to get the reaction's ID of Metacyc from a reconstruction with Aracyc as model."""

        no_match_list = []
        for reaction in self.blast_model.reactions:
            try:
                self.reactions_list += self.metacyc_matching_id_dict[reaction.name]
            except KeyError:
                if reaction.name in self.metacyc_matching_id_dict_reversed.keys():
                    self.reactions_list.append(reaction.name)
                else:
                    no_match_list.append(reaction + "\n")
                    print("No match for reaction :", reaction.id, " | ", reaction.name)
        print("Nb of reactions from Aracyc model : %i\n\Number of those reactions found in Metacyc : %i\n\Total of "
              "reactions not found : %i"
              % (len(self.blast_model.reactions), len(self.reactions_list), len(no_match_list)))
        no_match_list.append("------\nTotal no match : " + str(len(no_match_list)) + "\n------")
        utils.write_file(self.wd_log + self.name + "_error_reaction_blast_reconstruction.log", no_match_list)

    def _correct_gene_rule_reactions(self, reaction, verbose=True):
        """Function to correct the gene reaction rule in each reaction taken from Metacyc/Pathway Tools
        to make it fit the organism for which the model is reconstructed.
        
        Args:
            reaction: the reaction's gene reaction rule to change.
        Returns the corrected gene reaction rule.
        """

        reactions_file = utils.read_file_listed(self.wd_pgdb + "/reactions.dat")
        enzrxns_file = utils.read_file_listed(self.wd_pgdb + "/enzrxns.dat")
        proteins_file = utils.read_file_listed(self.wd_pgdb + "/proteins.dat")

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

    def merge(self, verbose=False):
        """Function to merge two models, one from an Aracyc model, the other one from Pathway Tools.
        
        ARGS:
            corres (str) -- the path to the file containing the correspondence 
            between short and long ID of Metacyc's reactions (see utils.metacyc_ids()).
            aracyc_model_path (str) -- the path to the aracyc based reconstructed model.
            metacyc_path (str) -- the path to the metacyc model. (=self.pwt_reactions/self.pwt_reactions_list)
            save_path (str) -- the path and name to the saving directory.
            WDlog (str) -- the working directory to store the log file containing 
            all the reactions not found in get_pwtools_reac().
            WD_pgdb (str) -- path to the .dat files for the organism (from the mpwt script).
            verbose (boolean) -- print or not enzyme matches.
        """

        self._get_pwt_reactions()
        self._get_aracyc_model_reactions()
        # Here we take away the reactions of the blast reconstruction because we will deep-copy them as they are
        # already correctly linked to the good genes and we just need to correct the gene reaction rule from the
        # pathway tools reconstruction.
        metacyc_reactions_set = set(self.pwt_reactions_list) - set.intersection(set(self.pwt_reactions_list),
                                                                                set(self.blast_reactions_list))
        new_model = copy.deepcopy(self.blast_model)
        list_fail = []
        for reaction in metacyc_reactions_set:
            try:
                if reaction[0] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                    reaction = "_" + reaction
                added_reactions = copy.deepcopy(self.metacyc_model.reactions.get_by_id(reaction))
                added_reactions_corrected = self._correct_gene_rule_reactions(added_reactions, verbose)
                new_model.add_reactions([added_reactions_corrected])
            except KeyError:
                list_fail.append(reaction)
        print("Nb of reactions not found in the metacyc model: ", len(list_fail))
        cobra.io.save_json_model(new_model, self.wd_log + "/../" + self.name + "_merged.json")
        print("Nb of reactions in the merged model : ", len(new_model.reactions))
