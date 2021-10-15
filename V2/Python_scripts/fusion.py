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


class Fusion:

    def __init__(self, _name, _metacyc_id_dict, _metacyc_reverse_id_dict, _wd_log, _wd_pgdb, _pwt_reactions_path,
                 _blast_sbml_path):
        self.name = _name
        self.metacyc_id_dict = _metacyc_id_dict
        self.metacyc_reverse_id_dict = _metacyc_reverse_id_dict
        self.wd_log = _wd_log
        self.wd_pgdb = _wd_pgdb
        self.pwt_reactions = utils.get_reactions_PT(self.pwt_reactions_path)
        self.blast_reactions = cobra.io.load_json_model(self.blast_sbml_path)
        self.reactions_list = []
        self.pwt_reactions_list = []
        self.blast_reactions_list = []

    def _get_pwt_reactions(self):
        """Function to get the reaction's ID of Metacyc from a Pathway Tools reconstruction."""

        no_match_list = []

        for reaction in self.pwt_reactions:
            try:
                self.pwt_reactions_list += self.metacyc_id_dict[reaction]
            except KeyError:
                if reaction in self.metacyc_reverse_id_dict.keys():
                    self.pwt_reactions_list.append(reaction)
                else:
                    no_match_list.append(reaction + "\n")
                    print("No match for reaction :", reaction.id, " | ", reaction.name)
        print("Number of reactions from PT model : %i\nNumber of those reactions found in Metacyc : %i\nTotal of "
              "reactions not found : %i "
              % (len(self.pwt_reactions), len(set(self.pwt_reactions_list)), len(no_match_list)))
        no_match_list.append("------\nTotal no match : " + str(len(no_match_list)) + "\n------")
        utils.write_file(self.wd_log, self.name + "_error_reaction_pwt.log", no_match_list)
        set(self.pwt_reactions_list)

    def _get_aracyc_model_reactions(self):
        """Function to get the reaction's ID of Metacyc from a reconstruction with Aracyc as model."""

        no_match_list = []
        for reaction in self.blast_reactions.reactions:
            try:
                self.reactions_list += self.metacyc_id_dict[reaction.name]
            except KeyError:
                if reaction.name in self.metacyc_reverse_id_dict.keys():
                    self.reactions_list.append(reaction.name)
                else:
                    no_match_list.append(reaction + "\n")
                    print("No match for reaction :", reaction.id, " | ", reaction.name)
        print("Nb of reactions from Aracyc model : %i\n\Number of those reactions found in Metacyc : %i\n\Total of "
              "reactions not found : %i"
              % (len(self.blast_reactions.reactions), len(self.reactions_list), len(no_match_list)))
        no_match_list.append("------\nTotal no match : " + str(len(no_match_list)) + "\n------")
        utils.write_file(self.wd_log, self.name + "_error_reaction_blast_reconstruction.log", no_match_list)
        set(self.reactions_list)

    def _correct_gene_rule_reactions(self, verbose=True):
        """Function to correct the gene reaction rule in each reaction taken from Metacyc/Pathway Tools
        to make it fit the organism for which the model is reconstructed."""

        reactions_file = utils.read_file(self.wd_pgdb + "/reactions.dat")
        enzrxns_file = utils.read_file(self.wd_pgdb + "/enzrxns.dat")
        proteins_file = utils.read_file(self.wd_pgdb + "/proteins.dat")

        # First step : gathering the ENZYME-REACTION fields in enzrxns (could be several or none for one ID).
        stop = False
        enzrxns = []
        for reaction_line in reactions_file:
            if "UNIQUE-ID" in reaction_line and "#" not in reaction_line: #TODO
                if stop:
                    break
                try:
                    unique_id = re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', reaction_line).group(0).rstrip()
                    if unique_id == reac.name or unique_id == metacyc_reverse_id_dict[reac.name]:
                        stop = True
                except AttributeError:
                    print("No UNIQUE-ID match for reactions.dat : ", reaction_line)
            if stop and "ENZYMATIC-REACTION " in reaction_line and "#" not in reaction_line:
                try:
                    enzrxns.append(
                        re.search('(?<=ENZYMATIC-REACTION - )[+-]*\w+(.*\w+)*(-*\w+)*', reaction_line).group(0).rstrip())
                except AttributeError:
                    print("No ENZYMATIC-REACTION match for reactions.dat : ", reaction_line)
        if verbose:
            print("%s : %i enzymatic reaction(s) found associated to this reaction." % (reac.name, len(enzrxns)))

        # Second step : getting the corresponding ENZYME for each ENZYME-REACTION.
        stop = False
        geneList = []
        if enzrxns:
            for enzrxn in enzrxns:
                for lineEnzrxn in enzrxns_file:
                    if "UNIQUE-ID" in lineEnzrxn and "#" not in lineEnzrxn:
                        if stop:
                            stop = False
                            break
                        try:
                            unique_id_rxn = re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', lineEnzrxn).group(
                                0).rstrip()
                            if unique_id_rxn == enzrxn:
                                stop = True
                        except AttributeError:
                            print("No UNIQUE-ID match for enzrxns.Dat : ", lineEnzrxn)
                    if stop and "ENZYME " in lineEnzrxn and "#" not in lineEnzrxn:
                        try:
                            enzyme = re.search('(?<=ENZYME - )[+-]*\w+(.*\w+)*(-*\w+)*', lineEnzrxn).group(0).rstrip()
                        except AttributeError:
                            print("No ENZYME match for enzrxns.dat : ", lineEnzrxn)
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
                            geneList.append(
                                re.search('(?<=GENE - )[+-]*\w+(.*\w+)*(-*\w+)*', lineProt).group(0).rstrip())
                        except AttributeError:
                            print("No GENE match for proteins.dat : ", lineProt)
            if verbose:
                print(unique_id, "\n", " or ".join(set(geneList)))
        else:
            pass
        reac.gene_reaction_rule = " or ".join(set(geneList))
        return reac


def pipeline_fusion(corres, aracyc_model_path, metacyc_path, save_path, WDlog, WD_pgdb, verbose=False):
    """Function to merge two models, one from an Aracyc model, the other one from Pathway Tools.
    
    ARGS:
        corres (str) -- the path to the file containing the correspondence 
        between short and long ID of Metacyc's reactions (see utils.metacyc_IDs()).
        aracyc_model_path (str) -- the path to the aracyc based reconstructed model.
        metacyc_path (str) -- the path to the metacyc model.
        save_path (str) -- the path and name to the saving directory.
        WDlog (str) -- the working directory to store the log file containing 
        all the reactions not found in get_pwtools_reac().
        WD_pgdb (str) -- path to the .dat files for the organism (from the mpwt script).
        verbose (boolean) -- print or not enzyme matches.
    """

    metacyc_id_dict, metacyc_reverse_id_dict = utils.corres_dico(corres)
    aracyc_model = cobra.io.load_json_model(aracyc_model_path)
    pwtools_reac_set = get_pwtools_reac(WD_pgdb + "reactions.dat", metacyc_id_dict, metacyc_reverse_id_dict, WDlog,
                                        aracyc_model.id)
    aracyc_reac_set = get_aracyc_model_reac(aracyc_model_path, metacyc_id_dict, metacyc_reverse_id_dict, WDlog,
                                            aracyc_model.id)
    # Here we take away the reactions of metacyc already present in the aracyc model otherwise it will be redundant
    # as we keep all the aracyc model reactions.
    metacyc_reac_set = pwtools_reac_set - set.intersection(pwtools_reac_set, aracyc_reac_set)

    # Then, addition of the metacyc reactions found in the Pathway Tools model to the aracyc model.
    new_model = copy.deepcopy(aracyc_model)

    reactionsFile = utils.read_file(WD_pgdb + "reactions.dat")
    enzrxnsFile = utils.read_file(WD_pgdb + "enzrxns.dat")
    proteinsFile = utils.read_file(WD_pgdb + "proteins.dat")
    metacyc = cobra.io.load_json_model(metacyc_path)
    list_fail = []
    for reac in metacyc_reac_set:
        try:
            if reac[0] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                reac = "_" + reac
            added_reac = copy.deepcopy(metacyc.reactions.get_by_id(reac))
            added_reac_cor = correct_gene_reac(added_reac, reactionsFile, enzrxnsFile, proteinsFile,
                                               metacyc_reverse_id_dict,
                                               verbose)
            new_model.add_reactions([added_reac_cor])
        except KeyError:
            list_fail.append(reac)
    print("Nb of reactions not found in the metacyc model: ", len(list_fail))
    cobra.io.save_json_model(new_model, save_path)
    print("Nb of reactions in the fusioned model : ", len(new_model.reactions))


if __name__ == "__main__":
    # Making of the file with the correspondence between short and long IDs, which will be needed here
    # utils.metacyc_IDs("/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/",
    #                   "/home/asa/INRAE/Logiciels/pathway-tools/metacyc.json")

    ##Tomato
    pipeline_fusion("/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/MetacycCorresIDs.csv",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/blasting/Tomato_Aracyc/Aracyc_genes_Tomato.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/metacyc.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/TomatoFusion.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/mpwt/output/SolLyPHFalse/")

    ##Kiwi
    pipeline_fusion("/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/MetacycCorresIDs.csv",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/blasting/Kiwi_Aracyc/Aracyc_genes_Kiwi.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/metacyc.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/KiwiFusion.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/mpwt/output/ActChPHFalse/")

    ##Cucumber
    pipeline_fusion("/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/MetacycCorresIDs.csv",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/blasting/Cucumber_Aracyc/Aracyc_genes_Cucumber.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/metacyc.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/CucumberFusion.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/mpwt/output/CucSaPHFalse/")

    ##Cherry
    pipeline_fusion("/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/MetacycCorresIDs.csv",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/blasting/Cherry_Aracyc/Aracyc_genes_Cherry.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/metacyc.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/CherryFusion.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/mpwt/output/PruAvPHFalse/")

    ##Camelina
    pipeline_fusion("/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/MetacycCorresIDs.csv",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/blasting/Camelina_Aracyc/Aracyc_genes_Camelina.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/metacyc.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/CamelinaFusion.json",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/Fusion/",
                    "/home/asa/INRAE/Work/FichiersRelancePipeline/mpwt/output/CamSaPHFalse/")
