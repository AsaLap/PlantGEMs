# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020 - Reworked on November 2023
"""This file is used for the gap filling of the previously reconstructed networks."""

import cobra
import copy
from pyasp.term import *
from meneco import run_meneco

import module
import utils


class Menecoing(module.Module):

    def __init__(self, _name, _main_directory, _repairnet, _draftnet, _seeds, _targets):
        """
        Args:
            _repairnet (str): name of the repair network, must be stored in the 'files' directory.
            _draftnet (str): name of the subject network to fill, must be stored in the 'meneco' directory.
            _seeds (str): name of the seeds file, must be stored in the 'files' directory.
            _targets (str): name of the targets file, must be stored in the 'files' directory.
        """

        super().__init__(_name, _main_directory)
        self.repairnet_path = self.main_directory + _repairnet
        self.draftnet_path = self.main_directory + _draftnet
        self.seeds_path = self.main_directory + _seeds
        self.targets_path = self.main_directory + _targets
        self.gap_filling_output = None

    def pipeline_gap_filling(self, enumeration=False, json=True):
        """The main function to make all the process.

        Args::
            seeds (str): the filename of the seeds.
            targets (str): the filename of the targets.
            enumeration (boolean): Meneco choice to list all the reactions found or not.
            json (boolean): Meneco choice of getting the result as a JSON (not working).
        """

        # Removed on the V2 version because it seems not to be needed
        # clean_draft = utils.clean_sbml(WD, draft)
        # clean_repair = utils.clean_sbml(WD, repair)
        self.gap_filling_output = run_meneco(draftnet=self.draftnet_path,
                            seeds=self.seeds_path,
                            targets=self.targets_path,
                            repairnet=self.repairnet_path,
                            enumeration=enumeration,
                            json_output=json)
        # /!\ This next function line will add every reaction without manual curation,
        # /!\ uncommenting it is not advised, as it is unspecific.
        # add_filled_reactions(WD, result['Union of cardinality minimal completions'], repair, clean_draft)

    def add_filled_reactions(WD, reacs, repair, draft, json=False):
        """Function to add the retrieved reactions to the model.

        Args::
            WD (str): the path to the working directory where there are the repair
            model and the draft.
            reacs (list of str): the list of reaction names.
            repair: the repair model file (SBML format).
            draft: the draft model file (SBML format).
            json (boolean): True if you want a save of the model in the JSON
            format instead of SBML (default).
        """

        repair_model = cobra.io.read_sbml_model(WD + repair)
        draft_model = cobra.io.read_sbml_model(WD + draft)
        old_size = len(draft_model.reactions)
        reac_to_blast = []
        for reac in reacs:
            reac = utils.cobra_compatibility(reac.rstrip())
            try:
                new_reac = copy.deepcopy(repair_model.reactions.get_by_id(reac))
                reac_to_blast.append([new_reac.id, new_reac.gene_reaction_rule])
                new_reac.gene_reaction_rule = ""
                draft_model.add_reactions([new_reac])
            except KeyError:
                print("No match for this reaction : ", reac)
        print("Number of reactions of the unrepaired model : %i\nNumber of reactions of the repaired model : %i" % (
        old_size, len(draft_model.reactions)))
        if json:
            cobra.io.save_json_model(draft_model, WD + "filled_" + draft.split(".")[0] + ".json")
        else:
            cobra.io.write_sbml_model(draft_model, WD + "filled_" + draft.split(".")[0] + ".sbml")
        utils.write_csv(WD, reac_to_blast, "rxns_to_blast_" + draft_model.id, "\t")


# def stats(res):
#     """NOT USED.
#     Function that transforms the dictionary of results from the utils.count_letters()
#     function into percentages, print and return them.
#
#     Args::
#         res: a dictionary from 'utils.count_letters()' (see count() for the structure).
#     Returns::
#         str_res (str): a string containing the computed information line by line.
#     """
#
#     total = 0
#     str_res = ""
#     for key in res.keys():
#         total += res[key]
#     for key in res.keys():
#         str_res += str(key) + "\t" + str(res[key]) + "\t" + str(round(res[key] * 100 / total, 2)) + "\n"
#         print("%s : %.0f = %.2f percent" % (key, res[key], res[key] * 100 / total))
#     return str_res


# def make_info(self):
#     """Function to organize all the previous ones to write a .txt file for an organism.
#     Args::
#         WD (str): the working directory.
#         DNA_file (str): the name of the fasta file DNA.
#         RNA_file (str): the name of the fasta file of only the amino_acids.
#         prot_file (str): the name of the fasta file of the proteins.
#         name (str): the name for the output .txt file.
#     """
#
#     letters_prot = "ACDEFGHIKLMNPQRSTVWY*"
#     letters_dna = "CGTAN"
#     str_results = "DNA stats :\n"
#     dna_file_cleaned = clean(WD + DNA_file)
#     results_dna = count(dna_file_cleaned, letters_dna)
#     str_results += stats(results_dna)
#     rna_file_cleaned = clean(WD + RNA_file)
#     results_rna = count(rna_file_cleaned, letters_dna)
#     str_results += "RNA stats :\n"
#     str_results += stats(results_rna)
#     prot_file_clean = clean(WD + prot_file)
#     results_prot = count(prot_file_clean, letters_prot)
#     str_results += "Prot stats :\n"
#     str_results += stats(results_prot)
#     utils.write_file(WD, name + "_stats.txt", str_results)
#     print(str_results)


if __name__ == "__main__":
    # vitis_json = cobra.io.load_json_model("/home/asa/INRAE/These/Bioinfo/Info/Tests/menecoing_2023_15_11/vitis_vinifera_merged.json")
    # cobra.io.write_sbml_model(vitis_json, "/home/asa/INRAE/These/Bioinfo/Info/Tests/menecoing_2023_15_11/vitis_vinifera_merged.sbml")
    vitis = Menecoing("vitis_vinifera",
                      "/home/asa/INRAE/These/Bioinfo/Info/Tests/menecoing_2023_15_11/",
                      "metacyc.sbml",
                      "vitis_vinifera_merged.sbml",
                      "seedsPlants.sbml",
                      "targetsTomato.sbml")
    vitis.pipeline_gap_filling()
    del vitis.gap_filling_output['Draft network file'], vitis.gap_filling_output['Seeds file'], vitis.gap_filling_output['Targets file'], vitis.gap_filling_output['Repair db file']
    utils.write_dict("/home/asa/INRAE/These/Bioinfo/Info/Tests/menecoing_2023_15_11/meneco_output.txt", vitis.gap_filling_output)

