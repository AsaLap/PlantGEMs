# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used to merge two differents GEMs, one coming from my own model based 
reconstruction and another one from a Pathway Tools reconstruction with the AuReMe's 
mpwt package."""

import cobra
import copy

import utils
import pathwayToolsPrep as PT


###For Pathway Tools models
def get_pwtools_reac(path, dico_matching, dico_matching_rev):
    no_match_cpt = 0
    reac_list = []
    model_reactions = utils.get_reactions_PT(path)
    for reac in model_reactions:
        try:
            reac_list += dico_matching[reac]
        except KeyError:
            try:
                dico_matching_rev[reac]
                reac_list.append(reac)
            except KeyError:
                no_match_cpt += 1
                print("No match for reaction :", reac) ###Store it into a log file
    print("Nb de réactions du modèle PT : %i\n\
Nb de ces réactions trouvées dans Metacyc : %i\n\
Total de réactions non trouvées : %i"
%(len(model_reactions), len(set(reac_list)), no_match_cpt))
    return set(reac_list)


###For Aracyc models
def get_aracyc_model_reac(path, dico_matching, dico_matching_rev):
    no_match_cpt = 0
    reac_list = []
    model = cobra.io.load_json_model(path)
    for reac in model.reactions:
        try:
            reac_list += dico_matching[reac.name]
        except KeyError:
            try:
                dico_matching_rev[reac.name]
                reac_list.append(reac.name)
            except KeyError:
                no_match_cpt += 1
                print("No match for reaction :", reac.id, " | ", reac.name)
    print("Nb de réactions du modèle Aracyc : %i\n\
Nb de ces réactions trouvées dans Metacyc : %i\n\
Total de réactions non trouvées : %i"
    %(len(model.reactions), len(reac_list), no_match_cpt))
    return set(reac_list)


def fusion(corres, pwtools_reac_path, aracyc_model_path, metacyc_path, save_path):
    matching = utils.read_file(corres)
    dico_matching = {}
    dico_matching_rev = {}
    for line in matching:
        if line:
            couple = line.rstrip().split("\t")
            if couple[0] in dico_matching.keys():
                dico_matching[couple[0]].append(couple[1])
            else:
                dico_matching[couple[0]] = [couple[1]]
            dico_matching_rev[couple[1]] = couple[0]

    pwtools_reac_set = get_pwtools_reac(pwtools_reac_path, dico_matching, dico_matching_rev)
    aracyc_reac_set = get_aracyc_model_reac(aracyc_model_path, dico_matching, dico_matching_rev)
    ###Here we take away the reactions of metacyc already present in the aracyc model otherwise 
    ###it will be redundant as we keep all the aracyc model reactions.
    metacyc_reac_set = pwtools_reac_set - set.intersection(pwtools_reac_set, aracyc_reac_set)
    
    ###Then, addition of the metacyc reactions found in the Pathway Tools model to the aracyc model
    new_model = copy.deepcopy(cobra.io.load_json_model(aracyc_model_path))

    metacyc = cobra.io.load_json_model(metacyc_path)
    list_fail = []
    for reac in metacyc_reac_set:
        try:
            if reac[0] in ['0','1','2','3','4','5','6','7','8','9']:
                reac = "_" + reac
            added_reac = copy.deepcopy(metacyc.reactions.get_by_id(reac))
            new_model.add_reactions([added_reac])
        except KeyError:
            list_fail.append(reac)
    print("Nb of reactions not found in the metacyc model: ", len(list_fail))
    cobra.io.save_json_model(new_model, save_path)
    print(len(new_model.reactions))


if __name__ == "__main__":
    fusion("/home/asa/INRAE/Work/fusion_test/correspondanceMetacycIds.tsv",
           "/home/asa/INRAE/Work/fusion_test/reactions.dat",
           "/home/asa/INRAE/Work/blasting_drafts/Tomato_Aracyc/Tomato.json",
           "/home/asa/INRAE/Work/fusion_test/metacyc.json",
           "/home/asa/INRAE/Work/fusion_test/test_fusion.json")