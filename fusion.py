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
import re

import utils


###For Pathway Tools models
def get_pwtools_reac(path, dico_matching, dico_matching_rev, WDlog):
    """Function to get the reaction's ID of Metacyc from a Pathway Tools reconstruction.
    
    ARGS:
        path (str) -- the path to the reactions.dat file.
        dico_matching -- the dictionary containing the correspondance between
        the short-ID and the long-ID of Metacyc's reactions (see the structure 
        in fusion()).
        dico_matching_rev -- same as above but key and values are reversed (see the structure 
        in fusion()).
        WDlog (str) -- the working directory to store the log file containing 
        all the reactions not found.
    RETURN:
        set(reac_list) (set of str) -- the list containing all the reaction's long ID.
    """
      
    no_match_list = []
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
                no_match_list.append(reac + "\n")
                print("No match for reaction :", reac) ###Store it into a log file
    utils.write_file(WDlog, "error_reaction.log", no_match_list)
    print("Number of reactions from PT model : %i\n\
Number of those reactions found in Metacyc : %i\n\
Total of reactions not found : %i"
%(len(model_reactions), len(set(reac_list)), len(no_match_list)))
    return set(reac_list)


###For Aracyc models
def get_aracyc_model_reac(path, dico_matching, dico_matching_rev):
    """Function to get the reaction's ID of Metacyc from an Aracyc model reconstruction.
    
    ARGS:
        path (str) -- the path to the reactions.dat file.
        dico_matching -- the dictionary containing the correspondance between
        the short-ID and the long-ID of Metacyc's reactions (see the construction 
        in utils.py).
        dico_matching_rev -- same as above but key and values are reversed (see the 
        construction in utils.py).
    RETURN:
        set(reac_list) (set of str) -- the list containing all the reaction's long ID.
    """
    
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
    print("Nb of reactions from Aracyc model : %i\n\
Number of those reactions found in Metacyc : %i\n\
Total of reactions not found : %i"
    %(len(model.reactions), len(reac_list), no_match_cpt))
    return set(reac_list)


def metacyc_correspondance(WD, path):
    """Function to make the correspondance file between short and long ID of Metacyc.
    
    ARGS:
    WD (str) -- the directory to save the correspondance file.
        path (str) -- the path to the metacyc model in JSON format.
    """
    
    data = utils.read_json(path)
    res = []
    print(len(data["reactions"]))
    for reac in data["reactions"]:
        long_ID = reac["name"].split("/")[0]
        ###Getting rid of the brackets in the name sometimes!
        reac_pattern = re.compile('([[].*[]])')
        tmp_ID = reac_pattern.sub("", long_ID)
        short_ID = tmp_ID
        ###Regexp to "clean" the metabolite's names
        meta_pattern = re.compile('(_CC[OI]-.*)|(^[_])|([_]\D$)')
        meta_list = []
        # print(len(reac["metabolites"].keys()))
        if len(reac["metabolites"].keys()) != 0:
            for metabolite in reac["metabolites"].keys():
                ###The metabolite are "cleaned" here
                metabolite = meta_pattern.sub("", metabolite)
                len_ID, len_meta = len(tmp_ID), len(metabolite)
                diff = len_ID - len_meta
                ###Small trick to get only the end of the ID removed and not the beginning
                ###(metabolite's names can be in the reaction's name)
                test_ID = tmp_ID[:diff - 1] + tmp_ID[diff - 1:].replace("-" + metabolite, "")
                if len(test_ID) < len(short_ID):
                    short_ID = test_ID
            res.append([short_ID, reac["name"]])
    utils.write_csv(WD, res, "MetacycCorresIDs", "\t")


def correct_gene_reac(reac, reactionsFile, enzrxnsFile, proteinsFile, dico_matching_rev, verbose = True):
    """Function to correct the gene reaction rule in each reaction taken from Metacyc/Pathway Tools
    to make it fit the organism for which the model is reconstructed.
    
    ARGS:
        reac -- the reaction to correct.
        reactionsFile -- the .dat file containing the associated enzymatic reaction(s).
        enzrxnsFile -- the .dat file containing the associated protein name.
        proteinsFile -- the .dat file containing the associated gene name.
        dico_matching_rev -- the dictionary of the correspondance between short IDs 
        and long IDs, in reverse (see construction in utils.py).
        verbose (boolean) -- print or not enzyme matches.
    RETURN:
        reac -- the reaction with the correct gene reaction rule.
    """
    
    ###First step : gathering the ENZYME-REACTION fields in enzrxns (could be several or none for one ID)
    stop = False
    global_stop = False
    enzrxns = []
    log_reactions = ""
    for lineReac in reactionsFile:
        if "UNIQUE-ID" in lineReac and not "#" in lineReac:
            if stop == True:
                stop = False
                break
            try:
                unique_id = re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', lineReac).group(0).rstrip()
                if unique_id == reac.name or unique_id == dico_matching_rev[reac.name]:
                    stop = True
            except AttributeError:
                print("No UNIQUE-ID match for reactions.dat : ", lineReac)
        if stop == True and "ENZYMATIC-REACTION " in lineReac and not "#" in lineReac:
            try:
                enzrxns.append(re.search('(?<=ENZYMATIC-REACTION - )[+-]*\w+(.*\w+)*(-*\w+)*', lineReac).group(0).rstrip())
            except AttributeError:
                print("No ENZYMATIC-REACTION match for reactions.dat : ", lineReac)
    if verbose:
        print("%s : %i enzymatic reaction(s) found associated to this reaction." %(reac.name, len(enzrxns)))
    
    ###Second step : getting the corresponding ENZYME for each ENZYME-REACTION
    stop = False
    geneList = []
    if enzrxns:
        for enzrxn in enzrxns:
            for lineEnzrxn in enzrxnsFile:
                if "UNIQUE-ID" in lineEnzrxn and not "#" in lineEnzrxn:
                    if stop == True:
                        stop = False
                        break
                    try:
                        unique_id_rxn = re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', lineEnzrxn).group(0).rstrip()
                        if unique_id_rxn == enzrxn:
                            stop = True
                    except AttributeError:
                        print("No UNIQUE-ID match for enzrxns.Dat : ", lineEnzrxn)
                if stop == True and "ENZYME " in lineEnzrxn and not "#" in lineEnzrxn:
                    try:
                        enzyme = re.search('(?<=ENZYME - )[+-]*\w+(.*\w+)*(-*\w+)*', lineEnzrxn).group(0).rstrip()
                    except AttributeError:
                        print("No ENZYME match for enzrxns.dat : ", lineEnzrxn)
            ###Third step into the second one : getting the corresponding GENE for each ENZYME
            ###and put it into geneList (which contains all that we're looking for)
            for lineProt in proteinsFile:
                if "UNIQUE-ID " in lineProt and not "#" in lineProt:
                    if stop == True:
                        stop = False
                        break
                    try:
                        unique_id_prot = re.search('(?<=UNIQUE-ID - )[+-]*\w+(.*\w+)*(-*\w+)*', lineProt).group(0).rstrip()
                        if unique_id_prot == enzyme:
                            stop = True
                    except AttributeError:
                        print("No UNIQUE-ID match for proteins.dat : ", lineProt)
                if stop == True and "GENE " in lineProt and not "#" in lineProt:
                    try:
                        geneList.append(re.search('(?<=GENE - )[+-]*\w+(.*\w+)*(-*\w+)*', lineProt).group(0).rstrip())
                    except AttributeError:
                        print("No GENE match for proteins.dat : ", lineProt)
        print(unique_id, "\n", " or ".join(set(geneList)))
    else:
        pass
    reac.gene_reaction_rule = " or ".join(set(geneList))
    return reac


def fusion(corres, aracyc_model_path, metacyc_path, save_path, WDlog, WD_pgdb):
    """Function to merge two models, one from an Aracyc model, the other one from Pathway Tools.
    
    ARGS:
        corres (str) -- the path to the file containing the correspondance 
        between short and long ID of Metacyc's reactions.
        aracyc_model_path (str) -- the path to the aracyc based reconstructed model.
        metacyc_path (str) -- the path to the metacyc model.
        save_path (str) -- the path and name to the saving directory.
        WDlog (str) -- the working directory to store the log file containing 
        all the reactions not found in get_pwtools_reac().
        WD_pgdb (str) -- path to the .dat files in ptools_local for the organism.
    """
    
    dico_matching, dico_matching_rev = utils.corres_dico(corres)
    pwtools_reac_set = get_pwtools_reac(WD_pgdb + "reactions.dat", dico_matching, dico_matching_rev, WDlog)
    aracyc_reac_set = get_aracyc_model_reac(aracyc_model_path, dico_matching, dico_matching_rev)
    ###Here we take away the reactions of metacyc already present in the aracyc model otherwise 
    ###it will be redundant as we keep all the aracyc model reactions.
    metacyc_reac_set = pwtools_reac_set - set.intersection(pwtools_reac_set, aracyc_reac_set)
    
    ###Then, addition of the metacyc reactions found in the Pathway Tools model to the aracyc model
    new_model = copy.deepcopy(cobra.io.load_json_model(aracyc_model_path))

    reactionsFile = utils.read_file(WD_pgdb + "reactions.dat")
    enzrxnsFile = utils.read_file(WD_pgdb + "enzrxns.dat")
    proteinsFile = utils.read_file(WD_pgdb + "proteins.dat")
    metacyc = cobra.io.load_json_model(metacyc_path)
    list_fail = []
    for reac in metacyc_reac_set:
        try:
            if reac[0] in ['0','1','2','3','4','5','6','7','8','9']:
                reac = "_" + reac
            added_reac = copy.deepcopy(metacyc.reactions.get_by_id(reac))
            added_reac = correct_gene_reac(added_reac, reactionsFile, enzrxnsFile, proteinsFile, dico_matching_rev)
            new_model.add_reactions([added_reac])
        except KeyError:
            list_fail.append(reac)
    print("Nb of reactions not found in the metacyc model: ", len(list_fail))
    cobra.io.save_json_model(new_model, save_path)
    print("Nb of reactions in the fusioned model : ", len(new_model.reactions))


if __name__ == "__main__":
    ###Tomato
    # fusion("/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv",
    #        "/home/asa/INRAE/Work/blasting_drafts/Tomato_Aracyc/Tomato.json",
    #        "/home/asa/INRAE/Work/Fusion/metacyc.json",
    #        "/home/asa/INRAE/Work/Fusion/TomatoFusionGenes.json",
    #        "/home/asa/INRAE/Work/Fusion/",
    #        "/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/sollyphfalsecyc/1.0/data/")

    ###Debug
    # WD_pgdb = "/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/sollyphfalsecyc/1.0/data/"
    # reactionsFile = reactionsFile = utils.read_file(WD_pgdb + "reactions.dat")
    # enzrxnsFile = utils.read_file(WD_pgdb + "enzrxns.dat")
    # proteinsFile = utils.read_file(WD_pgdb + "proteins.dat")
    # metacyc = cobra.io.load_json_model("/home/asa/INRAE/Work/Fusion/metacyc.json")
    # dico_matching, dico_matching_rev = utils.corres_dico("/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv")
    
    # reac = "RXN-14093"
    # added_reac = copy.deepcopy(metacyc.reactions.get_by_id(reac))
    # added_reac = correct_gene_reac(added_reac, reactionsFile, enzrxnsFile, proteinsFile, dico_matching_rev)

    ###Kiwi
    # fusion("/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv",
    #        "/home/asa/INRAE/Work/blasting_drafts/Kiwi_Aracyc/Kiwi.json",
    #        "/home/asa/INRAE/Work/Fusion/metacyc.json",
    #        "/home/asa/INRAE/Work/Fusion/KiwiFusionGenes.json",
    #        "/home/asa/INRAE/Work/Fusion/",
    #        "/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/actchphfalsecyc/1.0/data/")
 
    ###Cucumber
    # fusion("/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv",
    #        "/home/asa/INRAE/Work/blasting_drafts/Cucumber_Aracyc/Cucumber.json",
    #        "/home/asa/INRAE/Work/Fusion/metacyc.json",
    #        "/home/asa/INRAE/Work/Fusion/CucumberFusionGenes.json",
    #        "/home/asa/INRAE/Work/Fusion/",
    #        "/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/cucsaphfalsecyc/1.0/data/")

    ###Cherry
    # fusion("/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv",
    #        "/home/asa/INRAE/Work/blasting_drafts/Cherry_Aracyc/Cherry.json",
    #        "/home/asa/INRAE/Work/Fusion/metacyc.json",
    #        "/home/asa/INRAE/Work/Fusion/CherryFusionGenes.json",
    #        "/home/asa/INRAE/Work/Fusion/",
    #        "/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/pruavphfalsecyc/1.0/data/")

    ###Camelina
    # fusion("/home/asa/INRAE/Work/Fusion/MetacycCorresIDs.tsv",
    #        "/home/asa/INRAE/Work/blasting_drafts/Camelina_Aracyc/Camelina.json",
    #        "/home/asa/INRAE/Work/Fusion/metacyc.json",
    #        "/home/asa/INRAE/Work/Fusion/CamelinaFusionGenes.json",
    #        "/home/asa/INRAE/Work/Fusion/",
    #        "/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/camsaphfalsecyc/1.0/data/")
    
    # metacyc = cobra.io.load_json_model("/home/asa/INRAE/Work/Fusion/metacyc.json")
    # fileTest = utils.read_file("/home/asa/INRAE/Work/Fusion_repair/liste_reac.txt")
    # for reac in fileTest:
    #     old_reac = copy.deepcopy(metacyc.reactions.get_by_id(reac.rstrip()))
    #     correct_reac = correct_gene_reac("/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/sollyphfalsecyc/1.0/data/",
    #                     old_reac)