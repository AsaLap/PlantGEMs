# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used as a tool to gather data to make graph, 
do some calculations on my newly made drafts and so on..."""

import matplotlib.pyplot as plt
from statistics import mean
import csv
from upsetplot import plot
from upsetplot import from_memberships
import copy
import cobra
from supervenn import supervenn
import re

import blasting
import utils


def read_file(WD, file):
    """Function to load a file written under the previous model of save."""
    res = []
    f = open(WD + file, "r")
    for line in f.readlines():
        res.append(line.rstrip())
    return res


def graph_identity(data):
    """Function to gather the identity value of the blastp results
    (obtained with blasting.py) and prepare them for a plot.
    
    ARGS:
        data -- the dictionary containing the genes and the blastp result for each.
    RETURN:
        x -- list of x-axis values depending on the treshold value.
        y -- list of y-axis values depending on the treshold value.
    """
    
    res_plot = {}
    treshold = 5
    while treshold <= 100:
        count = 0
        for key in data.keys():
            for res in data[key]:
                if float(res.split(",")[6]) >= treshold:
                    count+=1
        res_plot[treshold] = count
        treshold+=5
    x, y = list(res_plot.keys()), list(res_plot.values())
    return x, y


def get_score(data):
    """Function that gathers the Score of all the blastp
    results and returns a list of these values."""
    
    list_value = []
    for key in data.keys():
        for res in data[key]:
            list_value.append(float(res.split(",")[7]))
    return list_value


def get_e_value(data):
    """Function that gathers the E-Value of all the blastp
    results and returns a list of these values."""
    
    list_value = []
    for key in data.keys():
        for res in data[key]:
            list_value.append(float(res.split(",")[8]))
    return list_value


def get_bit_score(data):
    """Function that gathers the Bit-Score of all the blastp
    results and returns a list of these values."""
    
    list_value = []
    for key in data.keys():
        for res in data[key]:
            list_value.append(float(res.split(",")[9]))
    return list_value


def get_all_scores(organism, data):
    """Function that gathers all the previous scores of all the blastp
    results and returns a list of these values, with the organism name as first value.
    This function is used to make graph with the R script 'regressions.R'"""
    
    list_value = []
    for key in data.keys():
        for res in data[key]:
            spl = res.split(",")
            list_value.append([
                organism,
                spl[2],
                float(spl[6]),
                float(spl[7]),
                float(spl[8]),
                float(spl[9])])
    return list_value


def make_upsetplot(WD, data, name):
    """Function to make an UpSetPlot.
    Need this three other functions : similarity_count(), get_clusters(), get_sub_clusters().
    
    ARGS:
        WD -- the working directory to save the result.
        data -- the dictionary containing the organisms as keys
        and the genes/reactions/others to treat for the UpSetPlot.
    """
    
    clusters = get_clusters(list(data.keys()))
    [clusters.insert(0,[key]) for key in data.keys()]
    count = []
    log = ""
    for c in clusters:
        others = list(data.keys())
        listInter = []
        for x in c:
            others.remove(x)
            listInter.append(set(data[x]))
        cluster_data, sim_count = similiraty_count(data, listInter, others)
        count.append(sim_count)
        for i in c:
            log += i + " "
        log += " :\n"
        for i in cluster_data:
            log += str(i) + "\n"
        log += "\n------\n\n"
    utils.write_file(WD, "logUpsetplot.txt", log)
    my_upsetplot = from_memberships(clusters, count)
    plot(my_upsetplot, show_counts = '%d', totals_plot_elements = 3)
    plt.suptitle("Unique intersection for each possible cluster")
    plt.savefig(WD + name + ".pdf")
    plt.show()


def similiraty_count(data, args, others):
    """Function which is part of the process to make the UpSetPlot,
    counting and retourning the similarities between the clusters."""
    
    cluster_set = set.intersection(*args)
    for i in others:
        cluster_set = cluster_set.difference(set(data[i]))
    return cluster_set, len(cluster_set)


def get_clusters(liste):
    """Function to create every individual cluster depending on 
    the number of organisms given to the UpSetPlot function."""
    
    res = []
    final_res = []
    for i in range(len(liste)-1):
        if i == 0:
            for x in liste:
                z = liste.index(x)
                for i in range(len(liste) - z - 1):
                    res.append([x, liste[z + i + 1]])
            [final_res.append(i) for i in res]
        else:
            res = get_sub_clusters(liste, res)
            [final_res.append(i) for i in res]
    return final_res


def get_sub_clusters(liste, res):
    """Subfunction of the clusters (algorithmic architecture)."""
    
    sub_res = []
    for y in res:
        z = liste.index(y[len(y)-1])
        for i in range(z + 1, len(liste)):
            x = copy.deepcopy(y)
            x.append(liste[i])
            sub_res.append(x)
    return sub_res


def data_venn(WD, model, name):
    """Function to gather the reactions in a model and write them into a CSV file.
    Used to make VENN graph"""
    
    list_id = []
    for reac in model.reactions:
        list_id.append([reac.id])
    utils.write_csv(WD, list_id, name + "_id_reac")


def help_treshold(WD, ini):
    """Function to help choosing the treshold value,
    personal use only, to combine with the R script 'tresholdSearch.R'"""
    
    listValues = []
    for i in range(101):
        test = 10 * i
        ###Putting the test values instead of default values
        draft = blasting.main(ini, 0, 100, 1, 0, test)
        listValues.append([test, len(draft.genes), len(draft.reactions)])
    listValues.insert(0,["Bit_Score", "Nb genes", "Nb reactions"])
    utils.write_csv(WD, listValues, "Treshold")
    

def get_reactions_PT(path):
    """Function to get the reactions in a reactions.dat file of Pathway Tools PGDB."""
    
    liste_Reac = []
    PTtomatoReac = open(path + "/1.0/data/reactions.dat","r")
    for line in PTtomatoReac:
        if "UNIQUE-ID" in line:
            try:
                liste_Reac.append(re.search('(?<=UNIQUE-ID - )\w+(.*\w+)*(-*\w+)*', line).group(0))
            except AttributeError:
                pass
    return liste_Reac


if __name__=="__main__":
    ###Files and working directory###
    WD = '/home/asa/INRAE/Work/Gap_filling/'
    #AraGem models
    # WDtom = '/home/asa/INRAE/Work/blasting_drafts/Tomato_Arabidopsis/'
    # WDkiw = '/home/asa/INRAE/Work/blasting_drafts/Kiwi_Arabidopsis/'
    # WDcuc = '/home/asa/INRAE/Work/blasting_drafts/Cucumber_Arabidopsis/'
    # WDche = '/home/asa/INRAE/Work/blasting_drafts/Cherry_Arabidopsis/'
    # WDcam = '/home/asa/INRAE/Work/blasting_drafts/Camelina_Arabidopsis/'
    #AraCyc models
    # WDtomCyc = '/home/asa/INRAE/Work/blasting_drafts/Tomato_Aracyc/'
    # WDkiwCyc = '/home/asa/INRAE/Work/blasting_drafts/Kiwi_Aracyc/'
    # WDcucCyc = '/home/asa/INRAE/Work/blasting_drafts/Cucumber_Aracyc/'
    # WDcheCyc = '/home/asa/INRAE/Work/blasting_drafts/Cherry_Aracyc/'
    # WDcamCyc = '/home/asa/INRAE/Work/blasting_drafts/Camelina_Aracyc/'
    
    ###The new models (AraGem)
    # tomatoDraft = cobra.io.read_sbml_model(WDtom + "Tomato.xml")
    # kiwiDraft = cobra.io.read_sbml_model(WDkiw + "Kiwi.xml")
    # cucumberDraft = cobra.io.read_sbml_model(WDcuc + "Cucumber.xml")
    # cherryDraft = cobra.io.read_sbml_model(WDche + "Cherry.xml")
    # camelinaDraft = cobra.io.read_sbml_model(WDcam + "Camelina.xml")
    
    ###The new models (AraCyc)
    # tomatoDraftCyc = cobra.io.load_json_model(WDtomCyc + "Tomato.json")
    # kiwiDraftCyc = cobra.io.load_json_model(WDkiwCyc + "Kiwi.json")
    # cucumberDraftCyc = cobra.io.load_json_model(WDcucCyc + "Cucumber.json")
    # cherryDraftCyc = cobra.io.load_json_model(WDcheCyc + "Cherry.json")
    # camelinaDraftCyc = cobra.io.load_json_model(WDcamCyc + "Camelina.json")

    ###Making plots###
    #Loading the data
    # tom = blasting.load_obj(WDtom + "resBlastp")
    # kiw = blasting.load_obj(WDkiw + "resBlastp")
    # cuc = blasting.load_obj(WDcuc + "resBlastp")
    # che = blasting.load_obj(WDche + "resBlastp")
    # cam = blasting.load_obj(WDcam + "resBlastp")
    
    ###Using UpSetplot
    ##AraGem
    # tomatoGenes = read_file(WDtom, "Tomato_id_reac.csv")
    # kiwiGenes = read_file(WDkiw, "Kiwi_id_reac.csv")
    # cucumberGenes = read_file(WDcuc, "Cucumber_id_reac.csv")
    # cherryGenes = read_file(WDche, "Cherry_id_reac.csv")
    # camelinaGenes = read_file(WDcam, "Camelina_id_reac.csv")
    ##AraCyc
    # tomatoGenesCyc = read_file(WDtomCyc, "TomatoCyc_id_reac.csv")
    # kiwiGenesCyc = read_file(WDkiwCyc, "KiwiCyc_id_reac.csv")
    # cucumberGenesCyc = read_file(WDcucCyc, "CucumberCyc_id_reac.csv")
    # cherryGenesCyc = read_file(WDcheCyc, "CherryCyc_id_reac.csv")
    # camelinaGenesCyc = read_file(WDcamCyc, "CamelinaCyc_id_reac.csv")
    
    #Getting the AraCyc reconstruction's reactions VS Pathway Tools:
    # print("Nombre réaction Aracyc perso : ", len(tomatoGenesCyc))
    
    # listeReacPTtomato = get_reactions_PT("/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/sollyphfalsecyc/")
    # print("Nombre réactions PT Tomate : ", len(listeReacPTtomato))
    # listeReacPTkiwi = get_reactions_PT("/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/actchphfalsecyc/")
    # print("Nombre réactions PT Kiwi : ", len(listeReacPTkiwi))
    # listeReacPTcucumber = get_reactions_PT("/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/cucsaphfalsecyc/")
    # print("Nombre réactions PT Concombre : ", len(listeReacPTcucumber))
    # listeReacPTcherry = get_reactions_PT("/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/pruavphfalsecyc/")
    # print("Nombre réactions PT Cerise : ", len(listeReacPTcherry))
    # listeReacPTcamelina = get_reactions_PT("/home/asa/INRAE/Logiciels/ptools-local/pgdbs/user/camsaphfalsecyc/")
    # print("Nombre réactions PT Cameline : ", len(listeReacPTcamelina))

    ###UpsetPlot on PT data
    # dicoUpset = {"Tomato" : listeReacPTtomato,
    #              "Kiwi" : listeReacPTkiwi,
    #              "Cucumber" : listeReacPTcucumber,
    #              "Cherry" : listeReacPTcherry,
    #              "Camelina" : listeReacPTcamelina}
    # make_upsetplot(WD, dicoUpset, "UpsetPlot_PT")

    ###Supervenn on PT data
    # sets_list = [set(listeReacPTtomato), set(listeReacPTkiwi), set(listeReacPTcucumber),
    #              set(listeReacPTcherry), set(listeReacPTcamelina)]
    # species_names = ["Tomato (Pathway Tools)", "Kiwi (Pathway Tools)", "Cucumber (Pathway Tools)",
    #                  "Cherry (Pathway Tools)", "Camelina (Pathway Tools)"]
    # plt.show(supervenn(sets_list, species_names, figsize=(20, 10), rotate_col_annotations=True,
    #       col_annotations_area_height=1.2, sets_ordering='minimize gaps',
    #       min_width_for_annotation=10))
    
    
    ###Supervenn Aracyc VS Pathway Tools Tomato
    # sets_list = [set(tomatoGenesCyc), set(listeReacPTtomato)]
    # species_names = ["Tomato AraCyc", "Tomato Pathway Tools"]
    # plt.show(supervenn(sets_list, species_names, figsize=(20, 10), rotate_col_annotations=True,
    #       col_annotations_area_height=1.2, sets_ordering='minimize gaps',
    #       min_width_for_annotation=10))
    
    ###Creating the dictionary to fit the function
    # dicoUpset = {"Tomato" : tomatoGenes,
    #              "Kiwi" : kiwiGenes,
    #              "Cucumber" : cucumberGenes,
    #              "Cherry" : cherryGenes,
    #              "Camelina" : camelinaGenes}
    # make_upsetplot(WD, dicoUpset, "upsetPlotAraGem")

    # dicoUpset = {"Tomato" : tomatoGenesCyc,
    #              "Kiwi" : kiwiGenesCyc,
    #              "Cucumber" : cucumberGenesCyc,
    #              "Cherry" : cherryGenesCyc,
    #              "Camelina" : camelinaGenesCyc}
    # make_upsetplot(WD, dicoUpset, "upsetPlotAraCyc")
    
    #Tried the supervenn plot, not convinced for this use
    # sets_list = [set(tomatoGenesCyc), set(kiwiGenesCyc), set(cucumberGenesCyc), set(cherryGenesCyc), set(camelinaGenesCyc)]
    # species_names = ["Tomato", "Kiwi", "Cucumber", "Cherry", "Camelina"]
    # plt.show(supervenn(sets_list, species_names, figsize=(20, 10), rotate_col_annotations=True,
    #       col_annotations_area_height=1.2, sets_ordering='minimize gaps',
    #       min_width_for_annotation=10))

    ###Get the reactions id for the Venn diagramm###
    ##First models (AraGem)
    # data_venn(WDtom, tomatoDraft, "Tomato")
    # data_venn(WDkiw, kiwiDraft, "Kiwi")
    # data_venn(WDcuc, cucumberDraft, "Cucumber")
    # data_venn(WDche, cherryDraft, "Cherry")
    # data_venn(WDcam, camelinaDraft, "Camelina")
    ##Second models (AraCyc)
    # data_venn(WDtomCyc, tomatoDraftCyc, "TomatoCyc")
    # data_venn(WDkiwCyc, kiwiDraftCyc, "KiwiCyc")
    # data_venn(WDcucCyc, cucumberDraftCyc, "CucumberCyc")
    # data_venn(WDcheCyc, cherryDraftCyc, "CherryCyc")
    # data_venn(WDcamCyc, camelinaDraftCyc, "CamelinaCyc")

    ###Help to choose the treshold###
    # help_treshold(WDtom, TomatoAracyc.ini)
    # help_treshold(WDkiw, KiwiAracyc.ini)
    # help_treshold(WDcuc, CucumberAracyc.ini)
    # help_treshold(WDche, CherryAracyc.ini)
    # help_treshold(WDcam, CamelinaAracyc.ini)

    #------Identity------#
    # tomato = graph_identity(tom)
    # kiwi = graph_identity(kiw)
    # cucumber = graph_identity(cuc)
    # cherry = graph_identity(che)
    # plt.plot(tomato[0], tomato[1], 'r', label = 'Tomato')
    # plt.plot(kiwi[0], kiwi[1], 'g', label = 'Kiwi')
    # plt.plot(cucumber[0], cucumber[1], 'b', label = 'Cucumber')
    # plt.plot(cherry[0], cherry[1], 'y', label = 'Cherry' )
    # plt.ylabel('Number of mappings')
    # plt.xlabel('Percentage of identity')
    # plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
    # plt.savefig(WD + '4plots.png', dpi=300)
    # plt.show()
    
    #------E-Value------#
    # tomato_e_value = get_e_value(tom)
    # tomato_e_value.insert(0, "Tomato")
    # kiwi_e_value = get_e_value(kiw)
    # kiwi_e_value.insert(0, "Kiwi")
    # cucumber_e_value = get_e_value(cuc)
    # cucumber_e_value.insert(0, "Cucumber")
    # cherry_e_value = get_e_value(che)
    # cherry_e_value.insert(0, "Cherry")
    # total_list = [tomato_e_value, kiwi_e_value, cucumber_e_value, cherry_e_value]
    # write_csv(WD, total_list)
    
    #------Score------#
    # tomato = get_score(tom)
    # print("Mean : %f, max : %f, min : %f" % (mean(tomato), max(tomato), min(tomato)))
    
    #------Bit-Score------#
    # tomato = get_bit_score(tom)
    # print("Mean : %f, max : %f, min : %f" % (mean(tomato), max(tomato), min(tomato)))
    
    #------All Scores------#
    # tomato = get_all_scores("Tomato", tom)
    # tomato.insert(0,["Organism", "Gene", "Identity", "Score", "E_Value", "Bit_Score"])
    # write_csv(WDtom, tomato, "tomato_values")
    
    # kiwi = get_all_scores("Kiwi", kiw)
    # kiwi.insert(0,["Organism", "Gene", "Identity", "Score", "E_Value", "Bit_Score"])
    # write_csv(WDkiw, kiwi, "kiwi_values")
    
    # cucumber = get_all_scores("Cucumber", cuc)
    # cucumber.insert(0,["Organism", "Gene", "Identity", "Score", "E_Value", "Bit_Score"])
    # write_csv(WDcuc, cucumber, "cucumber_values")
    
    # cherry = get_all_scores("Cherry", che)
    # cherry.insert(0,["Organism", "Gene", "Identity", "Score", "E_Value", "Bit_Score"])
    # write_csv(WDche, cherry, "cherry_values")
    
    # camelina = get_all_scores("Camelina", cam)
    # camelina.insert(0,["Organism", "Gene", "Identity", "Score", "E_Value", "Bit_Score"])
    # write_csv(WDcam, camelina, "camelina_values")
    
    #---Calculs savants---#    
    # x = 10
    # mean_score = 0
    # mean_e_value = 0
    # mean_bit_score = 0
    # score = 0
    # e_value = 0
    # bit_score = 0
    # count = 0
    # for i in total:
    #     mean_score += i[2]
    #     mean_e_value += i[3]
    #     mean_bit_score += i[4]
    #     if i[1] > x:
    #         score += i[2]
    #         e_value += i[3]
    #         bit_score += i[4]
    #         count += 1
    # mean_score /= len(total)
    # mean_e_value /= len(total)
    # mean_bit_score /= len(total)
    # print("- Means with identity > %i :\nScore : %f\nE-Value : %f\nBit-Score : %f" % (x, score/count, e_value/count, bit_score/count))
    # print("\n- Total means :\nScore : %f\nE-Value : %f\nBit-Score : %f" % (mean_score, mean_e_value, mean_bit_score))
    
    resTomato = {'Draft network file': '/home/asa/INRAE/Work/Gap_filling/clean_TomatoFusion.sbml', 'Seeds file': '/home/asa/INRAE/Work/Gap_filling/seedsPlants.sbml', 'Targets file': '/home/asa/INRAE/Work/Gap_filling/targetsTomato.sbml', 'Unproducible targets': ['MET_c', 'CYS_c', 'HIS_c', 'TRP_c', 'PRO_c', 'ILE_c', 'GLT_c', 'CIT_c', 'LYS_c', 'ARG_c', 'VAL_c', 'L__45__SELENOCYSTEINE_c', 'LEU_c', 'GLN_c'], 'Repair db file': '/home/asa/INRAE/Work/Gap_filling/metacyc.sbml', 'Unreconstructable targets': ['L__45__SELENOCYSTEINE_c'], 'Reconstructable targets': ['MET_c', 'CYS_c', 'HIS_c', 'TRP_c', 'PRO_c', 'ILE_c', 'GLT_c', 'CIT_c', 'LYS_c', 'ARG_c', 'LEU_c', 'VAL_c', 'GLN_c'], 'Essential reactions': {'MET_c': ['FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN'], 'CYS_c': ['FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN'], 'HIS_c': [], 'TRP_c': [], 'PRO_c': [], 'ILE_c': [], 'GLT_c': [], 'CIT_c': [], 'LYS_c': [], 'ARG_c': [], 'LEU_c': ['FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN'], 'VAL_c': [], 'GLN_c': []}, 'One minimal completion': ['FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN', 'RXN__45__1744'], 'Intersection of cardinality minimal completions': ['FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN'], 'Union of cardinality minimal completions': ['RXN__45__9929', 'SULFITE__45__REDUCT__45__RXN', 'RXN__45__1744', 'FUMARATE__45__REDUCTASE__45__NADH__45__RXN', 'FESO3OXI__45__RXN', 'RXN__45__9772']}
    resKiwi = {'Draft network file': '/home/asa/INRAE/Work/Gap_filling/clean_KiwiFusion.sbml', 'Seeds file': '/home/asa/INRAE/Work/Gap_filling/seedsPlants.sbml', 'Targets file': '/home/asa/INRAE/Work/Gap_filling/targetsTomato.sbml', 'Unproducible targets': ['CIT_c', 'ASN_c', 'L__45__SELENOCYSTEINE_c', 'TRP_c', 'LEU_c', 'HIS_c', 'PRO_c', 'GLY_c', 'MET_c', 'LYS_c', 'L__45__ASPARTATE_c', 'GLN_c', 'GLT_c', 'TYR_c', 'ARG_c', 'ILE_c', 'CYS_c', 'THR_c', 'VAL_c'], 'Repair db file': '/home/asa/INRAE/Work/Gap_filling/metacyc.sbml', 'Unreconstructable targets': ['L__45__SELENOCYSTEINE_c'], 'Reconstructable targets': ['MET_c', 'HIS_c', 'VAL_c', 'PRO_c', 'LYS_c', 'L__45__ASPARTATE_c', 'GLN_c', 'CIT_c', 'CYS_c', 'ASN_c', 'GLT_c', 'THR_c', 'TYR_c', 'GLY_c', 'TRP_c', 'ARG_c', 'ILE_c', 'LEU_c'], 'Essential reactions': {'MET_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'HIS_c': [], 'VAL_c': [], 'PRO_c': [], 'LYS_c': [], 'L__45__ASPARTATE_c': [], 'GLN_c': [], 'CIT_c': [], 'CYS_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'ASN_c': [], 'GLT_c': [], 'THR_c': [], 'TYR_c': [], 'GLY_c': [], 'TRP_c': [], 'ARG_c': [], 'ILE_c': [], 'LEU_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN']}, 'One minimal completion': ['RXN__45__8633', 'RXN__45__18377', 'KETOPANTOALDOLASE__45__RXN', '_2__45__DEHYDROPANTOATE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN'], 'Intersection of cardinality minimal completions': ['SULFITE__45__REDUCT__45__RXN', 'RXN__45__18377', 'FESO3OXI__45__RXN'], 'Union of cardinality minimal completions': ['DIMETHYLMALATE__45__DEHYDROGENASE__45__RXN', 'PANTOATE__45__4__45__DEHYDROGENASE__45__RXN', 'RXN__45__18377', 'RXN__45__8633', 'KETOPANTOALDOLASE__45__RXN', '_2__45__DEHYDROPANTOATE__45__REDUCT__45__RXN', 'RXN__45__10065', 'R__45__DEHYDROPANTOATE__45__DEHYDROGENASE__45__RXN', 'FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN']}
    resCucumber = {'Draft network file': '/home/asa/INRAE/Work/Gap_filling/clean_CucumberFusion.sbml', 'Seeds file': '/home/asa/INRAE/Work/Gap_filling/seedsPlants.sbml', 'Targets file': '/home/asa/INRAE/Work/Gap_filling/targetsTomato.sbml', 'Unproducible targets': ['CIT_c', 'GLT_c', 'L__45__SELENOCYSTEINE_c', 'HIS_c', 'MET_c', 'ARG_c', 'LEU_c', 'TRP_c', 'PRO_c', 'CYS_c', 'VAL_c', 'GLN_c', 'ILE_c', 'LYS_c'], 'Repair db file': '/home/asa/INRAE/Work/Gap_filling/metacyc.sbml', 'Unreconstructable targets': ['L__45__SELENOCYSTEINE_c'], 'Reconstructable targets': ['CIT_c', 'GLT_c', 'HIS_c', 'ARG_c', 'MET_c', 'LEU_c', 'TRP_c', 'PRO_c', 'CYS_c', 'VAL_c', 'GLN_c', 'ILE_c', 'LYS_c'], 'Essential reactions': {'CIT_c': [], 'GLT_c': [], 'HIS_c': [], 'ARG_c': [], 'MET_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'LEU_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN', 'DIHYDROXYISOVALDEHYDRAT__45__RXN'], 'TRP_c': [], 'PRO_c': [], 'CYS_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'VAL_c': ['DIHYDROXYISOVALDEHYDRAT__45__RXN'], 'GLN_c': [], 'ILE_c': ['DIHYDROXYMETVALDEHYDRAT__45__RXN'], 'LYS_c': []}, 'One minimal completion': ['SULFITE__45__REDUCT__45__RXN', 'DIHYDROXYISOVALDEHYDRAT__45__RXN', 'RXN__45__1744', 'FESO3OXI__45__RXN', 'DIHYDROXYMETVALDEHYDRAT__45__RXN', 'RXN__45__18377'], 'Intersection of cardinality minimal completions': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN', 'DIHYDROXYMETVALDEHYDRAT__45__RXN', 'DIHYDROXYISOVALDEHYDRAT__45__RXN'], 'Union of cardinality minimal completions': ['RXN__45__9929', 'SULFITE__45__REDUCT__45__RXN', 'DIHYDROXYISOVALDEHYDRAT__45__RXN', 'RXN__45__9772', 'FUMARATE__45__REDUCTASE__45__NADH__45__RXN', 'RXN__45__1744', 'FESO3OXI__45__RXN', 'DIHYDROXYMETVALDEHYDRAT__45__RXN', 'DIHYDRONEOPTERIN__45__MONO__45__P__45__DEPHOS__45__RXN', 'RXN__45__18377', 'RXN__45__18378']}
    resCherry = {'Draft network file': '/home/asa/INRAE/Work/Gap_filling/clean_CherryFusion.sbml', 'Seeds file': '/home/asa/INRAE/Work/Gap_filling/seedsPlants.sbml', 'Targets file': '/home/asa/INRAE/Work/Gap_filling/targetsTomato.sbml', 'Unproducible targets': ['GLT_c', 'THR_c', 'MET_c', 'L__45__SELENOCYSTEINE_c', 'TYR_c', 'L__45__ASPARTATE_c', 'PRO_c', 'HIS_c', 'GLN_c', 'LEU_c', 'VAL_c', 'GLY_c', 'TRP_c', 'ILE_c', 'ARG_c', 'CYS_c', 'ASN_c', 'CIT_c', 'LYS_c'], 'Repair db file': '/home/asa/INRAE/Work/Gap_filling/metacyc.sbml', 'Unreconstructable targets': ['L__45__SELENOCYSTEINE_c'], 'Reconstructable targets': ['TYR_c', 'ARG_c', 'LEU_c', 'L__45__ASPARTATE_c', 'VAL_c', 'GLY_c', 'GLT_c', 'CYS_c', 'PRO_c', 'ASN_c', 'CIT_c', 'TRP_c', 'MET_c', 'THR_c', 'HIS_c', 'LYS_c', 'GLN_c', 'ILE_c'], 'Essential reactions': {'TYR_c': [], 'ARG_c': [], 'LEU_c': ['FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN'], 'L__45__ASPARTATE_c': [], 'VAL_c': [], 'GLY_c': [], 'GLT_c': [], 'CYS_c': ['FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN'], 'PRO_c': [], 'ASN_c': [], 'CIT_c': [], 'TRP_c': [], 'MET_c': ['FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN'], 'THR_c': [], 'HIS_c': [], 'LYS_c': [], 'GLN_c': [], 'ILE_c': []}, 'One minimal completion': ['_2__45__DEHYDROPANTOATE__45__REDUCT__45__RXN', 'RXN__45__9772', '_1__46__4__46__1__46__21__45__RXN__45__L__45__ASPARTATE__47__NADP__47__WATER__47____47__OXALACETIC_ACID__47__AMMONIUM__47__NADPH__47__PROTON__46__62__46__', 'FESO3OXI__45__RXN', 'RXN__45__18377', 'SULFITE__45__REDUCT__45__RXN'], 'Intersection of cardinality minimal completions': ['FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN'], 'Union of cardinality minimal completions': ['GLYCOLALD__45__DEHYDROG__45__RXN', '_1__46__4__46__1__46__21__45__RXN__45__L__45__ASPARTATE__47__NADP__47__WATER__47____47__OXALACETIC_ACID__47__AMMONIUM__47__NADPH__47__PROTON__46__62__46__', 'RXN__45__18377', 'PANTOATE__45__4__45__DEHYDROGENASE__45__RXN', 'ASPARTASE__45__RXN', '_2__45__DEHYDROPANTOATE__45__REDUCT__45__RXN', 'RXN__45__9772', 'FUMARATE__45__REDUCTASE__45__NADH__45__RXN', 'DIMETHYLMALATE__45__DEHYDROGENASE__45__RXN', 'PYROXALTRANSAM__45__RXN', 'SULFITE__45__REDUCT__45__RXN', 'R__45__DEHYDROPANTOATE__45__DEHYDROGENASE__45__RXN', 'RXN__45__1744', 'RXN__45__9929', 'RXN__45__18378', 'KETOPANTOALDOLASE__45__RXN', 'RXN__45__10065', '_2__46__6__46__1__46__70__45__RXN', 'FESO3OXI__45__RXN', 'DIHYDRONEOPTERIN__45__MONO__45__P__45__DEPHOS__45__RXN', 'RXN__45__8633', '_1__46__4__46__1__46__21__45__RXN__45__L__45__ASPARTATE__47__NAD__47__WATER__47____47__OXALACETIC_ACID__47__AMMONIUM__47__NADH__47__PROTON__46__60__46__']}
    resCamelina = {'Draft network file': '/home/asa/INRAE/Work/Gap_filling/clean_CamelinaFusion.sbml', 'Seeds file': '/home/asa/INRAE/Work/Gap_filling/seedsPlants.sbml', 'Targets file': '/home/asa/INRAE/Work/Gap_filling/targetsTomato.sbml', 'Unproducible targets': ['ASN_c', 'TRP_c', 'PRO_c', 'TYR_c', 'GLN_c', 'CIT_c', 'L__45__ASPARTATE_c', 'MET_c', 'GLY_c', 'LEU_c', 'CYS_c', 'LYS_c', 'ARG_c', 'HIS_c', 'L__45__SELENOCYSTEINE_c', 'THR_c', 'VAL_c', 'GLT_c', 'ILE_c'], 'Repair db file': '/home/asa/INRAE/Work/Gap_filling/metacyc.sbml', 'Unreconstructable targets': ['L__45__SELENOCYSTEINE_c'], 'Reconstructable targets': ['CYS_c', 'THR_c', 'L__45__ASPARTATE_c', 'VAL_c', 'ILE_c', 'ASN_c', 'LYS_c', 'TRP_c', 'PRO_c', 'ARG_c', 'HIS_c', 'TYR_c', 'GLN_c', 'MET_c', 'GLT_c', 'GLY_c', 'LEU_c', 'CIT_c'], 'Essential reactions': {'CYS_c': ['FESO3OXI__45__RXN'], 'THR_c': [], 'L__45__ASPARTATE_c': [], 'VAL_c': [], 'ILE_c': [], 'ASN_c': [], 'LYS_c': [], 'TRP_c': [], 'PRO_c': [], 'ARG_c': [], 'HIS_c': [], 'TYR_c': [], 'GLN_c': [], 'MET_c': ['FESO3OXI__45__RXN'], 'GLT_c': [], 'GLY_c': [], 'LEU_c': ['FESO3OXI__45__RXN'], 'CIT_c': []}, 'One minimal completion': ['DIHYDRONEOPTERIN__45__MONO__45__P__45__DEPHOS__45__RXN', 'RXN__45__11108', 'FESO3OXI__45__RXN', 'FUMARATE__45__REDUCTASE__45__NADH__45__RXN'], 'Intersection of cardinality minimal completions': ['DIHYDRONEOPTERIN__45__MONO__45__P__45__DEPHOS__45__RXN', 'FESO3OXI__45__RXN'], 'Union of cardinality minimal completions': ['FUMARATE__45__REDUCTASE__45__NADH__45__RXN', 'DIHYDRONEOPTERIN__45__MONO__45__P__45__DEPHOS__45__RXN', 'RXN__45__15583', 'RXN__45__1744', 'FESO3OXI__45__RXN', 'RXN__45__11108', 'SULFITE__45__REDUCT__45__RXN']}

    resTomato['Union of cardinality minimal completions']
    dicoUpset = {
        "Tomato" : resTomato['Union of cardinality minimal completions'],
        "Kiwi" : resKiwi['Union of cardinality minimal completions'],
        "Cucumber" : resCucumber['Union of cardinality minimal completions'],
        "Cherry" : resCherry['Union of cardinality minimal completions'],
        "Camelina" : resCamelina['Union of cardinality minimal completions']}
    make_upsetplot(WD, dicoUpset, "UpsetPlot_Gap-Filling")