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
    for c in clusters:
        others = list(data.keys())
        listInter = []
        for x in c:
            others.remove(x)
            listInter.append(set(data[x]))
        count.append(similiraty_count(data, listInter, others))
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
    return len(cluster_set)


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
    WD = '/home/antoine/INRAE/Work/'
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
    
    # resTomato = {'Draft network file': '/home/antoine/INRAE/Work/Gap_filling/clean_TomatoFusion.sbml', 'Seeds file': '/home/antoine/INRAE/Work/Gap_filling/seedsPlants.sbml', 'Targets file': '/home/antoine/INRAE/Work/Gap_filling/targetsTomato.sbml', 'Unproducible targets': ['LEU_c', 'HIS_c', 'MET_c', 'L__45__SELENOCYSTEINE_c', 'ILE_c', 'TRP_c', 'CYS_c', 'ARG_c', 'GLN_c', 'PRO_c', 'VAL_c', 'CPD__45__10726_c', 'LYS_c', 'CIT_c', 'GLT_c'], 'Repair db file': '/home/antoine/INRAE/Work/Gap_filling/metacyc.sbml', 'Unreconstructable targets': ['CPD__45__10726_c', 'L__45__SELENOCYSTEINE_c'], 'Reconstructable targets': ['LEU_c', 'HIS_c', 'MET_c', 'ILE_c', 'TRP_c', 'CYS_c', 'GLN_c', 'PRO_c', 'VAL_c', 'LYS_c', 'ARG_c', 'CIT_c', 'GLT_c'], 'Essential reactions': {'LEU_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'HIS_c': [], 'MET_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'ILE_c': [], 'TRP_c': [], 'CYS_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'GLN_c': [], 'PRO_c': [], 'VAL_c': [], 'LYS_c': [], 'ARG_c': [], 'CIT_c': [], 'GLT_c': []}, 'One minimal completion': ['FUMARATE__45__REDUCTASE__45__NADH__45__RXN', 'SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'Intersection of cardinality minimal completions': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'Union of cardinality minimal completions': ['FESO3OXI__45__RXN', 'RXN__45__1744', 'FUMARATE__45__REDUCTASE__45__NADH__45__RXN', 'RXN__45__9772', 'RXN__45__9929', 'SULFITE__45__REDUCT__45__RXN']}
    # resKiwi = {'Draft network file': '/home/antoine/INRAE/Work/Gap_filling/clean_KiwiFusion.sbml', 'Seeds file': '/home/antoine/INRAE/Work/Gap_filling/seedsPlants.sbml', 'Targets file': '/home/antoine/INRAE/Work/Gap_filling/targetsKiwi.sbml', 'Unproducible targets': ['ARACHIDIC_ACID_c', 'L__45__ASPARTATE_c', 'LINOLEIC_ACID_c', 'ARABINOSE_c', 'GLY_c', 'CPD__45__12521_c', 'CPD__45__10726_c', 'PALMITATE_c', 'LYS_c', 'L__45__SELENOCYSTEINE_c', 'OLEATE__45__CPD_c', 'GLN_c', 'ARG_c', 'TETRACOSANOATE_c', 'THREO__45__DS__45__ISO__45__CITRATE_c', 'LINOLENIC_ACID_c', 'TRP_c', 'ILE_c', 'CPD__45__13293_c', 'GLT_c', 'MET_c', 'CYS_c', 'STEARIC_ACID_c', 'THR_c', 'CPD__45__12524_c', 'LEU_c', 'CIT_c', 'ASN_c', 'VAL_c', 'HIS_c', 'ASCORBATE_c', 'RHAMNOSE_c', 'TYR_c', 'XYLOSE_c', 'PRO_c'], 'Repair db file': '/home/antoine/INRAE/Work/Gap_filling/metacyc.sbml', 'Unreconstructable targets': ['ARACHIDIC_ACID_c', 'TETRACOSANOATE_c', 'STEARIC_ACID_c', 'RHAMNOSE_c', 'LINOLENIC_ACID_c', 'LINOLEIC_ACID_c', 'CPD__45__12524_c', 'CPD__45__12521_c', 'CPD__45__10726_c', 'PALMITATE_c', 'CPD__45__13293_c', 'L__45__SELENOCYSTEINE_c', 'OLEATE__45__CPD_c'], 'Reconstructable targets': ['L__45__ASPARTATE_c', 'ARABINOSE_c', 'GLY_c', 'LYS_c', 'GLN_c', 'ARG_c', 'THREO__45__DS__45__ISO__45__CITRATE_c', 'TRP_c', 'ILE_c', 'GLT_c', 'MET_c', 'CYS_c', 'THR_c', 'LEU_c', 'CIT_c', 'ASN_c', 'VAL_c', 'HIS_c', 'ASCORBATE_c', 'TYR_c', 'XYLOSE_c', 'PRO_c'], 'Essential reactions': {'L__45__ASPARTATE_c': [], 'ARABINOSE_c': ['L__45__ARABINITOL__45__4__45__DEHYDROGENASE__45__RXN'], 'GLY_c': [], 'LYS_c': [], 'GLN_c': [], 'ARG_c': [], 'THREO__45__DS__45__ISO__45__CITRATE_c': [], 'TRP_c': [], 'ILE_c': [], 'GLT_c': [], 'MET_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'CYS_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'THR_c': [], 'LEU_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'CIT_c': [], 'ASN_c': [], 'VAL_c': [], 'HIS_c': [], 'ASCORBATE_c': [], 'TYR_c': [], 'XYLOSE_c': [], 'PRO_c': []}, 'One minimal completion': ['SULFITE__45__REDUCT__45__RXN', 'L__45__ARABINITOL__45__4__45__DEHYDROGENASE__45__RXN', 'TRANS__45__RXN__45__10__45__L__45__ARABINOSE__47__PROTON__47____47__ARABINOSE__47__PROTON__46__37__46__', 'DEHYDRO__45__L__45__GULONATE__45__DECARBOXYLASE__45__RXN', 'RXN__45__14102', 'ATPSYN__45__RXN__91__CCO__45__PERI__45__BAC__45__CCO__45__CYTOSOL__93____45__ATP__47__WATER__47__PROTON__47____47__ADP__47__Pi__47__PROTON__46__58__46__', '_1__46__3__46__3__46__12__45__RXN', 'FESO3OXI__45__RXN', 'DIMETHYLMALATE__45__DEHYDROGENASE__45__RXN', 'R__45__DEHYDROPANTOATE__45__DEHYDROGENASE__45__RXN', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__L__45__ARABINOSE__47__PROTON__46__37__46__', '_3__45__DEHYDRO__45__L__45__GULONATE__45__2__45__DEHYDROGENASE__45__RXN__45__3__45__KETO__45__L__45__GULONATE__47__NADP__47____47__CPD__45__334__47__NADPH__47__PROTON__46__45__46__', 'PANTOATE__45__4__45__DEHYDROGENASE__45__RXN', 'RXN__45__11825', 'RXN__45__18377'], 'Intersection of cardinality minimal completions': ['SULFITE__45__REDUCT__45__RXN', 'L__45__ARABINITOL__45__4__45__DEHYDROGENASE__45__RXN', 'DEHYDRO__45__L__45__GULONATE__45__DECARBOXYLASE__45__RXN', 'RXN__45__14102', '_1__46__3__46__3__46__12__45__RXN', 'FESO3OXI__45__RXN', 'RXN__45__18377'], 'Union of cardinality minimal completions': ['ATPSYN__45__RXN__91__CCO__45__PM__45__BAC__45__NEG__93____45__ATP__47__WATER__47__PROTON__47____47__ADP__47__Pi__47__PROTON__46__48__46__', 'KETOPANTOALDOLASE__45__RXN', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__CPD__45__12046__47__PROTON__46__35__46__', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__L__45__ARABINOSE__47__PROTON__46__37__46__', 'TRANS__45__RXN__45__262', 'RXN__45__11825', 'SULFITE__45__REDUCT__45__RXN', 'ABC__45__2__45__RXN__45__ATP__47__CPD__45__12045__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__45__46__', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__BETA__45__L__45__ARABINOSE__47__PROTON__46__42__46__', 'ATPSYN__45__RXN__91__CCO__45__PERI__45__BAC__45__CCO__45__CYTOSOL__93____45__ATP__47__WATER__47__PROTON__47____47__ADP__47__Pi__47__PROTON__46__58__46__', 'TRANS__45__RXN__45__10__45__CPD__45__12046__47__PROTON__47____47__ARABINOSE__47__PROTON__46__35__46__', 'TRANS__45__RXN__45__10__45__BETA__45__L__45__ARABINOSE__47__PROTON__47____47__ARABINOSE__47__PROTON__46__42__46__', 'ABC__45__2__45__RXN__45__ATP__47__L__45__ARABINOSE__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__47__46__', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__CPD__45__15699__47__PROTON__46__35__46__', 'TRANS__45__RXN__45__10__45__CPD__45__15699__47__PROTON__47____47__ARABINOSE__47__PROTON__46__35__46__', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__ARABINOSE__47__PROTON__46__35__46__', '_2__45__DEHYDROPANTOATE__45__REDUCT__45__RXN', 'L__45__ARABINITOL__45__4__45__DEHYDROGENASE__45__RXN', 'RXN__45__8633', 'TRANS__45__RXN__45__10__45__CPD__45__12045__47__PROTON__47____47__ARABINOSE__47__PROTON__46__35__46__', 'TRANS__45__RXN0__45__571', 'L__45__XYLULOSE__45__REDUCTASE__45__RXN', '_1__46__3__46__3__46__12__45__RXN', 'ABC__45__2__45__RXN__45__ATP__47__CPD__45__12046__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__45__46__', 'FESO3OXI__45__RXN', 'DIMETHYLMALATE__45__DEHYDROGENASE__45__RXN', 'R__45__DEHYDROPANTOATE__45__DEHYDROGENASE__45__RXN', 'RXN__45__10065', 'TRANS__45__RXN__45__10__45__ARABINOSE__47__PROTON__47____47__ARABINOSE__47__PROTON__46__35__46__', 'TRANS__45__RXN__45__284', '_3__45__DEHYDRO__45__L__45__GULONATE__45__2__45__DEHYDROGENASE__45__RXN__45__3__45__KETO__45__L__45__GULONATE__47__NADP__47____47__CPD__45__334__47__NADPH__47__PROTON__46__45__46__', 'TRANS__45__RXN__45__300', 'RXN__45__18377', 'ABC__45__2__45__RXN__45__ATP__47__BETA__45__L__45__ARABINOSE__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__52__46__', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__CPD__45__12045__47__PROTON__46__35__46__', 'ABC__45__2__45__RXN__45__ATP__47__CPD__45__15699__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__45__46__', 'DEHYDRO__45__L__45__GULONATE__45__DECARBOXYLASE__45__RXN', 'RXN__45__14102', 'RXN0__45__703', 'TRANS__45__RXN__45__10__45__L__45__ARABINOSE__47__PROTON__47____47__ARABINOSE__47__PROTON__46__37__46__', 'PANTOATE__45__4__45__DEHYDROGENASE__45__RXN', 'ABC__45__2__45__RXN__45__ATP__47__ARABINOSE__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__45__46__']}
    # resCucumber = {'Draft network file': '/home/antoine/INRAE/Work/Gap_filling/clean_CucumberFusion.sbml', 'Seeds file': '/home/antoine/INRAE/Work/Gap_filling/seedsPlants.sbml', 'Targets file': '/home/antoine/INRAE/Work/Gap_filling/targetsCucumber.sbml', 'Unproducible targets': ['CPD__45__13293_c', 'STEARIC_ACID_c', 'RHAMNOSE_c', 'ARABINOSE_c', 'TETRACOSANOATE_c', 'CPD__45__12524_c', 'ILE_c', 'HIS_c', 'LINOLEIC_ACID_c', 'OLEATE__45__CPD_c', 'LYS_c', 'PRO_c', 'MET_c', 'ARACHIDIC_ACID_c', 'GLT_c', 'LEU_c', 'XYLOSE_c', 'CPD__45__12521_c', 'VAL_c', 'TRP_c', 'L__45__SELENOCYSTEINE_c', 'GLN_c', 'CYS_c', 'ASCORBATE_c', 'LINOLENIC_ACID_c', 'PALMITATE_c', 'CPD__45__10726_c', 'ARG_c', 'CIT_c'], 'Repair db file': '/home/antoine/INRAE/Work/Gap_filling/metacyc.sbml', 'Unreconstructable targets': ['ARACHIDIC_ACID_c', 'CPD__45__13293_c', 'STEARIC_ACID_c', 'L__45__SELENOCYSTEINE_c', 'RHAMNOSE_c', 'TETRACOSANOATE_c', 'LINOLEIC_ACID_c', 'OLEATE__45__CPD_c', 'LINOLENIC_ACID_c', 'CPD__45__12524_c', 'PALMITATE_c', 'CPD__45__10726_c'], 'Reconstructable targets': ['MET_c', 'CIT_c', 'ARABINOSE_c', 'HIS_c', 'GLT_c', 'GLN_c', 'LYS_c', 'LEU_c', 'XYLOSE_c', 'CYS_c', 'CPD__45__12521_c', 'ASCORBATE_c', 'ILE_c', 'ARG_c', 'VAL_c', 'PRO_c', 'TRP_c'], 'Essential reactions': {'MET_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'CIT_c': [], 'ARABINOSE_c': ['L__45__ARABINITOL__45__4__45__DEHYDROGENASE__45__RXN'], 'HIS_c': [], 'GLT_c': [], 'GLN_c': [], 'LYS_c': [], 'LEU_c': ['DIHYDROXYISOVALDEHYDRAT__45__RXN', 'SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'XYLOSE_c': [], 'CYS_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'CPD__45__12521_c': ['LUTEOLIN__45__7__45__O__45__GLUCORONOSYLTRANSFERASE__45__RXN', 'RXN__45__13942', 'DIHYDROXYISOVALDEHYDRAT__45__RXN', 'FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN'], 'ASCORBATE_c': [], 'ILE_c': ['DIHYDROXYMETVALDEHYDRAT__45__RXN'], 'ARG_c': [], 'VAL_c': ['DIHYDROXYISOVALDEHYDRAT__45__RXN'], 'PRO_c': [], 'TRP_c': []}, 'One minimal completion': ['LUTEOLIN__45__7__45__O__45__GLUCORONOSYLTRANSFERASE__45__RXN', 'DIHYDROXYISOVALDEHYDRAT__45__RXN', 'RXN__45__13942', 'L__45__ARABINITOL__45__4__45__DEHYDROGENASE__45__RXN', 'RXN__45__13935', 'RXN__45__1744', 'FESO3OXI__45__RXN', 'ABC__45__2__45__RXN__45__ATP__47__L__45__ARABINOSE__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__47__46__', 'ATPSYN__45__RXN__91__CCO__45__PM__45__BAC__45__NEG__93____45__ATP__47__WATER__47__PROTON__47____47__ADP__47__Pi__47__PROTON__46__48__46__', 'RXN__45__18377', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__L__45__ARABINOSE__47__PROTON__46__37__46__', '_1__46__3__46__3__46__12__45__RXN', 'RXN__45__16249', 'DIHYDROXYMETVALDEHYDRAT__45__RXN', 'RXN__45__11315__45__XYLITOL__47__NAD__47____47__D__45__Xylose__47__NADH__47__PROTON__46__34__46__', 'SULFITE__45__REDUCT__45__RXN', 'RXN__45__11825', 'RXN__45__14102'], 'Intersection of cardinality minimal completions': ['LUTEOLIN__45__7__45__O__45__GLUCORONOSYLTRANSFERASE__45__RXN', 'RXN__45__13942', 'DIHYDROXYISOVALDEHYDRAT__45__RXN', 'L__45__ARABINITOL__45__4__45__DEHYDROGENASE__45__RXN', 'RXN__45__13935', 'FESO3OXI__45__RXN', '_1__46__3__46__3__46__12__45__RXN', 'DIHYDROXYMETVALDEHYDRAT__45__RXN', 'SULFITE__45__REDUCT__45__RXN', 'RXN__45__14102'], 'Union of cardinality minimal completions': ['LUTEOLIN__45__7__45__O__45__GLUCORONOSYLTRANSFERASE__45__RXN', 'RXN__45__8773', 'TRANS__45__RXN__45__10__45__ARABINOSE__47__PROTON__47____47__ARABINOSE__47__PROTON__46__35__46__', 'TRANS__45__RXN__45__10__45__BETA__45__L__45__ARABINOSE__47__PROTON__47____47__ARABINOSE__47__PROTON__46__42__46__', 'L__45__ARABINITOL__45__4__45__DEHYDROGENASE__45__RXN', 'FESO3OXI__45__RXN', 'ABC__45__2__45__RXN__45__ATP__47__L__45__ARABINOSE__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__47__46__', 'RXN__45__11315', 'FUMARATE__45__REDUCTASE__45__NADH__45__RXN', 'ATPSYN__45__RXN__91__CCO__45__PM__45__BAC__45__NEG__93____45__ATP__47__WATER__47__PROTON__47____47__ADP__47__Pi__47__PROTON__46__48__46__', 'RXN__45__18377', 'RXN__45__11825', 'RXN__45__11759', 'ABC__45__2__45__RXN__45__ATP__47__BETA__45__L__45__ARABINOSE__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__52__46__', 'DIHYDROXYISOVALDEHYDRAT__45__RXN', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__BETA__45__L__45__ARABINOSE__47__PROTON__46__42__46__', 'ABC__45__2__45__RXN__45__ATP__47__ARABINOSE__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__45__46__', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__CPD__45__15699__47__PROTON__46__35__46__', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__CPD__45__12046__47__PROTON__46__35__46__', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__L__45__ARABINOSE__47__PROTON__46__37__46__', 'RXN__45__18378', 'L__45__XYLULOSE__45__REDUCTASE__45__RXN', 'SULFITE__45__REDUCT__45__RXN', 'TRANS__45__RXN0__45__571', 'ABC__45__2__45__RXN__45__ATP__47__CPD__45__12045__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__45__46__', 'ABC__45__2__45__RXN__45__ATP__47__CPD__45__12046__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__45__46__', 'ATPSYN__45__RXN__91__CCO__45__PERI__45__BAC__45__CCO__45__CYTOSOL__93____45__ATP__47__WATER__47__PROTON__47____47__ADP__47__Pi__47__PROTON__46__58__46__', 'TRANS__45__RXN__45__10__45__CPD__45__12045__47__PROTON__47____47__ARABINOSE__47__PROTON__46__35__46__', '_1__46__3__46__3__46__12__45__RXN', 'RXN__45__16249', 'DIHYDRONEOPTERIN__45__MONO__45__P__45__DEPHOS__45__RXN', 'DIHYDROXYMETVALDEHYDRAT__45__RXN', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__CPD__45__12045__47__PROTON__46__35__46__', 'RXN__45__13942', 'TRANS__45__RXN__45__300', 'RXN__45__9772', 'TRANS__45__RXN__45__284', 'RXN__45__9929', 'RXN__45__13935', 'XYLISOM__45__RXN__45__XYLOSE__47____47__D__45__XYLULOSE__46__19__46__', 'RXN__45__1744', 'TRANS__45__RXN__45__262', 'TRANS__45__RXN__45__10__45__CPD__45__12046__47__PROTON__47____47__ARABINOSE__47__PROTON__46__35__46__', 'TRANS__45__RXN__45__10__45__CPD__45__15699__47__PROTON__47____47__ARABINOSE__47__PROTON__46__35__46__', 'D__45__XYLULOSE__45__REDUCTASE__45__RXN', 'RXN__45__11315__45__XYLITOL__47__NAD__47____47__D__45__Xylose__47__NADH__47__PROTON__46__34__46__', 'ABC__45__2__45__RXN__45__ATP__47__CPD__45__15699__47__WATER__47____47__ADP__47__ARABINOSE__47__Pi__47__PROTON__46__45__46__', 'TRANS__45__RXN__45__40__45__CPD__45__15699__47__PROTON__47____47__ARABINOSE__47__PROTON__46__35__46__', 'TRANS__45__RXN__45__10__45__L__45__ARABINOSE__47__PROTON__47____47__ARABINOSE__47__PROTON__46__37__46__', 'RXN__45__14102']}
    # resCherry = {'Draft network file': '/home/antoine/INRAE/Work/Gap_filling/clean_CherryFusion.sbml', 'Seeds file': '/home/antoine/INRAE/Work/Gap_filling/seedsPlants.sbml', 'Targets file': '/home/antoine/INRAE/Work/Gap_filling/targetsCherry.sbml', 'Unproducible targets': ['CPD__45__7409_c', 'L__45__SELENOCYSTEINE_c', 'TYR_c', 'CPD1F__45__129_c', 'GLN_c', 'ASN_c', 'CPD__45__8491_c', 'CPD1F__45__118_c', 'THR_c', 'BENZYL__45__ALCOHOL_c', 'CPD1F__45__130_c', 'ILE_c', 'HEXANAL_c', 'CPD__45__10726_c', 'CPD__45__8490_c', 'TRP_c', 'HIS_c', 'VAL_c', 'LEU_c', 'PRO_c', 'CPD__45__12321_c', 'CYS_c', 'TRANS__45__2__45__HEXENAL_c', 'TRANS__45__2__45__HEXENOL_c', 'CPD__45__520_c', 'N__45__ACETYL__45__5__45__METHOXY__45__TRYPTAMINE_c', 'MET_c', 'BENZALDEHYDE_c', 'GLY_c', 'L__45__ASPARTATE_c', 'GLT_c', 'CPD__45__19100_c', 'CPD1F__45__119_c', 'LYS_c', 'CIS__45__3__45__HEXENAL_c', 'HEXANOL__45__CMPD_c', 'ARG_c'], 'Repair db file': '/home/antoine/INRAE/Work/Gap_filling/metacyc.sbml', 'Unreconstructable targets': ['CPD__45__7409_c', 'CPD__45__520_c', 'L__45__SELENOCYSTEINE_c', 'CPD1F__45__130_c', 'CPD__45__19100_c', 'CPD__45__10726_c', 'CPD1F__45__119_c', 'N__45__ACETYL__45__5__45__METHOXY__45__TRYPTAMINE_c', 'CIS__45__3__45__HEXENAL_c', 'CPD__45__8490_c', 'CPD__45__8491_c'], 'Reconstructable targets': ['TYR_c', 'CPD1F__45__129_c', 'GLN_c', 'ASN_c', 'CPD1F__45__118_c', 'THR_c', 'BENZYL__45__ALCOHOL_c', 'ILE_c', 'HEXANAL_c', 'VAL_c', 'HIS_c', 'TRP_c', 'LEU_c', 'PRO_c', 'CPD__45__12321_c', 'CYS_c', 'TRANS__45__2__45__HEXENAL_c', 'TRANS__45__2__45__HEXENOL_c', 'MET_c', 'BENZALDEHYDE_c', 'GLY_c', 'L__45__ASPARTATE_c', 'GLT_c', 'LYS_c', 'HEXANOL__45__CMPD_c', 'ARG_c'], 'Essential reactions': {'TYR_c': [], 'CPD1F__45__129_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'GLN_c': [], 'ASN_c': [], 'CPD1F__45__118_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN', 'RXN1F__45__148'], 'THR_c': [], 'BENZYL__45__ALCOHOL_c': ['BENZYL__45__ALC__45__DEHYDROGENASE__45__RXN'], 'ILE_c': [], 'HEXANAL_c': ['SULFITE__45__REDUCT__45__RXN', 'RXN__45__12565', 'FESO3OXI__45__RXN'], 'VAL_c': [], 'HIS_c': [], 'TRP_c': [], 'LEU_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'PRO_c': [], 'CPD__45__12321_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'CYS_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'TRANS__45__2__45__HEXENAL_c': ['SULFITE__45__REDUCT__45__RXN', 'RXN__45__12565', 'FESO3OXI__45__RXN'], 'TRANS__45__2__45__HEXENOL_c': ['SULFITE__45__REDUCT__45__RXN', 'RXN__45__12565', 'FESO3OXI__45__RXN'], 'MET_c': ['SULFITE__45__REDUCT__45__RXN', 'FESO3OXI__45__RXN'], 'BENZALDEHYDE_c': [], 'GLY_c': [], 'L__45__ASPARTATE_c': [], 'GLT_c': [], 'LYS_c': [], 'HEXANOL__45__CMPD_c': ['SULFITE__45__REDUCT__45__RXN', 'RXN__45__12565', 'HEXANOL__45__RXN', 'FESO3OXI__45__RXN'], 'ARG_c': []}, 'One minimal completion': ['FUMARATE__45__REDUCTASE__45__NADH__45__RXN', 'RXN__45__12559', 'ISOBUTYRYL__45__COA__45__MUTASE__45__RXN', 'RXN1F__45__148', 'RXN__45__9336', 'BENZYL__45__ALC__45__DEHYDROGENASE__45__RXN', 'RXN__45__19913', 'RXN__45__12570', 'HEXANOL__45__RXN', 'FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN', 'DIHYDRONEOPTERIN__45__MONO__45__P__45__DEPHOS__45__RXN', 'GLYCOLALD__45__DEHYDROG__45__RXN', 'RXN__45__12565', 'RXN__45__1348', '_2__45__DEHYDROPANTOATE__45__REDUCT__45__RXN', 'RXN__45__20677', 'RXN__45__12568__45__HEXANAL__47__CO__45__A__47__NADP__47____47__HEXANOYL__45__COA__47__NADPH__47__PROTON__46__45__46__'], 'Intersection of cardinality minimal completions': ['ISOBUTYRYL__45__COA__45__MUTASE__45__RXN', 'RXN1F__45__148', 'HEXANOL__45__RXN', 'FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN', 'RXN__45__12565', 'BENZYL__45__ALC__45__DEHYDROGENASE__45__RXN'], 'Union of cardinality minimal completions': ['FUMARATE__45__REDUCTASE__45__NADH__45__RXN', 'RXN__45__12559', 'DIMETHYLMALATE__45__DEHYDROGENASE__45__RXN', 'RXN__45__19913', 'RXN__45__12568__45__HEXANAL__47__CO__45__A__47__NAD__47____47__HEXANOYL__45__COA__47__NADH__47__PROTON__46__43__46__', 'RXN__45__9772', '_1__46__4__46__1__46__21__45__RXN__45__L__45__ASPARTATE__47__NADP__47__WATER__47____47__OXALACETIC_ACID__47__AMMONIUM__47__NADPH__47__PROTON__46__62__46__', 'RXN__45__9929', 'TRANSENOYLCOARED__45__RXN__45__HEXANOYL__45__COA__47__NADP__47____47__CPD0__45__2121__47__NADPH__47__PROTON__46__42__46__', 'RXN__45__12565', 'RXN__45__12299', 'ASPARTASE__45__RXN', '_2__45__DEHYDROPANTOATE__45__REDUCT__45__RXN', 'BENZALDEHYDE__45__DEHYDROGENASE__45__NAD__43____45__RXN', 'RXN__45__20677', 'RXN__45__12568__45__HEXANAL__47__CO__45__A__47__NADP__47____47__HEXANOYL__45__COA__47__NADPH__47__PROTON__46__45__46__', 'ISOBUTYRYL__45__COA__45__MUTASE__45__RXN', 'RXN__45__1744', 'RXN__45__18378', 'HEXANOL__45__RXN', 'GLYCOLALD__45__DEHYDROG__45__RXN', 'RXN__45__11270', '_1__46__4__46__1__46__21__45__RXN__45__L__45__ASPARTATE__47__NAD__47__WATER__47____47__OXALACETIC_ACID__47__AMMONIUM__47__NADH__47__PROTON__46__60__46__', '_2__46__6__46__1__46__70__45__RXN', 'RXN1F__45__148', 'RXN__45__7699__45__CPD__45__14915__47____47__CPD0__45__2121__47__WATER__46__27__46__', 'RXN__45__8633', 'PANTOATE__45__4__45__DEHYDROGENASE__45__RXN', 'FESO3OXI__45__RXN', 'SULFITE__45__REDUCT__45__RXN', 'R__45__DEHYDROPANTOATE__45__DEHYDROGENASE__45__RXN', 'RXN__45__18377', 'RXN__45__1348', 'BENZYL__45__ALC__45__DEHYDROGENASE__45__RXN', 'RXN__45__13279__45__CPD__45__14915__47__NADP__47____47__K__45__HEXANOYL__45__COA__47__NADPH__47__PROTON__46__44__46__', 'RXN__45__9336', 'RXN__45__12570', 'PYROXALTRANSAM__45__RXN', 'DIHYDRONEOPTERIN__45__MONO__45__P__45__DEPHOS__45__RXN', 'RXN__45__10065', 'RXN__45__11269', 'KETOPANTOALDOLASE__45__RXN']}
    # resCamelina = {'Draft network file': '/home/antoine/INRAE/Work/Gap_filling/clean_CamelinaFusion.sbml', 'Seeds file': '/home/antoine/INRAE/Work/Gap_filling/seedsPlants.sbml', 'Targets file': '/home/antoine/INRAE/Work/Gap_filling/targetsCamelina.sbml', 'Unproducible targets': ['CPD__45__7832_c', 'D__45__GALACTONATE_c', 'LEU_c', 'ERYTHRITOL_c', 'L__45__ASPARTATE_c', 'VAL_c', 'CYS_c', 'LINOLENIC_ACID_c', 'ARG_c', 'ILE_c', 'ADIPATE_c', 'TRP_c', 'TYR_c', 'L__45__THREONATE_c', 'LYS_c', 'MET_c', 'XYLITOL_c', 'PROPIONATE_c', 'PUTRESCINE_c', 'THR_c', 'GLY_c', 'CPD__45__10726_c', 'SINAPATE_c', 'PRO_c', 'CPD__45__373_c', 'ALPHA__45__MALTOSE_c', 'CPD__45__110_c', 'HIS_c', 'SPERMIDINE_c', 'ASCORBATE_c', 'ASN_c', 'SUC_c', 'CELLOBIOSE_c', 'L__45__SELENOCYSTEINE_c', 'LINOLEIC_ACID_c', 'GLN_c', 'GLT_c'], 'Repair db file': '/home/antoine/INRAE/Work/Gap_filling/metacyc.sbml', 'Unreconstructable targets': ['CPD__45__7832_c', 'LINOLENIC_ACID_c', 'ADIPATE_c', 'ERYTHRITOL_c', 'L__45__SELENOCYSTEINE_c', 'CPD__45__10726_c', 'LINOLEIC_ACID_c', 'SINAPATE_c'], 'Reconstructable targets': ['D__45__GALACTONATE_c', 'LEU_c', 'L__45__ASPARTATE_c', 'VAL_c', 'CYS_c', 'ARG_c', 'ILE_c', 'TRP_c', 'L__45__THREONATE_c', 'TYR_c', 'LYS_c', 'MET_c', 'XYLITOL_c', 'PROPIONATE_c', 'PUTRESCINE_c', 'THR_c', 'GLY_c', 'PRO_c', 'CPD__45__373_c', 'ALPHA__45__MALTOSE_c', 'CPD__45__110_c', 'HIS_c', 'SPERMIDINE_c', 'ASCORBATE_c', 'ASN_c', 'SUC_c', 'CELLOBIOSE_c', 'GLN_c', 'GLT_c'], 'Essential reactions': {'D__45__GALACTONATE_c': ['GALACTONOLACTONASE__45__RXN'], 'LEU_c': ['FESO3OXI__45__RXN'], 'L__45__ASPARTATE_c': [], 'VAL_c': [], 'CYS_c': ['FESO3OXI__45__RXN'], 'ARG_c': [], 'ILE_c': [], 'TRP_c': [], 'L__45__THREONATE_c': ['RXN__45__12874', 'RXN__45__12873'], 'TYR_c': [], 'LYS_c': [], 'MET_c': ['FESO3OXI__45__RXN'], 'XYLITOL_c': [], 'PROPIONATE_c': [], 'PUTRESCINE_c': [], 'THR_c': [], 'GLY_c': [], 'PRO_c': [], 'CPD__45__373_c': ['RIBOSE__45__1__45__DEHYDROGENASE__45__NADP__43____45__RXN'], 'ALPHA__45__MALTOSE_c': [], 'CPD__45__110_c': [], 'HIS_c': [], 'SPERMIDINE_c': [], 'ASCORBATE_c': [], 'ASN_c': [], 'SUC_c': [], 'CELLOBIOSE_c': [], 'GLN_c': [], 'GLT_c': []}, 'One minimal completion': ['RXN__45__7971', 'RXN0__45__268', 'RXN__45__12874', 'RXN__45__7790', 'RXN__45__11315__45__XYLITOL__47__NAD__47____47__D__45__Xylose__47__NADH__47__PROTON__46__34__46__', '_1__46__1__46__1__46__250__45__RXN', 'RXN__45__11108', 'RXN__45__1744', 'RXN__45__12873', 'CELLOBIOSE__45__PHOSPHORYLASE__45__RXN__45__CELLOBIOSE__47__Pi__47____47__GLC__45__1__45__P__47__Glucopyranose__46__37__46__', 'RXN__45__1981', '_1__46__3__46__3__46__12__45__RXN', 'GALACTONOLACTONASE__45__RXN', 'MALTOSE__45__SYNTHASE__45__RXN__45__GLC__45__1__45__P__47__WATER__47____47__ALPHA__45__MALTOSE__47__Pi__46__32__46__', 'FESO3OXI__45__RXN', 'GALACTODEHYDROG__45__RXN', 'RIBOSE__45__1__45__DEHYDROGENASE__45__NADP__43____45__RXN', 'DIHYDRONEOPTERIN__45__MONO__45__P__45__DEPHOS__45__RXN'], 'Intersection of cardinality minimal completions': ['RXN__45__12874', 'RXN__45__12873', '_1__46__3__46__3__46__12__45__RXN', 'RXN__45__1981', 'GALACTONOLACTONASE__45__RXN', 'FESO3OXI__45__RXN', 'RIBOSE__45__1__45__DEHYDROGENASE__45__NADP__43____45__RXN', 'DIHYDRONEOPTERIN__45__MONO__45__P__45__DEPHOS__45__RXN'], 'Union of cardinality minimal completions': ['RXN0__45__268', 'RXN__45__8773', '_3__45__DEHYDRO__45__L__45__GULONATE__45__2__45__DEHYDROGENASE__45__RXN__45__3__45__KETO__45__L__45__GULONATE__47__NADP__47____47__CPD__45__334__47__NADPH__47__PROTON__46__45__46__', 'LACCOA__45__RXN', 'RXN__45__1744', 'FUMARATE__45__REDUCTASE__45__NADH__45__RXN', 'RXN0__45__703', 'RXN__45__20926', '_1__46__3__46__3__46__12__45__RXN', 'DIHYDRONEOPTERIN__45__MONO__45__P__45__DEPHOS__45__RXN', 'PROPIONLACT__45__RXN', 'RXN__45__7790', 'CELLOBIOSE__45__PHOSPHORYLASE__45__RXN__45__CELLOBIOSE__47__Pi__47____47__GLC__45__1__45__P__47__ALPHA__45__GLUCOSE__46__37__46__', 'PROPANOATECOA__45__LIGASE__45__RXN', 'D__45__XYLULOSE__45__REDUCTASE__45__RXN', 'CELLOBIOSE__45__PHOSPHORYLASE__45__RXN__45__CELLOBIOSE__47__Pi__47____47__GLC__45__1__45__P__47__GLC__46__27__46__', 'CELLOBIOSE__45__PHOSPHORYLASE__45__RXN__45__CELLOBIOSE__47__Pi__47____47__GLC__45__1__45__P__47__Glucopyranose__46__37__46__', 'GALACTONOLACTONASE__45__RXN', 'KETOBUTFORMLY__45__RXN', 'PROPIONATE__45__COA__45__TRANSFERASE__45__RXN', 'FESO3OXI__45__RXN', 'GALACTODEHYDROG__45__RXN', 'SULFITE__45__REDUCT__45__RXN', 'RXN__45__7971', '_5__46__4__46__99__46__16__45__RXN', 'D__45__ARABINITOL__45__4__45__DEHYDROGENASE__45__RXN', 'RXN__45__11315__45__XYLITOL__47__NAD__47____47__D__45__Xylose__47__NADH__47__PROTON__46__34__46__', 'RXN__45__5041', 'RXN__45__12874', 'RXN__45__11315', 'RXN__45__7938', 'RXN__45__11108', 'RXN__45__12873', 'RXN__45__15583', 'GALACTODEHYDROG__45__RXN__45__ALPHA__45__D__45__GALACTOSE__47__NAD__47____47__D__45__GALACTONO__45__1__45__4__45__LACTONE__47__NADH__47__PROTON__46__59__46__', 'GALACTODEHYDROG__45__RXN__45__GALACTOSE__47__NAD__47____47__D__45__GALACTONO__45__1__45__4__45__LACTONE__47__NADH__47__PROTON__46__51__46__', 'MALTOSE__45__SYNTHASE__45__RXN__45__GLC__45__1__45__P__47__WATER__47____47__ALPHA__45__MALTOSE__47__Pi__46__32__46__', 'L__45__XYLULOSE__45__REDUCTASE__45__RXN', 'DEHYDRO__45__L__45__GULONATE__45__DECARBOXYLASE__45__RXN', 'RIBOSE__45__1__45__DEHYDROGENASE__45__NADP__43____45__RXN', '_1__46__1__46__1__46__250__45__RXN', 'RXN__45__7972', 'RXN__45__11825', 'RXN__45__1981']}

    # resTomato['Union of cardinality minimal completions']
    # dicoUpset = {
    #     "Tomato" : resTomato['Union of cardinality minimal completions'],
    #     "Kiwi" : resKiwi['Union of cardinality minimal completions'],
    #     "Cucumber" : resCucumber['Union of cardinality minimal completions'],
    #     "Cherry" : resCherry['Union of cardinality minimal completions'],
    #     "Camelina" : resCamelina['Union of cardinality minimal completions']}
    # make_upsetplot(WD, dicoUpset, "UpsetPlot_Gap-Filling")