# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used for analysing the reconstructed models."""

import cobra
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import utils

###Comptage sur les réseaux reconstruits
def count_genes_per_reac(model):
    nb_genes = 0
    dico_count = {"0" : 0, "1" : 0, "2" : 0, "3" : 0, "4" : 0, "5" : 0,
                  "6" : 0, "7" : 0, "8" : 0, "9" : 0, "10" : 0}
    for reac in model.reactions:
        if reac.gene_reaction_rule:
            nb_genes = len(reac.gene_reaction_rule.split(" or "))
        else:
            nb_genes = 0
        if nb_genes > 10:
            nb_genes = 10
        dico_count[str(nb_genes)] = add(dico_count[str(nb_genes)])
    return dico_count

def save_count_genes(WD, name, res):
    fic = []
    for key in res.keys():
        fic.append([int(key), res[key]])
    utils.write_csv(WD, fic, name)

###Analyse des fichiers de Sylvain
def read_pathways_csv(file):
    dico_res = {}
    res = utils.read_csv(file, "\t")
    for line in res:
        if len(line) == 3:
            line.append("Other")
        dico_res[line[0]] = {"Length" : int(line[1]), "Type" : line[2], "Precision1" : line[3]}
    return dico_res

def add(variable, x = 1):
    return variable + x

def count_pathways(dico_res, precision = None, size = 2):
    """Function to gather information about pathways
    ARGS:
        dico_res -- python dictionary resulting from the read_pathways_csv() function.
        precision (str) -- None by default makes the count on the first field of dico_res,
        if you want to look more closely, put a string corresponding to one type in the first field.
        size (int) -- the minimum for a pathway to be considered as is (2 by default).
    RETURN:
        dico_count -- python dictionary in which every count is stored.
    """
    
    dico_count = {}
    for key in dico_res.keys():
        if dico_res[key]["Length"] >= size:
            if not precision:
                if dico_res[key]["Type"] in list(dico_count.keys()):
                    dico_count[dico_res[key]["Type"]] = add(dico_count[dico_res[key]["Type"]])
                else:
                    dico_count[dico_res[key]["Type"]] = 1 #First appearance
            else:
                if dico_res[key]["Type"] == precision:
                    if dico_res[key]["Precision1"] in list(dico_count.keys()):
                        dico_count[dico_res[key]["Precision1"]] = add(dico_count[dico_res[key]["Precision1"]])
                    else:
                        dico_count[dico_res[key]["Precision1"]] = 1 #First appearance
    return dico_count


def make_dico_all(path, extension, precision, size):
    """Helpfull function to go faster making graphs for the pathways analysis."""
    
    dico_all = {"Tomato" : count_pathways(read_pathways_csv(path + "Tomato" + extension), precision, size), 
            "Kiwi" : count_pathways(read_pathways_csv(path + "Kiwi" + extension), precision, size),
            "Cucumber" : count_pathways(read_pathways_csv(path + "Cucumber" + extension), precision, size),
            "Cherry" : count_pathways(read_pathways_csv(path + "Cherry" + extension), precision, size),
            "Camelina" : count_pathways(read_pathways_csv(path + "Camelina" + extension), precision, size)}
    ###Add the keys that are not in every dictionary for the plot to work###
    all_keys = []
    for key in dico_all.keys():
        for type_key in dico_all[key].keys():
            all_keys.append(type_key)
    all_keys = set(all_keys)
    for bio_type in all_keys:
        for key in dico_all.keys():
            if bio_type not in list(dico_all[key].keys()):
                dico_all[key][bio_type] = 0
    return dico_all

def multi_barplot(dico_all, x_lab, y_lab, title):
    """Code coming mainly from matplotlib.org
    (https://matplotlib.org/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py)"""
    
    labels = list(dico_all[list(dico_all.keys())[0]].keys())
    tomato = []
    kiwi = []
    cucumber = []
    cherry = []
    camelina = []
    for label in labels:
        tomato.append(dico_all["Tomato"][label])
        kiwi.append(dico_all["Kiwi"][label])
        cucumber.append(dico_all["Cucumber"][label])
        cherry.append(dico_all["Cherry"][label])
        camelina.append(dico_all["Camelina"][label])

    x = np.arange(len(labels))  # the label locations
    width = 0.75  # the width of the bars
    size = 5

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - 2*width/size, tomato, width/size, label='Tomate')
    rects2 = ax.bar(x - width/size, kiwi, width/size, label='Kiwi')
    rects3 = ax.bar(x, cucumber, width/size, label='Concombre')
    rects4 = ax.bar(x + width/size, cherry, width/size, label='Cerise')
    rects5 = ax.bar(x + 2*width/size, camelina, width/size, label='Cameline')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    plt.setp(ax.get_xticklabels(), rotation=20, horizontalalignment='right')
    fig.tight_layout()

    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)
    autolabel(rects4)
    autolabel(rects5)
    plt.show()

if __name__=="__main__":
    ###Graph for the genes per reaction association
    # save_count_genes("/home/asa/INRAE/Work/", "distribution_genes_Tomato", count_genes_per_reac(cobra.io.load_json_model("/home/asa/INRAE/Work/final_Tomato.json")))
    # save_count_genes("/home/asa/INRAE/Work/", "distribution_genes_Kiwi", count_genes_per_reac(cobra.io.load_json_model("/home/asa/INRAE/Work/final_Kiwi.json")))
    # save_count_genes("/home/asa/INRAE/Work/", "distribution_genes_Cucumber", count_genes_per_reac(cobra.io.load_json_model("/home/asa/INRAE/Work/final_Cucumber.json")))
    # save_count_genes("/home/asa/INRAE/Work/", "distribution_genes_Cherry", count_genes_per_reac(cobra.io.load_json_model("/home/asa/INRAE/Work/final_Cherry.json")))
    # save_count_genes("/home/asa/INRAE/Work/", "distribution_genes_Camelina", count_genes_per_reac(cobra.io.load_json_model("/home/asa/INRAE/Work/final_Camelina.json")))
    # dico_all_genes = {
    #     "Tomato" : count_genes_per_reac(cobra.io.load_json_model("/home/asa/INRAE/Work/final_Tomato.json")),
    #     "Kiwi" : count_genes_per_reac(cobra.io.load_json_model("/home/asa/INRAE/Work/final_Kiwi.json")),
    #     "Cucumber" : count_genes_per_reac(cobra.io.load_json_model("/home/asa/INRAE/Work/final_Cucumber.json")),
    #     "Cherry" : count_genes_per_reac(cobra.io.load_json_model("/home/asa/INRAE/Work/final_Cherry.json")),
    #     "Camelina" : count_genes_per_reac(cobra.io.load_json_model("/home/asa/INRAE/Work/final_Camelina.json"))
    #     }
    # multi_barplot(dico_all_genes, "Nombre de gènes", "Nombre de réactions", "Nombre de réactions en fonction du nombre de gènes associés")
    
    ###Graph for the pathways's types (first level) distribution between the species
    # tomatoCount = count_pathways(read_pathways_csv("/home/asa/INRAE/Work/Analyse/perSpecies/Tomato-80.tsv"), size = 3)
    # kiwiCount = count_pathways(read_pathways_csv("/home/asa/INRAE/Work/Analyse/perSpecies/Kiwi-80.tsv"), size = 3)
    # cucumberCount = count_pathways(read_pathways_csv("/home/asa/INRAE/Work/Analyse/perSpecies/Cucumber-80.tsv"), size = 3)
    # cherryCount = count_pathways(read_pathways_csv("/home/asa/INRAE/Work/Analyse/perSpecies/Cherry-80.tsv"), size = 3)
    # camelinaCount = count_pathways(read_pathways_csv("/home/asa/INRAE/Work/Analyse/perSpecies/Camelina-80.tsv"), size = 3)
    # dico_all = {"Tomato" : tomatoCount,
    #         "Kiwi" : kiwiCount,
    #         "Cucumber" : cucumberCount,
    #         "Cherry" : cherryCount,
    #         "Camelina" : camelinaCount}
    # multi_barplot(dico_all, 'Nature des voies métaboliques', 'Nombre de voies métaboliques (taille mini = 3)', 'Nature et nombre de voies métaboliques par espèce')
    
    ###Graph for the pathways's types (Biomass only) distribution between the species
    ###Biosynthesis
    # dico_all = make_dico_all("/home/asa/INRAE/Work/Analyse/perSpecies/", "-80.tsv", "Biosynthesis", 4)
    # multi_barplot(dico_all, 'Type de biosynthèse',
    # 'Nombre de voies métaboliques (taille mini = 4)',
    # 'Type et nombre de voies métaboliques pour la biosynthèse')
    
    ###Degradation
    # dico_all = make_dico_all("/home/asa/INRAE/Work/Analyse/perSpecies/", "-80.tsv", "Degradation", 4)
    # multi_barplot(dico_all, 'Type de degradation',
    # 'Nombre de voies métaboliques (taille mini = 4)',
    # 'Type et nombre de voies métaboliques pour la dégradation')
    
    ###NRJ
    # dico_all = make_dico_all("/home/asa/INRAE/Work/Analyse/perSpecies/", "-80.tsv", "Energy-Metabolism", 2)
    # multi_barplot(dico_all, 'Type de métabolisme énergétique',
    # 'Nombre de voies métaboliques (taille mini = 2)',
    # "Type et nombre de voies métaboliques pour la production d'énergie")
    
    ###Activation, inactivation and interconversion
    # dico_all = make_dico_all("/home/asa/INRAE/Work/Analyse/perSpecies/", "-80.tsv", "Activation-Inactivation-Interconversion", 2)
    # multi_barplot(dico_all, "Type d'activation, désactivation et interconversion",
    # 'Nombre de voies métaboliques (taille mini = 2)',
    # "Type et nombre de voies métaboliques pour l'activation, la désactivation et l'interconversion")
    
    ###Macromolecule modification
    # dico_all = make_dico_all("/home/asa/INRAE/Work/Analyse/perSpecies/", "-80.tsv", "Macromolecule-Modification", 4)
    # multi_barplot(dico_all, "Type de modification des macromolécules",
    # 'Nombre de voies métaboliques (taille mini = 4)',
    # "Type et nombre de voies métaboliques pour la modification des macromolécules")
    
    ###Metabolic clusters
    # dico_all = make_dico_all("/home/asa/INRAE/Work/Analyse/perSpecies/", "-80.tsv", "Metabolic-Clusters", 4)
    # multi_barplot(dico_all, "Type de clusters métaboliques",
    # 'Nombre de voies métaboliques (taille mini = 4)',
    # "Type et nombre de voies métaboliques pour les clusters métaboliques")