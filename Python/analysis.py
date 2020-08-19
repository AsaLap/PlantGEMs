# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used for analysing the reconstructed models."""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import utils

###Analyse des fichiers de Sylvain
def read_Sylvain_csv(file):
    dico_res = {}
    res = utils.read_csv(file, "\t")
    for line in res:
        dico_res[line[0]] = {"Length" : int(line[1]), "Type" : line[2], "Comments" : line[3:]}
    return dico_res

def add(variable, x = 1):
    return variable + x

def count_firstLvl(dico_res, size = 2):
    len_sup2 = 0
    dico_count = {"Biosynthesis" : 0, "Degradation" : 0, "Energy-Metabolism" : 0, "Detoxification" : 0,
                  "Activation-Inactivation-Interconversion" : 0, "Macromolecule-Modification" : 0, "Metabolic-Clusters" : 0}
    for key in dico_res.keys():
        if dico_res[key]["Length"] > size:
            dico_count[dico_res[key]["Type"]] = add(dico_count[dico_res[key]["Type"]])
    return dico_count

tomatoCount = count_firstLvl(read_Sylvain_csv("/home/asa/INRAE/Work/Analyse/perSpecies/Tomato-80.tsv"), 0)
kiwiCount = count_firstLvl(read_Sylvain_csv("/home/asa/INRAE/Work/Analyse/perSpecies/Kiwi-80.tsv"), 0)
cucumberCount = count_firstLvl(read_Sylvain_csv("/home/asa/INRAE/Work/Analyse/perSpecies/Cucumber-80.tsv"), 0)
cherryCount = count_firstLvl(read_Sylvain_csv("/home/asa/INRAE/Work/Analyse/perSpecies/Cherry-80.tsv"), 0)
camelinaCount = count_firstLvl(read_Sylvain_csv("/home/asa/INRAE/Work/Analyse/perSpecies/Camelina-80.tsv"), 0)
dico_all = {"Tomato" : tomatoCount, 
            "Kiwi" : kiwiCount,
            "Cucumber" : cucumberCount,
            "Cherry" : cherryCount,
            "Camelina" : camelinaCount}



def multi_barplot(dico_all):
    """Code coming mainly from matplotlib.org
    (https://matplotlib.org/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py)"""
    
    labels = list(dico_all["Tomato"].keys())
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
    rects4 = ax.bar(x + width/size, cucumber, width/size, label='Cerise')
    rects5 = ax.bar(x + 2*width/size, cucumber, width/size, label='Cameline')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Scores')
    ax.set_title('Scores by group and gender')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()


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

    # fig.tight_layout()

    plt.show()

multi_barplot(dico_all)