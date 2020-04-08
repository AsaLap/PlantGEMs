# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Mars 2020

import matplotlib.pyplot as plt
from statistics import mean

import blasting


def graph_identity(data):
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


def graph_e_value(data):
    # list_e_value = []
    # for key in data.keys():
    #     for res in data[key]:
    #         list_e_value.append(float(res.split(",")[8]))
    # print(len(list_e_value))
    # print(list_e_value.index(min(list_e_value)), max(list_e_value), mean(list_e_value))
    res_plot = {}
    puissance = 10
    ##Python ne descend pas plus bas que 10**-320 !##
    while puissance <=320:
        treshold = 10**-puissance
        count = 0
        for key in data.keys():
            for res in data[key]:
                if float(res.split(",")[8]) >= treshold:
                    count+=1
        res_plot[treshold] = count
        puissance+=10
    x, y = list(res_plot.keys()), list(res_plot.values())
    return x, y


if __name__=="__main__":
    ###Files and working directory###
    WDtom = '/home/asa/INRAE/Work/Drafts/Tomato_Arabidopsis/'
    WDkiw = '/home/asa/INRAE/Work/Drafts/Kiwi_Arabidopsis/'
    WDcuc = '/home/asa/INRAE/Work/Drafts/Cucumber_Arabidopsis/'
    WDche = '/home/asa/INRAE/Work/Drafts/Cherry_Arabidopsis/'
    
    ###Making plots###
    #Loads for the plot
    tom = blasting.load_obj(WDtom + "resBlastp")
    kiw = blasting.load_obj(WDkiw + "resBlastp")
    cuc = blasting.load_obj(WDcuc + "resBlastp")
    # che = blasting.load_obj(WDche + "resBlastp")
    
    #------Identity------#
    # tomato = graph_identity(tom)
    # kiwi = graph_identity(kiw)
    # cucumber = graph_identity(cuc)
    # # cherry = graph_identity(che)
    # plt.plot(tomato[0], tomato[1], 'r', label = 'Tomato')
    # plt.plot(kiwi[0], kiwi[1], 'g', label = 'Kiwi')
    # plt.plot(cucumber[0], cucumber[1], 'b', label = 'Cucumber')
    # # plt.plot(cherry[0], cherry[1], 'y', label = 'Cherry' )
    # plt.ylabel('Number of mappings')
    # plt.xlabel('Percentage of identity')
    # plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
    # # plt.savefig(WD + '4plots.png', dpi=300)
    # plt.show()
    
    #------E-Value------#
    tomato = graph_e_value(tom)
    print(tomato[0])
    plt.plot(tomato[0], tomato[1], 'r', label = 'Tomato')
    plt.show()