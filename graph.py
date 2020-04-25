# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Mars 2020

import matplotlib.pyplot as plt
from statistics import mean
import csv

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


def get_score(data):
    list_value = []
    for key in data.keys():
        for res in data[key]:
            list_value.append(float(res.split(",")[7]))
    return list_value


def get_e_value(data):
    list_value = []
    for key in data.keys():
        for res in data[key]:
            list_value.append(float(res.split(",")[8]))
    return list_value


def get_bit_score(data):
    list_value = []
    for key in data.keys():
        for res in data[key]:
            list_value.append(float(res.split(",")[9]))
    return list_value


def get_all_scores(organism, data):
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


def write_csv(WD, list_value, name):
    with open(WD + name + '.csv', 'w', newline = '') as file:
        writer = csv.writer(file)
        for f in list_value:
            writer.writerow(f)


if __name__=="__main__":
    ###Files and working directory###
    WDtom = '/home/asa/INRAE/Work/Plant-GEMs/Drafts/Tomato_Arabidopsis/'
    WDkiw = '/home/asa/INRAE/Work/Plant-GEMs/Drafts/Kiwi_Arabidopsis/'
    WDcuc = '/home/asa/INRAE/Work/Plant-GEMs/Drafts/Cucumber_Arabidopsis/'
    WDche = '/home/asa/INRAE/Work/Plant-GEMs/Drafts/Cherry_Arabidopsis/'
    WDcam = '/home/asa/INRAE/Work/Plant-GEMs/Drafts/Camelina_Arabidopsis/'
    
    ###Making plots###
    #Loading the data
    tom = blasting.load_obj(WDtom + "resBlastp")
    kiw = blasting.load_obj(WDkiw + "resBlastp")
    cuc = blasting.load_obj(WDcuc + "resBlastp")
    che = blasting.load_obj(WDche + "resBlastp")
    cam = blasting.load_obj(WDcam + "resBlastp")
    
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