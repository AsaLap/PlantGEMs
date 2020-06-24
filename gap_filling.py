# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020
"""This file is used for the gap filling of the previously reconstructed models."""

import subprocess
import string
import os

import utils


def count(file, letters):
    res = {}
    for letter in letters:
        res[letter] = 0
        for i in file:
            if letter == i:
                res[letter] += 1
    return res


def clean(path):
    file = utils.read_file(path)
    res = ""
    for line in file:
        if ">" not in line:
            res += line.rstrip().upper()
    return res


def stats(res):
    total = 0
    for key in res.keys():
        total += res[key]
    for key in res.keys():
        print("%s : %.0f = %.2f percent" %(key, res[key], res[key] * 100 / total))


if __name__=="__main__":
    letterProt = "ACDEFGHIKLMNPQRSTVWY*"
    letterDNA = "CGTA"
    letterRNA = "CGUA"
    # total_cherry_prot = clean("/home/asa/INRAE/Work/Gap_filling/test.faa")
    total_cherry_prot = clean("/home/asa/INRAE/Work/Raw_Data/Tomato/ITAG/ITAG4.0_proteins.fasta")
    res = count(total_cherry_prot, letterProt)
    stats(res)

    total_cherry_CDS = clean("/home/asa/INRAE/Work/Raw_Data/Tomato/ITAG/ITAG4.0_CDS.fasta")
    res = count(total_cherry_CDS, letterDNA)
    stats(res)