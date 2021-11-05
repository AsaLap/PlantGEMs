# coding: utf8
# python 3.8
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux Aquitaine - DRC
# PhD : Genome scale metabolic reconstruction of several fruits
# 11/2021 - 10/2024

"""This file is the main, calling the others files to create an entirely new metabolic network."""

from blasting import Blasting

if __name__ == "__main__":
    # Tomate
    # tomate_blast = Blasting("tomate",
    #                         "/home/asa/INRAE/These/Reconstructions/Aracyc/aracyc.sbml",
    #                         "/home/asa/INRAE/These/Reconstructions/Aracyc/aracyc.fasta",
    #                         "/home/asa/INRAE/These/Reconstructions/Tomate/tomato.fasta")
    # tomate_blast.build()

    # Raisin
    vitis_blast = Blasting("raisin", "/home/asa/INRAE/Test/")
    # vitis_blast.build()
    # crashtest_blast = Blasting("tomate",
    #                            "/home/asa/INRAE/These/Reconstructions/Aracyc/aracyc.sbml",
    #                            "/home/asa/INRAE/These/Reconstructions/Aracyc/aracyc_crashtest.fasta",
    #                            "/home/asa/INRAE/These/Reconstructions/Crashtest/crashtest.fasta")
    # crashtest_blast.build()