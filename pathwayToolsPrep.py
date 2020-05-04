# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# 2020

def read_gff():
    print("TODO")


def read_eggNog():
    print("TODO")


def make_pf():
    print("TODO")


def make_dat():
    print("TODO")



if __name__=="__main__":
    ###Files and working directory###
    WDtom = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Tomato/'
    WDkiw = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Kiwi/'
    WDcuc = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Cucumber/'
    WDche = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Cherry/'
    WDcam = '/home/asa/INRAE/Work/Plant-GEMs/PathwayToolsData/Camelina/'
    
    ###Fasta files
    tomatoFasta = 'ITAG4.0_proteins.fasta'
    kiwiFasta = 'Actinidia_chinensis.Red5_PS1_1.69.0.pep.all.fa'
    cucumberFasta = 'Gy14_pep_v2.fa'
    cherryFasta = 'PRUAV_Regina_pep.fa'
    camelinaFasta = 'GCF_000633955.1_Cs_protein.faa'
    
    ###GFF files
    tomatoGFF = "ITAG4.0_gene_models.gff"
    kiwiGFF = "Actinidia_chinensis.Red5_PS1_1.69.0.44.gff3"
    cucumberGFF = "Gy14_gene_gff_v2.gff"
    cherryGFF = "PRUAV_Regina.gff3"
    camelinaGFF = "GCF_000633955.1_Cs_genomic.gff"
    
    ###eggNOG files
    tomatoEgg = "eggNOG_annotations.tsv"
    kiwiEgg = "eggNOG_annotations.tsv"
    cucumberEgg = "eggNOG_annotations.tsv"
    cherryEgg = "eggNOG_annotations.tsv"
    camelinaEgg = "eggNOG_annotations.tsv"