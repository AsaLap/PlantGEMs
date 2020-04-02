# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE de Bordeaux
# Mars 2020

import cobra
from os.path import join


#path for arabidopsis :
arabidopsis = '/home/asa/INRAE/Work/Drafts/Data/Arabidopsis/'
aragem = 'AraGEM.xml'
araFasta = 'genomic.in.fasta'

def get_genes_reactions(model):
    '''
        FUNCTION : browse model to get the genes and associated reactions
        INPUT : a model
        OUTPUT : a dictionary with the number as key, gene and associated reaction as value (list)
    '''
    dico_gene_reaction = {}
    for i in range(len(model.genes)):
        gene = model.genes[i].name
        try:
            reaction = model.reactions[i].name
        except IndexError:
            print("Reaction not found : ",i)
            reaction = "NA"
            pass
        dico_gene_reaction[gene] = reaction
    return dico_gene_reaction
    print("Metabolites : ", len(model.metabolites))
    print("Reactions : ", len(model.reactions))
    print("Genes : ", len(model.genes))
    print(model.genes)
    
if __name__=='__main__':
    aragemModel = cobra.io.read_sbml_model(join(arabidopsis, aragem))
    aragemDico = get_genes_reactions(aragemModel)
    print(len(aragemDico.keys()))
    # print(aragemDico)