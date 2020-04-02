# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE de Bordeaux
# Mars 2020

import cobra
from os.path import join
import subprocess

#path for WD :
WD = '/home/asa/INRAE/Work/Drafts/Data/Tests/'
aragem = 'AraGEM3.xml'
araFasta = 'genomic.in.fasta'
tomatoFasta = 'ITAG4.0_proteins.fasta'
aragemPATH = WD + aragem
aragemFastaPATH = WD + araFasta
tomatoFastaPATH = WD + tomatoFasta
testQueryPATH = WD + "testBlast.in.fasta"

class Model:
    def __init__(self, _name):
        self.name = _name
        self.dictionary = {}
        self.genes = []
        self.reactions = []
    
    def buil_model(self, _path):
        model = cobra.io.read_sbml_model(_path)
        for i in range(len(model.genes)):
            gene = model.genes[i].name
            try:
                reaction = model.reactions[i].name
            except IndexError:
                print("Reaction not found : ",i)
                reaction = "NA"
                pass
            if gene in self.dictionary.keys():
                self.dictionary[gene].append(reaction)
            else:
                self.dictionary[gene] = [reaction]
        self.genes = self.dictionary.keys()
        self.reactions = self.dictionary.values()
        
def blast_p(model, subjectPATH, queryPATH):
    '''
    "Query" is the sequence to search in several sequences
    of database/fasta file which is named "subject".
    The query here is each protein sequence of the model organism's GEM.
    '''
    fileFasta = open(queryPATH)
    queryFasta = fileFasta.read()
    for seq in queryFasta.split(">"):
        geneName = seq.split("\n")[0]
        f = open(WD+"Proteins/"+geneName+".txt", "w")
        f.write(">"+seq)
        f.close()
    fileFasta.close()
    dico_res = {}
    for gene in model.genes:
        ###The request used for Blastp
        requestBlastp = ["blastp", "-subject", subjectPATH, "-query", WD+"Proteins/"+"in_"+gene+".txt", "-outfmt", "10 delim=, qseqid sseqid qlen length slen nident pident score evalue bitscore"]
        ###The request passed via subprocess to the blastp command in bash
        dico_res[gene] = subprocess.run(requestBlastp, capture_output=True).stdout.decode('ascii').split("\n")[:-1]
        
        ###Debug
        # print(subprocess.run(requestBlastp, capture_output=True).stdout.decode('ascii').split("\n"))
        # print(requestBlastp)
        # print(dico_res[gene])
        
    
    
if __name__=='__main__':
    aragemModel = Model('AraGEM')
    aragemModel.buil_model(aragemPATH)
    print(len(aragemModel.genes))
    blast_p(aragemModel, aragemFastaPATH, aragemFastaPATH)
    # blast_p(aragemModel, aragemFastaPATH, testQueryPATH)