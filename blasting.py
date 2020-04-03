# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE de Bordeaux
# Mars 2020

import cobra
from os.path import join
import subprocess


def get_genes_reactions(_path):
    '''Browse an SBML model to get the genes and associated reactions.
    
    Args :
        _path -- the reference GEM in SBML format.
    Return : 
        genes_and_reactions -- a dictionary with the gene name as key and a dictionary
        as value, containing the reactions and their formulas,
        and the results of blastp.
    '''
    genes_and_reactions = {}
    model = cobra.io.read_sbml_model(_path)
    for i in range(len(model.genes)):
        gene = model.genes[i].name
        try:
            reac = model.reactions[i]
        except IndexError:
            print("Reaction not found : ",i)
            reac = "NA"
            pass
        genes_and_reactions[gene] = {reac.id:[reac.reaction]}
    return genes_and_reactions

def blast_p(dictionary, workDir, subjectFile, queryFile):
    """Runs multiple blastp between a subject file and each protein of the query file.
    
    ARGS : 
        dictionary -- the dictionary containing the genes and reactions,
        which will gather the results of the blastp.
        workDir -- the working directory where the files are, and where the proteins files
        will be temporarily stored.
        subjectPATH -- the subject file for the blastp.
        queryPATH -- the query file for the blastp corresponding to the BiModel object given.
    RETURN : 
        dictionary -- the input dictionary containing the result of the blastp command for each gene.
    """
    newDir = workDir + "Proteins/"
    queryDir = workDir + queryFile
    subjectDir = workDir + subjectFile
    subprocess.run(["mkdir", newDir])
    
    ###Creation of one file per protein of the query###
    print("\nCreating the small fasta files...")
    fileFasta = open(queryDir)
    queryFasta = fileFasta.read()
    for seq in queryFasta.split(">"):
        geneName = seq.split("\n")[0]
        f = open(newDir+geneName+".fa", "w")
        f.write(">"+seq)
        f.close()
    fileFasta.close()
    print("...done !")
    
    ###Blastp###
    print("\nLaunching the blastp !")
    i = 1
    x = len(dictionary.keys())
    for gene in dictionary.keys():
        if i%50==0:
            print("Protein %i out of %i" %(i,x))
        i+=1
        requestBlastp = [
            "blastp",
            "-subject",
            subjectDir,
            "-query",
            newDir+"in_"+gene+".fa",
            "-outfmt",
            "10 delim=, qseqid qlen sseqid slen length nident pident score evalue bitscore"]
        dictionary[gene]['BlastP'] = subprocess.run(requestBlastp, capture_output=True).stdout.decode('ascii').split("\n")[:-1]
        print(dictionary[gene])
    subprocess.run(["rm", "-rf", newDir])


if __name__=='__main__':
    ###Files and working directory###
    WD = '/home/asa/INRAE/Work/Drafts/Data/Tests/'
    modelGem = 'AraGEM3.xml'
    modelGemFasta = 'genomic.in.fasta'
    modelGemFastaTest = 'testBlast.in.fasta'
    tomatoFasta = 'ITAG4.0_proteins.fasta'
    tomatoTestFasta = 'TomatoTest.fasta'
    
    ###Pipeline###
    core_info = get_genes_reactions(WD+modelGem)
    blast_p(core_info, WD, tomatoFasta, modelGemFastaTest)