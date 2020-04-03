# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE de Bordeaux
# Mars 2020

import cobra
from os.path import join
import subprocess
import pickle


def get_genes_reactions(_path):
    '''Browse an SBML model to get the genes and associated reactions.
    
    Args :
        _path -- the reference GEM in SBML format.
    Return : 
        gene_dic -- a dictionary with the gene name as key and two dictionaries as values,
        one for the reactions (IDs and formula) and the other one for the blastp result.
    '''
    gene_dic = {}
    model = cobra.io.read_sbml_model(_path)
    for gene in model.genes:
        gene_dic[gene.name]={'Reactions':{}}
        for i in gene.reactions:
            gene_dic[gene.name]['Reactions'][i.id] = i.reaction
    return gene_dic

def blast_p(gene_dic, workDir, subjectFile, queryFile):
    """Runs multiple blastp between a subject file and each protein of the query file.
    
    ARGS : 
        gene_dic -- the dictionary containing the genes and reactions,
        which will gather the results of the blastp.
        workDir -- the working directory where the files are, and where the proteins files
        will be temporarily stored.
        subjectPATH -- the subject file for the blastp.
        queryPATH -- the query file for the blastp corresponding to the BiModel object given.
    RETURN : 
        gene_dic -- the input dictionary containing the result of the blastp command for each gene.
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
    x = len(gene_dic.keys())
    for gene in gene_dic.keys():
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
        gene_dic[gene]['BlastP'] = subprocess.run(requestBlastp, capture_output=True).stdout.decode('ascii').split("\n")[:-1]
    subprocess.run(["rm", "-rf", newDir])
    return gene_dic

def save_obj(obj, path):
    with open(path + '.pkl', 'wb+') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(path):
    with open(path + '.pkl', 'rb') as f:
        return pickle.load(f)

if __name__=='__main__':
    ###Files and working directory###
    WD = '/home/asa/INRAE/Work/Drafts/Data/Tests/'
    modelGem = 'AraGEM3.xml'
    modelGemFasta = 'genomic.in.fasta'
    modelGemFastaTest = 'testBlast.in.fasta'
    tomatoFasta = 'ITAG4.0_proteins.fasta'
    tomatoFastaTest = 'TomatoTest.fasta'
    
    ###Pipeline###
    core_info = get_genes_reactions(WD+modelGem)
    core_info = blast_p(core_info, WD, tomatoFastaTest, modelGemFasta)
    save_obj(core_info, WD + "dictionary")
    # test = load_obj(WD + "dictionary")