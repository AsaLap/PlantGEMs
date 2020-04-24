# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Mars 2020

import cobra
from os.path import join
import subprocess
import pickle
import time
import matplotlib.pyplot as plt
from statistics import mean
import copy

import graph


def save_obj(obj, path):
    """Saves the dictionary of Blastp results in a pickle file."""
    with open(path + '.pkl', 'wb+') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(path):
    """Loads a dictionary of Blastp results stored in a pickle file."""
    with open(path + '.pkl', 'rb') as f:
        return pickle.load(f)


def blast_run(workDir, model, queryFile, subjectFile):
    """Runs multiple blasts between a subject file and each protein of the query file.
    
    ARGS : 
        model -- the GEM of reference.
        workDir -- the working directory where the files are, and where the proteins files
        will be temporarily stored.
        subjectFile -- the subject file for the blast.
        queryFile -- the query file for the blast corresponding to the reference GEM's CDS fasta.
    RETURN : 
        blast_res -- dictionary containing the result of the blastp command for each gene of the model.
    """
    ###Concatenation of string to get the exact paths###
    newDir = workDir + "Proteins_tmp/"
    queryDir = workDir + queryFile
    subjectDir = workDir + subjectFile
    ###Creation of temporary directory containing one file per protein of the query###
    subprocess.run(["mkdir", newDir])
    print("\nCreating the individual CDS fasta files...")
    fileFasta = open(queryDir)
    queryFasta = fileFasta.read()
    for seq in queryFasta.split(">"):
        geneName = seq.split("\n")[0]
        f = open(newDir+geneName+".fa", "w")
        f.write(">"+seq)
        f.close()
    fileFasta.close()
    print("...done !")
    ###Blast###
    print("\nLaunching the blast !")
    i, x = 1, len(model.genes)
    total_time = lap_time = time.time()
    blast_res = {}
    for gene in model.genes:
        if i%10==0:
            print("Protein %i out of %i\nTime : %f s" %(i,x, time.time() - lap_time))
            lap_time = time.time()
        i+=1
        requestBlast = [
            "blastp",
            "-subject",
            subjectDir,
            "-query",
            newDir+"in_"+gene.name+".fa",
            "-outfmt",
            "10 delim=, qseqid qlen sseqid slen length nident pident score evalue bitscore"]
        blast_res[gene.name] = subprocess.run(requestBlast, capture_output=True).stdout.decode('ascii').split("\n")[:-1]
    subprocess.run(["rm", "-rf", newDir])
    print("Blast done !\nTotal time : %f s" %(time.time() - total_time))
    return blast_res


def select_genes(blast_res, treshold = 75, len_diff = 20, e_val = 1e-100):
    """Select the subject organism's genes with a good score.
    
    ARGS:
        blast_res -- the dictionary with the results of the blastp.
        treshold -- the treshold value of identity to select the subject genes.
    RETURN:
        dico_genes -- a dictionary with model gene as key and corresponding subject 
        key and coverage value as value.
    """
    dico_genes = {}
    for key in blast_res.keys():
        for res in blast_res[key]:
            spl = res.split(",")
            len_subject = float(spl[3])
            len_query = [float(spl[1]) * (100 - len_diff) / 100, float(spl[1]) * (100 + len_diff) / 100]
            if float(spl[6]) >= treshold and len_subject >= len_query[0] and len_subject <= len_query[1]:
            # length = float(spl[1]) * 100 / float(spl[3])
            # if float(spl[6]) >= treshold and length >= 100 - len_diff and length <= 100 + len_diff:
                try:
                    if float(spl[8]) <= e_val:
                        dico_genes[key].append(spl[2])
                except KeyError:
                    if float(spl[8]) <= e_val:
                        dico_genes[key] = [spl[2]]
    return dico_genes


def drafting(model, dico_genes, model_name):
    """Creates the new model for the subject organism.
    
    ARGS:
        model -- the COBRA model used for the reconstruction.
        dico_genes -- the dictionary containing the correspondance
        between genes of the model and the subject.
        model_name -- a name for the new model (string)
    RETURN:
        new_model -- the new COBRA model automatically generated.
    """
    new_model = cobra.Model(model_name)
    #Browsing the model reactions and associating the subject's genes to them
    for reac in model.reactions:
        to_add = []
        to_search = reac.gene_reaction_rule.split(" or ")
        for gene in to_search:
            try:
                to_add += dico_genes[gene]
            except KeyError:
                pass
        string_reaction_rule = " or ".join(to_add)
        if string_reaction_rule:
            x = copy.deepcopy(reac)
            x.gene_reaction_rule = string_reaction_rule
            new_model.add_reactions([x])
    return new_model


def pipeline(WD, ref_gem, queryFile, subjectFile, modelName, identity = 70, coverage = 20, EVal = 10**-100):
    model = cobra.io.read_sbml_model(WD + ref_gem)
    # blast_res = blast_run(WD, model, queryFile, subjectFile)
    # save_obj(blast_res, WD + "resBlastp")
    blast_res = load_obj(WD + "resBlastp")
    dico_genes = select_genes(blast_res, identity, coverage, EVal)
    new_model = drafting(model, dico_genes, modelName)
    new_model.add_metabolites(model.metabolites)
    
    ###Printing of verifications
    #Test blast_res, search for the genes with no matches with blastp
    no_results = []
    for key in blast_res.keys():
        if not blast_res[key]:
            no_results.append(key)
    print("The",len(no_results),"genes that have no matches : ", no_results)
    
    #Counting of different values
    nb_values = []
    for val in dico_genes.values():
        for i in val:
            nb_values.append(i)
    print(
        "Nb of genes in ref model : %s\n\
Nb of reactions in ref model : %s\n\
Nb of genes in the new model : %s\n\
Nb of reactions in the new model : %s\n\
Nb of genes in dico_genes : %s\n\
Nb of values in dico_genes : %s\
 (without doublons : %s)"
        %(len(model.genes),
        len(model.reactions),
        len(new_model.genes),
        len(new_model.reactions),
        len(dico_genes.keys()),
        len(nb_values),
        len(set(nb_values))))
    print("----------------------------------------")
    return new_model
    

def help_treshold(WD, model, modelFasta, subjectFasta, name, identity, coverage):
    """Function to help choosing the treshold value."""
    
    listValues = []
    for i in range(1,101):
        draft = pipeline(WD, model, modelFasta, subjectFasta, name, identity, coverage, 10**(-i))
        listValues.append([i, len(draft.genes), len(draft.reactions)])
    print(listValues)
    listValues.insert(0,["EValue (1e-x)", "Nb genes", "Nb reactions"])
    graph.write_csv(WD, listValues, name + "Treshold")


def data_venn(WD, model, name):
    list_id = []
    for reac in model.reactions:
        list_id.append([reac.id])
    graph.write_csv(WD, list_id, name + "_id_reac")

if __name__=='__main__':
    ###Files and working directory###
    WDtom = '/home/asa/INRAE/Work/Drafts/Tomato_Arabidopsis/'
    WDkiw = '/home/asa/INRAE/Work/Drafts/Kiwi_Arabidopsis/'
    WDcuc = '/home/asa/INRAE/Work/Drafts/Cucumber_Arabidopsis/'
    WDche = '/home/asa/INRAE/Work/Drafts/Cherry_Arabidopsis/'
    WDcam = '/home/asa/INRAE/Work/Drafts/Camelina_Arabidopsis/'
    WDara = '/home/asa/INRAE/Work/Drafts/Arabidopsis/'
    
    WDtests = '/home/asa/INRAE/Work/Drafts/Tests/'
    
    aragem = 'AraGEM3.xml'
    aragemFasta = 'genomic.in.fasta'
    tomatoFasta = 'ITAG4.0_proteins.fasta'
    kiwiFasta = 'Hongyang_pep_v2.0.fa'
    cucumberFasta = 'Gy14_pep_v2.fa'
    cherryFasta = 'PRUAV_Regina_CDS.fa'
    camelinaFasta = 'GCF_000633955.1_Cs_protein.faa'
    
    ###Main###
    # #For the tomato
    tomatoDraft = pipeline(WDtom, aragem, aragemFasta, tomatoFasta, "Tomato")
    # #For the kiwifruit
    # kiwiDraft = pipeline(WDkiw, aragem, aragemFasta, kiwiFasta, "Kiwi")
    # #For the cucumber
    # cucumberDraft = pipeline(WDcuc, aragem, aragemFasta, cucumberFasta, "Cucumber")
    # #For the cherry
    # cherryDraft = pipeline(WDche, aragem, aragemFasta, cherryFasta, "Cucumber")
    # #For the camelina
    # camelinaDraft = pipeline(WDcam, aragem, aragemFasta, camelinaFasta, "Camelina")
    
    ###Get the reactions id for the Venn diagramm###
    # data_venn(WDtom, tomatoDraft, "Tomato")
    # data_venn(WDkiw, kiwiDraft, "Kiwi")
    # data_venn(WDcuc, cucumberDraft, "Cucumber")
    # data_venn(WDche, cherryDraft, "Cherry")
    # data_venn(WDcam, camelinaDraft, "Camelina")
    # data_venn(WDara, cobra.io.read_sbml_model(WDara + aragem), "Arabidopsis")
    
    ###Help to choose the treshold###
    # help_treshold(WDtom, aragem, aragemFasta, tomatoFasta, "Tomato", 70, 20)
    # help_treshold(WDkiw, aragem, aragemFasta, kiwiFasta, "Kiwi", 70, 20)
    # help_treshold(WDcuc, aragem, aragemFasta, cucumberFasta, "Cucumber", 70, 20)
    # help_treshold(WDche, aragem, aragemFasta, cherryFasta, "Cherry", 70, 20)
    # help_treshold(WDcam, aragem, aragemFasta, camelinaFasta, "Camelina", 70, 20)