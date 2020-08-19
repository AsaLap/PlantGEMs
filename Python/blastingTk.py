# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020

"""Project Aborted"""

from tkinter import *
from tkinter import filedialog

import blasting

blastWindow = Tk()
blastWindow.title("Draft reconstruction with SBML model")
blastWindow.geometry("600x450")

#Working directory for the files and results
# wd = Label(blastWindow, text = "Working Directory :")
# WD = Entry(blastWindow)
# wd.grid(row = 1, column = 1, sticky = W)
# WD.grid(row = 1, column = 2)

#The model (browser)
file = filedialog.askopenfile(parent=blastWindow, mode='r', title='Choose the model file', filetypes=[("xml files", "*.xml"), ("sbml files", "*.sbml")])
if file != None:
    data = file.read()
    file.close()

#The model
model = Label(blastWindow, text = "COBRA compatible model filename :")
MODEL = Entry(blastWindow)
model.grid(row = 2, column = 1, sticky = W)
MODEL.grid(row = 2, column = 2)
modelFasta = Label(blastWindow, text = "Model fasta filename :")
modelFASTA = Entry(blastWindow)
modelFasta.grid(row = 3, column = 1, sticky = W)
modelFASTA.grid(row = 3, column = 2)

#The target
targetFasta = Label(blastWindow, text = "Target fasta filename :")
targetFASTA = Entry(blastWindow)
targetFasta.grid(row = 4, column = 1, sticky = W)
targetFASTA.grid(row = 4, column = 2)

#The parameters
#Name
modelName = Label(blastWindow, text = "The model name :")
modelNAME = Entry(blastWindow)
modelNAME.insert(0, "My_new_model")
modelName.grid(row = 5, column = 1, sticky = W)
modelNAME.grid(row = 5, column = 2)

#Blast
var1 = IntVar()
def blastOrNot():
    pass
blast = Checkbutton(blastWindow, text='Blastp', variable=var1, onvalue=1, offvalue=0, command=blastOrNot)
blast.grid(row = 6, column = 1, columnspan=2, sticky = W)

#Identity
identity = Label(blastWindow, text = "The percentage of identity :")
IDENTITY = Entry(blastWindow)
IDENTITY.insert(0, 50)
identity.grid(row = 7, column = 1, sticky = W)
IDENTITY.grid(row = 7, column = 2)

#Difference
diff = Label(blastWindow, text = "The percentage of difference :")
DIFF = Entry(blastWindow)
DIFF.insert(0, 30)
diff.grid(row = 8, column = 1, sticky = W)
DIFF.grid(row = 8, column = 2)

#E-Value
e_value = Label(blastWindow, text = "The maximum E-Value :")
EVAL = Entry(blastWindow)
EVAL.insert(0, 1e-100)
e_value.grid(row = 9, column = 1, sticky = W)
EVAL.grid(row = 9, column = 2)

#Coverage
coverage = Label(blastWindow, text = "The minimum coverage :")
COVERAGE = Entry(blastWindow)
COVERAGE.insert(0, 20)
coverage.grid(row = 10, column = 1, sticky = W)
COVERAGE.grid(row = 10, column = 2)

#Bit-Score
bitScore = Label(blastWindow, text = "The minimum Bit-Score :")
BITSCORE = Entry(blastWindow)
BITSCORE.insert(0, 300)
bitScore.grid(row = 11, column = 1, sticky = W)
BITSCORE.grid(row = 11, column = 2)

#Launch button
button = Button(blastWindow, text = "Go", command = blasting)
button.grid(row = 12, column = 1, sticky = NS)

#Quit button
button = Button(blastWindow, text = "Quit", command = blastWindow.destroy)
button.grid(row = 12, column = 2, sticky = NS)

blastWindow.mainloop()