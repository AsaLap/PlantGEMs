# coding: utf8
# python 3.8.2
# Antoine Laporte
# Université de Bordeaux - INRAE Bordeaux
# Reconstruction de réseaux métaboliques
# Mars - Aout 2020

from tkinter import *

import blasting

blastWindow = Tk()
blastWindow.title("Draft reconstruction with SBML model")
blastWindow.geometry("600x400")

#Working directory for the files and results
wd = Label(blastWindow, text = "Working Directory :")
WD = Entry(blastWindow)
wd.grid(row = 1, column = 1, sticky = W)
WD.grid(row = 1, column = 2)

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



button = Button(blastWindow, text = "Quit", command = blastWindow.destroy)
button.grid(row = 5, column = 1, columnspan=2, sticky = NS)

blastWindow.mainloop()