import cobra
from utils import list_reactions, make_upsetplot

WDtomCyc = '/home/asa/INRAE/StageMaster_2020/Work/FichiersRelancePipeline/blasting/Tomato_Aracyc/'
WDkiwCyc = '/home/asa/INRAE/StageMaster_2020/Work/FichiersRelancePipeline/blasting/Kiwi_Aracyc/'
WDcucCyc = '/home/asa/INRAE/StageMaster_2020/Work/FichiersRelancePipeline/blasting/Cucumber_Aracyc/'
WDcheCyc = '/home/asa/INRAE/StageMaster_2020/Work/FichiersRelancePipeline/blasting/Cherry_Aracyc/'
WDraisinCyc = '/home/asa/INRAE/These/Reconstructions/Raisin/raisin_aracyc/'

tomatoDraftCyc = cobra.io.load_json_model(WDtomCyc + "Tomato.json")
kiwiDraftCyc = cobra.io.load_json_model(WDkiwCyc + "Kiwi.json")
cucumberDraftCyc = cobra.io.load_json_model(WDcucCyc + "Cucumber.json")
cherryDraftCyc = cobra.io.load_json_model(WDcheCyc + "Cherry.json")
raisinDraftCyc = cobra.io.load_json_model(WDraisinCyc + "vitis_aracyc_blast.json")

dicoUpset = {"Tomato": list_reactions(tomatoDraftCyc),
             "Kiwi": list_reactions(kiwiDraftCyc),
             "Cucumber": list_reactions(cucumberDraftCyc),
             "Cherry": list_reactions(cherryDraftCyc),
             "Raisin": list_reactions(raisinDraftCyc)}
make_upsetplot("/home/asa/INRAE/These/Reconstructions/", "UpsetPlot 5 fruits", dicoUpset,
               "Total reactions intersections of reconstructed models")
