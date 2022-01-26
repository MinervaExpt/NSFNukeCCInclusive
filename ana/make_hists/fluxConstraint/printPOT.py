import ROOT
import os,sys
from ROOT import PlotUtils

infile= ROOT.TFile("/minerva/data2/users/anezkak/flux4Daisy_files_2022/CombinedPlaylists_Neutrinos_Tracker_nosys.root")

mcPOT = infile.Get("MCPOT").GetVal()

print("MC POT: "+str(mcPOT))

raw_input("Done")

