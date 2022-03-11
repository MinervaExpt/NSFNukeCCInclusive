import ROOT
import os,sys
from ROOT import PlotUtils

infile= ROOT.TFile("/pnfs/minerva/persistent/users/anezkak/default_analysis_loc/flux/AntiNu_sys/Hists_EventSelection_minervame6D_FluxConstraint_optim_sys_t99_z99_AntiNu.root")

mcPOT = infile.Get("MCPOT").GetVal()

print("MC POT: "+str(mcPOT))

raw_input("Done")

