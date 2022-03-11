import ROOT
import os,sys
from ROOT import PlotUtils

infile1= ROOT.TFile("/pnfs/minerva/persistent/users/anezkak/default_analysis_loc/Nu_flux/Hists_EventSelection_minervame1L_FluxConstraint_optimPetal_sys_t99_z99_Nu.root")
infile2= ROOT.TFile("/pnfs/minerva/persistent/users/anezkak/default_analysis_loc/Nu_flux/Hists_EventSelection_minervame1M_FluxConstraint_optimPetal_sys_t99_z99_Nu.root")
infile3= ROOT.TFile("/pnfs/minerva/persistent/users/anezkak/default_analysis_loc/Nu_flux/Hists_EventSelection_minervame1N_FluxConstraint_optimPetal_sys_t99_z99_Nu.root")
infile4= ROOT.TFile("/pnfs/minerva/persistent/users/anezkak/default_analysis_loc/Nu_flux/Hists_EventSelection_minervame1O_FluxConstraint_optimPetal_sys_t99_z99_Nu.root")
infile5= ROOT.TFile("/pnfs/minerva/persistent/users/anezkak/default_analysis_loc/Nu_flux//Hists_EventSelection_minervame1P_FluxConstraint_optimPetal_sys_t99_z99_Nu.root")

mcPOTL = infile1.Get("MCPOT").GetVal()
mcPOTM = infile2.Get("MCPOT").GetVal()
mcPOTN = infile3.Get("MCPOT").GetVal()
mcPOTO = infile4.Get("MCPOT").GetVal()
mcPOTP = infile5.Get("MCPOT").GetVal()

print("MC POT: "+str(mcPOTL))
print("MC POT: "+str(mcPOTM))
print("MC POT: "+str(mcPOTN))
print("MC POT: "+str(mcPOTO))
print("MC POT: "+str(mcPOTP))

mcPOT = mcPOTL + mcPOTM + mcPOTN + mcPOTO + mcPOTP
print("Water MC POT: "+str(mcPOT))

raw_input("Done")

