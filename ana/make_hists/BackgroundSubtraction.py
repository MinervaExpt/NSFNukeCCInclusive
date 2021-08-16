import os, sys
import ROOT
from ROOT import PlotUtils

# Bkg Subtration: subtract absolute number of bkg events
def BkgSubtractionMC(mc_hist, bkg_hist, var):
  #MC subtracted
  h_background_subtracted_mc = mc_hist.Clone("h_background_subtracted_mc_%s"%var)
  h_background_subtracted_mc.Add(bkg,-1)

  return h_background_subtracted_mc


## MAIN  

ROOT.TH1.AddDirectory(False)

infile = ROOT.TFile("/minerva/data/users/anezkak/ME6A_T3Fe/ROOT_files/Hists_EventSelection_Bkg_ML_ME6A_sys_t3_z26_AntiNu.root","READ")
# purity numerator = efficiency numerator
# for purity subtraction
infile2 = ROOT.TFile("/minerva/data/users/anezkak/ME6A_T3Fe/ROOT_files/Hists_Efficiency_ML_ME6A_sys_t3_z26_AntiNu.root","READ")

# Scale factor to scale MC to data
mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()
mcScale =  dataPOT/mcPOT

# files to write results in
out1 = ROOT.TFile("BackgroundSubtracted_EventSelection_bkgSubtract.root","RECREATE")
out2 = ROOT.TFile("BackgroundSubtracted_EventSelection_Purity.root","RECREATE")

vars = ["Enu", "x"]

for var in vars:
  # read in needed histograms
  reco = infile.Get("selected_mc_reco_%s"%var)
  bkg = infile.Get("selected_mc_reco_bkg_%s"%var)
  data = infile.Get("selected_data_reco_%s"%var)
  purityNum = infile2.Get("h_mc_%s"%var)

  ##########################################################
  # Bkg Subtration: subtract absolute number of bkg events
  h_background_subtracted_mc = BkgSubtractionMC(reco, bkg, var)

  h_background_mc_scale = bkg.Clone("h_background_mc_scale_%s"%var)
  h_background_mc_scale.Scale(mcScale)
  # subtract scaled bkg from data
  data.ClearAllErrorBands()
  data.AddMissingErrorBandsAndFillWithCV(h_background_mc_scale)
  h_background_subtracted_data = data.Clone("h_background_subtracted_data_%s"%var)
  h_background_subtracted_data.Add(h_background_mc_scale,-1) 

  # write all the histogramt to file
  out1.cd()
  h_background_subtracted_mc.Write()
  h_background_subtracted_data.Write()

  #######################################################
  # Purity subtraction
  out2.cd()
  purity = purityNum.Clone("purity_%s"%var)
  purity.Divide(purity,reco, 1.0, 1.0, "B")
  # binomial because numerator is the subset of the denominator, to make the stat error correct
 
  data.ClearAllErrorBands()
  data.AddMissingErrorBandsAndFillWithCV(reco)

  h_purity_background_subtracted_mc = reco.Clone("h_purity_background_subtracted_mc_%s"%var)
  h_purity_background_subtracted_data = data.Clone("h_purity_background_subtracted_data_%s"%var)

  h_purity_background_subtracted_mc.Multiply(h_purity_background_subtracted_mc, purity)
  h_purity_background_subtracted_data.Multiply(h_purity_background_subtracted_data, purity)

  h_purity_background_subtracted_mc.Write()
  h_purity_background_subtracted_data.Write()
  purity.Write()

out1.Close()
out2.Clear()
out1.Close()
out2.Close() 

'''
#MC subtracted
h_background_subtracted_mc = reco.Clone("h_background_subtracted_mc_%s"%var)
h_background_subtracted_mc.Add(bkg,-1)

# Data subtracted
# scale backgroud to data
h_background_mc_scale = bkg.Clone("h_background_mc_scale_%s"%var)
h_background_mc_scale.Scale(mcScale)
  
# subtract scaled bkg from data
data.ClearAllErrorBands()
data.AddMissingErrorBandsAndFillWithCV(h_background_mc_scale)
h_background_subtracted_data = data.Clone("h_background_subtracted_data_%s"%var)
h_background_subtracted_data.Add(h_background_mc_scale,-1) 

# write all the histogramt to file
out1.cd()
h_background_subtracted_mc.Write()
h_background_mc_scale.Write()
h_background_subtracted_data.Write()
'''

'''
# Purity based subtraction
out2.cd()
purity = purityNum.Clone("purity")
purity.Divide(purity,reco, 1.0, 1.0, "B")
# binomial because numerator is the subset of the denominator, to make the stat error correct
 
data.ClearAllErrorBands()
data.AddMissingErrorBandsAndFillWithCV(reco)

h_purity_background_subtracted_mc = reco.Clone("h_purity_background_subtracted_mc_%s"%var)
h_purity_background_subtracted_data = data.Clone("h_purity_background_subtracted_data_%s"%var)

h_purity_background_subtracted_mc.Multiply(h_purity_background_subtracted_mc, purity)
h_purity_background_subtracted_data.Multiply(h_purity_background_subtracted_data, purity)

h_purity_background_subtracted_mc.Write()
h_purity_background_subtracted_data.Write()
purity.Write()
'''


