import os, sys
import ROOT
from ROOT import PlotUtils
from ROOT import TParameter

# Bkg Subtration: subtract absolute number of bkg events
def BkgSubtractionMC(mc_hist, bkg_hist, petal, var):
  #MC subtracted
  h_background_subtracted_mc = mc_hist.Clone("h_bkg_subtracted_mc_daisy_%s_%s"%(petal,var))
  h_background_subtracted_mc.Add(bkg_hist,-1)

  return h_background_subtracted_mc


## MAIN  

ROOT.TH1.AddDirectory(False)

targetID = 99
targetZ = 99

pwd = sys.argv[1]
plist = sys.argv[2]

infile = ROOT.TFile(str(pwd)+"EventSelection_daisy_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ),"READ")

# Scale factor to scale MC to data
mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()
mcScale =  dataPOT/mcPOT

# files to write results in
out1 = ROOT.TFile("BkgSubtracted_EventSelection_daisy_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ),"RECREATE")

vars = ["Enu", "x"]

for var in vars:
  for petal in range(0,12):
    # read in needed histograms
    reco = infile.Get("selected_mc_reco_daisy_%s_%s"%(petal,var))
    bkg = infile.Get("selected_mc_reco_bkg_daisy_%s_%s"%(petal,var))
    data = infile.Get("selected_data_reco_daisy_%s_%s"%(petal,var))

    ##########################################################
    # Bkg Subtration: subtract absolute number of bkg events
    h_background_subtracted_mc = BkgSubtractionMC(reco, bkg, petal, var)

    h_background_mc_scale = bkg.Clone("h_bkg_mc_scale_%s"%var)
    h_background_mc_scale.Scale(mcScale)
    # subtract scaled bkg from data
    data.ClearAllErrorBands()
    data.AddMissingErrorBandsAndFillWithCV(reco)#h_background_mc_scale)
    h_background_subtracted_data = data.Clone("h_bkg_subtracted_data_daisy_%s_%s"%(petal,var))
    h_background_subtracted_data.Add(h_background_mc_scale,-1) 

    # write all the histogramt to file
    out1.cd()
    h_background_subtracted_mc.Write()
    h_background_subtracted_data.Write()

# write down POTs
dataPOTout = TParameter('double')("DataPOT", dataPOT)
mcPOTout = TParameter('double')("MCPOT", mcPOT)
dataPOTout.Write()
mcPOTout.Write()


out1.Close()
out1.Close()