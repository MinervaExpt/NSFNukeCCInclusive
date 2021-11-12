import os, sys
import ROOT
from ROOT import PlotUtils
from ROOT import TParameter

# Bkg Subtration: subtract absolute number of bkg events
def BkgSubtractionMC(mc_hist, bkg_hist, var):
  #MC subtracted
  h_background_subtracted_mc = mc_hist.Clone("h_bkg_subtracted_mc_%s"%var)
  h_background_subtracted_mc.Add(bkg_hist,-1)

  return h_background_subtracted_mc


## MAIN  

ROOT.TH1.AddDirectory(False)

infile = ROOT.TFile("Hists_EventSelection_Bkg_ML_ME6A_sys_t3_z26_AntiNu.root","READ")
infile2 = ROOT.TFile("Plastic_ScaleFactors_t3_z26_minervame6A.root", "READ")

# Scale factor to scale MC to data
mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()
mcScale =  dataPOT/mcPOT

# files to write results in
out1 = ROOT.TFile("BackgroundSubtracted_EventSelection_bkgSubtract_wPlastic.root","RECREATE")
#out2 = ROOT.TFile("BackgroundSubtracted_EventSelection_Purity.root","RECREATE")

vars = ["Enu", "x"]

for var in vars:
  # read in needed histograms
  reco = infile.Get("selected_mc_reco_%s"%var)
  bkg = infile.Get("selected_mc_reco_bkg_%s"%var)
  bkgPlasticUS = infile.Get("selected_mc_USplastic_%s"%var)
  bkgPlasticDS = infile.Get("selected_mc_DSplastic_%s"%var)
  bkgOther = infile.Get("selected_mc_other_%s"%var)
  data = infile.Get("selected_data_reco_%s"%var)
  #purityNum = infile.Get("selected_mc_reco_signal_%s"%var)

  # plastic scaling factors
  scaleUS = infile2.Get("scaleFactor_US_%s"%var)
  scaleDS = infile2.Get("scaleFactor_DS_%s"%var)

  ##########################################################
  # Bkg Subtration: subtract absolute number of bkg events

  #1 subtract other background 
  bkg_subtracted_mc = BkgSubtractionMC(reco, bkgOther, var)

  # multiply bkgPlasticUS and plastic DS by the corresponding scaling factor
  bkgPlasticUS_CHScaled = bkgPlasticUS.Clone("h_bkgPlasticUS_CHScaled_%s"%var)
  bkgPlasticDS_CHScaled = bkgPlasticDS.Clone("h_bkgPlasticDS_CHScaled_%s"%var)

  bkgPlasticUS_CHScaled.Multiply(bkgPlasticUS_CHScaled, scaleUS)
  bkgPlasticDS_CHScaled.Multiply(bkgPlasticDS_CHScaled, scaleDS)

  # subtract scaled US and DS plastic
  bkg_subtracted_mc = BkgSubtractionMC(bkg_subtracted_mc, bkgPlasticUS_CHScaled, var)
  bkg_subtracted_mc = BkgSubtractionMC(bkg_subtracted_mc, bkgPlasticDS_CHScaled, var)

  # save all background (other + plastic scaled US and DS)
  bkgAdded = bkgOther.Clone("h_bkgAdded_mc_%s"%var)
  bkgAdded.Add(bkgOther)
  bkgAdded.Add(bkgPlasticUS_CHScaled)
  bkgAdded.Add(bkgPlasticDS_CHScaled)

  # scale by MC scale (POT normalisation)
  otherBkg_mc_scale = bkgOther.Clone("h_bkgOther_mc_scaled_%s"%var)
  bkgPlasticUS_mc_scale = bkgPlasticUS_CHScaled.Clone("h_bkgPlasticUS_mc_scaled_%s"%var)
  bkgPlasticDS_mc_scale = bkgPlasticDS_CHScaled.Clone("h_bkgPlasticDS_mc_scaled_%s"%var)

  otherBkg_mc_scale.Scale(mcScale)
  bkgPlasticUS_mc_scale.Scale(mcScale)
  bkgPlasticDS_mc_scale.Scale(mcScale)

  # save POT scaled background
  bkgAdded_mc_scale = otherBkg_mc_scale.Clone("h_bkgAdded_mc_scaled_%s"%var)
  bkgAdded_mc_scale.Add(bkgPlasticUS_mc_scale)
  bkgAdded_mc_scale.Add(bkgPlasticDS_mc_scale)
  
  # -------------------------------------------------------------
  # SUBTRACT SCALED BKG FROM DATA
  # -------------------------------------------------------------
  data.ClearAllErrorBands()
  data.AddMissingErrorBandsAndFillWithCV(reco)

  bkg_subtracted_data = data.Clone("h_bkg_subtracted_data_%s"%var)
  bkg_subtracted_data.Add(bkgAdded_mc_scale,-1) 

  # write all the histograms to file
  out1.cd()

  # original
  bkg.Write()
  bkgOther.Write()
  bkgPlasticUS.Write()
  bkgPlasticDS.Write()

  # scaled by plastic scaling factors
  bkgPlasticUS_CHScaled.Write()
  bkgPlasticDS_CHScaled.Write()

  # MC scaled
  otherBkg_mc_scale.Write()
  bkgPlasticUS_mc_scale.Write()
  bkgPlasticDS_mc_scale.Write()

  # added backgrounds
  bkgAdded.Write()
  bkgAdded_mc_scale.Write()
  
  # subtracted
  bkg_subtracted_mc.Write()
  bkg_subtracted_data.Write()

# write down POTs
dataPOTout = TParameter(float)("DataPOT", dataPOT)
mcPOTout = TParameter(float)("MCPOT", mcPOT)
dataPOTout.Write()
mcPOTout.Write()

  #######################################################
  # Purity subtraction
  #out2.cd()
  #purity = purityNum.Clone("purity_%s"%var)
  #purity.Divide(purity,reco, 1.0, 1.0, "B")
  # binomial because numerator is the subset of the denominator, to make the stat error correct
 
  #data.ClearAllErrorBands()
  #data.AddMissingErrorBandsAndFillWithCV(reco)

  #h_purity_background_subtracted_mc = reco.Clone("h_purity_background_subtracted_mc_%s"%var)
  #h_purity_background_subtracted_data = data.Clone("h_purity_background_subtracted_data_%s"%var)

  #h_purity_background_subtracted_mc.Multiply(h_purity_background_subtracted_mc, purity)
  #h_purity_background_subtracted_data.Multiply(h_purity_background_subtracted_data, purity)

  #h_purity_background_subtracted_mc.Write()
  #h_purity_background_subtracted_data.Write()
  #purity.Write()

out1.Close()
#out2.Clear()
out1.Close()
raw_input("Done")
#out2.Close() 

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


