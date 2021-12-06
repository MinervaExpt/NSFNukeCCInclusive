import os, sys
import ROOT
from ROOT import PlotUtils
from ROOT import TParameter
import math as m

# Background subtraction using plastic scale factors
# Description:
# MC: Subtract absolute number of events predicted
# Data: Subtract POT scaled, plastic scaled number of events predicted from MC

# Bkg Subtration: subtract absolute number of bkg events
def BkgSubtractionMC(mc_hist, bkg_hist, var):
  #MC subtracted
  h_background_subtracted_mc = mc_hist.Clone("h_bkg_subtracted_mc_%s"%var)
  h_background_subtracted_mc.Add(bkg_hist,-1)

  return h_background_subtracted_mc

## MAIN  

ROOT.TH1.AddDirectory(False)

targetID = sys.argv[1] 
targetZ = sys.argv[2]

infile = ROOT.TFile("/minerva/data2/users/anezkak/ME6A_Targets/Hists_EventSelection_Bkg_ML_ME6A_sys_t%s_z%s_AntiNu.root"%(targetID, targetZ),"READ")
infileUntuned = ROOT.TFile.Open("/minerva/data2/users/anezkak/ME6A_Targets/Hists_PlasticBkg_sys_t%s_z%s_AntiNu.root"%(targetID, targetZ), "READ")#%(str(targetID),str(targetZ)))
scale = ROOT.TFile("/minerva/data2/users/anezkak/ME6A_Targets/Plastic_ScaleFactors_t%s_z%s_minervame6A.root"%(targetID, targetZ), "READ")

# Scale factor to scale MC to data
mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()
mcScale =  dataPOT/mcPOT

# files to write results in
out1 = ROOT.TFile("/minerva/data2/users/anezkak/ME6A_Targets/Hists_BkgSubtracted_EventSelection_sys_t%s_z%s_AntiNu.root"%(targetID, targetZ),"RECREATE")

vars = ["Enu", "x"]

for var in vars:
  # read in needed histograms
  reco = infile.Get("selected_mc_reco_%s"%var)
  bkg = infile.Get("selected_mc_reco_bkg_%s"%var)
  bkgPlasticUS = infile.Get("selected_mc_USplastic_%s"%var)
  bkgPlasticDS = infile.Get("selected_mc_DSplastic_%s"%var)
  bkgOther = infile.Get("selected_mc_other_%s"%var)
  data = infile.Get("selected_data_reco_%s"%var)

  # plastic scaling factors
  scaleUS = scale.Get("scaleFactor_US_%s"%var)
  scaleDS = scale.Get("scaleFactor_DS_%s"%var)

  hists_US_regUS = infileUntuned.Get("US_regUS_Iron_%s"%var)
  hists_DS_regDS = infileUntuned.Get("DS_regDS_Iron_%s"%var)

  integralUS = hists_US_regUS.Integral()
  statErrUS = 1/ m.sqrt(integralUS)
  for bin in range(scaleUS.GetNbinsX() + 1):
    scaleUS.SetBinError(bin,0+statErrUS)

  # Downstream: Remove original stat error and update it with 1/sqrt(N) in each bin
  integralDS = hists_DS_regDS.Integral()
  statErrDS = 1/ m.sqrt(integralDS)
  for bin in range(scaleDS.GetNbinsX() + 1):
    scaleDS.SetBinError(bin,0+statErrDS)

  ##########################################################
  # Bkg Subtraction: subtract absolute number of bkg events

  # subtract absolute number of bkg events predicted
  # all events at the same time
  bkg_subtracted_mc = BkgSubtractionMC(reco, bkg, var)
  
  # -------------------------------------------------------------
  # SUBTRACT SCALED BKG FROM DATA
  # -------------------------------------------------------------
  # multiply bkgPlasticUS and plastic DS by the corresponding scaling factor
  bkgPlasticUS_CHScaled = bkgPlasticUS.Clone("h_bkgPlasticUS_CHScaled_%s"%var)
  bkgPlasticDS_CHScaled = bkgPlasticDS.Clone("h_bkgPlasticDS_CHScaled_%s"%var)

  bkgPlasticUS_CHScaled.Multiply(bkgPlasticUS_CHScaled, scaleUS)
  bkgPlasticDS_CHScaled.Multiply(bkgPlasticDS_CHScaled, scaleDS)

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

  data.ClearAllErrorBands()
  data.AddMissingErrorBandsAndFillWithCV(reco)

  bkg_subtracted_data = data.Clone("h_bkg_subtracted_data_%s"%var)
  bkg_subtracted_data.Add(bkgAdded_mc_scale,-1) 

  bkg_subtracted_data_notConstrained = data.Clone("h_bkg_subtracted_data_notConstrained_%s"%var)
  bkg.Scale(mcScale)
  bkg_subtracted_data_notConstrained.Add(bkg,-1) 

  # write all the histograms to file
  out1.cd()
  
  # subtracted
  bkg_subtracted_mc.Write()
  bkg_subtracted_data.Write()
  bkg_subtracted_data_notConstrained.Write()

# write down POTs
dataPOTout = TParameter(float)("DataPOT", dataPOT)
mcPOTout = TParameter(float)("MCPOT", mcPOT)
dataPOTout.Write()
mcPOTout.Write()

out1.Close()
out1.Close()
raw_input("Done")