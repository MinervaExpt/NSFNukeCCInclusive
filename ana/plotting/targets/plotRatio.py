import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle

targetID = sys.argv[1] 
targetZ = sys.argv[2]

infile= ROOT.TFile("Hists_BkgSubtracted_EventSelection_sys_t%s_z%s_AntiNu.root"%(targetID, targetZ))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale =  dataPOT/mcPOT

mat = "Fe"
trueZ = "Iron"

if targetZ == 26:
    trueZ = "Iron"
    mat = "Fe"

if targetZ == 82:
    trueZ = "Lead"
    mat = "Pb"

if targetZ == 6:
    trueZ = "Carbon"
    mat = "C"

vars = ["Enu", "x"]

for var in vars:

    mc_hist = infile.Get("h_bkg_subtracted_mc_%s"%var)
    data_hist = infile.Get("h_bkg_subtracted_data_%s"%var)

    if var == "Enu":
        mc_hist.GetXaxis().SetTitle("Reconstructed Neutrino Energy (GeV)")
        mc_hist.GetYaxis().SetTitle("Events/GeV")

    if var == "x":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")

    '''
    if myvariable == "xfine":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")


    if myvariable == "xBrian":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")

    if myvariable == "x09":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
    ''' 
    mc_hist.GetXaxis().CenterTitle()
    mc_hist.GetYaxis().CenterTitle()
    #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
    #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
    #mc_hist.GetXaxis().SetRangeUser(0,0.001)

    #mnv.ApplyStyle(1)
    mnv.DrawDataMCRatio(data_hist, mc_hist, mcScale)

    mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)

    '''
    Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
    If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
    If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
    '''
    
    mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)
    mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)

    canvas1.Modified()
    canvas1.Print("ME6A_T%s%s_BkgSubtractedEventSelectionRatio_%s_sys.png"%(targetID,mat,var))

   


data_hist.SetDirectory(0)
mc_hist.SetDirectory(0)

raw_input("Done")
