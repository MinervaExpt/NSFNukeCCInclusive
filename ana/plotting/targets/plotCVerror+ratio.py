import ROOT
from ROOT import *
import os,sys
from ROOT import PlotUtils
import ctypes
import math as m

def plotCVerrorRatio(mnvPlotter, mc_hist, bkg_mc, data_hist, mcScale, mat, var):
    
    mc_hist.Scale(mcScale)
    mc_hist.Scale(mc_hist.GetNormBinWidth(), "width" )

    bkg_mc.Scale(mcScale)
    bkg_mc.Scale(bkg_mc.GetNormBinWidth(), "width" )

    data_hist.Scale(data_hist.GetNormBinWidth(), "width" )

    mc_hist.GetXaxis().CenterTitle()
    mc_hist.GetYaxis().CenterTitle()


    #Create a TCanvas on which to draw plots and split it into 2 panels
    overall = TCanvas("Data/MC for " + var, "plot", 800, 800)
    top = TPad("Overlay", "Overlay", 0, bottomFraction, 1, 1)
    bottom = TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction + margin)
    #Thou shalt Draw() new TPads lest they be blank!
    top.Draw()
    bottom.Draw()

    top.cd()

    if var == "Enu":
        mc_hist.GetYaxis().SetTitle("Events/GeV")
    if var == "x":
        mc_hist.GetYaxis().SetTitle("Events per unit x")
    
    mc_hist_err = mc_hist.GetCVHistoWithError(True, False).Clone() 
    # MC
    mc_hist.SetLineWidth(3)
    mc_hist.SetLineColor(2)
    # MC error
    mc_hist_err.SetFillColor(ROOT.kRed-10)
    mc_hist_err.SetFillStyle(1001)
    mc_hist_err.SetMarkerStyle(0)
    mc_hist_err.SetLineWidth(3)
    mc_hist_err.SetLineColor(2)
    # MC bkg
    bkg_mc.SetFillColor(14)
    bkg_mc.SetFillStyle(3005)
    bkg_mc.SetLineColor(1)
    bkg_mc.SetLineWidth(1)

    data_hist.SetMarkerSize(1)
    data_hist.SetMarkerStyle(20)

    mc_hist.Draw("HIST")
    mc_hist.SetMinimum(0.0001) # to get rid of the zero
    mc_hist.SetMaximum(1.4*mc_hist.GetMaximum())
    mc_hist_err.Draw("E2 SAME")
    bkg_mc.Draw("HIST SAME")
    mc_hist.Draw("HIST SAME")
    data_hist.Draw("X0 E1 SAME")

    mc_hist.GetYaxis().SetTitleOffset(1.1)

    legend = TLegend(0.45,0.65,0.8,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(data_hist, " Data", "lep")
    legend.AddEntry(mc_hist_err, " Simulation", "fl")
    legend.AddEntry(bkg_mc, " Background prediction", "f")
    legend.SetTextFont(62)
    legend.Draw()

    mnv.AddHistoTitle("Target %s %s: Event Selection"%(targetID, trueZ), 0.04, 1)
    mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.87)
    
    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.4)
    
    dataWithStatErrors = data_hist.GetCVHistoWithError().Clone()
    ratio = dataWithStatErrors.Clone() # stat

    ratio.Divide(ratio,mc_hist.GetCVHistoWithError()) # stat
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetTitleOffset(0.3)
    ratio.GetXaxis().CenterTitle()
    ratio.GetYaxis().CenterTitle()
    ratio.Draw("X0")

    ratio.SetMaximum(1.5)
    ratio.SetMinimum(0.5)
    ratio.GetXaxis().SetTitleSize(titleSize)
    ratio.GetXaxis().SetLabelSize(labelSize)
    ratio.GetYaxis().SetLabelSize(0.11)
    ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions
    ratio.GetYaxis().SetTitleSize(titleSize)

    # Copied from MINERvA-101-Cross-Section/backgroundStack.py
    # same as in void MnvPlotter::DrawDataMCRatio() in MnvPlotter
    # systematic error centered at y = 1
    sysError = mc_hist.GetTotalError(False, True, False) # False for stat error, True for as frac, False for covAreaNorm
    for whichBin in range(1, sysError.GetXaxis().GetNbins()+1):
        sysError.SetBinError(whichBin, max(sysError.GetBinContent(whichBin), 1e-9))
        sysError.SetBinContent(whichBin, 1)

    sysError.SetFillColor(ROOT.kRed-10)
    sysError.SetFillStyle(1001)
    sysError.SetMarkerStyle(0)
    sysError.Draw("E2 SAME")

    # line at 1
    line = TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
    line.SetLineColor(46)
    line.SetLineWidth(lineSize)
    line.SetLineStyle(9)
    line.Draw()

    ratio.Draw("X0 SAME")

    if var == "Enu":
        ratio.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")

    if var == "x":
        ratio.GetXaxis().SetTitle("Reconstructed Bjorken x")
    
    # ratio equivalent to this:
    #mnvPlotter.DrawDataMCRatio(data, totalMC, 1, True, True, 0.7, 1.5, "Data/Total MC" )

    overall.Print("T" + targetID + "_z" + targetZ + "_" + "eventSelection"+ "_" + var+".png")


# ---------------------------------------------------------------------------------------------------

gROOT.SetBatch() #Don't render histograms to a window.  Also gets filled areas correct.

bottomFraction = 0.15
margin = 0.116 #Tuned by hand
labelSize = 0.15
lineSize = 2
titleSize = 0.16

TH1.AddDirectory(False)

# Run like: python bkgStack+ratio.py 3 26
# To plot plots for target 3 iron
targetID = sys.argv[1] 
targetZ = sys.argv[2]

infile= ROOT.TFile("Hists_EventSelection_Bkg_ML_ME6A_sys_t%s_z%s_AntiNu.root"%(targetID, targetZ))
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
    # ----------------------------------------------------------------------------
    #                             HISTOGRAMS
    # ----------------------------------------------------------------------------
    # Read in histograms
    mc_hist = infile.Get("selected_mc_reco_%s"%var)
    data_hist = infile.Get("selected_data_reco_%s"%var)
    bkg_mc = infile.Get("selected_mc_reco_bkg_%s"%var)

    mcPOT = infile.Get("MCPOT").GetVal()
    dataPOT = infile.Get("DataPOT").GetVal()

    mcScale =  dataPOT/mcPOT

    plotCVerrorRatio(mnv, mc_hist, bkg_mc, data_hist, mcScale, mat, var)


   


raw_input("Done")
