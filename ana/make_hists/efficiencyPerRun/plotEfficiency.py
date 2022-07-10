# File:         plotEfficiency.py
# Brief:        Plots efficiency numerator, denominator and efficiency for ML and TBV vertex as a function
#               of MC run number.
# Usage:        python plotEfficiency.py rootfile mcLabel playlist
# Example Usage: python plotEfficiency.py Hists_Efficiency_ME6A_nosys_t12345_z00_AntiNu.root FullDet ME6A
# Author:       Anezka Klustova a.klustova20@imperial.ac.uk



import ROOT
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis, gPad, TLine
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
import os,sys
from ROOT import PlotUtils


rootfile = sys.argv[1]
mcLabel = sys.argv[2]  # FullDet or Nuke
print(str(mcLabel))
if str(mcLabel) == 'FullDet':
    print "Plotting Full Detector"
elif str(mcLabel) == "Nuke":
    print "Plotting Nuke"
else: 
    raise ValueError("Wrong MC Label: Choose FullDet or Nuke")

playlist = sys.argv[3]

infile = ROOT.TFile(str(rootfile))

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

ROOT.TH1.AddDirectory(False)
mcPOT = infile.Get("MCPOT").GetVal()

vars = ["mc_run"]
for var in vars:
    # numerator vs denominator
    num_hist_ML = infile.Get("selected_mc_MLreco_%s"%var) # selected signal
    num_hist_TBV = infile.Get("selected_mc_reco_%s"%var) # selected signal
    denom_hist = infile.Get("selected_truth_reco_%s"%var) # total signal

    numML = num_hist_ML.GetEntries()
    numTBV = num_hist_TBV.GetEntries()
    denom = denom_hist.GetEntries()

    num_hist_ML.SetNdivisions(505,"x")
    denom_hist.SetNdivisions(505,"x")
    num_hist_ML.SetNdivisions(505,"yz")
    denom_hist.SetNdivisions(505,"yz")

    print("Numerator ML: " + str(numML))
    print("Numerator TBV: " + str(numTBV))
    print("Denominator: " + str(denom))
    print("Efficiency ML: " + str(numML/denom))
    print("Efficiency TBV: " + str(numTBV/denom))

    num_hist_ML.SetMarkerColor(ROOT.kRed)
    num_hist_ML.SetMarkerStyle(ROOT.kFullCircle)
    num_hist_ML.SetMarkerSize(1)
    num_hist_TBV.SetMarkerColor(ROOT.kBlue)
    um_hist_TBV.SetMarkerStyle(ROOT.kFullCircle)
    num_hist_TBV.SetMarkerSize(1)

    num_hist_ML.GetXaxis().SetTitle("MC Run Number")
    num_hist_ML.GetYaxis().SetTitle("N Events")
    num_hist_ML.GetXaxis().CenterTitle()
    num_hist_ML.GetYaxis().CenterTitle()

    denom_hist.GetXaxis().SetTitle("MC Run Number")
    denom_hist.GetYaxis().SetTitle("N Events")
    denom_hist.GetXaxis().CenterTitle()
    denom_hist.GetYaxis().CenterTitle()

    denom_hist.SetMarkerStyle(ROOT.kFullCircle)
    num_hist_ML.SetMarkerSize(1)
    
    num_hist_ML.Draw('p')
    num_hist_ML.SetMaximum(num_hist_ML.GetMaximum()*1.5)
    num_hist_TBV.Draw('p SAME')

    legend = TLegend(0.55,0.7,0.80,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(num_hist_ML,  " Numerator ML"  , "p")
    legend.AddEntry(num_hist_TBV, " Numerator TBV", "p")
    legend.SetTextFont(42)
    legend.Draw()

    mnv.AddHistoTitle("Numerator %s %s"%(playlist,mcLabel), 0.05)
    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("%s_%s_Numerator_MLandTBV_%s.png"%(playlist, mcLabel, var))


    denom_hist.Draw('p')
    denom_hist.SetMaximum(denom_hist.GetMaximum()*1.5)
    mnv.AddHistoTitle("Denominator %s %s"%(playlist, mcLabel), 0.05)
    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("%s_%s_Denom_%s.png"%(playlist, mcLabel, var))

    num_hist_ML.Divide(num_hist_ML, denom_hist)
    num_hist_ML.Draw("p")
    num_hist_TBV.Divide(num_hist_TBV, denom_hist)
    num_hist_TBV.Draw('p SAME')
    num_hist_ML.SetMaximum(num_hist_ML.GetMaximum()*1.5)
    mnv.AddHistoTitle("Efficiency %s %s"%(playlist, mcLabel), 0.05)
    legend = TLegend(0.55,0.7,0.80,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(num_hist_ML,  " Efficiency ML"  , "p")
    legend.AddEntry(num_hist_TBV, " Efficiency TBV", "p")
    legend.SetTextFont(42)
    legend.Draw()
    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("%s_%s_Efficiency_MLandTBV_%s.png"%(playlist, mcLabel, var))



raw_input("Done")
