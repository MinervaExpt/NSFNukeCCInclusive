import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLine
from ROOT import TLegend
from ROOT.PlotUtils import MnvH1D
from array import array
import numpy as np


myfrw = PlotUtils.FluxReweighter(14,True,PlotUtils.FluxReweighter.minervame1A,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6, 500)
numu = myfrw.GetFluxReweighted(14)

tracker_file = ROOT.TFile("flux_tracker.root")

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

tracker = tracker_file.Get("flux")
rebinnedFlux = myfrw.GetRebinnedFluxGenerated(14,tracker)
rebinnedFlux.Scale(rebinnedFlux.GetBinWidth(1), "width")

rebinnedFlux.SetLineColor(ROOT.kViolet)


tracker.Draw("HIST")

rebinnedFlux.Draw("HIST SAME")


mnv.AddHistoTitle("Flux ME1A", 0.05, 1)
gStyle.SetErrorX(0)

#fix title
tracker.GetYaxis().CenterTitle(True)
tracker.GetYaxis().SetTitleOffset(1.3)
tracker.GetYaxis().SetTitleSize(0.05)
tracker.GetYaxis().SetLabelSize(0.05)
tracker.GetXaxis().CenterTitle(True)
tracker.GetXaxis().SetTitle("Neutrino Energy (GeV)")
#tracker.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
tracker.GetXaxis().SetTitleSize(0.05)
tracker.GetXaxis().SetLabelSize(0.05)
tracker.GetYaxis().SetTitle("#nu / m^{2} / P.O.T/ GeV")
#tracker.GetYaxis().SetTitle("#bar{#nu} / m^{2} / P.O.T/ GeV")

leg = ROOT.TLegend(0.60,0.65,0.85,0.9)
leg.SetTextSize(0.035)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetLineColor(0)
leg.AddEntry(tracker,"Tracker","l")
leg.AddEntry(rebinnedFlux,"FRW Tracker","l")

#leg.AddEntry(water,"Water","l")
leg.SetTextFont(42)
leg.Draw()

canvas1.Modified()
canvas1.Print("ME1A_Flux_FromGrid_sys.png")


# ----------------------------------------------------------------------------
# Data/MC ratio
print(type(tracker))
print(type(numu))
print(tracker.GetNbinsX())
print(numu.GetNbinsX())

# Cannot increase the size of the bins using rebin
# Rebin both to 100 for check


tracker.Divide(tracker, rebinnedFlux)
tracker.Draw("HIST")

#iron.Divide(iron,denom)
#lead.Divide(lead,denom)
#carbon.Divide(carbon,denom)
#tracker.Draw("HIST")
#iron.Draw("HIST SAME")
#lead.Draw("HIST SAME")
#carbon.Draw("HIST SAME")

tracker.SetMinimum(0)
tracker.SetMaximum(2)

mnv.AddHistoTitle("Flux ME1A Tracker", 0.05, 1)
gStyle.SetErrorX(0)

canvas1.Modified()
canvas1.Print("ME1A_TrackerFlux_ratio.png")

raw_input("Done")