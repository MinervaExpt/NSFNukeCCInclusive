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

path = "/minerva/data2/users/anezkak/flux4Daisy_files_2022/OptimisedEventLoopME1ABC_flux_grid_sys/"
tracker_file = ROOT.TFile(str(path)+"flux_tracker.root")
iron_file = ROOT.TFile(str(path)+"flux_iron.root")
lead_file = ROOT.TFile(str(path)+"flux_lead.root")
carbon_file = ROOT.TFile(str(path)+"flux_carbon.root")
#water_file = ROOT.TFile(str(path)+"flux_water_apothem.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

tracker = tracker_file.Get("flux")
iron = iron_file.Get("flux")
lead = lead_file.Get("flux")
carbon = carbon_file.Get("flux")
#water = water_file.Get("flux")

iron.SetLineColor(ROOT.kBlue)
lead.SetLineColor(ROOT.kGreen)
carbon.SetLineColor(ROOT.kRed)
numu.SetLineColor(ROOT.kViolet)
#water.SetLineColor(ROOT.kCyan)

tracker.Draw("HIST")
iron.Draw("HIST SAME")
lead.Draw("HIST SAME")
carbon.Draw("HIST SAME")
#water.Draw("HIST SAME")
numu.Draw("HIST SAME")


mnv.AddHistoTitle("Flux ME1A, ME1B, ME1C", 0.05, 1)
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
leg.AddEntry(numu,"FRW Tracker","l")
leg.AddEntry(iron,"Iron","l")
leg.AddEntry(lead,"Lead","l")
leg.AddEntry(carbon,"Carbon","l")
#leg.AddEntry(water,"Water","l")
leg.SetTextFont(42)
leg.Draw()

canvas1.Modified()
canvas1.Print("ME1ABC_Flux_FromGrid_sys.png")

'''
# ----------------------------------------------------------------------------
# Data/MC ratio
print(type(tracker))
print(type(numu))
print(tracker.GetNbinsX())
print(numu.GetNbinsX())

# Cannot increase the size of the bins using rebin
# Rebin both to 100 for check

bins_nuE = np.linspace(0, 100, 101)
nBins_nuE = len(bins_nuE)-1
referenceHist = MnvH1D( 'h_flux_reference' , 'h_flux_reference' , nBins_nuE , array('d',bins_nuE))

rebinned = tracker.Rebin(10)
rebinned.Scale(rebinned.GetNormBinWidth(), "width") # bin width normalise
rebinnedFlux = myfrw.GetRebinnedFluxGenerated(14,referenceHist)

print(rebinned.GetNbinsX())
print(type(rebinned))
print(rebinnedFlux.GetNbinsX())
print(type(rebinnedFlux))

rebinned.Divide(rebinned, rebinned)
rebinned.Draw("HIST")

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
'''
raw_input("Done")