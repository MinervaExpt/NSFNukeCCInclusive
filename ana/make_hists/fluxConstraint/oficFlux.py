import ROOT
import os,sys
from ROOT import PlotUtils

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

myfrw = PlotUtils.FluxReweighter(14,True,PlotUtils.FluxReweighter.minervame1A,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6, 500)
myfrw2 = PlotUtils.FluxReweighter(-14,True,PlotUtils.FluxReweighter.minervame6A,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6, 500)

numu = myfrw.GetFluxReweighted(14)
numubar = myfrw2.GetFluxReweighted(-14)

#fix title
numu.GetYaxis().CenterTitle(True)
numu.GetYaxis().SetTitleOffset(1.3)
numu.GetYaxis().SetTitleSize(0.05)
numu.GetYaxis().SetLabelSize(0.05)
numu.GetXaxis().CenterTitle(True)
numu.GetXaxis().SetTitle("Neutrino Energy (GeV)")
numu.GetXaxis().SetTitleSize(0.05)
numu.GetXaxis().SetLabelSize(0.05)
numu.GetYaxis().SetTitle("#nu / m^{2} / P.O.T/ GeV")

leg = ROOT.TLegend(0.55,0.7,0.7,0.8)
leg.SetTextSize(0.035)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.AddEntry(numu,"Neutrino ME1A","l")
leg.AddEntry(numubar,"Antineutrino ME6A","l")

numu.Draw("HIST")
numubar.SetLineColor(ROOT.kRed)
numubar.Draw("HIST SAME")
leg.Draw("SAME")
mnv.AddHistoTitle("Tracker Flux from FluxReweighter", 0.04, 1)

#canvas1.SetLogy(True)
numu.GetXaxis().SetRangeUser(0,20)

#numu.GetYaxis().SetRangeUser(1e-8,1e-3)

canvas1.Modified()
canvas1.Print("ME1AandME6A_oficFlux.png")