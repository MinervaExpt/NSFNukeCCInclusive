import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle


infile= ROOT.TFile("GENIEXSecExtract_CCInclusive_Tracker.root")
infile2 = ROOT.TFile("CrossSection_t99_z99_minervame6A2.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

flux = infile.Get("hFlux")
flux2 = infile2.Get('flux_Enu')

flux2.SetLineColor(ROOT.kRed)

#if var == "Enu":
#    mc_hist.GetXaxis().SetTitle("Reconstructed Neutrino Energy (GeV)")
#    mc_hist.GetYaxis().SetTitle("Events/GeV")


#mc_hist.GetXaxis().CenterTitle()
#mc_hist.GetYaxis().CenterTitle()
#print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
#print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
#mc_hist.GetXaxis().SetRangeUser(66,173)
flux.Scale(1.33E-46)

flux.Draw("HIST ")
flux2.Draw("HIST SAME ")
#flux.Divide(flux, flux2)
#flux.Draw("HIST ")
#canvas1.SetLogy()
print(flux.Integral())
# = 3.59945441751e+43
print(flux2.Integral())
canvas1.Modified()
canvas1.Print("GENIEXSecExtract_Flux.png")


flux.SetDirectory(0)

raw_input("Done")
