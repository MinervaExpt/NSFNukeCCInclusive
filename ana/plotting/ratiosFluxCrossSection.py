import ROOT
import os,sys
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
from ROOT import PlotUtils

ROOT.gROOT.SetBatch(True)

target = ROOT.TFile("/minerva/data2/users/anezkak/22-06-22-xSec_ratios_daisy/CrossSection_Daisy_t3_z06_minervame6A.root")
tracker = ROOT.TFile("/minerva/data2/users/anezkak/22-06-22-xSec_ratios_daisy/CrossSection_Daisy_t99_z99_minervame6A.root")
material = 'carbon'

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

ROOT.TH1.AddDirectory(False)

vars = ["Enu"]

for var in vars:
    #if var == "Enu":
    #target_data_hist = target.Get("crossSection_total_data_%s"%var)
    #target_simEventRate_hist = target.Get("simEventRate_crossSection_total_mc_%s"%var) # total signal
    #tracker_data_hist = tracker.Get("carbon_3_%s_std_xsec"%var)
    #tracker_simEventRate_hist = tracker.Get("carbon_3_%s_std_xsec"%var) 
    
    target_flux = target.Get("fluxRebinned_%s"%var)
    tracker_flux = tracker.Get("fluxRebinned_%s_%s"%(material,var))

    ratio = target_flux.Clone()
    ratio.Divide(target_flux,tracker_flux)

    if var == "Enu":
        ratio.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
        ratio.GetYaxis().SetTitle("Events/GeV")


    ratio.GetXaxis().CenterTitle()
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetXaxis().SetTitleSize(0.05)
    ratio.GetYaxis().SetTitleSize(0.05)
    ratio.GetYaxis().SetTitleOffset(1.6)

    # ----------------------------------------------------------------------------
    # Ratio

    ratio.GetYaxis().SetTitle("Flux %s/tracker"%material)
    ratio.GetXaxis().CenterTitle()
    ratio.GetYaxis().CenterTitle()
    ratio.Draw("HIST")

    ratio.Draw("HIST")

    gStyle.SetOptTitle(0)

    mnv.AddHistoTitle("ME6A Flux %s/tracker"%material, 0.05, 1)

    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("ME6A_Flux_xsec_ratio_%s_%s.png"%(material,var))


raw_input("Done")
