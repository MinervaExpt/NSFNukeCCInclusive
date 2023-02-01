import ROOT
import os,sys
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
from ROOT import PlotUtils

#targetID = sys.argv[1] 
#targetZ = sys.argv[2]

bkg = ROOT.TFile("Hists_BkgSubtracted_EventSelection_sys_t3_z26_AntiNu_TRUE.root")
unfold = ROOT.TFile("CrossSection_t3_z26_minervame6A.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()


ROOT.TH1.AddDirectory(False)

vars = ["Enu", "x"]

for var in vars:

    # numerator vs denominator
    bkg_sub = bkg.Get("h_bkg_subtracted_mc_%s"%var) # selected signal
    unfolded = unfold.Get("unfolded_mc_%s"%var) 

    ratio = bkg_sub.Clone()
    ratio.Divide(bkg_sub,unfolded)#, 1.0, 1.0, "B") # binomial because numerator is the subset of the denominator, to make the stat error correct
    sysError = ratio.Clone()
    ratio.SetLineColor(ROOT.kRed)

    if var == "Enu":
        ratio.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
        ratio.GetYaxis().SetTitle("Events/GeV")

    if var == "x":
        ratio.GetXaxis().SetTitle("Bjorken x")
        ratio.GetYaxis().SetTitle("Events (norm.)")

    ratio.GetXaxis().CenterTitle()
    ratio.GetYaxis().CenterTitle()

    # ----------------------------------------------------------------------------
    # Ratio

    ratio.GetYaxis().SetTitle("BkgSubMC True/Unfolded MC")
    ratio.GetXaxis().CenterTitle()
    ratio.GetYaxis().CenterTitle()
    ratio.Draw("HIST")

    '''
    sysError.SetFillColor(ROOT.kRed-10)
    sysError.SetFillStyle(1001)
    sysError.SetMarkerStyle(0)
    sysError.Draw("E2 SAME")
    sysError.SetMaximum(1)
    '''
    ratio.Draw("HIST SAME")
    ratio.SetMaximum(1.01)
    ratio.SetMinimum(0.99)

    gStyle.SetOptTitle(0)

    mnv.AddHistoTitle("ME6A Target 3 Iron", 0.05, 1)

    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("ME6A_Closure_BkgSubtractedTrue_UnfoldedMC_ratio_%s.png"%var)


raw_input("Done")
