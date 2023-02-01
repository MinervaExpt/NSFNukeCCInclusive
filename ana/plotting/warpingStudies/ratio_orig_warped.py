import ROOT
import os,sys
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
from ROOT import PlotUtils

#targetID = sys.argv[1] 
#targetZ = sys.argv[2]

new = ROOT.TFile("Hists_EventSelectionTracker_ML_ME6A_nosys_DeuteriumGeniePiTune_t99_z99_AntiNu.root")
orig = ROOT.TFile("Hists_EventSelectionTracker_ML_ME6A_nosys_t99_z99_AntiNu.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()


ROOT.TH1.AddDirectory(False)

vars = ["Enu", "x"]

for var in vars:

    # numerator vs denominator
    new_hist= new.Get("selected_mc_reco_%s"%var) # selected signal
    orig_hist = orig.Get("selected_mc_reco_%s"%var) 

    ratio = new_hist.Clone()
    ratio.Divide(new_hist,orig_hist)#, 1.0, 1.0, "B") # binomial because numerator is the subset of the denominator, to make the stat error correct
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
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetXaxis().SetTitleSize(0.05)
    ratio.GetYaxis().SetTitleSize(0.05)
    ratio.GetYaxis().SetTitleOffset(1.2)

    # ----------------------------------------------------------------------------
    # Ratio

    ratio.GetYaxis().SetTitle("Warped/Orig Ratio ")
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
    ratio.SetMaximum(1.2)
    ratio.SetMinimum(0.8)

    gStyle.SetOptTitle(0)

    mnv.AddHistoTitle("ME6A Tracker", 0.05, 1)
    mnv.AddPlotLabel("Integral Numerator: "+ str(new_hist.Integral()), 0.40, 0.87, 0.033, 12, 42)
    mnv.AddPlotLabel("Integral Denominator: "+ str(orig_hist.Integral()), 0.414, 0.82, 0.033, 12, 42)

    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("ME6A_Warped_Orig_ratio_%s.png"%var)


raw_input("Done")
