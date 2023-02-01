import ROOT
import os,sys
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
from ROOT import PlotUtils

#targetID = sys.argv[1] 
#targetZ = sys.argv[2]

efficiency = ROOT.TFile("Hists_EfficiencyTracker_ML_ME6A_sys_t99_z99_AntiNu.root")
evRate_file = ROOT.TFile("GENIEXSecExtract_CCInclusive_Tracker.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()


ROOT.TH1.AddDirectory(False)

vars = ["Enu", 'x']

for var in vars:

    eff_denom = efficiency.Get("h_truth_%s"%var) # total signal
    evRate = evRate_file.Get("dif_%s_xsec_evRate"%var) 

    eff_denom.GetCVHistoWithError()
    evRate.GetCVHistoWithError()

    ratio = eff_denom.Clone()
    ratio.Divide(eff_denom,evRate)#, 1.0, 1.0, "B") # binomial because numerator is the subset of the denominator, to make the stat error correct
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
    ratio.GetXaxis().SetTitleSize(0.05)
    ratio.GetYaxis().SetTitleSize(0.05)
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetXaxis().SetTitleSize(0.05)
    ratio.GetYaxis().SetTitleSize(0.05)
    ratio.GetYaxis().SetTitleOffset(1.6)

    # ----------------------------------------------------------------------------
    # Ratio

    ratio.GetYaxis().SetTitle("Eff Denom/GXSec Event Rate")
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

    mnv.AddHistoTitle("ME6A Tracker", 0.05, 1)
    mnv.AddPlotLabel("Integral Numerator: "+ str(eff_denom.Integral()), 0.40, 0.87, 0.033, 12, 42)
    mnv.AddPlotLabel("Integral Denominator: "+ str(evRate.Integral()), 0.414, 0.82, 0.033, 12, 42)

    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("ME6A_Closure_EffDenom_GENIEXSecEvRate_ratio_%s.png"%var)



raw_input("Done")
