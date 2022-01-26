import ROOT
import os,sys
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
from ROOT import PlotUtils

#targetID = sys.argv[1] 
#targetZ = sys.argv[2]

optim_file = ROOT.TFile("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/make_hists/fluxConstraint/Hists_EventSelection_minervame1A_FluxConstraint_optim_sys_t99_z99_Nu.root")
me_file = ROOT.TFile("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/make_hists/fluxConstraint/Hists_EventSelection_minervame1A_FluxConstraint_sys_t99_z99_Nu.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()


ROOT.TH1.AddDirectory(False)

# numerator vs denominator
optim = optim_file.Get("selected_mc_truth_Enu").GetCVHistoWithError() # selected signal
me = me_file.Get("selected_mc_truth_Enu").GetCVHistoWithError() # total selected events

ratio = optim.Clone()
ratio.Divide(me,optim)#, 1.0, 1.0, "B") # binomial because numerator is the subset of the denominator, to make the stat error correct
sysError = ratio.Clone()
ratio.SetLineColor(ROOT.kRed)


ratio.GetXaxis().SetTitle("Neutrino Energy (GeV)")

ratio.GetXaxis().CenterTitle()
ratio.GetYaxis().CenterTitle()

# ----------------------------------------------------------------------------
# Ratio

ratio.GetYaxis().SetTitle("optim/orig")
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
ratio.SetMaximum(1.1)
ratio.SetMinimum(0.8)

gStyle.SetOptTitle(0)

mnv.AddHistoTitle("ME1A Tracker Flux optim/orig", 0.05, 1)

canvas1.Modified()
canvas1.Update()
canvas1.Print("ME1A_ratio_optim_old_tracker.png")


raw_input("Done")
