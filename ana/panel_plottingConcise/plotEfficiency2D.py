import ROOT
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis, gPad, TLine
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
import os,sys
from ROOT import PlotUtils

dirpwd = sys.argv[1]
targetID = sys.argv[2] 
targetZ = sys.argv[3]
plist = sys.argv[4]

infile= ROOT.TFile(str(dirpwd)+"/Efficiency_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()
mnv.SetROOT6Palette(87) #kLightTenperature

ROOT.TH1.AddDirectory(False)
mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale =  dataPOT/mcPOT

mat = None
trueZ = None

if targetZ == "26":
    trueZ = "Iron"
    mat = "Fe"

if targetZ == "82":
    trueZ = "Lead"
    mat = "Pb"

if targetZ == "06":
    trueZ = "Carbon"
    mat = "C"

if targetZ == "99":
    trueZ = "Tracker"
    mat = "CH"

# numerator vs denominator
num_hist = infile.Get("h_mc_pZmu_pTmu") # selected signal
denom_hist = infile.Get("h_truth_pZmu_pTmu") # total signal

num = num_hist.GetEntries()
denom = denom_hist.GetEntries()
print("Numerator: " + str(num))
print("Denominator: " + str(denom))
print("Efficiency: " + str(num/denom))

ratio = num_hist.Clone()
ratio.Divide(num_hist,denom_hist, 1.0, 1.0, "B") # binomial because numerator is the subset of the denominator, to make the stat error correct

ratio.GetXaxis().SetTitle("Muon p_{||} (GeV/c)")
ratio.GetYaxis().SetTitle("Muon p_{t} (GeV/c)")

ratio.GetZaxis().SetTitle("Efficiency")
ratio.GetXaxis().CenterTitle()
ratio.GetYaxis().CenterTitle()
ratio.GetZaxis().CenterTitle()

ratio.GetXaxis().SetTitleFont(42)
ratio.GetXaxis().SetTitleSize(0.05)

ratio.GetYaxis().SetTitleOffset(1)
ratio.GetYaxis().SetTitleFont(42)
ratio.GetYaxis().SetTitleSize(0.05)

ratio.GetZaxis().SetTitleFont(42)
ratio.GetZaxis().SetTitleOffset(1.1)
ratio.GetZaxis().SetTitleSize(0.045)

ratio.Draw("COLZ")

ratio.SetMaximum(1)

gStyle.SetOptTitle(0)

legend = TLegend(0.56,0.70,0.80,0.89)
legend.SetBorderSize(0)
legend.SetTextSize(0.035)

if targetZ == "99":
     mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
else:
    mnv.AddHistoTitle("Target %s %s"%(targetID, trueZ), 0.05, 1)

canvas1.Modified()
canvas1.Update()
canvas1.Print("Efficiency2D_t%s_z%02s_%s.png"%(targetID, targetZ, plist))

print("DONE %s %s %02s"%(plist, targetID, targetZ))
