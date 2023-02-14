import ROOT
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis, gPad, TLine
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
import os,sys
from ROOT import PlotUtils

dirpwd = sys.argv[1]
plist = sys.argv[2]


infile= ROOT.TFile(str(dirpwd)+"/Efficiency_Daisy_%s_t99_z99_sys.root"%(plist))

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()
mnv.SetROOT6Palette(87) #kLightTenperature

ROOT.TH1.AddDirectory(False)
mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale =  dataPOT/mcPOT


trueZ = "Tracker"
mat = "CH"

var = "pZmu_pTmu"
for petal in range(0,12):

    # numerator vs denominator
    num_hist = infile.Get("h_mc2d_%s_pZmu_pTmu"%petal) # selected signal
    denom_hist = infile.Get("h_truth2d_%s_pZmu_pTmu"%petal) # total signal

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

    mnv.AddHistoTitle("Daisy %s: Petal %s"%(trueZ, petal), 0.05, 1)


    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("Efficiency2D_t99_z99_%s_%s_p%s.png"%(plist, var, petal))

print("DONE %s 99 99 Daisy"%(plist))
