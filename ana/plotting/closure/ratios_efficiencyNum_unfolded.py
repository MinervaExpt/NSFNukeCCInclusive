import ROOT
import os,sys
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
from ROOT import PlotUtils

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
targetID = sys.argv[2] 
targetZ = sys.argv[3]
plist = sys.argv[4]
targetIDSingle = sys.argv[5]

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

if targetZ == "99":
    infile= ROOT.TFile(str(dirpwd)+"/CrossSection_t%s_z%02s_%s.root"%(targetID, targetZ, plist))
else:
    infile= ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t%s_z%02s_%s.root"%(targetID, targetZ, plist))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()


ROOT.TH1.AddDirectory(False)

vars = ["Enu", "x"]

for var in vars:

    # numerator vs unfolded
    eff_num = infile.Get("efficiency_numerator_%s_%s"%(targetIDSingle,var)) # selected signal
    unfolded = infile.Get("unfolded_%s_mc_%s"%(targetIDSingle, var)) 

    ratio = eff_num.Clone()
    ratio.Divide(eff_num,unfolded)#, 1.0, 1.0, "B") # binomial because numerator is the subset of the denominator, to make the stat error correct
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

    ratio.GetYaxis().SetTitle("Eff. Num/Unfolded MC")
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
    ratio.SetMinimum(0.9)

    gStyle.SetOptTitle(0)

    if targetZ == "99":
        mnv.AddHistoTitle("%s (%s)"%(trueZ, plist), 0.05, 1)
    else:
        mnv.AddHistoTitle("Target %s %s (%s)"%(targetIDSingle, trueZ, plist), 0.04, 1)

    canvas1.Modified()
    canvas1.Update()
    if targetZ == "99":
        canvas1.Print("Closure_efficiencyNum_unfolded_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))
    else:
        canvas1.Print("Closure_Daisy_efficiencyNum_unfolded_t%s_z%02s_%s_%s.png"%(targetIDSingle, targetZ, var, plist))

print("DONE Closure test efficiencyNum/unfolded %s %s %02s"%(plist, targetID, targetZ))
