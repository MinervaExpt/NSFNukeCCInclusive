import ROOT
import os,sys
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
from ROOT import PlotUtils

newPWD = sys.argv[1] 
origPWD = sys.argv[2] 
targetID = sys.argv[3] 
targetZ = sys.argv[4]
plist = sys.argv[5]
warp = sys.argv[6]

new = ROOT.TFile(str(newPWD)+"/EventSelection/EventSelection_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))
orig = ROOT.TFile(str(origPWD)+"/EventSelection/EventSelection_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

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

ROOT.TH1.AddDirectory(False)

vars = ["Enu", "x", "pTmu1D", "pZmu1D", "ThetamuDeg"]

for var in vars:

    # numerator vs denominator
    new_hist= new.Get("selected_mc_reco_%s"%var) # selected signal
    orig_hist = orig.Get("selected_mc_reco_%s"%var) 

    new_hist_stat = new_hist.GetCVHistoWithStatError()
    orig_hist_stat = orig_hist.GetCVHistoWithStatError()
    ratio = new_hist_stat.Clone()
    ratio.Divide(ratio,orig_hist_stat)#, 1.0, 1.0, "B") # binomial because numerator is the subset of the denominator, to make the stat error correct
    sysError = ratio.Clone()
    ratio.SetLineColor(ROOT.kRed)

    if var == "Enu":
        ratio.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")
        ratio.GetYaxis().SetTitle("Events/GeV")

    if var == "x":
        ratio.GetXaxis().SetTitle("Reconstructed Bjorken x")
        ratio.GetYaxis().SetTitle("Events (norm.)")
    
    if var == "ThetamuDeg":
        ratio.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)")
        ratio.GetYaxis().SetTitle("Events/Deg")
    
    if var == "pTmu1D":
        ratio.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
        ratio.GetYaxis().SetTitle("Events/(GeV/c)")

    if var == "pZmu1D":
        ratio.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
        ratio.GetYaxis().SetTitle("Events/(GeV/c)")

    ratio.GetXaxis().CenterTitle()
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetXaxis().SetTitleSize(0.05)
    ratio.GetYaxis().SetTitleSize(0.05)
    ratio.GetYaxis().SetTitleOffset(1.2)

    # ----------------------------------------------------------------------------
    # Ratio

    ratio.GetYaxis().SetTitle("Warped/CV Ratio")
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
    ratio.SetMaximum(1.5)
    ratio.SetMinimum(0.5)

    gStyle.SetOptTitle(0)

    if targetZ == "99":
        mnv.AddHistoTitle("%s Warp: %s"%(warp, trueZ), 0.05, 1)
    else:
        mnv.AddHistoTitle("%s Warp: Target %s %s"%(warp, targetID, trueZ), 0.05, 1)
    mnv.AddPlotLabel("Integral Numerator: "+ str(new_hist.Integral()), 0.40, 0.87, 0.033, 12, 42)
    mnv.AddPlotLabel("Integral Denominator: "+ str(orig_hist.Integral()), 0.414, 0.82, 0.033, 12, 42)

    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("%sWarpCVRatio_t%s_z%02s_%s_%s.png"%(warp, targetID, targetZ, var, plist))


print("DONE %s %s %s %02s"%(warp, plist, targetID, targetZ))