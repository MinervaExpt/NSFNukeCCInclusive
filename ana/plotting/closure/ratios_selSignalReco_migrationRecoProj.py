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


signal = ROOT.TFile(str(dirpwd)+"/BackgroundSubtracted/BkgSubtracted_EventSelection_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))
mig = ROOT.TFile(str(dirpwd)+"/Migration/Migration_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))

print("BkgSubtracted_EventSelection_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))
print("Migration_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()


ROOT.TH1.AddDirectory(False)

vars = ["Enu", "x"]

for var in vars:

    # numerator vs denominator
    selsig = signal.Get("h_bkg_subtracted_mc_%s"%var) 
    migra = mig.Get("selected_Migration_%s"%var)
    migraX = migra.ProjectionX()

    print(type(migraX))
    
    ratio = selsig.Clone()
    ratio.Divide(selsig, migraX)#, 1.0, 1.0, "B") # binomial because numerator is the subset of the denominator, to make the stat error correct
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
    
    ratio.GetYaxis().SetTitle("Selected signal/Mig RecoProj")
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
        mnv.AddHistoTitle("Target %s %s (%s)"%(targetID, trueZ, plist), 0.04, 1)

    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("Closure_selSig_migRecoProj_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))


print("DONE Closure test selectedSignal/migrationRecoProjection %s %s %02s"%(plist, targetID, targetZ))

