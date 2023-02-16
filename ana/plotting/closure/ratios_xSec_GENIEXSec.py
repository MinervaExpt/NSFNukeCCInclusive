import ROOT
import os,sys
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
from ROOT import PlotUtils

dirpwd = sys.argv[1]
targetID = sys.argv[2] 
targetZ = sys.argv[3]
plist = sys.argv[4]

mat = None
trueZ = None

if targetZ == "26":
    trueZ = "Iron"
    truez = "iron"
    mat = "Fe"

if targetZ == "82":
    trueZ = "Lead"
    truez = "lead"
    mat = "Pb"

if targetZ == "06":
    trueZ = "Carbon"
    truez = "carbon"
    mat = "C"

if targetZ == "99":
    trueZ = "Tracker"
    truez = "tracker"
    mat = "CH"

if targetZ == "99":
    xSec_file= ROOT.TFile(str(dirpwd)+"/%s/CrossSection_t%s_z%02s_%s.root"%(plist, targetID, targetZ, plist))
    GENIEXSec_file = ROOT.TFile(str(dirpwd)+"/closure_%s/GENIEXSecExtract_CCInclusive_t%s_z%s.root"%(plist, targetID, targetZ))
else:
    xSec_file= ROOT.TFile(str(dirpwd)+"%s/CrossSection_Daisy_t%s_z%02s_%s.root"%(plist, targetID, targetZ, plist))
    GENIEXSec_file = ROOT.TFile(str(dirpwd)+"/closure_%s/GENIEXSecExtract_CCInclusive_t%s_z%s.root"%(plist, targetID, targetZ))


canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()


ROOT.TH1.AddDirectory(False)

vars = ["Enu", 'x']

for var in vars:

    xSec = xSec_file.Get("simEventRate_crossSection_mc_%s"%(var))
    GENIEXSec = GENIEXSec_file.Get("%s_%s_%s_std_xsec"%(truez,targetID, var))

    xSec.GetCVHistoWithError()
    GENIEXSec.GetCVHistoWithError()

    ratio = xSec.Clone()
    ratio.Divide(xSec,GENIEXSec)#, 1.0, 1.0, "B") # binomial because numerator is the subset of the denominator, to make the stat error correct
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

    ratio.GetYaxis().SetTitle("xSec/GENIEXSec")
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

    if targetZ == "99":
        mnv.AddHistoTitle("%s (%s)"%(trueZ, plist), 0.05, 1)
    else:
        mnv.AddHistoTitle("Target %s %s (%s)"%(targetID, trueZ, plist), 0.04, 1)

    canvas1.Modified()
    canvas1.Update()
    if targetZ == "99":
        canvas1.Print("Closure_xSec_GENIExSec_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))
    else:
        canvas1.Print("Closure_xSec_GENIExSec_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))




print("DONE Closure test xSec/GENIExSec %s %s %02s"%(plist, targetID, targetZ))
