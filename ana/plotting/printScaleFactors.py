import ROOT
import os,sys
from ROOT import PlotUtils
import math as m

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
outdir = sys.argv[2]
targetID = sys.argv[3] 
targetZ = sys.argv[4]
plist = sys.argv[5]

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

infile= ROOT.TFile(str(outdir)+"/Plastic_ScaleFactors_t%s_z%s_%s.root"%(targetID, targetZ,plist))
infileUntuned = ROOT.TFile.Open(str(dirpwd)+"/PlasticBkg_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

vars = ["planeDNN", "Enu", "x", "pTmu", "pZmu"]

for var in vars:
    upstream = infile.Get("scaleFactor_US_%s"%var)
    downstream = infile.Get("scaleFactor_DS_%s"%var)

    hists_US_regUS = infileUntuned.Get("US_regUS_%s_%s"%(trueZ,var))
    hists_DS_regDS = infileUntuned.Get("DS_regDS_%s_%s"%(trueZ, var))

    integralUS = hists_US_regUS.Integral()
    statErrUS = 1/ m.sqrt(integralUS)
    for bin in range(upstream.GetNbinsX() + 1):
        upstream.SetBinError(bin,statErrUS)

    # Downstream: Remove original stat error and update it with 1/sqrt(N) in each bin
    integralDS = hists_DS_regDS.Integral()
    statErrDS = 1/ m.sqrt(integralDS)
    for bin in range(downstream.GetNbinsX() + 1):
        downstream.SetBinError(bin,statErrDS)	

    bin = 1
    US_scale_factor = upstream.GetBinContent(bin)
    US_statErr = upstream.GetCVHistoWithStatError().GetBinError(bin)
    US_sysErr = upstream.GetCVHistoWithError(False).GetBinError(bin)

    DS_scale_factor = downstream.GetBinContent(bin)
    DS_statErr = downstream.GetCVHistoWithStatError().GetBinError(bin)
    DS_sysErr = downstream.GetCVHistoWithError(False).GetBinError(bin)


    print("Upstream (US):")
    print("Scale Factor = " + str(US_scale_factor))
    print("Statistical error = " + str(US_statErr))
    print("Systematic error = " + str(US_sysErr))

    print("-------------------------------------------------------------------------")
    print("Downstream (DS):")
    print("Scale Factor = " + str(DS_scale_factor))
    print("Statistical error = " + str(DS_statErr))
    print("Systematic error = " + str(DS_sysErr))

    # ---------------------------------------------------------------------------------
    # PLOTTING

    upstreamStat = upstream.GetCVHistoWithStatError().Clone()
    upstreamStat.SetLineColor(ROOT.kBlack)
    upstreamStat.Draw("X0")
    upstreamStat.SetMaximum(1.5)
    upstreamStat.SetMinimum(0.5)
    upstreamStat.SetMarkerSize(0.75)
    upstreamStat.SetMarkerStyle(1)
    sysError = upstream.GetCVHistoWithError(False).Clone() # systematic only
    sysError.SetFillColor(ROOT.kRed-10)
    sysError.SetFillStyle(1001)
    sysError.SetMarkerStyle(0)
    sysError.Draw("E2 SAME")

    if var == "Enu":
        upstreamStat.GetXaxis().SetTitle("Reconstructed Neutrino Energy (GeV)")

    if var == "x":
        upstreamStat.GetXaxis().SetTitle("Reconstructed Bjorken x")

    if var == "planeDNN":
        upstreamStat.GetXaxis().SetTitle("Plane Number")
    
    if var == "pTmu":
        upstreamStat.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
       
    if var == "pZmu":
        upstreamStat.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")

    upstreamStat.GetYaxis().SetTitle("Scaling Factor")
    upstreamStat.GetXaxis().CenterTitle()
    upstreamStat.GetYaxis().CenterTitle()

    upstreamStat.Draw("X0 SAME") # with stats
    mnv.AddHistoTitle("Target %s %s US Scale Factor"%(targetID, trueZ), 0.05, 1)


    canvas1.Modified()
    canvas1.Print("ScaleFactor_US_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))



    # -----------------------------------------------------------------------
    canvas2 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
    canvas2.cd()

    downstreamStat = downstream.GetCVHistoWithStatError().Clone()
    downstreamStat.SetLineColor(ROOT.kBlack)
    downstreamStat.Draw("X0")
    downstreamStat.SetMaximum(1.5)
    downstreamStat.SetMinimum(0.5)
    downstreamStat.SetMarkerSize(0.75)
    downstreamStat.SetMarkerStyle(1)
    sysError = downstream.GetCVHistoWithError(False).Clone() # systematic only
    sysError.SetFillColor(ROOT.kRed-10)
    sysError.SetFillStyle(1001)
    sysError.SetMarkerStyle(0)
    sysError.Draw("E2 SAME")

    if var == "Enu":
        downstreamStat.GetXaxis().SetTitle("Reconstructed Neutrino Energy (GeV)")

    if var == "x":
        downstreamStat.GetXaxis().SetTitle("Reconstructed Bjorken x")

    if var == "planeDNN":
        downstreamStat.GetXaxis().SetTitle("Plane Number")

    downstreamStat.GetYaxis().SetTitle("Scaling Factor")
    downstreamStat.GetXaxis().CenterTitle()
    downstreamStat.GetYaxis().CenterTitle()

    downstreamStat.Draw("X0 SAME") # with stats
    mnv.AddHistoTitle("Target %s %s DS Scale Factor"%(targetID, trueZ), 0.05, 1)

    canvas2.Modified()
    canvas2.Print("ScaleFactor_DS_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))

print("DONE %s %s %02s"%(plist, targetID, targetZ))

