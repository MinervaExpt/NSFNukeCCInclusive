import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLine

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
targetID = "99" #sys.argv[2] 
targetZ = "99" # sys.argv[3]
plist = "minervame5A6A6B6C6D6E6F6G6H6I6J" #sys.argv[4]
scale = 1 #sys.argv[5]

infile= ROOT.TFile(str(dirpwd)+"/BkgSubtracted_EventSelection_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale = dataPOT/mcPOT
if scale == "1":
    mcScale = 1

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

vars = ["x"]
#["Enu", "x", "pZmu1D", "pTmu1D", "ThetamuDeg"]

for var in vars:

    mc_hist = infile.Get("h_bkg_subtracted_mc_%s"%var)
    data_hist = infile.Get("h_bkg_subtracted_data_%s"%var)
    data_hist_notConstrained = None
    if targetZ != "99":
        data_hist_notConstrained = infile.Get("h_bkg_subtracted_data_notConstrained_%s"%var)

    if var == "Enu":
        mc_hist.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")
        mc_hist.GetYaxis().SetTitle("Events/GeV")
        data_hist.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")
        if targetZ != "99":
            data_hist_notConstrained.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")

    if var == "x":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
        data_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        if targetZ != "99":
            data_hist_notConstrained.GetXaxis().SetTitle("Reconstructed Bjorken x") 

    if var == "pTmu1D":
        mc_hist.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
        mc_hist.GetYaxis().SetTitle("Events/(GeV/c)")
        data_hist.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
        if targetZ != "99":
            data_hist_notConstrained.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)") 

    if var == "pZmu1D":
        mc_hist.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
        mc_hist.GetYaxis().SetTitle("Events/(GeV/c)")
        data_hist.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/)")
        if targetZ != "99":
            data_hist_notConstrained.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)") 
    
    if var == "ThetamuDeg":
        mc_hist.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)")
        mc_hist.GetYaxis().SetTitle("Events/Deg")
        data_hist.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)")
        if targetZ != "99":
            data_hist_notConstrained.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)") 


    '''
    if myvariable == "xfine":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")


    if myvariable == "xBrian":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")

    if myvariable == "x09":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
    ''' 
    mc_hist.GetXaxis().CenterTitle()
    mc_hist.GetYaxis().CenterTitle()
    #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
    #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
    #mc_hist.GetXaxis().SetRangeUser(0,0.001)

    data_hist_stat =  data_hist.GetCVHistoWithStatError() # stat error
    data_hist_total = data_hist.GetCVHistoWithError() # total error
    data_hist_sys = data_hist.GetCVHistoWithError(False) # sys error (bool is include stat)
    mc_hist_stat = mc_hist.GetCVHistoWithStatError() 
    
    # MC
    mc_hist.SetLineWidth(3)
    mc_hist.SetLineColor(2)
    # MC error
    mc_hist_stat.SetFillColor(ROOT.kRed-10)
    mc_hist_stat.SetFillStyle(1001)
    mc_hist_stat.SetMarkerStyle(0)
    mc_hist_stat.SetLineWidth(3)
    mc_hist_stat.SetLineColor(2)

    mc_hist.SetLineColor(ROOT.kRed)
    mc_hist.SetLineWidth(2)

    
    ratio = data_hist_stat.Clone()
    ratio.Divide(ratio,mc_hist_stat) # stat
    if var == "Enu":
        ratio.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")
    if var == "x":
        ratio.GetXaxis().SetTitle("Reconstructed Bjorken x")
    if var == "pTmu1D":
        ratio.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
    if var == "pZmu1D":
        ratio.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
    if var == "ThetamuDeg":
        ratio.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)")

    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetYaxis().CenterTitle()
    ratio.GetXaxis().CenterTitle()
    ratio.Draw("X0")
    if var == "Enu":
        ratio.SetMaximum(2.0)
        ratio.SetMinimum(0.0)
    else:
        ratio.SetMaximum(1.5)
        ratio.SetMinimum(0.5)
    ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions
    # Copied from MINERvA-101-Cross-Section/backgroundStack.py
    # same as in void MnvPlotter::DrawDataMCRatio() in MnvPlotter
    # systematic error centered at y = 1
    ratio_tot = data_hist_total.Clone()
    ratio_tot.Divide(ratio_tot, mc_hist_stat)

    # line at 1
    line = TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
    line.SetLineColor(46)
    line.SetLineWidth(2)
    line.SetLineStyle(9)
    line.Draw()
    
    ratio.Draw("X0 SAME E1")
    ratio_tot.Draw("E1 SAME")

    fit = ROOT.TF1( "fit", "pol2", 0, 20)  #plo6 is a 6th degree polynomial ypu can try other 
    ratio.Fit( "fit" )
    fit.Draw("X0 SAME")

    print("Last bin value")
    print(ratio.GetBinContent(ratio.GetNbinsX()))


    if targetZ == "99":
        mnv.AddHistoTitle("%s: Bkg Subtracted"%(trueZ), 0.05, 1)
    else:
        mnv.AddHistoTitle("Target %s %s: Bkg Subtracted"%(targetID, trueZ), 0.05, 1)
    mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.32)

    canvas1.SetLogx(False)
    if var == "x":
        canvas1.SetLogx()
    canvas1.Modified()
    canvas1.Print("DataMCratio_BkgSubtracted_EventSelection_t%s_z%02s_%s_%s_PolyFit.png"%(targetID, targetZ, var, plist))

data_hist.SetDirectory(0)
mc_hist.SetDirectory(0)

print("DONE %s %s %02s"%(plist, targetID, targetZ))
