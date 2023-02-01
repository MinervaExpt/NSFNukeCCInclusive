import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle


infile= ROOT.TFile("Hists_EventSelectionTracker_ML_ME6A_sys_t99_z99_AntiNu_v2.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale =  dataPOT/mcPOT

vars = ['planeDNN', 'vtxz']#"Enu", "x"]

for var in vars:

    mc_hist = infile.Get("selected_mc_sb_%s_NotTracker_true"%var)

    if var == "Enu":
        mc_hist.GetXaxis().SetTitle("Neutrino Energy (GeV)")
        mc_hist.GetYaxis().SetTitle("Events/GeV")

    if var == "x":
        mc_hist.GetXaxis().SetTitle("Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
    
    if var == "planeDNN":
        mc_hist.GetXaxis().SetTitle("Plane Number")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
        mc_hist.GetXaxis().SetRangeUser(0,200)
    
    if var == "ANNPlaneProb":
        mc_hist.GetXaxis().SetTitle("Plane Probability")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")

    if var == "vtxz":
        mc_hist.GetXaxis().SetTitle("Vertex Z (cm)")
        mc_hist.GetYaxis().SetTitle("Events/cm")
        mc_hist.GetXaxis().SetRangeUser(0,1000)


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
    #mc_hist.GetXaxis().SetRangeUser(66,173)

    #mnv.ApplyStyle(1)
    mnv.DrawMCWithErrorBand(mc_hist, mcScale)

    mnv.AddHistoTitle("Tracker", 0.05, 1)
    mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)

    canvas1.Modified()
    canvas1.Print("ME6A_Tracker_EventSelection_TrueOutOfTracker_%s_sys.png"%(var))


mc_hist.SetDirectory(0)

raw_input("Done")
