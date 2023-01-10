import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLegend


infile= ROOT.TFile("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/make_hists/planeDNNstudy/Hists_EventSelectionPlaneDNN_ME6A_nosys_allMaterial_AntiNu.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale =  dataPOT/mcPOT

vars = ['planeDNN']#"Enu", "x"]

for var in vars:

    mc_hist = infile.Get("selected_mc_reco_%s"%var)
    true_hist = infile.Get("selected_mc_truth_%s"%var)

    print(mc_hist.GetEntries()) 
    print(mc_hist.Integral()) 
    #print(true_hist.GetEntries())

    #print("Reco")
    #for bin in range(0, mc_hist.GetNbinsX()+1):
    #    print(mc_hist.GetBinContent(bin)*100./mc_hist.GetEntries())

    #print("True")
    #for bin in range(0, true_hist.GetNbinsX()+1):
    #    print(true_hist.GetBinContent(bin)*100./true_hist.GetEntries())


    
    if var == "planeDNN":
        mc_hist.GetXaxis().SetTitle("Plane Number")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
        mc_hist.GetXaxis().SetRangeUser(0,173)


    mc_hist.GetXaxis().CenterTitle()
    mc_hist.GetYaxis().CenterTitle()
    #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
    #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
    #mc_hist.GetXaxis().SetRangeUser(66,173)

    #mnv.ApplyStyle(1)
    mc_hist.Scale(mc_hist.GetNormBinWidth(),"width")
    mc_hist.Scale(mcScale)
    true_hist.Scale(mcScale)
    true_hist.Scale(true_hist.GetNormBinWidth(),"width")

    true_hist.SetLineColor(ROOT.kRed)
    mc_hist.Draw("HIST")
    true_hist.Draw("HIST SAME")
    #mc_hist.SetMaximum(1.7*true_hist.GetMaximum())

    #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)

    '''
    Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
    If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
    If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
    '''
    
    mnv.AddHistoTitle("True and (most probable) reconstructed plane number", 0.03, 1)
    #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.73, 0.88)

    legend = TLegend(0.48,0.7,0.62,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(true_hist, " True plane number", "fl")
    legend.AddEntry(mc_hist, " Reconstructed plane number", "fl")
    legend.SetTextFont(42)
    legend.Draw()

    canvas1.Modified()
    canvas1.Print("ME6A_EventSelection_%s_nosys_PlaneTrueReco.png"%(var))

 
 
mc_hist.SetDirectory(0)

raw_input("Done")
