import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLegend


infile= ROOT.TFile("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/make_hists/planeDNNstudy/Hists_EventSelectionPlaneDNN_ME6A_nosys_allMaterial_AntiNu.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

name = "Target 2 Iron"
tag = "target_2_iron"
name = "All material"
tag = "all_material"


mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale =  dataPOT/mcPOT

vars = ['planeDNN']#"Enu", "x"]

for var in vars:

    mc_hist = infile.Get("selected_mc_reco_%s"%var)
    data_hist = infile.Get("selected_data_reco_%s"%var)

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
        mc_hist.GetXaxis().SetRangeUser(7,173)


    mc_hist.GetXaxis().CenterTitle()
    mc_hist.GetYaxis().CenterTitle()
    #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
    #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
    #mc_hist.GetXaxis().SetRangeUser(66,173)

    #mnv.ApplyStyle(1)
    mc_hist.Scale(mc_hist.GetNormBinWidth(),"width")
    mc_hist.Scale(mcScale)
    #data_hist.Scale(mcScale)
    data_hist.Scale(data_hist.GetNormBinWidth(),"width")

    data_hist.SetLineColor(ROOT.kRed)
    mc_hist.Draw("HIST")
    data_hist.Draw("HIST SAME")
    #mc_hist.SetMaximum(1.7*true_hist.GetMaximum())

    #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)

    '''
    Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
    If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
    If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
    '''
    
    mnv.AddHistoTitle("%s"%name, 0.04, 1)
    mnv.AddPOTNormBox(dataPOT,mcPOT, 0.63, 0.68)

    legend = TLegend(0.40,0.7,0.55,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(data_hist, "Data reconstructed plane number", "fl")
    legend.AddEntry(mc_hist, " MC Reconstructed plane number ", "fl")
    legend.SetTextFont(42)
    legend.Draw()

    canvas1.Modified()
    canvas1.Print("ME6A_EventSelection_%s_nosys_Recoplane_%s.png"%(var, tag))

    ratio = mc_hist.Clone()
    ratio.Divide(data_hist, mc_hist)
    ratio.SetLineColor(ROOT.kRed)
    ratio.SetMaximum(1.5)
    ratio.SetMinimum(.5)
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.Draw("HIST")
    canvas1.Modified()

    legend = TLegend(0.40,0.7,0.55,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(ratio, "Reconstructed plane number", "fl")
    legend.SetTextFont(42)
    legend.Draw()

    mnv.AddHistoTitle("%s"%name, 0.04, 1)

    canvas1.Print("ME6A_EventSelection_%s_nosys_Recoplane_ratio_%s.png"%(var, tag))

 
 
mc_hist.SetDirectory(0)

raw_input("Done")
