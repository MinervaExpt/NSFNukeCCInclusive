import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLegend


infile= ROOT.TFile("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/make_hists/planeDNNstudy/Hists_EventSelectionPlaneDNN_ME6A_nosys_xdep03_t99_z99_AntiNu.root")
infile2= ROOT.TFile("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/make_hists/planeDNNstudy/Hists_EventSelectionPlaneDNN_ME6A_nosys_xdep07_t99_z99_AntiNu.root")
infile3= ROOT.TFile("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/make_hists/planeDNNstudy/Hists_EventSelectionPlaneDNN_ME6A_nosys_xdep07more_t99_z99_AntiNu.root")

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
    mc_hist2 = infile2.Get("selected_mc_reco_%s"%var)
    mc_hist3 = infile3.Get("selected_mc_reco_%s"%var)
    data_hist = infile.Get("selected_data_reco_%s"%var)
    data_hist2 = infile2.Get("selected_data_reco_%s"%var)
    data_hist3 = infile3.Get("selected_data_reco_%s"%var)

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
    mc_hist2.Scale(mc_hist2.GetNormBinWidth(),"width")
    mc_hist2.Scale(mcScale)
    mc_hist3.Scale(mc_hist3.GetNormBinWidth(),"width")
    mc_hist3.Scale(mcScale)
    #data_hist.Scale(mcScale)
    data_hist.Scale(data_hist.GetNormBinWidth(),"width")
    data_hist2.Scale(data_hist2.GetNormBinWidth(),"width")
    data_hist3.Scale(data_hist3.GetNormBinWidth(),"width")

    mc_hist.SetLineColor(46)
    data_hist.SetLineColor(46)
    data_hist.SetLineStyle(2)

    data_hist2.SetLineColor(30)
    mc_hist2.SetLineColor(30)
    data_hist2.SetLineStyle(2)

    mc_hist3.SetLineColor(39)
    data_hist3.SetLineColor(39)
    data_hist3.SetLineStyle(2)

    mc_hist.Draw("HIST")
    data_hist.Draw("HIST SAME")
    mc_hist2.Draw("HIST SAME")
    data_hist2.Draw("HIST SAME")
    mc_hist3.Draw("HIST SAME")
    data_hist3.Draw("HIST SAME")
    mc_hist.SetMaximum(1.2*data_hist.GetMaximum())

    #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)

    '''
    Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
    If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
    If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
    '''
    
    mnv.AddHistoTitle("%s"%(name), 0.04, 1)
    mnv.AddPOTNormBox(dataPOT,mcPOT, 0.63, 0.68)

    legend = TLegend(0.55,0.7,0.75,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(data_hist, "Data x < 0.3", "fl")
    legend.AddEntry(mc_hist, "MC x < 0.3", "fl")
    legend.AddEntry(data_hist2, "Data 0.3 < x < 0.7", "fl")
    legend.AddEntry(mc_hist2, "MC 0.3 < x < 0.7", "fl")
    legend.AddEntry(data_hist3, "Data x > 0.7", "fl")
    legend.AddEntry(mc_hist3, "MC x > 0.7", "fl")


    legend.SetTextFont(42)
    legend.Draw()

    #canvas1.SetLogy()
    canvas1.Modified()
    canvas1.Print("ME6A_EventSelection_%s_nosys_Recoplane_%s_all.png"%(var, tag))

    ratio = mc_hist.Clone()
    ratio.Divide(data_hist, mc_hist)
    ratio.SetLineColor(46)
    ratio.SetMaximum(2)
    ratio.SetMinimum(0)
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.Draw("HIST")
    canvas1.Modified()

    ratio2 = mc_hist2.Clone()
    ratio2.SetLineColor(30)
    ratio2.Divide(data_hist2, mc_hist2)
    ratio2.Draw("HIST SAME")

    ratio3 = mc_hist3.Clone()
    ratio3.SetLineColor(39)
    ratio3.Divide(data_hist3, mc_hist3)
    ratio3.Draw("HIST SAME")

    legend = TLegend(0.60,0.7,0.85,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(ratio, "x < 0.3", "fl")
    legend.AddEntry(ratio2, "0.3 < x < 0.7", "fl")
    legend.AddEntry(ratio3, "x > 0.7", "fl")

    legend.SetTextFont(42)
    legend.Draw()

    mnv.AddHistoTitle("Reconstructed plane number %s"%(name), 0.04, 1)

    canvas1.Print("ME6A_EventSelection_%s_nosys_Recoplane_ratio_%s_all.png"%(var, tag))

 
 
mc_hist.SetDirectory(0)

raw_input("Done")
