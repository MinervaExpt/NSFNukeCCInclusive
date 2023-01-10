import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLegend


infile= ROOT.TFile("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/make_hists/planeDNNstudy/Hists_EventSelectionPlaneDNN_ME6A_nosys_t2_z26_AntiNu.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()
name = "Target 2 Iron"
tag = "target_2_iron"
#name = "All material"
#tag = "all_material"

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale =  dataPOT/mcPOT

vars = ['planeDiff']#"Enu", "x"]

for var in vars:

    mc_hist = infile.Get("selected_mc_reco_%s"%var)
    true_hist = infile.Get("selected_mc_truth_%s"%var)
    data_hist = infile.Get("selected_data_reco_%s"%var)

    print(mc_hist.GetEntries()) 
    print(mc_hist.Integral()) 
    print(data_hist.GetEntries())

    print("Reconstructed(0) - reconstructed(1) plane diff")
    for bin in range(0, mc_hist.GetNbinsX()+1):
        print(mc_hist.GetBinContent(bin)*100./mc_hist.GetEntries())

    print("Reconstructed - true plane diff")
    for bin in range(0, true_hist.GetNbinsX()+1):
        print(true_hist.GetBinContent(bin)*100./true_hist.GetEntries())


    if var == "planeDiff":
        mc_hist.GetXaxis().SetTitle("Plane Difference")
        mc_hist.GetYaxis().SetTitle("Events %")
        mc_hist.SetMaximum(100 )

    if var == "x":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
    
    if var == "GetThetamuDeg":
        mc_hist.GetXaxis().SetTitle("Reconstructed Muon Angle")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
    
    if var == "planeDNN":
        mc_hist.GetXaxis().SetTitle("Reconstructed Plane Number")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
        mc_hist.GetXaxis().SetRangeUser(66,173)

    
    if var == "ANNPlaneProb":
        mc_hist.GetXaxis().SetTitle("Plane Probability")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")

    if var == "vtxz":
        mc_hist.GetXaxis().SetTitle("Reconstructed Vertex Z (cm)")
        mc_hist.GetYaxis().SetTitle("Events/cm")
        mc_hist.GetXaxis().SetRangeUser(599,842.2)

    if var == "vtxy":
        mc_hist.GetXaxis().SetTitle("Reconstructed Vertex Y (cm)")
        mc_hist.GetYaxis().SetTitle("Events/cm")

    if var == "vtxx":
        mc_hist.GetXaxis().SetTitle("Reconstructed Vertex X (cm)")
        mc_hist.GetYaxis().SetTitle("Events/cm")

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
    mc_hist.GetYaxis().SetTitleOffset(0.95)
    #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
    #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
    #mc_hist.GetXaxis().SetRangeUser(66,173)

    #mnv.ApplyStyle(1)
    #mc_hist.Scale(mc_hist.GetNormBinWidth(),"width")
    #true_hist.Scale(true_hist.GetNormBinWidth(),"width")
    
    #mc_hist.Scale(mcScale)
    #true_hist.Scale(mcScale)

    #data_hist.Scale(data_hist.GetNormBinWidth(),"width")

    mc_hist.Scale(100/mc_hist.GetEntries())
    true_hist.Scale(100/true_hist.GetEntries())
    data_hist.Scale(100/data_hist.GetEntries())

    data_hist.SetLineColor(ROOT.kRed)
    true_hist.SetLineColor(ROOT.kBlue)
    mc_hist.Draw("HIST")
    true_hist.Draw("HIST SAME")
    data_hist.Draw("HIST SAME")
    mc_hist.SetMaximum(100)
    
    #mc_hist.SetMaximum(1.2*data_hist.GetMaximum())

    #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)

    '''
    Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
    If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
    If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
    '''
    
    mnv.AddHistoTitle("%s"%name, 0.05, 1)
    #mnv.AddHistoTitle("Plane Reconstructed using ANN_segments, 0 - ANN_segments, 1", 0.023, 1)
    #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.65, 0.65)

    legend = TLegend(0.50 ,0.7,0.67,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.032)
    legend.AddEntry(mc_hist, " Simulation Reco(0)-Reco(1)", "fl")
    legend.AddEntry(true_hist, " Simulation Reco(0)-True", "fl")
    legend.AddEntry(data_hist, " Data Reco(0)-Reco(1)", "fl")
    legend.SetTextFont(42)
    legend.Draw()

    canvas1.Modified()
    canvas1.Print("ME6A_EventSelection_%s_nosys_PlaneStudy_%s.png"%(var, tag))

 
 
mc_hist.SetDirectory(0)

raw_input("Done")
