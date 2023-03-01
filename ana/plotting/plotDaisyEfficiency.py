import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLine
from ROOT import TLegend
from ROOT import gPad
import numpy as np

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
plist = sys.argv[2]
scale = sys.argv[3]

infile= ROOT.TFile(str(dirpwd)+"/Efficiency_Daisy_%s_t99_z99_sys.root"%(plist))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
canvas2 = ROOT.TCanvas()
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale = dataPOT/mcPOT
if scale == "1":
    mcScale = 1

mcColors = ROOT.MnvColors.GetColors(ROOT.MnvColors.kGlasbeyPalette)

vars = ["Enu", "x", "pZmu1D", "pTmu", "ThetamuDeg"]


for var in vars:

    mc_histos = [ ]
    denom_histos = [ ]   

    for petal in range(0,12):
        mc_hist = infile.Get("h_mc_daisy_%s_%s"%(petal, var))
        denom_hist = infile.Get("h_truth_daisy_%s_%s"%(petal, var))

        mc_hist.Scale(mc_hist.GetNormBinWidth(), "width")
        denom_hist.Scale(denom_hist.GetNormBinWidth(), "width")

        if var == "Enu":
            mc_hist.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
            mc_hist.GetYaxis().SetTitle("Events/GeV")
        elif var == "x":
            mc_hist.GetXaxis().SetTitle("Bjorken x")
            mc_hist.GetYaxis().SetTitle("Events per unit x")
        elif var == "vtxx_vtxy":
            mc_hist.GetXaxis().SetTitle("Vertex X (cm)")
            mc_hist.GetYaxis().SetTitle("Vertex Y (cm)")
        elif var == "pTmu":
            mc_hist.GetXaxis().SetTitle("Muon p_{T} (GeV/c)")
            mc_hist.GetYaxis().SetTitle("Events/(GeV/c)")
        elif var == "pZmu1D":
            mc_hist.GetXaxis().SetTitle("Muon p_{Z} (GeV/c)")
            mc_hist.GetYaxis().SetTitle("Events/(GeV/c)")
        elif var == "ThetamuDeg":
            mc_hist.GetXaxis().SetTitle("Muon Angle (Deg)")
            mc_hist.GetYaxis().SetTitle("Events/Deg")

        mc_hist.GetXaxis().CenterTitle()
        mc_hist.GetYaxis().CenterTitle()
        mc_hist.SetMaximum(mc_hist.GetMaximum()*1.2)


        mc_hist.SetLineColor(mcColors[petal])
        mc_hist.SetLineWidth(2)
        #denom_hist.SetMarkerColor(mcColors[petal])
        #denom_hist.SetMarkerStyle(ROOT.kFullCircle)
        #denom_hist.SetMarkerSize(1)

        denom_hist.SetLineColor(mcColors[petal])
        denom_hist.SetLineStyle(2)
        denom_hist.SetLineWidth(2)


        mc_histos.append(mc_hist)
        denom_histos.append(denom_hist)

    canvas1.cd()
    mc_histos[0].Draw("HIST")
    mc_histos[0].SetMaximum(denom_histos[0].GetMaximum()*1.2)
    denom_histos[0].Draw("HIST SAME")

    for petal in range(1,12):
        mc_histos[petal].Draw("HIST SAME")
        denom_histos[petal].Draw("HIST SAME") 


    mnv.AddHistoTitle("Efficiency: Daisy Tracker", 0.04, 1)
    #mnv.AddPOTNormBox(deomPOT, mcPOT, 0.5, 0.85)
    gStyle.SetErrorX(0)

    legend = TLegend(0.50,0.50,0.60,0.87)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    for petal in range(0,12):
        legend.AddEntry(mc_histos[petal], " Num. " + str(petal), "l")
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)
    legend.Draw()

    legend2 = TLegend(0.65,0.50,0.85,0.87)
    legend2.SetFillStyle(0)
    legend2.SetBorderSize(0)
    legend2.SetTextSize(0.035)
    for petal in range(0,12):
        legend2.AddEntry(denom_histos[petal], " Denom. " + str(petal), "l")
    legend2.SetTextFont(42)
    legend2.SetTextSize(0.03)
    legend2.Draw()

    canvas1.SetLogx(False)
    if var == "x":
        canvas1.SetLogx()
    canvas1.Modified()
    canvas1.Print("Efficiency_daisy_t99_z99_%s_%s_NumDenom.png"%(var, plist))

    
    canvas2.cd()
    ratios = [ ]
    for petal in range(0,12):
        ratio_hist = mc_histos[petal].Clone()
        ratio_hist.Divide(ratio_hist, denom_histos[petal])
        ratio_hist.SetLineColor(mcColors[petal])
        ratios.append(ratio_hist)
    
    ratios[0].Draw("HIST")
    ratios[0].SetMaximum(1)
    ratios[0].GetYaxis().SetTitle("Efficiency")

    for petal in range(1,12):
       ratios[petal].Draw("HIST SAME")
    
    mnv.AddHistoTitle("Efficiency: Daisy Tracker", 0.04, 1)
    #mnv.AddPOTNormBox(deomPOT, mcPOT, 0.5, 0.85)
    gStyle.SetErrorX(0)

    legend = TLegend(0.55,0.20,0.85,0.50)
    if var == "pTmu":
        legend = TLegend(0.55,0.60,0.85,0.90)
    if var == "ThetamuDeg":
        legend = TLegend(0.55,0.60,0.85,0.90)
    legend.SetNColumns(2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    for petal in range(0,12):
        legend.AddEntry(ratios[petal], " Petal " + str(petal), "l")
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    legend.Draw()

    #canvas2.SetLogy()
    canvas2.SetLogx(False)
    if var == "x":
        canvas2.SetLogx()
    canvas2.Modified()
    canvas2.Print("Efficiency_daisy_t99_z99_%s_%s.png"%(var, plist))
    
'''
mc_histos = [ ]
denom_histos = [ ]   

for petal in range(0,12):
    mc_hist = infile.Get("h_truth2d_%s_vtxx_vtxy"%(petal))
    #denom_hist = infile.Get("selected_denom_reco_daisy_%s_vtxx_vtxy"%(petal))

    mc_hist.GetXaxis().SetTitle("Vertex X (cm)")
    mc_hist.GetYaxis().SetTitle("Vertex Y (cm)")
    mc_hist.GetYaxis().SetTitleOffset(1)
    mc_hist.GetXaxis().SetTitleOffset(1.2)
    mc_hist.GetXaxis().CenterTitle()
    mc_hist.GetYaxis().CenterTitle()

    mc_hist.SetMarkerColor(mcColors[petal])
    mc_hist.SetLineColor(mcColors[petal])
    mc_hist.SetLineWidth(2)
    #denom_hist.SetMarkerColor(mcColors[petal])
    mc_hist.SetMarkerStyle(ROOT.kFullCircle)
    mc_hist.SetMarkerSize(0.3)

    mc_histos.append(mc_hist)
    #denom_histos.append(denom_hist)


mc_histos[0].Draw("")

for petal in range(1,12):
    mc_histos[petal].Draw("SAME")



#mnv.AddHistoTitle("Daisy Tracker (MC)", 0.05, 1)

mnv.AddPlotLabel("Daisy Tracker Efficiency Denominator", 0.31, 0.87, 0.033, 12, 42)


mnv.AddPOTNormBox(denomPOT, mcPOT, 0.5, 0.85)
gStyle.SetErrorX(0)

legend = TLegend(0.15,0.9,0.9,1.01)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetNColumns(6)
legend.SetTextSize(0.03)
for petal in range(0,12):
    legend.AddEntry(mc_histos[petal], " Petal " + str(petal), "pl")
legend.SetTextFont(42)
legend.Draw()
	
gPad.Update()
canvas1.Update()
canvas1.Modified()
#canvas1.Print("Efficiency_Denom_Daisy_Tracker_nosys_vtxx_vtxy.png")

mc_histos = [ ]
denom_histos = [ ]   

for petal in range(0,12):
    mc_hist = infile.Get("h_mc2d_%s_vtxx_vtxy"%(petal))
    #denom_hist = infile.Get("selected_denom_reco_daisy_%s_vtxx_vtxy"%(petal))

    mc_hist.GetXaxis().SetTitle("Vertex X (cm)")
    mc_hist.GetYaxis().SetTitle("Vertex Y (cm)")
    mc_hist.GetYaxis().SetTitleOffset(1)
    mc_hist.GetXaxis().SetTitleOffset(1.2)
    mc_hist.GetXaxis().CenterTitle()
    mc_hist.GetYaxis().CenterTitle()

    mc_hist.SetMarkerColor(mcColors[petal])
    mc_hist.SetLineColor(mcColors[petal])
    mc_hist.SetLineWidth(2)
    #denom_hist.SetMarkerColor(mcColors[petal])
    mc_hist.SetMarkerStyle(ROOT.kFullCircle)
    mc_hist.SetMarkerSize(0.3)

    mc_histos.append(mc_hist)
    #denom_histos.append(denom_hist)


mc_histos[0].Draw("")

for petal in range(1,12):
    mc_histos[petal].Draw("SAME")



#mnv.AddHistoTitle("Daisy Tracker (MC)", 0.05, 1)

mnv.AddPlotLabel("Daisy Tracker Efficiency Numerator", 0.31, 0.87, 0.033, 12, 42)


mnv.AddPOTNormBox(denomPOT, mcPOT, 0.5, 0.85)
gStyle.SetErrorX(0)

legend = TLegend(0.15,0.9,0.9,1.01)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetNColumns(6)
legend.SetTextSize(0.03)
for petal in range(0,12):
    legend.AddEntry(mc_histos[petal], " Petal " + str(petal), "pl")
legend.SetTextFont(42)
legend.Draw()
	
gPad.Update()
canvas1.Update()
canvas1.Modified()
#canvas1.Print("Efficiency_Num_Daisy_Tracker_nosys_vtxx_vtxy.png")
'''
ratio_hist.SetDirectory(0)
denom_hist.SetDirectory(0)
mc_hist.SetDirectory(0)
print("DONE %s 99 99 daisy"%(plist))
