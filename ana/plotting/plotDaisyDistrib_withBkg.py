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


infile= ROOT.TFile(str(dirpwd)+"EventSelection_daisy_%s_t99_z99_sys.root"%(plist))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale = dataPOT/mcPOT
if scale == "1":
    mcScale = 1

mcColors = ROOT.MnvColors.GetColors(ROOT.MnvColors.kGlasbeyPalette)

vars = ["Enu", "x", "pZmu1D", "pTmu1D", "ThetamuDeg"]

for var in vars:

    mc_histos = [ ]
    bkg_histos = [ ]
    data_histos = [ ]   

    for petal in range(0,12):
        mc_hist = infile.Get("selected_mc_reco_daisy_%s_%s"%(petal, var))
        bkg_hist = infile.Get("selected_mc_reco_bkg_daisy_%s_%s"%(petal, var))
        data_hist = infile.Get("selected_data_reco_daisy_%s_%s"%(petal, var))

        mc_hist.Scale(mc_hist.GetNormBinWidth(), "width")
        bkg_hist.Scale(bkg_hist.GetNormBinWidth(), "width")
        data_hist.Scale(data_hist.GetNormBinWidth(), "width")

        bkg_hist.Scale(mcScale)
        mc_hist.Scale(mcScale)

        if var == "Enu":
            mc_hist.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")
            mc_hist.GetYaxis().SetTitle("Events/GeV")
        elif var == "x":
            mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
            mc_hist.GetYaxis().SetTitle("Events per unit x")
        elif var == "vtxx_vtxy":
            mc_hist.GetXaxis().SetTitle("Reconstructed Vertex X (cm)")
            mc_hist.GetYaxis().SetTitle("Reconstructed Vertex Y (cm)")
        elif var == "pTmu1D":
            mc_hist.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
            mc_hist.GetYaxis().SetTitle("Events/(GeV/c)")
        elif var == "pZmu1D":
            mc_hist.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
            mc_hist.GetYaxis().SetTitle("Events/(GeV/c)")  
        elif var == "ThetamuDeg":
            mc_hist.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)")
            mc_hist.GetYaxis().SetTitle("Events/Deg")

        mc_hist.GetXaxis().CenterTitle()
        mc_hist.GetYaxis().CenterTitle()
        if var == "Enu":
            mc_hist.SetMaximum(mc_hist.GetMaximum()*1.25)
        if var == "x":
            mc_hist.SetMaximum(mc_hist.GetMaximum()*1.2)
        if var == "ThetamuDeg":
            mc_hist.SetMaximum(mc_hist.GetMaximum()*1.32)
        if var == "pZmu1D":
            mc_hist.SetMaximum(mc_hist.GetMaximum()*1.25)
        if var == "pTmu1D":
            mc_hist.SetMaximum(mc_hist.GetMaximum()*1.25)



        mc_hist.SetLineColor(mcColors[petal])
        mc_hist.SetLineWidth(2)

        bkg_hist.SetLineColor(mcColors[petal])
        bkg_hist.SetLineWidth(2)
        bkg_hist.SetLineStyle(7)

        data_hist.SetMarkerColor(mcColors[petal])
        data_hist.SetMarkerStyle(ROOT.kFullCircle)
        data_hist.SetMarkerSize(1)

        mc_histos.append(mc_hist)
        bkg_histos.append(bkg_hist)
        data_histos.append(data_hist)


    mc_histos[0].Draw("HIST")
    bkg_histos[0].Draw("HIST SAME")
    data_histos[0].Draw("HIST p SAME")

    for petal in range(1,12):
        mc_histos[petal].Draw("HIST SAME")
        bkg_histos[petal].Draw("HIST SAME")
        data_histos[petal].Draw("HIST p SAME") 


    mnv.AddHistoTitle("Daisy Tracker", 0.04, 1)
    if var == "x":
        mnv.AddPOTNormBox(dataPOT, mcPOT, 0.4, 0.3)
    elif var == "ThetamuDeg":
        mnv.AddPOTNormBox(dataPOT, mcPOT, 0.4, 0.3)
    else:
        mnv.AddPOTNormBox(dataPOT, mcPOT, 0.5, 0.85)
    gStyle.SetErrorX(0)

    legend = TLegend(0.65,0.50,0.85,0.87)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    for petal in range(0,12):
        legend.AddEntry(mc_histos[petal], " Petal " + str(petal), "l")
    legend.SetTextFont(42)
    legend.Draw()

    canvas1.SetLogx(False)
    if var == "x":
        canvas1.SetLogx()
    canvas1.Modified()
    canvas1.Print("EventSelection_daisy_t99_z99_%s_%s_wsimbkg.png"%(var, plist))
'''
mc_histos = [ ]
data_histos = [ ]   

for petal in range(0,12):
    mc_hist = infile.Get("selected_mc2d_reco_daisy_%s_vtxx_vtxy"%(petal))
    #data_hist = infile.Get("selected_data_reco_daisy_%s_vtxx_vtxy"%(petal))

    mc_hist.GetXaxis().SetTitle("Reconstructed Vertex X (cm)")
    mc_hist.GetYaxis().SetTitle("Reconstructed Vertex Y (cm)")
    mc_hist.GetYaxis().SetTitleOffset(1)
    mc_hist.GetXaxis().SetTitleOffset(1.2)
    mc_hist.GetXaxis().CenterTitle()
    mc_hist.GetYaxis().CenterTitle()

    mc_hist.SetMarkerColor(mcColors[petal])
    mc_hist.SetLineColor(mcColors[petal])
    mc_hist.SetLineWidth(2)
    #data_hist.SetMarkerColor(mcColors[petal])
    mc_hist.SetMarkerStyle(ROOT.kFullCircle)
    mc_hist.SetMarkerSize(0.3)

    mc_histos.append(mc_hist)
    #data_histos.append(data_hist)


mc_histos[0].Draw("")

for petal in range(1,12):
    mc_histos[petal].Draw("SAME")



#mnv.AddHistoTitle("Daisy Tracker (MC)", 0.05, 1)

mnv.AddPlotLabel("Daisy Tracker (MC)", 0.25, 0.87, 0.033, 12, 42)


mnv.AddPOTNormBox(dataPOT, mcPOT, 0.5, 0.85)
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
canvas1.Print("Daisy_Tracker_nosys_vtxx_vtxy.png")
'''

print("DONE %s 99 99 daisy"%(plist))