import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLine
from ROOT import TLegend
from ROOT import gPad
import numpy as np

dirpwd = sys.argv[1]
plist = sys.argv[2]

targetID = 99
targetZ = 99

infile= ROOT.TFile(str(dirpwd)+"/BkgSubtracted_EventSelection_daisy_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale = None
if len(sys.argv) > 2:
    mcScale = 1
else:
    mcScale =  dataPOT/mcPOT

mcColors = ROOT.MnvColors.GetColors(ROOT.MnvColors.kGlasbeyPalette)

vars = ["Enu", "x"]

for var in vars:

    mc_histos = [ ]
    data_histos = [ ]   

    for petal in range(0,12):
        mc_hist = infile.Get("h_bkg_subtracted_mc_daisy_%s_%s"%(petal, var))
        data_hist = infile.Get("h_bkg_subtracted_data_daisy_%s_%s"%(petal, var))

        mc_hist.Scale(mc_hist.GetNormBinWidth(), "width")
        data_hist.Scale(data_hist.GetNormBinWidth(), "width")
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
        mc_hist.GetXaxis().CenterTitle()
        mc_hist.GetYaxis().CenterTitle()
        mc_hist.SetMaximum(mc_hist.GetMaximum()*1.2)


        mc_hist.SetLineColor(mcColors[petal])
        mc_hist.SetLineWidth(2)
        data_hist.SetMarkerColor(mcColors[petal])
        data_hist.SetMarkerStyle(ROOT.kFullCircle)
        data_hist.SetMarkerSize(1)

        mc_histos.append(mc_hist)
        data_histos.append(data_hist)


    mc_histos[0].Draw("HIST")
    data_histos[0].Draw("HIST p SAME")

    for petal in range(1,12):
        mc_histos[petal].Draw("HIST SAME")
        data_histos[petal].Draw("HIST p SAME") 


    mnv.AddHistoTitle("Daisy Tracker: Bkg Subtracted", 0.04, 1)
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

    canvas1.Modified()
    canvas1.Print("BkgSubtracted_EventSelection_daisy_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))


mc_histos = [ ]
data_histos = [ ]   

'''
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
print("DONE %s %s %02s Daisy"%(plist, targetID, targetZ))
