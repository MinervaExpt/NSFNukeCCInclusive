import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import TLegend
from ROOT import THStack
from ROOT import gPad

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
plist = sys.argv[2]
scale = sys.argv[3]

targetID = 99
targetZ = 99

infile= ROOT.TFile(str(dirpwd)+"EventSelection_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale = dataPOT/mcPOT
if scale == "1":
    mcScale = 1

vars = ["Enu", "x", "pTmu", "pZmu"]

for var in vars:
    total = infile.Get("selected_mc_reco_%s"%var)
    signal = infile.Get("selected_mc_reco_signal_%s"%var)
    # All background categories
    NotTracker = infile.Get("selected_mc_reco_NotTracker_%s"%var)
    WrongSign = infile.Get("selected_mc_reco_WrongSign_%s"%var)
    NC = infile.Get("selected_mc_reco_NC_%s"%var)
    NotEmu = infile.Get("selected_mc_reco_NotEmu_%s"%var)

    data_hist = infile.Get("selected_data_reco_%s"%var)

    # Background percentage
    signal_per = round(signal.Integral()*100/total.Integral(),1)
    NotTracker_per = round(NotTracker.Integral()*100/total.Integral(),1)
    WrongSign_per = round(WrongSign.Integral()*100/total.Integral(),1)
    NC_per = round(NC.Integral()*100/total.Integral(),1)
    NotEmu_per = round(NotEmu.Integral()*100/total.Integral(),1) #100 - signal_per - US_hist_per - DS_hist_per - other_hist_per - WrongSign_per - NC_per
    print(NotEmu_per)

    signal.Scale(signal.GetBinWidth(1), "width")
    NotTracker.Scale(NotTracker.GetBinWidth(1), "width")
    WrongSign.Scale(WrongSign.GetBinWidth(1), "width")
    NC.Scale(NC.GetBinWidth(1), "width")
    NotEmu.Scale(NotEmu.GetBinWidth(1), "width")

    signal.Scale(mcScale)
    NotTracker.Scale(mcScale)
    WrongSign.Scale(mcScale)
    NC.Scale(mcScale)
    NotEmu.Scale(mcScale)

    signal.SetLineColor(ROOT.kGray+2)
    NotTracker.SetLineColor(30)
    WrongSign.SetLineColor(ROOT.kRed-6)
    NC.SetLineColor(ROOT.kBlue-5)
    NotEmu.SetLineColor(ROOT.kGray)

    signal.SetFillColor(ROOT.kGray+2)
    signal.SetFillStyle(3002)
    NotTracker.SetFillColor(30)
    WrongSign.SetFillColor(ROOT.kRed-6)
    NC.SetFillColor(ROOT.kBlue-5)
    NotEmu.SetFillColor(ROOT.kGray)

    signal.SetLineWidth(1)
    NotTracker.SetLineWidth(2)
    WrongSign.SetLineWidth(2)
    NC.SetLineWidth(2)
    NotEmu.SetLineWidth(2)


    data_hist.Scale(data_hist.GetBinWidth(1), "width")
    data_hist.SetMarkerStyle(20)
    data_hist.SetMarkerSize(1)
    data_hist.SetMarkerColor(1)
    data_hist.SetLineWidth(1)
    data_hist.SetLineStyle(1)
    data_hist.SetLineColor(1)
    data_hist.SetMaximum(data_hist.GetMaximum()*1.3)
    data_hist.Draw("HIST p E1 X0") # for error bars, suppressed error bars along X


    stack = ROOT.THStack("stack","stack")
    if NotEmu.GetEntries() != 0:
        stack.Add(NotEmu)
    if NC.GetEntries() != 0:
        stack.Add(NC)
    stack.Add(WrongSign)
    stack.Add(NotTracker)
    stack.Add(signal)
    stack.Draw("SAME HIST")

    if var == "Enu":
        data_hist.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")
        data_hist.GetYaxis().SetTitle("Events/GeV")

    if var == "x":
        data_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        data_hist.GetYaxis().SetTitle("Events (norm.)")
    
    if var == "pTmu":
        data_hist.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
        data_hist.GetYaxis().SetTitle("Events/(GeV/c)")

    if var == "pZmu":
        data_hist.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
        data_hist.GetYaxis().SetTitle("Events/(GeV/c)")

    gPad.RedrawAxis()
    gPad.Update()

    data_hist.GetXaxis().CenterTitle()
    data_hist.GetYaxis().CenterTitle()

    data_hist.Draw("SAME HIST p E1 X0") # for error bars, suppressed error bars along X

    mnv.AddHistoTitle("Tracker", 0.05, 1)

    legend = TLegend(0.40,0.60,0.80,0.85)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(data_hist, " Data", "ep")
    legend.AddEntry(signal, " Signal ("+ str(signal_per) +"%)", "f")
    if WrongSign.GetEntries() != 0:
        legend.AddEntry(WrongSign, " Wrong sign ("+ str(WrongSign_per) +"%)", "fl")
    if NC.GetEntries() != 0:
        legend.AddEntry(NC, " Neutral current ("+ str(NC_per) +"%)", "fl")
    if NotTracker.GetEntries() != 0:
        legend.AddEntry(NotTracker, " Outside tracker ("+ str(NotTracker_per) +"%)", "fl")
    if NotEmu.GetEntries() != 0:
        legend.AddEntry(NotEmu, " Outside muon energy ("+ str(NotEmu_per) +"%)", "fl")
    legend.SetTextFont(42)
    legend.Draw()

    canvas1.SetLogy(False)
    canvas1.Modified()
    canvas1.Print("EventSelection_Bkg_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))

    #data_hist.SetMaximum(data_hist.GetMaximum()*10)
    #data_hist.SetMinimum(0.1)
    #canvas1.SetLogy()
    #canvas1.Print("ME6A_Tracker_StackBackground_%s_log.png"%var)



print("DONE %s %s %02s"%(plist, targetID, targetZ))