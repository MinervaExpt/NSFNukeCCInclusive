import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import TLegend
from ROOT import THStack

infile= ROOT.TFile("Hists_EventSelectionTracker_ML_ME6A_nosys_t99_z99_AntiNu.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale =  dataPOT/mcPOT

myvariable = sys.argv[1]

signal = infile.Get("selected_mc_sb_%s_Signal"%myvariable)
# All background categories
NotTracker = infile.Get("selected_mc_sb_%s_NotTracker"%myvariable)
WrongSign = infile.Get("selected_mc_sb_%s_WrongSign"%myvariable)
NC = infile.Get("selected_mc_sb_%s_NC"%myvariable)
NotEmu = infile.Get("selected_mc_sb_%s_NotEmu"%myvariable)

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

signal.SetLineColor(30)
NotTracker.SetLineColor(39)
WrongSign.SetLineColor(41)
NC.SetLineColor(32)
NotEmu.SetLineColor(38)

signal.SetFillColor(0)
NotTracker.SetFillColor(0)
WrongSign.SetFillColor(0)
NC.SetFillColor(0)
NotEmu.SetFillColor(0)

signal.SetLineWidth(2)
NotTracker.SetLineWidth(2)
WrongSign.SetLineWidth(2)
NC.SetLineWidth(2)
NotEmu.SetLineWidth(2)

data_hist = infile.Get("selected_data_reco_%s"%myvariable)
data_hist.Scale(data_hist.GetBinWidth(1), "width")
data_hist.SetMarkerStyle(20)
data_hist.SetMarkerSize(1)
data_hist.SetMarkerColor(1)
data_hist.SetLineWidth(1)
data_hist.SetLineStyle(1)
data_hist.SetLineColor(1)
data_hist.SetFillColor(0)
data_hist.Draw("HIST p E1 X0") # for error bars, suppressed error bars along X

NotTracker.Draw("HIST SAME")
signal.Draw("HIST SAME")
WrongSign.Draw("HIST SAME")
if NC.GetEntries() != 0:   
    NC.Draw("HIST SAME")
if NotEmu.GetEntries() != 0:
    NotEmu.Draw("HIST SAME")

data_hist.SetMinimum(0.1)
data_hist.SetMaximum(data_hist.GetMaximum()*1.2)


if myvariable == "Enu":
    data_hist.GetXaxis().SetTitle("Reconstructed Neutrino Energy [GeV]")
    data_hist.GetYaxis().SetTitle("Events/GeV")

if myvariable == "x":
    data_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
    data_hist.GetYaxis().SetTitle("Events (norm.)")


data_hist.GetXaxis().CenterTitle()
data_hist.GetYaxis().CenterTitle()
data_hist.Draw("SAME HIST p E1 X0")

mnv.AddHistoTitle("Tracker", 0.05, 1)

legend = TLegend(0.55,0.64,0.80,0.89)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextSize(0.035)
legend.AddEntry(data_hist, " Data", "ep")
legend.AddEntry(signal, " Signal", "fl")
if WrongSign.GetEntries() != 0:
    legend.AddEntry(WrongSign, " Wrong sign", "fl")
if NC.GetEntries() != 0:
    legend.AddEntry(NC, " Neutral current", "fl")
if NotTracker.GetEntries() != 0:
    legend.AddEntry(NotTracker, " Outside tracker", "fl")
if NotEmu.GetEntries() != 0:
    legend.AddEntry(NotEmu, " Outside muon energy", "fl")
legend.SetTextFont(62)
legend.Draw()

#canvas1.SetLogy()
canvas1.Modified()
canvas1.Print("ME6A_Tracker_NOTStackBackground_%s_sys.png"%myvariable)

raw_input("Done")