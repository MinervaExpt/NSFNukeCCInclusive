import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLine
from ROOT import TLegend

path = "/minerva/app/users/anezkak/MAT_GitHub/opt/lib/data/flux_daisy/targets_2345_temp/flux_with_errors/"
tracker_file = ROOT.TFile(str(path)+"flux_tracker.root")
iron_file = ROOT.TFile(str(path)+"flux_iron.root ")
lead_file = ROOT.TFile(str(path)+"flux_lead.root ")
carbon_file = ROOT.TFile(str(path)+"flux_carbon.root ")
water_file = ROOT.TFile(str(path)+"flux_water.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

tracker = tracker_file.Get("flux")
iron = iron_file.Get("flux")
lead = lead_file.Get("flux")
carbon = carbon_file.Get("flux")
water = water_file.Get("flux")

iron.SetLineColor(ROOT.kBlue)
lead.SetLineColor(ROOT.kGreen)
carbon.SetLineColor(ROOT.kRed)
water.SetLineColor(ROOT.kCyan)

#fix title
tracker.GetYaxis().CenterTitle(True)
tracker.GetYaxis().SetTitleOffset(1.3)
tracker.GetYaxis().SetTitleSize(0.05)
tracker.GetYaxis().SetLabelSize(0.05)
tracker.GetXaxis().CenterTitle(True)
tracker.GetXaxis().SetTitle("Neutrino Energy (GeV)")
#tracker.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
tracker.GetXaxis().SetTitleSize(0.05)
tracker.GetXaxis().SetLabelSize(0.05)
tracker.GetYaxis().SetTitle("#nu / m^{2} / P.O.T/ GeV")
#tracker.GetYaxis().SetTitle("#bar{#nu} / m^{2} / P.O.T/ GeV")

tracker.Draw("HIST")
iron.Draw("HIST SAME")
lead.Draw("HIST SAME")
carbon.Draw("HIST SAME")
#water.Draw("HIST SAME")

mnv.AddHistoTitle("Flux Neutrinos", 0.05, 1)
gStyle.SetErrorX(0)

legend = TLegend(0.65,0.65,0.85,0.89)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextSize(0.035)
#legend.AddEntry(reco, " Simulation", "l")
#legend.AddEntry(data_hist, " Data", "ep")
legend.AddEntry(tracker, " Tracker (C) ", "l")
legend.AddEntry(iron, " Iron", "l")
legend.AddEntry(lead, " Lead", "l")
legend.AddEntry(carbon, " Carbon", "l")
legend.AddEntry(water, " Water (O)", "l")
legend.SetTextFont(42)
legend.Draw()

canvas1.Modified()
canvas1.Print("Neutrino_Fluxes_Jeffrey.png")


# ----------------------------------------------------------------------------
# Data/MC ratio

#tracker.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
tracker.GetXaxis().SetTitle("Neutrino Energy (GeV)")
tracker.GetYaxis().SetTitle("Ratio to Tracker")
tracker.GetXaxis().CenterTitle()
tracker.GetYaxis().CenterTitle()

denom = tracker.Clone()

iron.Divide(iron,denom)
lead.Divide(lead,denom)
carbon.Divide(carbon,denom)
water.Divide(water,denom)
tracker.Divide(tracker,denom)
tracker.Draw("HIST")
iron.Draw("HIST SAME")
lead.Draw("HIST SAME")
carbon.Draw("HIST SAME")
water.Draw("HIST SAME")

tracker.SetMinimum(0)
tracker.SetMaximum(2)

mnv.AddHistoTitle("Flux Neutrinos", 0.05, 1)
gStyle.SetErrorX(0)

legend = TLegend(0.65,0.65,0.85,0.89)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextSize(0.035)
#legend.AddEntry(reco, " Simulation", "l")
#legend.AddEntry(data_hist, " Data", "ep")
legend.AddEntry(tracker, " Tracker (C) ", "l")
legend.AddEntry(iron, " Iron", "l")
legend.AddEntry(lead, " Lead", "l")
legend.AddEntry(carbon, " Carbon", "l")
legend.AddEntry(water, " Water (O)", "l")
legend.SetTextFont(42)
legend.Draw()

canvas1.Modified()
canvas1.Print("Neutrino_Fluxes_ratio_Jeffrey.png")


raw_input("Done")
