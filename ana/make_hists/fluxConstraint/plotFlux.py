import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLine
from ROOT import TLegend

path = "/minerva/data2/users/anezkak/fluxFinal/Antineutrinos/"
tracker_file = ROOT.TFile(str(path)+"flux_tracker.root")
iron_file = ROOT.TFile(str(path)+"flux_iron.root")
lead_file = ROOT.TFile(str(path)+"flux_lead.root")
carbon_file = ROOT.TFile(str(path)+"flux_carbon.root")
#water_file = ROOT.TFile(str(path)+"flux_water_apothem.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()


infile = ROOT.TFile(str(path)+"CombinedPlaylists_AntiNu_sys_Enu_vs_Petal.root ")
mcPOT = infile.Get("MCPOT").GetVal()
#mcPOTwater = 1.80406844465e+21

tracker = tracker_file.Get("flux")
iron = iron_file.Get("flux")
lead = lead_file.Get("flux")
carbon = carbon_file.Get("flux")
#water = water_file.Get("flux")

iron.SetLineColor(38)
lead.SetLineColor(30)
carbon.SetLineColor(46)
tracker.SetLineColor(12)

iron.SetLineWidth(2)
lead.SetLineWidth(2)
carbon.SetLineWidth(2)
tracker.SetLineWidth(2)
#water.SetLineColor(ROOT.kCyan)

#fix title
tracker.GetYaxis().CenterTitle(True)
tracker.GetYaxis().SetTitleOffset(1.3)
tracker.GetYaxis().SetTitleSize(0.05)
tracker.GetYaxis().SetLabelSize(0.05)
tracker.GetXaxis().CenterTitle(True)
#tracker.GetXaxis().SetTitle("Neutrino Energy (GeV)")
tracker.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
tracker.GetXaxis().SetTitleSize(0.05)
tracker.GetXaxis().SetLabelSize(0.05)
#tracker.GetYaxis().SetTitle("#nu / m^{2} / P.O.T/ GeV")
tracker.GetYaxis().SetTitle("#bar{#nu} / m^{2} / P.O.T/ GeV")

tracker.Draw("HIST C")
iron.Draw("HIST C SAME")
lead.Draw("HIST C SAME")
carbon.Draw("HIST C SAME")
#water.Draw("HIST SAME")

tracker.SetMaximum(tracker.GetMaximum()*1.15)

#mnv.AddPlotLabel("MC POT "+ "{:.2e}".format(mcPOT), 0.29, 0.87, 0.04, 12, 42)
#mnv.AddPlotLabel("MC POT "+ "{:.2e}".format(mcPOT) + " (Water "+ "{:.2e}".format(mcPOTwater) + ")", 0.38, 0.87, 0.033, 12, 42)
#mnv.AddPlotLabel("MC POT "+ "{:.2e}".format(mcPOT), 0.28, 0.87, 0.033, 12, 42)

#mnv.AddHistoTitle("Neutrino Flux", 0.05, 1)
mnv.AddHistoTitle("Antineutrino Flux", 0.05, 1)
gStyle.SetErrorX(0)

legend = TLegend(0.60,0.65,0.85,0.89)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextSize(0.040)
#legend.AddEntry(reco, " Simulation", "l")
#legend.AddEntry(data_hist, " Data", "ep")
legend.AddEntry(tracker, " Tracker (C) ", "l")
legend.AddEntry(iron, " Iron", "l")
legend.AddEntry(lead, " Lead", "l")
legend.AddEntry(carbon, " Carbon", "l")
#legend.AddEntry(water, " Water (O)", "l")
legend.SetTextFont(42)
legend.Draw()

canvas1.Modified()
canvas1.Print("AllAntiNu_Fluxes_sys.png")


# ----------------------------------------------------------------------------
# Data/MC ratio

tracker.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
#tracker.GetXaxis().SetTitle("Neutrino Energy (GeV)")
tracker.GetYaxis().SetTitle("Ratio to Tracker")
tracker.GetXaxis().CenterTitle()
tracker.GetYaxis().CenterTitle()

denom = tracker.Clone()

iron.Divide(iron,denom)
lead.Divide(lead,denom)
carbon.Divide(carbon,denom)
#water.Divide(water,denom)
tracker.Divide(tracker,denom)
'''
iron.Smooth(10)
lead.Smooth(10)
carbon.Smooth(10)
water.Smooth(10)
'''
tracker.Draw("HIST C")
iron.Draw("HIST  C SAME")
lead.Draw("HIST C SAME")
carbon.Draw("HIST C SAME")
#water.Draw("HIST SAME")

tracker.SetMinimum(0.7)
tracker.SetMaximum(1.3)

#mnv.AddHistoTitle("Neutrino Flux", 0.05, 1)
mnv.AddHistoTitle("Antineutrino Flux", 0.05, 1)
gStyle.SetErrorX(0)

legend = TLegend(0.35,0.75,0.82,0.89)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextSize(0.04)
legend.SetNColumns(2)
#legend.AddEntry(reco, " Simulation", "l")
#legend.AddEntry(data_hist, " Data", "ep")
legend.AddEntry(tracker, " Tracker (C) ", "l")
legend.AddEntry(iron, " Iron", "l")
legend.AddEntry(lead, " Lead", "l")
legend.AddEntry(carbon, " Carbon", "l")
#legend.AddEntry(water, " Water (O)", "l")
legend.SetTextFont(42)
legend.Draw()

#mnv.AddPlotLabel("MC POT "+ "{:.2e}".format(mcPOT), 0.29, 0.20, 0.04, 12, 42)
#mnv.AddPlotLabel("MC POT "+ "{:.2e}".format(mcPOT) + " (Water "+"{:.2e}".format(mcPOTwater) + ")", 0.38, 0.87, 0.033, 12, 42)
#mnv.AddPlotLabel("MC POT "+ "{:.2e}".format(mcPOT) , 0.28, 0.87, 0.033, 12, 42)


canvas1.Modified()
canvas1.Print("AllAntiNu_Fluxes_sys_ratio_close_up.png")


raw_input("Done")
