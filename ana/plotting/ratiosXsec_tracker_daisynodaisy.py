import ROOT
import os,sys
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
from ROOT import PlotUtils

ROOT.gROOT.SetBatch(True)

nodaisy = ROOT.TFile("rootFiles/CrossSection_t99_z99_minervame6A.root")
infile = ROOT.TFile("rootFiles/CrossSection_Daisy_t99_z99_minervame6A.root")

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

ROOT.TH1.AddDirectory(False)

vars = ["Enu","x"]

for var in vars:
    if var=="Enu":
        data_hist_C = infile.Get("crossSection_total_carbon_data_%s"%(var))
        data_hist_Fe = infile.Get("crossSection_total_iron_data_%s"%(var))
        data_hist_Pb = infile.Get("crossSection_total_lead_data_%s"%(var))

        data_hist = nodaisy.Get("crossSection_total_data_%s"%(var))

        mc_hist_C = infile.Get("simEventRate_crossSection_total_carbon_mc_%s"%(var))
        mc_hist_Fe = infile.Get("simEventRate_crossSection_total_iron_mc_%s"%(var))
        mc_hist_Pb = infile.Get("simEventRate_crossSection_total_lead_mc_%s"%(var))

        mc_hist = nodaisy.Get("simEventRate_crossSection_total_mc_%s"%(var))

    else: 
        data_hist_C = infile.Get("crossSection_carbon_data_%s"%(var))
        data_hist_Fe = infile.Get("crossSection_iron_data_%s"%(var))
        data_hist_Pb = infile.Get("crossSection_lead_data_%s"%(var))

        data_hist = nodaisy.Get("crossSection_data_%s"%(var))

        mc_hist_C = infile.Get("simEventRate_crossSection_carbon_mc_%s"%(var))
        mc_hist_Fe = infile.Get("simEventRate_crossSection_iron_mc_%s"%(var))
        mc_hist_Pb = infile.Get("simEventRate_crossSection_lead_mc_%s"%(var))

        mc_hist = nodaisy.Get("simEventRate_crossSection_mc_%s"%(var))


    ratio_data_C =  data_hist_C.Clone()
    ratio_data_C.Divide(data_hist_C,data_hist)

    ratio_data_Fe =  data_hist_Fe.Clone()
    ratio_data_Fe.Divide(data_hist_Fe,data_hist)

    ratio_data_Pb =  data_hist_Pb.Clone()
    ratio_data_Pb.Divide(data_hist_Pb,data_hist)

    ratio_mc_C =  mc_hist_C.Clone()
    ratio_mc_C.Divide(mc_hist_C,mc_hist)

    ratio_mc_Fe =  mc_hist_Fe.Clone()
    ratio_mc_Fe.Divide(mc_hist_Fe,mc_hist)

    ratio_mc_Pb =  mc_hist_Pb.Clone()
    ratio_mc_Pb.Divide(mc_hist_Pb,mc_hist)

    if var == "Enu":
        ratio_mc_C.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
        ratio_mc_C.GetYaxis().SetTitle("Events/GeV")
        ratio_data_C.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
        ratio_data_C.GetYaxis().SetTitle("Events/GeV")

    if var == "x":
        ratio_mc_C.GetXaxis().SetTitle("Bjorken x")
        ratio_mc_C.GetYaxis().SetTitle("Events (norm.)")
        ratio_data_C.GetXaxis().SetTitle("Bjorken x")
        ratio_data_C.GetYaxis().SetTitle("Events (norm.)")

    ratio_mc_C.GetXaxis().CenterTitle()
    ratio_mc_C.GetYaxis().CenterTitle()
    ratio_mc_C.GetYaxis().SetNdivisions(505)
    ratio_mc_C.GetXaxis().SetTitleSize(0.05)
    ratio_mc_C.GetYaxis().SetTitleSize(0.05)
    ratio_mc_C.GetYaxis().SetTitleOffset(1.6)

    # ----------------------------------------------------------------------------
    # Ratio MC

    ratio_mc_C.GetYaxis().SetTitle("xsec daisy tracker/tracker")
    ratio_mc_C.GetXaxis().CenterTitle()
    ratio_mc_C.GetYaxis().CenterTitle()

    ratio_mc_C.SetLineColor(46)
    ratio_mc_Fe.SetLineColor(38)
    ratio_mc_Pb.SetLineColor(30)

    ratio_mc_stat_C = ratio_mc_C.GetCVHistoWithStatError() 
    ratio_mc_stat_Fe = ratio_mc_Fe.GetCVHistoWithStatError() 
    ratio_mc_stat_Pb = ratio_mc_Pb.GetCVHistoWithStatError() 

    ratio_mc_stat_C.SetMarkerColor(46)
    ratio_mc_stat_C.SetFillColor(46)
    ratio_mc_stat_C.SetFillStyle(3003)
    ratio_mc_stat_Fe.SetMarkerColor(38)
    ratio_mc_stat_Fe.SetFillColor(38)
    ratio_mc_stat_Fe.SetFillStyle(3003)
    ratio_mc_stat_Pb.SetMarkerColor(38)
    ratio_mc_stat_Pb.SetFillColor(30)
    ratio_mc_stat_Pb.SetFillStyle(3003)


    ratio_mc_C.Draw("HIST")
    ratio_mc_stat_C.Draw("E2 SAME")
    ratio_mc_Fe.Draw("HIST SAME")
    ratio_mc_stat_Fe.Draw("E2 SAME")
    ratio_mc_Pb.Draw("HIST SAME")
    ratio_mc_stat_Pb.Draw("E2 SAME")

    #ratio_data_C.Draw("HIST SAME")
    #ratio_data_Fe.Draw("HIST SAME")
    #ratio_data_Pb.Draw("HIST SAME")

    ratio_mc_C.SetMaximum(1.1)
    ratio_mc_C.SetMinimum(0.9)

    gStyle.SetOptTitle(0)

    legend = TLegend(0.20,0.7,0.40,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(ratio_mc_stat_C, " Simulation C", "fl")
    legend.AddEntry(ratio_mc_stat_Fe, " Simulation Fe", "fl")
    legend.AddEntry(ratio_mc_stat_Pb, " Simulation Pb", "fl")
    legend.SetTextFont(42)
    legend.Draw()

    mnv.AddHistoTitle("MC xsec daisy tracker/tracker" , 0.05, 1)

    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("ME6A_xSec_ratio_MC_daisy-nodaisy_%s.png"%(var))

    # ----------------------------------------------------------------------------
    ## Ratio data

    ratio_data_C.GetYaxis().SetTitle("xsec daisy tracker/tracker")
    ratio_data_C.GetXaxis().CenterTitle()
    ratio_data_C.GetYaxis().CenterTitle()

    ratio_data_C.SetLineColor(46)
    ratio_data_Fe.SetLineColor(38)
    ratio_data_Pb.SetLineColor(30)

    ratio_data_stat_C = ratio_data_C.GetCVHistoWithStatError() 
    ratio_data_stat_Fe = ratio_data_Fe.GetCVHistoWithStatError() 
    ratio_data_stat_Pb = ratio_data_Pb.GetCVHistoWithStatError() 

    ratio_data_total_C = ratio_data_C.GetCVHistoWithError() # total error
    ratio_data_total_Fe = ratio_data_Fe.GetCVHistoWithError() # total error
    ratio_data_total_Pb = ratio_data_Pb.GetCVHistoWithError() # total error 

    ratio_data_stat_C.SetMarkerColor(46)
    ratio_data_stat_C.SetFillColor(46)
    ratio_data_stat_C.SetFillStyle(3003)
    ratio_data_stat_Fe.SetMarkerColor(38)
    ratio_data_stat_Fe.SetFillColor(38)
    ratio_data_stat_Fe.SetFillStyle(3003)
    ratio_data_stat_Pb.SetMarkerColor(38)
    ratio_data_stat_Pb.SetFillColor(30)
    ratio_data_stat_Pb.SetFillStyle(3003)

    ratio_data_total_C.SetMarkerColor(46)
    ratio_data_total_C.SetFillColor(46)
    ratio_data_total_C.SetFillStyle(3003)
    ratio_data_total_Fe.SetMarkerColor(38)
    ratio_data_total_Fe.SetFillColor(38)
    ratio_data_total_Fe.SetFillStyle(3003)
    ratio_data_total_Pb.SetMarkerColor(38)
    ratio_data_total_Pb.SetFillColor(30)
    ratio_data_total_Pb.SetFillStyle(3003)


    ratio_data_C.Draw("HIST")
    ratio_data_total_C.Draw("E2 SAME")
    ratio_data_Fe.Draw("HIST SAME")
    ratio_data_total_Fe.Draw("E2 SAME")
    ratio_data_Pb.Draw("HIST SAME")
    ratio_data_total_Pb.Draw("E2 SAME")

    #ratio_data_C.Draw("HIST SAME")
    #ratio_data_Fe.Draw("HIST SAME")
    #ratio_data_Pb.Draw("HIST SAME")

    ratio_data_C.SetMaximum(1.1)
    ratio_data_C.SetMinimum(0.9)

    gStyle.SetOptTitle(0)

    legend = TLegend(0.20,0.7,0.40,0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(ratio_data_C, " Data C", "fl")
    legend.AddEntry(ratio_data_Fe, " Data Fe", "fl")
    legend.AddEntry(ratio_data_Pb, " Data Pb", "fl")
    legend.SetTextFont(42)
    legend.Draw()

    mnv.AddHistoTitle("Data xsec daisy tracker/tracker" , 0.05, 1)

    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("ME6A_xSec_ratio_Data_daisy-nodaisy_%s_v2.png"%(var))


raw_input("Done")
