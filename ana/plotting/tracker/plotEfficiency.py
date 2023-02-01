import ROOT
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis, gPad, TLine
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
import os,sys
from ROOT import PlotUtils


infile = ROOT.TFile("Hists_Efficiency_ML_ME6A_sys_t99_z99_AntiNu.root")

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

ROOT.TH1.AddDirectory(False)
mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

vars = ["Enu", "x"]
for var in vars:
    # numerator vs denominator
    num_hist = infile.Get("h_mc_%s"%var) # selected signal
    denom_hist = infile.Get("h_truth_%s"%var) # total signal
    num = num_hist.GetEntries()
    denom = denom_hist.GetEntries()

    num_hist.Scale(0.5)
    denom_hist.Scale(0.5)
    print("Numerator: " + str(num))
    print("Denominator: " + str(denom))
    print("Efficinecy: " + str(num/denom))

    ratio = num_hist.Clone()
    ratio.Divide(num_hist,denom_hist, 1.0, 1.0, "B") # binomial because numerator is the subset of the denominator, to make the stat error correct
    sysError = ratio.Clone()
    ratio.SetLineColor(ROOT.kRed)

    if var == "Enu":
        ratio.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
        denom_hist.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
        denom_hist.GetYaxis().SetTitle("Events/GeV")
        num_hist.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
        num_hist.GetYaxis().SetTitle("Events/GeV")
    if var == "x":   
        ratio.GetXaxis().SetTitle("Bjorken x")
        denom_hist.GetXaxis().SetTitle("Bjorken x")
        denom_hist.GetYaxis().SetTitle("Events (norm.)")
        num_hist.GetXaxis().SetTitle("Bjorken x")
        num_hist.GetYaxis().SetTitle("Events (norm.)")

    ratio.GetYaxis().SetTitle("Efficiency")
    ratio.GetXaxis().CenterTitle()
    ratio.GetYaxis().CenterTitle()
    ratio.Draw("HIST")

    sysError.SetFillColor(ROOT.kRed-10)
    sysError.SetFillStyle(1001)
    sysError.SetMarkerStyle(0)
    sysError.Draw("E2 SAME")
    sysError.SetMaximum(1)

    ratio.Draw("HIST SAME")
    ratio.SetMaximum(1)

    gStyle.SetOptTitle(0)

    legend = TLegend(0.56,0.70,0.80,0.89)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)

    mnv.AddHistoTitle("Target Tracker", 0.05, 1)

    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("ME6A_Tracker_Eficiency_ratio_%s_sys.png"%(var))

    # Numerator and Denominator

    num_hist.Scale(num_hist.GetBinWidth(1), "width")
    denom_hist.Scale(denom_hist.GetBinWidth(1), "width")

    denom_hist.SetLineColor(ROOT.kRed+2)
    denom_hist.Draw("HIST")
    sysErr_denom = denom_hist.Clone()
    sysErr_denom.SetFillColor(ROOT.kRed-10)
    sysErr_denom.SetFillStyle(1001)
    sysErr_denom.SetMarkerStyle(0)
    sysErr_denom.SetLineWidth(3)
    sysErr_denom.SetLineColor(2)
    sysErr_denom.Draw("E2 SAME")
    denom_hist.Draw("HIST SAME")

    num_hist.SetLineColor(ROOT.kGray+2)
    num_hist.Draw("SAME HIST")
    sysErr_num = num_hist.Clone()
    sysErr_num.SetFillColor(ROOT.kGray)
    sysErr_num.SetFillStyle(1001)
    sysErr_num.SetMarkerStyle(0)
    sysErr_num.SetLineWidth(3)
    sysErr_num.SetLineColor(ROOT.kGray+2)
    sysErr_num.Draw("E2 SAME")
    num_hist.Draw("HIST SAME")

    denom_hist.GetXaxis().CenterTitle()
    denom_hist.GetYaxis().CenterTitle()

    legend = TLegend(0.56,0.70,0.80,0.89)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    #legend.AddEntry(reco, " Simulation", "l")
    legend.AddEntry(sysErr_num, " Numerator", "fl")
    legend.AddEntry(sysErr_denom, " Denominator", "fl")
    legend.Draw()

    mnv.AddHistoTitle("Efficiency: Tracker", 0.05, 1)

    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("ME6A_Tracker_Eficiency_NumDenum_%s_sys.png"%(var))

    # ----------------------------------------------------------------------------
    # Fractional Error Groups
    mnv2 = PlotUtils.MnvPlotter()
    mnv2.ApplyStyle(7)

    mnv2.error_summary_group_map.clear()
    mnv2.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_Resolution")
    mnv2.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_MINOS")
    mnv2.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_MINERvA")

    mnv2.error_summary_group_map["MINOS"].push_back("MINOS_Reconstruction_Efficiency")

    mnv2.error_summary_group_map["Muon Angle"].push_back("BeamAngleX")
    mnv2.error_summary_group_map["Muon Angle"].push_back("BeamAngleY")

    mnv2.error_summary_group_map["Flux"].push_back("Flux")

    
    mnv2.error_summary_group_map["Interaction Model"].push_back("RPA_HighQ2")
    mnv2.error_summary_group_map["Interaction Model"].push_back("RPA_LowQ2")
    mnv2.error_summary_group_map["Interaction Model"].push_back("Low_Recoil_2p2h_Tune")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_AhtBY")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_BhtBY")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_CCQEPauliSupViaKF")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_CV1uBY")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_CV2uBY")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_EtaNCEL")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_MaCCQE")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_MaNCEL")
    #mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_MaRES")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_MvRES")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_NormDISCC")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_NormNCRES")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvn1pi")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvn2pi")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvp1pi")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvp2pi")
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_VecFFCCQEshape")


    mnv2.error_summary_group_map["FSI"].push_back("GENIE_AGKYxF1pi")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_RDecBR1gamma")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_Theta_Delta2Npi")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_FrAbs_N")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_FrAbs_pi")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_FrCEx_N")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_FrCEx_pi")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_FrElas_N")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_FrElas_pi")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_FrInel_N")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_FrPiProd_N")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_FrPiProd_pi")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_MFP_N")
    mnv2.error_summary_group_map["FSI"].push_back("GENIE_MFP_pi")


    mnv2.error_summary_group_map["Particle Response"].push_back("response_em")
    mnv2.error_summary_group_map["Particle Response"].push_back("response_low_proton")
    mnv2.error_summary_group_map["Particle Response"].push_back("response_mid_proton")
    mnv2.error_summary_group_map["Particle Response"].push_back("response_high_proton")
    mnv2.error_summary_group_map["Particle Response"].push_back("response_meson")
    mnv2.error_summary_group_map["Particle Response"].push_back("response_other")

    mnv2.error_summary_group_map["GEANT4 Hadrons"].push_back("GEANT_Neutron")
    mnv2.error_summary_group_map["GEANT4 Hadrons"].push_back("GEANT_Proton")
    mnv2.error_summary_group_map["GEANT4 Hadrons"].push_back("GEANT_Pion")

    mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_C")
    mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_Fe")
    mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_Pb")
    mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_CH")
    mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_H2O")

    # Change colour
    mnv2.error_color_map["Target Mass"] = 15
    mnv2.error_color_map["Muon Efficiency"] = 46
    mnv2.error_color_map["Particle Response"] = 65
    
    # Numerator
    mnv2.DrawErrorSummary(num_hist, "TL", True, True, 0.0, False, "",True);
    # last boolean decides whether frac or not

    keys = canvas1.GetListOfPrimitives();
    for k in keys:
        if(k.ClassName().find("TH1")!=-1):
            if var == "Enu":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.02) # to 0.35 for Enu
            if var == "x":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.05)

        if(k.ClassName().find("Legend")!=-1):
            if var == "Enu":
                k.SetNColumns(2)
                k.SetX2(45) #Enu
                k.SetY1(0.25) #Enu

            if var == "x":
                k.SetNColumns(2)
                k.SetX2(2.0) #x
                k.SetY1(0.23) #x

    mnv2.AddHistoTitle("Efficiency Numerator: Tracker", 0.04, 1)

    canvas1.Modified()
    canvas1.Print("ME6A_Tracker_Efficiency_FracErrNum_%s_sys.png"%(var))

    # Denominator
    mnv2.DrawErrorSummary(denom_hist, "TL", True, True, 0.0, False, "",True);
    # last boolean decides whether frac or not

    keys = canvas1.GetListOfPrimitives();
    for k in keys:
        if(k.ClassName().find("TH1")!=-1):
            if var == "Enu":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.02) # to 0.35 for Enu
            if var == "x":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.05)

        if(k.ClassName().find("Legend")!=-1):
            if var == "Enu":
                k.SetNColumns(2)
                k.SetX2(45) #Enu
                k.SetY1(0.25) #Enu

            if var == "x":
                k.SetNColumns(2)
                k.SetX2(2.0) #x
                k.SetY1(0.23) #x

    mnv2.AddHistoTitle("Efficiency Denominator: Tracker", 0.04, 1)

    canvas1.Modified()
    canvas1.Print("ME6A_Tracker_EventSelection_Efficiency_FracErrDenom_%s_sys.png"%(var))

    # Ratio
    mnv2.DrawErrorSummary(ratio, "TL", True, True, 0.0, False, "",True);
    # last boolean decides whether frac or not

    keys = canvas1.GetListOfPrimitives();
    for k in keys:
        if(k.ClassName().find("TH1")!=-1):
            if var == "Enu":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.02) # to 0.35 for Enu
            if var == "x":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.05)

        if(k.ClassName().find("Legend")!=-1):
            if var == "Enu":
                k.SetNColumns(2)
                k.SetX2(45) #Enu
                k.SetY1(0.025) #Enu

            if var == "x":
                k.SetNColumns(2)
                k.SetX2(2.0) #x
                k.SetY1(0.024) #x

    mnv2.AddHistoTitle("Efficiency Tracker", 0.04, 1)

    canvas1.Modified()
    canvas1.Print("ME6A_Tracker_EventSelection_Efficiency_FracErr_%s_sys.png"%(var))




raw_input("Done")
