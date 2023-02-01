import ROOT
import os,sys
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
from ROOT import PlotUtils

#targetID = sys.argv[1] 
#targetZ = sys.argv[2]

infile= ROOT.TFile("Hists_EventSelectionTracker_ML_ME6A_sys_t99_z99_AntiNu.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

ROOT.TH1.AddDirectory(False)

mcScale =  dataPOT/mcPOT

vars = ["Enu", "x"]

for var in vars:

    # numerator vs denominator
    num_hist = infile.Get("selected_mc_reco_signal_%s"%var) # selected signal
    denom_hist = infile.Get("selected_mc_reco_%s"%var) # total selected events
    denom_hist.SetLineColor(ROOT.kRed)

    ratio = num_hist.Clone()
    ratio.Divide(num_hist,denom_hist, 1.0, 1.0, "B") # binomial because numerator is the subset of the denominator, to make the stat error correct
    sysError = ratio.Clone()
    ratio.SetLineColor(ROOT.kRed)

    if var == "Enu":
        ratio.GetXaxis().SetTitle("Reconstructed Neutrino Energy (GeV)")
        denom_hist.GetXaxis().SetTitle("Reconstructed Neutrino Energy (GeV)")
        denom_hist.GetYaxis().SetTitle("Events/GeV")
        num_hist.GetXaxis().SetTitle("Reconstructed Neutrino Energy (GeV)")
        num_hist.GetYaxis().SetTitle("Events/GeV")

    if var == "x":
        ratio.GetXaxis().SetTitle("Reconstructed Bjorken x")
        denom_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        denom_hist.GetYaxis().SetTitle("Events (norm.)")
        num_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        num_hist.GetYaxis().SetTitle("Events (norm.)")

    denom_hist.GetXaxis().CenterTitle()
    denom_hist.GetYaxis().CenterTitle()
    
    denom_hist.Draw("HIST")
    sysErr_denom = denom_hist.Clone()
    sysErr_denom.SetFillColor(ROOT.kRed-10)
    sysErr_denom.SetFillStyle(1001)
    sysErr_denom.SetMarkerStyle(0)
    sysErr_denom.SetLineWidth(3)
    sysErr_denom.SetLineColor(2)
    sysErr_denom.Draw("E2 SAME")
    denom_hist.Draw("HIST SAME")

    num_hist.Draw("HIST SAME")
    sysErr_num = num_hist.Clone()
    sysErr_num.SetFillColor(ROOT.kGray)
    sysErr_num.SetFillStyle(1001)
    sysErr_num.SetMarkerStyle(0)
    sysErr_num.SetLineWidth(3)
    sysErr_num.SetLineColor(ROOT.kGray+2)
    sysErr_num.Draw("E2 SAME")
    num_hist.Draw("HIST SAME")

    legend = TLegend(0.56,0.70,0.80,0.89)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    #legend.AddEntry(reco, " Simulation", "l")
    legend.AddEntry(sysErr_num, " Numerator", "fl")
    legend.AddEntry(sysErr_denom, " Denominator", "fl")
    legend.Draw()

    mnv.AddHistoTitle("Purity: Tracker", 0.05, 1)
    canvas1.Modified()
    canvas1.Update()

    canvas1.Print("ME6A_Tracker_Purity_NumDenum_%s_sys.png"%var)

    # ----------------------------------------------------------------------------
    # Ratio

    ratio.GetYaxis().SetTitle("Purity")
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

    mnv.AddHistoTitle("Tracker", 0.05, 1)

    canvas1.Modified()
    canvas1.Update()
    canvas1.Print("ME6A_Tracker_Purity_ratio_%s_sys.png"%var)

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
    mnv2.error_summary_group_map["Interaction Model"].push_back("GENIE_MaRES")
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


    mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass")
    
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
                k.SetY1(0.2) #Enu

            if var == "x":
                k.SetNColumns(2)
                k.SetX2(2.0) #x
                k.SetY1(0.16) #x

    mnv2.AddHistoTitle("Purity Numerator: Tracker", 0.04, 1)

    canvas1.Modified()
    canvas1.Print("ME6A_Tracker_Purity_FracErrNum_%s_sys.png"%(var))

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
                k.SetY1(0.2) #Enu

            if var == "x":
                k.SetNColumns(2)
                k.SetX2(2.0) #x
                k.SetY1(0.16) #x

    mnv2.AddHistoTitle("Purity Denominator: Tracker", 0.04, 1)

    canvas1.Modified()
    canvas1.Print("ME6A_Tracker_Purity_FracErrDenom_%s_sys.png"%(var))

    # Ratio
    mnv2.DrawErrorSummary(ratio, "TL", True, True, 0.0, False, "",True)
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
                k.SetY1(0.02) #Enu

            if var == "x":
                k.SetNColumns(2)
                k.SetX2(2.0) #x
                k.SetY1(0.0012) #x

    mnv2.AddHistoTitle("Purity Tracker", 0.04, 1)

    canvas1.Modified()
    canvas1.Print("ME6A_Tracker_Purity_FracErr_%s_sys.png"%(var))


raw_input("Done")
