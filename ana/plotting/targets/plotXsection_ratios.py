import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLine


infile= ROOT.TFile("CrossSection_t3_z26_minervame6A.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale =  dataPOT/mcPOT

vars = ["Enu", "x"]

steps = ["unfolded_effCorrected_3"]#, "total_unfolded_effCorrected"]
#steps = ['crossSection', 'crossSection_total']

for step in steps:

    for var in vars:

        data_hist = infile.Get("%s_data_%s"%(step, var))

        if step == "unfolded":
            mc_hist = infile.Get("%s_mc_%s"%(step, var))
        elif step == "unfolded_effCorrected":
            mc_hist = infile.Get("%s_mc_%s"%(step, var))
        elif step == "unfolded_effCorrected_3":
            mc_hist = infile.Get("%s_mc_%s"%(step, var))
        elif step == "total_unfolded_effCorrected":
            mc_hist = infile.Get("%s_mc_%s"%(step, var))
        else:
            mc_hist = infile.Get("simEventRate_%s_mc_%s"%(step, var))
    

        if var == "Enu":
            data_hist.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
            if step == "crossSection":
                data_hist.GetYaxis().SetTitle("d#sigma/dE_{#bar{#nu}} (10^{-39} cm^{2}/GeV/CH)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                data_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                data_hist.GetYaxis().SetTitle("Events/GeV")

        if var == "x":
            data_hist.GetXaxis().SetTitle("Bjorken x")
            if step == "crossSection":
                data_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/GeV/CH)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                data_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                data_hist.GetYaxis().SetTitle("Events (norm.)")

        # Data/MC ratio UPDATED
        if step == "unfolded":
            mc_hist.Scale(mcScale)
        elif step == "unfolded_effCorrected":
            mc_hist.Scale(mcScale)
        else:
            mc_hist.Scale(1)
  
        gStyle.SetErrorX(0)
        data_hist_stat =  data_hist.GetCVHistoWithStatError().Clone() # stat error
        data_hist_total = data_hist.GetCVHistoWithError() # total error
        data_hist_sys = data_hist.GetCVHistoWithError(True) # sys error
        # in cross-section, we ALWAYS SHOW ONLY STAT ERROR ON MC
        mc_hist_stat = mc_hist.GetCVHistoWithStatError()

        ratio = data_hist_stat.Clone()
        ratio.Divide(ratio,mc_hist_stat) # stat
        ratio.GetYaxis().SetTitle("Data/MC")
        ratio.GetYaxis().CenterTitle()
        ratio.GetXaxis().CenterTitle()
        ratio.Draw("X0")
        #ratio.SetMaximum(1.5)
        #ratio.SetMinimum(0.5)
        ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions
        # Copied from MINERvA-101-Cross-Section/backgroundStack.py
        # same as in void MnvPlotter::DrawDataMCRatio() in MnvPlotter
        # systematic error centered at y = 1
        ratio_tot = data_hist_total.Clone()
        ratio_tot.Divide(ratio_tot, mc_hist_stat)
        #sysError = data_hist_sys.Clone()
        #mc_hist.GetTotalError(False, True, False) # False for stat error, True for as frac, False for covAreaNorm
        #for whichBin in range(1, sysError.GetXaxis().GetNbins()+1):
        #    sysError.SetBinError(whichBin, max(sysError.GetBinContent(whichBin), 1e-9))
        #    sysError.SetBinContent(whichBin, 1)

        #sysError.SetFillColor(ROOT.kRed-10)
        #sysError.SetFillStyle(1001)
        #sysError.SetMarkerStyle(0)
        #sysError.Draw("E2 SAME")
        #sysError.Draw("E4 SAME")
        
        # line at 1
        line = TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
        line.SetLineColor(46)
        line.SetLineWidth(2)
        line.SetLineStyle(9)
        line.Draw()
        
        ratio.Draw("X0 SAME E1")
        ratio_tot.Draw("E1 SAME")
    
        '''
        Stat error for data and sys error for MC (for pre-background subtraction). 
        After background subtraction both stat+sys error on data and sys error on MC"
        '''
        #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
        mnv.AddHistoTitle("T3 Iron %s"%step, 0.05, 1)
        if step == "unfolded":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        elif step == "unfolded_effCorrected":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)

        #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
        if var == "ANNPlaneProb":
            canvas1.SetLogy()
        canvas1.Modified()
        canvas1.Print("ME6A_T3Iron_%s_%s_ratio.png"%(step, var))
    
        # ----------------------------------------------------------------------------
        # Event selection with simulated background

        #mnv.DrawDataMCWithErrorBand(data_hist, mc_hist, mcScale, "TR", False, bkg_mc, None, False, True, False)
        #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
        #mnv.AddHistoTitle("Tracker", 0.05, 1)
        #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        #canvas1.Print("ME6A_Tracker_EventSelection_SimBkg_%s_sys.png"%(var))

        '''
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
        
        #mnv2.ApplyStyle(7)
        mnv2.DrawErrorSummary(mc_hist, "TL", True, True, 0.0, False, "",True);
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

        mnv2.AddHistoTitle("Tracker %s MC"%step, 0.05, 1)

        canvas1.Modified()
        #canvas1.Print("ME6A_Tracker_EventSelection_FracErr_%s_sys.png"%(var))

        # ---------------------------------------------
        # Data fractional error
        mnv2.DrawErrorSummary(data_hist, "TL", True, True, 0.0, False, "",True);
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
                    k.SetY1(0.27) #Enu

                if var == "x":
                    k.SetNColumns(2)
                    k.SetX2(2.0) #x
                    k.SetY1(0.16) #x

        mnv2.AddHistoTitle("Tracker %s Data"%step, 0.05, 1)

        canvas1.Modified()
        #canvas1.Print("ME6A_Tracker_BackgroundEventSelection_FracErr_%s_sys.png"%(var))
        '''

canvas2 = ROOT.TCanvas()
flux = infile.Get("fluxRebinned_Enu")
flux.Draw("HIST")
flux.GetXaxis().CenterTitle()
flux.GetYaxis().CenterTitle()
flux.GetYaxis().CenterTitle(True)
flux.GetYaxis().SetTitleOffset(1.3)
flux.GetYaxis().SetTitleSize(0.05)
flux.GetYaxis().SetLabelSize(0.05)
flux.GetXaxis().CenterTitle(True)
#flux.GetXaxis().SetTitle("Neutrino Energy (GeV)")
flux.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
flux.GetXaxis().SetTitleSize(0.05)
flux.GetXaxis().SetLabelSize(0.05)
flux.GetYaxis().SetTitle("#bar{#nu} / m^{2} / P.O.T/ GeV")
mnv.AddHistoTitle("ME6A Flux", 0.05, 1)
canvas2.Print("ME6A_FluxRebinned.png")

        

data_hist.SetDirectory(0)
mc_hist.SetDirectory(0)

raw_input("Done")
