import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLegend


infile= ROOT.TFile("CrossSection_t99_z99_minervame6A_NEW.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale =  dataPOT/mcPOT

vars = ["Enu", "x"]

steps = ['unfolded','unfolded_effCorrected', 'crossSection', 'crossSection_total']
#steps = ['crossSection', 'crossSection_total']

for step in steps:

    for var in vars:

        data_hist = infile.Get("%s_data_%s"%(step, var))

        if step == "unfolded":
            mc_hist = infile.Get("efficiency_numerator_%s"%(var))
        elif step == "unfolded_effCorrected":
            mc_hist = infile.Get("simEventRate_%s"%(var))
        else:
            mc_hist = infile.Get("simEventRate_%s_mc_%s"%(step, var))

        if var == "Enu":
            mc_hist.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dE_{#bar{#nu}} (10^{-39} cm^{2}/GeV/CH)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events/GeV")

        if var == "x":
            mc_hist.GetXaxis().SetTitle("Bjorken x")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/CH)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events (norm.)")

        if step == "crossSection":
            mc_hist.Scale(1E39)
            data_hist.Scale(1E39)
        
        if step == "crossSection_total":
            mc_hist.Scale(1E38)
            data_hist.Scale(1E38)

        mc_hist.GetXaxis().CenterTitle()
        mc_hist.GetYaxis().CenterTitle()
        #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
        #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
        #mc_hist.GetXaxis().SetRangeUser(66,173)

        #mnv.ApplyStyle(1)
        if step == "unfolded":
            mc_hist.Scale(mcScale)
        elif step == "unfolded_effCorrected":
            mc_hist.Scale(mcScale)
        else:
            mc_hist.Scale(1)
  
        mc_hist.Scale(mc_hist.GetNormBinWidth(), "width")
        data_hist.Scale(mc_hist.GetNormBinWidth(), "width")
        data_hist_stat =  data_hist.GetCVHistoWithStatError() # stat error
        data_hist_total = data_hist.GetCVHistoWithError() # total error
        data_hist_sys = data_hist.GetCVHistoWithError(True) # sys error
        mc_hist_stat = mc_hist.GetCVHistoWithStatError() 
        
        # MC
        mc_hist.SetLineWidth(3)
        mc_hist.SetLineColor(2)
        # MC error
        mc_hist_stat.SetFillColor(ROOT.kRed-10)
        mc_hist_stat.SetFillStyle(1001)
        mc_hist_stat.SetMarkerStyle(0)
        mc_hist_stat.SetLineWidth(3)
        mc_hist_stat.SetLineColor(2)

        mc_hist.SetLineColor(ROOT.kRed)
        mc_hist.SetLineWidth(2)
        mc_hist.SetMaximum(data_hist.GetMaximum()*1.5)

        #gStyle.SetErrorX()
        mc_hist.Draw("HIST")
        #gStyle.SetErrorX(0.5)
        mc_hist_stat.Draw("E2 SAME")
        #gStyle.SetErrorX(0)
        mc_hist.Draw("HIST SAME")
        data_hist.Draw("SAME E1 X0")
        #data_hist_stat.Draw("E2 SAME")
        data_hist_total.Draw("E1 SAME X0")

        ''' void MnvPlotter::DrawDataMCWithErrorBand(
        const MnvH1D* dataHist,
        const MnvH1D* mcHist,
        const Double_t mcScale /*= 1.0*/,
        const std::string& legPos /*= "L"*/,
        const bool useHistTitles /*=false*/,
        const MnvH1D* bkgdHist /*= NULL*/,
        const MnvH1D* dataBkgdHist /*= NULL*/,
        const bool covAreaNormalize/*= false, Area Normalize considerations for covariance matrix*/,
        const bool statPlusSys /* = false */,
        const bool isSmooth /* = false */)
        '''
        #if step == 'crossSection' or 'crossSection_total':
        #   chi2 = mnv.Chi2DataMC(data_hist, mc_hist,mcScale,False,True )
        #    ndf = mc_hist.GetNbinsX()

        #print(chi2)
        if step == "unfolded":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        elif step == "unfolded_effCorrected":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)

        '''
        Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
        If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
        If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
        '''
        if step == 'crossSection':
            mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
        
        if step == 'crossSection_total':
            mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
        
        mnv.AddHistoTitle("Tracker %s"%step, 0.05, 1)

        mc_hist.GetYaxis().SetTitleOffset(0.99)

        legend = TLegend(0.55,0.7,0.80,0.89)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.035)
        legend.AddEntry(mc_hist_stat, " Simulation", "fl")
        legend.AddEntry(data_hist, " Data", "lep")
        #legend.SetTextFont(42)
        legend.Draw()

        canvas1.Modified()
        canvas1.Print("ME6A_Tracker_%s_%s.png"%(step,var))

        ''''
        # ----------------------------------------------------------------------------
        # Data/MC ratio
        if step == "unfolded":
            mnv.DrawDataMCRatio(data_hist, mc_hist, mcScale)
        elif step == "unfolded_effCorrected":
            mnv.DrawDataMCRatio(data_hist, mc_hist, mcScale)
        else:
            mnv.DrawDataMCRatio(data_hist, mc_hist, 1)
        
        Stat error for data and sys error for MC (for pre-background subtraction). 
        After background subtraction both stat+sys error on data and sys error on MC"
        
        #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
        mnv.AddHistoTitle("Tracker %s"%step, 0.05, 1)
        if step == "unfolded":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        elif step == "unfolded_effCorrected":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)

        #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
        if var == "ANNPlaneProb":
            canvas1.SetLogy()
        canvas1.Modified()
        canvas1.Print("ME6A_Tracker_%s_%s_ratio.png"%(step, var))
        '''
        # ----------------------------------------------------------------------------
        # Event selection with simulated background

        #mnv.DrawDataMCWithErrorBand(data_hist, mc_hist, mcScale, "TR", False, bkg_mc, None, False, True, False)
        #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
        #mnv.AddHistoTitle("Tracker", 0.05, 1)
        #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        #canvas1.Print("ME6A_Tracker_EventSelection_SimBkg_%s_sys.png"%(var))

        
        # ----------------------------------------------------------------------------
        # Fractional Error Groups
        mnv2 = PlotUtils.MnvPlotter()
        mnv2.ApplyStyle(7)

        mnv2.error_summary_group_map.clear()
        mnv2.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_Resolution")
        mnv2.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_MINOS")
        mnv2.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_MINERvA")

        mnv2.error_summary_group_map["Muon Efficiency"].push_back("MINOS_Reconstruction_Efficiency")

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


        mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_C")
        mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_Fe")
        mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_Pb")
        mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_CH")
        mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_H2O")

        '''
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
        '''
        # ---------------------------------------------
        # Data fractional error
        data_hist.GetXaxis().CenterTitle()
        data_hist.GetYaxis().CenterTitle()

        if var == "Enu":
            data_hist.GetXaxis().SetTitle("Antineutrino Energy (GeV)")

        if var == "x":
            data_hist.GetXaxis().SetTitle("Bjorken x")

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
                    k.SetY1(0.3) #Enu

                if var == "x":
                    k.SetNColumns(2)
                    k.SetX2(2.0) #x
                    k.SetY1(0.065) #x

        mnv2.AddHistoTitle("Tracker %s (Data)"%step, 0.04, 1)

        canvas1.Modified()
        canvas1.Print("ME6A_Tracker_%s_FracErr_%s_sys.png"%(step, var))
        

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
