import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
from ROOT import TLegend

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
targetID = sys.argv[2] 
targetZ = sys.argv[3]
plist = sys.argv[4]
scale = sys.argv[5]

mat = None
trueZ = None

if targetZ == "26":
    trueZ = "Iron"
    mat = "Fe"

if targetZ == "82":
    trueZ = "Lead"
    mat = "Pb"

if targetZ == "06":
    trueZ = "Carbon"
    mat = "C"

if targetZ == "99":
    trueZ = "Tracker"
    mat = "CH"

if targetZ == "99":
    infile= ROOT.TFile(str(dirpwd)+"/CrossSection_t%s_z%02s_%s.root"%(targetID, targetZ, plist))
else:
    infile= ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t%s_z%02s_%s.root"%(targetID, targetZ, plist))
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale = dataPOT/mcPOT
if scale == "1":
    mcScale = 1

vars = ["Enu", "x", "pZmu1D", "pTmu"]

#steps = ['unfolded','unfolded_effCorrected', 'crossSection', 'crossSection_total']
steps = ['total_unfolded_effCorrected','crossSection', 'crossSection_total']

for step in steps:

    for var in vars:

        data_hist = infile.Get("%s_data_%s"%(step, var))

        if step == "unfolded":
            mc_hist = infile.Get("total_efficiency_numerator_%s"%(var))
        elif step == "total_unfolded_effCorrected":
            mc_hist = infile.Get("total_simEventRate_%s"%(var))
            mc_hist_QE = infile.Get("total_simEventRate_QE_%s"%(var))
            mc_hist_RES = infile.Get("total_simEventRate_RES_%s"%(var))
            mc_hist_DIS = infile.Get("total_simEventRate_DIS_%s"%(var))
            mc_hist_Other = infile.Get("total_simEventRate_Other_%s"%(var))
            mc_hist_2p2h = infile.Get("total_simEventRate_2p2h_%s"%(var))
        else:
            mc_hist = infile.Get("simEventRate_%s_mc_%s"%(step, var))
            mc_hist_QE = infile.Get("simEventRate_QE_%s_mc_%s"%(step, var))
            mc_hist_RES = infile.Get("simEventRate_RES_%s_mc_%s"%(step, var))
            mc_hist_DIS = infile.Get("simEventRate_DIS_%s_mc_%s"%(step, var))
            mc_hist_Other = infile.Get("simEventRate_Other_%s_mc_%s"%(step, var))
            mc_hist_2p2h = infile.Get("simEventRate_2p2h_%s_mc_%s"%(step, var))
    

        if var == "Enu":
            mc_hist.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dE_{#bar{#nu}} (10^{-39} cm^{2}/GeV/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3}/GeV")

        if var == "x":
            mc_hist.GetXaxis().SetTitle("Bjorken x")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")

        if var == "pTmu":
            mc_hist.GetXaxis().SetTitle("Muon p_{T} (GeV/c)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")
       
        if var == "pZmu1D":
            mc_hist.GetXaxis().SetTitle("Muon p_{Z} (GeV/c)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")

        if var == "ThetamuDeg":
            mc_hist.GetXaxis().SetTitle("Muon Angle (Deg)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")

        if step == "total_unfolded_effCorrected":
            mc_hist.Scale(1/1E3)
            mc_hist_QE.Scale(1/1E3)
            mc_hist_RES.Scale(1/1E3)
            mc_hist_DIS.Scale(1/1E3)
            mc_hist_Other.Scale(1/1E3)
            mc_hist_2p2h.Scale(1/1E3)
            data_hist.Scale(1/1E3)

        if step == "crossSection":
            mc_hist.Scale(1E39)
            mc_hist_QE.Scale(1E39)
            mc_hist_RES.Scale(1E39)
            mc_hist_DIS.Scale(1E39)
            mc_hist_Other.Scale(1E39)
            mc_hist_2p2h.Scale(1E39)
            data_hist.Scale(1E39)
        
        if step == "crossSection_total":
            mc_hist.Scale(1E38)
            mc_hist_QE.Scale(1E38)
            mc_hist_RES.Scale(1E38)
            mc_hist_DIS.Scale(1E38)
            mc_hist_Other.Scale(1E38)
            mc_hist_2p2h.Scale(1E38)
            data_hist.Scale(1E38)

        mc_hist.GetXaxis().CenterTitle()
        mc_hist.GetYaxis().CenterTitle()
        #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
        #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
        #mc_hist.GetXaxis().SetRangeUser(66,173)

        #mnv.ApplyStyle(1)
        if step == "unfolded":
            mc_hist.Scale(mcScale)
        elif step == "total_unfolded_effCorrected":
            mc_hist.Scale(mcScale)
            mc_hist_QE.Scale(mcScale)
            mc_hist_RES.Scale(mcScale)
            mc_hist_DIS.Scale(mcScale)
            mc_hist_Other.Scale(mcScale)
            mc_hist_2p2h.Scale(mcScale)
        else:
            mc_hist.Scale(1)
            mc_hist_QE.Scale(1)
            mc_hist_RES.Scale(1)
            mc_hist_DIS.Scale(1)
            mc_hist_Other.Scale(1)
            mc_hist_2p2h.Scale(1)

        #if var == "Enu":
        #    mc_hist.GetXaxis().SetRangeUser(2, 20)
        
        mc_hist.Scale(mc_hist.GetNormBinWidth(), "width")
        mc_hist_QE.Scale(mc_hist_QE.GetNormBinWidth(), "width")
        mc_hist_RES.Scale(mc_hist_RES.GetNormBinWidth(), "width")
        mc_hist_DIS.Scale(mc_hist_DIS.GetNormBinWidth(), "width")
        mc_hist_Other.Scale(mc_hist_Other.GetNormBinWidth(), "width")
        mc_hist_2p2h.Scale(mc_hist_2p2h.GetNormBinWidth(), "width")
        data_hist.Scale(mc_hist.GetNormBinWidth(), "width")
        data_hist_stat =  data_hist.GetCVHistoWithStatError() # stat error
        data_hist_total = data_hist.GetCVHistoWithError() # total error
        data_hist_sys = data_hist.GetCVHistoWithError(False) # sys error (bool is include stat)
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
        if step == "crossSection_total":
            mc_hist.GetYaxis().SetRangeUser(0, data_hist.GetMaximum()*1.2)
        elif var == "x":
            mc_hist.SetMaximum(data_hist.GetMaximum()*1.5)
        else:
            mc_hist.SetMaximum(data_hist.GetMaximum()*1.3)

        # Int channels
        mc_hist_QE.SetLineWidth(3)
        mc_hist_QE.SetLineColor(38)
        mc_hist_RES.SetLineColor(30)
        mc_hist_DIS.SetLineColor(40)
        mc_hist_Other.SetLineColor(14)
        mc_hist_2p2h.SetLineColor(41)

        #gStyle.SetErrorX()
        mc_hist.Draw("HIST")
        #gStyle.SetErrorX(0.5)
        mc_hist_stat.Draw("E2 SAME")
        #gStyle.SetErrorX(0)
        mc_hist.Draw("HIST SAME")
        mc_hist_QE.Draw("HIST SAME")
        mc_hist_RES.Draw("HIST SAME")
        mc_hist_DIS.Draw("HIST SAME")
        mc_hist_Other.Draw("HIST SAME")
        mc_hist_2p2h.Draw("HIST SAME")
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

        #chi2 = mnv.Chi2DataMC(data_hist, mc_hist,mcScale,False,True )
        #ndf = mc_hist.GetNbinsX()
        #print(chi2)
        if step == "unfolded":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        elif step == "total_unfolded_effCorrected":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        elif step == "crossSection_total":
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.82)
        else:
            mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)

        '''
        Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
        If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
        If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
        '''
        if step == 'crossSection':
            mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
        
        elif step == 'crossSection_total':
            mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TR", 0.035, 0.0, True, False)
        
        if targetZ == "99":
            mnv.AddHistoTitle("%s"%(trueZ), 0.05, 1)
        else:
            mnv.AddHistoTitle("Target %s %s: %s"%(targetID, trueZ, step), 0.04, 1)

        mc_hist.GetYaxis().SetTitleOffset(0.96)

        if step == "crossSection_total":
            legend = TLegend(0.20,0.45,0.50,0.89)
        else:
            legend = TLegend(0.55,0.55,0.80,0.89)
        if var == "x":
            legend = TLegend(0.60,0.55,0.8,0.89)
            if step == "total_unfolded_effCorrected":
                legend = TLegend(0.45,0.7,0.8,0.89)
                legend.SetNColumns(2)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.035)
        legend.AddEntry(mc_hist_stat, " Simulation", "fl")
        legend.AddEntry(mc_hist_QE, " QE", "fl")
        legend.AddEntry(mc_hist_RES, " RES", "fl")
        legend.AddEntry(mc_hist_2p2h, " 2p2h", "fl")
        legend.AddEntry(mc_hist_DIS, " DIS", "fl")
        legend.AddEntry(mc_hist_Other, " Other", "fl")
        legend.AddEntry(data_hist, " Data", "lep")
        #legend.SetTextFont(42)
        legend.Draw()

        canvas1.SetLogx(False)
        if var == "x":
            canvas1.SetLogx()
        canvas1.Modified()
        if targetZ == "99":
            canvas1.Print("%s_t%s_z%02s_%s_%s.png"%(step, targetID, targetZ, var, plist))
        else:
            canvas1.Print("%s_Daisy_t%s_z%02s_%s_%s.png"%(step, targetID, targetZ, var, plist))

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

        mnv2.error_summary_group_map["Hadrons"].push_back("response_em")
        mnv2.error_summary_group_map["Hadrons"].push_back("response_low_proton")
        mnv2.error_summary_group_map["Hadrons"].push_back("response_mid_proton")
        mnv2.error_summary_group_map["Hadrons"].push_back("response_high_proton")
        mnv2.error_summary_group_map["Hadrons"].push_back("response_meson")
        mnv2.error_summary_group_map["Hadrons"].push_back("response_other")

        mnv2.error_summary_group_map["Hadrons"].push_back("GEANT_Neutron")
        mnv2.error_summary_group_map["Hadrons"].push_back("GEANT_Proton")
        mnv2.error_summary_group_map["Hadrons"].push_back("GEANT_Pion")

        mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_C")
        mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_Fe")
        mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_Pb")
        mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_CH")
        mnv2.error_summary_group_map["Target Mass"].push_back("Target_Mass_H2O")
    
        # Change colour
        mnv2.error_color_map["Target Mass"] = 15
        mnv2.error_color_map["Muon Efficiency"] = 46
        mnv2.error_color_map["Particle Response"] = 65
        
        # ---------------------------------------------
        # Data fractional error
        if var == "Enu":
            data_hist.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
            #data_hist.GetXaxis().SetRangeUser(2, 20)

        if var == "x":
            data_hist.GetXaxis().SetTitle("Bjorken x")

        if var == "pZmu1D":
            data_hist.GetXaxis().SetTitle("Muon p_{Z} (GeV/c)")
            #data_hist.GetXaxis().SetRangeUser(2, 20)

        if var == "pTmu":
            data_hist.GetXaxis().SetTitle("Muon p_{T} (GeV/c)")

        if var == "pTmu":
            data_hist.GetXaxis().SetTitle("Muon Angle (Deg)")

        mnv2.DrawErrorSummary(data_hist, "TL", True, True, 0.0, False, "",True);
        # last boolean decides whether frac or not

        keys = canvas1.GetListOfPrimitives();
        for k in keys:
            if(k.ClassName().find("TH1")!=-1):
                if var == "Enu":
                    k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.03) # to 0.35 for Enu
                if var == "x":
                    k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.05)
                if var == "pTmu":
                    k.GetYaxis().SetRangeUser(0,k.GetMaximum())
                if var == "pZmu1D":
                    k.GetYaxis().SetRangeUser(0,k.GetMaximum())
                if var == "ThetamuDeg":
                    k.GetYaxis().SetRangeUser(0,k.GetMaximum())

            if(k.ClassName().find("Legend")!=-1):
                if var == "Enu":
                    k.SetNColumns(2)
                    k.SetX2(20) #Enu
                    k.SetY1(0.15) #Enu
                    if len(sys.argv) > 4:
                        k.SetY1(0.1) #Enu
                    if step == "total_unfolded_effCorrected":
                        if targetZ == "06":
                            k.SetY1(0.25) #Enu
                            if len(sys.argv) > 4:
                                k.SetY1(0.08) #Enu
                        if targetZ == "82":
                            if len(sys.argv) > 4:
                                k.SetY1(0.06) #Enu
                        if targetZ == "99":
                            k.SetY1(0.04) #Enu
                    if step == "crossSection_total":
                        k.SetY1(0.35) #Enu
                        if targetZ == "06":
                            k.SetY1(0.1) #Enu
                        if targetZ == "26":
                            k.SetY1(0.12) #Enu
                        if targetZ == "82":
                            k.SetY1(0.1) #Enu
                        if targetZ == "99":
                            k.SetY1(0.08) #Enu
                    if step == "crossSection":
                        if targetZ == "06":
                            k.SetY1(0.25) #Enu
                            if len(sys.argv) > 4:
                                k.SetY1(0.1) #Enu
                        if targetZ == "99":
                            k.SetY1(0.07) #Enu
                        else:
                            k.SetY1(0.15) #Enu
                            if len(sys.argv) > 4:
                                k.SetY1(0.1) #Enu
                                if targetZ == "82":
                                    k.SetY1(0.08) #Enu

                if var == "x":
                    k.SetNColumns(2)
                    k.SetX2(0.5) #x
                    k.SetY1(0.12) #x
                    if step == "total_unfolded_effCorrected":
                        if targetZ == "06":
                            k.SetY1(0.052) #x
                        if targetZ == "26":
                            k.SetY1(0.11) #x
                            if len(sys.argv) > 4:
                                k.SetY1(0.09) #x
                        if targetZ == "82":
                            k.SetY1(0.05) #x
                        if targetZ == "99":
                            k.SetY1(0.027) #x
                    if step == "crossSection":
                        if targetZ == "06" or "82":
                            k.SetY1(0.07) #x
                        if targetZ == "26":
                            k.SetY1(0.105) #x
                        if targetZ == "99":
                            k.SetY1(0.065) #x

                if var == "pTmu":
                    k.SetNColumns(2)
                    k.SetX2(4.5) #x
                    if targetZ == "06":
                        k.SetY1(0.45) #x
                    if targetZ == "26":
                        k.SetY1(0.3) #x
                    if targetZ == "82":
                        k.SetY1(0.25) #x
                    if targetZ == "99":
                        k.SetY1(0.065) #x
                if var == "pZmu1D":
                    k.SetNColumns(2)
                    k.SetX2(20) #x
                if var == "ThetamuDeg":
                    k.SetNColumns(2)
                    k.SetX2(17) #x
                    if step == "total_unfolded_effCorrected":
                        if targetZ == "06":
                            k.SetY1(0.08) #x
                        if targetZ == "26":
                            k.SetY1(0.065) #x
                        if targetZ == "82":
                            k.SetY1(0.06) #x
                        if targetZ == "99":
                            k.SetY1(0.017) #x
                    if step == "crossSection":
                        if targetZ == "06":
                            k.SetY1(0.1) #x
                        if targetZ == "26":
                            k.SetY1(0.085) #x
                        if targetZ == "82":
                            k.SetY1(0.085) #x
                        if targetZ == "99":
                            k.SetY1(0.06) #x

        if targetZ == "99":
            mnv2.AddHistoTitle("%s"%(trueZ), 0.05, 1)
        else:
            mnv2.AddHistoTitle("Target %s %s: %s"%(targetID, trueZ, step), 0.04, 1)
        if var == "Enu":
            #mnv2.WritePreliminary(0.37, 0.5, 0.035, True)
            if targetZ == "82":
                mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.61, 0.033, 12, 42)
            else:
                mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.57, 0.033, 12, 42)

        if var == "x":
            #mnv2.WritePreliminary(0.37, 0.53, 0.035, True)
            if step == "total_unfolded_effCorrected":
                if targetZ == "06":
                    mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.65, 0.033, 12, 42)
                if targetZ == "26":
                    mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.65, 0.033, 12, 42)
                if targetZ == "82":
                    mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.65, 0.033, 12, 42)
                if targetZ == "99":
                    mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.65, 0.033, 12, 42)

            else:
                if targetZ == "06":
                    mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.62, 0.033, 12, 42)
                if targetZ == "26":
                    mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.65, 0.033, 12, 42)
                if targetZ == "82":
                    mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.63, 0.033, 12, 42)
                if targetZ == "99":
                    mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.63, 0.033, 12, 42)
        
        if var == "pTmu":
            #mnv2.WritePreliminary(0.37, 0.53, 0.035, True)
            if step == "total_unfolded_effCorrected":
                mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.62, 0.033, 12, 42)
            else:
                mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.30, 0.65, 0.033, 12, 42)
        
        if var == "pZmu1D":
            #mnv2.WritePreliminary(0.37, 0.53, 0.035, True)
            if step == "total_unfolded_effCorrected":
                mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.35, 0.033, 12, 42)
            else:
                mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.70, 0.35, 0.033, 12, 42)
        
        if var == "ThetamuDeg":
            #mnv2.WritePreliminary(0.37, 0.53, 0.035, True)
            if step == "total_unfolded_effCorrected":
                mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.40, 0.62, 0.033, 12, 42)
            else:
                mnv2.AddPlotLabel("Data POT "+ "{:.2e}".format(dataPOT), 0.35, 0.65, 0.033, 12, 42)


        canvas1.Modified()
        if targetZ == "99":
            canvas1.Print("%s_FracErr_t%s_z%02s_%s_%s.png"%(step, targetID, targetZ, var, plist))  
        else:
            canvas1.Print("%s_Daisy_FracErr_t%s_z%02s_%s_%s.png"%(step, targetID, targetZ, var, plist))        

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
#canvas2.Print("ME6A_FluxRebinned.png")

        

data_hist.SetDirectory(0)
mc_hist.SetDirectory(0)

print("DONE %s %s %02s"%(plist, targetID, targetZ))
