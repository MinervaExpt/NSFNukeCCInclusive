import ROOT
import os,sys
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TGaxis
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TColor, TLegend
from ROOT import PlotUtils

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
targetZ = sys.argv[2]
plist = sys.argv[3]
intType = sys.argv[4]
#material = 'iron t235'

if targetZ == "26":
  material = "iron"
  mat = "Fe"
  targetID = "235"

if targetZ == "82":
  material = "lead"
  mat = "Pb"
  targetID = "2345"

if targetZ == "06":
  material = "carbon"
  mat = "C"
  targetID = "3"


target = ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t%s_z%02s_%s.root"%(targetID, targetZ, plist))
tracker = ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t99_z99_%s.root"%(plist))

print("CrossSection_Daisy_t%s_z%02s_%s.root"%(targetID, targetZ, plist))

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = target.Get("MCPOT").GetVal()
dataPOT = tracker.Get("DataPOT").GetVal()

ROOT.TH1.AddDirectory(False)

vars = ["Enu","x", "pZmu1D", "pTmu1D"]

for var in vars:
    if var=="Enu":
        target_data_hist = target.Get("crossSection_total_data_%s"%var)
        target_simEventRate_hist = target.Get("simEventRate_crossSection_total_mc_%s"%var) # total signal
        target_simEventRate_hist_QE = target.Get("simEventRate_QE_crossSection_total_mc_%s"%var) # QE
        target_simEventRate_hist_RES = target.Get("simEventRate_RES_crossSection_total_mc_%s"%var)
        target_simEventRate_hist_DIS = target.Get("simEventRate_DIS_crossSection_total_mc_%s"%var) 
        target_simEventRate_hist_Other = target.Get("simEventRate_Other_crossSection_total_mc_%s"%var)
        target_simEventRate_hist_2p2h = target.Get("simEventRate_2p2h_crossSection_total_mc_%s"%var) 
        tracker_data_hist = tracker.Get("crossSection_total_%s_data_%s"%(material,var))
        tracker_simEventRate_hist = tracker.Get("simEventRate_crossSection_total_%s_mc_%s"%(material,var))
        tracker_simEventRate_hist_QE = tracker.Get("simEventRate_QE_crossSection_total_%s_mc_%s"%(material,var))
        tracker_simEventRate_hist_RES = tracker.Get("simEventRate_RES_crossSection_total_%s_mc_%s"%(material,var))
        tracker_simEventRate_hist_DIS = tracker.Get("simEventRate_DIS_crossSection_total_%s_mc_%s"%(material,var))
        tracker_simEventRate_hist_Other = tracker.Get("simEventRate_Other_crossSection_total_%s_mc_%s"%(material,var))
        tracker_simEventRate_hist_2p2h = tracker.Get("simEventRate_2p2h_crossSection_total_%s_mc_%s"%(material,var))

    else: 
        target_data_hist = target.Get("crossSection_data_%s"%var)
        target_simEventRate_hist = target.Get("simEventRate_crossSection_mc_%s"%var) # total signal
        target_simEventRate_hist_QE = target.Get("simEventRate_QE_crossSection_mc_%s"%var)
        target_simEventRate_hist_RES = target.Get("simEventRate_RES_crossSection_mc_%s"%var)
        target_simEventRate_hist_DIS = target.Get("simEventRate_DIS_crossSection_mc_%s"%var)
        target_simEventRate_hist_Other = target.Get("simEventRate_Other_crossSection_mc_%s"%var)
        target_simEventRate_hist_2p2h = target.Get("simEventRate_2p2h_crossSection_mc_%s"%var)
        tracker_data_hist = tracker.Get("crossSection_%s_data_%s"%(material,var))
        tracker_simEventRate_hist = tracker.Get("simEventRate_crossSection_%s_mc_%s"%(material,var))
        tracker_simEventRate_hist_QE = tracker.Get("simEventRate_QE_crossSection_%s_mc_%s"%(material,var))
        tracker_simEventRate_hist_RES = tracker.Get("simEventRate_RES_crossSection_%s_mc_%s"%(material,var))
        tracker_simEventRate_hist_DIS = tracker.Get("simEventRate_DIS_crossSection_%s_mc_%s"%(material,var))
        tracker_simEventRate_hist_Other = tracker.Get("simEventRate_Other_crossSection_%s_mc_%s"%(material,var))
        tracker_simEventRate_hist_2p2h = tracker.Get("simEventRate_2p2h_crossSection_%s_mc_%s"%(material,var))


    ratio_data = target_data_hist.Clone()
    ratio_data.Divide(target_data_hist,tracker_data_hist) 
    ratio_mc = target_simEventRate_hist.Clone()
    ratio_mc.Divide(target_simEventRate_hist,tracker_simEventRate_hist)

    # int types
    ratio_mc_QE = target_simEventRate_hist_QE.Clone()
    ratio_mc_QE.Divide(target_simEventRate_hist_QE,tracker_simEventRate_hist_QE)

    ratio_mc_RES = target_simEventRate_hist_RES.Clone()
    ratio_mc_RES.Divide(target_simEventRate_hist_RES,tracker_simEventRate_hist_RES)

    ratio_mc_DIS = target_simEventRate_hist_DIS.Clone()
    ratio_mc_DIS.Divide(target_simEventRate_hist_DIS,tracker_simEventRate_hist_DIS)

    ratio_mc_Other = target_simEventRate_hist_Other.Clone()
    ratio_mc_Other.Divide(target_simEventRate_hist_Other,tracker_simEventRate_hist_Other)

    ratio_mc_2p2h = target_simEventRate_hist_2p2h.Clone()
    ratio_mc_2p2h.Divide(target_simEventRate_hist_2p2h,tracker_simEventRate_hist_2p2h)

    if var == "Enu":
        ratio_mc.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
        ratio_mc.GetYaxis().SetTitle("Events/GeV")
        #ratio_mc.GetXaxis().SetRangeUser(2, 20)

    if var == "x":
        ratio_mc.GetXaxis().SetTitle("Bjorken x")
        ratio_mc.GetYaxis().SetTitle("Events (norm.)")

    if var == "pTmu1D":
        ratio_mc.GetXaxis().SetTitle("Muon p_{T} (GeV/c)")
        ratio_mc.GetYaxis().SetTitle("Events/(GeV/c)")
       
    if var == "pZmu1D":
        ratio_mc.GetXaxis().SetTitle("Muon p_{Z} (GeV/c)")
        ratio_mc.GetYaxis().SetTitle("Events/(GeV/c)")

    if var == "ThetamuDeg":
        ratio_mc.GetXaxis().SetTitle("Muon Angle (Deg)")
        ratio_mc.GetYaxis().SetTitle("Events/Deg")

    ratio_mc.GetXaxis().CenterTitle()
    ratio_mc.GetYaxis().CenterTitle()
    ratio_mc.GetYaxis().SetNdivisions(505)
    ratio_mc.GetXaxis().SetTitleSize(0.05)
    ratio_mc.GetYaxis().SetTitleSize(0.05)
    ratio_mc.GetYaxis().SetTitleOffset(1.6)

    # ----------------------------------------------------------------------------
    # Ratio

    if var == "Enu":
        ratio_mc.GetYaxis().SetTitle("#sigma^{%s}/#sigma^{CH}"%mat)
        ratio_mc.GetYaxis().SetTitleOffset(1.1)

    else:
        ratio_mc.GetYaxis().SetTitle("#frac{d#sigma^{%s}}{dx_{#bar{#nu}}} / #frac{d#sigma^{CH}}{dx_{#bar{#nu}}} "%mat)
        ratio_mc.GetYaxis().SetTitleOffset(1.2)
    ratio_mc.GetXaxis().CenterTitle()
    ratio_mc.GetYaxis().CenterTitle()
    ratio_mc.Draw("HIST")

    ratio_data_stat = ratio_data.GetCVHistoWithStatError() # stat error
    ratio_data_total = ratio_data.GetCVHistoWithError() # total error

    ratio_mc_stat = ratio_mc.GetCVHistoWithStatError() 

    # MC
    ratio_mc.SetLineWidth(3)
    ratio_mc.SetLineColor(2)
    ratio_mc.SetLineColor(ROOT.kRed)
    ratio_mc.SetLineWidth(2)
    # MC error
    ratio_mc_stat.SetFillColor(ROOT.kRed-10)
    ratio_mc_stat.SetFillStyle(1001)
    ratio_mc_stat.SetMarkerStyle(0)
    ratio_mc_stat.SetLineWidth(3)
    ratio_mc_stat.SetLineColor(2)

    # int types
    ratio_mc_QE.SetLineWidth(3)
    ratio_mc_QE.SetLineColor(38)
    ratio_mc_RES.SetLineWidth(3)
    ratio_mc_RES.SetLineColor(30)
    ratio_mc_DIS.SetLineColor(9)
    ratio_mc_Other.SetLineColor(14)
    ratio_mc_2p2h.SetLineColor(41)

    ratio_mc_stat.Draw("E2 SAME")
    ratio_mc.Draw("HIST SAME")
    ratio_data_stat.Draw("SAME E1 X0")
    ratio_data_total.Draw("E1 SAME X0")
    if intType == "0":
        ratio_mc_QE.Draw("HIST SAME")
        ratio_mc_RES.Draw("HIST SAME")
        ratio_mc_DIS.Draw("HIST SAME")
        ratio_mc_Other.Draw("HIST SAME")
        ratio_mc_2p2h.Draw("HIST SAME")

    #ratio_mc.SetMaximum(ratio_data.GetMaximum()*1.2)
    #ratio_mc.SetMinimum(ratio_data.GetMinimum()*0.8)

    ratio_mc.SetMaximum(1.5)
    ratio_mc.SetMinimum(0.5)

    gStyle.SetOptTitle(0)
    mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)

    mnv.AddHistoTitle("xSec %s/tracker"%material, 0.05, 1)
    if intType == "0":
        legend = TLegend(0.85,0.55,0.99,0.89)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.035)
        legend.AddEntry(ratio_mc_QE, " QE", "l")
        legend.AddEntry(ratio_mc_RES, " RES", "l")
        legend.AddEntry(ratio_mc_2p2h, " 2p2h", "l")
        legend.AddEntry(ratio_mc_DIS, " DIS", "l")
        legend.AddEntry(ratio_mc_Other, " Other", "l")
        legend.SetTextFont(42)
        legend.Draw()

    canvas1.SetLogx(False)
    if var == "x":
        canvas1.SetLogx()
    canvas1.Modified()
    canvas1.Update()
    if intType == "0":
        canvas1.Print("xSec_ratio_Daisy_t%s_z%02s_%s_%s_intType.png"%(targetID, targetZ, var, plist))
    else:
        canvas1.Print("xSec_ratio_Daisy_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))

 # ----------------------------------------------------------------------------
    # Fractional Error Groups
    mnv2 = PlotUtils.MnvPlotter()
    mnv2.ApplyStyle(7)

    if var == "Enu":
        ratio_data.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
        ratio_data.GetYaxis().SetTitle("Events/GeV")
        #ratio_data.GetXaxis().SetRangeUser(2, 20)


    if var == "x":
        ratio_data.GetXaxis().SetTitle("Bjorken x")
        ratio_data.GetYaxis().SetTitle("Events (norm.)")

    if var == "pTmu1D":
        ratio_data.GetXaxis().SetTitle("Reconstructed Muon p_{T} (GeV/c)")
        ratio_data.GetYaxis().SetTitle("Events/(GeV/c)")
       
    if var == "pZmu1D":
        ratio_data.GetXaxis().SetTitle("Reconstructed Muon p_{Z} (GeV/c)")
        ratio_data.GetYaxis().SetTitle("Events/(GeV/c)")

    if var == "ThetamuDeg":
        ratio_data.GetXaxis().SetTitle("Reconstructed Muon Angle (Deg)")
        ratio_data.GetYaxis().SetTitle("Events/Deg")

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
    

    #mnv2.ApplyStyle(7)
    mnv2.DrawErrorSummary(ratio_data, "TL", True, True, 0.0, False, "",True);
    # last boolean decides whether frac or not

    keys = canvas1.GetListOfPrimitives();
    for k in keys:
        if(k.ClassName().find("TH1")!=-1):
            if var == "Enu":
                #k.GetYaxis().SetRangeUser(0, 0.35)
                k.GetMaximum()*1.0000001 # to 0.35 for Enu
            if var == "x":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.005)
            if var == "pTmu1D":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum())
            if var == "pZmu1D":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum())
            if var == "ThetamuDeg":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum())

        if(k.ClassName().find("Legend")!=-1):
            if var == "Enu":
                k.SetNColumns(2)
                k.SetX2(20) #Enu
                k.SetY1(0.1) #Enu
                if targetZ =="06":
                    k.SetY1(0.22) #Enu
                    if len(sys.argv) > 4:
                        k.SetY1(0.08) #Enu
                else:
                    k.SetY1(0.15) #Enu
                    if len(sys.argv) > 4:
                        if targetZ =="26":
                            k.SetY1(0.09) #Enu
                        if targetZ =="82":
                            k.SetY1(0.06) #Enu

            if var == "x":
                k.SetNColumns(2)
                k.SetX2(0.5) #x
                k.SetY1(0.035) #x
                if targetZ =="06":
                    k.SetY1(0.055) #x
                    if len(sys.argv) > 4:
                        k.SetY1(0.035) #x
                if targetZ =="26":
                    k.SetY1(0.09) #x
                    if len(sys.argv) > 4:
                        k.SetY1(0.075) #x
                if targetZ =="82":
                    k.SetY1(0.05)
                    if len(sys.argv) > 4:
                        k.SetY1(0.045) #x

            if var == "pTmu1D":
                k.SetNColumns(2)
                k.SetX2(4.5) #Enu
                if targetZ =="06":
                    k.SetY1(0.45) #x
                if targetZ =="26":
                    k.SetY1(0.3) #x
                if targetZ =="82":
                    k.SetY1(0.25) #x
            if var == "pZmu1D":
                k.SetNColumns(2)
                k.SetX2(20) #Enu
            if var == "ThetamuDeg":
                k.SetNColumns(2)
                k.SetX2(17) #Enu
                if targetZ =="06":
                    k.SetY1(0.08) #x
                if targetZ =="26":
                    k.SetY1(0.06) #x
                if targetZ =="82":
                    k.SetY1(0.062) #x


    mnv2.AddHistoTitle("xSec %s/tracker"%material, 0.04, 1)

    canvas1.Modified()
    canvas1.Print("xSec_ratio_Daisy_FracErr_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))


print("DONE %s %s %02s"%(plist, targetID, targetZ))
