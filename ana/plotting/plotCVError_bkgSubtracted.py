import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle

dirpwd = sys.argv[1]
targetID = sys.argv[2] 
targetZ = sys.argv[3]
plist = sys.argv[4]

infile= ROOT.TFile(str(dirpwd)+"/BkgSubtracted_EventSelection_%s_t%s_z%02s_sys.root"%(plist, targetID, targetZ))

canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale = None
if len(sys.argv) > 4:
    mcScale = 1
else:
    mcScale =  dataPOT/mcPOT

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

vars = ["Enu", "x"]

for var in vars:

    mc_hist = infile.Get("h_bkg_subtracted_mc_%s"%var)
    data_hist = infile.Get("h_bkg_subtracted_data_%s"%var)
    data_hist_notConstrained = None
    if targetZ != "99":
        data_hist_notConstrained = infile.Get("h_bkg_subtracted_data_notConstrained_%s"%var)

    if var == "Enu":
        mc_hist.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")
        mc_hist.GetYaxis().SetTitle("Events/GeV")
        data_hist.GetXaxis().SetTitle("Reconstructed Antineutrino Energy (GeV)")
        if targetZ != "99":
            data_hist_notConstrained.GetXaxis().SetTitle("Reconstructed Neutrino Energy (GeV)")

    if var == "x":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
        data_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        if targetZ != "99":
            data_hist_notConstrained.GetXaxis().SetTitle("Reconstructed Bjorken x") 

    '''
    if myvariable == "xfine":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")


    if myvariable == "xBrian":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")

    if myvariable == "x09":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x")
        mc_hist.GetYaxis().SetTitle("Events (norm.)")
    ''' 
    mc_hist.GetXaxis().CenterTitle()
    mc_hist.GetYaxis().CenterTitle()
    #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
    #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
    #mc_hist.GetXaxis().SetRangeUser(0,0.001)

    #mnv.ApplyStyle(1)
    mnv.DrawDataMCWithErrorBand(data_hist, mc_hist, mcScale, "TR", False, None, None, False, True, False)

    #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)

    '''
    Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
    If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
    If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
    '''
    if targetZ == "99":
        mnv.AddHistoTitle("%s: Bkg Subtracted"%(trueZ), 0.05, 1)
    else:
        mnv.AddHistoTitle("Target %s %s: Bkg Subtracted"%(targetID, trueZ), 0.05, 1)
    mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)

    canvas1.Modified()
    canvas1.Print("BkgSubtracted_EventSelection_t%s_z%02s_%s_%s.png"%(targetID, targetZ, var, plist))

    # ----------------------------------------------------------------------------
    # Background subtracted with data not constrained
    '''
    mnv.DrawDataMCWithErrorBand(data_hist, mc_hist, mcScale, "TR", False, None, None, False, True, False)
    mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
    mnv.AddHistoTitle("Target %s %s: Background subtracted (not constrained)"%(targetID, trueZ), 0.035, 1)
    mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
    canvas1.Print("ME6A_T%s%s_BkgSubtracted_NotConstrained_EventSelection_%s_sys.png"%(targetID,mat,var))
    '''
    # ----------------------------------------------------------------------------
    # Fractional Error Groups: MC
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
    
    #mnv2.ApplyStyle(7)
    mnv2.DrawErrorSummary(mc_hist, "TL", True, True, 0.0, False, "",True);
    # last boolean decides whether frac or not

    keys = canvas1.GetListOfPrimitives();
    for k in keys:
        if(k.ClassName().find("TH1")!=-1):
            if var == "Enu":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()) # to 0.35 for Enu
            if var == "x":
                k.GetYaxis().SetRangeUser(0,k.GetMaximum()*1.05)

        if(k.ClassName().find("Legend")!=-1):
            if var == "Enu":
                k.SetNColumns(2)
                k.SetX2(45) #Enu
                k.SetY1(0.25) #Enu
                if targetZ == "99":
                    k.SetX2(45) #Enu
                    k.SetY1(0.2) #Enu
            if var == "x":
                k.SetNColumns(2)
                k.SetX2(2.0) #x
                k.SetY1(0.16) #x
                if targetZ == "99":
                    k.SetX2(2.0) #x
                    k.SetY1(0.17) #x

    if targetZ == "99":
        mnv2.AddHistoTitle("%s: Bkg Subtracted (MC)"%(trueZ), 0.05, 1)
    else:
        mnv2.AddHistoTitle("Target %s %s: Bkg Subtracted (MC)"%(targetID, trueZ), 0.04, 1)

    canvas1.Modified()
    canvas1.Print("BkgSubtracted_EventSelection_t%s_z%02s_%s_%s_FracErr_MC.png"%(targetID, targetZ, var, plist))

    # ----------------------------------------------------------------------------
    # Fractional Error Groups: Data

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
                k.SetY1(0.25) #Enu
                if len(sys.argv) > 4:
                    k.SetY1(0.15) #Enu
                if targetZ == "99":
                    k.SetX2(45) #Enu
                    k.SetY1(0.05) #Enu
                    if len(sys.argv) > 4:
                        k.SetY1(0.03) #Enu

            if var == "x":
                k.SetNColumns(2)
                k.SetX2(2.0) #x
                k.SetY1(0.040) #x
                if len(sys.argv) > 4:
                    k.SetY1(0.03) #x
                if targetZ == "99":
                    k.SetX2(2.0) #x
                    k.SetY1(0.007) #x
                    if len(sys.argv) > 4:
                        k.SetY1(0.005) #x



    if targetZ == "99":
        mnv2.AddHistoTitle("%s: Bkg Subtracted (Data)"%(trueZ), 0.05, 1)
    else:
        mnv2.AddHistoTitle("Target %s %s: Bkg Subtracted (Data)"%(targetID, trueZ), 0.04, 1)

    canvas1.Modified()
    canvas1.Print("BkgSubtracted_EventSelection_t%s_z%02s_%s_%s_FracErr_Data.png"%(targetID, targetZ, var, plist))

    # ----------------------------------------------------------------------------
    # Fractional Error Groups: Data (NOT cosntrained)

    # Data fractional error
    if targetZ != "99":
        mnv2.DrawErrorSummary(data_hist_notConstrained, "TL", True, True, 0.0, False, "",True);
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
                    if len(sys.argv) > 4:
                        k.SetY1(0.15) #Enu

                if var == "x":
                    k.SetNColumns(2)
                    k.SetX2(2.0) #x
                    k.SetY1(0.045) #x
                    if len(sys.argv) > 4:
                        k.SetY1(0.03) #x

        mnv2.AddHistoTitle("Target %s %s: Bkg Subtracted (Data, not constr.)"%(targetID, trueZ), 0.03, 1)

        canvas1.Modified()
        canvas1.Print("BkgSubtracted_EventSelection_t%s_z%02s_%s_%s_FracErr_DataNotConstrained.png"%(targetID, targetZ, var, plist))

data_hist.SetDirectory(0)
mc_hist.SetDirectory(0)

print("DONE %s %s %02s"%(plist, targetID, targetZ))