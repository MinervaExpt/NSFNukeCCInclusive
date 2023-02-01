import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import TMath
import math as m
_fnsf = ROOT.TFile("Hists_PlasticBkg_sys_t3_z26_AntiNu.root","READ")
_fnsf2 = ROOT.TFile("Plastic_ScaleFactors_t3_z26_minervame6A.root")

myvariable = sys.argv[1]

mc_name ="DS_regDS_Iron_Enu"
scale = _fnsf2.Get("scaleFactor_DS_Enu")
mc_new = _fnsf.Get(mc_name) # just cv with sys
integral = mc_new.Integral()
statErr = 1/ m.sqrt(integral)
for bin in range(scale.GetNbinsX() + 1):
	scale.SetBinError(bin,0+statErr)	
	print(scale.GetBinError(bin))

mc_new.Multiply(mc_new, scale)

if myvariable == "Enu":
	mc_new.GetXaxis().SetTitle("Reconstructed Neutrino Energy (GeV)")

if myvariable == "x":
	mc_new.GetXaxis().SetTitle("Reconstructed Bjorken x")

canvas1 = ROOT.TCanvas()
mnv = PlotUtils.MnvPlotter()

mnv.ApplyStyle(7)

mnv.error_summary_group_map.clear()
mnv.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_Resolution")
mnv.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_MINOS")
mnv.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_MINERvA")

mnv.error_summary_group_map["MINOS"].push_back("MINOS_Reconstruction_Efficiency")

mnv.error_summary_group_map["Muon Angle"].push_back("BeamAngleX")
mnv.error_summary_group_map["Muon Angle"].push_back("BeamAngleY")

mnv.error_summary_group_map["Flux"].push_back("Flux")

mnv.error_summary_group_map["Interaction Model"].push_back("RPA_HighQ2")
mnv.error_summary_group_map["Interaction Model"].push_back("RPA_LowQ2")
mnv.error_summary_group_map["Interaction Model"].push_back("Low_Recoil_2p2h_Tune")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_AhtBY")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_BhtBY")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_CCQEPauliSupViaKF")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_CV1uBY")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_CV2uBY")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_EtaNCEL")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_MaCCQE")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_MaNCEL")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_MaRES")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_MvRES")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_NormDISCC")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_NormNCRES")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvn1pi")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvn2pi")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvp1pi")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvp2pi")
mnv.error_summary_group_map["Interaction Model"].push_back("GENIE_VecFFCCQEshape")


mnv.error_summary_group_map["FSI"].push_back("GENIE_AGKYxF1pi")
mnv.error_summary_group_map["FSI"].push_back("GENIE_RDecBR1gamma")
mnv.error_summary_group_map["FSI"].push_back("GENIE_Theta_Delta2Npi")
mnv.error_summary_group_map["FSI"].push_back("GENIE_FrAbs_N")
mnv.error_summary_group_map["FSI"].push_back("GENIE_FrAbs_pi")
mnv.error_summary_group_map["FSI"].push_back("GENIE_FrCEx_N")
mnv.error_summary_group_map["FSI"].push_back("GENIE_FrCEx_pi")
mnv.error_summary_group_map["FSI"].push_back("GENIE_FrElas_N")
mnv.error_summary_group_map["FSI"].push_back("GENIE_FrElas_pi")
mnv.error_summary_group_map["FSI"].push_back("GENIE_FrInel_N")
mnv.error_summary_group_map["FSI"].push_back("GENIE_FrPiProd_N")
mnv.error_summary_group_map["FSI"].push_back("GENIE_FrPiProd_pi")
mnv.error_summary_group_map["FSI"].push_back("GENIE_MFP_N")
mnv.error_summary_group_map["FSI"].push_back("GENIE_MFP_pi")


mnv.error_summary_group_map["Target Mass"].push_back("Target_Mass")

mnv.DrawErrorSummary(mc_new, "TL", True, True, 0.0, False, "",True);
# last boolean decides whether frac or not

keys = canvas1.GetListOfPrimitives();
for k in keys:
	if(k.ClassName().find("Legend")!=-1):
		if myvariable == "Enu":
			k.SetNColumns(2)
			k.SetX2(50) #Enu
        	k.SetY1(0.20) #Enu

		if myvariable == "x":
			k.SetNColumns(2)
			k.SetX2(1.52) #x
			k.SetY1(0.18) #x

		#k.SetX2(40) #Emu
        #k.SetY1(0.28) #Emu
	if(k.ClassName().find("TH1")!=-1):
		if myvariable == "Enu":
			k.GetYaxis().SetRangeUser(0,0.3) # to 0.35 for Enu
		if myvariable == "x":
			k.GetYaxis().SetRangeUser(0,0.35)

mnv.AddHistoTitle("Tuned Iron Downstream ", 0.05, 1)

canvas1.Modified()
canvas1.Print("ME6A_T3Fe_TunedScaledDan_MC_DS_%s.png"%myvariable)

raw_input("Done")
