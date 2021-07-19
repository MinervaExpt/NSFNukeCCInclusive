import ROOT
import os,sys
from PlotUtils import MnvH1D,MnvH2D,MnvPlotter
_fnsf = ROOT.TFile("/minerva/data/users/zdar/MADHists/histsEventsDIS/histsEventsDIS/Hists_EventSelection_t1_z26_Nu_v1_.root","READ")
#_fnsf = ROOT.TFile("/minerva/data/users/afilkins/DIS_ME_NSF/2020-7-20_effvalidation/Hists_Efficiency_t1_z26_Nu_v1_.root","READ")
mc_name = "selected_mc_reco_x"
mc_new = _fnsf.Get(mc_name)


mnv = MnvPlotter()
canvas1 = ROOT.TCanvas()
mnv.error_summary_group_map.clear()
mnv.error_summary_group_map["Muon"].push_back("Muon_Energy_Resolution");
mnv.error_summary_group_map["Muon"].push_back("Muon_Energy_MINOS");
mnv.error_summary_group_map["Muon"].push_back("Muon_Energy_MINERvA");
mnv.error_summary_group_map["MINOS"].push_back("MINOS_Reconstruction_Efficiency");
mnv.error_summary_group_map["Det Res"].push_back("BeamAngleX");
mnv.error_summary_group_map["Det Res"].push_back("BeamAngleY");
mnv.error_summary_group_map["Flux"].push_back("Flux");
mnv.error_summary_group_map["RPA"].push_back("RPA_HighQ2");
mnv.error_summary_group_map["RPA"].push_back("RPA_LowQ2");
mnv.error_summary_group_map["GENIE"].push_back("Low_Recoil_2p2h_Tune");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrAbs_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrAbs_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrCEx_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrCEx_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrElas_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrElas_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrInel_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrInel_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrPiProd_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrPiProd_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_MFP_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_MFP_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_AGKYxF1pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_AhtBY");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_BhtBY");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_CCQEPauliSupViaKF");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_CV1uBY");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_CV2uBY");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_EtaNCEL");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_MaCCQE");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_NormCCQE");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_MaCCQEshape");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_MaNCEL");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_MaRES");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_MvRES");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_NormDISCC");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_NormCCRES");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_NormNCRES");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_RDecBR1gamma");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_Rvn1pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_Rvn2pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_Rvn3pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_Rvp1pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_Rvp2pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_Theta_Delta2Npi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_VecFFCCQEshape");
#mnv.DrawErrorSummary(mc_new, "TR", 0, 1, 0.0, 0, "GENIE");
#mnv.DrawErrorSummary(mc_new, "TR", 0, 1, 0.0, 0);
mnv.DrawErrorSummary(mc_new);
#mnv.DrawDataMCWithErrorBand(mc_new);
#mnv.DrawErrorSummary(mc_new, "TR", "false", "true", 0.0)
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "GENIE")
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "GENIE")
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "Flux")
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "GENIE_VecFFCCQEshape")
canvas1.Print("abc.png")
input("Press enter to end program")
