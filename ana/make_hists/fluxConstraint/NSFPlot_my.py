import ROOT
import os,sys
from PlotUtils import MnvH1D,MnvH2D,MnvPlotter
_fnsf = ROOT.TFile("/minerva/data2/users/anezkak/flux4Daisy_files_2022/OptimisedEventLoopME5A6ABCDEFGHJ_flux_grid_sys/flux_tracker.root","READ")
#_fnsf = ROOT.TFile("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/make_hists/fluxConstraint/flux_tracker.root","READ")

#_fnsf = ROOT.TFile("/minerva/data/users/afilkins/DIS_ME_NSF/2020-7-20_effvalidation/Hists_Efficiency_t1_z26_Nu_v1_.root","READ")
mc_name = "flux"
mc_new = _fnsf.Get(mc_name)
mc_new.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
canvas1 = ROOT.TCanvas()
mnv = MnvPlotter()

mnv.ApplyStyle(7)

mnv.error_summary_group_map.clear()

#mnv.error_summary_group_map["Flux"].push_back("Flux");


#mnv.DrawErrorSummary(mc_new, "TR", 0, 1, 0.0, 0, "GENIE");
#mnv.DrawErrorSummary(mc_new, "TR", 0, 1, 0.0, 0);
mnv.DrawErrorSummary(mc_new, "TL", True, True, 0.0, False, "",True);
"""
bool MnvPlotter::DrawErrorSummary(
        MnvH1D* h,
        const std::string& legPos  /* = "TR"    */,
        const bool   includeStat     /* = true    */,
        const bool   solidLinesOnly  /* = true    */,
        const double ignoreThreshold /* = 0.00001 */,
        const bool covAreaNormalize/* = false*/,
        const std::string& errorGroupName /* = "" */,
        const bool  asfrac  /* = false */,
        const std::string &Ytitle,
        bool ignoreUngrouped
        )
"""
keys = canvas1.GetListOfPrimitives();
for k in keys:
	if(k.ClassName().find("Legend")!=-1):
		k.SetNColumns(2)
		#k.SetX2(1.4) #x
		#k.SetY1(0.18) #x
		k.SetX2(44) #Enu
        	k.SetY1(0.2) #Enu

		#k.SetX2(40) #Emu
                #k.SetY1(0.28) #Emu
	if(k.ClassName().find("TH1")!=-1):
		k.GetYaxis().SetRangeUser(0,0.25)
canvas1.Modified()
#mnv.DrawDataMCWithErrorBand(mc_new);
#mnv.DrawErrorSummary(mc_new, "TR", "false", "true", 0.0)
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "GENIE")
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "GENIE")
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "Flux")
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "GENIE_VecFFCCQEshape")
mnv.AddHistoTitle("Tracker Flux", 0.04, 1)

canvas1.Print("FracUnc_AntiNuExceptME6I_Flux_Enu_Tracker.png")
input("Press enter to end program")
