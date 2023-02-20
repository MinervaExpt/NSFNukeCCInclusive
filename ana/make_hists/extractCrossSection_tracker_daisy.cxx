#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TArrayD.h"
#include "TParameter.h"

#include <iostream>
#include <stdlib.h>

#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/MacroUtil.h" 
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
        
#include "Math/Factory.h"     
#include "Math/Functor.h"     
#include "Math/Minimizer.h"   

#include "MinervaUnfold/MnvUnfold.h"

#include "../../NUKECCSRC/include/Binning.h"
#include "../../NUKECCSRC/include/UtilsNSF.h"
#include "../../NUKECCSRC/include/CVUniverse.h"
#include "../include/Variable.h"

#include "../include/systematics/Systematics.h"

#include "../include/plotting_functions.h"

#ifndef __CINT__
#endif

using namespace NUKECC_ANA;
using namespace PlotUtils;
using namespace std;

//Plot a step in cross section extraction.
void Plot(PlotUtils::MnvH1D& hist, const std::string& stepName, const std::string& prefix)
{
  TCanvas can(stepName.c_str());
  hist.GetCVHistoWithError().Clone()->Draw();
  can.Print((prefix + "_" + stepName + ".png").c_str());

  //Uncertainty summary
  PlotUtils::MnvPlotter plotter;
  plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
  plotter.axis_maximum = 0.4;

  plotter.DrawErrorSummary(&hist);
  can.Print((prefix + "_" + stepName + "_uncertaintySummary.png").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Other");
  can.Print((prefix + "_" + stepName + "_otherUncertainties.png").c_str());
}

//Unfolding function from Aaron Bercelle
//TODO: Trim it down a little?  Remove that static?
PlotUtils::MnvH1D* UnfoldHist( PlotUtils::MnvH1D* h_folded, PlotUtils::MnvH2D* h_migration, int num_iter )
{
  static MinervaUnfold::MnvUnfold unfold;
  PlotUtils::MnvH1D* h_unfolded = nullptr;

  //bool bUnfolded = false;

  TMatrixD dummyCovMatrix;
  if(!unfold.UnfoldHisto( h_unfolded, dummyCovMatrix, h_migration, h_folded, RooUnfold::kBayes, num_iter, true, false ))
    return nullptr;

  /////////////////////////////////////////////////////////////////////////////////////////  
  //No idea if this is still needed
  //Probably.  This gets your stat unfolding covariance matrix
  TMatrixD unfoldingCovMatrixOrig; 
  int correctNbins;
  int matrixRows;  
  TH1D* hUnfoldedDummy  = new TH1D(h_unfolded->GetCVHistoWithStatError());
  TH1D* hRecoDummy      = new TH1D(h_migration->ProjectionX()->GetCVHistoWithStatError());
  TH1D* hTruthDummy     = new TH1D(h_migration->ProjectionY()->GetCVHistoWithStatError());
  TH1D* hBGSubDataDummy = new TH1D(h_folded->GetCVHistoWithStatError());
  TH2D* hMigrationDummy = new TH2D(h_migration->GetCVHistoWithStatError());

  // Perform unfolding
  unfold.UnfoldHisto(hUnfoldedDummy, unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy,RooUnfold::kBayes, num_iter);//Stupid RooUnfold.  This is dummy, we don't need iterations

  correctNbins=hUnfoldedDummy->fN;
  matrixRows=unfoldingCovMatrixOrig.GetNrows();
  if(correctNbins!=matrixRows){
    std::cout << "****************************************************************************" << std::endl;
    std::cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << std::endl;
    std::cout << "****************************************************************************" << std::endl;
    // It looks like this, since the extra last two bins don't have any content
    unfoldingCovMatrixOrig.ResizeTo(correctNbins, correctNbins);
  }

  for(int i=0; i<unfoldingCovMatrixOrig.GetNrows(); ++i) unfoldingCovMatrixOrig(i,i)=0;
  // Release memory
  delete hUnfoldedDummy;
  delete hMigrationDummy;
  delete hRecoDummy;
  delete hTruthDummy;
  delete hBGSubDataDummy;
  h_unfolded->PushCovMatrix("unfoldingCov",unfoldingCovMatrixOrig);
  unfoldingCovMatrixOrig.Print();

  /////////////////////////////////////////////////////////////////////////////////////////  
  return h_unfolded;
}

//  Write covariance and correlation marices
// Inspired by Gonzalo Diaz's function
void WriteCovAndCorrMatrices(const auto&var, PlotUtils::MnvH1D* unfolded_histo, TFile* fout, const auto&prefix)
{
  fout->cd();

  // Unfolding covariance matrix
  // It’s preferably that, once extracted, you can check unfold_cov_matrix. 
  // If your unfolding has been done correctly, it will be a symmetric matrix with no entries in the diagonal.
  // Get syst. error matrix
  // 1st bool: Show as fractional error?
  // 2nd bool: Area-norm covariance?
  TMatrixD unfold_cov_matrix = unfolded_histo->GetSysErrorMatrix("unfoldingCov", false, false);
  TH2D* h_unfold_cov_matrix =  new TH2D(unfold_cov_matrix);
  h_unfold_cov_matrix->Write(Form("unfold_cov_matrix_%s_%s", prefix.c_str(), var.c_str() ));

  // Statistical error matrix
  // bool: Show as fractional error?
  TMatrixD stat_err_matrix = unfolded_histo->GetStatErrorMatrix(false);
  TH2D* h_stat_err_matrix =  new TH2D(stat_err_matrix);
  h_stat_err_matrix->Write(Form("stat_err_matrix_%s_%s", prefix.c_str(), var.c_str() ));

  // Statistical covariance matrix
  // = sum of the unfolding covariance and the stat-only error; not fully diagonal
  // includes the statistical errors of each bins AND the effect of the unfolding
  // between different bins, which in fact adds non-diagonal terms.
  TMatrixD stat_cov_matrix = unfold_cov_matrix + stat_err_matrix;
  TH2D* h_stat_cov_matrix =  new TH2D(stat_cov_matrix);
  h_stat_cov_matrix->Write(Form("stat_cov_matrix_%s_%s", prefix.c_str(), var.c_str() ));
  
  // Statistical correlation matrix
  const int size = stat_cov_matrix.GetNrows();
  TMatrixD stat_corr_matrix(size,size);
  for ( int x = 0; x < size; ++x ) { //loop over columns
    for ( int y = 0; y < size; ++y ) { // loop over rows
      stat_corr_matrix[x][y] = (stat_cov_matrix[x][x] == 0.0 || stat_cov_matrix[y][y] == 0.0 ) ?
                                0.0 : stat_cov_matrix[x][y] / std::sqrt(stat_cov_matrix[x][x]*stat_cov_matrix[y][y]);
    }
  }
  TH2D* h_stat_corr_matrix =  new TH2D(stat_corr_matrix);
  h_stat_corr_matrix->Write(Form("stat_corr_matrix_%s_%s", prefix.c_str(), var.c_str() ));

  // TOTAL covariance and correlations
  
  // Write total error matrix
  // 1st bool: Include stat. error?
  // 2nd bool: Show as fractional error?
  // 3rd bool: Area-norm covariance?
  TH2D* h_cov_matrix = new TH2D(unfolded_histo->GetTotalErrorMatrix(true, false, false));
  h_cov_matrix->Write(Form("total_cov_matrix_%s_%s", prefix.c_str(), var.c_str() ));

  // Write total correlation matrix
  // 1st bool: Area-norm covariance?
  // 2nd bool: Include stat. error?
  TH2D* h_corr_matrix = new TH2D(unfolded_histo->GetTotalCorrelationMatrixTH2D(false, true));
  h_corr_matrix->Write(Form("total_corr_matrix_%s_%s", prefix.c_str(), var.c_str() ));
  
  delete h_unfold_cov_matrix;
  delete h_stat_err_matrix;
  delete h_stat_cov_matrix;
  delete h_stat_corr_matrix;
  delete h_cov_matrix;
  delete h_corr_matrix;

}

double GetTotalScatteringCenters(int targetZ, bool isMC)
{
  // TARGET INFO
  PlotUtils::TargetUtils targetInfo;
  double Nucleons;

  // Target 1 is generally excluded due to rock muon contamination (in the inclusive analysis)
  if(targetZ == 6){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ); // Target 3
  }
  
  if(targetZ == 26){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 2, targetZ, isMC ) // Target 2
             + targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ) // Target 3
             + targetInfo.GetPassiveTargetNNucleons( 5, targetZ, isMC );// Target 5
  }

  if(targetZ == 82){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 2, targetZ, isMC ) // Target 2
             + targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ) // Target 3
             + targetInfo.GetPassiveTargetNNucleons( 4, targetZ, isMC ) // Target 4
             + targetInfo.GetPassiveTargetNNucleons( 5, targetZ, isMC );// Target 5
  }
  if(targetZ > 90 ){
    Nucleons = targetInfo.GetTrackerNNucleons(5980, 8422, isMC, 850);
    //double TargetUtils::GetTrackerNNucleons( double minZ, double maxZ, bool isMC, double apothem /* = 850. */ ) const

  }

  return Nucleons;     
}


//The final step of cross section extraction: normalize by flux, bin width, POT, and number of targets
// Differential cross-section
PlotUtils::MnvH1D* normalize(PlotUtils::MnvH1D* efficiencyCorrected, PlotUtils::MnvH1D* fluxIntegral, const double nNucleons, const double POT)
{
  efficiencyCorrected->Divide(efficiencyCorrected, fluxIntegral);
  efficiencyCorrected->Scale(1./nNucleons/POT);
  efficiencyCorrected->Scale(1.e4); //Flux histogram is in m^-2, but convention is to report cm^2
  efficiencyCorrected->Scale(1., "width");

  return efficiencyCorrected;
}

// Total cross-section
PlotUtils::MnvH1D* normalizeTotal(PlotUtils::MnvH1D* efficiencyCorrected, PlotUtils::MnvH1D* fluxRebinned, const double nNucleons, const double POT)
{
  efficiencyCorrected->Divide(efficiencyCorrected, fluxRebinned);
  efficiencyCorrected->Scale(1./nNucleons/POT);
  efficiencyCorrected->Scale(1.e4); //Flux histogram is in m^-2, but convention is to report cm^2

  return efficiencyCorrected;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------
// MAIN FUNCTION
// -------------------------------------------------------------------------------------------------------------------------------------------------

//int UnfoldIterations( const std::string& var );

int main(int argc, char * argv[]){
  TH1::AddDirectory(false);

  if(argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------\
-----"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./runEventLoop Path_to_Output_file \n\n"<<
      "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
      "\t-Target_number\t \t = \t Number of target, e.g. 1-5, or tracker 99 \n"<<
      "\t-Material_atomic_number\t =\t Atomic number of material, e.g. 6, 26, 82, or tracker 99  \n" <<
      "\t-Playlist\t \t =\t eg. minervame1A"<< std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    return 0;
  }

  string outdir = argv[1];
  int targetID = atoi(argv[2]);
  int targetZ = atoi(argv[3]);
  const string playlist= argv[4];
  const std::string plist_string(playlist);

  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(-14);


  TString eventLoop, migr, efficiency;

  //bool RunCodeWithSystematics = false;
  // Read in background subtracted event rate, Migration and Efficiency
  if(RunCodeWithSystematics){
        eventLoop = Form("%s/BackgroundSubtracted/BkgSubtracted_EventSelection_daisy_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetID, targetZ);
        migr = Form("%s/Migration/Migration_Daisy_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetID, targetZ);
        efficiency = Form("%s/Efficiency/Efficiency_Daisy_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetID, targetZ);
    }
  else{
      eventLoop = Form("%s/BackgroundSubtracted/BkgSubtracted_EventSelection_daisy_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetID, targetZ);
      migr = Form("%s/Migration/Migration_Daisy_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetID, targetZ);
      efficiency = Form("%s/Efficiency/Efficiency_Daisy_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetID, targetZ);
  }

  TFile *fEventLoop = new TFile( eventLoop,"read" );
  TFile *fMigration = new TFile( migr,"read" );
  TFile *fEfficiency = new TFile( efficiency,"read" );
  
  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF("minervame6A");
  NukeCC_Binning  *binsDef = new NukeCC_Binning();
  HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist("minervame6A");

  // to iterate over variables
  std::vector<string> vars;
  vars.push_back("Enu");
  vars.push_back("x");
  vars.push_back("pTmu");
  vars.push_back("pZmu");
  
  //string var = "Enu";

  // to iterate over mc and data
  std::vector<string> prefixes;
  prefixes.push_back("mc");
  prefixes.push_back("data");
  
  // output file
  TFile *fUnfold = new TFile( Form("%s/CrossSection_Daisy_t%02d_z%02d_%s.root", outdir.c_str(), targetID, targetZ, plist_string.c_str() ), "recreate" );
   std::cout<<"Hi"<<std::endl;

  // read in the POT information
  TParameter<double> *mcPOT = (TParameter<double>*)fEventLoop->Get("MCPOT");
  TParameter<double> *dataPOT = (TParameter<double>*)fEventLoop->Get("DataPOT");
  double mcpot = mcPOT->GetVal();
  double datapot = dataPOT->GetVal();

  fUnfold->cd();
  // write the POT information into the output file
  auto dataPOTOut = new TParameter<double>("DataPOT", datapot);
  auto mcPOTOut = new TParameter<double>("MCPOT", mcpot);
  dataPOTOut->Write();
  mcPOTOut->Write(); 

  double dataMCScale = datapot/mcpot;
  cout<<"Data POT = "<<datapot<<" MC POT = "<<mcpot<<"   Data/MC POT = "<<dataMCScale<<endl;

  for(const auto&var: vars){

    // read in histogram for each petal (variable and prefix)
    // unfold
    // efficiency correct
    for (int petal=0; petal<12; petal++){

      // Var specific ingredients
      auto migration  = dynamic_cast<MnvH2D*>(fMigration->Get( Form("selected_Migration_daisy_%d_%s", petal, var.c_str())));
      auto effNum = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_mc_daisy_%d_%s", petal, var.c_str())));
      auto effDenom = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_daisy_%d_%s", petal, var.c_str())));
      auto effDenom_QE = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_daisy_QE_%d_%s", petal, var.c_str())));
      auto effDenom_RES = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_daisy_RES_%d_%s", petal, var.c_str())));
      auto effDenom_DIS = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_daisy_DIS_%d_%s", petal, var.c_str())));
      auto effDenom_Other = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_daisy_Other_%d_%s", petal, var.c_str())));
      auto effDenom_2p2h = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_daisy_2p2h_%d_%s", petal, var.c_str())));
      auto simEventRate = effDenom->Clone();
      auto simEventRate_QE = effDenom_QE->Clone();
      auto simEventRate_RES = effDenom_RES->Clone();
      auto simEventRate_DIS = effDenom_DIS->Clone();
      auto simEventRate_Other = effDenom_Other->Clone();
      auto simEventRate_2p2h = effDenom_2p2h->Clone();
      fUnfold->cd();

      simEventRate->Clone()->Write(Form("simEventRate_daisy_%d_%s", petal, var.c_str() ));
      simEventRate_QE->Clone()->Write(Form("simEventRate_QE_daisy_%d_%s", petal, var.c_str() ));
      simEventRate_RES->Clone()->Write(Form("simEventRate_RES_daisy_%d_%s", petal, var.c_str() ));
      simEventRate_DIS->Clone()->Write(Form("simEventRate_DIS_daisy_%d_%s", petal, var.c_str() ));
      simEventRate_Other->Clone()->Write(Form("simEventRate_Other_daisy_%d_%s", petal, var.c_str() ));
      simEventRate_2p2h->Clone()->Write(Form("simEventRate_2p2h_daisy_%d_%s", petal, var.c_str() ));

      effNum->Clone()->Write(Form("efficiency_numerator_daisy_%d_%s", petal, var.c_str() ));

      // EFFICIENCY
      auto efficiency = effNum->Clone();
      efficiency->Divide(efficiency, effDenom, 1.0, 1.0, "B"); //Only the 2 parameter version of MnvH1D::Divide() //handles systematics correctly.
      //Plot(*effNum, "efficiency", prefix);
      efficiency->Clone()->Write(Form("efficiency_daisy_%d_%s", petal, var.c_str() )); 

      // Var and prefix specific ingredients
      for(const auto& prefix: prefixes){

        auto bkgSubtracted =  dynamic_cast<MnvH1D*>(fEventLoop->Get( Form("h_bkg_subtracted_%s_daisy_%d_%s", prefix.c_str(), petal, var.c_str())));

        cout<<"I am here"<<endl;
        //d'Aogstini unfolding
        std::cout << "===================================== \n" << std::flush;
        std::cout << "============= UNFOLDING =============\n" << std::flush;
        int nIterations=3;    
        std::cout << "Number of iterations: " << nIterations << "\n" << std::flush;

        // UNFOLDING
        auto unfolded = UnfoldHist(bkgSubtracted, migration, nIterations);
        if(!unfolded) throw std::runtime_error(std::string("Failed to unfold ") + bkgSubtracted->GetName() + " using " + migration->GetName());
        //Plot(*unfolded_mc, "unfolded", prefix);
        unfolded->Clone()->Write(Form("unfolded_%s_daisy_%d_%s", prefix.c_str(), petal, var.c_str() ));
        // histogram unfolded data stores the diagonals, "unfoldinCov" contains the covariant values
        // in unfoldingCov, the diagonals are zero, but there should be off-diagonal content

        //unfolded_data->GetStatErrorMatrix(false) + unfolded_data->GetSysCorrelationMatrix("unfoldingCov", false);
        //unfolded_data->Clone()->Write(Form("unfolded_data_covMatrix_%s", var.c_str() ));

        std::cout << "Survived writing the unfolded histogram.\n" << std::flush;
        bool isMC = false;
        if (prefix=="mc"){
          isMC = true;
        }

        WriteCovAndCorrMatrices(var, unfolded, fUnfold, prefix);
        std::cout << "Wrote covariance and correlation matrices.\n" << std::flush;

        // CROSS-SECTION EXTRACTION

        // Unfolded efficiency corrected
        unfolded->Divide(unfolded, efficiency);
        //Plot(*unfolded_data, "efficiencyCorrected", prefix);
        unfolded->Clone()->Write(Form("unfolded_effCorrected_%s_daisy_%d_%s", prefix.c_str(), petal, var.c_str() )); 
        std::cout << "Efficiency corrected.\n" << std::flush;
      }
    }

    // ---------------------------------------------------------------------------
    // Combine efficiency corrected unfolded petal distributions for tracker here
    // ---------------------------------------------------------------------------

    for(const auto& prefix: prefixes){  
      bool isMC = false;
      if (prefix=="mc"){
        isMC = true;
      }
    

      // ---------------------------------------------------------------------
      // Flux reweighter information, get reweighted daisy sum according to a material
      // ---------------------------------------------------------------------
      int n_flux_universes =  PlotUtils::MinervaUniverse::GetNFluxUniverses();
      int nu_pdg = PlotUtils::MinervaUniverse::GetAnalysisNuPDG();
      const bool use_nue_constraint = true;
      bool useMuonCorrelations = true;
      const std::string project_dir = "targets_2345_jointNueIMD";
      double min_energy = 0;
      double max_energy = 120;

      auto& frw = PlotUtils::flux_reweighter("minervame6A", nu_pdg, use_nue_constraint, n_flux_universes);

      //map of petal distributions
      std::map<int, MnvH1D*> daisy_petal_hists;
      std::map<int, MnvH1D*> daisy_petal_hists_QE;
      std::map<int, MnvH1D*> daisy_petal_hists_RES;
      std::map<int, MnvH1D*> daisy_petal_hists_DIS;
      std::map<int, MnvH1D*> daisy_petal_hists_Other;
      std::map<int, MnvH1D*> daisy_petal_hists_2p2h;
      MnvH1D* h_eff_corr_tracker_to_carbon;
      MnvH1D* h_eff_corr_tracker_to_carbon_QE;
      MnvH1D* h_eff_corr_tracker_to_carbon_RES;
      MnvH1D* h_eff_corr_tracker_to_carbon_DIS;
      MnvH1D* h_eff_corr_tracker_to_carbon_Other;
      MnvH1D* h_eff_corr_tracker_to_carbon_2p2h;
      MnvH1D* h_eff_corr_tracker_to_iron;
      MnvH1D* h_eff_corr_tracker_to_iron_QE;
      MnvH1D* h_eff_corr_tracker_to_iron_RES;
      MnvH1D* h_eff_corr_tracker_to_iron_DIS;
      MnvH1D* h_eff_corr_tracker_to_iron_Other;
      MnvH1D* h_eff_corr_tracker_to_iron_2p2h;
      MnvH1D* h_eff_corr_tracker_to_lead;
      MnvH1D* h_eff_corr_tracker_to_lead_QE;
      MnvH1D* h_eff_corr_tracker_to_lead_RES;
      MnvH1D* h_eff_corr_tracker_to_lead_DIS;
      MnvH1D* h_eff_corr_tracker_to_lead_Other;
      MnvH1D* h_eff_corr_tracker_to_lead_2p2h;

      // Efficiency corrected distribution to targets
      if(prefix == "mc"){
        for (int petal=0; petal<12; petal++){
          // MC after unfolding is actually efficiency denominator (otherwise we introduce additional stat uncertainty from unfolding)
          daisy_petal_hists[petal] = dynamic_cast<MnvH1D*>(fUnfold->Get(Form("simEventRate_daisy_%d_%s", petal, var.c_str() )));
          daisy_petal_hists_QE[petal] = dynamic_cast<MnvH1D*>(fUnfold->Get(Form("simEventRate_QE_daisy_%d_%s", petal, var.c_str() )));
          daisy_petal_hists_RES[petal] = dynamic_cast<MnvH1D*>(fUnfold->Get(Form("simEventRate_RES_daisy_%d_%s", petal, var.c_str() )));
          daisy_petal_hists_DIS[petal] = dynamic_cast<MnvH1D*>(fUnfold->Get(Form("simEventRate_DIS_daisy_%d_%s", petal, var.c_str() )));
          daisy_petal_hists_Other[petal] = dynamic_cast<MnvH1D*>(fUnfold->Get(Form("simEventRate_Other_daisy_%d_%s", petal, var.c_str() )));
          daisy_petal_hists_2p2h[petal] = dynamic_cast<MnvH1D*>(fUnfold->Get(Form("simEventRate_2p2h_daisy_%d_%s", petal, var.c_str() )));
        }
        // Efficiency corrected distribution to targets
        //Carbon
        h_eff_corr_tracker_to_carbon = frw.GetReweightedDaisySum(nu_pdg, "carbon", daisy_petal_hists, project_dir );
        h_eff_corr_tracker_to_carbon_QE = frw.GetReweightedDaisySum(nu_pdg, "carbon", daisy_petal_hists_QE, project_dir );
        h_eff_corr_tracker_to_carbon_RES = frw.GetReweightedDaisySum(nu_pdg, "carbon", daisy_petal_hists_RES, project_dir );
        h_eff_corr_tracker_to_carbon_DIS = frw.GetReweightedDaisySum(nu_pdg, "carbon", daisy_petal_hists_DIS, project_dir );
        h_eff_corr_tracker_to_carbon_Other = frw.GetReweightedDaisySum(nu_pdg, "carbon", daisy_petal_hists_Other, project_dir );
        h_eff_corr_tracker_to_carbon_2p2h = frw.GetReweightedDaisySum(nu_pdg, "carbon", daisy_petal_hists_2p2h, project_dir );
        //Iron
        h_eff_corr_tracker_to_iron = frw.GetReweightedDaisySum(nu_pdg, "iron", daisy_petal_hists, project_dir );
        h_eff_corr_tracker_to_iron_QE = frw.GetReweightedDaisySum(nu_pdg, "iron", daisy_petal_hists_QE, project_dir );
        h_eff_corr_tracker_to_iron_RES = frw.GetReweightedDaisySum(nu_pdg, "iron", daisy_petal_hists_RES, project_dir );
        h_eff_corr_tracker_to_iron_DIS = frw.GetReweightedDaisySum(nu_pdg, "iron", daisy_petal_hists_DIS, project_dir );
        h_eff_corr_tracker_to_iron_Other = frw.GetReweightedDaisySum(nu_pdg, "iron", daisy_petal_hists_Other, project_dir );
        h_eff_corr_tracker_to_iron_2p2h = frw.GetReweightedDaisySum(nu_pdg, "iron", daisy_petal_hists_2p2h, project_dir );
        //Lead
        h_eff_corr_tracker_to_lead = frw.GetReweightedDaisySum(nu_pdg, "lead", daisy_petal_hists, project_dir );
        h_eff_corr_tracker_to_lead_QE = frw.GetReweightedDaisySum(nu_pdg, "lead", daisy_petal_hists_QE, project_dir );
        h_eff_corr_tracker_to_lead_RES = frw.GetReweightedDaisySum(nu_pdg, "lead", daisy_petal_hists_RES, project_dir );
        h_eff_corr_tracker_to_lead_DIS = frw.GetReweightedDaisySum(nu_pdg, "lead", daisy_petal_hists_DIS, project_dir );
        h_eff_corr_tracker_to_lead_Other = frw.GetReweightedDaisySum(nu_pdg, "lead", daisy_petal_hists_Other, project_dir );
        h_eff_corr_tracker_to_lead_2p2h = frw.GetReweightedDaisySum(nu_pdg, "lead", daisy_petal_hists_2p2h, project_dir );
      }
      
      else{
        for (int petal=0; petal<12; petal++){
          daisy_petal_hists[petal] = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("unfolded_effCorrected_%s_daisy_%d_%s", prefix.c_str(), petal, var.c_str() )));
        }
        // Efficiency corrected distribution to targets
        h_eff_corr_tracker_to_carbon = frw.GetReweightedDaisySum(nu_pdg, "carbon", daisy_petal_hists, project_dir );
        h_eff_corr_tracker_to_iron = frw.GetReweightedDaisySum(nu_pdg, "iron", daisy_petal_hists, project_dir );
        h_eff_corr_tracker_to_lead = frw.GetReweightedDaisySum(nu_pdg, "lead", daisy_petal_hists, project_dir );
      }

      //std::cout << "Here.\n" << std::flush;

      //auto h_eff_corr_tracker_to_lead = frw.GetReweightedDaisySum(nu_pdg, "lead", daisy_petal_hists, project_dir );
      //auto h_eff_corr_tracker_to_iron = frw.GetReweightedDaisySum(nu_pdg, "iron", daisy_petal_hists, project_dir );

      // Write efficiency corrected distributions to targets to the file
      fUnfold->cd();
      // carbon
      h_eff_corr_tracker_to_carbon->Clone()->Write(Form("h_eff_corr_tracker_to_carbon_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_carbon_QE->Clone()->Write(Form("h_eff_corr_tracker_to_carbon_QE_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_carbon_RES->Clone()->Write(Form("h_eff_corr_tracker_to_carbon_RES_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_carbon_DIS->Clone()->Write(Form("h_eff_corr_tracker_to_carbon_DIS_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_carbon_Other->Clone()->Write(Form("h_eff_corr_tracker_to_carbon_Other_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_carbon_2p2h->Clone()->Write(Form("h_eff_corr_tracker_to_carbon_2p2h_%s_%s", prefix.c_str(),var.c_str() ));
      // iron
      h_eff_corr_tracker_to_iron->Clone()->Write(Form("h_eff_corr_tracker_to_iron_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_iron_QE->Clone()->Write(Form("h_eff_corr_tracker_to_iron_QE_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_iron_RES->Clone()->Write(Form("h_eff_corr_tracker_to_iron_RES_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_iron_DIS->Clone()->Write(Form("h_eff_corr_tracker_to_iron_DIS_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_iron_Other->Clone()->Write(Form("h_eff_corr_tracker_to_iron_Other_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_iron_2p2h->Clone()->Write(Form("h_eff_corr_tracker_to_iron_2p2h_%s_%s", prefix.c_str(),var.c_str() ));
      //lead
      h_eff_corr_tracker_to_lead->Clone()->Write(Form("h_eff_corr_tracker_to_lead_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_lead_QE->Clone()->Write(Form("h_eff_corr_tracker_to_lead_QE_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_lead_RES->Clone()->Write(Form("h_eff_corr_tracker_to_lead_RES_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_lead_DIS->Clone()->Write(Form("h_eff_corr_tracker_to_lead_DIS_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_lead_Other->Clone()->Write(Form("h_eff_corr_tracker_to_lead_Other_%s_%s", prefix.c_str(),var.c_str() ));
      h_eff_corr_tracker_to_lead_2p2h->Clone()->Write(Form("h_eff_corr_tracker_to_lead_2p2h_%s_%s", prefix.c_str(),var.c_str() ));


      std::map<int, MnvH1D*> h_eff_corr_tracker_to_material;
      h_eff_corr_tracker_to_material[0] = h_eff_corr_tracker_to_carbon;
      h_eff_corr_tracker_to_material[1] = h_eff_corr_tracker_to_iron;
      h_eff_corr_tracker_to_material[2] = h_eff_corr_tracker_to_lead;

      // To iterate over materials now
      std::vector<string> materials;
      materials.push_back("carbon");
      materials.push_back("iron");
      materials.push_back("lead");


      int iter = 0;
      for(const auto&mat: materials){
        //std::cout<< mat.c_str() << std::endl;
        std::cout<< materials[iter].c_str() << std::endl;
        // Integrated target flux
        //GetIntegratedFluxReweighted arguments
        //int nuPDG
        //std::string tar_mat
        //MnvHistoType* template_hist
        //double min_energy
        //double max_energy
        //std::string project_dir
        auto fluxIntegral = frw.GetIntegratedTargetFlux(nu_pdg, mat.c_str(), h_eff_corr_tracker_to_material[iter], min_energy, max_energy, project_dir);

        auto flux = frw.GetTargetFluxMnvH1D(nu_pdg, mat.c_str(), project_dir);
        auto fluxRebinned = frw.GetRebinnedFluxReweighted_FromInputFlux(flux, h_eff_corr_tracker_to_material[iter]);

        if (prefix == "mc"){
          fUnfold->cd();
          flux->Clone()->Write(Form("flux_%s_%s", mat.c_str(), var.c_str()));
          fluxIntegral->Clone()->Write(Form("flux_integral_%s_%s", mat.c_str(), var.c_str()));
          fluxRebinned->Clone()->Write(Form("fluxRebinned_%s_%s", mat.c_str(), var.c_str()));
        }
      

        
        PlotUtils::TargetUtils targetInfo;
        //double nNucleons = targetInfo.GetPassiveTargetNNucleons( targetID, targetZ, isMC );
        double nNucleons = GetTotalScatteringCenters(targetZ, isMC);

        std::cout << prefix + " number of nucleons" << std::endl;
        std::cout << nNucleons << std::endl;

        double pot = datapot;

        std::cout << prefix + " POT" << std::endl;
        std::cout << pot << std::endl;
        
        if (prefix=="mc"){
        // simulated event rate from efficiency denominator -> cross-section (to check with GENIE xSec)
        // MC cross-section calculated from the efficiency numerator
        // Unfolding inflates the statistical uncertainty and will give a reader the wrong impression of the size of the MC sample.
          if(plist_string.size() > 11){
          pot = datapot; // if COMBINED histos, then scale by dataPOT!! because MCs already scaled by mcScale when added
          }
          else{
            pot = mcpot;
          }
          if (var=="Enu"){ // calculated simulated total cross-section for Enu
            auto simEventRate_cross_total = h_eff_corr_tracker_to_material[iter]->Clone();
            simEventRate_cross_total = normalizeTotal(simEventRate_cross_total, fluxRebinned, nNucleons, pot);
            fUnfold->cd();
            simEventRate_cross_total->Clone()->Write(Form("simEventRate_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));

            if(mat=="carbon"){
              // QE
              auto simEventRate_QE_cross_total = h_eff_corr_tracker_to_carbon_QE->Clone();
              simEventRate_QE_cross_total = normalizeTotal(simEventRate_QE_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_QE_cross_total->Clone()->Write(Form("simEventRate_QE_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // RES
              auto simEventRate_RES_cross_total = h_eff_corr_tracker_to_carbon_RES->Clone();
              simEventRate_RES_cross_total = normalizeTotal(simEventRate_RES_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_RES_cross_total->Clone()->Write(Form("simEventRate_RES_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // DIS
              auto simEventRate_DIS_cross_total = h_eff_corr_tracker_to_carbon_DIS->Clone();
              simEventRate_DIS_cross_total = normalizeTotal(simEventRate_DIS_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_DIS_cross_total->Clone()->Write(Form("simEventRate_DIS_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // Other
              auto simEventRate_Other_cross_total = h_eff_corr_tracker_to_carbon_Other->Clone();
              simEventRate_Other_cross_total = normalizeTotal(simEventRate_Other_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_Other_cross_total->Clone()->Write(Form("simEventRate_Other_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // 2p2h
              auto simEventRate_2p2h_cross_total = h_eff_corr_tracker_to_carbon_2p2h->Clone();
              simEventRate_2p2h_cross_total = normalizeTotal(simEventRate_2p2h_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_2p2h_cross_total->Clone()->Write(Form("simEventRate_2p2h_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
            }
            if(mat=="iron"){
              // QE
              auto simEventRate_QE_cross_total = h_eff_corr_tracker_to_iron_QE->Clone();
              simEventRate_QE_cross_total = normalizeTotal(simEventRate_QE_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_QE_cross_total->Clone()->Write(Form("simEventRate_QE_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // RES
              auto simEventRate_RES_cross_total = h_eff_corr_tracker_to_iron_RES->Clone();
              simEventRate_RES_cross_total = normalizeTotal(simEventRate_RES_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_RES_cross_total->Clone()->Write(Form("simEventRate_RES_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // DIS
              auto simEventRate_DIS_cross_total = h_eff_corr_tracker_to_iron_DIS->Clone();
              simEventRate_DIS_cross_total = normalizeTotal(simEventRate_DIS_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_DIS_cross_total->Clone()->Write(Form("simEventRate_DIS_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // Other
              auto simEventRate_Other_cross_total = h_eff_corr_tracker_to_iron_Other->Clone();
              simEventRate_Other_cross_total = normalizeTotal(simEventRate_Other_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_Other_cross_total->Clone()->Write(Form("simEventRate_Other_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // 2p2h
              auto simEventRate_2p2h_cross_total = h_eff_corr_tracker_to_iron_2p2h->Clone();
              simEventRate_2p2h_cross_total = normalizeTotal(simEventRate_2p2h_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_2p2h_cross_total->Clone()->Write(Form("simEventRate_2p2h_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));


            }
            if(mat=="lead"){
              auto simEventRate_QE_cross_total = h_eff_corr_tracker_to_lead_QE->Clone();
              simEventRate_QE_cross_total = normalizeTotal(simEventRate_QE_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_QE_cross_total->Clone()->Write(Form("simEventRate_QE_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // RES
              auto simEventRate_RES_cross_total = h_eff_corr_tracker_to_lead_RES->Clone();
              simEventRate_RES_cross_total = normalizeTotal(simEventRate_RES_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_RES_cross_total->Clone()->Write(Form("simEventRate_RES_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // DIS
              auto simEventRate_DIS_cross_total = h_eff_corr_tracker_to_lead_DIS->Clone();
              simEventRate_DIS_cross_total = normalizeTotal(simEventRate_DIS_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_DIS_cross_total->Clone()->Write(Form("simEventRate_DIS_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // Other
              auto simEventRate_Other_cross_total = h_eff_corr_tracker_to_lead_Other->Clone();
              simEventRate_Other_cross_total = normalizeTotal(simEventRate_Other_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_Other_cross_total->Clone()->Write(Form("simEventRate_Other_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
              // 2p2h
              auto simEventRate_2p2h_cross_total = h_eff_corr_tracker_to_lead_2p2h->Clone();
              simEventRate_2p2h_cross_total = normalizeTotal(simEventRate_2p2h_cross_total, fluxRebinned, nNucleons, pot);
              fUnfold->cd();
              simEventRate_2p2h_cross_total->Clone()->Write(Form("simEventRate_2p2h_crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));

            }

            //auto crossSection_total_mc = effNum->Clone();
            //crossSection_total_mc= normalizeTotal(crossSection_total_mc, fluxRebinned, nNucleons, pot);
            //crossSection_total_mc->Clone()->Write(Form("crossSection_total_%s_%s", prefix.c_str(), var.c_str() ));
            //std::cout << "MC total cross-section  DONE.\n" << std::flush;

          }
          auto simEventRate_cross = h_eff_corr_tracker_to_material[iter]->Clone();
          simEventRate_cross = normalize(simEventRate_cross , fluxIntegral, nNucleons, pot);
          fUnfold->cd();
          simEventRate_cross->Clone()->Write(Form("simEventRate_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));

          if(mat=="carbon"){
            // QE
            auto simEventRate_QE_cross = h_eff_corr_tracker_to_carbon_QE->Clone();
            simEventRate_QE_cross = normalize(simEventRate_QE_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_QE_cross->Clone()->Write(Form("simEventRate_QE_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // RES
            auto simEventRate_RES_cross = h_eff_corr_tracker_to_carbon_RES->Clone();
            simEventRate_RES_cross = normalize(simEventRate_RES_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_RES_cross->Clone()->Write(Form("simEventRate_RES_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // DIS
            auto simEventRate_DIS_cross = h_eff_corr_tracker_to_carbon_DIS->Clone();
            simEventRate_DIS_cross = normalize(simEventRate_DIS_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_DIS_cross->Clone()->Write(Form("simEventRate_DIS_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // Other
            auto simEventRate_Other_cross = h_eff_corr_tracker_to_carbon_Other->Clone();
            simEventRate_Other_cross = normalize(simEventRate_Other_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_Other_cross->Clone()->Write(Form("simEventRate_Other_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // 2p2h
            auto simEventRate_2p2h_cross = h_eff_corr_tracker_to_carbon_2p2h->Clone();
            simEventRate_2p2h_cross = normalize(simEventRate_2p2h_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_2p2h_cross->Clone()->Write(Form("simEventRate_2p2h_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
          }
          if(mat=="iron"){
            // QE
            auto simEventRate_QE_cross = h_eff_corr_tracker_to_iron_QE->Clone();
            simEventRate_QE_cross = normalize(simEventRate_QE_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_QE_cross->Clone()->Write(Form("simEventRate_QE_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // RES
            auto simEventRate_RES_cross = h_eff_corr_tracker_to_iron_RES->Clone();
            simEventRate_RES_cross = normalize(simEventRate_RES_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_RES_cross->Clone()->Write(Form("simEventRate_RES_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // DIS
            auto simEventRate_DIS_cross = h_eff_corr_tracker_to_iron_DIS->Clone();
            simEventRate_DIS_cross = normalize(simEventRate_DIS_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_DIS_cross->Clone()->Write(Form("simEventRate_DIS_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // Other
            auto simEventRate_Other_cross = h_eff_corr_tracker_to_iron_Other->Clone();
            simEventRate_Other_cross = normalize(simEventRate_Other_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_Other_cross->Clone()->Write(Form("simEventRate_Other_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // 2p2h
            auto simEventRate_2p2h_cross = h_eff_corr_tracker_to_iron_2p2h->Clone();
            simEventRate_2p2h_cross = normalize(simEventRate_2p2h_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_2p2h_cross->Clone()->Write(Form("simEventRate_2p2h_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));

          }
          if(mat=="lead"){
            // QE
            auto simEventRate_QE_cross = h_eff_corr_tracker_to_lead_QE->Clone();
            simEventRate_QE_cross = normalize(simEventRate_QE_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_QE_cross->Clone()->Write(Form("simEventRate_QE_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // RES
            auto simEventRate_RES_cross = h_eff_corr_tracker_to_lead_RES->Clone();
            simEventRate_RES_cross = normalize(simEventRate_RES_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_RES_cross->Clone()->Write(Form("simEventRate_RES_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // DIS
            auto simEventRate_DIS_cross = h_eff_corr_tracker_to_lead_DIS->Clone();
            simEventRate_DIS_cross = normalize(simEventRate_DIS_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_DIS_cross->Clone()->Write(Form("simEventRate_DIS_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // Other
            auto simEventRate_Other_cross = h_eff_corr_tracker_to_lead_Other->Clone();
            simEventRate_Other_cross = normalize(simEventRate_Other_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_Other_cross->Clone()->Write(Form("simEventRate_Other_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));
            // 2p2h
            auto simEventRate_2p2h_cross = h_eff_corr_tracker_to_lead_2p2h->Clone();
            simEventRate_2p2h_cross = normalize(simEventRate_2p2h_cross , fluxIntegral, nNucleons, pot);
            fUnfold->cd();
            simEventRate_2p2h_cross->Clone()->Write(Form("simEventRate_2p2h_crossSection_%s_%s_%s",mat.c_str(), prefix.c_str(), var.c_str() ));

          }
          //auto crossSection_mc = effNum->Clone();
          //crossSection_mc = normalize(crossSection_mc, fluxIntegral, nNucleons, pot);
          //crossSection_mc->Clone()->Write(Form("crossSection_%s_%s", prefix.c_str(), var.c_str() ));
          
          std::cout << "MC differential cross-section  DONE.\n" << std::flush;
        }

        // Total CrossSection
        if (prefix == "data"){
          if (var=="Enu"){
            auto efficiencyCorrected = h_eff_corr_tracker_to_material[iter]->Clone();
            auto crossSection_total = normalizeTotal( efficiencyCorrected, fluxRebinned, nNucleons, pot);
            //Plot(*crossSection_mc, "crossSection_mc", prefix);
            fUnfold->cd();
            crossSection_total->Clone()->Write(Form("crossSection_total_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
            
            std::cout << "Data total cross-section  DONE.\n" << std::flush;
          }
          // Differential CrossSection
          auto crossSection = normalize(h_eff_corr_tracker_to_material[iter], fluxIntegral, nNucleons, pot);
          //Plot(*crossSection_mc, "crossSection_mc", prefix);
          fUnfold->cd();
          crossSection->Clone()->Write(Form("crossSection_%s_%s_%s", mat.c_str(), prefix.c_str(), var.c_str() ));
          
          std::cout << "Data differential cross-section  DONE.\n" << std::flush;
          //WriteCovAndCorrMatrices(var, crossSection, fUnfold, prefix);
          //std::cout << "Wrote covariance and correlation matrices.\n" << std::flush;
        } 
        iter++;
      }
    }
  }

  std::cout << "DONE" << std::endl;


}
