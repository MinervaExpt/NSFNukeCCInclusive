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
      "\t-Material_atomic_number\t =\t Atomic number of material, e.g. 6, 26, 82, or tracker 99  \n" << std::endl;
      //"\t-Playlist\t \t =\t eg. minervame1A"<< std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    return 0;
  }

  string outdir = argv[1];
  int targetID = atoi(argv[2]);
  int targetZ = atoi(argv[3]);

  const std::string plist_string("minervame6A");

  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);

  TString eventLoop, migr, efficiency;

  // Read in background subtracted event rate, Migration and Efficiency
  if(RunCodeWithSystematics){
        eventLoop = Form("%s/Hists_BkgSubtracted_EventSelection_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
        migr = Form("%s/Hists_Migration_ML_ME6A_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
        efficiency = Form("%s/Hists_Efficiency_ML_ME6A_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
    }
  else{
      eventLoop = Form("%s/Hists_BkgSubtracted_EventSelection_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
      migr = Form("%s/Hists_Migration_ML_ME6A_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
      efficiency = Form("%s/Hists_Efficiency_ML_ME6A_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
  }

  TFile *fEventLoop = new TFile( eventLoop,"read" );
  TFile *fMigration = new TFile( migr,"read" );
  TFile *fEfficiency = new TFile( efficiency,"read" );

  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
  NukeCC_Binning  *binsDef = new NukeCC_Binning();
  HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);

  // to iterate over variables
  std::vector<string> vars;
  vars.push_back("Enu");
  vars.push_back("x");
  
  //string var = "Enu";

  // to iterate over mc and data
  std::vector<string> prefixes;
  prefixes.push_back("mc");
  prefixes.push_back("data");

  // output file
  TFile *fUnfold = new TFile( Form("%s/CrossSection_t%d_z%d_%s_NEW.root", outdir.c_str(), targetID, targetZ, plist_string.c_str() ), "recreate" );
  
  // read in the POT information
  TParameter<double> *mcPOT = (TParameter<double>*)fMigration->Get("MCPOT");
  TParameter<double> *dataPOT = (TParameter<double>*)fMigration->Get("DataPOT");
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

    // Var specific ingredients
    auto migration  = dynamic_cast<MnvH2D*> (fMigration->Get( Form("selected_Migration_%s", var.c_str())));
    auto effNum = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_mc_%s", var.c_str())));
    auto effDenom = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_%s", var.c_str())));
    auto simEventRate = effDenom->Clone();

    effNum->Clone()->Write(Form("efficiency_numerator_%s", var.c_str() ));
    // EFFICIENCY
    auto efficiency = effNum->Clone();
    efficiency->Divide(efficiency, effDenom, 1.0, 1.0, "B"); //Only the 2 parameter version of MnvH1D::Divide() //handles systematics correctly.
    //Plot(*effNum, "efficiency", prefix);
    efficiency->Clone()->Write(Form("efficiency_%s", var.c_str() )); 

    simEventRate->Clone()->Write(Form("simEventRate_%s", var.c_str() ));

    // Var and prefix specific ingredients
    for(const auto& prefix: prefixes){

      auto bkgSubtracted =  dynamic_cast<MnvH1D*>(fEventLoop->Get( Form("h_bkg_subtracted_%s_%s", prefix.c_str(), var.c_str())));

      cout<<"I am here"<<endl;
      //d'Aogstini unfolding
      int nIterations=2;  

      // UNFOLDING
      auto unfolded = UnfoldHist(bkgSubtracted, migration, nIterations);
      if(!unfolded) throw std::runtime_error(std::string("Failed to unfold ") + bkgSubtracted->GetName() + " using " + migration->GetName());
      //Plot(*unfolded_mc, "unfolded", prefix);
      unfolded->Clone()->Write(Form("unfolded_%s_%s", prefix.c_str(),var.c_str() ));
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
      unfolded->Clone()->Write(Form("unfolded_effCorrected_%s_%s", prefix.c_str(),var.c_str() )); 
      std::cout << "Efficiency corrected.\n" << std::flush;

      // FLUX INFO
      auto& frw = PlotUtils::flux_reweighter(plist_string,-14, true, 100);
      auto flux = frw.GetFluxReweighted(-14);
      auto fluxIntegral = frw.GetIntegratedFluxReweighted(-14, unfolded, 0, 120, false);
      auto& frw2 = PlotUtils::flux_reweighter(plist_string,-14, true, 100);
      auto fluxRebinned = frw2.GetRebinnedFluxReweighted(-14, unfolded);
      if (prefix == "mc"){
        flux->Clone()->Write(Form("flux_%s", var.c_str()));
        fluxRebinned->Clone()->Write(Form("fluxRebinned_%s", var.c_str()));
      }
      /*
      bool isMC = false;
      if (prefix=="mc"){
        isMC = true;
      }
      */
      PlotUtils::TargetUtils targetInfo;
      double nNucleons = targetInfo.GetPassiveTargetNNucleons( targetID, targetZ, isMC );
      //double nNucleons = GetTotalScatteringCenters(targetZ, isMC);

      std::cout << prefix + " number of nucleons" << std::endl;
      std::cout << nNucleons << std::endl;

      double pot = datapot;

      std::cout << prefix + " POT" << std::endl;
      std::cout << pot << std::endl;
      
      if (prefix=="mc"){
      // simulated event rate from efficiency denominator -> cross-section (to check with GENIE xSec)
      // MC cross-section calculated from the efficiency numerator
      // Unfolding inflates the statistical uncertainty and will give a reader the wrong impression of the size of the MC sample.
        pot = mcpot;
        if (var=="Enu"){ // calculated simulated total cross-section for Enu
          auto simEventRate_cross_total = simEventRate->Clone();
          simEventRate_cross_total = normalizeTotal(simEventRate_cross_total, fluxRebinned, nNucleons, pot);
          simEventRate_cross_total->Clone()->Write(Form("simEventRate_crossSection_total_%s_%s", prefix.c_str(), var.c_str() ));

          //auto crossSection_total_mc = effNum->Clone();
          //crossSection_total_mc= normalizeTotal(crossSection_total_mc, fluxRebinned, nNucleons, pot);
          //crossSection_total_mc->Clone()->Write(Form("crossSection_total_%s_%s", prefix.c_str(), var.c_str() ));
          //std::cout << "MC total cross-section  DONE.\n" << std::flush;

        }
        auto simEventRate_cross = simEventRate->Clone();
        simEventRate_cross = normalize(simEventRate_cross, fluxIntegral, nNucleons, pot);
        simEventRate_cross->Clone()->Write(Form("simEventRate_crossSection_%s_%s", prefix.c_str(), var.c_str() ));

        //auto crossSection_mc = effNum->Clone();
        //crossSection_mc = normalize(crossSection_mc, fluxIntegral, nNucleons, pot);
        //crossSection_mc->Clone()->Write(Form("crossSection_%s_%s", prefix.c_str(), var.c_str() ));
        
        std::cout << "MC differential cross-section  DONE.\n" << std::flush;
      }

      // Total CrossSection
      if (prefix == "data"){
        if (var=="Enu"){
          auto efficiencyCorrected = unfolded->Clone();
          auto crossSection_total = normalizeTotal(efficiencyCorrected, fluxRebinned, nNucleons, pot);
          //Plot(*crossSection_mc, "crossSection_mc", prefix);
          crossSection_total->Clone()->Write(Form("crossSection_total_%s_%s", prefix.c_str(), var.c_str() ));
          
          std::cout << "Data total cross-section  DONE.\n" << std::flush;
        }
        // Differential CrossSection
        auto crossSection = normalize(unfolded, fluxIntegral, nNucleons, pot);
        //Plot(*crossSection_mc, "crossSection_mc", prefix);
        crossSection->Clone()->Write(Form("crossSection_%s_%s", prefix.c_str(), var.c_str() ));
        
        std::cout << "Data differential cross-section  DONE.\n" << std::flush;
        //WriteCovAndCorrMatrices(var, crossSection, fUnfold, prefix);
        //std::cout << "Wrote covariance and correlation matrices.\n" << std::flush;
      } 
    }
  }

  std::cout << "DONE" << std::endl;


}
