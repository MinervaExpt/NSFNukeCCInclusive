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
  if(!unfold.UnfoldHisto( h_unfolded, dummyCovMatrix, h_migration, h_folded, RooUnfold::kBayes, num_iter, true, true ))
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

std::vector<int> whichTargetIDs(int targetZ)
{
  std::vector<int> targetIDs;

 if( targetZ == 99){
    targetIDs.push_back(99);
  }
  
  else{
    if(targetZ==26){
      targetIDs.push_back(2);
      targetIDs.push_back(3);
      targetIDs.push_back(5);
    }
    if(targetZ==82){
      targetIDs.push_back(2);
      targetIDs.push_back(3);
      targetIDs.push_back(4);
      targetIDs.push_back(5);
    }
    if(targetZ==6){
      targetIDs.push_back(3);
    }
  }
  return targetIDs;

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
      "\t-Material_atomic_number\t =\t Atomic number of material, e.g. 6, 26, 82, or tracker 99  \n" <<
      "\t-Playlist\t \t =\t eg. minervame1A \n" <<
      "\t-Target_number\t \t = \t Number of target, e.g. 1-5, or tracker 99  \n" <<
      "\t-Tracker_flux\t =\t Use FRW instead of target specific flux \n"<< std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    return 0;
  }

  const string playlist= argv[3];
  const std::string plist_string(playlist);

  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF("minervame6A");
  NukeCC_Binning  *binsDef = new NukeCC_Binning();
  HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist("minervame6A");

  // to iterate over variables
  std::vector<string> vars;
  vars.push_back("Enu");
  vars.push_back("x");
  vars.push_back("pTmu1D");
  vars.push_back("pZmu1D");
  vars.push_back("ThetamuDeg");

  // to iterate over mc and data
  std::vector<string> prefixes;
  prefixes.push_back("mc");
  prefixes.push_back("data");

  std::vector<int> targetIDs;

  string outdir = argv[1];
  int targetZ = atoi(argv[2]);
  int targetID; 

  if (argc == 5){
    targetID = atoi(argv[4]);
    targetIDs.push_back(targetID);
  }
  // To combine different targets
  else{
    targetIDs = whichTargetIDs(targetZ);
    if (targetIDs.size() == 1)
    {
      targetID = targetIDs[0];
    }
    else{
      std::stringstream ss_result;
      std::copy(targetIDs.begin(), targetIDs.end(), std::ostream_iterator<int>(ss_result, ""));
      string s_result = ss_result.str();
      targetID = std::stoi(s_result);
    }

  }

  bool trackerFlux;
  if (argc == 6){
    trackerFlux = argv[5];
  }
  else{
    trackerFlux= false;
  }

  // output file
  TFile *fUnfold;
  if (trackerFlux == true){
    fUnfold = new TFile( Form("%s/CrossSection_t%d_z%02d_%s.root", outdir.c_str(), targetID, targetZ, plist_string.c_str() ), "recreate" );
  }
  else{
    fUnfold = new TFile( Form("%s/CrossSection_Daisy_t%d_z%02d_%s.root", outdir.c_str(), targetID, targetZ, plist_string.c_str() ), "recreate" );
  }

  TString eventLoop, migr, efficiency;
  double mcpot, datapot;

  
  //=========================================
  // Targets combining Loop
  //=========================================
  for(int t = 0; t < targetIDs.size(); t++){
    
    // Read in background subtracted event rate, Migration and Efficiency
    if(RunCodeWithSystematics){
          eventLoop = Form("%s/BackgroundSubtracted/BkgSubtracted_EventSelection_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetIDs[t], targetZ);
          migr = Form("%s/Migration/Migration_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetIDs[t], targetZ);
          efficiency = Form("%s/Efficiency/Efficiency_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetIDs[t], targetZ);
      }
    else{
        eventLoop = Form("%s/BackgroundSubtracted/BkgSubtracted_EventSelection_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetIDs[t], targetZ);
        migr = Form("%s/Migration/Migration_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetIDs[t], targetZ);
        efficiency = Form("%s/Efficiency/Efficiency_%s_t%d_z%02d_sys.root", outdir.c_str(), plist_string.c_str(), targetIDs[t], targetZ);
    }

    TFile *fEventLoop = new TFile( eventLoop,"read" );
    TFile *fMigration = new TFile( migr,"read" );
    TFile *fEfficiency = new TFile( efficiency,"read" );
    
    // read in the POT information
    if (t == 0){
      TParameter<double> *mcPOT = (TParameter<double>*)fMigration->Get("MCPOT");
      TParameter<double> *dataPOT = (TParameter<double>*)fMigration->Get("DataPOT");
      mcpot = mcPOT->GetVal();
      datapot = dataPOT->GetVal();
    }

    fUnfold->cd();

    for(const auto&var: vars){

      // Var specific ingredients
      auto migration  = dynamic_cast<MnvH2D*> (fMigration->Get( Form("selected_Migration_%s", var.c_str())));
      auto effNum = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_mc_%s", var.c_str())));
      auto effDenom = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_%s", var.c_str())));
      auto effDenom_QE = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_QE_%s", var.c_str())));
      auto effDenom_RES = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_RES_%s", var.c_str())));
      auto effDenom_DIS = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_DIS_%s", var.c_str())));
      auto effDenom_Other = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_Other_%s", var.c_str())));
      auto effDenom_2p2h = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_2p2h_%s", var.c_str())));
      auto simEventRate = effDenom->Clone();
      auto simEventRate_QE = effDenom_QE->Clone();
      auto simEventRate_RES = effDenom_RES->Clone();
      auto simEventRate_DIS = effDenom_DIS->Clone();
      auto simEventRate_Other = effDenom_Other->Clone();
      auto simEventRate_2p2h = effDenom_2p2h->Clone();
      
      fUnfold->cd();
      simEventRate->Clone()->Write(Form("simEventRate_%d_%s", targetIDs[t],var.c_str() ));
      simEventRate_QE->Clone()->Write(Form("simEventRate_QE_%d_%s", targetIDs[t],var.c_str() ));
      simEventRate_RES->Clone()->Write(Form("simEventRate_RES_%d_%s", targetIDs[t],var.c_str() ));
      simEventRate_DIS->Clone()->Write(Form("simEventRate_DIS_%d_%s", targetIDs[t],var.c_str() ));
      simEventRate_Other->Clone()->Write(Form("simEventRate_Other_%d_%s", targetIDs[t],var.c_str() ));
      simEventRate_2p2h->Clone()->Write(Form("simEventRate_2p2h_%d_%s", targetIDs[t],var.c_str() ));

      effNum->Clone()->Write(Form("efficiency_numerator_%d_%s", targetIDs[t], var.c_str() ));
      // EFFICIENCY
      effNum->Divide(effNum, effDenom, 1.0, 1.0, "B"); //Only the 2 parameter version of MnvH1D::Divide() //handles systematics correctly.
      //Plot(*effNum, "efficiency", prefix);
      effNum->Clone()->Write(Form("efficiency_%d_%s", targetIDs[t], var.c_str() )); 

      //  Data or MC
      for(const auto& prefix: prefixes){
        bool isMC = false;
        if (prefix=="mc"){
          isMC = true;
        }

        std::cout << "Reading in: " << prefix<<".\n" << std::flush;


        auto bkgSubtracted = dynamic_cast<MnvH1D*>(fEventLoop->Get( Form("h_bkg_subtracted_%s_%s", prefix.c_str(), var.c_str())));

        cout<<"I am here"<<endl;
        //d'Aogstini unfolding
        std::cout << "===================================== \n" << std::flush;
        std::cout << "============= UNFOLDING =============\n" << std::flush;
        int nIterations=3;    
        std::cout << "Number of iterations: " << nIterations << "\n" << std::flush;

        // UNFOLDING
        auto unfolded = UnfoldHist(bkgSubtracted, migration, nIterations);
        std::cout<< "Hello"<<std::endl;
        if(!unfolded) throw std::runtime_error(std::string("Failed to unfold ") + bkgSubtracted->GetName() + " using " + migration->GetName());
        //Plot(*unfolded_mc, "unfolded", prefix);
        unfolded->Clone()->Write(Form("unfolded_%d_%s_%s", targetIDs[t], prefix.c_str(),var.c_str() ));

        std::cout << "Survived writing the unfolded histogram.\n" << std::flush;

 // CROSS-SECTION EXTRACTION

        // Unfolded efficiency corrected
        unfolded->Divide(unfolded, effNum);
        //Plot(*unfolded_data, "efficiencyCorrected", prefix);
        unfolded->Clone()->Write(Form("unfolded_effCorrected_%d_%s_%s", targetIDs[t], prefix.c_str(),var.c_str() )); 
        std::cout << "Efficiency corrected.\n" << std::flush;

        if (t == targetIDs.size() -1){
          std::cout << "All files with given targetZ unfolded and efficiency corrected.\n" << std::flush;
        }
      }
    }
  }

  for(const auto&var: vars){

    std::cout << "Reading in effNum, simEventRate from individual targets to combine.\n" << std::flush;
    auto effNum_combined = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("efficiency_numerator_%d_%s", targetIDs[0],var.c_str() )));
    auto simEventRate_combined = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_%d_%s", targetIDs[0],var.c_str() )));
    auto simEventRate_QE_combined = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_QE_%d_%s", targetIDs[0],var.c_str() )));
    auto simEventRate_RES_combined = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_RES_%d_%s", targetIDs[0],var.c_str() )));
    auto simEventRate_DIS_combined = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_DIS_%d_%s", targetIDs[0],var.c_str() )));
    auto simEventRate_Other_combined = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_Other_%d_%s", targetIDs[0],var.c_str() )));
    auto simEventRate_2p2h_combined = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_2p2h_%d_%s", targetIDs[0],var.c_str() )));
      
    for(int t = 1; t < targetIDs.size(); t++){
      auto effNum = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("efficiency_numerator_%d_%s", targetIDs[t],var.c_str() )));
      auto simEventRate = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_%d_%s", targetIDs[t],var.c_str() )));
      auto simEventRate_QE = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_QE_%d_%s", targetIDs[t],var.c_str() )));
      auto simEventRate_RES = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_RES_%d_%s", targetIDs[t],var.c_str() )));
      auto simEventRate_DIS = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_DIS_%d_%s", targetIDs[t],var.c_str() )));
      auto simEventRate_Other = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_Other_%d_%s", targetIDs[t],var.c_str() )));
      auto simEventRate_2p2h = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("simEventRate_2p2h_%d_%s", targetIDs[t],var.c_str() )));

      effNum_combined->Add(effNum);
      simEventRate_combined->Add(simEventRate);
      simEventRate_QE_combined->Add(simEventRate_QE);
      simEventRate_RES_combined->Add(simEventRate_RES);
      simEventRate_DIS_combined->Add(simEventRate_DIS);
      simEventRate_Other_combined->Add(simEventRate_Other);
      simEventRate_2p2h_combined->Add(simEventRate_2p2h);
    }

    effNum_combined->Clone()->Write(Form("total_efficiency_numerator_%s",var.c_str() )); 
    simEventRate_combined->Clone()->Write(Form("total_simEventRate_%s",var.c_str() )); 
    simEventRate_QE_combined->Clone()->Write(Form("total_simEventRate_QE_%s",var.c_str() ));
    simEventRate_RES_combined->Clone()->Write(Form("total_simEventRate_RES_%s",var.c_str() )); 
    simEventRate_DIS_combined->Clone()->Write(Form("total_simEventRate_DIS_%s",var.c_str() )); 
    simEventRate_DIS_combined->Clone()->Write(Form("total_simEventRate_Other_%s",var.c_str() )); 
    simEventRate_2p2h_combined->Clone()->Write(Form("total_simEventRate_2p2h_%s",var.c_str() ));  

    for(const auto& prefix: prefixes){
      bool isMC = false;
      if (prefix=="mc"){
        isMC = true;
      }

      std::cout << "Reading in unfolded mc and data from individual targets to combine.\n" << std::flush;
      auto unfolded_combined = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("unfolded_effCorrected_%d_%s_%s", targetIDs[0], prefix.c_str(),var.c_str() )));
      
      for(int t = 1; t < targetIDs.size(); t++){
        auto unfolded = dynamic_cast<MnvH1D*>(fUnfold->Get( Form("unfolded_effCorrected_%d_%s_%s", targetIDs[t], prefix.c_str(),var.c_str() )));
        
        unfolded_combined->Add(unfolded);
      }

      unfolded_combined->Clone()->Write(Form("total_unfolded_effCorrected_%s_%s", prefix.c_str(),var.c_str() )); 


      // FLUX INFO
      // ---------------------------------------------------------------------
      // Flux reweighter information
      // ---------------------------------------------------------------------
      int n_flux_universes =  PlotUtils::MinervaUniverse::GetNFluxUniverses();
      int nu_pdg = PlotUtils::MinervaUniverse::GetAnalysisNuPDG();
      const bool use_nue_constraint = true;
      bool useMuonCorrelations = true;
      const std::string project_dir = "targets_2345_jointNueIMD";
      double min_energy = 0;
      double max_energy = 120;

      auto& frw = PlotUtils::flux_reweighter("minervame6A", nu_pdg, use_nue_constraint, n_flux_universes);
      auto& frw2 = PlotUtils::flux_reweighter("minervame6A", nu_pdg, use_nue_constraint, n_flux_universes);

     
      MnvH1D* flux;
      MnvH1D* fluxIntegral;
      MnvH1D* fluxRebinned;

      string material; 
      
      if(targetZ == 6) material = "carbon";
      else if(targetZ == 26) material = "iron";
      else if(targetZ == 82) material = "lead";
      else material = "scintillator";

      if(trackerFlux == true){
        flux = frw.GetFluxReweighted(-14);
        fluxIntegral = frw.GetIntegratedFluxReweighted(-14, unfolded_combined, 0, 120, false);
        fluxRebinned = frw2.GetRebinnedFluxReweighted(-14, unfolded_combined);
        std::cout<< "Using FRW flux!" << std::endl;
      }

      else{
        fluxIntegral = frw.GetIntegratedTargetFlux(nu_pdg, material, unfolded_combined, min_energy, max_energy, project_dir);
        flux = frw2.GetTargetFluxMnvH1D(nu_pdg, material, project_dir);
        fluxRebinned = frw2.GetRebinnedFluxReweighted_FromInputFlux(flux, unfolded_combined); // issue here
        std::cout<< "Using target generated flux!" << std::endl;
      }

      if (prefix == "mc"){
        fUnfold->cd();
        flux->Clone()->Write(Form("flux_%s", var.c_str()));
        fluxRebinned->Clone()->Write(Form("fluxRebinned_%s", var.c_str()));
      }

      std::cout << "FluxReweighter information generated.\n" << std::flush;

      PlotUtils::TargetUtils targetInfo;
      double nNucleons = targetInfo.GetPassiveTargetNNucleons( targetID, targetZ, isMC );
      //double nNucleons = GetTotalScatteringCenters(targetZ, isMC);

      std::cout << prefix + " number of nucleons" << std::endl;
      std::cout << nNucleons << std::endl;

      double pot = datapot;
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
          auto simEventRate_cross_total = simEventRate_combined->Clone();
          simEventRate_cross_total = normalizeTotal(simEventRate_cross_total, fluxRebinned, nNucleons, pot);
          fUnfold->cd();
          simEventRate_cross_total->Clone()->Write(Form("simEventRate_crossSection_total_%s_%s", prefix.c_str(), var.c_str() ));

          // QE
          auto simEventRate_QE_cross_total = simEventRate_QE_combined->Clone();
          simEventRate_QE_cross_total = normalizeTotal(simEventRate_QE_cross_total, fluxRebinned, nNucleons, pot);
          fUnfold->cd();
          simEventRate_QE_cross_total->Clone()->Write(Form("simEventRate_QE_crossSection_total_%s_%s", prefix.c_str(), var.c_str() ));

          //RES
          auto simEventRate_RES_cross_total = simEventRate_RES_combined->Clone();
          simEventRate_RES_cross_total = normalizeTotal(simEventRate_RES_cross_total, fluxRebinned, nNucleons, pot);
          fUnfold->cd();
          simEventRate_RES_cross_total->Clone()->Write(Form("simEventRate_RES_crossSection_total_%s_%s", prefix.c_str(), var.c_str() ));
          
          //DIS
          auto simEventRate_DIS_cross_total = simEventRate_DIS_combined->Clone();
          simEventRate_DIS_cross_total = normalizeTotal(simEventRate_DIS_cross_total, fluxRebinned, nNucleons, pot);
          fUnfold->cd();
          simEventRate_DIS_cross_total->Clone()->Write(Form("simEventRate_DIS_crossSection_total_%s_%s", prefix.c_str(), var.c_str() ));
          
          //Other
          auto simEventRate_Other_cross_total = simEventRate_Other_combined->Clone();
          simEventRate_Other_cross_total = normalizeTotal(simEventRate_Other_cross_total, fluxRebinned, nNucleons, pot);
          fUnfold->cd();
          simEventRate_Other_cross_total->Clone()->Write(Form("simEventRate_Other_crossSection_total_%s_%s", prefix.c_str(), var.c_str() ));

          //2p2h
          auto simEventRate_2p2h_cross_total = simEventRate_2p2h_combined->Clone();
          simEventRate_2p2h_cross_total = normalizeTotal(simEventRate_2p2h_cross_total, fluxRebinned, nNucleons, pot);
          fUnfold->cd();
          simEventRate_2p2h_cross_total->Clone()->Write(Form("simEventRate_2p2h_crossSection_total_%s_%s", prefix.c_str(), var.c_str() ));
        
          std::cout << "MC total cross-section DONE.\n" << std::flush;

        }
        auto simEventRate_cross_dif = simEventRate_combined->Clone();
        simEventRate_cross_dif = normalize(simEventRate_cross_dif, fluxIntegral, nNucleons, pot);
        fUnfold->cd();
        simEventRate_cross_dif->Clone()->Write(Form("simEventRate_crossSection_%s_%s", prefix.c_str(), var.c_str() ));

        //QE
        auto simEventRate_QE_cross_dif = simEventRate_QE_combined->Clone();
        simEventRate_QE_cross_dif = normalize(simEventRate_QE_cross_dif, fluxIntegral, nNucleons, pot);
        fUnfold->cd();
        simEventRate_QE_cross_dif->Clone()->Write(Form("simEventRate_QE_crossSection_%s_%s", prefix.c_str(), var.c_str() ));

        //RES
        auto simEventRate_RES_cross_dif = simEventRate_RES_combined->Clone();
        simEventRate_RES_cross_dif = normalize(simEventRate_RES_cross_dif, fluxIntegral, nNucleons, pot);
        fUnfold->cd();
        simEventRate_RES_cross_dif->Clone()->Write(Form("simEventRate_RES_crossSection_%s_%s", prefix.c_str(), var.c_str() ));

        //DIS
        auto simEventRate_DIS_cross_dif = simEventRate_DIS_combined->Clone();
        simEventRate_DIS_cross_dif = normalize(simEventRate_DIS_cross_dif, fluxIntegral, nNucleons, pot);
        fUnfold->cd();
        simEventRate_DIS_cross_dif->Clone()->Write(Form("simEventRate_DIS_crossSection_%s_%s", prefix.c_str(), var.c_str() ));

        //Other
        auto simEventRate_Other_cross_dif = simEventRate_Other_combined->Clone();
        simEventRate_Other_cross_dif = normalize(simEventRate_Other_cross_dif, fluxIntegral, nNucleons, pot);
        fUnfold->cd();
        simEventRate_Other_cross_dif->Clone()->Write(Form("simEventRate_Other_crossSection_%s_%s", prefix.c_str(), var.c_str() ));

        //2p2h
        auto simEventRate_2p2h_cross_dif = simEventRate_2p2h_combined->Clone();
        simEventRate_2p2h_cross_dif = normalize(simEventRate_2p2h_cross_dif, fluxIntegral, nNucleons, pot);
        fUnfold->cd();
        simEventRate_2p2h_cross_dif->Clone()->Write(Form("simEventRate_2p2h_crossSection_%s_%s", prefix.c_str(), var.c_str() ));

        std::cout << "MC differential cross-section DONE.\n" << std::flush;
      }

      std::cout << prefix + " POT" << std::endl;
      std::cout << pot << std::endl;

      if (prefix=="data"){
        // Total CrossSection
        if (var=="Enu"){
          auto efficiencyCorrected = unfolded_combined->Clone();
          auto crossSection_total = normalizeTotal(efficiencyCorrected, fluxRebinned, nNucleons, pot);
          //Plot(*crossSection_mc, "crossSection_mc", prefix);
          fUnfold->cd();
          crossSection_total->Clone()->Write(Form("crossSection_total_%s_%s", prefix.c_str(), var.c_str() ));
          std::cout << "Total cross-section DONE.\n" << std::flush;
        }
        // Differential CrossSection
        auto crossSection = normalize(unfolded_combined, fluxIntegral, nNucleons, pot);
        //Plot(*crossSection_mc, "crossSection_mc", prefix);
        fUnfold->cd();
        crossSection->Clone()->Write(Form("crossSection_%s_%s", prefix.c_str(), var.c_str() ));

        WriteCovAndCorrMatrices(var, unfolded_combined, fUnfold, prefix);
        std::cout << "Wrote covariance and correlation cross-section matrices.\n" << std::flush;
      }
    }
  }

  fUnfold->cd();
  // write the POT information into the output file
  auto dataPOTOut = new TParameter<double>("DataPOT", datapot);
  auto mcPOTOut = new TParameter<double>("MCPOT", mcpot);
  dataPOTOut->Write();
  mcPOTOut->Write(); 

  double dataMCScale = datapot/mcpot;
  cout<<"Data POT = "<<datapot<<" MC POT = "<<mcpot<<"   Data/MC POT = "<<dataMCScale<<endl;

  std::cout << "DONE" << std::endl;


}
