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

// Differential cross-section
PlotUtils::MnvH1D* normalizeTotal(PlotUtils::MnvH1D* efficiencyCorrected, PlotUtils::MnvH1D* fluxRebinned, const double nNucleons, const double POT)
{
  efficiencyCorrected->Divide(efficiencyCorrected, fluxRebinned);

  efficiencyCorrected->Scale(1./nNucleons/POT);
  efficiencyCorrected->Scale(1.e4); //Flux histogram is in m^-2, but convention is to report cm^2

  return efficiencyCorrected;
}

 
//int UnfoldIterations( const std::string& var );

int main(int argc, char * argv[]){
  //ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if(argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------\
-----"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n"<<
       "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<< std::endl;
      //"\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
       //"\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" <<
      //"\t-Playlist\t \t =\t eg. minervame1A"<< std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    return 0;
  }

  string outdir=argv[1];
  int targetID = 99;//atoi(argv[2]);
  int targetZ = 99;//atoi(argv[3]);
  //const string playlist= argv[4];
  const std::string plist_string("minervame6A");

  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);

  TString eventLoop, migr, efficiency;
  if(RunCodeWithSystematics){
        eventLoop = Form("%s/BackgroundSubtracted_EventSelection_bkgSubtract_tracker.root", outdir.c_str());
        //hname1 = Form("%s/Hists_EventSelectionTracker_ML_ME6A_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
        migr = Form("%s/Hists_MigrationTracker_ML_ME6A_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
        efficiency = Form("%s/Hists_EfficiencyTracker_ML_ME6A_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
    }
    else{
        eventLoop = Form("%s/BackgroundSubtracted_EventSelection_bkgSubtract_tracker.root", outdir.c_str());
        //hname1 = Form("%s/Hists_EventSelectionTracker_ML_ME6A_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
        migr = Form("%s/Hists_MigrationTracker_ML_ME6A_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
        efficiency = Form("%s/Hists_EfficiencyTracker_ML_ME6A_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ);
    }

    TFile *fEventLoop = new TFile( eventLoop,"read" );
    TFile *fMigration = new TFile( migr,"read" );
    TFile *fEfficiency = new TFile( efficiency,"read" );

    NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
    NukeCC_Binning  *binsDef = new NukeCC_Binning();
    HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);

    //vector<string> vars;
    //vars.clear();
    //vars.push_back("Enu");
    //vars.push_back("x");

    string var = "Enu";

    auto bkgSubtracted_mc =  dynamic_cast<MnvH1D*>(fEventLoop->Get( Form("h_background_subtracted_mc_%s", var.c_str())));
    auto bkgSubtracted_data =  dynamic_cast<MnvH1D*>(fEventLoop->Get( Form("h_background_subtracted_data_%s", var.c_str()))); 
    auto migration  = dynamic_cast<MnvH2D*> (fMigration->Get( Form("selected_Migration_%s", var.c_str())));
    
    auto effNum = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_mc_%s", var.c_str())));
    auto effDenom = dynamic_cast<MnvH1D*>(fEfficiency->Get( Form("h_truth_%s", var.c_str())));

    cout<<"I am here"<<endl;
    
    TFile *fUnfold = new TFile( Form("%s/CrossSection_t%d_z%d_%s.root", outdir.c_str(), targetID, targetZ, plist_string.c_str() ), "recreate" );
    
    TParameter<double> *mcPOT = (TParameter<double>*)fMigration->Get("MCPOT");
    TParameter<double> *dataPOT = (TParameter<double>*)fMigration->Get("DataPOT");
    double mcpot = mcPOT->GetVal();
    double datapot = dataPOT->GetVal();

    double dataMCScale = datapot/mcpot;
    cout<<"Data POT = "<<datapot<<" MC POT = "<<mcpot<<"   Data/MC POT = "<<dataMCScale<<endl;

    //d'Aogstini unfolding
    int nIterations=15;
    string prefix;
    prefix = "mc";

    // MC
    auto unfolded_mc = UnfoldHist(bkgSubtracted_mc, migration, nIterations);
    if(!unfolded_mc) throw std::runtime_error(std::string("Failed to unfold ") + bkgSubtracted_mc->GetName() + " using " + migration->GetName());
    //Plot(*unfolded_mc, "unfolded", prefix);
    unfolded_mc->Clone()->Write(Form("unfolded_mc_%s", var.c_str() ));
     // histogram unfolded data stores the diagonals
    // "unfoldinCov" contains the covariant values
    // in unfoldingCov, the diagonals are zero, but there should be off-diagonal content
    
    prefix.clear();
    prefix = "data";

    // Data
    auto unfolded_data = UnfoldHist(bkgSubtracted_data, migration, nIterations);
    if(!unfolded_data) throw std::runtime_error(std::string("Failed to unfold ") + bkgSubtracted_data->GetName() + " using " + migration->GetName());
    //Plot(*unfolded_data, "unfolded", prefix);
    unfolded_data->Clone()->Write(Form("unfolded_data_%s", var.c_str() )); 
    // Get unfolding-only correlation matrix
    // False in GetStatErrorMatrix is fractional
    // False in GetSysCorrelationMatrix means that the covariance is not area-normalised

    //unfolded_data->GetStatErrorMatrix(false) + unfolded_data->GetSysCorrelationMatrix("unfoldingCov", false);
    //unfolded_data->Clone()->Write(Form("unfolded_data_covMatrix_%s", var.c_str() ));

    std::cout << "Survived writing the unfolded histogram.\n" << std::flush;

    fUnfold->cd();
    auto dataPOTOut = new TParameter<double>("DataPOT", datapot);
    auto mcPOTOut = new TParameter<double>("MCPOT", mcpot);
    dataPOTOut->Write();
    mcPOTOut->Write(); 

    // Cross-section extraction

    // Efficiency
    effNum->Divide(effNum, effDenom, 1.0, 1.0, "B"); //Only the 2 parameter version of MnvH1D::Divide() //handles systematics correctly.
    //Plot(*effNum, "efficiency", prefix);
    effNum->Clone()->Write(Form("efficiency_%s", var.c_str() )); 

    // Unfolded MC efficiency corrected
    prefix.clear();
    prefix = "data";
    unfolded_mc->Divide(unfolded_mc, effNum);
    //Plot(*unfolded_data, "efficiencyCorrected", prefix);
    unfolded_mc->Clone()->Write(Form("unfolded_effCorrected_mc_%s", var.c_str() )); 
    
    prefix.clear();
    prefix = "data";

    // Unfolded data efficiency corrected
    unfolded_data->Divide(unfolded_data, effNum);
    //Plot(*unfolded_data, "efficiencyCorrected", prefix);
    unfolded_data->Clone()->Write(Form("unfolded_effCorrected_data_%s", var.c_str() ));

    auto& frw = PlotUtils::flux_reweighter(plist_string,-14, true, 100);
    auto flux = frw.GetFluxReweighted(-14);
    flux->Clone()->Write(Form("flux_%s", var.c_str()));
    auto fluxIntegral = frw.GetIntegratedFluxReweighted(-14, unfolded_data, 2, 50, false);

    auto fluxRebinned = frw.GetRebinnedFluxReweighted(-14, unfolded_mc);

    PlotUtils::TargetUtils targetInfo;
    double nNucleons = targetInfo.GetTrackerNNucleons(5980, 8422, false, 850);
    //double TargetUtils::GetTrackerNNucleons( double minZ, double maxZ, bool isMC, double apothem /* = 850. */ ) const

    // Total CrossSection
    prefix.clear();
    prefix = "MC";
    auto efficiencyCorrected_mc = unfolded_mc->Clone();
    auto crossSection_total_mc = normalizeTotal(efficiencyCorrected_mc, fluxRebinned, nNucleons, datapot);
    //Plot(*crossSection_mc, "crossSection_mc", prefix);
    crossSection_total_mc->Clone()->Write(Form("crossSection_total_mc_%s", var.c_str() ));

    // Total CrossSection
    prefix.clear();
    prefix = "data";
    auto efficiencyCorrected_data = unfolded_data->Clone();
    auto crossSection_total_data = normalizeTotal(efficiencyCorrected_data, fluxRebinned, nNucleons, datapot);
    //Plot(*crossSection_data, "crossSection_data", prefix);
    crossSection_total_data->Clone()->Write(Form("crossSection_total_data_%s", var.c_str() ));
    
    
    // Differential CrossSection
    prefix.clear();
    prefix = "MC";
    auto crossSection_mc = normalize(unfolded_mc, fluxIntegral, nNucleons, datapot);
    //Plot(*crossSection_mc, "crossSection_mc", prefix);
    crossSection_mc->Clone()->Write(Form("crossSection_mc_%s", var.c_str() ));


    // Data differential cross-section
    prefix.clear();
    prefix = "data";

    auto crossSection_data = normalize(unfolded_data, fluxIntegral, nNucleons, datapot);
    //Plot(*crossSection_data, "crossSection_data", prefix);
    crossSection_data->Clone()->Write(Form("crossSection_data_%s", var.c_str() ));

    std::cout << "DONE" << std::endl;


}






