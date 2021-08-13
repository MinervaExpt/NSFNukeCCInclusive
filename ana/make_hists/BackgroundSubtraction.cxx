//#include "include/CommonIncludes.h"
#include "../../NUKECCSRC/ana_common/include/CommonIncludes.h"
#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "../include/VariableRun.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "../../NUKECCSRC/ana_common/include/LateralSystematics.h"
#include <iostream>
#include <stdlib.h>
//#include "Cintex/Cintex.h"
#include "../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Cuts.h"
#include "TParameter.h"

#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/TargetMassSystematics.h"

// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we 
// shield it.
#ifndef __CINT__
#include "../include/plotting_functions.h"
#endif
#include "PlotUtils/MacroUtil.h" 
//using namespace globalV;
using namespace NUKECC_ANA;

int main(){
  TH1::AddDirectory(false);

  // file with background filled
  TFile f1(Form("/minerva/data/users/anezkak/ME6A_T3Fe/Hists_EventSelection_Bkg_ME6A_t3_z26_AntiNu.root"));

  MnvH1D* reco = (MnvH1D*)f1.Get("selected_mc_reco_Enu");
  MnvH1D* bkg = (MnvH1D*)f1.Get("selected_mc_reco_bkg_Enu");
  MnvH1D* data = (MnvH1D*)f1.Get("selected_data_reco_Enu");

  TFile f2(Form("/minerva/data/users/anezkak/ME6A_T3Fe/Efficiency/Hists_Efficiency_t3_z26_AntiNu.root")); // purity numerator = efficiency numerator
  MnvH1D* purityNum = (MnvH1D*)f2.Get("h_mc_Enu");

  reco->GetCVHistoWithError(true);
  bkg->GetCVHistoWithError(true);
  purityNum->GetCVHistoWithError(true);
  TFile *out1 = new TFile("BackgroundSubtracted_EventSelection_bkgSubtract.root","RECREATE");
  TFile *out2 = new TFile("BackgroundSubtracted_EventSelection_Purity.root","RECREATE");

  TParameter<double> *mcPOT = (TParameter<double>*)f1.Get("MCPOT");
  TParameter<double> *dataPOT = (TParameter<double>*)f1.Get("DataPOT");
  double mcpot = mcPOT->GetVal();
  double datapot = dataPOT->GetVal();

  double POTNormalization = datapot/mcpot;

  // background subtraction - method #1 
  MnvH1D* h_background_subtracted_mc;
  MnvH1D* h_background_subtracted_data;
  MnvH1D* h_background_mc_scale; 

  // MC subtracted
  h_background_subtracted_mc = (MnvH1D*)reco->Clone("h_background_subtracted_mc");
  h_background_subtracted_mc->Add(bkg,-1); 

  // Data subtracted
  // scale backgroud to data
  h_background_mc_scale = (MnvH1D*)bkg->Clone("h_background_mc_scale");
  h_background_mc_scale->Scale(POTNormalization);
  
  // subtract scaled bkg
  data->ClearAllErrorBands(); 
  data->AddMissingErrorBandsAndFillWithCV(*h_background_mc_scale);
  h_background_subtracted_data = (MnvH1D*)data->Clone("h_background_subtracted_data");
  h_background_subtracted_data->Add(h_background_mc_scale,-1);  

  // write all the histogramt to file
  out1->cd();
  h_background_subtracted_mc->Write();
  h_background_mc_scale->Write();
  h_background_subtracted_data->Write();  
 
  // Purity background subtraction
  out2->cd();
  MnvH1D* purity;
  purity = (MnvH1D*)purityNum->Clone("purity");
  purity->Divide(purity,reco, 1.0, 1.0, "B"); // binomial because numerator is the subset of the denominator, to make the stat error correct

  MnvH1D* h_purity_background_subtracted_mc;
  MnvH1D* h_purity_background_subtracted_data;

  data->ClearAllErrorBands(); 
  data->AddMissingErrorBandsAndFillWithCV(*reco);

  h_purity_background_subtracted_mc= (MnvH1D*)reco->Clone("h_purity_background_subtracted_mc");
  h_purity_background_subtracted_data = (MnvH1D*)data->Clone("h_purity_background_subtracted_data");

  h_purity_background_subtracted_mc->Multiply(h_purity_background_subtracted_mc, purity);
  h_purity_background_subtracted_data->Multiply(h_purity_background_subtracted_data, purity);

  h_purity_background_subtracted_mc->Write();
  h_purity_background_subtracted_data->Write();
  purity->Write();
  

  out1->Clear();
  out1->Close();
  out2->Clear();
  out2->Close();
} 


