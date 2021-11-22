#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TArrayD.h" 
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
#include "../../../NUKECCSRC/include/CVUniverse.h"
#include "../../include/Variable_plasticSB.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../../NUKECCSRC/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "../../../NUKECCSRC/include/LateralSystematics.h"
#include <iostream>
#include <stdlib.h>
//#include "Cintex/Cintex.h"
#include "../../../NUKECCSRC/include/NukeCCUtilsNSF.h"

// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we
// shield it.
#ifndef __CINT__
#include "../../include/plotting_functions.h"
#endif
using namespace std;
using namespace NUKECC_ANA;
using namespace PlotUtils;

//----Helper Functions-----------------------------------------------
void PlotBGStuff(MnvPlotter mnvPlotter, vector<MnvH1D*> histos_mc, const string var, int targetID, int targetZ, double dataMCScale, bool plotUS, bool plotDS, string playlist, string Name);
void PlotBGRatioStuff(MnvPlotter mnvPlotter, MnvH1D* histos_mc, TObject hlist, const string var, int targetID, int targetZ, double dataMCScale, double dataPOT, double mcPOT, bool plotUS, bool plotDS, string playlist, string Name);
void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name);
void GetMinAndMaxAxisRatio( TH1D* histData, TH1D* histMC, double dataMCScale, double &plotMin, double &plotMax );
//void PlotBGRatioStuff(MnvPlotter mnvPlotter, MnvH1D* histos_mc, MnvH1D* histos_mc_tuned, const string var, int targetID, int targetZ, double dataMCScale, double dataPOT, double mcPOT, bool plotUS, bool plotDS, string playlist, string Name);
//void PlotChi2Stat( MnvPlotter &mnvPlotter, /*vector<MnvH1D*> TObjArray*/ MnvH1D* histos_mc, MnvH1D* histos_mc_tuned, double dataMCScale, const string var, int targetZ);
//============================================================================================================================
// Main
int main(int argc, char * argv[]){
   //ROOT::Cintex::Cintex::Enable();
   TH1::AddDirectory(false);
	
   if(argc==1){
     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
     std::cout<<"MACROS HELP:\n\n"<<
       "\t-./runEventLoop Path_to_Input_file Target_number Material_atomic_number Playlist\n\n"<<
       "\t-Path_to_Input_file\t =\t Path to the directory where the input ROOT file is called \n"<<
       "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" <<
       "\t-Playlist\t \t =\t eg. minervame1A"<< std::endl;
     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
     return 0;
   }

   string outdir=argv[1];
   int targetID = atoi(argv[2]);
   int targetZ = atoi(argv[3]);
   const string playlist= argv[4];
   
   //ROOT::Cintex::Cintex::Enable();
   TH1::AddDirectory(false);
   
   //bool RunCodeWithSystematics = false;
   bool RunCodeWithSystematics = true;
   TString histFileName, SF_file;

   cout<<"******************************************************************************************"<<endl;
   cout<<"                   I am making plots for both before and after tuning!!                   "<<endl;
   cout<<"******************************************************************************************"<<endl;

   if(RunCodeWithSystematics){
       //histFileName = Form("%s/Hists_PlasticBackgd_with_SYS_t%d_z%02d_Nu_%s.root", outdir.c_str(), targetID, targetZ, getenv("NUKECC_TAG") ); 
       histFileName = Form("%s/Hists_PlasticBackgd_with_SYS_FullDet_q2WfromBranch_ME1A_targetsCombined_t%d_z%02d_Nu.root", outdir.c_str(), targetID, targetZ ); 
     }
   else{
       histFileName = Form("%s/Hists_PlasticBackgd_without_SYS_t%d_z%02d_Nu_%s.root", outdir.c_str(), targetID, targetZ, getenv("NUKECC_TAG"));
     } 
   
   cout<<histFileName<<endl;
  

  MnvPlotter mnvPlotter(kNukeCCStyle);
   //MnvPlotter mnvPlotter(kCCQENuStyle);

   string trueZ, mat;
   if( targetZ == 26 ){
      trueZ = "Iron";
      mat = "Fe";
    }
   if( targetZ == 82 ){
       trueZ = "Lead";
       mat = "Pb";
    }
   if( targetZ == 6 ){
        trueZ = "Carbon";
        mat = "C";
    }

   vector<string> vars;
   vars.clear();
   //vars.push_back("planeDNN");
   vars.push_back("Emu");
   vars.push_back("Ehad");
   vars.push_back("Enu");
   vars.push_back("x");
   vars.push_back("y");
   vars.push_back("planeDNN");
   //vars.push_back("W");
   //vars.push_back("Q2");

   TFile *f1 = new TFile( histFileName,"read" );

   TParameter<double> *mcPOT = (TParameter<double>*)f1->Get("MCPOT");
   TParameter<double> *dataPOT = (TParameter<double>*)f1->Get("DataPOT");
   double mcpot = mcPOT->GetVal();
   double datapot = dataPOT->GetVal(); 
   double dataMCScale = datapot/mcpot;
   cout<<"MCPOT = "<<mcpot<<"DataPOT = "<< datapot << "Scale = " << dataMCScale<<endl;

// Add mc hists to obj array, add the smallest to the bottom
   TObjArray hlistSidebandUS, hlistSidebandDS, hlistFiducial, hlistRegionUS, hlistRegionDS;

   vector< MnvH1D*> hists_US_data, hists_US_mat, hists_US_other, hists_US_regUS, hists_US_regDS;
   vector< MnvH1D*> hists_DS_data, hists_DS_mat, hists_DS_other, hists_DS_regUS, hists_DS_regDS;
   vector< MnvH1D*> histos_mc_US, histos_mc_DS;

   for( vector<string>::iterator i = vars.begin(); i != vars.end(); ++i ){
   const string& var = *i;

   f1->cd();
   hists_US_data.push_back(  (MnvH1D*)f1->Get(Form("selected_data_reco_US_%s", var.c_str())));
   hists_US_mat.push_back(   (MnvH1D*)f1->Get(Form("US_%s_%s_%s", mat.c_str(), trueZ.c_str(), var.c_str())));
   hists_US_other.push_back( (MnvH1D*)f1->Get(Form("US_other_%s_%s", trueZ.c_str(), var.c_str())));
   hists_US_regUS.push_back( (MnvH1D*)f1->Get(Form("US_regUS_%s_%s", trueZ.c_str(), var.c_str())));
   hists_US_regDS.push_back( (MnvH1D*)f1->Get(Form("US_regDS_%s_%s", trueZ.c_str(), var.c_str())));

   hists_DS_data.push_back(  (MnvH1D*)f1->Get(Form("selected_data_reco_DS_%s", var.c_str())));
   hists_DS_mat.push_back(   (MnvH1D*)f1->Get(Form("DS_%s_%s_%s", mat.c_str(), trueZ.c_str(), var.c_str())));
   hists_DS_other.push_back( (MnvH1D*)f1->Get(Form("DS_other_%s_%s", trueZ.c_str(), var.c_str())));
   hists_DS_regUS.push_back( (MnvH1D*)f1->Get(Form("DS_regUS_%s_%s", trueZ.c_str(), var.c_str())));
   hists_DS_regDS.push_back( (MnvH1D*)f1->Get(Form("DS_regDS_%s_%s", trueZ.c_str(), var.c_str())));

{
   MnvH1D* hists_US_data   = (MnvH1D*)f1->Get(Form("selected_data_reco_US_%s", var.c_str()));
   MnvH1D* hists_US_mat    = (MnvH1D*)f1->Get(Form("US_%s_%s_%s", mat.c_str(), trueZ.c_str(), var.c_str()));
   MnvH1D* hists_US_other  = (MnvH1D*)f1->Get(Form("US_other_%s_%s", trueZ.c_str(), var.c_str()));
   MnvH1D* hists_US_regUS  = (MnvH1D*)f1->Get(Form("US_regUS_%s_%s", trueZ.c_str(), var.c_str()));
   MnvH1D* hists_US_regDS  = (MnvH1D*)f1->Get(Form("US_regDS_%s_%s", trueZ.c_str(), var.c_str()));
 
  //save them into TObjArray
  TObjArray hlistUntunedUS;
  hlistUntunedUS.Add(hists_US_mat);
  hlistUntunedUS.Add(hists_US_other);
  hlistUntunedUS.Add(hists_US_regUS);
  hlistUntunedUS.Add(hists_US_regDS);
   
   MnvH1D* hists_DS_data   = (MnvH1D*)f1->Get(Form("selected_data_reco_DS_%s", var.c_str()));
   MnvH1D* hists_DS_mat    = (MnvH1D*)f1->Get(Form("DS_%s_%s_%s", mat.c_str(), trueZ.c_str(), var.c_str()));
   MnvH1D* hists_DS_other  = (MnvH1D*)f1->Get(Form("DS_other_%s_%s", trueZ.c_str(), var.c_str()));
   MnvH1D* hists_DS_regUS  = (MnvH1D*)f1->Get(Form("DS_regUS_%s_%s", trueZ.c_str(), var.c_str()));
   MnvH1D* hists_DS_regDS  = (MnvH1D*)f1->Get(Form("DS_regDS_%s_%s", trueZ.c_str(), var.c_str()));

  //save them into TObjArray
  TObjArray hlistUntunedDS;
  hlistUntunedDS.Add(hists_DS_mat);
  hlistUntunedDS.Add(hists_DS_other);
  hlistUntunedDS.Add(hists_DS_regUS);
  hlistUntunedDS.Add(hists_DS_regDS);
  
  PlotChi2Stat( mnvPlotter, hists_US_data, hlistUntunedUS, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "US_Untuned");
  PlotChi2Stat( mnvPlotter, hists_DS_data, hlistUntunedDS, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "DS_Untuned");

}
   }

for( unsigned int i = 0; i != vars.size(); ++i ){
  const string& var = vars[i];
  cout<<" variable in the loop = "<<var<<endl;

  histos_mc_US.push_back( hists_US_data[i] );
  histos_mc_US.push_back( hists_US_mat[i] );
  histos_mc_US.push_back( hists_US_other[i] );
  histos_mc_US.push_back( hists_US_regUS[i] );
  histos_mc_US.push_back( hists_US_regDS[i] );
      
  histos_mc_DS.push_back( hists_DS_data[i] );
  histos_mc_DS.push_back( hists_DS_mat[i] );
  histos_mc_DS.push_back( hists_DS_other[i] );
  histos_mc_DS.push_back( hists_DS_regUS[i] );
  histos_mc_DS.push_back( hists_DS_regDS[i] );
 
  PlotBGStuff(mnvPlotter, histos_mc_US, vars[i], targetID, targetZ, dataMCScale, true, false, playlist, "Untuned" );  //Filling US sidebands
  PlotBGStuff(mnvPlotter, histos_mc_DS, vars[i], targetID, targetZ, dataMCScale, false, true, playlist, "UnTuned");   //Filling DS sidebands


    //clear vector before going to the next variable 
    histos_mc_US.clear();
    histos_mc_DS.clear();
 
  }

// Now I am going to make tuned plots!!!! 
  vector< MnvH1D*> hists_US_mat_tuned, hists_US_other_tuned, hists_US_regUS_tuned, hists_US_regDS_tuned;
  vector< MnvH1D*> hists_DS_mat_tuned, hists_DS_other_tuned, hists_DS_regUS_tuned, hists_DS_regDS_tuned;
  vector< MnvH1D*> histos_mc_US_tuned, histos_mc_DS_tuned;
  vector< MnvH1D*> hists_US_data_tuned, hists_DS_data_tuned;
    bool plotUS = true;
    bool plotDS = true;
  for( vector<string>::iterator i = vars.begin(); i != vars.end(); ++i ){  
    const string& var = *i;
    TString histFileNameUS = Form("%s/TunedPlasticSidebands_US_t%d_z%02d_%s_%s.root", outdir.c_str(), targetID, targetZ, var.c_str(), playlist.c_str());
    TString histFileNameDS = Form("%s/TunedPlasticSidebands_DS_t%d_z%02d_%s_%s.root", outdir.c_str(), targetID, targetZ, var.c_str(), playlist.c_str());
 
    TFile *f_us = new TFile( histFileNameUS,"read" );
    TFile *f_ds = new TFile( histFileNameDS,"read" );

    f_us->cd();
    hists_US_data_tuned.push_back(  (MnvH1D*)f_us->Get("tuned_data"));
    hists_US_mat_tuned.push_back(   (MnvH1D*)f_us->Get("tuned_mc_signal"));
    hists_US_other_tuned.push_back( (MnvH1D*)f_us->Get("tuned_mc_other"));
    hists_US_regUS_tuned.push_back( (MnvH1D*)f_us->Get("tuned_mc_us"));
    hists_US_regDS_tuned.push_back( (MnvH1D*)f_us->Get("tuned_mc_ds"));

    hists_DS_data_tuned.push_back(  (MnvH1D*)f_ds->Get("tuned_data"));
    hists_DS_mat_tuned.push_back(   (MnvH1D*)f_ds->Get("tuned_mc_signal"));
    hists_DS_other_tuned.push_back( (MnvH1D*)f_ds->Get("tuned_mc_other"));
    hists_DS_regUS_tuned.push_back( (MnvH1D*)f_ds->Get("tuned_mc_us"));
    hists_DS_regDS_tuned.push_back( (MnvH1D*)f_ds->Get("tuned_mc_ds"));

{ //For ratio plotting
    f_us->cd();
    MnvH1D* hists_US_data_tuned   =  (MnvH1D*)f_us->Get("tuned_data");
    MnvH1D* hists_US_mat_tuned    =  (MnvH1D*)f_us->Get("tuned_mc_signal");
    MnvH1D* hists_US_other_tuned  =  (MnvH1D*)f_us->Get("tuned_mc_other");
    MnvH1D* hists_US_regUS_tuned  =  (MnvH1D*)f_us->Get("tuned_mc_us");
    MnvH1D* hists_US_regDS_tuned  =  (MnvH1D*)f_us->Get("tuned_mc_ds");
  
  //save them into TObjArray
  TObjArray hlisttunedUS;
  hlisttunedUS.Add(hists_US_mat_tuned);
  hlisttunedUS.Add(hists_US_other_tuned);
  hlisttunedUS.Add(hists_US_regUS_tuned);
  hlisttunedUS.Add(hists_US_regDS_tuned);
   
    f_ds->cd();
    MnvH1D*  hists_DS_data_tuned    = (MnvH1D*)f_ds->Get("tuned_data");
    MnvH1D*  hists_DS_mat_tuned     = (MnvH1D*)f_ds->Get("tuned_mc_signal");
    MnvH1D*  hists_DS_other_tuned   = (MnvH1D*)f_ds->Get("tuned_mc_other");
    MnvH1D*  hists_DS_regUS_tuned   = (MnvH1D*)f_ds->Get("tuned_mc_us");
    MnvH1D*  hists_DS_regDS_tuned   = (MnvH1D*)f_ds->Get("tuned_mc_ds");

  //save them into TObjArray
  TObjArray hlisttunedDS;
  hlisttunedDS.Add(hists_DS_mat_tuned);
  hlisttunedDS.Add(hists_DS_other_tuned);
  hlisttunedDS.Add(hists_DS_regUS_tuned);
  hlisttunedDS.Add(hists_DS_regDS_tuned);

/////CHi2 plotting

  PlotChi2Stat( mnvPlotter, hists_US_data_tuned, hlisttunedUS, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "US_tuned");
  PlotChi2Stat( mnvPlotter, hists_DS_data_tuned, hlisttunedDS, dataMCScale, var.c_str(), targetID, targetZ, false, true, playlist, "DS_tuned");


}//For ratio plotting

  }



  for( unsigned int i = 0; i != vars.size(); ++i ){  
    const string& var = vars[i];
    cout<<" variable in the after tuning loop = "<<var<<endl;
  
    histos_mc_US_tuned.push_back( hists_US_data_tuned[i] );
    histos_mc_US_tuned.push_back( hists_US_mat_tuned[i] );
    histos_mc_US_tuned.push_back( hists_US_other_tuned[i] );
    histos_mc_US_tuned.push_back( hists_US_regUS_tuned[i] );
    histos_mc_US_tuned.push_back( hists_US_regDS_tuned[i] );
      
    histos_mc_DS_tuned.push_back( hists_DS_data_tuned[i] );
    histos_mc_DS_tuned.push_back( hists_DS_mat_tuned[i] );
    histos_mc_DS_tuned.push_back( hists_DS_other_tuned[i] );
    histos_mc_DS_tuned.push_back( hists_DS_regUS_tuned[i] );
    histos_mc_DS_tuned.push_back( hists_DS_regDS_tuned[i] );

  
    PlotBGStuff(mnvPlotter, histos_mc_US_tuned, vars[i], targetID, targetZ, dataMCScale, true, false, playlist, "Tuned" );   //Filling US sidebands
    PlotBGStuff(mnvPlotter, histos_mc_DS_tuned, vars[i], targetID, targetZ, dataMCScale, false, true, playlist, "Tuned" );   //Filling DS sidebands
    
    //clear vector before going to the next variable 
      histos_mc_US_tuned.clear();
      histos_mc_DS_tuned.clear();
            
    }

}// End of Main

void PlotBGStuff(MnvPlotter mnvPlotter, vector<MnvH1D*> histos_mc, const string var, int targetID, int targetZ, double dataMCScale, bool plotUS, bool plotDS, string playlist, string Name){
 
  MnvH1D* data_sideband       = (MnvH1D*)histos_mc[0]->Clone("data");
  MnvH1D* mc_sideband_signal  = (MnvH1D*)histos_mc[1]->Clone("mc_signal");
  MnvH1D* mc_sideband_Other   = (MnvH1D*)histos_mc[2]->Clone("mc_other");
  MnvH1D* mc_sideband_regUS   = (MnvH1D*)histos_mc[3]->Clone("mc_us");
  MnvH1D* mc_sideband_regDS   = (MnvH1D*)histos_mc[4]->Clone("mc_ds");
  //MnvH1D* mc_sideband_regDS   = (MnvH1D*)histos_mc[4]->Clone("mc_us");

  string suffix = "";
  PlotUtils::HistFolio<PlotUtils::MnvH1D> hlist_grouped(Form("%s_%s", suffix.c_str(), var.c_str()), mc_sideband_signal);
    hlist_grouped.AddComponentHist("Signal", mc_sideband_signal);
    hlist_grouped.AddComponentHist("Other",  mc_sideband_Other);
    hlist_grouped.AddComponentHist("regUS", mc_sideband_regUS);
    hlist_grouped.AddComponentHist("regDS", mc_sideband_regDS);

  string name;
  if(plotUS) name = Form("US_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
  if(plotDS) name = Form("DS_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
 

  // Plotting sidebands  
  //PlotStacked(data_sideband, hlist, dataMCScale, var.c_str(), name, "");
  PlotStacked(data_sideband, hlist_grouped.GetHistArray(), dataMCScale, var.c_str(), name, " ");


}

void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name ){
  
  //first declared some plotting parameters
  mnvPlotter.legend_offset_x = .03;
  mnvPlotter.width_xspace_per_letter = .33;
  /*
  string targetString;  
  if( targetZ == 6 )
    targetString = "Carbon";
  else if( targetZ == 26 )
    targetString = "Iron";
  if( targetZ == 82 )
    targetString = "Lead";
 */
   
  //just get the CV for now 
  MnvH1D* data_cv = (MnvH1D*)dataHisto->Clone("mnvh1d_data_histo");
  cout << "data_cv = " << data_cv->Integral() << endl; 
  //get total MC
  MnvH1D* totalMC = (MnvH1D*)mchistos[0]->Clone("mnvh1d_mc_histo");
  totalMC->UseCurrentStyle();

  /*string cname = (string)canname;
  if( targetZ == 6 && string::npos != cname.find("_US_") && string::npos != cname.find("_planeDNN_") ) totalMC->GetXaxis()->SetRangeUser(22., 28.);
  else if( targetZ == 6 && string::npos != cname.find("_DS_") && string::npos != cname.find("_planeDNN_") ) totalMC->GetXaxis()->SetRangeUser(34., 40.);
  */
  for( int i = 1; i < mchistos.GetLast()+1; i++ ){
    MnvH1D* mc_histo_add = (MnvH1D*)mchistos[i]->Clone("mnvh1d_mc_histo_add");
    mc_histo_add->SetFillColor(0);
    mc_histo_add->SetFillStyle(0);
    
    totalMC->Add(mc_histo_add);
  }
  
  cout << "total MC = " << totalMC->Integral() << endl; 
  if( data_cv->GetSumw2N() == 0 ) 
    data_cv->Sumw2();
  if( totalMC->GetSumw2N() == 0 ){ 
    totalMC->Sumw2();
  }


  
  int ndfStat = 1;
  double chi2Sys = mnvPlotter.Chi2DataMC( data_cv, totalMC, ndfStat, dataMCScale );
  ndfStat -= 1;
  
  //TString labelStat = Form("Stat. + Sys. Error, #chi^{2}/ndf = %3.2f/%d = %3.2f", chi2Sys, ndfStat, chi2Sys/(Double_t)ndfStat);
  TString labelStat = Form("Stat. Error, #chi^{2}/ndf = %3.2f/%d = %3.2f", chi2Sys, ndfStat, chi2Sys/(Double_t)ndfStat);
  
 // string region = "Fiducial";
 // if( string::npos != cname.find("_US_") ) region = "Upstream Plastic Background";
 // if( string::npos != cname.find("_DS_") ) region = "Downstream Plastic Background";
  
 // string status = "After Tuning";
 // if( string::npos != cname.find("BeforeTuning_") ) status = "Before Tuning";

 // if( writeChi2 && !tune_fiducial ) fprintf(Chi2Table, "%s %s %s %.2f \n", var.c_str(), region.c_str(), status.c_str(), chi2Sys/(Double_t)ndfStat ); 
 // if( writeChi2 && tune_fiducial  ) fprintf(Chi2Table, "%s Fiducial %s, With Sys. %.2f \n", var.c_str(), status.c_str(), chi2Sys/(Double_t)ndfStat ); 

  double plotMin = -1., plotMax = -1.;
  GetMinAndMaxAxisRatio( data_cv, totalMC, dataMCScale, plotMin, plotMax );
  
  //plot data/mc with stat errors
  {
  TString cName;
  if(plotUS) cName = Form("Ratio_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
  if(plotDS) cName = Form("Ratio_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
    //TString cName = Form("BGPlotRatio_%s",var.c_str());
    //TString Yaxis = "Data/MC (" + status + ")"; 
    TString Yaxis = "Data/MC"; 
    std::string x_label = var.c_str();
    TCanvas c( cName, cName, 1200, 800 );
    cout << "plotting " << cName << endl;
    //mnvPlotter.DrawDataMCRatio( data_cv, totalMC, dataMCScale, true, true, plotMin, plotMax, Yaxis );
    mnvPlotter.DrawDataMCRatio( data_cv, totalMC, dataMCScale, true, true, 0.6, 1.4, x_label.c_str() );
    mnvPlotter.AddPlotLabel( labelStat, .46, .875, .05 );
    
    if( WRITE_PRELIMINARY ) mnvPlotter.WritePreliminary("BC");
    
    //if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("%s of %s", region.c_str(), targetString.c_str() ), .06);
    mnvPlotter.MultiPrint( &c );
  }
  
  {    
    TString cName;
   if(plotUS) cName = Form("FracError_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
   if(plotDS) cName = Form("FracError_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
    TCanvas can( cName, cName, 1200,800 );
    
    cout << "Fill color: "<< totalMC->GetFillColor() << endl;
 
    can.SetFillColor(0);
    can.SetFillStyle(4050);
    totalMC->UseCurrentStyle();
    
  string legLoc  = "TL";
  string prelimLoc = "TR";
  
    cout << "plotting " << cName << endl;

    //PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
    plotter.DrawErrorSummary( totalMC, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
    
    
    if(WRITE_PRELIMINARY)
      mnvPlotter.WritePreliminary(prelimLoc);
    
    mnvPlotter.MultiPrint( &can, can.GetName() );
  }
  //scale to data
  //totalMC->Scale(dataMCScale);
}

void GetMinAndMaxAxisRatio( TH1D* histData, TH1D* histMC, double dataMCScale, double &plotMin, double &plotMax ){
  
  TH1D* ratio = (TH1D*)histData->Clone("ratio_plot");
  TH1D* totalMC = (TH1D*)histMC->Clone("mc_clone");
  totalMC->Scale(dataMCScale);
  ratio->Divide( histData, totalMC);
  
  int firstBin = ratio->GetXaxis()->GetFirst();
  int lastBin = ratio->GetXaxis()->GetLast();
  
  std::map<double, int> bins;
  for( int b = firstBin; b <= lastBin; b++ ){
    bins.insert( std::pair<double,int>(ratio->GetBinContent(b), b) ); 
  }
  
  vector<int> bins_in_order;
  for( std::map<double,int>::iterator it = bins.begin(); it != bins.end(); ++it ){
    bins_in_order.push_back( it->second );
  }
  
  int last_bin = bins_in_order.size() - 1;
  plotMax = ratio->GetBinContent( bins_in_order[last_bin] ) + ratio->GetBinError( bins_in_order[last_bin] ) + 0.3;
  if( plotMax > 4.0 ) plotMax = ratio->GetBinContent( bins_in_order[last_bin-1] ) + ratio->GetBinError( bins_in_order[last_bin-1] ) + 0.1;
  plotMin = ratio->GetBinContent( bins_in_order[0] ) - ratio->GetBinError( bins_in_order[0] ) - 0.1;
  if( plotMin <= -0.1 ) plotMin = ratio->GetBinContent( bins_in_order[1] ) - ratio->GetBinError( bins_in_order[1] ) - 0.1;
  
}





