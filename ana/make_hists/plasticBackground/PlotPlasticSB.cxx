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

#include "../../include/systematics/Systematics.h"

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
void PlotStackedSideband(MnvPlotter mnvPlotter, vector<MnvH1D*> histos_mc, const string var, int targetID, int targetZ, double dataMCScale, bool plotUS, bool plotDS, string playlist, string Name);
void PlotRatiowithChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name);
void PlotFracUncertainty( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name);
void GetMinAndMaxAxisRatio( TH1D* histData, TH1D* histMC, double dataMCScale, double &plotMin, double &plotMax );
//============================================================================================================================
// Main
int main(int argc, char * argv[]){
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
  
  TString histFileName, SF_file;

  cout<<"******************************************************************************************"<<endl;
  cout<<"                   I am making plots for both before and after tuning!!                   "<<endl;
  cout<<"******************************************************************************************"<<endl;

  if(RunCodeWithSystematics){
    histFileName = Form("%s/Hists_PlasticBkg_sys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ ); 
  }

  else{
    histFileName = Form("%s/Hists_PlasticBkg_nosys_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ, getenv("NUKECC_TAG"));
  } 
   
  cout<<histFileName<<endl;

  MnvPlotter mnvPlotter(kCCQENuStyle);

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
  vars.push_back("Enu");
  vars.push_back("x");
  vars.push_back("planeDNN");

  TFile *f1 = new TFile( histFileName,"read" );

  TParameter<double> *mcPOT = (TParameter<double>*)f1->Get("MCPOT");
  TParameter<double> *dataPOT = (TParameter<double>*)f1->Get("DataPOT");
  double mcpot = mcPOT->GetVal();
  double datapot = dataPOT->GetVal(); 
  double dataMCScale = datapot/mcpot;
  cout<<"MCPOT = "<<mcpot<<"DataPOT = "<< datapot << "Scale = " << dataMCScale<<endl;

  //---------------------------------------------------------------------------------------------------
  // ------------------------ UNTUNED PLOTS -------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------

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
    
    PlotRatiowithChi2Stat( mnvPlotter, hists_US_data, hlistUntunedUS, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "US_untuned");
    PlotFracUncertainty( mnvPlotter, hists_US_data, hlistUntunedUS, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "US_untuned");

    PlotRatiowithChi2Stat( mnvPlotter, hists_DS_data, hlistUntunedDS, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "DS_untuned");
    PlotFracUncertainty( mnvPlotter, hists_DS_data, hlistUntunedDS, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "DS_untuned");

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
  
    PlotStackedSideband(mnvPlotter, histos_mc_US, vars[i], targetID, targetZ, dataMCScale, true, false, playlist, "US_untuned" );  //Filling US sidebands
    PlotStackedSideband(mnvPlotter, histos_mc_DS, vars[i], targetID, targetZ, dataMCScale, false, true, playlist, "DS_untuned");   //Filling DS sidebands


    //clear vector before going to the next variable 
    histos_mc_US.clear();
    histos_mc_DS.clear();
 
  }

  //---------------------------------------------------------------------------------------------------
  // ------------------------ TUNED PLOTS -------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------
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

    //For ratio plotting
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

    // Tuned Ratio with Chi2
    // Tuned Fractional uncertainty

    PlotRatiowithChi2Stat( mnvPlotter, hists_US_data_tuned, hlisttunedUS, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "US_tuned");
    PlotFracUncertainty( mnvPlotter, hists_US_data_tuned, hlisttunedUS, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "US_tuned");
    
    PlotRatiowithChi2Stat( mnvPlotter, hists_DS_data_tuned, hlisttunedDS, dataMCScale, var.c_str(), targetID, targetZ, false, true, playlist, "DS_tuned");
    PlotFracUncertainty( mnvPlotter, hists_DS_data_tuned, hlisttunedDS, dataMCScale, var.c_str(), targetID, targetZ, false, true, playlist, "DS_tuned");

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

    PlotStackedSideband(mnvPlotter, histos_mc_US_tuned, vars[i], targetID, targetZ, dataMCScale, true, false, playlist, "US_tuned" );   //Filling US sidebands
    PlotStackedSideband(mnvPlotter, histos_mc_DS_tuned, vars[i], targetID, targetZ, dataMCScale, false, true, playlist, "DS_tuned" );   //Filling DS sidebands
    
    //clear vector before going to the next variable 
    histos_mc_US_tuned.clear();
    histos_mc_DS_tuned.clear();
            
  }

}// End of Main


void PlotStackedSideband(MnvPlotter mnvPlotter, vector<MnvH1D*> histos_mc, const string var, int targetID, int targetZ, double dataMCScale, bool plotUS, bool plotDS, string playlist, string Name){
 
  MnvH1D* data_sideband       = (MnvH1D*)histos_mc[0]->Clone("data");
  MnvH1D* mc_sideband_signal  = (MnvH1D*)histos_mc[1]->Clone("mc_signal");
  MnvH1D* mc_sideband_Other   = (MnvH1D*)histos_mc[2]->Clone("mc_other");
  MnvH1D* mc_sideband_regUS   = (MnvH1D*)histos_mc[3]->Clone("mc_us");
  MnvH1D* mc_sideband_regDS   = (MnvH1D*)histos_mc[4]->Clone("mc_ds");

  data_sideband->GetYaxis()->CenterTitle();
  mc_sideband_signal->GetYaxis()->CenterTitle();
  mc_sideband_Other->GetYaxis()->CenterTitle();
  mc_sideband_regUS->GetYaxis()->CenterTitle();
  mc_sideband_regDS->GetYaxis()->CenterTitle();

  string suffix = "";
  PlotUtils::HistFolio<PlotUtils::MnvH1D> hlist_grouped(Form("%s_%s", suffix.c_str(), var.c_str()), mc_sideband_signal);
  hlist_grouped.AddComponentHist("Signal", mc_sideband_signal);
  hlist_grouped.AddComponentHist("Other",  mc_sideband_Other);
  hlist_grouped.AddComponentHist("Upstream", mc_sideband_regUS);
  hlist_grouped.AddComponentHist("Downstream", mc_sideband_regDS);

  //hlist_grouped.GetHistArray()
  TObjArray array = *(TObjArray*)hlist_grouped.GetHistArray().Clone("mc");

  string outfile_name;
  if(plotUS) outfile_name = Form("StackedBreakdown_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
  if(plotDS) outfile_name = Form("StackedBreakdown_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
  TCanvas cE("c1", "c1");

  string x_label;
  if (var == "Enu") {x_label = "Neutrino Energy [GeV]";}
  if (var == "x") {x_label= "Bjorken x";}
  if (var == "planeDNN") {x_label = "Plane Number";}

  string y_label;
  if (var == "Enu") {y_label = "Events/GeV";}
  if (var == "x") {y_label= "Events (norm.)";}
  if (var == "planeDNN") {y_label = "Events (norm.)";}

  mnvPlotter.SetLegendNColumns(1);
  mnvPlotter.DrawDataStackedMC(data_sideband, &array, dataMCScale, "TL", "Data", -2, -2, 3001, x_label.c_str(), y_label.c_str());
  mnvPlotter.SetLegendNColumns(1);
  //mnvPlotter.SetLegendNColumns(1);
                                                                    // base color, offset colot, fill style
  if (Name == "DS_tuned") mnvPlotter.AddHistoTitle("Tuned Downstream ", 0.05);
  if (Name == "DS_untuned") mnvPlotter.AddHistoTitle("Untuned Downstream ", 0.05);
  if (Name == "US_tuned") mnvPlotter.AddHistoTitle("Tuned Upstream ", 0.05);
  if (Name == "US_untuned") mnvPlotter.AddHistoTitle("Untuned Upstream ", 0.05);
  if( WRITE_PRELIMINARY ) mnvPlotter.WritePreliminary("TR");
  //mnv_plotter.WriteNorm("Abs-Normalized", 0.8, 0.6);
  mnvPlotter.MultiPrint(&cE, outfile_name, "png");

}

void PlotRatiowithChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name ){
  //first declared some plotting parameters
  mnvPlotter.legend_offset_x = .03;
  mnvPlotter.width_xspace_per_letter = .33;
   
  //just get the CV for now 
  MnvH1D* data_cv = (MnvH1D*)dataHisto->Clone("mnvh1d_data_histo");
  cout << "data_cv = " << data_cv->Integral() << endl; 
  //get total MC
  MnvH1D* totalMC = (MnvH1D*)mchistos[0]->Clone("mnvh1d_mc_histo");
  totalMC->UseCurrentStyle();

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

  // Calculate chi2 stats
  int ndfStat = 1;
  double chi2Sys = mnvPlotter.Chi2DataMC( data_cv, totalMC, ndfStat, dataMCScale );
  ndfStat -= 1;
  
  TString labelStat;
  if (RunCodeWithSystematics){
    labelStat = Form("Stat. + Sys. Error, #chi^{2}/ndf = %3.2f/%d = %3.2f", chi2Sys, ndfStat, chi2Sys/(Double_t)ndfStat);
  }
  else {
    labelStat = Form("Stat. Error, #chi^{2}/ndf = %3.2f/%d = %3.2f", chi2Sys, ndfStat, chi2Sys/(Double_t)ndfStat);
  }

  double plotMin = -1., plotMax = -1.;
  GetMinAndMaxAxisRatio( data_cv, totalMC, dataMCScale, plotMin, plotMax );
  
  //plot Data CV/total MC ratio with stat/sys error
  std::string  cName = Form("Ratio_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());

  TCanvas c("c1", "c1", 1200, 800 );
  cout << "plotting " << cName << endl;
  if (var == "Enu") {totalMC->GetXaxis()->SetTitle("Neutrino Energy [GeV]");}
  if (var == "x") {totalMC->GetXaxis()->SetTitle("Bjorken x");}
  if (var == "planeDNN") {totalMC->GetXaxis()->SetTitle("Plane Number");}
  totalMC->GetXaxis()->CenterTitle();
  totalMC->GetYaxis()->CenterTitle();
  mnvPlotter.DrawDataMCRatio( data_cv, totalMC, dataMCScale, true, true, 0.5, 1.5, "Data/Total MC" );
  mnvPlotter.AddPlotLabel( labelStat, .41, .875, .04 );
  
  if( WRITE_PRELIMINARY ) mnvPlotter.WritePreliminary("BR");

  if (Name == "DS_tuned") mnvPlotter.AddHistoTitle("Tuned Downstream ", 0.05);
  if (Name == "DS_untuned") mnvPlotter.AddHistoTitle("Untuned Downstream ", 0.05);
  if (Name == "US_tuned") mnvPlotter.AddHistoTitle("Tuned Upstream ", 0.05);
  if (Name == "US_untuned") mnvPlotter.AddHistoTitle("Untuned Upstream ", 0.05);

  mnvPlotter.MultiPrint( &c, cName, "png");
  
}


void PlotFracUncertainty( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name ){
  //first declared some plotting parameters
  mnvPlotter.legend_offset_x = .03;
  mnvPlotter.width_xspace_per_letter = .33;

  //just get the CV for now 
  MnvH1D* data_cv = (MnvH1D*)dataHisto->Clone("mnvh1d_data_histo");
  cout << "data_cv = " << data_cv->Integral() << endl; 
  //get total MC
  MnvH1D* totalMC = (MnvH1D*)mchistos[0]->Clone("mnvh1d_mc_histo");

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

  std::string cName;
  if(plotUS) cName = Form("FracError_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
  if(plotDS) cName = Form("FracError_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());

  TCanvas c("c1", "c1", 1200, 800 );
  cout << "plotting " << cName << endl;
  if (var == "Enu") {totalMC->GetXaxis()->SetTitle("Neutrino Energy [GeV]");}
  if (var == "x") {totalMC->GetXaxis()->SetTitle("Bjorken x");}
  if (var == "planeDNN") {totalMC->GetXaxis()->SetTitle("Plane Number");}
  totalMC->GetXaxis()->CenterTitle();
  totalMC->GetYaxis()->CenterTitle();
  
  cout << "Fill color: "<< totalMC->GetFillColor() << endl;
  c.SetFillColor(0);
  c.SetFillStyle(4050);
    
  string legLoc  = "TL";

  cout << "plotting " << cName << endl;

  mnvPlotter.ApplyStyle(kCCQENuStyle);
  //mnvPlotter.ApplyStyle(kDefaultStyle);
  //mnvPlotter.ApplyStyle(kNukeCCStyle);
  //mnvPlotter.ApplyStyle(7);

  // Use my own grouping consistent with the rest
  // Should eventually put it in the mnvPlotter
  mnvPlotter.error_summary_group_map.clear();
  mnvPlotter.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_Resolution");
  mnvPlotter.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_MINOS");
  mnvPlotter.error_summary_group_map["Muon Energy"].push_back("Muon_Energy_MINERvA");
  
  mnvPlotter.error_summary_group_map["MINOS"].push_back("MINOS_Reconstruction_Efficiency");

  mnvPlotter.error_summary_group_map["Muon Angle"].push_back("BeamAngleX");
  mnvPlotter.error_summary_group_map["Muon Angle"].push_back("BeamAngleY");

  mnvPlotter.error_summary_group_map["Flux"].push_back("Flux");

  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("RPA_HighQ2");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("RPA_LowQ2");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("Low_Recoil_2p2h_Tune");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_AhtBY");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_BhtBY");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_CCQEPauliSupViaKF");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_CV1uBY");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_CV2uBY");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_EtaNCEL");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_MaCCQE");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_MaNCEL");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_MaRES");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_MvRES");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_NormDISCC");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_NormNCRES");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvn1pi");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvn2pi");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvp1pi");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_Rvp2pi");
  mnvPlotter.error_summary_group_map["Interaction Model"].push_back("GENIE_VecFFCCQEshape");

  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_AGKYxF1pi");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_RDecBR1gamma");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_Theta_Delta2Npi");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_FrAbs_N");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_FrAbs_pi");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_FrCEx_N");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_FrCEx_pi");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_FrElas_N");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_FrElas_pi");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_FrInel_N");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_FrPiProd_N");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_FrPiProd_pi");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_MFP_N");
  mnvPlotter.error_summary_group_map["FSI"].push_back("GENIE_MFP_pi");

  mnvPlotter.error_summary_group_map["Target Mass"].push_back("Target_Mass");

  totalMC->GetYaxis()->CenterTitle();

  mnvPlotter.SetLegendNColumns(2);
  mnvPlotter.DrawErrorSummary( totalMC, legLoc, true, true, 0.0001, true, "", true, var.c_str(), false); 
  // cov area normalised doesn't make a difference for fractional error

  //bool MnvPlotter::DrawErrorSummary(
  //        MnvH1D* h,
  //        const std::string& legPos  /* = "TR"    */,
  //        const bool   includeStat     /* = true    */,
  //        const bool   solidLinesOnly  /* = true    */,
  //        const double ignoreThreshold /* = 0.00001 */,
  //        const bool covAreaNormalize /* = false*/,
  //        const std::string& errorGroupName /* = "" */,
  //        const bool  asfrac  /* = false */,
  //        const std::string &Ytitle,
  //        bool ignoreUngrouped
  //        )

  mnvPlotter.SetLegendNColumns(2);
    
  //if( WRITE_PRELIMINARY ) mnvPlotter.WritePreliminary("TR");
  if(ADD_HISTO_TITLES){
    if (Name == "DS_tuned") {mnvPlotter.AddHistoTitle("Tuned Downstream ", 0.05);}
    if (Name == "DS_untuned") {mnvPlotter.AddHistoTitle("Untuned Downstream ", 0.05);}
    if (Name == "US_tuned") {mnvPlotter.AddHistoTitle("Tuned Upstream ", 0.05);}
    if (Name == "US_untuned") {mnvPlotter.AddHistoTitle("Untuned Upstream ", 0.05);}
  }
  mnvPlotter.MultiPrint( &c, cName, "png" );
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





