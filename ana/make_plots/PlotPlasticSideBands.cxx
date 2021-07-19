#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
#include "include/NukeCCUtilsNSF.h"
//#include "syst_common.h"
#include "TArrayD.h"
#include <iostream>
                              
#include "Math/Factory.h"     
#include "Math/Functor.h"     
#include "Math/Minimizer.h"   

#ifndef __CINT__
#endif

using namespace NUKECC_ANA;
using namespace PlotUtils;
using namespace std;

//plotting parameters
double statX = .68, statY = .44;
string legPos = "TL", normPos="TR";
string legLoc  = "TL";
string prelimLoc = "TR";
double infoX = .68;
double infoY = .78;

const bool makeDISPlots = true; //set to true to make DIS plots!
const bool makeLowQ2Plots = false; //set to true to make DIS plots!
const bool makeLowWPlots = false; //set to true to make DIS plots!
const bool onlyCalculateCV = true; //set to false to do the tuning in systematics universe - set to true if you're not calculating systematics
const bool writeSysCode = false; //set to true to write the systematics scale factors code
const bool writeChi2 = true;// set to true to make file with chi2s
const bool tune_fiducial = false; //set to true to see effects in passive targets
FILE *Chi2Table;
string outdir;

//----Helper Functions-----------------------------------------------
//void SetHistosAxis( MnvH1D* histo, string var, const bool makeDISPlots );
//void SetStackedHistosMaxAxis( MnvPlotter &mnvPlotter, int targetZ, string var, const bool makeDISPlots, int opt );
void DrawMeSomeDataMCStackedHistos( MnvPlotter mnvPlotter, MnvH1D* histo_data, TObjArray hlist, TString canname, string plotname, string outdir, string var, int targetZ, const bool makeDISPlots );
void DrawTunedBG( MnvPlotter mnvPlotter, vector<MnvH1D*> histos, const string var, string plot_title, TString canname, int targetZ,  MnvH1D* scaleFactorHist, bool fixUS, bool fixDS, string playlist, string outdir);
void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, string var, int targetZ, TString canname );
void GetMinAndMaxAxisRatio( TH1D* histData, TH1D* histMC,  double &plotMin, double &plotMax );
//void TuneFiducial( MnvPlotter &mnvPlotter, vector<MnvH1D*> histos_mc_tuning, int targetZ, string var, TString canname, string plot_title, double dataMCScale, string playlist );
//void SubtractFiducial( MnvPlotter &mnvPlotter, vector<MnvH1D*> histos_mc_tuning, int targetZ, string var, TString canname, string plot_title, double dataMCScale, string playlist );
void DrawSomeMCStackedHistos( MnvPlotter mnvPlotter, TObjArray hlist, double scale, TString canname, string plotname, string var, int targetZ, const bool makeDISPlots );
void TuneScintBackground(MnvH1D* mnvh1d_data, MnvH1D* mnvh1d_mc_signal, MnvH1D* mnvh1d_mc_other, MnvH1D* mnvh1d_mc_background_US, MnvH1D* mnvh1d_mc_background_DS,MnvH1D* mc_signal, MnvH1D* mc_other, MnvH1D* mc_background_US, MnvH1D* mc_background_DS, MnvH1D* tmp_scale_factor,  bool fixUS, bool fixDS, string channel, bool onlyCV , bool writeScaleFactor );
vector<double> CalcScaleFactorMinimizer( TH1D* h1d_data, TH1D* h1d_mc_signal, TH1D* h1d_mc_other, TH1D* h1d_mc_backgroundUS, TH1D* h1d_mc_backgroundDS, const std::string& errName, bool fixUS, bool fixDS);
double getChi2( const double * par );
      TH1D* m_histo_data;
      TH1D* m_histo_sig_fix;
      TH1D* m_histo_other_fix;
      TH1D* m_histo_plastic_fix;
      TH1D* m_histo_plastic_float;
      TH1D* m_histo_plastic_us_float;
      TH1D* m_histo_plastic_ds_float;
//void AddHistosToAVectorUS( vector<MnvH1D*> &v_histos, vector<MnvH1D*> histos, int targetZ ); 
//void AddHistosToAVectorDS( vector<MnvH1D*> &v_histos, vector<MnvH1D*> histos, int targetZ ); 

//void DrawFracMCHistos( MnvPlotter mnvPlotter, TObjArray hlist, TString canname, string plotname, string var, int targetZ, const bool makeDISPlots );
//-------------------------------------------------------------------
int main(int argc, char * argv[]){
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if(argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./runEventLoop Path_to_input_files Path_to_Output_plots Target_number Material_atomic_number Playlist\n\n"<<
      "\t-Path_to_input_file\t =\t Path to the directory with histogram and scale factors to plot \n"<<
       "\t-Path_to_Output_plots\t =\t Path to the directory where the plots will be created \n"<<
      //"\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" <<
      "\t-Playlist\t \t =\t eg. minervame1A"<< std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    return 0;
  }


  string indir=argv[1];
  string outdir=argv[2];
  //  int targetID = atoi(argv[2]);
  int targetZ = atoi(argv[3]);
  const string playlist= argv[4];

  TString hname1 = Form("/minerva/data/users/fakbar/NukeHists/v1_MAD/Hists_EventSelection_t1_z26_Nu_v1_.root");//Form("%s/Hists_EventSelection_t1_z%d_Nu_v1_.root",indir.c_str(),targetZ); //or change that to whatever we start calling the background histogram file. 
  cout<<hname1<<endl;

  vector< MnvH1D* > hists_sideband_US, hists_sideband_DS;
  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(playlist);

  MnvPlotter mnvPlotter(kNukeCCStyle);
  TFile *f1 = new TFile( hname1,"read" );
  f1->cd();
  
  if( writeChi2 ) Chi2Table = fopen( Form("%s/Chi2Sideband_z%02d.txt", outdir.c_str(), targetZ), "w" );
  
  MnvH1D *planeDNN_US_data  = (MnvH1D*)f1->Get("planeDNN_Data_US");
  MnvH1D *planeDNN_US_Fe  = (MnvH1D*)f1->Get("planeDNN_US_Fe");
  MnvH1D *planeDNN_US_Pb  = (MnvH1D*)f1->Get("planeDNN_US_Pb");
  MnvH1D *planeDNN_US_C  = (MnvH1D*)f1->Get("planeDNN_US_C");
  MnvH1D *planeDNN_US_other  = (MnvH1D*)f1->Get("planeDNN_US_other");
  MnvH1D *planeDNN_US_regUS  = (MnvH1D*)f1->Get("planeDNN_US_regUS");
  MnvH1D *planeDNN_US_regDS  = (MnvH1D*)f1->Get("planeDNN_US_regDS");
  
  MnvH1D *planeDNN_DS_data = (MnvH1D*)f1->Get("planeDNN_Data_DS");
  MnvH1D *planeDNN_DS_Fe = (MnvH1D*)f1->Get("planeDNN_DS_Fe");
  MnvH1D *planeDNN_DS_Pb = (MnvH1D*)f1->Get("planeDNN_DS_Pb");
  MnvH1D *planeDNN_DS_C = (MnvH1D*)f1->Get("planeDNN_DS_C");
  MnvH1D *planeDNN_DS_other  = (MnvH1D*)f1->Get("planeDNN_DS_other");
  MnvH1D *planeDNN_DS_regUS  = (MnvH1D*)f1->Get("planeDNN_DS_regUS");
  MnvH1D *planeDNN_DS_regDS  = (MnvH1D*)f1->Get("planeDNN_DS_regDS");
  
  MnvH1D *planeDNN_tgt_data = (MnvH1D*)f1->Get("planeDNN_Data_tgt");
  MnvH1D *planeDNN_tgt_Fe = (MnvH1D*)f1->Get("planeDNN_tgt_Fe");
  MnvH1D *planeDNN_tgt_Pb = (MnvH1D*)f1->Get("planeDNN_tgt_Pb");
  MnvH1D *planeDNN_tgt_C = (MnvH1D*)f1->Get("planeDNN_tgt_C");
  MnvH1D *planeDNN_tgt_other  = (MnvH1D*)f1->Get("planeDNN_tgt_other");
  MnvH1D *planeDNN_tgt_regUS  = (MnvH1D*)f1->Get("planeDNN_tgt_regUS");
  MnvH1D *planeDNN_tgt_regDS  = (MnvH1D*)f1->Get("planeDNN_tgt_regDS");
  
  double datapot = utils->getPOTData(f1);
  double mcpot = utils->getPOTMC(f1);
  double dataMCpot = datapot/mcpot;
  cout<<"Data POT = "<<datapot<<"   MC POT = "<<mcpot<<"   Data/MC POT = "<<datapot/mcpot<<endl;
  
  string trueZ;
  if( targetZ == 26 ) trueZ = "Iron";
  if( targetZ == 82 ) trueZ = "Lead";
  if( targetZ == 6 ) trueZ = "Carbon";

  //Scale all plots by data/MC POT ratio
  planeDNN_US_Fe->Scale(dataMCpot);
  planeDNN_US_Pb->Scale(dataMCpot);
  planeDNN_US_C->Scale(dataMCpot);
  planeDNN_US_regUS->Scale(dataMCpot);
  planeDNN_US_regDS->Scale(dataMCpot);
  planeDNN_US_other->Scale(dataMCpot);
  planeDNN_DS_Fe->Scale(dataMCpot);  
  planeDNN_DS_Pb->Scale(dataMCpot);  
  planeDNN_DS_C->Scale(dataMCpot);  
  planeDNN_DS_regUS->Scale(dataMCpot);  
  planeDNN_DS_regDS->Scale(dataMCpot);  
  planeDNN_DS_other->Scale(dataMCpot);  
  planeDNN_tgt_Fe->Scale(dataMCpot);  
  planeDNN_tgt_Pb->Scale(dataMCpot);  
  planeDNN_tgt_C->Scale(dataMCpot);  
  planeDNN_tgt_regUS->Scale(dataMCpot);  
  planeDNN_tgt_regDS->Scale(dataMCpot);  
  planeDNN_tgt_other->Scale(dataMCpot);  
  //-------------------------------------------------------------------------------plain plots

  //Upstream sideband
  THStack *us = new THStack("us","Stacked 1D US SB histograms");
  
  planeDNN_US_Fe->SetFillColor(kRed);
  planeDNN_US_Fe->SetMarkerColor(kRed);
  planeDNN_US_Fe->SetTitle("Iron");
  us->Add(planeDNN_US_Fe);
  
  planeDNN_US_Pb->SetFillColor(kBlue);
  planeDNN_US_Pb->SetMarkerColor(kBlue);
  planeDNN_US_Pb->SetTitle("Lead");
  us->Add(planeDNN_US_Pb);
  
  planeDNN_US_C->SetFillColor(kCyan+2);
  planeDNN_US_C->SetMarkerColor(kCyan+2);
  planeDNN_US_C->SetTitle("Carbon");
  us->Add(planeDNN_US_C);
  
  planeDNN_US_regUS->SetFillColor(kMagenta);
  planeDNN_US_regUS->SetMarkerColor(kMagenta);
  planeDNN_US_regUS->SetTitle("Upstream");
  us->Add(planeDNN_US_regUS);
  
  planeDNN_US_regDS->SetFillColor(kOrange+1);
  planeDNN_US_regDS->SetMarkerColor(kOrange+1);
  planeDNN_US_regDS->SetTitle("Downstream");
  us->Add(planeDNN_US_regDS);
  
  planeDNN_US_other->SetFillColor(kMagenta+2);
  planeDNN_US_other->SetMarkerColor(kMagenta+2);
  planeDNN_US_other->SetTitle("Other");
  us->Add(planeDNN_US_other);
  
  TCanvas *USstacked = new TCanvas("USstacked","stacked hists for Upstream SB",10,10,1200,800);
  USstacked->SetGrid(1);
  planeDNN_US_data->Draw("p");   
  us->Draw("hist same");
  
  planeDNN_US_data->GetXaxis()->SetTitle( "PlaneDNN"  );
  
  gPad->SetGrid(1);
  gPad->Print(Form("%s/US_planeDNN_%s.png",outdir.c_str(),trueZ.c_str()),"png");
  gPad->Update();
  
  USstacked->Update();
      
  THStack *ds = new THStack("ds","Stacked 1D DS SB histograms");
  
  planeDNN_DS_Fe->SetFillColor(kRed);
  planeDNN_DS_Fe->SetLineColor(kRed);
  planeDNN_DS_Fe->SetTitle("Iron");
  ds->Add(planeDNN_DS_Fe);
  
  planeDNN_DS_Pb->SetFillColor(kBlue);
  planeDNN_DS_Pb->SetLineColor(kBlue);
  planeDNN_DS_Pb->SetTitle("Lead");
  ds->Add(planeDNN_DS_Pb);
  
  planeDNN_DS_C->SetFillColor(kCyan+2);
  planeDNN_DS_C->SetLineColor(kCyan+2);
  planeDNN_DS_C->SetTitle("Carbon");
  ds->Add(planeDNN_DS_C);
  
  planeDNN_DS_regUS->SetFillColor(kMagenta);
  planeDNN_DS_regUS->SetLineColor(kMagenta);
  planeDNN_DS_regUS->SetTitle("Upstream");
  ds->Add(planeDNN_DS_regUS);
  
  planeDNN_DS_regDS->SetFillColor(kOrange+1);
  planeDNN_DS_regDS->SetLineColor(kOrange+1);
  planeDNN_DS_regDS->SetTitle("Downstream");
  ds->Add(planeDNN_DS_regDS);
  
  planeDNN_DS_other->SetFillColor(kMagenta+2);
  planeDNN_DS_other->SetLineColor(kMagenta+2);
  planeDNN_DS_other->SetTitle("Other");
  ds->Add(planeDNN_DS_other);
  
  TCanvas *DSstacked = new TCanvas("DSstacked","stacked hists for Downstream SB",10,10,1200,800);
  DSstacked->SetGrid(1);
  planeDNN_DS_data->Draw("p");   
  ds->Draw("hist same");
  
  planeDNN_DS_data->GetXaxis()->SetTitle( "PlaneDNN"  );
  
  gPad->SetGrid(1);
  gPad->Print(Form("%s/DS_planeDNN_%s.png",outdir.c_str(),trueZ.c_str()),"png");
  gPad->Update();
  
  DSstacked->Update();

      
  //target region
  THStack *tgt = new THStack("tgt","Stacked 1D Target region histograms");
  
  planeDNN_tgt_Fe->SetFillColor(kRed);
  planeDNN_tgt_Fe->SetMarkerColor(kRed);
  planeDNN_US_Fe->SetTitle("Iron");
  ds->Add(planeDNN_tgt_Fe);
  
  planeDNN_tgt_Pb->SetFillColor(kBlue);
  planeDNN_tgt_Pb->SetMarkerColor(kBlue);
  planeDNN_tgt_Pb->SetTitle("Lead");
  ds->Add(planeDNN_tgt_Pb);
  
  planeDNN_tgt_C->SetFillColor(kCyan+2);
  planeDNN_tgt_C->SetMarkerColor(kCyan+2);
  planeDNN_tgt_C->SetTitle("Carbon");
  ds->Add(planeDNN_tgt_C);
  
  planeDNN_tgt_regUS->SetFillColor(kMagenta);
  planeDNN_tgt_regUS->SetMarkerColor(kMagenta);
  planeDNN_tgt_regUS->SetTitle("Upstream");
  ds->Add(planeDNN_tgt_regUS);
  
  planeDNN_tgt_regDS->SetFillColor(kOrange+1);
  planeDNN_tgt_regDS->SetMarkerColor(kOrange+1);
  planeDNN_tgt_regDS->SetTitle("Downstream");
  ds->Add(planeDNN_tgt_regDS);
  
  planeDNN_tgt_other->SetFillColor(kMagenta+2);
  planeDNN_tgt_other->SetMarkerColor(kMagenta+2);
  planeDNN_tgt_other->SetTitle("Other");
  ds->Add(planeDNN_tgt_other);
  
  TCanvas *TGTstacked = new TCanvas("TGTstacked","stacked hists for Target region",10,10, 1200,800);
  TGTstacked->SetGrid(1);
  planeDNN_tgt_data->Draw("p");   
  tgt->Draw("hist same");
  
  planeDNN_tgt_data->GetXaxis()->SetTitle( "PlaneDNN"  );
  
  gPad->SetGrid(1);
  gPad->Print(Form("%s/Tgt_planeDNN_%s.png",outdir.c_str(),trueZ.c_str()),"png");
  gPad->Update();
  
  TGTstacked->Update();

      
  //Legend
  TCanvas *c = new TCanvas("c","Legend",600,500);
  planeDNN_US_data->SetLineColor(kBlack);
  planeDNN_US_Fe->SetLineColor(kRed);
  planeDNN_US_Pb->SetLineColor(kBlue);
  planeDNN_US_C->SetLineColor(kCyan+2);
  planeDNN_US_other->SetLineColor(kMagenta+2);
  planeDNN_US_regUS->SetLineColor(kMagenta);
  planeDNN_US_regDS->SetLineColor(kOrange+1);
  
  TLegend* leg = new TLegend(0.25,0.3,0.8,0.9);
  leg->AddEntry(planeDNN_US_data,"Data", "lp");
  leg->AddEntry(planeDNN_US_C,"Carbon", "f");
  leg->AddEntry(planeDNN_US_Fe,"Iron", "f");
  leg->AddEntry(planeDNN_US_Pb,"Lead", "f");
  leg->AddEntry(planeDNN_US_regUS,"US", "f");
  leg->AddEntry(planeDNN_US_regDS,"DS", "f");
  leg->AddEntry(planeDNN_US_other,"Other", "f");
  leg->Draw();
  gPad->Print(Form("%s/legend.png",outdir.c_str()),"png");
  gPad->Update();
  
  c->Update();
  

  TObjArray hlistUS, hlistDS, hlisttgt;
  //Add hists to lists
  if(targetZ==26)hlistUS.Add(planeDNN_US_Fe );
  else if(targetZ==82)hlistUS.Add(planeDNN_US_Pb );
  else if (targetZ==6) hlistUS.Add(planeDNN_US_C );
  else {
    cout<<"Unknown target"<<endl;
    return 1;}
  hlistUS.Add(planeDNN_US_other );
  hlistUS.Add(planeDNN_US_regUS );
  hlistUS.Add(planeDNN_US_regDS );

  if(targetZ==26)      hlistDS.Add(planeDNN_DS_Fe );
  else if(targetZ==82) hlistDS.Add(planeDNN_DS_Pb );
  else if (targetZ==6) hlistDS.Add(planeDNN_DS_C );
  else {
    cout<<"Unknown target"<<endl;
    return 1;}
  hlistDS.Add(planeDNN_DS_other );
  hlistDS.Add(planeDNN_DS_regUS );
  hlistDS.Add(planeDNN_DS_regDS );

  if(targetZ==26)hlisttgt.Add(planeDNN_tgt_Fe );
  else if(targetZ==82)hlisttgt.Add(planeDNN_tgt_Pb );
  else if (targetZ==6) hlisttgt.Add(planeDNN_tgt_C );
  else {
    cout<<"Unknown target"<<endl;
    return 1;}
  hlisttgt.Add(planeDNN_tgt_other );
  hlisttgt.Add(planeDNN_tgt_regUS );
  hlisttgt.Add(planeDNN_tgt_regDS );


  //=========Stack plots before tuning============

  TString canname = Form( "SidebandUSPosStackPlot_planeDNN_z%d",targetZ ); //File and canvas name("png", "root", "C", "pdf")
  string plotname = "Upstream Region - "+ trueZ; // Histogram title
  DrawMeSomeDataMCStackedHistos( mnvPlotter, planeDNN_US_data, hlistUS, canname, plotname, outdir, "Enu", targetZ, makeDISPlots );


  canname = Form( "SidebandDSPosStackPlot_planeDNN_z%d",targetZ ); //File and canvas name("png", "root", "C", "pdf")
  plotname = "Downstream Region - "+ trueZ; // Histogram title
  DrawMeSomeDataMCStackedHistos( mnvPlotter, planeDNN_DS_data, hlistDS, canname, plotname, outdir, "planeDNN", targetZ, makeDISPlots );

  
  canname = Form( "TargetPosStackPlot_planeDNN_z%d",targetZ ); //File and canvas name("png", "root", "C", "pdf")
  plotname = "Target Region - "+ trueZ; // Histogram title
  DrawMeSomeDataMCStackedHistos( mnvPlotter, planeDNN_tgt_data, hlisttgt, canname, plotname, outdir, "planeDNN", targetZ, makeDISPlots );

  //*********************************************************STARTING TUNING STUFF
  
  vector<MnvH1D*> histos_US, histos_DS,  histos_mc_fid;
  
  histos_US.push_back( planeDNN_US_data );
  if( targetZ == 6 ) histos_US.push_back(planeDNN_US_C );
  if( targetZ == 26 ) histos_US.push_back(planeDNN_US_Fe );
  if( targetZ == 82 ) histos_US.push_back(planeDNN_US_Pb );
  histos_US.push_back(planeDNN_US_other );
  histos_US.push_back(planeDNN_US_regUS );
  histos_US.push_back(planeDNN_US_regDS );
  
  histos_DS.push_back( planeDNN_DS_data );
  if( targetZ == 6 ) histos_DS.push_back(planeDNN_DS_C );
  if( targetZ == 26 ) histos_DS.push_back(planeDNN_DS_Fe );
  if( targetZ == 82 ) histos_DS.push_back(planeDNN_DS_Pb );
  histos_DS.push_back(planeDNN_DS_other );
  histos_DS.push_back(planeDNN_DS_regUS );
  histos_DS.push_back(planeDNN_DS_regDS );
      

      
  /* histos_mc_tuning_US.push_back( planeDNN_US_data );
  if( targetZ == 6 ) histos_mc_tuning_US.push_back(planeDNN_US_C );
  if( targetZ == 26 ) histos_mc_tuning_US.push_back(planeDNN_US_Fe );
  if( targetZ == 82 ) histos_mc_tuning_US.push_back(planeDNN_US_Pb );
  histos_mc_tuning_US.push_back(planeDNN_US_other );
  histos_mc_tuning_US.push_back(planeDNN_US_regUS );
  histos_mc_tuning_US.push_back(planeDNN_US_regDS );
  
  histos_mc_tuning_DS.push_back( planeDNN_DS_data );
  if( targetZ == 6 ) histos_mc_tuning_DS.push_back(planeDNN_DS_C );
  if( targetZ == 26 ) histos_mc_tuning_DS.push_back(planeDNN_DS_Fe );
  if( targetZ == 82 ) histos_mc_tuning_DS.push_back(planeDNN_DS_Pb );
  histos_mc_tuning_DS.push_back(planeDNN_DS_other );
  histos_mc_tuning_DS.push_back(planeDNN_DS_regUS );
  histos_mc_tuning_DS.push_back(planeDNN_DS_regDS );
  */
  
  //All Selected DIS - Upstream Region
  //  histos_US.push_back( planeDNN_US_data );
  //  AddHistosToAVectorUS( histos_US, hists_sideband_US, targetZ );
  //  cout<<"Hist data = "; hists_data_sideband_US->Print("ALL");
  
  //All Selected DIS - Downstream Region
  //  histos_DS.push_back( planeDNN_DS_data );
  //  AddHistosToAVectorDS( histos_DS, hists_sideband_DS, targetZ );
      
  //create a rootfile that will hold all the scale factors
  TFile *scaleFactor= new TFile(Form("%s/CHScaleFactors_z%d_%s_%s.root", indir.c_str(), targetZ, playlist.c_str(), getenv("NUKECC_TAG") ), "read") ;
  
  //=========needs work=======
  //if (tune_fiducial) presumably some other file 
      
      
  // Loop over all variables and make the plots
  //  for( unsigned int ivar = 0; ivar != vars.size(); ++ivar ){
  
  //    TString canname;
  //    if(!tune_fiducial){

  MnvH1D *scaleFactor_US  = (MnvH1D*)scaleFactor->Get("scaleFactor_US_planeDNN");
  MnvH1D *scaleFactor_DS  = (MnvH1D*)scaleFactor->Get("scaleFactor_DS_planeDNN");
      
  //All Selected DIS - Upstream Region
  //====================================
  
  //      histos_mc_tuning_US.push_back(  planeDNN_US_data );
  //      AddHistosToAVectorUS( histos_mc_tuning_US, hists_sideband_US, targetZ );
  
  //     cout<<"Hist data 2 = "; hists_data_sideband_US[ivar]->Print("ALL");
  //pass the vectors of histos for tuning
  string plot_title = "Tuned Upstream Plastic Background -  " + trueZ;
  
  //pass the vectors of histos for tuning
  canname = Form( "AfterTuning_US_%s_%d",  "planeDNN", targetZ );
  DrawTunedBG( mnvPlotter, histos_US, "planeDNN", plot_title, canname, targetZ, scaleFactor_US, false, true, playlist, outdir);
  
  //All Selected DIS - Downstream Region
  //=====================================
  //      histos_mc_tuning_DS.push_back( planeDNN_DS_data );
  //      AddHistosToAVectorDS( histos_mc_tuning_DS, hists_sideband_DS, targetZ );
  
  //pass the vectors of histos for tuning
  plot_title = "Tuned Downstream Plastic Background -  " + trueZ;
  
  //pass the vectors of histos for tuning
  canname = Form( "AfterTuning_DS_%s_%d",  "planeDNN", targetZ );
  DrawTunedBG( mnvPlotter, histos_DS,  "planeDNN", plot_title, canname, targetZ, scaleFactor_DS, true, false, playlist, outdir);
      

  
      /*    } else {
	    
      //All Selected DIS - Downstream Region
      histos_mc_tuning_fid.push_back( hists_data_fiducial[ivar] );
      AddHistosToAVector( histos_mc_tuning_fid, hists_fiducial[ivar], targetZ );
      
      //pass the vectors of histos for background fitting and tuning
      string plot_title = "Tuned Fiducial -  " + trueZ;
      canname = Form( "%s_%s",  vars[ivar].c_str(), suffix.c_str() );
      TuneFiducial( mnvPlotter, histos_mc_tuning_fid, targetZ, vars[ivar], canname, plot_title, dataMCScale, playlist );
      
      //pass the vectors of histos for background fitting and tuning
      plot_title = "BG Subtracted -  " + trueZ;
      canname = Form( "%s_%s",  vars[ivar].c_str(), suffix.c_str() );
      SubtractFiducial( mnvPlotter, histos_mc_tuning_fid, targetZ, vars[ivar], canname, plot_title, dataMCScale, playlist );
      
      //clear vector before going to the next variable
      histos_mc_tuning_fid.clear(); //second index is DS scale
      }
      */
      //  }//loop over all variables
      
      //write scale factor histograms 
      //if( ! tune_fiducial ) scaleFactor->Write();
      //if( ! tune_fiducial ) scaleFactor->Close();
      
      //close text file containing chi^2/ndf values
      if( writeChi2 )fclose(Chi2Table);
      
      //*********************************************************END TUNING STUFF}
      return 0;
}

void DrawMeSomeDataMCStackedHistos( MnvPlotter mnvPlotter, MnvH1D* histo_data, TObjArray hlist, TString canname, string plotname, string dir, string var, int targetZ, const bool makeDISPlots ){

  //first declared some plotting parameters
  mnvPlotter.legend_offset_x = .03;
  mnvPlotter.width_xspace_per_letter = .33;

  cout << "canname ***** " << canname << endl;
  TCanvas can( canname, canname, 1200,800 );
  mnvPlotter.DrawDataStackedMC( histo_data, &hlist, 1, legLoc, "Data", -1, -1, -1 ); //use my colors and styles

  mnvPlotter.WriteNorm("Abs-Normalized", infoX, infoY);
  mnvPlotter.WriteNorm("Stat. Errors Only", infoX, infoY-.04);
  mnvPlotter.MultiPrint( &can, Form("%s/%s",dir.c_str(),can.GetName()) );

}

void  DrawTunedBG/*CalcScaleFactorAndTuneBG*/( MnvPlotter mnvPlotter, vector<MnvH1D*> histos, const string var, string plot_title, TString canname, int targetZ,  MnvH1D* scaleFactorHist, bool fixUS, bool fixDS, string playlist, string outdir ){
  
  //only want to save scale factors histograms
  string channelstr = "Inclusive";
 

  //clone for before tuned mc templates
  MnvH1D* data_sideband_untuned      = (MnvH1D*)histos[0]->Clone("untuned_data");
  MnvH1D* untuned_mc_sideband_signal = (MnvH1D*)histos[1]->Clone("untuned_mc_signal");
  MnvH1D* untuned_mc_sideband_other  = (MnvH1D*)histos[2]->Clone("untuned_mc_other");
  MnvH1D* untuned_mc_sideband_us     = (MnvH1D*)histos[3]->Clone("untuned_mc_us");
  MnvH1D* untuned_mc_sideband_ds     = (MnvH1D*)histos[4]->Clone("untuned_mc_ds");
  //clone the tuned mc templates
  MnvH1D* data_sideband            = (MnvH1D*)histos[0]->Clone("tuned_data");
  MnvH1D* tuned_mc_sideband_signal = (MnvH1D*)histos[1]->Clone("tuned_mc_signal");
  MnvH1D* tuned_mc_sideband_other  = (MnvH1D*)histos[2]->Clone("tuned_mc_other");
  MnvH1D* tuned_mc_sideband_us     = (MnvH1D*)histos[3]->Clone("tuned_mc_us");
  MnvH1D* tuned_mc_sideband_ds     = (MnvH1D*)histos[4]->Clone("tuned_mc_ds");
 
  //untunned scintillator background before adding it to the stacks 
  TObjArray hlistUntuned;
  hlistUntuned.Add(untuned_mc_sideband_signal);
  hlistUntuned.Add(untuned_mc_sideband_other);
  hlistUntuned.Add(untuned_mc_sideband_us);
  hlistUntuned.Add(untuned_mc_sideband_ds);
  
  //---- PLOT DATA/MC RATIO and CHI^2
  string targetString;
  if( targetZ == 6 ) targetString = "t03_z06";
  if( targetZ == 26 ) targetString = "t01_z26";
  if( targetZ == 82 ) targetString = "t01_z82";
  
  string cname = (string)canname;
  cout <<"CNAME = "<<cname<<endl;// Defining name for the plot  "AfterTuning_....".
  string regID;
  if( string::npos != cname.find("_US_") ){
    regID = "US";
  }
  if( string::npos != cname.find("_DS_") ){
    regID = "DS";
  }
  
  TString cName = "BeforeTuning_" + regID + "_" + var + "_" + targetString;// + "_" + NukeUtils::Get().GetPlotSuffix().Data(); 
  cout<<"cName = "<<cName<<endl;
  PlotChi2Stat( mnvPlotter, data_sideband_untuned, hlistUntuned, var, targetZ, cName ); 
 
  TuneScintBackground( data_sideband_untuned, untuned_mc_sideband_signal, untuned_mc_sideband_other, untuned_mc_sideband_us, untuned_mc_sideband_ds, tuned_mc_sideband_signal, tuned_mc_sideband_other, tuned_mc_sideband_us, tuned_mc_sideband_ds, scaleFactorHist,  fixUS, fixDS, channelstr, onlyCalculateCV, writeSysCode ); //tune all universes as well

 
  //tuned the scintillator background before adding it to the stacks 
  TObjArray hlistTuned;
  hlistTuned.Add(tuned_mc_sideband_signal);
  hlistTuned.Add(tuned_mc_sideband_other);
  hlistTuned.Add(tuned_mc_sideband_us);
  hlistTuned.Add(tuned_mc_sideband_ds);
  
  //--- PLOT THE TUNED BACKGROUND --------
  // now let's draw the PLOT!
  //cout<<"Plot tuned background"<<"\t Outdir:"<<outdir.c_str()<<endl;
  DrawMeSomeDataMCStackedHistos( mnvPlotter, data_sideband, hlistTuned, canname, plot_title, outdir, var, targetZ, makeDISPlots ); 
  
  //---- PLOT DATA/MC RATIO and CHI^2
  PlotChi2Stat( mnvPlotter, data_sideband, hlistTuned,  var, targetZ, canname ); 
  
}

void  DrawSomeMCStackedHistos( MnvPlotter mnvPlotter, TObjArray hlist, double scale, TString canname, string plotname, string var, int targetZ, const bool makeDISPlots ){
  
  //first declared some plotting parameters
  mnvPlotter.legend_offset_x = .03;
  mnvPlotter.width_xspace_per_letter = .33;
  
  TCanvas can( canname, canname, 1200,800 );
  mnvPlotter.DrawStackedMC( &hlist, scale, legLoc, -1, -1, -1 ); //use my colors and styles
  
  mnvPlotter.WriteNorm("Abs-Normalized", infoX, infoY);
  mnvPlotter.WriteNorm("Stat. Errors Only", infoX, infoY-.04);
  mnvPlotter.MultiPrint( &can, can.GetName() );

}

void  PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos,  string var, int targetZ, TString canname ){
  
  //first declared some plotting parameters
  mnvPlotter.legend_offset_x = .03;
  mnvPlotter.width_xspace_per_letter = .33;
  
  string targetString;  
  if( targetZ == 6 )
    targetString = "Carbon";
  else if( targetZ == 26 )
    targetString = "Iron";
  if( targetZ == 82 )
    targetString = "Lead";
 
   
  //just get the CV for now 
  MnvH1D* data_cv = (MnvH1D*)dataHisto->Clone("mnvh1d_data_histo");
  cout << "data_cv = " << data_cv->Integral() << endl; 
  //get total MC
  MnvH1D* totalMC = (MnvH1D*)mchistos[0]->Clone("mnvh1d_mc_histo");
  totalMC->SetFillColor(0);
  totalMC->SetFillStyle(0);

  string cname = (string)canname; //"AfterTuning_..." plots
  if( targetZ == 6 && string::npos != cname.find("_US_") && string::npos != cname.find("_planeDNN_") ) totalMC->GetXaxis()->SetRangeUser(22., 28.);
  else if( targetZ == 6 && string::npos != cname.find("_DS_") && string::npos != cname.find("_planeDNN_") ) totalMC->GetXaxis()->SetRangeUser(34., 40.);
  
  for( int i = 1; i < mchistos.GetLast()+1; i++ ){
    MnvH1D* mc_histo_add = (MnvH1D*)mchistos[i]->Clone("mnvh1d_mc_histo_add");
    mc_histo_add->SetFillColor(0);
    mc_histo_add->SetFillStyle(0);
    totalMC->Add(mc_histo_add);
  }
  
//I do not know whats happening here...  
  cout << "total MC = " << totalMC->Integral() << endl; 
  if( data_cv->GetSumw2N() == 0 ) 
    data_cv->Sumw2();
  if( totalMC->GetSumw2N() == 0 ){ 
    totalMC->Sumw2();
  }
  
  int ndfStat = 1;
  double chi2Sys = mnvPlotter.Chi2DataMC( data_cv, totalMC, ndfStat);
  ndfStat -= 1;
  
  TString labelStat = Form("Stat. + Sys. Error, #chi^{2}/ndf = %3.2f/%d = %3.2f", chi2Sys, ndfStat, chi2Sys/(Double_t)ndfStat);
  cout<<"Label stat = "<<labelStat<<endl;
  
  string region = "Fiducial";
  if( string::npos != cname.find("_US_") ) region = "Upstream Plastic Background";
  if( string::npos != cname.find("_DS_") ) region = "Downstream Plastic Background";
  
  string status = "After Tuning";
  if( string::npos != cname.find("BeforeTuning_") ) status = "Before Tuning";

  if( writeChi2 && !tune_fiducial ) fprintf(Chi2Table, "%s %s %s %.2f \n", var.c_str(), region.c_str(), status.c_str(), chi2Sys/(Double_t)ndfStat ); 
  if( writeChi2 && tune_fiducial  ) fprintf(Chi2Table, "%s Fiducial %s, With Sys. %.2f \n", var.c_str(), status.c_str(), chi2Sys/(Double_t)ndfStat ); 

  double plotMin = -1., plotMax = -1.;
  GetMinAndMaxAxisRatio( data_cv, totalMC, plotMin, plotMax );
}

void  GetMinAndMaxAxisRatio( TH1D* histData, TH1D* histMC, double &plotMin, double &plotMax ){
  
  TH1D* ratio = (TH1D*)histData->Clone("ratio_plot");
  TH1D* totalMC = (TH1D*)histMC->Clone("mc_clone");
  //totalMC->Scale(dataMCScale);
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

/*
void  TuneFiducial( MnvPlotter &mnvPlotter, vector<MnvH1D*> histos_mc_tuning, int targetZ, string var, TString canname, string plot_title, double dataMCScale, string playlist  )
{

  //clone the tuned mc templates
  MnvH1D* data_fiducial              = (MnvH1D*)histos_mc_tuning[0]->Clone("tuned_data");
  MnvH1D* untuned_mc_fiducial_signal = (MnvH1D*)histos_mc_tuning[1]->Clone("untuned_mc_signal");
  MnvH1D* untuned_mc_fiducial_other  = (MnvH1D*)histos_mc_tuning[2]->Clone("untuned_mc_other");
  MnvH1D* tuned_mc_fiducial_us       = (MnvH1D*)histos_mc_tuning[3]->Clone("tuned_mc_us");
  MnvH1D* tuned_mc_fiducial_ds       = (MnvH1D*)histos_mc_tuning[4]->Clone("tuned_mc_ds");
  MnvH1D* untuned_mc_fiducial_us     = (MnvH1D*)histos_mc_tuning[3]->Clone("untuned_mc_us");
  MnvH1D* untuned_mc_fiducial_ds     = (MnvH1D*)histos_mc_tuning[4]->Clone("untuned_mc_ds");
  MnvH1D* multiply_mc_fiducial_us     = (MnvH1D*)histos_mc_tuning[3]->Clone("multiply_mc_us");
  MnvH1D* multiply_mc_fiducial_ds     = (MnvH1D*)histos_mc_tuning[4]->Clone("multiply_mc_ds");
  
  //save the untuned into TObjArray
  //save them into TObjArray
  TObjArray hlistUntuned;
  //hlistUntuned.Add(untuned_mc_fiducial_signal);
  hlistUntuned.Add(untuned_mc_fiducial_other);
  hlistUntuned.Add(untuned_mc_fiducial_us);
  hlistUntuned.Add(untuned_mc_fiducial_ds);
  
  //tune the US and DS plastic
  tuned_mc_fiducial_us = TunePlasticBackground( multiply_mc_fiducial_us, "US", targetZ, playlist, var );
  tuned_mc_fiducial_ds = TunePlasticBackground( multiply_mc_fiducial_ds, "DS", targetZ, playlist, var );
 
  //save them into TObjArray
  TObjArray hlistTuned;
  hlistTuned.Add(untuned_mc_fiducial_signal);
  hlistTuned.Add(untuned_mc_fiducial_other);
  hlistTuned.Add(tuned_mc_fiducial_us);
  hlistTuned.Add(tuned_mc_fiducial_ds);
  
  //--- PLOT THE TUNED FIDUCIAL --------
  // now let's draw the PLOTS
  TString cName = "BeforeTuning_Fiducial_" + canname;
  cout << "plotting **** " << cName << endl; 
  DrawMeSomeDataMCStackedHistos( mnvPlotter, data_fiducial, hlistUntuned, dataMCScale, cName, plot_title, outdir, var, targetZ, makeDISPlots ); 
  cName = "AfterTuning_Fiducial_" + canname;
  cout << "plotting ***** " << cName << endl; 
  DrawMeSomeDataMCStackedHistos( mnvPlotter, data_fiducial, hlistTuned, dataMCScale, cName, plot_title, outdir, var, targetZ, makeDISPlots ); 
  //---- PLOT DATA/MC RATIO and CHI^2
  cName = "BeforeTuning_Fiducial_" + canname;
  cout << "plotting " << canname << endl; 
  PlotChi2Stat( mnvPlotter, data_fiducial, hlistUntuned, dataMCScale, var, targetZ, cName );
  cName = "AfterTuning_Fiducial_" + canname; 
  cout << "plotting " << canname << endl; 
  PlotChi2Stat( mnvPlotter, data_fiducial, hlistTuned, dataMCScale, var, targetZ, cName );
  
}

void  SubtractFiducial( MnvPlotter &mnvPlotter, vector<MnvH1D*> histos_mc_tuning, int targetZ, string var, TString canname, string plot_title, double dataMCScale, string playlist  )
{

  //clone the tuned mc templates
  MnvH1D* data_fiducial              = (MnvH1D*)histos_mc_tuning[0]->Clone("tuned_data");
  MnvH1D* untuned_mc_fiducial_signal = (MnvH1D*)histos_mc_tuning[1]->Clone("untuned_mc_signal");
  MnvH1D* untuned_mc_fiducial_other  = (MnvH1D*)histos_mc_tuning[2]->Clone("untuned_mc_other");
  MnvH1D* tuned_mc_fiducial_us       = (MnvH1D*)histos_mc_tuning[3]->Clone("tuned_mc_us");
  MnvH1D* tuned_mc_fiducial_ds       = (MnvH1D*)histos_mc_tuning[4]->Clone("tuned_mc_ds");
  MnvH1D* untuned_mc_fiducial_us     = (MnvH1D*)histos_mc_tuning[3]->Clone("untuned_mc_us");
  MnvH1D* untuned_mc_fiducial_ds     = (MnvH1D*)histos_mc_tuning[4]->Clone("untuned_mc_ds");

  //save the untuned into TObjArray
  //save them into TObjArray
  TObjArray hlistUntuned;
  hlistUntuned.Add(untuned_mc_fiducial_signal);
  hlistUntuned.Add(untuned_mc_fiducial_other);
  hlistUntuned.Add(untuned_mc_fiducial_us);
  hlistUntuned.Add(untuned_mc_fiducial_ds);
  
  //tune the US and DS plastic
  tuned_mc_fiducial_us = TunePlasticBackground( untuned_mc_fiducial_us, "US", targetZ, playlist, var );
  tuned_mc_fiducial_ds = TunePlasticBackground( untuned_mc_fiducial_ds, "DS", targetZ, playlist, var );

  data_fiducial->Add( tuned_mc_fiducial_us, -dataMCScale ); 
  data_fiducial->Add( tuned_mc_fiducial_ds, -dataMCScale );
 
  //save them into TObjArray
  TObjArray hlistTuned;
  hlistTuned.Add(untuned_mc_fiducial_signal);
  hlistTuned.Add(untuned_mc_fiducial_other);
  //hlistTuned.Add(tuned_mc_fiducial_us);
  //hlistTuned.Add(tuned_mc_fiducial_ds);
  
  //--- PLOT THE TUNED FIDUCIAL --------
  // now let's draw the PLOTS
  cout << "cName fiducial = " << canname << endl;
  //TString cName = "BeforeTuning_BGSubtracted_" + canname;
  //DrawMeSomeDataMCStackedHistos( mnvPlotter, data_fiducial, hlistUntuned, dataMCScale, cName, plot_title, outdir, var, targetZ, makeDISPlots ); 
  TString cName = "AfterTuning_BGSubtracted_" + canname;
  cout << "plotting *********** " << cName << endl; 
 
  DrawMeSomeDataMCStackedHistos( mnvPlotter, data_fiducial, hlistTuned, dataMCScale, cName, plot_title, outdir, var, targetZ, makeDISPlots ); 
  //---- PLOT DATA/MC RATIO and CHI^2
  cName = "BeforeTuning_BGSubtracted_" + canname;
  cout << "plotting " << canname << endl; 
  PlotChi2Stat( mnvPlotter, data_fiducial, hlistUntuned, dataMCScale, var, targetZ, cName );
  cName = "AfterTuning_BGSubtracted_" + canname; 
  cout << "plotting " << canname << endl; 
  PlotChi2Stat( mnvPlotter, data_fiducial, hlistTuned, dataMCScale, var, targetZ, cName );
  
}
*/
MnvH1D*  TunePlasticBackground( MnvH1D *histo, string region, int targetZ, string playlist, string var ){
    
    TFile *scaleFactor = new TFile( Form("%s/CHScaleFactors_z%d_%s_%s.root", getenv("FILES"), targetZ, playlist.c_str(), getenv("NUKECC_TAG") ), "read" );
    assert(scaleFactor);
    
    MnvH1D* h_scaleFactors = (MnvH1D*)scaleFactor->Get( Form("scaleFactor_%s_%s", region.c_str(), var.c_str() ) ); //region should be "US" or "DS"
    cout << "name of scale factors = " << scaleFactor->GetName() << endl;
    int nonzerobin = h_scaleFactors->FindFirstBinAbove(0.);
    MnvH1D* tunedHisto = (MnvH1D*)histo->Clone("tunedHisto");
    MnvH1D* scaleHisto = (MnvH1D*)histo->Clone("scaleHisto");
    scaleHisto->Reset();
    cout<<" I cloned the histogram "<<endl;
    double cv_val = h_scaleFactors->GetBinContent(nonzerobin);
    cout<<"I've got the scale factor stuff for SB "<<region.c_str()<<"  the CV value is "<<cv_val<<endl;
    
    for(int i=0;i<scaleHisto->GetNbinsX()+2;i++) scaleHisto->SetBinContent(i,cv_val);
  
    tunedHisto->Multiply(tunedHisto,scaleHisto);
    
    return tunedHisto;
    
}

// ===========================================================
// you're entering background tune zone 
// ============================================================

void  TuneScintBackground(MnvH1D* mnvh1d_data,
			  MnvH1D* mnvh1d_mc_signal,
			  MnvH1D* mnvh1d_mc_other,
			  MnvH1D* mnvh1d_mc_background_US,
			  MnvH1D* mnvh1d_mc_background_DS,
			  MnvH1D* mc_signal,
			  MnvH1D* mc_other,
			  MnvH1D* mc_background_US,
			  MnvH1D* mc_background_DS,
			  MnvH1D* scale_factor,
			  bool fixUS,
			  bool fixDS,
			  string channel /*= "DIS"*/,
			  bool onlyCV /*=false*/,
			  bool writeScaleFactor /*=false*/ )
{
    
    //Get the CV histo
    TH1D* h1d_data_cv = new TH1D(mnvh1d_data->GetCVHistoWithStatError());
    TH1D* h1d_signal_cv = new TH1D(mnvh1d_mc_signal->GetCVHistoWithStatError());
    TH1D* h1d_other_cv = new TH1D(mnvh1d_mc_other->GetCVHistoWithStatError());
    TH1D* h1d_backgroundUS_cv = new TH1D(mnvh1d_mc_background_US->GetCVHistoWithStatError());
    TH1D* h1d_backgroundDS_cv = new TH1D(mnvh1d_mc_background_DS->GetCVHistoWithStatError());
    
    std::cout << "Total data: " << h1d_data_cv->Integral() << std::endl;
    
 
    if( !fixUS && fixDS ){
      ((TH1D*)mc_background_US)->Multiply((TH1D*)scale_factor);
      //((TH1D*)tmp_scale_factor)->Divide( ((TH1D*)tmp_scale_factor), ((TH1D*)tmp_scale_factor), results_cv[0], 1.0);
    }
    else if( !fixDS && fixUS ){
      ((TH1D*) mc_background_DS)->Multiply((TH1D*)scale_factor);
        //((TH1D*)tmp_scale_factor)->Divide( ((TH1D*)tmp_scale_factor), ((TH1D*)tmp_scale_factor), results_cv[0], 1.0);
    }
    
    if( onlyCV ) return;
  
}

double  getChi2( const double * par )
{
    
    double scale = par[0];
    
    // Sideband is bins from plane number 11 to 65
    // Assume no error on the MC, and use the stat errors on the data
    int binLow = m_histo_data->FindBin( 11. ); // same bins on both sets of histograms
    int binHigh = m_histo_data->FindBin( 65. );
    
    double chi2 = 0.0;
    // clone and scale the floating histograms
    TH1D * temp_plastic_float = new TH1D( *(TH1D*)m_histo_plastic_float->Clone("temp_plastic_float") );
    temp_plastic_float->Scale( scale );
    
    double tot_dat = 0.0;
    double tot_sig_fix = 0.0;
    double tot_other_fix = 0.0;
    double tot_plastic_fix = 0.0;
    double tot_plastic_float = 0.0;
    
    for( int b = binLow; b <= binHigh; ++b ) {
        double data = m_histo_data->GetBinContent(b);
        double err = m_histo_data->GetBinError(b);
        double sig_fix = m_histo_sig_fix->GetBinContent(b); // fixed part of MC, should be small
        double other_fix = m_histo_other_fix->GetBinContent(b); // fixed part of MC, should be small
        double plastic_fix = m_histo_plastic_fix->GetBinContent(b); // fixed part of MC, should be small
        double plastic_float = temp_plastic_float->GetBinContent(b); // let the background floats in the fitting
        
        if( data != 0.0 ){
            chi2 += pow( (plastic_float + plastic_fix + other_fix + sig_fix - data)/err, 2 ); // track bin chi2
            tot_dat += data;
            tot_sig_fix += sig_fix;
            tot_other_fix += other_fix;
            tot_plastic_fix += plastic_fix;
            tot_plastic_float += plastic_float;
        }
    }
    
    double unbinned_chi2 = pow( (tot_sig_fix+tot_other_fix+tot_plastic_fix+tot_plastic_float-tot_dat)/sqrt(tot_dat), 2 );
/*    
    bool floatUS = false;
    string histoName = m_histo_plastic_float->GetName();
    if( string::npos != histoName.find("plasticUS") ) floatUS = true;
    
    if( floatUS ){
        cout << "nominal event after fit ( data, signal, other, plastic US, plastic DS ) : " << std::fixed << std::setprecision(2) << tot_dat << " &  "
        << std::fixed << std::setprecision(2) <<  tot_sig_fix << " & "
        << std::fixed << std::setprecision(2) <<  tot_other_fix << " & "
        << std::fixed << std::setprecision(2) <<  tot_plastic_float << " & "
        << std::fixed << std::setprecision(2) <<  tot_plastic_fix << " & " << endl;
        cout << "nominal event after fit ( data, signal, other, plastic US, plastic DS ) : " << std::fixed << tot_dat << " & "
        << std::fixed << std::setprecision(2) <<  tot_sig_fix << " & "
        << std::fixed << std::setprecision(2) <<  tot_other_fix << " & "
        << std::fixed << std::setprecision(2) <<  tot_plastic_fix << " & "
        << std::fixed << std::setprecision(2) <<  tot_plastic_float << " & " << endl;

    }
*/    
    return unbinned_chi2;
    
}
/*
void AddHistosToAVectorUS( vector<MnvH1D*> &v_histos, vector<MnvH1D*> histos, int targetZ ){

  if( targetZ == 6 ) v_histos.push_back(planeDNN_US_C );
  if( targetZ == 26 ) v_histos.push_back(planeDNN_US_Fe );
  if( targetZ == 82 ) v_histos.push_back(planeDNN_US_Pb );
  v_histos.push_back(planeDNN_US_other );
  v_histos.push_back(planeDNN_US_regUS );
  v_histos.push_back(planeDNN_US_regDS );

}
void AddHistosToAVectorDS( vector<MnvH1D*> &v_histos, vector<MnvH1D*> histos, int targetZ ){

  if( targetZ == 6 ) v_histos.push_back(planeDNN_DS_C );
  if( targetZ == 26 ) v_histos.push_back(planeDNN_DS_Fe );
  if( targetZ == 82 ) v_histos.push_back(planeDNN_DS_Pb );
  v_histos.push_back(planeDNN_DS_other );
  v_histos.push_back(planeDNN_DS_regUS );
  v_histos.push_back(planeDNN_DS_regDS );

}
*/
