#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
//#include "syst_common.h"
#include "TArrayD.h"
#include <iostream>
                              
#include "Math/Factory.h"     
#include "Math/Functor.h"     
#include "Math/Minimizer.h"   

#include "RooUnfold/RooUnfoldBayes.h"
#include "RooUnfold/RooUnfoldSvd.h"
#include "RooUnfold/RooUnfoldTUnfold.h"
#include "RooUnfold/RooUnfoldInvert.h"
#include "RooUnfold/RooUnfoldBinByBin.h"
#include "RooUnfold/RooUnfoldResponse.h"
#include "RooUnfold/RooUnfold.h"
#include "MinervaUnfold/MnvResponse.h"
#include "MinervaUnfold/MnvUnfold.h"

#include "PlotUtils/MacroUtil.h" 
//#include "include/CVUniverse_faiza.h"
#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "../include/Variable.h"
//#include "../../include/Variable_faiza.h"
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
#include "TParameter.h"
#include "../include/plotting_functions.h"
#ifndef __CINT__
#endif

using namespace NUKECC_ANA;
using namespace PlotUtils;
using namespace std;

void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name);
void GetMinAndMaxAxisRatio( TH1D* histData, TH1D* histMC, double dataMCScale, double &plotMin, double &plotMax );



//const bool do_fits = true;//Run this code with this boolean as true first and then run it again with boolean as false before plotiing
const bool do_fits = false;//second

TH1D* GetVertErrorBandUniverseHist(MnvH1D* mnvh1d, const std::string& errName, int universe); 
void ConstrainPhysBackground(MnvH1D* mnvh1d_trans_data,
                                   MnvH1D* mnvh1d_mc_signal_in_trans,
                                   MnvH1D* mnvh1d_mc_trans_in_trans,
                                   MnvH1D* mnvh1d_mc_contin_in_trans,
                                   MnvH1D* mnvh1d_contin_data,
                                   MnvH1D* mnvh1d_mc_signal_in_contin,
                                   MnvH1D* mnvh1d_mc_trans_in_contin,
                                   MnvH1D* mnvh1d_mc_contin_in_contin,
                                   MnvH1D* mc_trans_in_trans,
                                   MnvH1D* mc_contin_in_trans,
                                   MnvH1D* mc_trans_in_contin,
                                   MnvH1D* mc_contin_in_contin,
                                   MnvH1D* tmp_scale_factor_trans,
                                   MnvH1D* tmp_scale_factor_contin,
                                   double mc_scale, bool onlyCV /*=false*/, bool writeScaleFactor /*=false*/ );

vector<double> CalcPhysScaleFactorMinimizer( TH1D* h1d_data_trans_cv,
                                                  TH1D* h1d_signal_trans_cv,
                                                  TH1D* h1d_trans_trans_cv,
                                                  TH1D* h1d_contin_trans_cv,
                                                  TH1D* h1d_data_contin_cv,
                                                  TH1D* h1d_signal_contin_cv,
                                                  TH1D* h1d_trans_contin_cv,
                                                  TH1D* h1d_contin_contin_cv,
                                                  double mc_scale,
                                                  const std::string& errName);
double getPhysChi2( const double * par);
      TH1D* m_histo_trans_data;
      TH1D* m_histo_sig_trans;
      TH1D* m_histo_trans_trans;
      TH1D* m_histo_contin_trans;
      TH1D* m_histo_contin_data;
      TH1D* m_histo_sig_contin;
      TH1D* m_histo_trans_contin;
      TH1D* m_histo_contin_contin;

//bool UnfoldHisto( MnvH1D* hUnfolded, MnvH1D *hFolded, const std::string& var, int targetID, int targetZ, string playist, HelicityType::t_HelicityType helicity, bool accCorr /*=false*/, bool addSysErrors /*= true*/, Unfold::t_UnfoldMethod method );
bool UnfoldHisto( MnvH1D* hUnfolded, MnvH1D *hFolded, const std::string& var, string outdir, int targetID, int targetZ, string playist, HelicityType::t_HelicityType helicity, bool accCorr /*=false*/);
 
int UnfoldIterations( const std::string& var );



MnvH1D* TunePhysBackground( MnvH1D *histo, string SB, string outdir, int targetID, int targetZ, string playlist, string var, bool cv_only = false );
MnvH1D* TunePlasticBackground( MnvH1D *tunedHisto, string region, string outdir, int targetID, int targetZ, string playlist, string var );
//MnvH1D* ReBinPhysSBTuning( MnvH1D* h_scaleFactors, string SB, int targetID, int targetZ, string playlist, string var, bool cv_only = false, NukeCC_Binning* binsDef, MnvH1D* h_sb=NULL);
MnvH1D* ReBinPhysSBTuning( MnvH1D* h_scaleFactors, string SB, int targetID, int targetZ, string playlist, string var, bool cv_only = false, MnvH1D* h_sb=NULL);

int main(int argc, char * argv[]){
  //ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if(argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------\
-----"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n"<<
       "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
      "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" <<
      "\t-Playlist\t \t =\t eg. minervame1A"<< std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    return 0;
  }

  string outdir=argv[1];
  int targetID = atoi(argv[2]);
  int targetZ = atoi(argv[3]);
  const string playlist= argv[4];

  TString hname, hname1, hname2;
   if(RunCodeWithSystematics){
  //hname= Form("%s/Hists_PhysicsBackgd_without_SYS_t%d_z%02d_Nu_v1.root", outdir.c_str(), targetID, targetZ);
  hname1= Form("%s/Hists_EventSelection_t%d_z%02d_Nu.root", outdir.c_str(), targetID, targetZ);
  hname2= Form("%s/Hists_Migration_t%d_z%02d_Nu.root", outdir.c_str(), targetID, targetZ);
     }
   else{
  hname1= Form("%s/Hists_EventSelection_t%d_z%02d_Nu.root", outdir.c_str(), targetID, targetZ);
     }
  TFile *f1 = new TFile( hname1,"read" );
  //TFile *f2 = new TFile( hname2,"read" );
  cout<<hname1<<endl;

  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(playlist);
  NukeCC_Binning  *binsDef = new NukeCC_Binning();
 
  MnvPlotter mnvPlotter(kNukeCCStyle);

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

  const string suffix = Form( "t%02d_z%02d_%s", targetID, targetZ, trueZ.c_str() );

  vector<string> vars;
  vars.clear();
  vars.push_back("Emu");
  vars.push_back("Enu");
  vars.push_back("Ehad");
  vars.push_back("x");
  vars.push_back("y");

   HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(playlist);
   //FileType::t_FileType type = utils->GetFileTypeString(std::string fType)
   cout<<"I am here"<<endl;
   
   TFile *fUnfold = new TFile( Form("%s/Unfolded_t%d_z%d_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str() ), "recreate" );
   assert(fUnfold);
  
  TParameter<double> *mcPOT = (TParameter<double>*)f1->Get("MCPOT");
  TParameter<double> *dataPOT = (TParameter<double>*)f1->Get("DataPOT");
  double mcpot = mcPOT->GetVal();
  double datapot = dataPOT->GetVal();

  double dataMCScale = datapot/mcpot;
  cout<<"Data POT = "<<datapot<<"   MC POT = "<<mcpot<<"   Data/MC POT = "<<dataMCScale<<endl;

    // get the hists
    std::vector<MnvH1D*> histsMC, histsData, hists_acccorr, hists_dis;
    std::vector<MnvH1D*> histsMCBG, histsDataBG, mcBG, dataBG, histsMCWrongTarget, histsMCNonDIS;

    vector< MnvH1D*> hists_trans_data, hists_contin_data;
    vector< MnvH1D* > hists_true_signal_in_trans,  hists_true_signal_in_contin;
    vector< MnvH1D* > hists_true_trans_in_trans,   hists_true_trans_in_contin, hists_trans_in_target;
    vector< MnvH1D* > hists_true_contin_in_contin, hists_true_contin_in_trans, hists_contin_in_target;
    vector< MnvH1D* > hists_true_signal_in_signal;
    vector< MnvH1D* > hists_trans_mc_CH_US,  hists_trans_mc_CH_DS;
    vector< MnvH1D* > hists_contin_mc_CH_US, hists_contin_mc_CH_DS;
    vector< MnvH1D* > tmp_scale_factor_trans, tmp_scale_factor_contin;
 
    vector< MnvH1D* >true_trans_in_trans_tuned;
    vector< MnvH1D* >true_contin_in_trans_tuned;
    vector< MnvH1D* >true_trans_in_contin_tuned;
    vector< MnvH1D* >true_contin_in_contin_tuned;
    vector< MnvH1D* >tuned_trans_mc_US, tuned_trans_mc_DS, tuned_contin_mc_US, tuned_contin_mc_DS; 

    vector<MnvH1D*> hists_true_wrong_target_USCH, hists_true_wrong_target_DSCH, hists_true_wrong_target_Other;
    vector<MnvH1D*> hists_true_dis_wrong_target_USCH, hists_true_dis_wrong_target_DSCH, hists_true_dis_wrong_target_Other;
    std::vector<MnvH1D*> HistsMC, HistsData, mcSignals, mcHistsBG, mcHists_nondis, FullSysts;
    std::vector<MnvH1D*> hists_true_sig, hists_trans_SB, hists_contin_SB;
  
   for( vector<string>::iterator i = vars.begin(); i != vars.end(); ++i )
    {
       const string& var = *i;
       f1->cd();

        //Get the sample
        MnvH1D *mnv1  = (MnvH1D*) f1->Get( Form("h_mc_%s",  var.c_str()));
        
        cout <<FileType::kAny<< " sample integral = " << mnv1->Integral() << endl;
        if( !mnv1 )
            cout << "ERROR: could not get MC MnvH1D for variable: " << var << endl;
        
        MnvH1D *mnv2  = (MnvH1D*) f1->Get( Form("h_data_%s",  var.c_str()));
        
        cout <<FileType::kAny<< " sample integral = " << mnv2->Integral() << endl;
        if( !mnv2 )
            cout << "ERROR: could not get MC MnvH1D for variable: " << var << endl;
        
        HistsMC.push_back( mnv1 );
        HistsData.push_back( mnv2 );
        
        MnvH1D *fullsyst = (MnvH1D*)f1->Get( Form("h_mc_%s",  var.c_str()));; //Hist we subtract true DIS BG and extracted CH BG from
        FullSysts.push_back( fullsyst );
       
         cout<<"15ii9"<<endl;
        MnvH1D *true_signal = (MnvH1D*)f1->Get( Form("true_in_target_%s",  var.c_str()));; //Hist we subtract true DIS BG and extracted CH BG from
        mcSignals.push_back(true_signal);
        // ===================================================
        //           PHYSICS BACKGROUND
        // ===================================================
        // Signal
        // Transition
        MnvH1D *mnv_true_trans_bg_target  = (MnvH1D*) f1->Get( Form("true_trans_in_target_%s",  var.c_str()));
        hists_trans_SB.push_back( mnv_true_trans_bg_target );
        // Continum 
        MnvH1D *mnv_true_contin_bg_target = (MnvH1D*) f1->Get( Form("true_contin_in_target_%s", var.c_str()));
        hists_contin_SB.push_back( mnv_true_contin_bg_target );
        
        // ===================================================
        //           PLASTIC BACKGROUND
        // ===================================================
        // --- US Plastic
        MnvH1D *hists_wrong_target_USCH  = (MnvH1D*) f1->Get( Form("h_USCH_%s",  var.c_str()));
        hists_true_wrong_target_USCH.push_back(hists_wrong_target_USCH);
        // --- DS Plastic
        MnvH1D *hists_wrong_target_DSCH  = (MnvH1D*) f1->Get( Form("h_DSCH_%s",  var.c_str()));
        hists_true_wrong_target_DSCH.push_back(hists_wrong_target_DSCH);
        // --- other
        MnvH1D *hists_wrong_target_other  = (MnvH1D*) f1->Get( Form("h_other_%s",  var.c_str()));
        hists_true_wrong_target_Other.push_back(hists_wrong_target_other);
        
    }


    for( unsigned int ivar = 0; ivar != vars.size(); ++ivar )
    {
        cout<<"199 "<<ivar<<endl;
        f1->cd();
        const string var = vars[ivar];
        //TString hname = utils->GetHistName( "signal", FileType::kAny, var, targetID, targetZ );
        //MnvH1D *mnvD  = dynamic_cast<MnvH1D*> (f1->Get( Form("h_mc_%s",  var.c_str())));
        MnvH1D *mnvD  = dynamic_cast<MnvH1D*>( HistsMC[ivar]->Clone("signal")) ;
        cout<<"what is mnvD"<<mnvD<<endl;
        //Add Missing Error Bands
        //if(type==FileType::kData)
        if(FileType::kAny == FileType::kData)
            mnvD->AddMissingErrorBandsAndFillWithCV( *hists_true_wrong_target_USCH[ivar] );
        //======================================================================
        // Wrong Target Background
        //======================================================================
        
        MnvH1D *mnvMCtrue_trans_BG = dynamic_cast<MnvH1D*> (hists_trans_SB[ivar]->Clone("tmptransSB"));
        MnvH1D *mnvMCtrue_contin_BG = dynamic_cast<MnvH1D*> (hists_contin_SB[ivar]->Clone("tmpcontinSB"));
        
        cout<<"You got kidin me"<<endl;    
        MnvH1D *mnvMCtrue_trans_BG_tuned = dynamic_cast<MnvH1D*> (hists_trans_SB[ivar]->Clone("tmptransSBTuned"));
        MnvH1D *mnvMCtrue_contin_BG_tuned = dynamic_cast<MnvH1D*> (hists_contin_SB[ivar]->Clone("tmpcontinSBTuned"));

        cout<<"Going to do physics tuning now"<<endl;
        mnvMCtrue_trans_BG_tuned  = TunePhysBackground( mnvMCtrue_trans_BG,  "trans",  outdir.c_str(), targetID, targetZ, playlist, var);
        cout<<"I'm done tuning the Transition SB"<<endl;
        mnvMCtrue_contin_BG_tuned = TunePhysBackground( mnvMCtrue_contin_BG, "contin", outdir.c_str(), targetID, targetZ, playlist, var);

        cout<<"I'm done tuning the physics background"<<endl;
        
        //======================================================================
        // Wrong Target Background
        //============= Inclusive ==============================================
        MnvH1D *mnvMCtrue_wrong_target_USCH = dynamic_cast<MnvH1D*> (hists_true_wrong_target_USCH[ivar]->Clone("tmpTrueWrongTargetUSCH"));
        MnvH1D *mnvMCtrue_wrong_target_DSCH = dynamic_cast<MnvH1D*> (hists_true_wrong_target_DSCH[ivar]->Clone("tmpTrueWrongTargetDSCH"));
        MnvH1D *mnvMCtrue_wrong_target_Other = dynamic_cast<MnvH1D*> (hists_true_wrong_target_Other[ivar]->Clone("tmpTrueWrongTargetOther"));
        //======================================================================
        // Tuned Wrong Target Background (FOR DATA)
        //============= Inclusive ==============================================
        MnvH1D *mnvMCtrue_wrong_target_USCH_tuned = dynamic_cast<MnvH1D*> (hists_true_wrong_target_USCH[ivar]->Clone("tmpTrueWrongTargetUSCHTuned"));
        MnvH1D *mnvMCtrue_wrong_target_DSCH_tuned = dynamic_cast<MnvH1D*> (hists_true_wrong_target_DSCH[ivar]->Clone("tmpTrueWrongTargetDSCHTuned"));
        MnvH1D *mnvMCtrue_wrong_target_Other_tuned = dynamic_cast<MnvH1D*> (hists_true_wrong_target_Other[ivar]->Clone("tmpTrueWrongTargetOtherTuned")); // this is not tuned using any scale factor because it's a small background,  but I have it here just for consistency
        //=====================================================================
        
        // ==============  TUNE THE MC BACKGROUND! ============================
        // Tune the US and DS plastic based on the sideband constraint
        if(targetID<10){
            mnvMCtrue_wrong_target_USCH_tuned = TunePlasticBackground( mnvMCtrue_wrong_target_USCH, "US", outdir.c_str(), targetID, targetZ, playlist, var ); //first index is USCH scale
            mnvMCtrue_wrong_target_DSCH_tuned = TunePlasticBackground( mnvMCtrue_wrong_target_DSCH, "DS", outdir.c_str(), targetID, targetZ, playlist, var ); //first index is DSCH scale
        }
        
        MnvH1D* true_kinematic_bg_in_target_Data;
        MnvH1D* true_kinematic_bg_in_target_MC;
        MnvH1D* true_dis_wrong_target_Data;
        MnvH1D*  true_dis_wrong_target_MC;
        
        //if(FileType::kAny==FileType::kData){
            true_kinematic_bg_in_target_Data = (MnvH1D*) mnvMCtrue_trans_BG_tuned->Clone("total_kin_bg");
            true_kinematic_bg_in_target_Data->Add(mnvMCtrue_contin_BG_tuned);
             cout<<"I summed my Data backgrounds!"<<true_kinematic_bg_in_target_Data->Integral()<<endl;
            
            true_dis_wrong_target_Data = (MnvH1D*)mnvMCtrue_wrong_target_USCH_tuned->Clone("total_plastic_bg");
            true_dis_wrong_target_Data->Add(mnvMCtrue_wrong_target_DSCH_tuned);
            true_dis_wrong_target_Data->Add(mnvMCtrue_wrong_target_Other_tuned);
            cout<<"I got the Data tuned background histograms "<<true_dis_wrong_target_Data->Integral()<<endl;

        //}
        //if(FileType::kAny==FileType::kMC || FileType::kAny==FileType::kNukeOnlyMC){
            true_kinematic_bg_in_target_MC = (MnvH1D*) mnvMCtrue_trans_BG->Clone("total_kin_bg");
            true_kinematic_bg_in_target_MC->Add(mnvMCtrue_contin_BG);
            cout<<"I summed my MC backgrounds!"<<true_kinematic_bg_in_target_MC->Integral()<<endl;
            
            true_dis_wrong_target_MC = (MnvH1D*)mnvMCtrue_wrong_target_USCH->Clone("total_plastic_bg");
            true_dis_wrong_target_MC->Add(mnvMCtrue_wrong_target_DSCH);
            true_dis_wrong_target_MC->Add(mnvMCtrue_wrong_target_Other);
            cout<<"I got the MC untuned background histograms "<<true_dis_wrong_target_MC->Integral()<<endl;
            
        MnvH1D *mcNonDIS = dynamic_cast<MnvH1D*>(true_kinematic_bg_in_target_MC->Clone("mc_Non_DIS"));        
        MnvH1D *mcWrongTarget = dynamic_cast<MnvH1D*>(true_dis_wrong_target_MC->Clone("mc_wrong_target"));        
       // }
        
        if( !mnvD )
            cout << "ERROR: could not get energy MnvH1D for variable: " << var << endl;
        
        cout << "how many sample? type? " << mnvD->Integral() << ", " << FileType::kAny << endl;
        
       // if(FileType::kAny==FileType::kData){
            //mnvD->AddMissingErrorBandsAndFillWithCV(*true_kinematic_bg_in_target);
            //mnvD->AddMissingErrorBandsAndFillWithCV( *FullSysts[ivar] );
       // }
        
        //Subtract all non-signal events   
        if( !true_kinematic_bg_in_target_Data )
            cout << "ERROR: could not get MC MnvH1D for variable: " << var << endl;
        
   
        ////TString hName = utils->GetHistName( "signal", FileType::kAny, var, targetID, targetZ );
        //TString s1  = dynamic_cast<MnvH1D*> (f1->Get( Form("h_mc_%s",  var.c_str())));
        //MnvH1D *signal  = dynamic_cast<MnvH1D*> (f1->Get( Form("h_mc_%s",  var.c_str())));
        MnvH1D *signalMC  = dynamic_cast<MnvH1D*>( HistsMC[ivar]->Clone("signal")) ;
        MnvH1D *signalData  = dynamic_cast<MnvH1D*>( HistsData[ivar]->Clone("signal")) ;
        //MnvH1D *signal = (MnvH1D*)mnvD->Clone(hName);
        cout<<" The scale is "<<dataMCScale<<endl;
        ////cout << "now how many physics bg ? type? " << true_kinematic_bg_in_target->Integral() << ", " << FileType::kAny << endl;
        ////cout << "now how many plastic bg ? type? " << true_dis_wrong_target->Integral() << ", " << FileType::kAny << endl;
        
        signalMC->Add( true_kinematic_bg_in_target_MC, -1.0 );
        signalMC->Add( true_dis_wrong_target_MC, -1.0 );
        

        MnvH1D * mcBG = dynamic_cast<MnvH1D*>(true_kinematic_bg_in_target_MC->Clone("mc_background"));
        mcBG->Add(true_dis_wrong_target_MC, 1.0 );



 
        signalData->Add( true_kinematic_bg_in_target_Data, -1.0 * dataMCScale );
        signalData->Add( true_dis_wrong_target_Data, -1.0 * dataMCScale );
        
        MnvH1D * dataBG = dynamic_cast<MnvH1D*>(true_kinematic_bg_in_target_Data->Clone("data_background"));
        dataBG->Add(true_dis_wrong_target_Data, 1.0);
        dataBG->Scale(dataMCScale);
        

       
       /* std::vector<std::string> signal_vertNames = signal->GetVertErrorBandNames();
        cout<<var<<" The Signal Error Bands "<<signal_vertNames.size()<<endl;
        
        cout << "now how many sample? type? " << mnvD->Integral() << ", " << FileType::kAny << endl;
        cout << "now how many signal? type? " << signal->Integral() << ", " << FileType::kAny << endl;
        */for( int bin=1; bin < mnvD->GetNbinsX()+1; bin++ ){
        cout << "bin content, bin error? type? " << bin << ", " << mnvD->GetBinContent(bin) << ", " << mnvD->GetBinError(bin) << ", " << FileType::kAny << endl;
        cout << "bin content, bin error? type? " << bin << ", " << signalMC->GetBinContent(bin) << ", " << signalMC->GetBinError(bin) << ", " << FileType::kAny << endl;
        cout << "bin content, bin error? type? " << bin << ", " << signalData->GetBinContent(bin) << ", " << signalData->GetBinError(bin) << ", " << FileType::kAny << endl;
        }
        histsMC.push_back( signalMC );
        histsData.push_back( signalData );
        histsDataBG.push_back(dataBG);
        histsMCBG.push_back(mcBG);  
        histsMCWrongTarget.push_back(mcWrongTarget);  
        histsMCNonDIS.push_back(mcNonDIS);  
        
    }
    cout << "Create Unfolded" << endl;
    //Create unfolded histograms
    std::vector<MnvH1D*> histsUnfold, histsUnfoldData, histsUnfold_accCorr, histsUnfold_dis, histsUnfoldSig;
//    fUnfold->cd();


    int j=0;
    for( unsigned int ivar = 0; ivar != vars.size(); ++ivar )
    {
        
        const string var = vars[ivar];
        
        vector<double> binsLowEdge = binsDef->GetBins( var );
        
        ///TString hnameUnfold  = utils->GetHistName( "unfolded_signal",  FileType::kAny, var, targetID, targetZ );
        TString hnameUnfold  = utils->GetHistName( "unfolded_signal",  FileType::kAny, var);
        //TString hnameUnfold  = utils->GetHistName( "unfolded",  FileType::kAny, var);
        cout<<"Unfolded hist name "<<hnameUnfold<<endl;
        MnvH1D *mnv = new MnvH1D( hnameUnfold, hnameUnfold, binsLowEdge.size()-1, &binsLowEdge[0] );
        
        //TString hnameUnfoldData  = utils->GetHistName( "unfolded_data",  FileType::kAny, var, targetID, targetZ );
        TString hnameUnfoldData  = utils->GetHistName( "unfolded_data",  FileType::kAny, var);
        //TString hnameUnfold  = utils->GetHistName( "unfolded",  FileType::kAny, var);
        cout<<"Unfolded hist name "<<hnameUnfoldData<<endl;
        MnvH1D *mnv1 = new MnvH1D( hnameUnfoldData, hnameUnfoldData, binsLowEdge.size()-1, &binsLowEdge[0] );
        ////mnv->AddMissingErrorBandsAndFillWithCV(*hists[ivar]);
        ////cout<<"Made the histo "<<endl;
        ////std::vector<std::string> unfold_vertNames = mnv->GetVertErrorBandNames();
        ////cout<<var<<" The Unfold Error Bands "<<unfold_vertNames.size()<<endl;
        
        
        ////if( !mnv )
        ////    cout << "ERROR: could not create Unfold MnvH1D for variable: " << var << endl;
        histsUnfold.push_back( mnv );
        histsUnfoldData.push_back( mnv1 );
        
        hnameUnfold  = utils->GetHistName( "unfolded_signal",  FileType::kAny, var, targetID, targetZ );
        //hnameUnfold  = utils->GetHistName( "unfolded_signal",  FileType::kAny, var );
        cout<<"Unfolded hist name "<<hnameUnfold<<endl;
        MnvH1D *truesignal = new MnvH1D( hnameUnfold, hnameUnfold, binsLowEdge.size()-1, &binsLowEdge[0] );
        ////truesignal->AddMissingErrorBandsAndFillWithCV(*hists[ivar]);
        cout<<"Made the histo "<<endl;
        ////std::vector<std::string> unfold_vertNames_truesig = truesignal->GetVertErrorBandNames();
        ////cout<<var<<" The Unfold Error Bands "<<unfold_vertNames_truesig.size()<<endl;
        
        
        ////if( !truesignal )
        ////    cout << "ERROR: could not create Unfold MnvH1D for variable: " << var << endl;
        histsUnfoldSig.push_back( truesignal );
        
        j++;
    }
    
    
    //there are fewer unfold hists than hists
    int unfoldIdx = 0;



    
    for( unsigned int i = 0; i != vars.size(); ++i )
    {
        cout << "****************** what variables do I have? " << vars[i] << endl;    
        //set the number of iterations (if you don't want the code to use the defaults)
        //if(vars[i]=="Enu") utils->SetUnfoldIter( 5 );
	//if(vars[i]=="x") utils->SetUnfoldIter( 8 );
        cout<<"I set my iterations "<<endl;
        
        const string var = vars[i];
        //if( ! NukeUtils::Get().UnfoldVar( var ) )
            //continue;
        cout<<"This is a variable to unfold  "<<endl;
        MnvH1D *mnvSig = histsMC[i];
        MnvH1D *mnvSigData = histsData[i];
        MnvH1D *mnvTrueSig = mcSignals[i];
        MnvH1D *mnvUnfold = histsUnfold[i];
        MnvH1D *mnvUnfoldData = histsUnfoldData[i];
        MnvH1D *mnvUnfoldTrueSig = histsUnfoldSig[i];
         
        MnvH1D *mnvBGData = histsDataBG[i];
        MnvH1D *mnvBGMC   = histsMCBG[i];
        MnvH1D *mnvMCWrongTarget   = histsMCWrongTarget[i];
        MnvH1D *mnvMCNonDIS   = histsMCNonDIS[i];

         
        TString histSignal  = "signal_"+ var;
        TString histData  = "data_"+ var;
        TString BGData  = "BG_data_"+ var;
        TString BGMC  = "BG_mc_"+ var;
        TString MCWrongTarget  = "BG_WrongTarget_"+ var;
        TString MCNonDIS  = "BG_NonDIS_"+ var;
        ////if( !mnvSig || !mnvUnfold || !mnvTrueSig || !mnvUnfoldTrueSig )
        ////    continue;
        cout<<"I have my histograms "<<endl;
        
        cout << "before unfolding now how many signal? type? " << mnvSig->Integral() << ", " << FileType::kAny << endl;
        cout << "before unfolding now how many unfolded? type? " << mnvUnfold->Integral() << ", " << FileType::kAny << endl;


        cout<<"Please print Data before unfolding"<<endl;


        cout << "before unfolding now how many Data? type? " << mnvSigData->Integral() << ", " << FileType::kAny << endl;
        cout << "before unfolding now how many Data? type? " << mnvUnfoldData->Integral() << ", " << FileType::kAny << endl;
         
        // Do the unfolding
        //        bool accCorr = false;
        //bool unfoldWorked = UnfoldHisto( mnvUnfold, mnvSig, var, targetID, targetZ, helicity, false, addSys );
        bool unfoldWorkedMC = UnfoldHisto( mnvUnfold, mnvSig, var, outdir.c_str(), targetID, targetZ, playlist, helicity, false );
        if( !unfoldWorkedMC )
            Warning( "UnfoldHists", Form("  Could not unfold for variable : %s", var.c_str() ) );
        
        bool unfoldWorkedData = UnfoldHisto( mnvUnfoldData, mnvSigData, var, outdir.c_str(), targetID, targetZ, playlist, helicity, false );
        if( !unfoldWorkedMC )
            Warning( "UnfoldHists", Form("  Could not unfold for variable : %s", var.c_str() ) );
        //unfoldWorked = UnfoldHisto( mnvUnfoldTrueSig, mnvTrueSig, var, targetID, targetZ, helicity, false, addSys );
        ////unfoldWorked = UnfoldHisto( mnvUnfoldTrueSig, mnvTrueSig, var, outdir.c_str(), targetID, targetZ, playlist, helicity, false );
        ///if( !unfoldWorked )
        ////    Warning( "UnfoldHists true signal", Form("  Could not unfold for variable : %s", var.c_str() ) );
      
        cout << "after unfolding now how many unfolded? type? " << mnvUnfold->Integral() << ", " << FileType::kAny << endl;
        cout << "after unfolding now how many unfoldedData? type? " << mnvUnfoldData->Integral() << ", " << FileType::kAny << endl;
        
        ////cout << "after unfolding now how many signal? type?   " << mnvTrueSig->Integral() << ", " << FileType::kAny << endl;
        cout << "after unfolding now how many unfolded? type? " << mnvUnfoldTrueSig->Integral() << ", " << FileType::kAny << endl;
	//set the norm bin width to be something more reasonable instead of the 2.5 GeV in Enu, and 0.045  in x that are the default.  
	// ====Note this will cause bigger differences with Anne's results which say they are normalized per 1 GeV, but are really normalized per 2.5 GeV====
	if(var=="Enu" || var=="Emu" || var=="Q2" || var=="W" || var=="ThetaMu"){
	  ////mnvTrueSig->SetNormBinWidth(1.0);
	  mnvSig->SetNormBinWidth(1.0);
	  mnvUnfold->SetNormBinWidth(1.0);
	  mnvUnfoldTrueSig->SetNormBinWidth(1.0);
	  
          mnvSigData->SetNormBinWidth(1.0);
	  mnvUnfoldData->SetNormBinWidth(1.0);
	} 

	if(var=="x" || var=="y"){
	  ///mnvTrueSig->SetNormBinWidth(0.1);
	  mnvSig->SetNormBinWidth(0.1);
	  mnvUnfold->SetNormBinWidth(0.1);
	  mnvUnfoldTrueSig->SetNormBinWidth(0.1);
	  
          mnvSigData->SetNormBinWidth(0.1);
	  mnvUnfoldData->SetNormBinWidth(0.1);
	} 

        //save memory
        TH1::AddDirectory( kTRUE );
        fUnfold->cd();
        ///mnvTrueSig->Write();
        mnvSig->Write(histSignal);
        mnvUnfold->Write();
       
        mnvBGData->Write(BGData);
        mnvBGMC->Write(BGMC);       
        mnvMCWrongTarget->Write(MCWrongTarget);       
        mnvMCNonDIS->Write(MCNonDIS);       
        
        mnvSigData->Write(histData);
        mnvUnfoldData->Write();
        //mnvUnfoldTrueSig->Write();
        delete mnvUnfold;
        
        cout<<"I deleted my unfolded histogram "<<endl;
        
        ++unfoldIdx;
    } //end loop over vars
    
    cout<<"274"<<endl;
   
    
    //cleanup
    fUnfold->Write();
    delete fUnfold;
    
    delete f1;
    //if(fBG)
    //delete fBG;
    
    return 0;

}


MnvH1D* TunePhysBackground( MnvH1D *histo, string SB, string outdir, int targetID, int targetZ, string playlist, string var, bool cv_only ){
    cout<<"I'm tuning the Physics background now "<<endl;
    int temp_tgtZ = 0;
    if(targetID > 10)
        temp_tgtZ = 99;
    if(targetID< 10)
        temp_tgtZ = targetZ;
    
    TFile *scaleFactor = new TFile( Form("%s/PhysBGScaleFactors_t%d_z%d_%s.root", outdir.c_str(), targetID, temp_tgtZ, playlist.c_str() ), "read" );
    cout<<"I'm going to use the following file "<<Form("%s/PhysBGScaleFactor_z%d_%s.root", outdir.c_str(), temp_tgtZ, playlist.c_str() )<<endl;
    assert(scaleFactor);
    cout<<"Got the scale factor file"<<endl;
    scaleFactor->ls();
    
    MnvH1D* h_scaleFactors = (MnvH1D*)scaleFactor->Get( Form("scale_factor_%s_%s", SB.c_str(), var.c_str() ) ); 
    //MnvH1D* h_scaleFactors = (MnvH1D*)scaleFactor->Get( Form("scale_factor_%s_Emu", SB.c_str()) ); 
    cout<<"The scale factor histogram exists"<<endl;
    cout << "name of scale factors = " << scaleFactor->GetName() << endl;
    int nonzerobin = h_scaleFactors->FindFirstBinAbove(0.);
    MnvH1D* tunedHisto = (MnvH1D*)histo->Clone("tunedHisto");
    MnvH1D* scaleHisto = (MnvH1D*)histo->Clone("scaleHisto");
    scaleHisto->Reset();
    cout<<" I cloned the histogram "<<endl;
    double cv_val = h_scaleFactors->GetBinContent(nonzerobin);
    cout<<"I've got the scale factor stuff for SB (it is here) "<<SB.c_str()<<"  the CV value is "<<cv_val<<endl;
    
    for(int i=0;i<scaleHisto->GetNbinsX()+2;i++) scaleHisto->SetBinContent(i,cv_val);
////////new addition, you have to use this as this was used before
  /*  std::vector<std::string> vertNames = histo->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k ){
        MnvVertErrorBand *errBand = h_scaleFactors->GetVertErrorBand( vertNames[k] );
        const int universes = errBand->GetNHists();
        std::vector<TH1D*> vert_hists;
        for(int u=0;u<universes;++u){
            TH1D* tmp_scale = new TH1D(*errBand->GetHist( u ));
            TH1D* tmp_template = new TH1D(scaleHisto->GetCVHistoWithStatError());
            tmp_template->SetName(Form("Scale_universe_%d",u));
            for(int i=0;i<scaleHisto->GetNbinsX()+2;i++) tmp_template->SetBinContent(i,tmp_scale->GetBinContent(nonzerobin));
            vert_hists.push_back(tmp_template);
        }
        scaleHisto->AddVertErrorBand( vertNames[k],vert_hists);
    }*/
/////////////////////////////////////////////////////////////////////
    tunedHisto->Multiply(tunedHisto,scaleHisto);
    
    return tunedHisto;
  
}

// ===========================================================
// you're entering background constraint zone 
// ============================================================

void ConstrainPhysBackground(MnvH1D* mnvh1d_trans_data,
                                        MnvH1D* mnvh1d_mc_signal_in_trans,
                                        MnvH1D* mnvh1d_mc_trans_in_trans,
                                        MnvH1D* mnvh1d_mc_contin_in_trans,
                                        MnvH1D* mnvh1d_contin_data,
                                        MnvH1D* mnvh1d_mc_signal_in_contin,
                                        MnvH1D* mnvh1d_mc_trans_in_contin,
                                        MnvH1D* mnvh1d_mc_contin_in_contin,
                                        MnvH1D* mc_trans_in_trans,
                                        MnvH1D* mc_contin_in_trans,
                                        MnvH1D* mc_trans_in_contin,
                                        MnvH1D* mc_contin_in_contin,
                                        MnvH1D* tmp_scale_factor_trans,
                                        MnvH1D* tmp_scale_factor_contin,
                                        double mc_scale, bool onlyCV /*=false*/, bool writeScaleFactor /*=false*/ )
{
    std::cout<<"I'm going to constrain the physics background"<<std::endl;
    
    mnvh1d_mc_signal_in_trans->Scale(mc_scale);
    mnvh1d_mc_trans_in_trans->Scale(mc_scale);
    mnvh1d_mc_contin_in_trans->Scale(mc_scale);
    mnvh1d_mc_signal_in_contin->Scale(mc_scale);
    mnvh1d_mc_trans_in_contin->Scale(mc_scale);
    mnvh1d_mc_contin_in_contin->Scale(mc_scale);
    
    mc_trans_in_trans->Scale(mc_scale);
    mc_contin_in_trans->Scale(mc_scale);
    mc_trans_in_contin->Scale(mc_scale);
    mc_contin_in_contin->Scale(mc_scale);
    
    //Get the CV histo
    TH1D* h1d_data_trans_cv = new TH1D(mnvh1d_trans_data->GetCVHistoWithStatError());
    TH1D* h1d_signal_trans_cv = new TH1D(mnvh1d_mc_signal_in_trans->GetCVHistoWithStatError());
    TH1D* h1d_trans_trans_cv = new TH1D(mnvh1d_mc_trans_in_trans->GetCVHistoWithStatError());
    TH1D* h1d_contin_trans_cv = new TH1D(mnvh1d_mc_contin_in_trans->GetCVHistoWithStatError());
    TH1D* h1d_data_contin_cv = new TH1D(mnvh1d_contin_data->GetCVHistoWithStatError());
    TH1D* h1d_signal_contin_cv = new TH1D(mnvh1d_mc_signal_in_contin->GetCVHistoWithStatError());
    TH1D* h1d_trans_contin_cv = new TH1D(mnvh1d_mc_trans_in_contin->GetCVHistoWithStatError());
    TH1D* h1d_contin_contin_cv = new TH1D(mnvh1d_mc_contin_in_contin->GetCVHistoWithStatError());
    
    
    
    std::cout << "Total Trans data:  " << h1d_data_trans_cv->Integral() << std::endl;
    std::cout << "Total Contin data: " << h1d_data_contin_cv->Integral() << std::endl;
    
    vector<double> scales_cv, results_cv;
    // fit tells how much the MC should be scaled to match data
    results_cv = CalcPhysScaleFactorMinimizer(h1d_data_trans_cv,
                                              h1d_signal_trans_cv,
                                              h1d_trans_trans_cv,
                                              h1d_contin_trans_cv,
                                              h1d_data_contin_cv,
                                              h1d_signal_contin_cv,
                                              h1d_trans_contin_cv,
                                              h1d_contin_contin_cv,
                                              mc_scale,
                                              "CV");
    
    //scale only the cv here, universe histograms below
    //scale only the floated histogram;
   
    ((TH1D*)mc_trans_in_trans)->Scale(results_cv[0]);
    ((TH1D*) mc_trans_in_contin)->Scale(results_cv[0]);
    ((TH1D*)tmp_scale_factor_trans)->Divide( ((TH1D*)tmp_scale_factor_trans), ((TH1D*)tmp_scale_factor_trans), results_cv[0], 1.0);
    
   
    ((TH1D*)mc_contin_in_trans)->Scale(results_cv[1]);
    ((TH1D*)mc_contin_in_contin)->Scale(results_cv[1]);
    ((TH1D*)tmp_scale_factor_contin)->Divide( ((TH1D*)tmp_scale_factor_contin), ((TH1D*)tmp_scale_factor_contin), results_cv[1], 1.0);
    
    if( onlyCV ) return;
 
    // vertical error band
    std::vector<std::string> vertNames = mc_trans_in_trans->GetVertErrorBandNames();
    for (std::vector<std::string>::iterator name = vertNames.begin();
         name != vertNames.end(); ++name) {
        std::string errName = *name;
        const unsigned int n_universe = mc_trans_in_trans->GetVertErrorBand(errName)->GetNHists();
        
        //Tune CV in each universe
        if( mc_trans_in_trans->HasVertErrorBand(errName)) {
            MnvVertErrorBand *vertError = mc_trans_in_trans->GetVertErrorBand(errName);
            vertError->Scale(results_cv[0],"",false);
            MnvVertErrorBand *scaleVertError = tmp_scale_factor_trans->GetVertErrorBand(errName);
            scaleVertError->TH1D::Divide( ((TH1D*)tmp_scale_factor_trans), ((TH1D*)tmp_scale_factor_trans), results_cv[0], 1.0);
        }
        if( mc_contin_in_trans->HasVertErrorBand(errName)) {
            MnvVertErrorBand *vertError = mc_contin_in_trans->GetVertErrorBand(errName);
            vertError->Scale(results_cv[1],"",false);
            MnvVertErrorBand *scaleVertError = tmp_scale_factor_contin->GetVertErrorBand(errName);
            scaleVertError->TH1D::Divide( ((TH1D*)tmp_scale_factor_contin), ((TH1D*)tmp_scale_factor_contin), results_cv[1], 1.0);
        }
        if( mc_trans_in_contin->HasVertErrorBand(errName)) {
            MnvVertErrorBand *vertError = mc_trans_in_contin->GetVertErrorBand(errName);
            vertError->Scale(results_cv[0],"",false);
            MnvVertErrorBand *scaleVertError = tmp_scale_factor_trans->GetVertErrorBand(errName);
            scaleVertError->TH1D::Divide( ((TH1D*)tmp_scale_factor_trans), ((TH1D*)tmp_scale_factor_trans), results_cv[0], 1.0);
        }
        if( mc_contin_in_contin->HasVertErrorBand(errName)) {
            MnvVertErrorBand *vertError = mc_contin_in_contin->GetVertErrorBand(errName);
            vertError->Scale(results_cv[1],"",false);
            MnvVertErrorBand *scaleVertError = tmp_scale_factor_contin->GetVertErrorBand(errName);
            scaleVertError->TH1D::Divide( ((TH1D*)tmp_scale_factor_contin), ((TH1D*)tmp_scale_factor_contin), results_cv[1], 1.0);
        }
      
        for (unsigned int i = 0; i < n_universe; ++i) {
            
            std::cout << "\t Uncertainty band name:   " << errName << " universe " << i << std::endl;
            std::cout << "\t\tTotal data:      " << h1d_data_trans_cv->Integral() << std::endl;
            
            TH1D* h1d_signal_trans_universe = GetVertErrorBandUniverseHist(mnvh1d_mc_signal_in_trans, errName, i);
            TH1D* h1d_trans_trans_universe  = GetVertErrorBandUniverseHist(mc_trans_in_trans, errName, i);
            TH1D* h1d_contin_trans_universe = GetVertErrorBandUniverseHist(mc_contin_in_trans, errName, i);
            
            TH1D* h1d_signal_contin_universe = GetVertErrorBandUniverseHist(mnvh1d_mc_signal_in_contin, errName, i);
            TH1D* h1d_trans_contin_universe  = GetVertErrorBandUniverseHist(mc_trans_in_contin, errName, i);
            TH1D* h1d_contin_contin_universe = GetVertErrorBandUniverseHist(mc_contin_in_contin, errName, i);
            
            vector<double> results_universe = CalcPhysScaleFactorMinimizer(h1d_data_trans_cv,
                                                                       h1d_signal_trans_universe,
                                                                       h1d_trans_trans_universe,
                                                                       h1d_contin_trans_universe,
                                                                       h1d_data_contin_cv,
                                                                       h1d_signal_contin_universe,
                                                                       h1d_trans_contin_universe,
                                                                       h1d_contin_contin_universe,
                                                                       mc_scale,
                                                                       errName);
            
            
            assert(results_universe[0] > 0.0);
           
          
            if( mnvh1d_mc_signal_in_contin->HasVertErrorBand(errName)) {
                
                TH1D* universe_trans_trans       = GetVertErrorBandUniverseHist(mc_trans_in_trans, errName, i);
                TH1D* universe_scale_histo_trans = GetVertErrorBandUniverseHist(tmp_scale_factor_trans,   errName, i);
                TH1D* universe_scalefactor_trans = GetVertErrorBandUniverseHist(tmp_scale_factor_trans,   errName, i);
                universe_trans_trans->Scale(results_universe[0]);
                universe_scalefactor_trans->Divide(universe_scale_histo_trans,universe_scale_histo_trans,results_universe[0],1.0);
                
               
                TH1D* universe_contin_contin = GetVertErrorBandUniverseHist(mc_contin_in_contin, errName, i);
                TH1D* universe_scale_histo_contin = GetVertErrorBandUniverseHist(tmp_scale_factor_contin,   errName, i);
                TH1D* universe_scalefactor_contin = GetVertErrorBandUniverseHist(tmp_scale_factor_contin,   errName, i);
                universe_contin_contin->Scale(results_universe[1]);
                universe_scalefactor_contin->Divide(universe_scale_histo_contin,universe_scale_histo_contin,results_universe[1],1.0);
                
                TH1D* universe_trans_contin       = GetVertErrorBandUniverseHist(mc_trans_in_contin, errName, i);
                universe_trans_contin->Scale(results_universe[0]);
               
                TH1D* universe_contin_trans = GetVertErrorBandUniverseHist(mc_contin_in_trans, errName, i);
                universe_contin_trans->Scale(results_universe[1]);
                
                std::cout.precision(8);
                std::cout << "\t\t " << errName << " " << i << " check scale US: " << results_universe[0] << std::endl;
                
            }// end of upstream
            
        }// end of universes
        
    }// end of vertical error band names
    
}

TH1D* GetVertErrorBandUniverseHist(MnvH1D* mnvh1d,
                                              const std::string& errName,
                                              int universe)
{
    //cout << "histo, errName, uni = " << mnvh1d << ", " << errName << ", " << universe << endl;
    return mnvh1d->GetVertErrorBand(errName)->GetHists().at(universe);
}


vector<double> CalcPhysScaleFactorMinimizer( TH1D* h1d_data_trans_cv,
                                                       TH1D* h1d_signal_trans_cv,
                                                       TH1D* h1d_trans_trans_cv,
                                                       TH1D* h1d_contin_trans_cv,
                                                       TH1D* h1d_data_contin_cv,
                                                       TH1D* h1d_signal_contin_cv,
                                                       TH1D* h1d_trans_contin_cv,
                                                       TH1D* h1d_contin_contin_cv,
                                                       double mc_scale,
                                                       const std::string& errName)
{
    
    std::cout<<"I'm inside the physics scale minimizer" <<std::endl;
    
    //    cout << "nominal event before fit ( data, signal, other, plastic US, plastic DS ) : " << std::fixed << std::setprecision(2) << h1d_data->Integral() << " & "
    //    << std::fixed << std::setprecision(2) <<  h1d_mc_signal->Integral() << " & "
    //    << std::fixed << std::setprecision(2) <<  h1d_mc_other->Integral() << " & "
    //    << std::fixed << std::setprecision(2) <<  h1d_mc_backgroundUS->Integral() << " & "
    //    << std::fixed << std::setprecision(2) <<  h1d_mc_backgroundDS->Integral() << "\\\\" << endl;
    
    m_histo_trans_data   = h1d_data_trans_cv;
    m_histo_sig_trans    = h1d_signal_trans_cv;
    m_histo_trans_trans  = h1d_trans_trans_cv;
    m_histo_contin_trans = h1d_contin_trans_cv;
    m_histo_contin_data   = h1d_data_contin_cv;
    m_histo_sig_contin    = h1d_signal_contin_cv;
    m_histo_trans_contin  = h1d_trans_contin_cv;
    m_histo_contin_contin = h1d_contin_contin_cv;
    
    std::cout<<"I've got my histograms now" <<std::endl;
    
    
    // initialize stuff
    double scale_factor_trans  = 1.0;
    double scale_factor_contin = 1.0;
    double minChi2 = -1.0;
    
    vector<double> results;
    results.push_back(scale_factor_trans);
    results.push_back(scale_factor_contin);
    results.push_back(minChi2);
    
    // Pick some numbers to be the tolerance and stuff
    ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2");
    fitter->SetMaxFunctionCalls(1000000);
    fitter->SetMaxIterations(100000);
    fitter->SetTolerance(0.01);
    fitter->SetPrintLevel(1);
    
    std::cout<<"Made a Minimizer" <<std::endl;
    
    
    //BkgFitter::getChi2 has 1 parameter, minimizes chi2 for joint fit
    // initialize scale to 1.0 with step size of 0.001
    //    fitter->SetVariable( 0, "scale", 1.0, 0.001 );
   
//    ROOT::Math::Functor lf( this, &getPhysChi2, 2 );
      ROOT::Math::Functor lf( & getPhysChi2, 2 );
    ROOT::Math::Functor functor( lf, 2 );
    fitter->SetFunction( functor );
    fitter->SetLimitedVariable( 0, "trans", 3.0, 0.01, 0.1, 5.0 );
    fitter->SetLimitedVariable( 1, "contin", 3.0, 0.01, 0.1, 5.0 );
    //    fitter->setFitRange(2.0,25.0);

    std::cout<<"Defined my functor " <<std::endl;
    
    // Go!
    fitter->Minimize();
    
    std::cout<<"I minimized!" <<std::endl;
    if( fitter->Status() != 0 ) {
        std::cout << "Womp womp womp... The background fit failed!" << std::endl;
        return results;
    }
    
    const double *bestfit = fitter->X();
    minChi2 = getPhysChi2( bestfit );
    results[0] = bestfit[0]; //scale_factor_trans;
    results[1] = bestfit[1]; //scale_factor_contin;
    results[2] = minChi2;
    
    std::cout<<"The results are "<<results[0]<<"  "<<results[1]<<"  "<<results[2]<<std::endl;
    return results;
    
}


double getPhysChi2( const double * par )
{
   
    double scale_trans  = par[0];
    double scale_contin = par[1];
    
    // Sideband is bins from plane number 11 to 65
    // Assume no error on the MC, and use the stat errors on the data
    int binLow  = m_histo_trans_data->FindBin( 2.0  ); // same bins on both sets of histograms
    int binHigh = m_histo_trans_data->FindBin( 25.0 );
    
    double chi2 = 0.0;
    // clone and scale the floating histograms
    TH1D * temp_trans_float = new TH1D( *(TH1D*)m_histo_trans_trans->Clone("temp_trans_float") );
    TH1D * temp_contin_trans_float = new TH1D( *(TH1D*)m_histo_contin_trans->Clone("temp_contin_trans_float") );
    
    temp_trans_float->Scale( scale_trans );
    temp_contin_trans_float->Scale( scale_contin );
    temp_trans_float->Add(temp_contin_trans_float);
    temp_trans_float->Add(m_histo_sig_trans);
    
    TH1D * temp_contin_float = new TH1D( *(TH1D*)m_histo_contin_contin->Clone("temp_contin_float") );
    TH1D * temp_trans_contin_float = new TH1D( *(TH1D*)m_histo_trans_contin->Clone("temp_trans_contin_float") );
    temp_contin_float->Scale( scale_contin );
    temp_trans_contin_float->Scale( scale_trans );
    temp_contin_float->Add(temp_trans_contin_float);
    temp_contin_float->Add(m_histo_sig_contin);
    
    double tot_trans_data   = 0.0;
    double tot_trans_float  = 0.0;
    double tot_contin_data  = 0.0;
    double tot_contin_float = 0.0;
    
    for( int b = binLow; b <= binHigh; ++b ) {
        double trans_data   = m_histo_trans_data->GetBinContent(b);
        double trans_float  = temp_trans_float->GetBinContent(b);
        
        if( trans_data != 0.0 ){
            tot_trans_data  += trans_data;
            tot_trans_float += trans_float;
            
            chi2 += pow( (trans_float - trans_data)/sqrt(trans_data),2);
        }
    }
    
    for( int b = binLow; b <= binHigh; ++b ) {
        double contin_data   = m_histo_contin_data->GetBinContent(b);
        double contin_float  = temp_contin_float->GetBinContent(b);
        
        if( contin_data != 0.0 ){
            tot_contin_data  += contin_data;
            tot_contin_float += contin_float;
            
            chi2 += pow( (contin_float - contin_data)/sqrt(contin_data),2);
        }
    }
    
    double unbinned_chi2 = pow( (tot_contin_float-tot_contin_data)/sqrt(tot_contin_data), 2 );
    unbinned_chi2 += pow( (tot_trans_float-tot_trans_data)/sqrt(tot_trans_data), 2 );
    
    
//    cout<<"The chi2 integrated is "<<unbinned_chi2<<"  and the unintegrated is "<<chi2<<endl;
    return unbinned_chi2;
    
     //return chi2;
    
}

MnvH1D*  TunePlasticBackground( MnvH1D *histo, string region, string outdir, int targetID, int targetZ, string playlist, string var ){
    cout<<"OutDir = "<<outdir<<"  Z = "<<targetZ<<"   playlist = "<<playlist<<"   var = "<<var<<endl;  
    TFile *scaleFactor = new TFile( Form("%s/CHScaleFactors_t%d_z%d_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str() ), "read" );
    assert(scaleFactor);
    
    MnvH1D* h_scaleFactors = (MnvH1D*)scaleFactor->Get( Form("scaleFactor_%s_%s", region.c_str(), var.c_str() ) ); //region should be "US" or "DS"
    cout << "name of scale factors = " << scaleFactor->GetName() << endl;
    int nonzerobin = h_scaleFactors->FindFirstBinAbove(0.);
    MnvH1D* tunedHisto = (MnvH1D*)histo->Clone("tunedHisto");
    MnvH1D* scaleHisto = (MnvH1D*)histo->Clone("scaleHisto");
    scaleHisto->Reset();
    cout<<" I cloned the histogram "<<endl;
    double cv_val = h_scaleFactors->GetBinContent(nonzerobin);
    cout<<"I've got the scale factor stuff for "<<region.c_str()<<"  the CV value is "<<cv_val<<endl;
    
    for(int i=0;i<scaleHisto->GetNbinsX()+2;i++) scaleHisto->SetBinContent(i,cv_val);
  
    tunedHisto->Multiply(tunedHisto,scaleHisto);
    cout<<"done tuning plastic SB"<<endl; 
    return tunedHisto;
    
}


//bool UnfoldHisto( MnvH1D* hUnfolded, MnvH1D *hFolded, const std::string& var, int targetID, int targetZ, string playist, HelicityType::t_HelicityType helicity, bool accCorr /*=false*/, bool addSysErrors /*= true*/, Unfold::t_UnfoldMethod method )
bool UnfoldHisto( MnvH1D* hUnfolded, MnvH1D *hFolded, const std::string& var, string outdir, int targetID, int targetZ, string playist, HelicityType::t_HelicityType helicity, bool accCorr /*=false*/)
{
    
    TH1::AddDirectory(kFALSE);
    hUnfolded->Reset();
    hUnfolded->ClearAllErrorBands();
  
    NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF();
    
    bool unfoldOK = false;
    //bool errsOK   = false;
    //const int errIter = 4;
    //f2->cd();
    
    bool forceDisk = true;
    //TString histFileName = utils->GetHistFileName( "Migration", FileType::kAny, targetID, targetZ, helicity, forceDisk );
    TString migration = Form("%s/Hists_Migration_t%d_z%02d_Nu.root", outdir.c_str(), targetID, targetZ);
    //TString histFileName = utils->GetHistFileName( "Migration", FileType::kAny, targetID, targetZ, helicity);
    //TFile fMigration( migration, "READ" );
    TFile *fMigration = new TFile( migration,"read" );

    cout << "Start Unfolding" << endl;
    MnvH2D *hResp  = dynamic_cast<MnvH2D*> (fMigration->Get( Form("selected_Migration_%s_%s",  var.c_str(), var.c_str())));
    //MnvH2D *hResp  = dynamic_cast<MnvH2D*> (fMigration->Get( Form("selected_Migration_Emu_Ehad")));
    //TString hNameResp = utils->GetHistName( ( accCorr ? "migration_accCorr" : "migration" ), FileType::kAny, var, targetID, targetZ );
    //MnvH2D *hResp = dynamic_cast<MnvH2D*>( fMigration.Get( hNameResp ) );
    //MnvH2D *hResp = dynamic_cast<MnvH2D*>( fMigration.Get( hNameResp ) );
    if( 0 == hResp )
    {
        Error( "NukeUtils::UnfoldHisto", Form("Could not find response migration histogram for var %s, targetID %d, targetZ %d", var.c_str(), targetID, targetZ ) );
        return false;
    }
    cout<<"Ok, I have my migration matrix now"<<endl;
//    if( hFolded->HasVertErrorBand( ErrorName::PlasticBG ) && !hResp->HasVertErrorBand( ErrorName::PlasticBG ) ){
        cout<<"I added a plsticBG error band to my response matrix "<<endl;
    ////////    hResp->PopVertErrorBand( ErrorName::PlasticBG );
    ///////    hResp->AddVertErrorBandAndFillWithCV( ErrorName::PlasticBG, 2 );
//    }
    
    //int nIter = UnfoldIterations( var );
    int nIter = 8;
    /*if("Enu" = var)
       {
        nIter = 5;
       }
    else nIter = 8;
*/
    cout << "Number of iterations " << nIter << endl;
//    hResp->PopVertErrorBand(ErrorName::Hadronization);
//    if(addSysErrors){
//        if( hResp->HasVertErrorBand("CCQEPauliSupViaFK") )
//            hResp->PopVertErrorBand("CCQEPauliSupViaFK");
//        if( hFolded->HasVertErrorBand("CCQEPauliSupViaFK") )
//            hFolded->PopVertErrorBand("CCQEPauliSupViaFK");
//    }
//    cout<<"Ok, I have my CCQE error bands"<<endl;
    
    if(nIter == 0){
        TString histName = hUnfolded->GetName();
        hUnfolded = (MnvH1D*) hFolded->Clone(histName);
        return true;
    }
    cout<<"Ok, I have my folded distribution now"<<endl;
    std::vector<std::string> mnv_names = hFolded->GetVertErrorBandNames();
    
    //            for(int i=0;i<mnv_names.size();i++)
    //                std::cout<<"folded name "<<mnv_names[i]<<std::endl;
    cout<<"Ok, I have my migration matrix now"<<endl;
    
    hResp->AddMissingErrorBandsAndFillWithCV(*hFolded);
   /* 
    if(hResp->HasLatErrorBand(ErrorName::MuonAngle)) hResp->PopLatErrorBand(ErrorName::MuonAngle);
    if(hResp->HasLatErrorBand(ErrorName::MuonE)) hResp->PopLatErrorBand(ErrorName::MuonE);
*/
    std::vector<std::string> mnv_names_resp = hResp->GetVertErrorBandNames();
            for(int i=0;i<mnv_names_resp.size();i++)
                std::cout<<"migration name "<<mnv_names_resp[i]<<std::endl;
    
    
    
    cout<<"number of error bands "<<mnv_names.size()<<"  "<<mnv_names_resp.size()<<endl;
    
    //    MnvH2D *hResp_temp = dynamic_cast<MnvH2D*>( hResp );
    
    //    if(mnv_names.size() == 0)
    //        hResp_temp = new MnvH2D( hResp->GetCVHistoWithStatError() );
    //    else
    //        MnvH2D *hResp_temp = dynamic_cast<MnvH2D*>( hResp );
    
    
    TMatrixD covmx;
    
    //unfoldOK = MinervaUnfold::MnvUnfold::Get().UnfoldHisto( hUnfolded, covmx, hResp, hFolded, RooUnfold::kBayes, nIter, addSysErrors);
    unfoldOK = MinervaUnfold::MnvUnfold::Get().UnfoldHisto( hUnfolded, covmx, hResp, hFolded, RooUnfold::kBayes, nIter );
//       unfoldOK = MinervaUnfold::MnvUnfold::Get().UnfoldHisto( hUnfolded, hResp, hFolded, RooUnfold::kBayes, nIter, true, true );
    
    cout<<"I actually ran the unfolding now... "<<endl;
    
    //change this for data to get the stat covariance matrix
    
    //for the data, assign a statistical error based on errIter iterations
 /*   if(!addSysErrors){
     errsOK = MinervaUnfold::MnvUnfold::Get().UnfoldHisto( hUnfoldErr, hResp, hFolded, RooUnfold::kBayes, errIter, addSysErrors, addSysErrors );
     if(errsOK){
     for(int bin = 1; bin < hUnfolded->GetNbinsX()+1; ++bin){
     //cout << hFolded->GetBinError(bin) << " " << hUnfoldErr->GetBinError(bin) << endl;
     
     hUnfolded->SetBinError(bin,hUnfoldErr->GetBinError(bin) );
     
     }
     
     }
     
     else
     Error( "NukeUtils::UnfoldHisto", Form("Could not unfold stat errors for var %s, targetID %d, targetZ %d", var.c_str(), targetID, targetZ ) );
     
     
     }//end of !addSysErrors */
    cout<<"I finished unfolding ok "<<endl;
    return unfoldOK;
}
/*
int UnfoldIterations( const std::string& var )
{

    if(fUnfoldIter != - 1)
        return fUnfoldIter;

    if(
       "Emu" == var ||
       "Ehad" == var ||
       "Enu" == var ||
       "ThetaMu" == var
       )
        return 2;

    else if( "x" == var )
        return 2;

    return 0;
}
*/






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
  
  TString labelStat = Form("Stat. + Sys. Error, #chi^{2}/ndf = %3.2f/%d = %3.2f", chi2Sys, ndfStat, chi2Sys/(Double_t)ndfStat);
  
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
    cout << "data_cv =" << data_cv->GetBinContent(2) << endl;
    cout << "mc_cv =" << totalMC->GetBinContent(2) << endl;
    //mnvPlotter.DrawDataMCRatio( data_cv, totalMC, dataMCScale, true, true, plotMin, plotMax, Yaxis );
    mnvPlotter.DrawDataMCRatio( data_cv, totalMC, dataMCScale, true, true, 0., 2., x_label.c_str(), Yaxis );
    mnvPlotter.AddPlotLabel( labelStat, .46, .875, .05 );
    cout<<"Is this correct number = "<<dataMCScale<<endl;
 
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
  std::string x_label = var.c_str(); 
    cout << "plotting " << cName << endl;
    //totalMC->PopVertErrorBand(ErrorName::PlasticBG);
    //cout<<"Canname:" <<canname<<"\t totalMC Plastic_SB? "<< totalMC->HasUncorrError(ErrorName::CHSB)<<endl;
    //mnvPlotter.DrawDataMCErrorSummary( totalMC, data_cv, /*cname,*/ legLoc, false, true, 0.00001, false ); //use my colors and styles
    PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
    //PlotUtils::MnvPlotter plotter(kCCQENuStyle);
    //plotter.DrawErrorSummary(totalMC, legLoc);
    //mnvPlotter.DrawErrorSummary( totalMC, legLoc, false, true, 0.0001, true, "", "", "", false); 
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


//MnvH1D* ReBinPhysSBTuning( MnvH1D* h_scaleFactors, string SB, int targetID, int targetZ, string playlist, string var, bool cv_only, NukeCC_Binning *binsDef,  MnvH1D* h_sb){
MnvH1D* ReBinPhysSBTuning( MnvH1D* h_scaleFactors, string SB, int targetID, int targetZ, string playlist, string var, bool cv_only,  MnvH1D* h_sb){
    cout<<"I'm saving the other binnings for the sidebands now "<<endl;
    int temp_tgtZ = 0;
    if(targetID > 10)
        temp_tgtZ = 99;
    if(targetID< 10)
        temp_tgtZ = targetZ;
    NukeCC_Binning  *binsDef = new NukeCC_Binning();     
    int nonzerobin = h_scaleFactors->FindFirstBinAbove(0.);

    vector<double> binsLowEdge = binsDef->GetEnergyBins( var, true ); //Do I really want this? The sideband histograms are using SidebandBins (only relevant for Q2 and W) probably?
    vector<double> SBbins = binsDef->GetSidebandBins( var ); 
    TString hName = Form("scale_factor_%s_%s", SB.c_str(),var.c_str());
    MnvH1D* tunedHisto = new MnvH1D( hName, hName, SBbins.size()-1, &SBbins[0] ); //binsLowEdge.size()-1, &binsLowEdge[0] ); 
    MnvH1D* scaleHisto = (MnvH1D*)tunedHisto->Clone(Form("scale_factor_%s_%s", SB.c_str(),var.c_str()));
    scaleHisto->Reset();
    scaleHisto->ClearAllErrorBands();
    //cout<<" I cloned the histogram "<<endl;
    double cv_val = h_scaleFactors->GetBinContent(nonzerobin);
    cout<<"I've got the scale factor stuff for SB "<<SB.c_str()<<"  the CV value is "<<cv_val<<endl;

    
    for(int i=0;i<scaleHisto->GetNbinsX()+2;i++) scaleHisto->SetBinContent(i,cv_val);
    
    std::vector<std::string> vertNames = h_scaleFactors->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k ){
        MnvVertErrorBand *errBand = h_scaleFactors->GetVertErrorBand( vertNames[k] );
        const int universes = errBand->GetNHists();
        std::vector<TH1D*> vert_hists;
        for(int u=0;u<universes;++u){
            TH1D* tmp_scale = new TH1D(*errBand->GetHist( u ));
            TH1D* tmp_template = new TH1D(scaleHisto->GetCVHistoWithStatError());
            tmp_template->SetName(Form("Scale_universe_%d",u));
            for(int i=0;i<scaleHisto->GetNbinsX()+2;i++) tmp_template->SetBinContent(i,tmp_scale->GetBinContent(nonzerobin));
            vert_hists.push_back(tmp_template);
        }
        scaleHisto->AddVertErrorBand( vertNames[k],vert_hists);
    }
    // Uncorrelated errors
    if(h_sb != NULL){
      int offset=0;
      if(binsLowEdge[0] != SBbins[0]){
	cout<<"The EnergyBins and SidebandBins differ."<<endl;
	/*for(int j=0; j<=h_sb->GetNbinsX()+2; ++j){
	  if (binsLowEdge[0] == SBbins[j]){
	    offset=j;
	    cout<<"Var: "<<var.c_str()<<"\t offset: "<<j<< "Value: "<<SBbins[j]<<endl;
	    break;
	  }
	}//end loop over bins
	if(offset==0){
            Error( "NukeUtils::ReBinPhysSBTuning", "Inconsistent binnings between EnergyBins and SidebandBins. Check that SidebandBins contains EnergyBins" );
            throw 1;
	}
	*/
      }//end if case for different binnings
    cout<<"What about here"<<endl;
      std::vector<std::string> uncorrNames = h_sb->GetUncorrErrorNames();
      scaleHisto->AddMissingErrorBandsAndFillWithCV(*h_sb);
      
      for (std::vector<std::string>::iterator name = uncorrNames.begin(); name != uncorrNames.end(); ++name){
	std::string errName = *name;
	TH1D* sferr= scaleHisto->GetUncorrError(errName);
	TH1D* tmp_err= h_sb->GetUncorrError(errName); 
	

	for( int i=0; i<=sferr->GetNbinsX()+2;++i){//keep the fractional error from the untuned sideband
	  sferr->SetBinContent(i,cv_val); //make sure CV is consistent
	  double mccontent=tmp_err->GetBinContent(i);
	  double mcerr=tmp_err->GetBinError(i);
	  if(mccontent==0.0) sferr->SetBinError(i, 0.0);
	  else sferr->SetBinError(i, cv_val * mcerr/mccontent);
	  
	} 

      }//end loop over uncorr errors
    }// end if h_sb not NULL



    cout<<"Scale histogram CV "<<scaleHisto->GetBinContent(2)<<endl;
    return scaleHisto;
//    scaleHisto->Write();
//    scaleFactor->Close();
}
