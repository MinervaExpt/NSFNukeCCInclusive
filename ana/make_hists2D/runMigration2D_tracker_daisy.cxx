//==============================================================================
// Loop entries, make cuts, fill histograms.
// * Uses the New Systematics Framework and "Universe" objects.
// * loop universes, make cuts and fill histograms with the correct lateral
// shifts and weights for each universe.
// * TChain --> PlotUtils::ChainWrapper.
// * MnvHXD --> PlotUtils::HistWrapper.
// * Genie, flux, non-resonant pion, and some detector systematics calculated.
//==============================================================================

#include "../../NUKECCSRC/include/CommonIncludes.h"
#include "../../NUKECCSRC/include/CVUniverse.h"
#include "../include/VariableDaisy.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../NUKECCSRC/include/Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include <iostream>
#include <stdlib.h>
#include "../../NUKECCSRC/include/UtilsNSF.h"
#include "../../NUKECCSRC/include/Cuts.h"
#include "TParameter.h"
#include "PlotUtils/MnvH2D.h"


#include "../include/systematics/Systematics.h"

// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we 
// shield it.
#ifndef __CINT__
#include "../include/plotting_functions.h"
#endif
#include "PlotUtils/MacroUtil.h" 
//using namespace globalV;
using namespace NUKECC_ANA;
//======================================================================
typedef VarLoop::Variable Var;
typedef Var2DLoop::Variable2D Var2D;
typedef PlotUtils::MnvH2D MnvH2D;



void FillVariable(TFile* fout, PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d, bool isMC, int targetID=1, int targetZ=26, const string playlist="minervame1A", bool doDIS=true);


//============================================================================================================================
// Main
//============================================================================================================================

int main(int argc, char *argv[]){
	 //ROOT::Cintex::Cintex::Enable();
	 TH1::AddDirectory(false);

	 if(argc==1){
	     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
	     std::cout<<"MACROS HELP:\n\n"<<
	       "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n"<<
	       "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
	       "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
	       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" << std::endl;
	     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
	     return 0;
	 }

	 TString dir(argv[1]);
	 //int targetID = atoi(argv[2]);
	 //int targetZ  = atoi(argv[3]);
  int targetID = 99; int targetZ = 99; 

		bool doDIS = false; 
	 // TString dir(argv[1]);
	 // int targetID = 1;
	 // int targetZ = 26;
	 // const string playlist= argv[4];
   	const string playlist= argv[2];

    const std::string plist_string(playlist);    

    //const std::string mc_file_list("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/include/playlists/try.txt"); 
    //const std::string data_file_list("/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/include/playlists/MasterAnaDev_data_AnaTuple_run00022246_Playlist.root"); 

    const std::string mc_file_list(Form("../include/playlists/MasterAnaDev_MC_%s.txt", plist_string.c_str()));
    const std::string data_file_list(Form("../include/playlists/MasterAnaDev_Data_%s.txt",plist_string.c_str()));
    const std::string reco_tree_name("MasterAnaDev");

    //const std::string mc_file_list("../include/playlists/NukeCC_MC_minervame6A_MuonKludged.txt");
    //const std::string data_file_list("../include/playlists/NukeCC_Data_minervame6A_MuonKludged.txt");
    //const std::string reco_tree_name("NukeCC");
  
    const bool wants_truth = false;
    //const bool is_grid = false;
    // is grid removed after update of MAT 07/12/2021

	 PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth);

	 util.PrintMacroConfiguration("main");

	 //=========================================
	 // Systematics
	 //=========================================
	 //std::map<std::string, std::vector<CVUniverse*> > error_bands =
	 //  GetErrorBands(util.m_mc);

	 PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
	 PlotUtils::MinervaUniverse::SetNuEConstraint(true);
	 PlotUtils::MinervaUniverse::SetAnalysisNuPDG(-14);
	 PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
   PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
   PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(true);
   PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);
   // Defined for MnvHadronReweighter (GEANT Hadron sytematics)
  //Tracker or nuke (what clusters are accepted for reconstruction)
  PlotUtils::MinervaUniverse::SetReadoutVolume("Nuke");
  //Neutron CV reweight is on by default (recommended you keep this on)
  PlotUtils::MinervaUniverse::SetMHRWeightNeutronCVReweight(true);
  //Elastics are on by default (recommended you keep this on)
  PlotUtils::MinervaUniverse::SetMHRWeightElastics(true);

   
	    
	 NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
	 NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
	 NukeCC_Binning  *binsDef = new NukeCC_Binning();
	  
	 PlotUtils::ChainWrapper* chainData = util.m_data;
	 PlotUtils::ChainWrapper* chainMC = util.m_mc;
	 HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);
	 double DataPot=  util.m_data_pot; 
	 double MCPot=  util.m_mc_pot; 
	 //double total_pot_data,total_pot_mc;
	 //utils->getPOT(total_pot_data,total_pot_mc);  
	 double  MCscale=DataPot/MCPot;
	 //double  MCscale=1.0;
	 
	 std::cout<<" MCScale= "<<MCscale<<std::endl; 
	 std::vector<Var*> variablesMC;
	 std::vector<Var2D*> variables2DMC;

         //TString histFileName = utils->GetHistFileName( "Migration", FileType::kAny, targetID, targetZ );
	TString histFileName;
  if(RunCodeWithSystematics){
    histFileName += Form("/Migration2D_Daisy_%s_t%d_z%02d_sys.root", plist_string.c_str(), targetID, targetZ);
  }

  else{
    histFileName += Form("/Migration2D_Daisy_%s_t%d_z%02d_nosys.root", plist_string.c_str(), targetID, targetZ);
  }  
   	   
	TFile* fout = new TFile(dir.Append(histFileName),"RECREATE");
  	   

   FillVariable(fout, chainMC, helicity, utils, cutter, binsDef, variablesMC, variables2DMC, true, targetID, targetZ, plist_string, doDIS);

	     
         	 
	 //for (auto v : variables2DMC) {
      //for(int petal=0; petal<4; petal++){
       // v-> mresp.SyncCVHistos();
      //}
   //}


	 //For 2D variable

	 //for (auto v : variables2DMC) {
    //fout.cd();
    //for(int petal=0; petal<4; petal++){
    //  mresp[petal].hist.Write()
    //}
	 // v->WriteAllHistogramsToFileMig(fout,true);
	 //}

  fout->cd();
  auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  dataPOTOut->Write();
  mcPOTOut->Write(); 
	 
  std::cout << "DONE" << std::endl;


}//End Main


//============================================================================================================================
// FillVariable 
//============================================================================================================================
    
void FillVariable(TFile* fout, PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d, bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){

   // std::map< std::string, std::vector<CVUniverse*> > error_bands = utils->GetErrorBands(chain)
   
   std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
   
   std::map<std::string, std::vector<CVUniverse*> >::iterator itr;
   
   std::map<const std::string, int> error_name;
   for(itr = error_bands.begin(); itr != error_bands.end(); ++itr) error_name.insert(pair<std::string, const int>((itr->second)[0]->ShortName(), (itr->second).size())); 

   std::map<std::string, const int>::iterator itr_m;

   std::vector<double> pTbin, pZbin;

    
    pTbin = binsDef->GetEnergyBins("muonPt"); 
    pZbin = binsDef->GetEnergyBins("muonPz"); 

   Var* pTmu = new Var("pTmu", "pTmu", pTbin, &CVUniverse::GetMuonPt, &CVUniverse::GetlepPtTrue);
   Var* pZmu = new Var("pZmu", "pZmu", pZbin, &CVUniverse::GetMuonPz, &CVUniverse::GetlepPzTrue);

   //For 2D variable

   Var2D* pZmu_pTmu = new Var2D(*pZmu, *pTmu);
  
   variables2d = {pZmu_pTmu};

   for (auto v : variables2d) v->InitializeAllHistograms(error_bands);

   //Migration starts here!
   MinervaUnfold::MnvResponse* Response[12];

  for (auto v : variables2d){
    for(int petal=0; petal<12; petal++){
      MinervaUnfold::MnvResponse* resp = Response[petal];
      Response[petal] = v->SetupResponse(resp, petal, error_name);
    }
  }


    
  int reco0=0;
  int reco1=0;
  int reco2=0;
  int reco3=0;
  int reco4=0; 
  int reco5=0; 
  int reco6=0;  
  int reco7 = 0;
  int reco8=0;
   
   CVUniverse *dataverse = new CVUniverse(chain,0);

   
    
   //=========================================
   // Entry Loop
   //=========================================

   std::cout<<"# of entries = "<<chain->GetEntries()<<std::endl;

   for(int i=0; i<chain->GetEntries(); ++i){
     if(i%500000==0) std::cout << (i/1000) << "k " << std::endl;
     //=========================================
     // For every systematic, loop over the universes, and fill the
     // appropriate histogram in the MnvH1D
     //=========================================

       for (auto band : error_bands){
	       int unv_count = 0;
               std::vector<CVUniverse*> error_band_universes = band.second;
	       for (auto universe : error_band_universes){
		 // Tell the Event which entry in the TChain it's looking at
		 //=========================================
		 // CUTS in each universe
		 //========================================
	 
		 universe->SetEntry(i);
     reco0++;
     // reco cuts

     if( universe->GetInt("muon_corrected_p") == -999 ) continue; // additional cut to get rid of an issue
     if(!cutter->PassReco(universe,helicity)) continue;
     reco1++;

	   if(!cutter->TrackerOnly(universe)) continue;
     reco2++;

     if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
     reco4++;

      if (!cutter->PassTrueCC(universe, helicity)) continue; //true CC, true antinu
      // NO xy separation,  NO APOTHEM CUT
      reco5++;

	   if(!cutter->TrackerOnlyTrue(universe)) continue; // true tracker  
     reco6++; 

     int petal = TargetUtils::Get().GetDaisyPetal(universe->GetVertexXMy()*10, universe->GetVertexYMy()*10 ); // multiply from cm to mm


		 for (auto v : variables2d){
	    //if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(universe)) continue;
	    if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(universe)) continue;
      if( v->GetName()=="pZmu_pTmu") reco7++;

      //if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassTrueMuEnergyCut(universe)) continue;
	    //if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassTrueThetaCut(universe)) continue;
      // NO TRUE angle cut, efficiency corrected
	    if( v->GetName()=="pZmu_pTmu") reco8++;   
      
      //v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
      
      //Migration stuff
      //std::cout<<"Tried to fill"<<std::endl;
      v->FillResponse(Response[petal], v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),universe->ShortName(),universe->GetWeight(),unv_count); 
      //std::cout<<"Filled one"<<std::endl;
		 }

		 unv_count++;
    } // End band's universe loop
     }
  }//End entries loop
 std::cout<<"Filling DONE"<<std::endl;
 //for (auto v : variables2d) v->getResponseObjects(Response, migrationH, h_reco, h_truth, error_bands);
 //PlotUtils::Hist2DWrapper<CVUniverse> mresp;
 for (auto v : variables2d){
  for(int petal=0; petal<12; petal++){
    fout->cd();
    Response[petal]->GetMigrationMatrix()->Write();
	  Response[petal]->GetTruth2D()->Write();
	  Response[petal]->GetReco2D()->Write();
    
    /*MnvH2D *migrationH = NULL;
    MnvH2D *h_reco = NULL;
    MnvH2D *h_truth = NULL;
    bool check = Response[petal]->GetMigrationObjects( migrationH, h_reco, h_truth );
    const bool clear_bands = false;  
    mresp = PlotUtils::Hist2DWrapper<CVUniverse> (migrationH, error_bands, clear_bands);
    mresp.SyncCVHistos();
    fout->cd();
    mresp.hist->Write();
    h_reco->Write();
    h_truth->Write();


    std::cout<<"Petal " << petal <<std::endl;
    
    delete Response[petal];
    delete migrationH;
    delete h_reco;
    delete h_truth;
    */
    //mresp[petal] = v->getResponseObjects(Response[petal], error_bands);
  }
 }
 
 for(auto band : error_bands){
    std::vector<CVUniverse*> band_universes = band.second;
    for(unsigned int i_universe = 0; i_universe < band_universes.size(); ++i_universe)  
      delete band_universes[i_universe];
 } 

    

   std::cout<<"**********************************"<<std::endl;
      std::cout<<" Summary "<<std::endl;
   std::cout<<"**********************************"<<std::endl;
     std::cout<<"Migration Matrix "<<std::endl;
    std::cout << "No cuts = " << reco0 << std::endl;
    std::cout << "Reco Cut = " << reco1 << std::endl;
    std::cout << "Tracker Cut = " << reco2 << std::endl;
    std::cout << "Plane prob. cut = " << reco4 << std::endl;
    std::cout << "Truth cut (fiducial, CC, antinu) = " << reco5 << std::endl;
    std::cout << "True Tracker cut  = "<< reco6 << std::endl;
    std::cout << "Muon kinematics cut  = " << reco7 << std::endl;
    std::cout << "True muon kinematics cut (energy)  = " << reco8 << std::endl;

}
