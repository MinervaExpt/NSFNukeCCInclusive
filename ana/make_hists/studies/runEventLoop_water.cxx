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
#include "../include/VariableRun.h"  
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

//void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType
//helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef,
//std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID=1,
//int targetZ=26, const string playlist="minervame1A", bool doDIS=true);

//=============================================================================
//=============================================================================
// MAIN FUNCTION
//=============================================================================
//=============================================================================

void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType
helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef,
std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID=1,
int targetZ=26, const string playlist="minervame1A", bool doDIS=true);

int main(int argc, char *argv[]){
  //ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if(argc==1){
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "MACROS HELP:\n\n" <<
    "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n" <<
    "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
    "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n" <<
    "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    return 0;
  }

  TString dir(argv[1]);
  //int targetID = atoi(argv[2]);
  //int targetZ = atoi(argv[3]);
  int targetID = 0; int targetZ = 0;  

  bool doDIS=false;
  const string playlist= argv[2];


  // MasterAnaDev tuples?
  //const std::string mc_file_list("../include/playlists/OfficialMAD_minervaME6B_MCNukeOnly_merged.txt");
  //const std::string data_file_list("../include/playlists/OfficialMAD_minervaME6B_Data_merged.txt");
  //const std::string reco_tree_name("MasterAnaDev");

  // NukeCC Tuples ?
  const std::string plist_string(playlist);    

  const std::string mc_file_list(Form("../include/playlists/MasterAnaDev_MC_%s.txt", plist_string.c_str()));
  const std::string data_file_list(Form("../include/playlists/MasterAnaDev_Data_%s.txt",plist_string.c_str()));
  const std::string reco_tree_name("MasterAnaDev");
  
  const bool wants_truth = false;
  //const bool is_grid = false;
  // is grid removed after update of MAT 07/12/2021

  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth);

  util.PrintMacroConfiguration("main");

  // SYSTEMATICS

  //std::map<std::string, std::vector<CVUniverse*> > error_bands =
  //GetErrorBands(util.m_mc);

  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(-14);
  PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100); 
  PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);
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
  double MCscale=DataPot/MCPot;
 
  std::cout << "MC Scale = " << MCscale << std::endl; 
  std::cout << "Data POT: " << DataPot << std::endl;
  std::cout << "MC POT: " << MCPot << std::endl;

  std::vector<Var*> variablesMC,variablesData; 
  std::vector<Var2D*> variables2DMC,variables2DData; 

  TString histFileName;
  if(RunCodeWithSystematics){
    histFileName += Form("/EventSelection_%s_t%d_z%02d_sys.root", plist_string.c_str(), targetID, targetZ);
  }

  else{
    histFileName += Form("/EventSelection_%s_t%d_z%02d_nosys.root", plist_string.c_str(), targetID, targetZ);
  } 

  //TString histFileName = utils->GetHistFileName( "EventSelection_ML_ME6A", FileType::kAny, targetID, targetZ, helicity ); 

  //Works good for the grid submission
  //TString histFileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ );

  TFile fout(dir.Append(histFileName),"RECREATE");	
   
  // MC 
  std::cout << "Processing MC and filling histograms" << std::endl;

  FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true,targetID, targetZ, plist_string,doDIS);     
  for (auto v : variablesMC) v->m_selected_mc_reco.SyncCVHistos();
  for (auto v : variables2DMC) v->m_selected_mc_reco.SyncCVHistos();
   
  // DATA
  std::cout << "Processing Data and filling histograms" << std::endl;

  FillVariable(chainData, helicity, utils, cutter,binsDef,variablesData,variables2DData,false,targetID, targetZ, plist_string,doDIS);
  for (auto v : variablesData) v->m_selected_data_reco.SyncCVHistos();
  for (auto v : variables2DData) v->m_selected_data_reco.SyncCVHistos();

  // WRITE HISTOGRAMS TO FILE

  // 1D variables
  for (auto v : variablesMC) {
    v->WriteAllHistogramsToFile(fout, true);
  }

  for (auto v : variablesData) {
    v->WriteAllHistogramsToFile(fout, false);
  }

  // 2D Variables
  for (auto v : variables2DMC) {
    v->WriteAllHistogramsToFile(fout,true);
  }

  for (auto v : variables2DData) {
    v->WriteAllHistogramsToFile(fout,false);
  }
	
  //Writing POT to the HistFile
  fout.cd();
  auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  dataPOTOut->Write();
  mcPOTOut->Write(); 

  std::cout << "DONE" << std::endl;

}//End Main


//=============================================================================
//=============================================================================
// OTHER FUNCTIONS
//=============================================================================
//============================================================================

// Fill Variables
   
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
  
  std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
  std::cout<< "error_bands.size() = " << error_bands.size()<<std::endl;
  std::cout<<"Number of Universes set is = "<<  MinervaUniverse::GetNFluxUniverses()<<std::endl;
  
  std::vector<double> ThetaMuBin,Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin,xbinBrian;
  std::vector<double> x09bin, xfinebin;
  std::vector<double> ANNPlaneProbBin;
  std::vector<double> vtxxbin, vtxybin, vtxzbin;
  std::vector<double> planeDNNbin; 
  std::vector<double> pTbin, pZbin;
  
  if (doDIS){
    Enubin = binsDef->GetDISBins("Enu"); 
    Emubin = binsDef->GetDISBins("Emu"); 
    Ehadbin = binsDef->GetDISBins("Ehad");
    Q2bin = binsDef->GetDISBins("Q2");
    Wbin = binsDef->GetDISBins("W");
    xbin = binsDef->GetDISBins("x");
    ybin = binsDef->GetDISBins("y");
    ThetaMuBin = binsDef->GetDISBins("ThetaMu");
  }
  else{
    Enubin = binsDef->GetEnergyBins("Enu"); 
    Emubin = binsDef->GetEnergyBins("Emu"); 
    Ehadbin = binsDef->GetEnergyBins("Ehad");
    Q2bin = binsDef->GetEnergyBins("Q2");
    Wbin = binsDef->GetEnergyBins("W");
    xbin = binsDef->GetEnergyBins("x");
    x09bin = binsDef->GetEnergyBins("x09");
    xfinebin = binsDef->GetEnergyBins("xfine");
    xbinBrian    = binsDef->GetEnergyBins("xBrian");
    ybin = binsDef->GetEnergyBins("y");
    ThetaMuBin = binsDef->GetEnergyBins("ThetaMu");
    ANNPlaneProbBin = binsDef->GetEnergyBins("ANNPlaneProb");
    vtxxbin = binsDef->GetEnergyBins("vtxx"); 
    vtxybin = binsDef->GetEnergyBins("vtxy");
    vtxzbin = binsDef->GetEnergyBins("vtxz");
    planeDNNbin = binsDef->GetEnergyBins("planeDNN");
    pTbin = binsDef->GetEnergyBins("muonPt"); 
    pZbin = binsDef->GetEnergyBins("muonPz"); 
  }
  //Q2bin = binsDef->GetSidebandBins("Q2");
  //Wbin = binsDef->GetSidebandBins("W");

  // 1D Variables
  Var* thetaMu = new Var("ThetamuDeg", "ThetamuDeg", ThetaMuBin, &CVUniverse::GetThetamuDeg, &CVUniverse::GetThetamuTrueDeg);
  Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
  Var* ehad = new Var("Ehad", "Ehad (GeV)", Ehadbin, &CVUniverse::GetEhadGeV, &CVUniverse::GetEhadTrueGeV);
  Var* Q2 = new Var("Q2", "Q2 (GeV^2)", Q2bin, &CVUniverse::GetQ2RecoGeV, &CVUniverse::GetQ2TrueGeV);
  Var* W = new Var("W", "W (GeV)", Wbin, &CVUniverse::GetWRecoGeV, &CVUniverse::GetWTrueGeV);
  Var* emu = new Var("Emu", "Emu (GeV)", Emubin, &CVUniverse::GetMuonEGeV, &CVUniverse::GetMuonETrueGeV);
  Var* x = new Var("x", "x", xbin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
  Var* x09 = new Var("x09", "x09", x09bin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
  Var* xfine = new Var("xfine", "xfine", xfinebin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
  Var* xBrian = new Var("xBrian", "xBrian", xbinBrian, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
  Var* y = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);

  Var* pTmu = new Var("pTmu", "pTmu", pTbin, &CVUniverse::GetMuonPt, &CVUniverse::GetlepPtTrue);
  Var* pZmu = new Var("pZmu", "pZmu", pZbin, &CVUniverse::GetMuonPz, &CVUniverse::GetlepPzTrue);
  Var* vtxx = new Var("vtxx", "Vertex X", vtxxbin, &CVUniverse::GetVertexXMy, &CVUniverse::GetVertexXTrueMy);
  Var* vtxy = new Var("vtxy", "Vertex Y", vtxybin, &CVUniverse::GetVertexYMy, &CVUniverse::GetVertexYTrueMy);
  Var* vtxz = new Var("vtxz", "Vertex Z", vtxzbin, &CVUniverse::GetVertexZMy, &CVUniverse::GetVertexZTrueMy);
  
  Var *ANNPlaneProb = new Var("ANNPlaneProb", "ANNPlaneProb", ANNPlaneProbBin, &CVUniverse::GetANNPlaneProb, &CVUniverse::GetANNPlaneProb);
  Var* planeDNN = new Var("planeDNN", "planeDNN", planeDNNbin, &CVUniverse::GetplaneDNNReco, &CVUniverse::GetplaneDNNTrue);

  //variables = {emu, ehad, enu, thetaMu, x, x09, xfine, xBrian, y, Q2, W, vtxx, vtxy, vtxz, ANNPlaneProb, planeDNN, pTmu, pZmu}; //{enu,ehad}; 
  variables = {enu,  x,vtxz, planeDNN, Q2, pTmu, pZmu, thetaMu}; //{enu,ehad}; 

   // 2D Variables 
  Var2D* pZmu_pTmu = new Var2D(*pZmu, *pTmu);
  Var2D* W_Q2 = new Var2D(*W, *Q2);
  Var2D* enu_ehad = new Var2D(*enu, *ehad);
  Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
  Var2D* x_y = new Var2D(*x, *y);  // y var
  Var2D* x_Q2 = new Var2D(*x, *Q2);  // y var
  
  variables2d = {pZmu_pTmu };
   
  for (auto v : variables2d) v->InitializeAllHistograms(error_bands);
  for (auto v : variables) v->InitializeAllHistograms(error_bands);

  int reco0=0;
  int reco1=0;
  int reco2=0;
  int reco3=0;
  int reco4=0; 
  int reco5=0; 
  int reco6=0; 
  int Signal = 0;
  int Bkg = 0;
  int NotEmu = 0;
  int NotTracker = 0 ;
  int NC = 0;
  int WrongSign = 0;

  //int neutralcur = 0;
  //int wrongsign = 0;
  
  CVUniverse *dataverse = new CVUniverse(chain,0);
  

    //=========================================
    // Entry Loop
    //=========================================

  std::cout<<"# of entries = "<<chain->GetEntries()<<std::endl;
  for(int i=0; i<chain->GetEntries(); ++i){
    if(i%5000==0) std::cout << (i/1000) << "k " << std::endl;
      //=========================================
      // For every systematic, loop over the universes, and fill the
      // appropriate histogram in the MnvH1D
      //=========================================
      if(isMC){     
        for (auto band : error_bands){
          std::vector<CVUniverse*> error_band_universes = band.second;
          for (auto universe : error_band_universes){
            // Tell the Event which entry in the TChain it's looking at
            //=========================================
            // CUTS in each universe
            //========================================
            universe->SetEntry(i);
            reco0++;

            if( universe->GetInt("muon_corrected_p") == -999 ) continue; // additional cut to get rid of an issue
            if(!cutter->PassReco(universe,helicity)) continue;
            reco1++;
            
            if(!cutter->Water(universe)) continue;
            reco2++;
            
            if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
            //if( universe->GetVecElem("ANN_plane_probs",0) < 0.2 ) continue;	   
            reco3++;

            for (auto v : variables2d){
              if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(universe)) continue;
              if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(universe)) continue;
              v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
            }

            for (auto v : variables){
              if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(universe)) continue;
              if( v->GetName()=="Enu") reco4++;
            
              if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(universe))continue;
              if (v->GetName()=="Enu") reco5++;
      
              v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              //v->m_selected_mc_sb.GetComponentHist("MC")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
            }
          } // End band's universe loop
        }// End Band loop
      }

      else{

      dataverse->SetEntry(i);
      reco0++;
  
      if( dataverse->GetInt("muon_corrected_p") == -999 ) continue; // additional cut to get rid of an issue
      if(!cutter->PassReco(dataverse,helicity)) continue;
      reco1++;
      
      if(!cutter->Water(dataverse)) continue;
      reco2++;

      if( dataverse->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
      reco3++;
            

      for (auto v : variables2d){
        if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(dataverse)) continue;
        if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(dataverse)) continue;	     
    
        v->m_selected_data_reco.hist->Fill(v->GetRecoValueX(*dataverse), v->GetRecoValueY(*dataverse));
      }
    
      for (auto v : variables){
        if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(dataverse)) continue; 
        if (v->GetName()=="Enu") reco4++;
        if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(dataverse))continue;
        if (v->GetName()=="Enu") reco5++;
      
        v->m_selected_data_reco.hist->Fill(v->GetRecoValue(*dataverse, 0));	      }       
    }
  }//End entries loop



  for(auto band : error_bands){
    std::vector<CVUniverse*> band_universes = band.second;
    for(unsigned int i_universe = 0; i_universe < band_universes.size(); ++i_universe)
    delete band_universes[i_universe];
  } 
    
  delete dataverse;

  int univ_norm = error_bands.size();

  // Printing summary
  std::cout << "**********************************" << std::endl;
  std::cout << "Printing the ";
    isMC? std::cout << "MC ": std::cout << "Data ";
  std::cout << "Summary " << std::endl;
  std::cout << "No cuts = " << reco0 << std::endl;
  std::cout << "Reco Cut = " << reco1 << std::endl;
  std::cout << "Water Cut = " << reco2 << std::endl;
  std::cout << "Plane prob. cut = " << reco3 << std::endl;
  std::cout << "Muon Energy cut  = "<< reco4 << std::endl;
  std::cout << "Muon theta cut  = " << reco5 << std::endl;
  std::cout << "**********************************" << std::endl;
  
  //return variables;
}
//=============================================================================
