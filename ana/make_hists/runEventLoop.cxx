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
#include "../include/VariableRunBkg.h" 
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
  int targetID = atoi(argv[2]);
  int targetZ = atoi(argv[3]);
  const string playlist= argv[4];
   
  bool doDIS=false;

  // MasterAnaDev tuples?
  //const std::string mc_file_list("../include/playlists/OfficialMAD_minervaME6B_MCNukeOnly_merged.txt");
  //const std::string data_file_list("../include/playlists/OfficialMAD_minervaME6B_Data_merged.txt");
  //const std::string reco_tree_name("MasterAnaDev");

  // NukeCC Tuples ?
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

  //Works good for the grid submission
  //TString histFileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ );

  TFile fout(dir.Append(histFileName),"RECREATE");	
   
  // MC 
  std::cout << "Processing MC and filling histograms" << std::endl;

  FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true,targetID, targetZ, plist_string,doDIS);     
  for (auto v : variablesMC) {
    v->m_selected_mc_reco.SyncCVHistos();
    // Interaction types
    v->m_selected_mc_reco_QE.SyncCVHistos();
    v->m_selected_mc_reco_RES.SyncCVHistos();
    v->m_selected_mc_reco_DIS.SyncCVHistos();
    v->m_selected_mc_reco_2p2h.SyncCVHistos();
    v->m_selected_mc_reco_OtherIT.SyncCVHistos();
    // Background breakdown
    v->m_selected_mc_reco_bkg.SyncCVHistos();
    v->m_selected_mc_reco_signal.SyncCVHistos();
    v->m_selected_mc_USplastic.SyncCVHistos();
    v->m_selected_mc_DSplastic.SyncCVHistos();
    v->m_selected_mc_other.SyncCVHistos();
    v->m_selected_mc_WrongSign.SyncCVHistos();
    v->m_selected_mc_NC.SyncCVHistos();
    v->m_selected_mc_NotEmu.SyncCVHistos();
  }

  for (auto v : variables2DMC) {
    v->m_selected_mc_reco.SyncCVHistos();
    // Interaction types
    v->m_selected_mc_reco_QE.SyncCVHistos();
    v->m_selected_mc_reco_RES.SyncCVHistos();
    v->m_selected_mc_reco_DIS.SyncCVHistos();
    v->m_selected_mc_reco_2p2h.SyncCVHistos();
    v->m_selected_mc_reco_OtherIT.SyncCVHistos();
    // Background breakdown
    v->m_selected_mc_reco_bkg.SyncCVHistos();
    v->m_selected_mc_reco_signal.SyncCVHistos();
    v->m_selected_mc_USplastic.SyncCVHistos();
    v->m_selected_mc_DSplastic.SyncCVHistos();
    v->m_selected_mc_other.SyncCVHistos();
    v->m_selected_mc_WrongSign.SyncCVHistos();
    v->m_selected_mc_NC.SyncCVHistos();
    //v->m_selected_mc_NotEmu.SyncCVHistos();
  }
   
  // DATA
  std::cout << "Processing Data and filling histograms" << std::endl;

  FillVariable(chainData, helicity, utils, cutter,binsDef,variablesData,variables2DData,false,targetID, targetZ, plist_string,doDIS);
  for (auto v : variablesData) v->m_selected_data_reco.SyncCVHistos();
  for (auto v : variablesData) v->m_selected_data_reco_sb.SyncCVHistos();
  for (auto v : variables2DData) v->m_selected_data_reco.SyncCVHistos();

  // WRITE HISTOGRAMS TO FILE

  // 1D variables
  for (auto v : variablesMC) {
    v->WriteAllHistogramsToFile(fout, true);
  }

  for (auto v : variablesData) {
    v->WriteAllHistogramsToFile(fout, false);
  }

  // 1D Plotting
  //for(int i=0; i < variablesMC.size();i++){
    //PlotCVAndError(variablesData[i]->m_selected_data_reco.hist,variablesMC[i]->m_selected_mc_reco.hist,variablesMC[i]->GetName(),MCscale); 
    //PlotErrorSummary(variablesMC[i]->m_selected_mc_reco.hist, variablesMC[i]->GetName());
    //PlotStacked(variablesData[i]->m_selected_data_reco_sb.hist,variablesMC[i]->m_selected_mc_sb.GetHistArray(),MCscale, variablesMC[i]->m_selected_mc_sb.GetName(), variablesMC[i]->m_selected_mc_sb.GetName());
  //}//End 1D plotting 
 
  // 2D Variables
  for (auto v : variables2DMC) {
    v->WriteAllHistogramsToFile(fout,true);
  }

  for (auto v : variables2DData) {
    v->WriteAllHistogramsToFile(fout,false);
  }
 
  // 2D Plotting
  //for(int i=0; i< variables2DMC.size();i++){
    //Plot2D(variables2DMC[i]->m_selected_mc_reco.hist, variables2DMC[i]->GetName(), variables2DMC[i]->GetNameX(), variables2DMC[i]->GetNameY());
    //Plot2D(variables2DData[i]->m_selected_data_reco.hist, variables2DData[i]->GetName(), variables2DData[i]->GetNameX(),variables2DData[i]->GetNameY());   
  //}//End 2D plotting
	
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
//=============================================================================

// Fill Variables
   
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
  
  std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
  
  std::vector<double> ThetaMuBin,Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin,xbinBrian;
  std::vector<double> x09bin, xfinebin;
  std::vector<double> ANNPlaneProbBin;
  std::vector<double> vtxzbin;
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
  Var* vtxz = new Var("vtxz", "Vertex Z", vtxzbin, &CVUniverse::GetVertexZMy, &CVUniverse::GetVertexZTrueMy);
  Var *ANNPlaneProb = new Var("ANNPlaneProb", "ANNPlaneProb", ANNPlaneProbBin, &CVUniverse::GetANNPlaneProb, &CVUniverse::GetANNPlaneProb);
  Var* planeDNN = new Var("planeDNN", "planeDNN", planeDNNbin, &CVUniverse::GetplaneDNNReco, &CVUniverse::GetplaneDNNTrue);

  variables = {emu, enu, x, vtxz, planeDNN, pTmu, pZmu, thetaMu}; //{enu,ehad}; 

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

  std::vector<int> targetIDs;

  if( targetID > 10){
    targetIDs.push_back(14);
    targetIDs.push_back(24);
    targetIDs.push_back(34);
    targetIDs.push_back(44);
    targetIDs.push_back(54);
    targetIDs.push_back(64);
    targetIDs.push_back(74);
    targetIDs.push_back(84);
    targetIDs.push_back(94);
  }
  
  else{
    if(targetZ==26){
      targetIDs.push_back(1);
      targetIDs.push_back(2);
      targetIDs.push_back(3);
      targetIDs.push_back(5);
    }
    if(targetZ==82){
      targetIDs.push_back(1);
      targetIDs.push_back(2);
      targetIDs.push_back(3);
      targetIDs.push_back(4);
      targetIDs.push_back(5);
    }
    if(targetZ==6){
      targetIDs.push_back(3);
    }
  }

  int reco0=0;
  int reco1=0;
  int reco2=0;
  int reco3=0;
  int reco4=0; 
  int reco5=0; 
  int reco6=0; 
  int Signal = 0;
  int NotEmu = 0;
  int WrongMaterialOrTarget = 0;
  int plastic = 0;
  int Bkg = 0;
  int WrongSign = 0; // just CC antinu
  int NC = 0; // both nu and antinu
  int otherneutrinotype = 0;
  int US=0;
  int DS=0;
  int other=0;

  int QE = 0;
  int RES = 0;
  int DIS = 0;
  int npnh = 0;
  int otherIT = 0 ;

  int Signal2d = 0;
  //int NotEmu2d = 0;
  int WrongMaterialOrTarget2d = 0;
  int plastic2d = 0;
  int Bkg2d = 0;
  int WrongSign2d = 0; // just CC antinu
  int NC2d = 0; // both nu and antinu
  int otherneutrinotype2d = 0;
  int US2d=0;
  int DS2d=0;
  int other2d=0;
  
  CVUniverse *dataverse = new CVUniverse(chain,0);
    
  //=========================================
  // Targets combining Loop
  //=========================================

  //for(int t = 0; t < targetIDs.size(); t++){

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
              
              if(!cutter->IsInMaterial(universe,targetID,targetZ, /*anyTrakerMod*/false)) continue;
              //if(!cutter->IsInMaterial(universe,targetIDs[t],targetZ, /*anyTrakerMod*/false)) continue;
              reco2++;
              
              if(targetID<10 && universe->GetInt("MasterAnaDev_ANN_targetID") != targetID) continue;
              //if(targetIDs[t]<10 && universe->GetInt("NukeCC_targetID") != targetIDs[t]) continue;
              reco3++;
              
              if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
              //if( universe->GetVecElem("ANN_plane_probs",0) < 0.2 ) continue;	   
              reco4++;

              // 2D SELECTION

              for (auto v : variables2d){
                //if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(universe)) continue;
                if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(universe)) continue;	     
                
                v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 

                // Interaction breakdown
                if (universe->GetInt("mc_intType") == 1){ // QE
                  v->m_selected_mc_reco_QE.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                }
                else if (universe->GetInt("mc_intType") == 2){ // RES
                  v->m_selected_mc_reco_RES.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                }
                else if (universe->GetInt("mc_intType") == 3){ // DIS 
                  v->m_selected_mc_reco_DIS.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                }
                else if (universe->GetInt("mc_intType") == 8){ // 2p2h
                  v->m_selected_mc_reco_2p2h.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                }
                else{ // other
                  v->m_selected_mc_reco_OtherIT.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                }

                // Background breakdown
                // Signal
                if( 1 == universe->GetInt("mc_current") &&  -14 == universe->GetInt("mc_incoming") ){
                  if(cutter->IsInTrueMaterial(universe,targetID, targetZ,false)) { // true fiducial z distance
                    v->m_selected_mc_reco_signal.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                    if (v->GetName()=="pZmu_pTmu")  Signal2d++;
                  }
                  else{ 
                    //Background: different material !
                    v->m_selected_mc_reco_bkg.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                    if (v->GetName()=="pZmu_pTmu")  Bkg2d++;
                    
                    if(universe->GetInt("truth_targetID") == 0 ){ // targetID is plastic == same as picking everything that is not target 1,2,3,4,5 or water
                      v->m_selected_mc_plastic.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                      if (v->GetName()=="pZmu_pTmu")  plastic2d++;
                      
                      // DS
                      if(universe->GetVecElem("ANN_vtx", 2) < universe->GetVecElem("mc_vtx", 2)){ // DS
                        v->m_selected_mc_DSplastic.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                        if (v->GetName()=="pZmu_pTmu")  DS2d++;
                      }
                      // US
                      else if(universe->GetVecElem("ANN_vtx", 2) > universe->GetVecElem("mc_vtx", 2)){ // US
                      v->m_selected_mc_USplastic.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                      if (v->GetName()=="pZmu_pTmu")  US2d++;
                      }
                    }
                    else{
                      v->m_selected_mc_other.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                      if (v->GetName()=="pZmu_pTmu")  other2d++;
                    }
                  }
                }
                // Background: wrong sign events and neutral current events
                else { 
                  v->m_selected_mc_reco_bkg.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                  if (v->GetName()=="pZmu_pTmu")  Bkg2d++;

                  if( 1 != universe->GetInt("mc_current")){ // NC neutrino and antineutrino
                    v->m_selected_mc_NC.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                    if (v->GetName()=="pZmu_pTmu")  NC2d++;
                  }

                  else if( -14 != universe->GetInt("mc_incoming") ){ // CC neutrino only
                    v->m_selected_mc_WrongSign.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
                    if (v->GetName()=="pZmu_pTmu")  WrongSign2d++;
                  }
                }
	            }

              // 1D SELECTION
	            for (auto v : variables){
	              if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(universe)) continue;
	              if( v->GetName()=="Enu") reco5++;
             
	              if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(universe))continue;
	              if (v->GetName()=="Enu") reco6++;
	     
                v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                //v->m_selected_mc_sb.GetComponentHist("MC")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());

                // Interaction breakdown
                if (universe->GetInt("mc_intType") == 1){ // QE
                  v->m_selected_mc_reco_QE.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  QE++;
                }
                else if (universe->GetInt("mc_intType") == 2){ // RES
                  v->m_selected_mc_reco_RES.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  RES++;
                }
                else if (universe->GetInt("mc_intType") == 3){ // DIS 
                  v->m_selected_mc_reco_DIS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  DIS++;
                }
                else if (universe->GetInt("mc_intType") == 8){ // 2p2h
                  v->m_selected_mc_reco_2p2h.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  npnh++;
                }
                else{ // other
                  v->m_selected_mc_reco_OtherIT.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  otherIT++;
                }

                // Background breakdown

                // Signal
                if( 1 == universe->GetInt("mc_current") &&  -14 == universe->GetInt("mc_incoming") ){
                  if(cutter->IsInTrueMaterial(universe,targetID, targetZ,false)) { // true fiducial z distance
                    if( v->GetName()!="Emu")  if(cutter->PassTrueMuEnergyCut(universe)){
                      v->m_selected_mc_reco_signal.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                      if (v->GetName()=="Enu")  Signal++;
                    }
                    else{
                      // Background: out of muon energy range !
                      v->m_selected_mc_reco_bkg.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                      if (v->GetName()=="Enu") Bkg++;
                      v->m_selected_mc_NotEmu.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                      if (v->GetName()=="Enu")  NotEmu++;
                    }
                  }
                  else{ 
                    //Background: different material !
                    v->m_selected_mc_reco_bkg.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                    if (v->GetName()=="Enu") Bkg++;
                    if (v->GetName()=="Enu") WrongMaterialOrTarget++; 
                    
                    if(universe->GetInt("truth_targetID") == 0 ){ // targetID is plastic == same as picking everything that is not target 1,2,3,4,5 or water
                      v->m_selected_mc_plastic.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                      if (v->GetName()=="Enu") plastic++;
                      
                      // DS
                      if(universe->GetVecElem("ANN_vtx", 2) < universe->GetVecElem("mc_vtx", 2)){ // DS
                        v->m_selected_mc_DSplastic.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                        if (v->GetName()=="Enu") DS++;
                      }
                      // US
                      else if(universe->GetVecElem("ANN_vtx", 2) > universe->GetVecElem("mc_vtx", 2)){ // US
                      v->m_selected_mc_USplastic.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                      if (v->GetName()=="Enu") US++;
                      }
                    }
                    else{
                      v->m_selected_mc_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                      if (v->GetName()=="Enu") other++;
                    }
                  }
                }
                // Background: wrong sign events and neutral current events
                else { 
                  v->m_selected_mc_reco_bkg.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  if (v->GetName()=="Enu") Bkg++;

                  if( 1 != universe->GetInt("mc_current")){ // NC neutrino and antineutrino
                    v->m_selected_mc_NC.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                    if (v->GetName()=="Enu") NC++;
                    if (-14 != universe->GetInt("mc_incoming") && 14 != universe->GetInt("mc_incoming")){
                      // if not muon antineutrino or muon neutrino -> print
                      std::cout << universe->GetInt("mc_incoming") << std::endl;
                      otherneutrinotype++;
                    }
                  }

                  else if( -14 != universe->GetInt("mc_incoming") ){ // CC neutrino only
                    v->m_selected_mc_WrongSign.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                    if (v->GetName()=="Enu") WrongSign++;
                  }
                }
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
        
        if(!cutter->IsInMaterial(dataverse,targetID,targetZ, false)) continue;
        //if(!cutter->IsInMaterial(dataverse,targetIDs[t],targetZ, false)) continue;
        reco2++;
        
        if(targetID<10 && dataverse->GetInt("MasterAnaDev_ANN_targetID") != targetID) continue;
        //if(targetIDs[t]<10 && dataverse->GetInt("NukeCC_targetID") != targetIDs[t]) continue;
        reco3++;
	 
        if( dataverse->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
        //if( dataverse->GetVecElem("ANN_plane_probs",0) < 0.2 ) continue;	    
        reco4++;
              
        /* GLOBAL VS LOCAL CUT*/
        //if(!cutter->PassMuEnergyCut(dataverse)) continue;
        //reco5++;
          
        //if(!cutter->PassThetaCut(dataverse)) continue; 
        //reco6++;

 	      for (auto v : variables2d){
	        if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(dataverse)) continue;
	        if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(dataverse)) continue;	     
	    
	        v->m_selected_data_reco.hist->Fill(v->GetRecoValueX(*dataverse), v->GetRecoValueY(*dataverse));
	      }
	   
	      for (auto v : variables){
	        if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(dataverse)) continue; 
	        if (v->GetName()=="Enu") reco5++;
	        if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(dataverse))continue;
	        if (v->GetName()=="Enu") reco6++;
	     
	        v->m_selected_data_reco.hist->Fill(v->GetRecoValue(*dataverse, 0));
	        v->m_selected_data_reco_sb.hist->Fill(v->GetRecoValue(*dataverse, 0));
	      }       
      }

    }//End entries loop
  //}//Target loop closed


  for(auto band : error_bands){
    std::vector<CVUniverse*> band_universes = band.second;
    for(unsigned int i_universe = 0; i_universe < band_universes.size(); ++i_universe)
    delete band_universes[i_universe];
  } 
    
  delete dataverse;

  // Printing summary
  std::cout << "**********************************" << std::endl;
  std::cout << "Printing the ";
    isMC? std::cout << "MC ": std::cout << "Data ";
  std::cout << "Summary " << std::endl;
  std::cout << "No cuts = " << reco0 << std::endl;
  std::cout << "Reco Cut = " << reco1 << std::endl;
  std::cout << "Material Cut = " << reco2 << std::endl;
  std::cout << "TargetID Cuts = " << reco3 << std::endl;
  std::cout << "Plane prob. cut = " << reco4 << std::endl;
  //std::cout<<" Neutral Current cuts = " << neutralcur << std::endl;
  //std::cout<<" Wrong Sign Cuts = " << wrongsign << std::endl;
  //std::cout<<" In Hexagon Cuts = "<<reco4<<std::endl;
  //std::cout<<" Muon Kinematics Cuts = "<<reco4<<std::endl;
  std::cout << "Muon Energy cut  = "<< reco5 << std::endl;
  std::cout << "Muon theta cut  = " << reco6 << std::endl;
  std::cout << "**********************************" << std::endl;
  if (isMC) {
    std::cout<<"Interaction Type "<<std::endl;
    std::cout << "QE = " << QE << std::endl;
    std::cout << "RES = " << RES << std::endl;
    std::cout << "DIS= " << DIS << std::endl;
    std::cout << "2p2h = " << npnh << std::endl;
    std::cout << "Other = " << otherIT << std::endl;
    std::cout << "**********************************" << std::endl;

    std::cout << "Signal and Background " << std::endl;
    std::cout << "Signal  = "<< Signal<< std::endl;
    std::cout << "All background = "<< Bkg << std::endl;
    std::cout << "Wrong target or material  = "<< WrongMaterialOrTarget << std::endl;
    std::cout << "Out of that is NOT Plastic = " << other << std::endl;
    std::cout << "Out of that is Plastic = " << plastic << std::endl;
    std::cout << "Out of that is US plastic  = " << US << std::endl;
    std::cout << "Out of that is DS plastic  = " << DS << std::endl;
    std::cout << "Neutral current (nu+antinu) = "<< NC<< std::endl;
    std::cout << "Other neutrino types = " << otherneutrinotype << std::endl;
    std::cout << "Wrong sign (CC) = "<< WrongSign << std::endl;
    std::cout << "Not muon energy = "<< NotEmu << std::endl;
    std::cout << "**********************************" << std::endl;

    std::cout << "Signal and Background 2D " << std::endl;
    std::cout << "Signal  = "<< Signal2d<< std::endl;
    std::cout << "All background = "<< Bkg << std::endl;
    std::cout << "Wrong target or material  = "<< WrongMaterialOrTarget2d << std::endl;
    std::cout << "Out of that is NOT Plastic = " << other2d << std::endl;
    std::cout << "Out of that is Plastic = " << plastic2d << std::endl;
    std::cout << "Out of that is US plastic  = " << US2d << std::endl;
    std::cout << "Out of that is DS plastic  = " << DS2d << std::endl;
    std::cout << "Neutral current (nu+antinu) = "<< NC2d<< std::endl;
    std::cout << "Other neutrino types = " << otherneutrinotype2d << std::endl;
    std::cout << "Wrong sign (CC) = "<< WrongSign2d << std::endl;
    //std::cout << "Not muon energy = "<< NotEmu2d << std::endl;
    std::cout << "**********************************" << std::endl;
  }
  //return variables;
}
//=============================================================================

