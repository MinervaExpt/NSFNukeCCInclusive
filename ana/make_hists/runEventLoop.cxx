//==============================================================================
// Loop entries, make cuts, fill histograms.
// * Uses the New Systematics Framework and "Universe" objects.
// * loop universes, make cuts and fill histograms with the correct lateral
// shifts and weights for each universe.
// * TChain --> PlotUtils::ChainWrapper.
// * MnvHXD --> PlotUtils::HistWrapper.
// * Genie, flux, non-resonant pion, and some detector systematics calculated.
//==============================================================================

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
//#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
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
   
  bool doDIS=false;

  // MasterAnaDev tuples?
  //const std::string mc_file_list("../include/playlists/OfficialMAD_minervaME6B_MCNukeOnly_merged.txt");
  //const std::string data_file_list("../include/playlists/OfficialMAD_minervaME6B_Data_merged.txt");
  //const std::string reco_tree_name("MasterAnaDev");

  // NukeCC Tuples ?
  const std::string mc_file_list("../include/playlists/NukeCC_MC_minervame6A_MuonKludged.txt");
  const std::string data_file_list("../include/playlists/NukeCC_Data_minervame6A_MuonKludged.txt");
  const std::string reco_tree_name("NukeCC");
  
  const std::string plist_string("minervame6A");
  const bool wants_truth = false;
  const bool is_grid = false;

  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth, is_grid);

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

  TString histFileName = utils->GetHistFileName( "EventSelection_ML_ME6A", FileType::kAny, targetID, targetZ, helicity ); 

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

// Error Bands
std::map<std::string, std::vector<CVUniverse*> > GetErrorBands(PlotUtils::ChainWrapper* chain) {
  
  typedef std::map<std::string, std::vector<CVUniverse*> > SystMap;
  
  // return map
  SystMap error_bands;

  // CV 
  error_bands[std::string("CV")].push_back(new CVUniverse(chain));

  if(RunCodeWithSystematics) {
  
    //========================================================================
    // FLUX
    //========================================================================
    int n_flux_universes = 100 ;
    SystMap flux_systematics = 
    PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain,n_flux_universes);
  
    error_bands.insert(flux_systematics.begin(), flux_systematics.end());

    //========================================================================
    // GENIE
    //========================================================================
    SystMap genie_systematics = PlotUtils::GetGenieSystematicsMap<CVUniverse>(chain);
    // change that true to a switch on do_nonrespi_tune
    error_bands.insert(genie_systematics.begin(), genie_systematics.end());
    // includes standard Genie + Rvx1pi (final state pion normalization) + FaCCQE

    //========================================================================
    // MnvTunes
    //======================================================================== 
  
    // RPA
    SystMap RPA_systematics = PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
    error_bands.insert(RPA_systematics.begin(),RPA_systematics.end());

    // 2p2h
    SystMap _2p2h_systematics = PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
    error_bands.insert(_2p2h_systematics.begin(),_2p2h_systematics.end());

  
    //========================================================================
    //  Muons
    //========================================================================

    // Muon momentum
    // MINOS
    SystMap muonminosP_systematics = PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
    error_bands.insert(muonminosP_systematics.begin(), muonminosP_systematics.end());
  
    // Muon momentum
    // MINERvA
    SystMap muonminervaP_systematics = PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
    error_bands.insert(muonminervaP_systematics.begin(), muonminervaP_systematics.end());

    //Muon Momentum Resolution
    SystMap muonP_resolutions = PlotUtils::GetMuonResolutionSystematicsMap<CVUniverse>(chain);
    error_bands.insert(muonP_resolutions.begin(),muonP_resolutions.end());
  
    //MINOS Efficiency
    SystMap minos_efficiency = PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
    error_bands.insert(minos_efficiency.begin(),minos_efficiency.end());

    //========================================================================
    // Detector
    //========================================================================
  
    SystMap angle_systematics = PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);
    ///////////////////////NSFDefaults::beamThetaX_Err,NSFDefaults::beamThetaY_Err);
    error_bands.insert(angle_systematics.begin(), angle_systematics.end());


    // Target Mass Systematics
    SystMap targetMass_systematics = PlotUtils::GetTargetMassSystematicsMap<CVUniverse>(chain);
    error_bands.insert(targetMass_systematics.begin(), targetMass_systematics.end());

  }

  return error_bands;
}


// Fill Variables
   
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
  
  std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
  
  std::vector<double> ThetaMuBin,Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin;
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
  Var* thetaMu = new Var("GetThetamuDeg", "GetThetamuDeg (Degree)", ThetaMuBin, &CVUniverse::GetThetamuDeg, &CVUniverse::GetThetamuTrueDeg);
  Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
  Var* ehad = new Var("Ehad", "Ehad (GeV)", Ehadbin, &CVUniverse::GetEhadGeV, &CVUniverse::GetEhadTrueGeV);
  Var* Q2 = new Var("Q2", "Q2 (GeV^2)", Q2bin, &CVUniverse::GetQ2RecoGeV, &CVUniverse::GetQ2TrueGeV);
  Var* W = new Var("W", "W (GeV)", Wbin, &CVUniverse::GetWRecoGeV, &CVUniverse::GetWTrueGeV);
  Var* emu = new Var("Emu", "Emu (GeV)", Emubin, &CVUniverse::GetMuonEGeV, &CVUniverse::GetMuonETrueGeV);
  Var* x = new Var("x", "x", xbin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
  Var* y = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);

  Var* pTmu = new Var("pTmu", "pTmu", pTbin, &CVUniverse::GetMuonPt, &CVUniverse::GetlepPtTrue);
  Var* pZmu = new Var("pZmu", "pZmu", pZbin, &CVUniverse::GetMuonPz, &CVUniverse::GetlepPzTrue);
  Var* vtxz = new Var("vtxz", "Vertex Z", vtxzbin, &CVUniverse::GetVertexZMy, &CVUniverse::GetVertexZTrueMy);
  Var *ANNPlaneProb = new Var("ANNPlaneProb", "ANNPlaneProb", ANNPlaneProbBin, &CVUniverse::GetANNPlaneProb, &CVUniverse::GetANNPlaneProb);
  Var* planeDNN = new Var("planeDNN", "planeDNN", planeDNNbin, &CVUniverse::GetplaneDNNReco, &CVUniverse::GetplaneDNNTrue);

  variables = {emu, ehad, enu, thetaMu, x, y, Q2, W, vtxz, ANNPlaneProb, planeDNN, pTmu, pZmu}; //{enu,ehad}; 

  // 2D Variables 
  Var2D* pTmu_pZmu = new Var2D(*pTmu, *pZmu);
  Var2D* W_Q2 = new Var2D(*W, *Q2);
  Var2D* enu_ehad = new Var2D(*enu, *ehad);
  Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
  Var2D* x_y = new Var2D(*x, *y);  // y var
  Var2D* x_Q2 = new Var2D(*x, *Q2);  // y var
  
  variables2d = {emu_ehad,enu_ehad, x_y, W_Q2, pTmu_pZmu };//{enu_ehad, Q2_W};
   
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
  //int neutralcur = 0;
  //int wrongsign = 0;
  
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

              if(!cutter->PassReco(universe,helicity)) continue;
              reco1++;
              
              if(!cutter->IsInMaterial(universe,targetID,targetZ, /*anyTrakerMod*/false)) continue;
              //if(!cutter->IsInMaterial(universe,targetIDs[t],targetZ, /*anyTrakerMod*/false)) continue;
              reco2++;
              
              if(targetID<10 && universe->GetInt("NukeCC_targetID") != targetID) continue;
              //if(targetIDs[t]<10 && universe->GetInt("NukeCC_targetID") != targetIDs[t]) continue;
              reco3++;
              
              if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
              //if( universe->GetVecElem("ANN_plane_probs",0) < 0.2 ) continue;	   
              reco4++;
              
              
              /* GLOBAL VS LOCAL CUT*/
              //if(!cutter->PassMuEnergyCut(universe)) continue;
              //reco5++;

              //if(!cutter->PassThetaCut(universe)) continue;
              //reco6++;

              //pure signal cuts on MC to estimate bkg (counts)    
              //if( 1 != universe->GetInt("mc_current") ) continue; 
              //neutralcur++;
              //if(-14 != universe->GetInt("mc_incoming") ) continue;
              //wrongsign++;

              for (auto v : variables2d){
	              if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(universe)) continue;
	              if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(universe)) continue;
	              v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
	            }
	
	            for (auto v : variables){
	              if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(universe)) continue;
	              if( v->GetName()=="Enu") reco5++;
             
	              if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(universe))continue;
	              if (v->GetName()=="Enu") reco6++;
	     
                v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                //v->m_selected_mc_sb.GetComponentHist("MC")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              
                /*
                // Signal
                if( 1 == universe->GetInt("mc_current") &&  -14 == universe->GetInt("mc_incoming") ){
                  v->m_selected_mc_sb.GetComponentHist("Sig")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                }
                // Background
                else { 
                  v->m_selected_mc_sb.GetComponentHist("Bkg")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  if( 1 != universe->GetInt("mc_current")){
                  v->m_selected_mc_sb.GetComponentHist("NC")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  }
                  if( -14 != universe->GetInt("mc_incoming") ){
                  v->m_selected_mc_sb.GetComponentHist("WRGsign")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  }
                }*/
              }
            } // End band's universe loop
          }// End Band loop
        }

        else{

        dataverse->SetEntry(i);
        reco0++;
	 
	      if(!cutter->PassReco(dataverse,helicity)) continue;
	      reco1++;
        
        if(!cutter->IsInMaterial(dataverse,targetID,targetZ, false)) continue;
        //if(!cutter->IsInMaterial(dataverse,targetIDs[t],targetZ, false)) continue;
        reco2++;
        
        if(targetID<10 && dataverse->GetInt("NukeCC_targetID") != targetID) continue;
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
  
  //return variables;
}
//=============================================================================

