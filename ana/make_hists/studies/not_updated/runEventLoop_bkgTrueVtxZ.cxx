// HOW MANY WRONG SIGN EVENTS AND NEUTRAL CURRENT EVENTS ARE IN MY PASS RECO MC

//==============================================================================
// Loop entries, make cuts, fill histograms.
// * Uses the New Systematics Framework and "Universe" objects.
// * loop universes, make cuts and fill histograms with the correct lateral
// shifts and weights for each universe.
// * TChain --> PlotUtils::ChainWrapper.
// * MnvHXD --> PlotUtils::HistWrapper.
// * Genie, flux, non-resonant pion, and some detector systematics calculated.
//==============================================================================

#include "include/CommonIncludes.h"
#include "include/CVUniverse.h"
#include "VariableRunBkgMat.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "include/LateralSystematics.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>

#include "Cintex/Cintex.h"
#include "include/NukeCCUtilsNSF.h"
#include "include/NukeCC_Cuts.h"
#include "PlotUtils/MnvPlotter.h" 
#include "TParameter.h"

#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"

// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we 
// shield it.
#ifndef __CINT__
#include "../../include/plotting_functions.h"
#endif
#include "PlotUtils/MacroUtil.h" 
//using namespace globalV;
using namespace NUKECC_ANA;

//======================================================================

typedef VarLoop::Variable Var;
typedef Var2DLoop::Variable2D Var2D;

void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType
helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef,
std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID=1,
int targetZ=26, const string playlist="minervame1A", bool doDIS=true);

//=============================================================================
//=============================================================================
// MAIN FUNCTION
//=============================================================================
//=============================================================================

int main(int argc, char *argv[]){
  ROOT::Cintex::Cintex::Enable();
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
  const std::string mc_file_list("../../include/playlists/NukeCC_MC_minervame6A_MuonKludged.txt");
  const std::string data_file_list("../../include/playlists/NukeCC_Data_minervame6A_MuonKludged.txt");
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
  cout<<"Helicity"<<helicity<<endl;
  
  double DataPot=  util.m_data_pot; 
  double MCPot=  util.m_mc_pot;  
  double MCscale=DataPot/MCPot;
 
  std::cout << "MC Scale = " << MCscale << std::endl; 
  std::cout << "Data POT: " << DataPot << std::endl;
  std::cout << "MC POT: " << MCPot << std::endl;

  std::vector<Var*> variablesMC,variablesData; 
  std::vector<Var2D*> variables2DMC,variables2DData; 

  TString histFileName = utils->GetHistFileName( "EventSelection_Bkg_Materials_MLNukeCC_ME6A", FileType::kAny, targetID, targetZ, helicity );   
  TFile fout(dir.Append(histFileName),"RECREATE");	
   
  // MC 
  std::cout << "Processing MC and filling histograms" << std::endl;

  FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true,targetID, targetZ, plist_string,doDIS);     
  for (auto v : variablesMC) v->m_selected_mc_reco.SyncCVHistos();
  for (auto v : variablesMC) v->m_selected_mc_reco_signal.SyncCVHistos();
  for (auto v : variablesMC) v->m_selected_mc_reco_bkg.SyncCVHistos();
  for (auto v : variablesMC) v->m_selected_mc_reco_bkg_NC.SyncCVHistos();
  for (auto v : variablesMC) v->m_selected_mc_reco_bkg_WrongSign.SyncCVHistos();
  for (auto v : variablesMC) v->m_selected_mc_reco_bkg_fidcut.SyncCVHistos();
  for (auto v : variablesMC) v->m_selected_mc_reco_bkg_muon.SyncCVHistos();
  for (auto v : variablesMC) v->m_selected_mc_reco_bkg_water.SyncCVHistos();
  for (auto v : variables2DMC) v->m_selected_mc_reco.SyncCVHistos();
  for (auto v : variables2DMC) v->m_selected_mc_reco_Sg.SyncCVHistos();
  for (auto v : variables2DMC) v->m_selected_mc_reco_Fe.SyncCVHistos();
  for (auto v : variables2DMC) v->m_selected_mc_reco_Pb.SyncCVHistos();
  for (auto v : variables2DMC) v->m_selected_mc_reco_C.SyncCVHistos();
  for (auto v : variables2DMC) v->m_selected_mc_reco_other.SyncCVHistos();
  for (auto v : variables2DMC) v->m_selected_mc_reco_otherFe.SyncCVHistos();
  for (auto v : variables2DMC) v->m_selected_mc_reco_otherPb.SyncCVHistos();
  for (auto v : variables2DMC) v->m_selected_mc_reco_otherC.SyncCVHistos();
  for (auto v : variables2DMC) v->m_selected_mc_reco_Plastic.SyncCVHistos();
   
  // DATA
  std::cout << "Processing Data and filling histograms" << std::endl;

  FillVariable(chainData, helicity, utils, cutter,binsDef,variablesData,variables2DData,false,targetID, targetZ, plist_string,doDIS);
  for (auto v : variablesData) v->m_selected_data_reco.SyncCVHistos();
  //for (auto v : variablesData) v->m_selected_data_reco_sb.SyncCVHistos();
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
  //for(int i=0; i< variablesMC.size();i++){
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

// struct to feed in
struct SliceID{
  int run;
  int subrun;
  int gate;
  int slice;
};

// ARACHNE LINKS
std::string arachne(SliceID &id, const bool isData = false, const bool useRodriges = false){
  return std::string("http://minerva05.fnal.gov") + (useRodriges?"/rodriges":"") + std::string("/Arachne/arachne.html?det=") + (isData?"MV":"SIM_minerva") +
                       "&recoVer=v21r1p1&run=" + std::to_string(id.run) +
                       "&subrun=" + std::to_string(id.subrun) +
                       "&gate=" + std::to_string(id.gate + !isData) +
                       "&slice=" + std::to_string(id.slice);
};


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

  }

  return error_bands;
}


// Fill Variables
   
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
  
  std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
  
  std::vector<double> ThetaMuBin,Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin, vtxxbin, vtxybin, vtxzbin;
  std::vector<double> planeDNNbin; 

  if (doDIS){
    Enubin = binsDef->GetDISBins("Enu"); 
    Emubin = binsDef->GetDISBins("Emu"); 
    Ehadbin = binsDef->GetDISBins("Ehad");
    Q2bin = binsDef->GetDISBins("Q2");
    Wbin = binsDef->GetDISBins("W");
    xbin = binsDef->GetDISBins("x");
    ybin = binsDef->GetDISBins("y");
    ThetaMuBin = binsDef->GetDISBins("ThetaMu");

    vtxxbin = binsDef->GetDISBins("vtxx"); 
    vtxybin = binsDef->GetDISBins("vtxy");
    vtxzbin = binsDef->GetDISBins("vtxz");
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
    planeDNNbin = binsDef->GetEnergyBins("planeDNN");
    vtxxbin = binsDef->GetEnergyBins("vtxx"); 
    vtxybin = binsDef->GetEnergyBins("vtxy"); 
    vtxzbin = binsDef->GetEnergyBins("vtxz");
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
    // Get XYZ coordinates of the vertex position
  Var* vtxx = new Var("vtxx", "Vertex X", vtxxbin, &CVUniverse::GetVertexXMy, &CVUniverse::GetVertexXTrueMy);
  Var* vtxy = new Var("vtxy", "Vertex Y", vtxybin, &CVUniverse::GetVertexYMy, &CVUniverse::GetVertexYTrueMy);
  Var* vtxz = new Var("vtxz", "Vertex Z", vtxzbin, &CVUniverse::GetVertexZMy, &CVUniverse::GetVertexZTrueMy);
    Var* planeDNN = new Var("planeDNN", "planeDNN", planeDNNbin, &CVUniverse::GetplaneDNNReco, &CVUniverse::GetplaneDNNTrue);

  variables = {vtxx, vtxy, vtxz, enu, x}; 
  //{emu, ehad, enu, thetaMu, x, y}; //{enu,ehad}; 

  // 2D Variables 
  Var2D* W_Q2 = new Var2D(*W, *Q2);
  Var2D* enu_ehad = new Var2D(*enu, *ehad);
  Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
  Var2D* x_y = new Var2D(*x, *y);  // y var
  Var2D* x_Q2 = new Var2D(*x, *Q2);  // y var

  Var2D* vtxx_vtxy = new Var2D(*vtxx, *vtxy);
  
  variables2d = {vtxx_vtxy};//{enu_ehad, Q2_W};
   
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
  int neutralcur = 0;
  int wrongsign = 0;
  int signal = 0;
  int bkg = 0;
  int materialbkg = 0;
  int notmuonkin = 0;
  int notmuonangle = 0;
  int fidcut = 0;

  int fe = 0;
  int pb = 0;
  int carbon = 0;
  int ch = 0;

  int otherTarget = 0;
  int otherTargetFe = 0;
  int otherTargetPb = 0;
  int otherTargetC = 0;
  int otherTargetPlastic = 0;

  int t1fe = 0;
  int t2fe = 0;
  int t4fe = 0;
  int t5fe = 0;

  int t1pb = 0;
  int t2pb = 0;
  int t4pb = 0;
  int t5pb = 0;

  int water = 0;

  int plastH = 0;
  int plastC = 0;
  int plastO = 0;
  int plastAl = 0;
  int plastSi = 0;
  int plastTi = 0;
  int plastN = 0;
  int plastHe = 0;
  int plastCl = 0;

  
  CVUniverse *dataverse = new CVUniverse(chain,0);

  // File with Links
  std::ofstream MyFile("links.txt");
  std::vector<int> v = {};
    
  //=========================================
  // Targets combining Loop
  //=========================================

  ////for(int t = 0; t < targetIDs.size(); t++){

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
              reco4++;
              
              
              /* GLOBAL VS LOCAL CUT*/
              if(!cutter->PassMuEnergyCut(universe)) continue;
              reco5++;

              if(!cutter->PassThetaCut(universe)) continue;
              reco6++;

              // COUNTER
              //if( 1 == universe->GetInt("mc_current") &&  -14 == universe->GetInt("mc_incoming") ){
              if(cutter->PassTruth(universe,helicity) ){ // CC, antineutrino, fiducial cuts
                if (cutter->PassTrueMuEnergyCut(universe) && cutter->PassTrueThetaCut(universe)){
                  if(cutter->IsInTrueMaterial(universe,targetID, targetZ) && universe->GetInt("truth_targetID") == targetID){
                    signal++;
                  }
                  // material + other target bkg
                  else{
                    bkg++;
                    materialbkg++;
                    // bkg in t3
                    if(universe->GetInt("truth_targetID") == targetID){
                      if(cutter->IsInTrueMaterial(universe,targetID, 26)){
                        fe++ ;
                      }
                      else if(cutter->IsInTrueMaterial(universe,targetID, 82) ){
                        pb++;
                      }
                      else if(cutter->IsInTrueMaterial(universe,targetID, 6)){
                        carbon++;

                        // event comign from Carbon: LINK
                        struct SliceID Carb;
                        Carb.run = universe->GetMCRunN();
                        Carb.subrun = universe->GetMCSubRunN();
                        Carb.gate = universe->GetMCGateN();
                        Carb.slice = universe->GetMCSliceN();

                        MyFile << "Carbon Signal event in T3 # " << carbon <<std::endl;
                        MyFile << arachne(Carb) << std::endl;
                      }
                      else{
                        ch++;
                        if (universe->GetInt("truth_targetZ") != 0){
                        //std::cout << "Target Z = " << universe->GetInt("truth_targetZ")<< std::endl;
                        std::cout << "Target ID = " << universe->GetInt("truth_targetID")<< std::endl;
                        }
                      }
                    } //end of bkg in T3
                  
                    // bkg from outside target3
                    else{
                      otherTarget++;
                      if(universe->GetInt("truth_targetZ") == 26){
                        otherTargetFe++;
                        if(universe->GetInt("truth_targetID") == 1) {t1fe++;};
                        if(universe->GetInt("truth_targetID") == 2) {t2fe++;};
                        if(universe->GetInt("truth_targetID") == 4) {t4fe++;};
                        if(universe->GetInt("truth_targetID") == 5) {t5fe++;};
                      }
                      else if(universe->GetInt("truth_targetZ") == 82){
                        otherTargetPb++;
                        if(universe->GetInt("truth_targetID") == 1) {t1pb++;};
                        if(universe->GetInt("truth_targetID") == 2) {t2pb++;};
                        if(universe->GetInt("truth_targetID") == 4) {t4pb++;};
                        if(universe->GetInt("truth_targetID") == 5) {
                          t5pb++;
                          struct SliceID PbT5;
                          PbT5.run = universe->GetMCRunN();
                          PbT5.subrun = universe->GetMCSubRunN();
                          PbT5.gate = universe->GetMCGateN();
                          PbT5.slice = universe->GetMCSliceN();

                          MyFile << "Pb T5 CC antinu event # " << t5pb <<std::endl;
                          MyFile << arachne(PbT5) << std::endl;
                        };
                      }
                      else if(universe->GetInt("truth_targetZ") == 12){
                        otherTargetC++;
                      }
                      else{
                        if(517 < universe->GetVertexZTrueMy() &&  universe->GetVertexZTrueMy() < 544) {
                          water++;
                          /*if (universe->GetInt("mc_targetNucleus") == 2212){plastH++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000060120){plastC++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000080160){plastO++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000130270){plastAl++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000140280){plastSi++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000220480){plastTi++;};*/
                        }
                        else{ 
                          otherTargetPlastic++;
                          // UNIQUE values of background I include in the plastic
                          v.push_back(universe->GetInt("mc_targetNucleus"));
                          // count different 
                          if (universe->GetInt("mc_targetNucleus") == 2212){plastH++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000060120){plastC++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000080160){plastO++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000130270){plastAl++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000140280){plastSi++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000220480){plastTi++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000070140){plastN++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000020040){plastHe++;};
                          if (universe->GetInt("mc_targetNucleus") == 1000170350){plastCl++;};

                          //if (universe->GetInt("truth_targetZ") != 0){
                          //std::cout << "Target Z = " << universe->GetInt("truth_targetZ")<< std::endl;
                          //std::cout << "Target ID = " << universe->GetInt("truth_targetID")<< std::endl;
                          //std::cout << "MC nuclei Z = " << universe->GetInt("mc_targetNucleus")<< std::endl;
                          //}
                        }
                      }
                    } // end of bkg from outside t3
                  } // end of material+other target bkg
                } // end of muon kinematics signal

                // muon kinematics bkg
                else{
                  bkg++;
                  notmuonkin++;
                  if (!cutter->PassTrueThetaCut(universe)){notmuonangle++;}
                }
              } // end of CC antineutrino 
              // bkg: not CC/antineutrino
              else {
                bkg++;
                if( 1 != universe->GetInt("mc_current")){
                  neutralcur++;

                  struct SliceID NC;
                  NC.run = universe->GetMCRunN();
                  NC.subrun = universe->GetMCSubRunN();
                  NC.gate = universe->GetMCGateN();
                  NC.slice = universe->GetMCSliceN();

                  //std::cout << "Neutral Current Event # " << neutralcur <<std::endl;
                  //std::cout << arachne(NC) << std::endl;

                  MyFile << "Neutral Current Event # " << neutralcur <<std::endl;
                  MyFile << arachne(NC) << std::endl;
                }

                else if( -14 != universe->GetInt("mc_incoming") ){wrongsign++;}

                else {fidcut++;} // NOT to double count!!
              } 


              for (auto v : variables2d){
                //Reco 2D
	              v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 

                //if( 1 == universe->GetInt("mc_current") &&  -14 == universe->GetInt("mc_incoming") ){
                if(cutter->PassTruth(universe,helicity) ){ // CC, antineutrino, fiducial cuts
                  if (cutter->PassTrueMuEnergyCut(universe) && cutter->PassTrueThetaCut(universe)){
                    // just filling in different materials...
                    if(cutter->IsInTrueMaterial(universe,targetID, targetZ) && universe->GetInt("truth_targetID") == targetID){
                      v->m_selected_mc_reco_Sg.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
                    }
                    // material + other target bkg

                    // TRUE INFORMATION
                    else{
                      // bkg in t3
                      if(universe->GetInt("truth_targetID") == targetID){
                        if(cutter->IsInTrueMaterial(universe,targetID, 26)){
                          v->m_selected_mc_reco_Fe.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight()); 
                        }
                        else if(cutter->IsInTrueMaterial(universe,targetID, 82) ){
                          v->m_selected_mc_reco_Pb.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight()); 
                        }
                        else if(cutter->IsInTrueMaterial(universe,targetID, 6)){
                          v->m_selected_mc_reco_C.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight()); 

                        }
                        else{
                          if (universe->GetInt("truth_targetZ") != 0){
                            std::cout << "Target Z = " << universe->GetInt("truth_targetZ")<< std::endl;
                          }
                        }
                      } //end of bkg in T3
                    
                      // bkg from outside target3
                      else{
                        // all other target/plastic
                        v->m_selected_mc_reco_other.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight()); 
                          
                        if(universe->GetInt("truth_targetZ") == 26){
                          v->m_selected_mc_reco_otherFe.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight()); 
                        }
                        else if(universe->GetInt("truth_targetZ") == 82){
                          v->m_selected_mc_reco_otherPb.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight()); 
                        }
                        else if(universe->GetInt("truth_targetZ") == 12){
                          v->m_selected_mc_reco_otherC.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight()); 
                        }
                        else{
                          // all the other bkg (plastic)
                          v->m_selected_mc_reco_Plastic.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight()); 
                          if (universe->GetInt("truth_targetZ") != 0){
                            std::cout << "Target Z = " << universe->GetInt("truth_targetZ")<< std::endl;
                          }
                        }
                      } // end of bkg from outside t3
                    } //
                  } // here would be else: events not passing true muon kinematics cuts
                } // here would be else: neutral current events and wrong sign event
              } // end of 2D variable for loop
	
	            for (auto v : variables){
                // All Pass Reco	     
                v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());

                // signal vs background 
                //if( 1 == universe->GetInt("mc_current") &&  -14 == universe->GetInt("mc_incoming") ){
                if(cutter->PassTruth(universe,helicity) ){ // CC, antineutrino, fiducial cuts
                  if (cutter->PassTrueMuEnergyCut(universe) && cutter->PassTrueThetaCut(universe)){
                    if(cutter->IsInTrueMaterial(universe,targetID, targetZ) && universe->GetInt("truth_targetID") == targetID){
                      // signal definition
                      v->m_selected_mc_reco_signal.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                      //signal++;
                    }
                    // material + other target bkg
                    else{
                      // All background reco
                      v->m_selected_mc_reco_bkg.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                      //bkg++;

                      // All material True Value
                      v->m_selected_mc_sb.GetComponentHist("MatBkg")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                      //materialbkg++;
                      // bkg in t3
                      if(universe->GetInt("truth_targetID") == targetID){
                        if(cutter->IsInTrueMaterial(universe,targetID, 26)){
                          v->m_selected_mc_sb.GetComponentHist("FeT3")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                          //fe++ ;
                        }
                        else if(cutter->IsInTrueMaterial(universe,targetID, 82) ){
                          v->m_selected_mc_sb.GetComponentHist("PbT3")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                          //pb++;
                        }
                        else if(cutter->IsInTrueMaterial(universe,targetID, 6)){
                          v->m_selected_mc_sb.GetComponentHist("CT3")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                          //carbon++;
                        }
                        else{
                          //ch++;
                          if (universe->GetInt("truth_targetZ") != 0){
                          std::cout << "Target Z = " << universe->GetInt("truth_targetZ")<< std::endl;
                          }
                        }
                      } //end of bkg in T3
                    
                      // bkg from outside target3
                      else{
                        if(universe->GetInt("truth_targetZ") == 26){
                          v->m_selected_mc_sb.GetComponentHist("otherFe")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                          //otherTargetFe++;
                        }
                        else if(universe->GetInt("truth_targetZ") == 82){
                        v->m_selected_mc_sb.GetComponentHist("otherPb")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                          //otherTargetPb++;
                        }
                        else if(universe->GetInt("truth_targetZ") == 12){
                          v->m_selected_mc_sb.GetComponentHist("otherC")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                          //otherTargetC++;
                        }
                        else{
                          // water target!!
                          if(517 < universe->GetVertexZTrueMy() &&  universe->GetVertexZTrueMy() < 544){
                            v->m_selected_mc_reco_bkg_water.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                          }
                          else{
                            //otherTargetPlastic++;
                            v->m_selected_mc_sb.GetComponentHist("Plastic")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());

                            // different sources of plastic!
                            if (universe->GetInt("mc_targetNucleus") == 2212){ // H
                              v->m_selected_mc_plastic.GetComponentHist("H")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                            };
                            if (universe->GetInt("mc_targetNucleus") == 1000060120){ // C
                              v->m_selected_mc_plastic.GetComponentHist("C")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                            };
                            if (universe->GetInt("mc_targetNucleus") == 1000080160){ //O
                              v->m_selected_mc_plastic.GetComponentHist("O")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                            };
                            if (universe->GetInt("mc_targetNucleus") == 1000130270){ // Al
                              v->m_selected_mc_plastic.GetComponentHist("Al")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                            };
                            if (universe->GetInt("mc_targetNucleus") == 1000140280){ // Si
                              v->m_selected_mc_plastic.GetComponentHist("Si")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                            };
                            if (universe->GetInt("mc_targetNucleus") == 1000220480){ // Ti
                              v->m_selected_mc_plastic.GetComponentHist("Ti")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                            };
                            if (universe->GetInt("mc_targetNucleus") == 1000070140){ // N
                              v->m_selected_mc_plastic.GetComponentHist("N")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                            };
                            if (universe->GetInt("mc_targetNucleus") == 1000020040){ // He
                              v->m_selected_mc_plastic.GetComponentHist("He")->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
                            };
                          }
                        }
                      } // end of bkg from outside t3
                    } // end of material+other target bkg
                  } // end of muon kinematics signal

                  // muon kinematics bkg
                  else{
                    // All background reco
                    v->m_selected_mc_reco_bkg.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                    //bkg++;

                    v->m_selected_mc_reco_bkg_muon.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                    //notmuonkin++;
                    if (!cutter->PassTrueThetaCut(universe)){/*notmuonangle++;*/}
                  }
                } // end of CC antineutrino 

                // bkg: not CC/antineutrino
                else {
                  // All background reco
                  v->m_selected_mc_reco_bkg.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  //bkg++;

                  if( 1 != universe->GetInt("mc_current")){
                    v->m_selected_mc_reco_bkg_NC.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                    //neutralcur++;
                  }
                  else if( -14 != universe->GetInt("mc_incoming")){
                    v->m_selected_mc_reco_bkg_WrongSign.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  }
                  else {
                    v->m_selected_mc_reco_bkg_fidcut.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                  } // NOT to double count!!
                } // end of bkg not CC/antineutrino
              } // end of variable for loop
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
        reco4++;
              
        /* GLOBAL VS LOCAL CUT*/
        if(!cutter->PassMuEnergyCut(dataverse)) continue;
        reco5++;
          
        if(!cutter->PassThetaCut(dataverse)) continue; 
        reco6++;

 	      for (auto v : variables2d){	    
	        v->m_selected_data_reco.hist->Fill(v->GetRecoValueX(*dataverse), v->GetRecoValueY(*dataverse));
	      }
	   
	      for (auto v : variables){
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
  std::cout << "Plane probability cut = " << reco4 <<std::endl;
  std::cout << "Muon Energy cut  = "<< reco5 << std::endl;
  std::cout << "Muon theta cut  = " << reco6 << std::endl;
  std::cout << "**********************************" << std::endl;
  std::cout << "Signal vs background " << std::endl;
  std::cout << "Signal  = "<< signal << std::endl;
  std::cout << "All background = " << bkg << std::endl;
  std::cout << "NC  = "<< neutralcur << std::endl;
  std::cout << "Wrong sign  = "<< wrongsign << std::endl;
  std::cout << "Fiducial cut  = "<< fidcut << std::endl;
  std::cout << "Material background  = "<< materialbkg << std::endl;
  std::cout << "Not true muon kinematics  = "<< notmuonkin << std::endl;
  std::cout << "From that not true angle = " << notmuonangle << std::endl;
  std::cout << "**********************************" << std::endl;
  // Other backround category in terms of materials
  std::cout << "**Material + Other target background** " << std::endl;
  std::cout << "**Background from T3** " << std::endl;
  std::cout << "Fe T3 bkg (should be 0) = " << fe << std::endl;
  std::cout << "Pb T3 = " << pb << std::endl;
  std::cout << "C T3 = " << carbon << std::endl;
  std::cout << "Plastic T3 (should be 0) = " << ch << std::endl;

  std::cout << "**Other target**" << std::endl;
  std::cout << "Other target Total = " << otherTarget << std::endl;
  std::cout << "Other target Fe = " << otherTargetFe << std::endl;
  std::cout << "Other target Pb = " << otherTargetPb << std::endl;
  std::cout << "Other target C = " << otherTargetC << std::endl;
  std::cout << "Other target Plastic = " << otherTargetPlastic << std::endl;
  std::cout << "Water = " << water << std::endl;

  std::cout << "**********************************" << std::endl;
  std::cout << "Other Plastic" << std::endl;
  std::cout << "Plastic H = " << plastH << std::endl;
  std::cout << "Plastic C = " << plastC << std::endl;
  std::cout << "Plastic O = " << plastO << std::endl;
  std::cout << "Plastic Al = " << plastAl << std::endl;
  std::cout << "Plastic Si = " << plastSi << std::endl;
  std::cout << "Plastic Ti = " << plastTi << std::endl;
  std::cout << "Plastic N = " << plastN << std::endl;
  std::cout << "Plastic He = " << plastHe << std::endl;
   std::cout << "Plastic Cl135 = " << plastCl << std::endl;

  std::cout << "**********************************" << std::endl;
  std::cout << "Iron bkg separated into targets" << std::endl;
  std::cout << "Fe T1 = " << t1fe << std::endl;
  std::cout << "Fe T2 = " << t2fe << std::endl;
  std::cout << "Fe T4 = " << t4fe << std::endl;
  std::cout << "Fe T5 = " << t5fe << std::endl;

  std::cout << "**********************************" << std::endl;
  std::cout << "Lead bkg separated into targets" << std::endl;
  std::cout << "Pb T1 = " << t1pb << std::endl;
  std::cout << "Pb T2 = " << t2pb << std::endl;
  std::cout << "Pb T4 = " << t4pb << std::endl;
  std::cout << "Pb T5 = " << t5pb << std::endl;
  
  //return variables;
  MyFile.close();

  sort( v.begin(), v.end() );
  v.erase( unique( v.begin(), v.end() ), v.end() );
                      
  // Displaying the vector after applying std::unique
  for (int i=0; i<v.size();i++){
    cout << v[i] << endl;
  }
}
//=============================================================================

