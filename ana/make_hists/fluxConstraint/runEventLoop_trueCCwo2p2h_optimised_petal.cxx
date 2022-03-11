//==============================================================================
// Loop entries, make cuts, fill histograms.
// * Uses the New Systematics Framework and "Universe" objects.
// * loop universes, make cuts and fill histograms with the correct lateral
// shifts and weights for each universe.
// * TChain --> PlotUtils::ChainWrapper.
// * MnvHXD --> PlotUtils::HistWrapper.
// * Genie, flux, non-resonant pion, and some detector systematics calculated.
//==============================================================================

#include "../../../NUKECCSRC/include/CommonIncludes.h"
#include "../../../NUKECCSRC/include/CVUniverse.h"
#include "VariableAll_Petal.h"  
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../../NUKECCSRC/include/Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include <iostream>
#include <stdlib.h>
#include "../../../NUKECCSRC/include/UtilsNSF.h"
#include "../../../NUKECCSRC/include/Cuts.h"
#include "TParameter.h"

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
    "\t-./runEventLoop Path_to_Output_file Playlist\n\n" <<
    "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
    "\t-Playlist\t \t =\t eg. minervame1A  \n" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    return 0;
  }

  TString dir(argv[1]);
  const string playlist= argv[2];
  int targetID = 99;
  int targetZ = 99;
   
  bool doDIS=false;
  //bool RunCodeWithSystematics = false;
  //bool RunCodeWithSysFluxUniverses = false;

  const std::string plist_string(playlist);
  // Tuples ?
  const std::string mc_file_list(Form("../../include/playlists/DanTuples/MasterAnaDev_MC_%s_MuonKludged.txt", plist_string.c_str()));
  const std::string reco_tree_name("MasterAnaDev");
  
  const bool wants_truth = true; // read in truth!

  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, plist_string, wants_truth);

  util.PrintMacroConfiguration("main");

  // SYSTEMATICS

  //std::map<std::string, std::vector<CVUniverse*> > error_bands =
  //GetErrorBands(util.m_mc);

  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(-14);
  PlotUtils::MinervaUniverse::SetNonResPiReweight(false);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(1);  
  PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);


  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
  NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
  NukeCC_Binning  *binsDef = new NukeCC_Binning();
  
  //PlotUtils::ChainWrapper* chainData = util.m_data;
  //PlotUtils::ChainWrapper* chainMC = util.m_mc;
  PlotUtils::ChainWrapper* chainTruth = util.m_truth; // reading in the Truth tree
  HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);
    
  //double DataPot=  util.m_data_pot; 
  double MCPot=  util.m_mc_pot;  
  //double MCscale=DataPot/MCPot;

  //std::cout << "MC Scale = " << MCscale << std::endl; 
  //std::cout << "Data POT: " << DataPot << std::endl;
  std::cout << "MC POT: " << MCPot << std::endl;

  std::vector<Var*> variablesMC;//,variablesData; 
  std::vector<Var2D*> variables2DMC;//,variables2DData; 


  TString histFileName;
  if(RunCodeWithSystematics){
    histFileName = utils->GetHistFileName( Form("EventSelection_%s_FluxConstraint_optimPetal_sys",plist_string.c_str()), FileType::kAny, targetID, targetZ, helicity ); 
  }

  else{
    histFileName = utils->GetHistFileName( Form("EventSelection_%s_FluxConstraint_nosys",plist_string.c_str()), FileType::kAny, targetID, targetZ, helicity ); 
  } 

  //Works good for the grid submission
  //TString histFileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ );

  TFile fout(dir.Append(histFileName),"RECREATE");	
  
  // MC 
  std::cout << "Processing MC and filling histograms" << std::endl;

  FillVariable(chainTruth, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true,targetID, targetZ, plist_string,doDIS);
    
  for (auto v : variablesMC) {
    v->m_selected_mc_truth_trackerC.SyncCVHistos();
    v->m_selected_mc_truth_waterO.SyncCVHistos();
    v->m_selected_mc_truth_t1fe.SyncCVHistos();
    v->m_selected_mc_truth_t1pb.SyncCVHistos();
    v->m_selected_mc_truth_t2fe.SyncCVHistos();
    v->m_selected_mc_truth_t2pb.SyncCVHistos();
    v->m_selected_mc_truth_t3fe.SyncCVHistos();
    v->m_selected_mc_truth_t3pb.SyncCVHistos();
    v->m_selected_mc_truth_t3c.SyncCVHistos();
    v->m_selected_mc_truth_t4pb.SyncCVHistos();
    v->m_selected_mc_truth_t5fe.SyncCVHistos();
    v->m_selected_mc_truth_t5pb.SyncCVHistos();
    v->m_selected_mc_truth_t25fe.SyncCVHistos();
    v->m_selected_mc_truth_t25pb.SyncCVHistos();
    v->m_selected_mc_truth_t15fe.SyncCVHistos();
    v->m_selected_mc_truth_t15pb.SyncCVHistos();
  }


  for (auto v : variables2DMC) {
    v->m_selected_mc_truth_trackerC.SyncCVHistos();
    v->m_selected_mc_truth_t3c.SyncCVHistos();
    v->m_selected_mc_truth_t25fe.SyncCVHistos();
    v->m_selected_mc_truth_t15fe.SyncCVHistos();
    /*
    v->m_selected_mc_truth_waterO.SyncCVHistos();
    v->m_selected_mc_truth_t1fe.SyncCVHistos();
    v->m_selected_mc_truth_t1pb.SyncCVHistos();
    v->m_selected_mc_truth_t2fe.SyncCVHistos();
    v->m_selected_mc_truth_t2pb.SyncCVHistos();
    v->m_selected_mc_truth_t3fe.SyncCVHistos();
    v->m_selected_mc_truth_t3pb.SyncCVHistos();
    v->m_selected_mc_truth_t3c.SyncCVHistos();
    v->m_selected_mc_truth_t4pb.SyncCVHistos();
    v->m_selected_mc_truth_t5fe.SyncCVHistos();
    v->m_selected_mc_truth_t5pb.SyncCVHistos();
    v->m_selected_mc_truth_t25fe.SyncCVHistos();
    v->m_selected_mc_truth_t25pb.SyncCVHistos();
    v->m_selected_mc_truth_t15fe.SyncCVHistos();
    v->m_selected_mc_truth_t15pb.SyncCVHistos();
    */
  }

  // WRITE HISTOGRAMS TO FILE

  // 1D variables
  for (auto v : variablesMC) {
    v->WriteAllHistogramsToFile(fout, true);
  }

  // 2D Variables
  for (auto v : variables2DMC) {
    v->WriteAllHistogramsToFile(fout,true);
  }
  
  //Writing POT to the HistFile
  fout.cd();
  //auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  //dataPOTOut->Write();
  mcPOTOut->Write(); 
  fout.Close();

  std::cout << "DONE " << plist_string << std::endl;


}//End Main


//=============================================================================
//=============================================================================
// OTHER FUNCTIONS
//=============================================================================
//=============================================================================

std::map<std::string, std::vector<CVUniverse*> > GetErrorBands(PlotUtils::ChainWrapper* chain) {
  typedef std::map<std::string, std::vector<CVUniverse*> > SystMap;
  SystMap error_bands;
  // CV
  error_bands[std::string("CV")].push_back(new CVUniverse(chain));

  if(RunCodeWithSystematics){
    SystMap flux_systematics = PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain,CVUniverse::GetNFluxUniverses());
    error_bands.insert(flux_systematics.begin(), flux_systematics.end());
  }
  return error_bands;
}
// Based on Jeffrey Kleykamp's method
// Do not need other systematics, just the flux
// measuring the number of neutrinos which doesn't change
//  as a function of interaction model since you're dividing by the xsec in the end

// Group universes that change by vertical shift == flux universes
// Meaning: if 500 flux universes, can just apply one cut rather than run the vent 500 times
std::vector<std::vector<CVUniverse*>> groupCompatibleUniverses(const std::map<std::string, std::vector<CVUniverse*>> bands) {
  std::vector<std::vector<CVUniverse*>> groupedUnivs;
  std::vector<CVUniverse*> vertical;
  for(const auto& band: bands){
    if(band.first == "cv") vertical.insert(vertical.begin(), band.second.begin(), band.second.end());
    else{
      for(const auto univ: band.second){
        if(univ->IsVerticalOnly()) vertical.push_back(univ);
        else groupedUnivs.push_back(std::vector<CVUniverse*>{univ});
      }
    }
  }
  groupedUnivs.insert(groupedUnivs.begin(), vertical);
  return groupedUnivs;
}


// Fill Variables
   
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
  
  std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
  for(auto band :error_bands ){std::cout<<"Checking Universe this is universe with name : " << band.first<<std::endl;}
  std::vector<std::vector<CVUniverse*>> error_bands_GROUP = groupCompatibleUniverses(error_bands);
  std::cout<< "error_bands.size() = " << error_bands.size()<<std::endl;
  std::cout<< "Error_Band_GROUPS.size() = " << error_bands_GROUP.size()<<std::endl;
  std::cout<<"Number of Universes set is = "<<    MinervaUniverse::GetNFluxUniverses()<<std::endl;
  
  std::vector<double> Enubin, petalbin;
  std::vector<double> vtxzbin, vtxybin, vtxxbin;
  Enubin = binsDef->GetEnergyBins("EnuFlux"); // 0.1 GeV binning from 0 to 100 GeV
  petalbin = binsDef->GetEnergyBins("petal");
  vtxxbin = binsDef->GetEnergyBins("vtxx"); 
  vtxybin = binsDef->GetEnergyBins("vtxy");
  //vtxzbin = binsDef->GetEnergyBins("vtxz");

  // 1D Variables
  Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
  Var* petal = new Var("Petal", "Petal", petalbin, &CVUniverse::GetDaisyPetalTrue, &CVUniverse::GetDaisyPetalTrue);
  Var* vtxx = new Var("vtxx", "Vertex X", vtxxbin, &CVUniverse::GetVertexXMy, &CVUniverse::GetVertexXTrueMy);
  Var* vtxy = new Var("vtxy", "Vertex Y", vtxybin, &CVUniverse::GetVertexYMy, &CVUniverse::GetVertexYTrueMy);
  //Var* vtxz = new Var("vtxz", "Vertex Z", vtxzbin, &CVUniverse::GetVertexZMy, &CVUniverse::GetVertexZTrueMy);
  //Var* vtxx = new Var("vtxx", "Vertex X", vtxxbin, &CVUniverse::GetVertexXMy, &CVUniverse::GetVertexXTrueMy);
  //Var* vtxy = new Var("vtxy", "Vertex Y", vtxybin, &CVUniverse::GetVertexYMy, &CVUniverse::GetVertexYTrueMy);

  variables = {enu, petal};//, vtxz, vtxx, vtxy};

  Var2D* enu_petal = new Var2D(*enu, *petal);
  Var2D* petal_enu = new Var2D(*petal, *enu);
  Var2D* vtxx_vtxy = new Var2D(*vtxx, *vtxy);
  variables2d = {vtxx_vtxy, enu_petal, petal_enu};

  for (auto v : variables) v->InitializeAllHistograms(error_bands);
  for (auto v : variables2d) v->InitializeAllHistograms(error_bands);


  int all =0 ;
  int truth = 0 ;
  int exclude2p2h = 0;
  int trackerC = 0;
  int waterO = 0;
  int t1fe = 0;
  int t1pb = 0;
  int t2fe = 0;
  int t2pb = 0;
  int t3fe = 0;
  int t3pb = 0;
  int t3c = 0;
  int t4pb = 0;
  int t5fe = 0;
  int t5pb = 0;

  
  //CVUniverse *dataverse = new CVUniverse(chain,0);

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
        for (auto band_GROUP : error_bands_GROUP){

          // Tell the Event which entry in the TChain it's looking at
          //=========================================
          // CUTS in each universe
          //========================================
          band_GROUP.front()->SetEntry(i);
          all++;

          //if(!cutter->PassTruth(universe,helicity)) continue; // CC + antinu + fiducial (apothem)
          if( ! cutter->PassTrueCC(band_GROUP.front(),helicity)) continue;
          if( ! cutter->InHexagonTrue(band_GROUP.front(), 850.) ) continue;
          truth++;

          // exclude if true 2p2h
          if (8 == band_GROUP.front()->GetInt("mc_intType") ) continue;
          exclude2p2h++;

          //if(!cutter->IsInMaterial(universe,targetID,targetZ, /*anyTrakerMod*/false)) continue;
          //if (!(band_GROUP.front()->GetVecElem("mc_vtx",2) >= 5890 && band_GROUP.front()->GetVecElem("mc_vtx",2) <= 8467)) continue;
          //5890.0, 8467.0
          //material++;

          //if(! (band_GROUP.front()->GetInt("mc_targetZ") == 6)) continue;
          //carbon++;

          // TRACKER
          if (cutter->TrackerOnlyTrue4Flux(band_GROUP.front())){
            //5890.0, 8467.0
            if( (band_GROUP.front()->GetInt("mc_targetZ") == 6)){
              trackerC++;
              for (auto universe : band_GROUP){
                for (auto v : variables) v->m_selected_mc_truth_trackerC.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                for (auto v : variables2d) v->m_selected_mc_truth_trackerC.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              }
            }    
          }

          // WATER
          if (cutter->WaterTrue(band_GROUP.front())){
            //5170.0, 5440.0
            if( (band_GROUP.front()->GetInt("mc_targetZ") == 8)){
              waterO++;
              for (auto universe : band_GROUP){
                for (auto v : variables) v->m_selected_mc_truth_waterO.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                //for (auto v : variables2d) v->m_selected_mc_truth_waterO.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              }
            }    
          } 

          // TARGET 1 Fe
          if(cutter->IsInTrueMaterialFlux(band_GROUP.front(), 1, 26, /*anyTrakerMod*/false)){
            t1fe++;
            for (auto universe : band_GROUP){
              for (auto v : variables) {
                v->m_selected_mc_truth_t1fe.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t15fe.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
              }
              for (auto v : variables2d) v->m_selected_mc_truth_t15fe.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
            }
          } 

          // TARGET 1 Pb
          if(cutter->IsInTrueMaterialFlux(band_GROUP.front(), 1, 82, /*anyTrakerMod*/false)){
            t1pb++;
            for (auto universe : band_GROUP){
              for (auto v : variables) {
                v->m_selected_mc_truth_t1pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t15pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux()); 
              }
              //for (auto v : variables2d) v->m_selected_mc_truth_t15pb.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());

            }
          } 

          // TARGET 2 Fe
          if(cutter->IsInTrueMaterialFlux(band_GROUP.front(), 2, 26, /*anyTrakerMod*/false)){
            t2fe++;
            for (auto universe : band_GROUP){
              for (auto v : variables) {
                v->m_selected_mc_truth_t2fe.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t15fe.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t25fe.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
              }
              for (auto v : variables2d){
                v->m_selected_mc_truth_t15fe.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t25fe.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              }        
            }
          } 

          // TARGET 2 Pb
          if(cutter->IsInTrueMaterialFlux(band_GROUP.front(), 2, 82, /*anyTrakerMod*/false)){
            t2pb++;
            for (auto universe : band_GROUP){
              for (auto v : variables) {
                v->m_selected_mc_truth_t2pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t15pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t25pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
              }
              //for (auto v : variables2d){
              //  v->m_selected_mc_truth_t15pb.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              //  v->m_selected_mc_truth_t25pb.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              //} 
            }
          }  

          // TARGET 3 Fe
          if(cutter->IsInTrueMaterialFlux(band_GROUP.front(), 3, 26, /*anyTrakerMod*/false)){
            t3fe++;
            for (auto universe : band_GROUP){
              for (auto v : variables){
                v->m_selected_mc_truth_t3fe.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t15fe.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t25fe.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
              }
              for (auto v : variables2d){
                v->m_selected_mc_truth_t15fe.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t25fe.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              } 
            }
          }

          // TARGET 3 Pb
          if(cutter->IsInTrueMaterialFlux(band_GROUP.front(), 3, 82, /*anyTrakerMod*/false)){
            t3pb++;
            for (auto universe : band_GROUP){
              for (auto v : variables) {
                v->m_selected_mc_truth_t3pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t15pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t25pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
              }
              //for (auto v : variables2d){
              //  v->m_selected_mc_truth_t15pb.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              //  v->m_selected_mc_truth_t25pb.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              //}
            }
          }

          // TARGET 3 C
          if(cutter->IsInTrueMaterialFlux(band_GROUP.front(), 3, 6, /*anyTrakerMod*/false)){
            t3c++;
            for (auto universe : band_GROUP){
              for (auto v : variables) v->m_selected_mc_truth_t3c.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
              for (auto v : variables2d) v->m_selected_mc_truth_t3c.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());

            }            
          }

          // TARGET 4 Pb
          if(cutter->IsInTrueMaterialFlux(band_GROUP.front(), 4, 82, /*anyTrakerMod*/false)){
            t4pb++;
            for (auto universe : band_GROUP){
              for (auto v : variables) {
                v->m_selected_mc_truth_t4pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t15pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t25pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
              }
              //for (auto v : variables2d){
              //  v->m_selected_mc_truth_t15pb.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              //  v->m_selected_mc_truth_t25pb.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              //}
            }
          }

          // TARGET 5 Fe
          if(cutter->IsInTrueMaterialFlux(band_GROUP.front(), 5, 26, /*anyTrakerMod*/false)){
            t5fe++;
            for (auto universe : band_GROUP){
              for (auto v : variables) {
                v->m_selected_mc_truth_t5fe.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t15fe.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t25fe.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
              }
              for (auto v : variables2d){
                v->m_selected_mc_truth_t15fe.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t25fe.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              }
            }
          }

          // TARGET 5 Pb
          if(cutter->IsInTrueMaterialFlux(band_GROUP.front(), 5, 82, /*anyTrakerMod*/false)){
            t5pb++;
            for (auto universe : band_GROUP){
              for (auto v : variables) {
                v->m_selected_mc_truth_t5pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t15pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                v->m_selected_mc_truth_t25pb.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
              }
              //for (auto v : variables2d){
              //  v->m_selected_mc_truth_t15pb.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              //  v->m_selected_mc_truth_t25pb.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
              //}
            }
          }

        }// End Band loop
      }

      else{
        std::cout << "Not MC!!" << std::endl;
      }       
    }//End entries loop



  for(auto band : error_bands){
    std::vector<CVUniverse*> band_universes = band.second;
    for(unsigned int i_universe = 0; i_universe < band_universes.size(); ++i_universe)
    delete band_universes[i_universe];
  } 

  int univ_norm = error_bands_GROUP.size();

  // Printing summary
  std::cout << "**********************************" << std::endl;
  std::cout << "Printing the ";
    isMC? std::cout << "MC ": std::cout << "Data ";
  std::cout << "Summary " << std::endl;
  std::cout << "No cuts = " << all/univ_norm << std::endl;
  std::cout << "True CC antinu in the fiducial volume = " << truth/univ_norm << std::endl;
  std::cout << "After excluding 2p2h = " << exclude2p2h/univ_norm<< std::endl;
  std::cout << "Tracker C = " << trackerC/univ_norm << std::endl;
  std::cout << "Water O = " << waterO/univ_norm << std::endl;
  std::cout << "T1 Fe = " << t1fe/univ_norm << std::endl;
  std::cout << "T1 Pb = " << t1pb/univ_norm << std::endl;
  std::cout << "T2 Fe = " << t2fe/univ_norm << std::endl;
  std::cout << "T2 Pb = " << t2pb/univ_norm << std::endl;
  std::cout << "T3 Fe = " << t3fe/univ_norm << std::endl;
  std::cout << "T3 Pb = " << t3pb/univ_norm << std::endl;
  std::cout << "T3 C = " << t3c/univ_norm << std::endl;
  std::cout << "T4 Pb = " << t4pb/univ_norm << std::endl;
  std::cout << "T5 Fe = " << t5fe/univ_norm << std::endl;
  std::cout << "T5 Pb = " << t5pb/univ_norm << std::endl;
  std::cout << "All Fe (1-5) = " << (t1fe + t2fe + t3fe +t5fe)/univ_norm  << std::endl;
  std::cout << "All Fe (2-5) = " << (t2fe + t3fe +t5fe)/univ_norm  << std::endl;
  std::cout << "All Pb (1-5) = " << (t1pb + t2pb + t3pb + t4pb +t5pb)/univ_norm  << std::endl;
  std::cout << "All Pb (2-5) = " <<  (t2pb + t3pb + t4pb +t5pb)/univ_norm  << std::endl;
  std::cout << "**********************************" << std::endl;
  
  //return variables;
}
//=============================================================================

