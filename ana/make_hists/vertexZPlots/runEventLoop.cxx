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
#include "Variable.h"  
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

#include "../../include/systematics/Systematics.h"

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
    "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n" <<
    "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
    "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n" <<
    "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    return 0;
  }

  TString dir(argv[1]);
  int targetID = 99;
  int targetZ = 99;

  bool doDIS=false;

  // MasterAnaDev tuples?
  //const std::string mc_file_list("../include/playlists/OfficialMAD_minervaME6B_MCNukeOnly_merged.txt");
  //const std::string data_file_list("../include/playlists/OfficialMAD_minervaME6B_Data_merged.txt");
  //const std::string reco_tree_name("MasterAnaDev");

  // NukeCC Tuples ?
  const std::string mc_file_list("../../include/playlists/minervame6B_mc_ML.txt"); 
  const std::string data_file_list("../../include/playlists/minervame5A_data_TBV.txt");
  const std::string reco_tree_name("MasterAnaDev");
  
  const std::string plist_string("minervame6A");
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
    histFileName = utils->GetHistFileName( "EventSelection_ML_ME6A_sys", FileType::kAny, targetID, targetZ, helicity ); 
  }

  else{
    histFileName = utils->GetHistFileName( "EventSelection_ML_Target3_ME6A_nosys", FileType::kAny, targetID, targetZ, helicity ); 
  };

  //TString histFileName = utils->GetHistFileName( "EventSelection_ML_ME6A", FileType::kAny, targetID, targetZ, helicity ); 

  //Works good for the grid submission
  //TString histFileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ );

  TFile fout(dir.Append(histFileName),"RECREATE");	
   
  // MC 
  std::cout << "Processing MC and filling histograms" << std::endl;

  /*FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true,targetID, targetZ, plist_string,doDIS);     
  for (auto v : variablesMC) {
    v->m_selected_mc_reco_water.SyncCVHistos();
    v->m_selected_mc_reco_carbon.SyncCVHistos();
    v->m_selected_mc_reco_lead.SyncCVHistos();  
    v->m_selected_mc_reco_iron.SyncCVHistos();
    v->m_selected_mc_reco_scintillator.SyncCVHistos();

  }*/
  //for (auto v : variables2DMC) v->m_selected_mc_reco.SyncCVHistos();
   
  // DATA
  std::cout << "Processing Data and filling histograms" << std::endl;

  FillVariable(chainData, helicity, utils, cutter,binsDef,variablesData,variables2DData,false,targetID, targetZ, plist_string,doDIS);
  for (auto v : variablesData) {
    v->m_selected_data_reco.SyncCVHistos();
    v->m_selected_data_reco_water.SyncCVHistos();
    v->m_selected_data_reco_carbon.SyncCVHistos();
    v->m_selected_data_reco_lead.SyncCVHistos();  
    v->m_selected_data_reco_iron.SyncCVHistos();
    v->m_selected_data_reco_scintillator.SyncCVHistos();
  }
  //for (auto v : variables2DData) v->m_selected_data_reco.SyncCVHistos();

  // WRITE HISTOGRAMS TO FILE

  // 1D variables
  //for (auto v : variablesMC) {
 //   v->WriteAllHistogramsToFile(fout, true);
  //}

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
  //for (auto v : variables2DMC) {
  //  v->WriteAllHistogramsToFile(fout,true);
  //}

  //for (auto v : variables2DData) {
  //  v->WriteAllHistogramsToFile(fout,false);
  //}
 
  // 2D Plotting
  //for(int i=0; i< variables2DMC.size();i++){
    //Plot2D(variables2DMC[i]->m_selected_mc_reco.hist, variables2DMC[i]->GetName(), variables2DMC[i]->GetNameX(), variables2DMC[i]->GetNameY());
    //Plot2D(variables2DData[i]->m_selected_data_reco.hist, variables2DData[i]->GetName(), variables2DData[i]->GetNameX(),variables2DData[i]->GetNameY());   
  //}//End 2D plotting
	
  //Writing POT to the HistFile
  fout.cd();
  auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  //auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  dataPOTOut->Write();
  //mcPOTOut->Write(); 

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
  
  std::vector<double> vtxzbin;

  if (doDIS){
  }
  else{
    vtxzbin = binsDef->GetEnergyBins("vtxz_all"); 
  
  }
  //Q2bin = binsDef->GetSidebandBins("Q2");
  //Wbin = binsDef->GetSidebandBins("W");

  // 1D Variables
  Var* vtxz = new Var("vtxz", "Vertex Z (cm)", vtxzbin, &CVUniverse::GetVertexZMy, &CVUniverse::GetVertexZTrueMy);

  variables = {vtxz}; //{enu,ehad};   

  // 2D Variables 
  //Var2D* pTmu_pZmu = new Var2D(*pTmu, *pZmu);
  
  //variables2d = {pTmu_pZmu };//{enu_ehad, Q2_W};
   
  //for (auto v : variables2d) v->InitializeAllHistograms(error_bands);
  for (auto v : variables) v->InitializeAllHistograms(error_bands);

  int reco0=0;
  int reco1=0;
  int reco2=0;
  int reco3=0;
  int reco4=0;
  int scintillator=0;
  int water=0;
  int iron=0;
  int lead=0;
  int carbon=0; 
  
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
            if( ! cutter->PassTrueCC(universe,helicity)) continue;
            if( ! cutter->InHexagonTrue(universe, 850.) ) continue;
            //if(!cutter->PassReco(universe,helicity)) continue;   
            reco1++;

            //if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
            //if( universe->GetVecElem("ANN_plane_probs",0) < 0.2 ) continue;	   
            reco2++;

            if(!cutter->PassMuEnergyCut(universe)) continue;
            reco3++;

            if(!cutter->PassThetaCut(universe))continue;
            reco4++;

            // WATER
            if(cutter->WaterTrue(universe)){
            //if(universe->GetInt("MasterAnaDev_targetID") ==2){
              water++;
              for (auto v : variables){
              v->m_selected_mc_reco_water.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              }
            }
            
            // CARBON
            else if(cutter->IsInTrueMaterial(universe,3,6, /*anyTrakerMod*/false)){
              carbon++;
              for (auto v : variables){
                v->m_selected_mc_reco_carbon.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              }
            }
          
            //IRON T1235
            else if( (cutter->IsInTrueMaterial(universe,1,26, /*anyTrakerMod*/false)) || (cutter->IsInTrueMaterial(universe,2,26, /*anyTrakerMod*/false))
                || (cutter->IsInTrueMaterial(universe,3,26, /*anyTrakerMod*/false)) ||  (cutter->IsInTrueMaterial(universe,5,26, /*anyTrakerMod*/false))){
            //else if( (cutter->IsInTrueMaterial(universe,3,26, /*anyTrakerMod*/false))){
              iron++;
              for (auto v : variables){
                v->m_selected_mc_reco_iron.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              }
            }
            
            // LEAD T12345
            else if( (cutter->IsInTrueMaterial(universe,1,82, /*anyTrakerMod*/false)) || (cutter->IsInTrueMaterial(universe,2,82, /*anyTrakerMod*/false)) 
                || (cutter->IsInTrueMaterial(universe,3,82, /*anyTrakerMod*/false)) || (cutter->IsInTrueMaterial(universe,4,82, /*anyTrakerMod*/false)) 
                ||  (cutter->IsInTrueMaterial(universe,5,82, /*anyTrakerMod*/false)) ){
            //else if((cutter->IsInTrueMaterial(universe,3,82, /*anyTrakerMod*/false))  ){
              lead++;
              for (auto v : variables){
                v->m_selected_mc_reco_lead.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());   
              }
            }

            else if(cutter->TrackerOnly(dataverse)){
              for (auto v : variables){
                scintillator++;
                v->m_selected_mc_reco_scintillator.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              }
            }

          } // End band's universe loops
        }// End Band loop
      }

      else{

      dataverse->SetEntry(i);
      reco0++;
  
      if(!cutter->PassReco(dataverse,helicity)) continue;
      reco1++;
      
      //if( dataverse->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
      //if( dataverse->GetVecElem("ANN_plane_probs",0) < 0.2 ) continue;	    
      reco2++;

      if(!cutter->PassMuEnergyCut(dataverse)) continue;   
      reco3++;

      if(!cutter->PassThetaCut(dataverse))continue;
      reco4++;
    

      // WATER
      if(cutter->Water(dataverse)){
      //if(universe->GetInt("MasterAnaDev_targetID") ==2){
        water++;
        for (auto v : variables){
        v->m_selected_data_reco_water.hist->Fill(v->GetRecoValue(*dataverse, 0));	
        }
      }
      
      // CARBON
      else if( (cutter->IsInMaterial(dataverse,3,6, /*anyTrakerMod*/false) && dataverse->GetInt("MasterAnaDev_targetID") == 3)){
        carbon++;
        for (auto v : variables){
          v->m_selected_data_reco_carbon.hist->Fill(v->GetRecoValue(*dataverse, 0));	

        }
      }
    
      //IRON Ts235
      else if( (cutter->IsInMaterial(dataverse,2,26, /*anyTrakerMod*/false)  && dataverse->GetInt("MasterAnaDev_targetID") == 2)
          || (cutter->IsInMaterial(dataverse,3,26, /*anyTrakerMod*/false) && dataverse->GetInt("MasterAnaDev_targetID") == 3) ||  (cutter->IsInMaterial(dataverse,5,26, /*anyTrakerMod*/false) && dataverse->GetInt("MasterAnaDev_targetID") == 5)){
      //else if( (cutter->IsInTrueMaterial(universe,3,26, /*anyTrakerMod*/false))){
        iron++;
        for (auto v : variables){
          v->m_selected_data_reco_iron.hist->Fill(v->GetRecoValue(*dataverse, 0));	

        }
      }
      
      // LEAD T2345
      else if( (cutter->IsInMaterial(dataverse,2,82, /*anyTrakerMod*/false) && dataverse->GetInt("MasterAnaDev_targetID") == 2) 
          || (cutter->IsInMaterial(dataverse,3,82, /*anyTrakerMod*/false)  && dataverse->GetInt("MasterAnaDev_targetID") == 3) || (cutter->IsInMaterial(dataverse,4,82, /*anyTrakerMod*/false) && dataverse->GetInt("MasterAnaDev_targetID") == 4) 
          ||  (cutter->IsInMaterial(dataverse,5,82, /*anyTrakerMod*/false)  && dataverse->GetInt("MasterAnaDev_targetID") == 5) ){
      //else if((cutter->IsInTrueMaterial(universe,3,82, /*anyTrakerMod*/false))  ){
        lead++;
        for (auto v : variables){
          v->m_selected_data_reco_lead.hist->Fill(v->GetRecoValue(*dataverse, 0));	

        }
      }

      else if(cutter->TrackerOnly(dataverse)){
        for (auto v : variables){
          scintillator++;
          v->m_selected_data_reco_scintillator.hist->Fill(v->GetRecoValue(*dataverse, 0));	

        }
      }
          
    }

  }//End entries loop


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
  std::cout << "Plane prob. cut = " << reco2 << std::endl;
  std::cout << "Muon Energy cut  = "<< reco3 << std::endl;
  std::cout << "Muon theta cut  = " << reco4 << std::endl;
  std::cout<<" Water = " << water << std::endl;
  std::cout<<" Carbon = " << carbon << std::endl;
  std::cout<<" Iron = "<< iron <<std::endl;
  std::cout<<" Lead = " << lead << std::endl; 
  std::cout<< "Scintillator = "<<scintillator<<std::endl;

  std::cout << "**********************************" << std::endl;
  
  //return variables;
}
//=============================================================================

