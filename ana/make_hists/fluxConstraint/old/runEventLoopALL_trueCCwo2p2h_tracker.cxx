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
  const std::string mc_file_list(Form("../../include/playlists/MasterAnaDev_MC_%s_MuonKludged1.txt", plist_string.c_str()));
  const std::string reco_tree_name("MasterAnaDev");
  
  const bool wants_truth = true; // read in truth!

  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, plist_string, wants_truth);

  util.PrintMacroConfiguration("main");

  // SYSTEMATICS

  //std::map<std::string, std::vector<CVUniverse*> > error_bands =
  //GetErrorBands(util.m_mc);

  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
  PlotUtils::MinervaUniverse::SetNonResPiReweight(false);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(10);
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
    histFileName = utils->GetHistFileName( Form("EventSelection_%s_FluxConstraint_sys",plist_string.c_str()), FileType::kAny, targetID, targetZ, helicity ); 
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
  for (auto v : variablesMC) v->m_selected_mc_truth.SyncCVHistos();

  // WRITE HISTOGRAMS TO FILE

  // 1D variables
  for (auto v : variablesMC) {
    v->WriteAllHistogramsToFile(fout, true);
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


// Fill Variables
   
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
  
  std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
  
  std::vector<double> Enubin;
  std::vector<double> vtxzbin, vtxybin, vtxxbin;
  Enubin = binsDef->GetEnergyBins("EnuFlux"); // 0.1 GeV binning from 0 to 100 GeV
  //vtxxbin = binsDef->GetEnergyBins("vtxx"); 
  //vtxybin = binsDef->GetEnergyBins("vtxy");
  //vtxzbin = binsDef->GetEnergyBins("vtxz");

  // 1D Variables
  Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
  //Var* vtxz = new Var("vtxz", "Vertex Z", vtxzbin, &CVUniverse::GetVertexZMy, &CVUniverse::GetVertexZTrueMy);
  //Var* vtxx = new Var("vtxx", "Vertex X", vtxxbin, &CVUniverse::GetVertexXMy, &CVUniverse::GetVertexXTrueMy);
  //Var* vtxy = new Var("vtxy", "Vertex Y", vtxybin, &CVUniverse::GetVertexYMy, &CVUniverse::GetVertexYTrueMy);

  variables = {enu};//, vtxz, vtxx, vtxy};


  for (auto v : variables) v->InitializeAllHistograms(error_bands);


  int all =0 ;
  int truth = 0 ;
  int material = 0;
  int exclude2p2h = 0;
  int carbon = 0;
  
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
        for (auto band : error_bands){
          std::vector<CVUniverse*> error_band_universes = band.second;
          for (auto universe : error_band_universes){
            // Tell the Event which entry in the TChain it's looking at
            //=========================================
            // CUTS in each universe
            //========================================
            universe->SetEntry(i);
            all++;

            //if(!cutter->PassTruth(universe,helicity)) continue; // CC + antinu + fiducial (apothem)
            if( ! cutter->PassTrueCC(universe,helicity)) continue;
            if( ! cutter->InHexagonTrue(universe, 850.) ) continue;
            truth++;

            // exclude if true 2p2h
            if (8 == universe->GetInt("mc_intType") ) continue;
            exclude2p2h++;

            //if(!cutter->IsInMaterial(universe,targetID,targetZ, /*anyTrakerMod*/false)) continue;
            if (!(universe->GetVecElem("mc_vtx",2) >= 5890 && universe->GetVecElem("mc_vtx",2) <= 8467)) continue;
            //5890.0, 8467.0
            material++;

            if(! (universe->GetInt("mc_targetZ") == 6)) continue;
            carbon++;
            

            for (auto v : variables){
              v->m_selected_mc_truth.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
              // True Value and True Weight!!
            }           
          
          } // End band's universe loop
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


  // Printing summary
  std::cout << "**********************************" << std::endl;
  std::cout << "Printing the ";
    isMC? std::cout << "MC ": std::cout << "Data ";
  std::cout << "Summary " << std::endl;
  std::cout << "No cuts = " << all << std::endl;
  std::cout << "True CC antinu in the fiducial volume = " << truth << std::endl;
  std::cout << "After excluding 2p2h = " << exclude2p2h << std::endl;
  std::cout << "Tracker Cut = " << material << std::endl;
  std::cout << "Carbon only = " << carbon << std::endl;
  std::cout << "**********************************" << std::endl;
  
  //return variables;
}
//=============================================================================

