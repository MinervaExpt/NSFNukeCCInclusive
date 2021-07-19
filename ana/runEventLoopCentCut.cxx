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
#include "include/CVUniverse.h"
#include "../include/Variable.h"
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
#include <stdlib.h>
#include "Cintex/Cintex.h"
#include "include/NukeCCUtilsNSF.h"
#include "include/NukeCC_Cuts.h"
#include "include/NukeCC_CentCuts.h"
#include "PlotUtils/Cutter.h"
//#include "../includes/GetVariables.h"
#include <iomanip>
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
//#include "PlotUtils/MnvTuneSystematics.h"



// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we
// shield it.
#ifndef __CINT__
#include "../include/plotting_functions.h"
//#ifndef __CINT__//for PlotUtils::cuts_t<>
#endif
#include "PlotUtils/MacroUtil.h" 
using namespace reco;
using namespace NUKECC_ANA;

typedef VarLoop::Variable Var;
typedef Var2DLoop::Variable2D Var2D;


void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils ,const cut_map_t passes_cuts ,NukeCC_Binning  *binsDef,PlotUtils::Cutter<CVUniverse>& selection,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID=1, int targetZ=26, const string playlist="minervame1A");

//============================================================================================================================
// GetErrorBands 
//============================================================================================================================

  std::map<std::string, std::vector<CVUniverse*> > GetErrorBands(
  PlotUtils::ChainWrapper* chain) {
 
  typedef std::map<std::string, std::vector<CVUniverse*> > SystMap;

  SystMap error_bands;

  // CV
  error_bands[std::string("CV")].push_back(new CVUniverse(chain));

  if(RunCodeWithSystematics) {
  	int n_flux_universes =100 ;
  	SystMap flux_systematics = 
      	PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain,n_flux_universes);
  	error_bands.insert(flux_systematics.begin(), flux_systematics.end());

	//GENIE
	SystMap genie_systematics = 
	PlotUtils::GetGenieSystematicsMap<CVUniverse>(chain);// change that true to a switch on do_nonrespi_tune
	error_bands.insert(genie_systematics.begin(), genie_systematics.end());
	    
	//Muon Angle systematics...
	SystMap angle_systematics = 
	PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);//NSFDefaults::beamThetaX_Err and beamThetaY_Err);
	error_bands.insert(angle_systematics.begin(), angle_systematics.end());
	  
	//Muon P Systematics...
	//0.01 ---->Why>
	SystMap muonminosP_systematics = PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
	error_bands.insert(muonminosP_systematics.begin(), muonminosP_systematics.end());
	  
	SystMap muonminervaP_systematics = PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
	error_bands.insert(muonminervaP_systematics.begin(), muonminervaP_systematics.end()); 

	//MuonMomentum Resolution Systematics....
	SystMap muonP_resolutions = PlotUtils::GetMuonResolutionSystematicsMap<CVUniverse>(chain);
	error_bands.insert(muonP_resolutions.begin(),muonP_resolutions.end());
	  
	//Minos Efficiency Systematics....
	SystMap minos_efficiency = PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
	error_bands.insert(minos_efficiency.begin(),minos_efficiency.end());

	//2p2h 
	SystMap _2p2h_systematics = PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
	error_bands.insert(_2p2h_systematics.begin(),_2p2h_systematics.end());

	//RPA
	SystMap RPA_systematics = PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
	error_bands.insert(RPA_systematics.begin(),RPA_systematics.end());  
  }

  return error_bands;
}

//============================================================================================================================
// FillVariable 
//============================================================================================================================
    
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , const cut_map_t passes_cuts ,NukeCC_Binning  *binsDef ,PlotUtils::cuts_t<CVUniverse>& sidebands,PlotUtils::Cutter<CVUniverse>& selection, std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist){
   // std::map< std::string, std::vector<CVUniverse*> > error_bands = utils->GetErrorBands(chain)
   
    reco::PassMuEnergyCut<CVUniverse, PlotUtils::detail::empty> energyCut;
    reco::PassThetaCut<CVUniverse, PlotUtils::detail::empty> thetaCut;
    reco::IsDIS<CVUniverse, PlotUtils::detail::empty> DIScut;
     //reco::PassMuEnergyCut energyCut; 
    std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
   
    std::map<std::string, std::vector<CVUniverse*> >::iterator itr;
   
    std::map<const std::string, int> error_name;
    for(itr = error_bands.begin(); itr != error_bands.end(); ++itr) error_name.insert(pair<std::string, const int>((itr->second)[0]->ShortName(), (itr->second).size())); 

    std::map<std::string, const int>::iterator itr_m;
   
    std::vector<double> Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin;

   //if (doDIS){
     Enubin  = binsDef->GetDISBins("Enu"); 
     Emubin  = binsDef->GetDISBins("Emu"); 
     Ehadbin = binsDef->GetDISBins("Ehad");
     xbin    = binsDef->GetDISBins("x");
     ybin    = binsDef->GetDISBins("y");
  //   }
  // else{
  //   Enubin  = binsDef->GetEnergyBins("Enu"); 
  //   Emubin  = binsDef->GetEnergyBins("Emu"); 
  //   Ehadbin = binsDef->GetEnergyBins("Ehad");
  //   xbin    = binsDef->GetEnergyBins("x");
  //   ybin    = binsDef->GetEnergyBins("y");
  // }
  
   Q2bin = binsDef->GetSidebandBins("Q2");
   Wbin  = binsDef->GetSidebandBins("W");

   //For 1D varriable

   Var* enu  = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
   Var* ehad = new Var("Ehad", "Ehad (GeV)", Ehadbin, &CVUniverse::GetEhadGeV, &CVUniverse::GetEhadTrueGeV);
   Var* Q2   = new Var("Q2", "Q2 (GeV^2)", Q2bin, &CVUniverse::GetQ2Reco, &CVUniverse::GetQ2True);
   Var* W    = new Var("W", "W (GeV)", Wbin, &CVUniverse::GetWReco, &CVUniverse::GetWTrue);
   Var* emu  = new Var("Emu", "Emu (GeV)", Emubin, &CVUniverse::GetMuonEGeV, &CVUniverse::GetMuonETrueGeV);
   Var* x    = new Var("x", "x", xbin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
   Var* y    = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);
   
    //reco::IsDIS discut()
   //std::vector<Var*> variables = {enu,ehad}; 
   variables = {emu,enu};//{enu,ehad}; 
   

     for (auto v : variables) {
         	std::cout << "Started with this mask: " << v->ignoreTheseCuts << "\n";
           	v->MatchCutsToVars(sidebands); //N.B.: Only sidebands can be ignored.
    //Double-check MatchCutsToVars
    std::cout << "Ignored these Cuts for Variable " << v->GetName() << " using mask " << v->ignoreTheseCuts << ":\n";
    for(size_t whichCut = 0; whichCut < v->ignoreTheseCuts.size(); ++whichCut)	
    {
      if(v->ignoreTheseCuts[whichCut]) std::cout << sidebands[whichCut]->getName() << "\n";	
    }
  }         

     
   //For 2D variable

   Var2D* Q2_W     = new Var2D(*Q2, *W);
   Var2D* enu_ehad = new Var2D(*enu, *ehad);
   Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var

   variables2d = {emu_ehad};//{enu_ehad, Q2_W};
   
   for (auto v : variables2d) v->InitializeAllHistograms(error_bands);
   for (auto v : variables) v->InitializeAllHistograms(error_bands);

   //Migration starts here!
   for (auto v : variables2d) v->SetupResponse(error_name);
   for (auto v : variables)   v->SetupResponse1D(error_name);
   
   int VarCut=0;
   int reco0=0;
   int reco00=0;
   int reco1=0;
   int reco2=0;
   int reco3=0;
   int reco4=0; 
   int allcuts=0;
   int DIS=0;
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
     PlotUtils::detail::empty event;
      
      bool vert_universe_passes_cuts = false;
      bool vert_universe_checked_cuts = false;
      //int weightForCuts = 1;
     reco00++;
     if(isMC){     
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
   //bool passes_cuts = true;
   //cut_map_t passes_cuts = selection.isMCSelected(*universe, event, universe->GetWeight()); //PassesCuts(*universe);
     cut_map_t passes_cuts = selection.isMCSelected(*universe, event, universe->GetWeight()); //PassesCuts(*universe);
             
      for (auto v : variables2d){
                    
      if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")       //   if(!cutter->PassMuEnergyCut(universe)) continue;
      if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu") // if(!cutter->PassThetaCut(universe)) continue;
        v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());

      }

      for (auto v : variables){
                   	       	
      //if(passes_cuts.all()){	
      if(v->IgnoreMyVars(passes_cuts).all()){	
      PlotUtils::detail::empty IDontNeedEvent;     
      //if( v->GetName()!="Emu") if(!energyCut.passesCut(*universe, IDontNeedEvent))continue;
      //if( v->GetName()!="ThetaMu")  if(!thetaCut.passesCut(*universe, IDontNeedEvent))continue;
          VarCut++;
      if(DIScut.passesCut(*universe, IDontNeedEvent)){
      //if (PassesCuts(*universe)){      
          DIS++;  
      v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight()); 
      //v->m_selected_mc_sb.GetComponentHist("MC")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
      //v->m_selected_mc_sb.GetComponentHist("DIS")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight()); 
      }//DIS cut close
   }
 }//Close 1D variable loop
} 
      unv_count++;

   } // End band's universe loop
       
 }//isMC

     else{

       dataverse->SetEntry(i);
             
     //PlotUtils::detail::empty event;
     //if (PassesCuts(*dataverse)){
     //tdead++;
     cut_map_t passes_cuts = selection.isDataSelected(*dataverse, event); //PassesCuts(*universe);
     // reco3++;
     //bool passes_cuts = true;
     //if(!selection.isDataSelected(*dataverse, event).all()) continue;
      for (auto v : variables2d){
      if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu") // if(!cutter->PassMuEnergyCut(dataverse)) continue;
      if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu") // if(!cutter->PassThetaCut(dataverse)) continue;	     
        v->m_selected_data_reco.univHist(dataverse)->Fill(v->GetRecoValueX(*dataverse), v->GetRecoValueY(*dataverse));
	v->m_selected_data_reco.hist->Fill(v->GetRecoValueX(*dataverse), v->GetRecoValueY(*dataverse));
      }

      for (auto v : variables){
            
      if(v->IgnoreMyVars(passes_cuts).all()){
      //if(passes_cuts.all()){
                 
      PlotUtils::detail::empty IDontNeedEvent;     
      //if( v->GetName()!="Emu") if(!energyCut.passesCut(*dataverse, IDontNeedEvent))continue;
      //if( v->GetName()!="ThetaMu")  if(!thetaCut.passesCut(*dataverse, IDontNeedEvent))continue;
            VarCut++;
      if(DIScut.passesCut(*dataverse, IDontNeedEvent)){
      DIS++;  
      //if (PassesCuts(*dataverse)){      
      v->m_selected_data_reco.hist->Fill(v->GetRecoValue(*dataverse, 0));
      v->m_selected_data_reco_sb.hist->Fill(v->GetRecoValue(*dataverse, 0));
     }      
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
 // std::cout<<"**********************************"<<std::endl;
 std::cout<<"Printing the ";
 isMC? std::cout<<"MC": std::cout<<"Data";

 //std::cout<<" Summary "<<std::endl;

 std::cout<<" No1 cuts = "<<reco00<<std::endl;
 std::cout<<" VarCut = "<<VarCut<<std::endl;
 std::cout<<" DIS = "<<DIS<<std::endl;
/* std::cout<<" Reco Cut = "<<reco1<<std::endl;
 std::cout<<" Material Cut = "<<reco2<<std::endl;
 std::cout<<" Plane/targetID Cuts = "<<reco3<<std::endl;
 std::cout<<" Muon Kinematics Cuts = "<<reco4<<std::endl;
 std::cout<<"**********************************"<<std::endl;
 //return variables;
*/
}

//============================================================================================================================
// Main
//============================================================================================================================

int main(int argc, char *argv[]){
	 ROOT::Cintex::Cintex::Enable();
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
	 int targetID = atoi(argv[2]);
	 int targetZ  = atoi(argv[3]);

		
	   
	 bool doDIS = true;
	 bool isSYS = true;
	 const std::string mc_file_list("../include/playlists/NukeCC_minervame1A_MC_Inextinguishable_merged.txt");
	 const std::string data_file_list("../include/playlists/NukeCC_minervame1A_DATA_Inextinguishable_merged.txt");
	 const std::string plist_string("minervame1A");
	 const std::string reco_tree_name("NukeCC");
	 const bool wants_truth = false;
	 const bool is_grid = false;

	 // PlotUtils::MacroUtil util(reco_tree_name, file_list, plist_string, wants_truth,
	 //                         is_grid);
	 PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth, is_grid);

	 util.PrintMacroConfiguration("main");

	 //=========================================
	 // Systematics
	 //=========================================
	 //std::map<std::string, std::vector<CVUniverse*> > error_bands =
	 //  GetErrorBands(util.m_mc);

         PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
         PlotUtils::MinervaUniverse::SetNuEConstraint(true);
       //  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
         PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
       //  PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
         PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);
         PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);

	 NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
//	 NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
	 NukeCC_Binning  *binsDef = new NukeCC_Binning();
	  const cut_map_t passes_cuts;
	 PlotUtils::ChainWrapper* chainData = util.m_data;
	 PlotUtils::ChainWrapper* chainMC = util.m_mc;
	 HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);
	 double DataPot=  util.m_data_pot; 
	 double MCPot=  util.m_mc_pot; 
	 //double total_pot_data,total_pot_mc;
	 //utils->getPOT(total_pot_data,total_pot_mc);  
	 double  MCscale=DataPot/MCPot;
	 //double  MCscale=1.0;
	 //PlotUtils::cuts_t<CVUniverse> sidebands;

    
	 std::cout<<" MCScale= "<<MCscale<<std::endl; 
	 std::vector<Var*> variablesMC,variablesData; 
	 std::vector<Var2D*> variables2DMC,variables2DData; 

         PlotUtils::cuts_t<CVUniverse> sidebands = reco::getCCInclusiveCuts<CVUniverse>(targetID,targetZ,helicity,false), preCuts;
         PlotUtils::constraints_t<CVUniverse> signal, phaseSpace;






     //  PlotUtils::Cutter<CVUniverse> selectionCriteria(reco::getCCInclusiveCuts<CVUniverse>(targetID,targetZ,helicity,false), std::move(sidebands), std::move(signal), std::move(phaseSpace));
         PlotUtils::Cutter<CVUniverse> selectionCriteria(std::move(preCuts), std::move(sidebands), std::move(signal), std::move(phaseSpace));
         

	 TString histFileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ, helicity ); 
	   
	 TFile fout(dir.Append(histFileName),"RECREATE");	
	   
	 FillVariable(chainMC, helicity, utils,  passes_cuts, binsDef, sidebands,selectionCriteria, variablesMC, variables2DMC, true, targetID, targetZ, plist_string);
            for (auto v : variablesMC)   v-> m_selected_mc_reco.SyncCVHistos();
         	 
         for (auto v : variables2DMC) v-> m_selected_mc_reco.SyncCVHistos();
	 // For 1D variables 
         // std::cout << "\nCut summary for Data:\n" << selectionCriteria << "\n";
           
         std::cout << std::fixed << selectionCriteria << "\n";
         selectionCriteria.resetStats();
	 FillVariable(chainData, helicity, utils, passes_cuts, binsDef,sidebands,selectionCriteria, variablesData, variables2DData, false, targetID, targetZ, plist_string);


	
          std::cout << "\nCut summary for Data:\n"  << selectionCriteria << "\n";
// reco::summarize(std::cout, selectionCriteria, util.GetDataEntries());
       //  for(auto& cut: selectionCriteria) cut->resetStats();


	 for (auto v : variablesData)   v-> m_selected_data_reco.SyncCVHistos();
	 for (auto v : variablesData)   v-> m_selected_data_reco_sb.SyncCVHistos();
	 for (auto v : variables2DData) v-> m_selected_data_reco.SyncCVHistos();
	 
          



	 for (auto v : variablesMC) {
	   v->WriteAllHistogramsToFile(fout, true);
	 }

	 for (auto v : variablesData) {
	   v->WriteAllHistogramsToFile(fout, false);
	 }

	 // Plotting If you want for 1D
	  
	 for(int i=0; i< variablesMC.size();i++){
		 PlotCVAndError(variablesData[i]->m_selected_data_reco.hist,variablesMC[i]->m_selected_mc_reco.hist,variablesMC[i]->GetName(),MCscale);
		       
	 PlotErrorSummary(variablesMC[i]->m_selected_mc_reco.hist, variablesMC[i]->GetName());
	//	 PlotStacked(variablesData[i]->m_selected_data_reco_sb.hist,variablesMC[i]->m_selected_mc_sb.GetHistArray(),MCscale, variablesMC[i]->m_selected_mc_sb.GetName(), variablesMC[i]->m_selected_mc_sb.GetName());
	 }//End 1D plotting 
	 

	 //For 2D variable

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
	 
	 //Plotting in 2D
	   
//	 for(int i=0; i< variables2DMC.size();i++){
//	   Plot2D(variables2DMC[i]->mresp.hist, variables2DMC[i]->GetName(), variables2DMC[i]->GetNameX(), variables2DMC[i]->GetNameY()); //Plotting line that I somehow cannot delete without producing memory errors, but no one else can reproduce. --ANF 2020.4.6
	 //Plot2D(variables2DData[i]->m_selected_data_reco.hist, variables2DData[i]->GetName(), variables2DData[i]->GetNameX(),variables2DData[i]->GetNameY());
	     
//	 }//End 2D plotting

}//End Main

