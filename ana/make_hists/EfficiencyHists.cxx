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
#include "../include/Variable.h" 
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


void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID=1, int targetZ=26, const string playlist="minervame1A", bool doDIS=true);

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
   int targetID = atoi(argv[2]);
   int targetZ = atoi(argv[3]);

	
  // TString dir(argv[1]);
  // int targetID = 1;
  // int targetZ = 26;
  // const string playlist= argv[4];
   
  const std::string mc_file_list("../include/playlists/NukeCC_MC_minervame6A_MuonKludged.txt");
  const std::string data_file_list("../include/playlists/NukeCC_Data_minervame6A_MuonKludged.txt");
  const std::string reco_tree_name("NukeCC");
  
  bool doDIS=false;
  const std::string plist_string("minervame6A");
  const bool wants_truth = true;
  //const bool is_grid = false;
  // is grid removed after update of MAT 07/12/2021

  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, plist_string, wants_truth);
   //PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth, is_grid);

    util.PrintMacroConfiguration("main");   

  //=========================================
  // Systematics
  //=========================================
 // std::map<std::string, std::vector<CVUniverse*> > error_bands =
   //   GetErrorBands(util.m_mc);

   PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
   PlotUtils::MinervaUniverse::SetNuEConstraint(true);
   PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
   PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
   PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
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
  
    PlotUtils::ChainWrapper* chainTruth = util.m_truth;
    PlotUtils::ChainWrapper* chainMC = util.m_mc;
    HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);
    cout<<"Helicity"<<helicity<<endl;

    double DataPot=  util.m_data_pot; 
    double MCPot=  util.m_mc_pot; 
 //  double total_pot_data,total_pot_mc;
  // utils->getPOT(total_pot_data,total_pot_mc);  
   double  MCscale=DataPot/MCPot;
  // double  MCscale=1.0;
 
    std::cout<<" MCScale= "<<MCscale<<std::endl; 
   std::vector<Var*> variablesMC,variablesTruth; 
   std::vector<Var2D*> variables2DMC,variables2DTruth;

  TString histFileName;
  if(RunCodeWithSystematics){
    histFileName = utils->GetHistFileName( "Efficiency_ML_ME6A_sys", FileType::kAny, targetID, targetZ, helicity ); 
  }

  else{
    histFileName = utils->GetHistFileName( "Efficiency_ML_ME6A_nosys", FileType::kAny, targetID, targetZ, helicity ); 
  } 
   
     
   TFile fout(dir.Append(histFileName),"RECREATE");	
   
   // For 1D variables 
   FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true, targetID, targetZ, plist_string,doDIS);
       
   for (auto v : variablesMC) v->m_selected_mc_reco.SyncCVHistos();
   for (auto v : variables2DMC) v->m_selected_mc_reco.SyncCVHistos();
   
   FillVariable(chainTruth, helicity, utils, cutter,binsDef,variablesTruth,variables2DTruth,false, targetID, targetZ, plist_string,doDIS);
   
   for (auto v : variablesTruth) v->m_selected_truth_reco.SyncCVHistos();
   for (auto v : variables2DTruth) v->m_selected_truth_reco.SyncCVHistos();
 
   for (auto v : variablesMC) {
     v->WriteAllHistogramsToFileEff(fout, true);
   }


   for (auto v : variablesTruth) {
     v->WriteAllHistogramsToFileEff(fout, false);
   }

   // Plotting If you want for 1D
   /*   for(int i=0; i< variablesMC.size();i++){
     PlotCVAndError(variablesData[i]->data_reco.hist,variablesMC[i]->mc_reco.hist,variablesMC[i]->GetName(),MCscale);
       
     PlotErrorSummary(variablesMC[i]->mc_reco.hist, variablesMC[i]->GetName());
     PlotStacked(variablesData[i]->data_reco_sb.hist,variablesMC[i]->mc_sb.GetHistArray(),MCscale, variablesMC[i]->mc_sb.GetName(), variablesMC[i]->mc_sb.GetName());
   }//End 1D plotting 
   */
   for (auto v : variables2DMC) {
     v->WriteAllHistogramsToFileEff(fout,true);
   }

   for (auto v : variables2DTruth) {
     v->WriteAllHistogramsToFileEff(fout,false);
   }

  fout.cd();
  auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  dataPOTOut->Write();
  mcPOTOut->Write();   


 
   //Plotting in 2D
   
   //for(int i=0; i< variables2DMC.size();i++){
     //Plot2D(variables2DMC[i]->m_selected_mc_reco.hist, variables2DMC[i]->GetName(), variables2DMC[i]->GetNameX(), variables2DMC[i]->GetNameY()); //Plotting line that I somehow cannot delete without producing memory errors, but no one else can reproduce. --ANF 2020.4.6
     //Plot2D(variables2DData[i]->data_reco.hist, variables2DData[i]->GetName(), variables2DData[i]->GetNameX(),variables2DData[i]->GetNameY());
     
   //}//End 2D plotting
  std::cout << "DONE" << std::endl;

}//End Main

   
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID, int targetZ, const string playlist, bool doDIS){
 // std::map< std::string, std::vector<CVUniverse*> > error_bands = utils->GetErrorBands(chain);
    
  std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
  
  std::vector<double> ThetaMuBin,Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin,xbinBrian;
  std::vector<double> x09bin, xfinebin;
  //std::vector<double> ANNPlaneProbBin;
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
    //ANNPlaneProbBin = binsDef->GetEnergyBins("ANNPlaneProb");
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
  Var* x09 = new Var("x09", "x09", x09bin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
  Var* xfine = new Var("xfine", "xfine", xfinebin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
  Var* xBrian = new Var("xBrian", "xBrian", xbinBrian, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
  Var* y = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);

  Var* pTmu = new Var("pTmu", "pTmu", pTbin, &CVUniverse::GetMuonPt, &CVUniverse::GetlepPtTrue);
  Var* pZmu = new Var("pZmu", "pZmu", pZbin, &CVUniverse::GetMuonPz, &CVUniverse::GetlepPzTrue);
  Var* vtxz = new Var("vtxz", "Vertex Z", vtxzbin, &CVUniverse::GetVertexZMy, &CVUniverse::GetVertexZTrueMy);
  //Var *ANNPlaneProb = new Var("ANNPlaneProb", "ANNPlaneProb", ANNPlaneProbBin, &CVUniverse::GetANNPlaneProb, &CVUniverse::GetANNPlaneProb);
  Var* planeDNN = new Var("planeDNN", "planeDNN", planeDNNbin, &CVUniverse::GetplaneDNNReco, &CVUniverse::GetplaneDNNTrue);

  variables = {emu, enu, x, x09, xfine, xBrian, vtxz, planeDNN, pTmu, pZmu}; //{enu,ehad}; 

  // 2D Variables 
  Var2D* pTmu_pZmu = new Var2D(*pTmu, *pZmu);
  Var2D* W_Q2 = new Var2D(*W, *Q2);
  Var2D* enu_ehad = new Var2D(*enu, *ehad);
  Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
  Var2D* x_y = new Var2D(*x, *y);  // y var
  Var2D* x_Q2 = new Var2D(*x, *Q2);  // y var
  
  variables2d = {pTmu_pZmu };
   
  //smakefor (auto v : variables2d) v->InitializeAllHistograms(error_bands);
  for (auto v : variables2d) v->InitializeAllHistograms(error_bands);
  for (auto v : variables) v->InitializeAllHistograms(error_bands);

   int mc0=0;
   int mc1=0;
   int mc2=0;
   int mc3=0;
   int mc4=0; 
   int mc5=0; 
   int mc6=0; 
   int mc7=0; 
   int mc8=0; 
   int mc9=0; 
   int mc10=0; 
   int mc11=0; 
   int mc12=0; 
   int mc13=0; 
   int mc14=0; 
   int mc15=0; 
   int mc16=0; 
   int mc17=0; 
     int mc_truth0=0.0;
     int mc_truth1=0.0;
     int mc_truth2=0.0; 
     int mc_truth3=0.0; 
     int mc_truth4=0.0; 
     int mc_truth5=0.0; 
     int mc_truth6=0.0; 
     int mc_truth7=0.0; 
   int allcuts=0;
   
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
       std::vector<CVUniverse*> error_band_universes = band.second;
       for (auto universe : error_band_universes){
	 // Tell the Event which entry in the TChain it's looking at
	 universe->SetEntry(i);
	 //=========================================
	 // CUTS in each universe
	 //=========================================
	 
	if(isMC){ // NUMERATOR
/////////////////////////////////////////////////////////////////
	
	   universe->SetEntry(i);
     mc0++;

    if(!cutter->PassReco(universe,helicity)) continue;  
    mc1++;
              
     if(!cutter->IsInMaterial(universe,targetID,targetZ, /*anyTrakerMod*/false)) continue;
     //if(!cutter->IsInMaterial(universe,targetIDs[t],targetZ, /*anyTrakerMod*/false)) continue;
     mc2++;
              
    if(targetID<10 && universe->GetInt("NukeCC_targetID") != targetID) continue;
    //if(targetIDs[t]<10 && universe->GetInt("NukeCC_targetID") != targetIDs[t]) continue;
      mc3++;
              
      if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
      //if( universe->GetVecElem("ANN_plane_probs",0) < 0.2 ) continue;	   
      mc4++;

      if (!cutter->PassTrueCC(universe, helicity)) continue; //true CC, true antinu
      // NO xy separation,  NO APOTHEM CUT
      mc5++;

	    if(!cutter->IsInTrueMaterial(universe,targetID, targetZ,false)) continue; // true target + material
      mc6++;
    
	   for (auto v : variables2d){
	     if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(universe)) continue;
	     if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(universe)) continue;

       if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassTrueMuEnergyCut(universe)) continue;
       // NO TRUE angle cut, efficiency corrected

        v->m_selected_mc_reco.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight()); 
      }

	   for (auto v : variables){
	     if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(universe)) continue;
	     if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(universe))continue;
	     if( v->GetName()=="Enu") mc7++;

       if( v->GetName()!="Emu")   if(!cutter->PassTrueMuEnergyCut(universe)) continue;
	     // NO TRUE angle cut, efficiency corrected
	     if( v->GetName()=="Enu") mc8++;
	       v->m_selected_mc_reco.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetWeight());
	     //v->m_selected_mc_sb.GetComponentHist("MC")->Fill(v->GetTrueValue(*universe, 0), universe->GetWeight());
	   }
	 } else if(!isMC ){ // DENOMINATOR
	   if (!cutter->PassTruth(universe, helicity)) continue; // True fiducial z distance+apothem, true CC, true antinu
      mc_truth0++;

	    if(!cutter->IsInTrueMaterial(universe,targetID, targetZ,false)) continue; // true target + material
      mc_truth1++;

	   for (auto v : variables2d){
	    if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassTrueThetaCut(universe)) continue;	     
	    if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassTrueMuEnergyCut(universe)) continue;
	    v->m_selected_truth_reco.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
	   }
	   for (auto v : variables){
	     if( v->GetName()!="ThetaMu") if(!cutter->PassTrueThetaCut(universe))continue;
       if( v->GetName()=="Enu") mc_truth2++;
	     if( v->GetName()!="Emu")   if(!cutter->PassTrueMuEnergyCut(universe)) continue;
       if( v->GetName()=="Enu") mc_truth3++;
      v->m_selected_truth_reco.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
	   }
	 }
       } // End band's universe loop
     }// End Band loop
   }//End entries loop

   std::cout<<"**********************************"<<std::endl;
      std::cout<<" Summary "<<std::endl;
   std::cout<<"**********************************"<<std::endl;
     std::cout<<"Printing the Numerator Summary "<<std::endl;
     std::cout<<" Reco cut0 = "<<mc0<<std::endl;
     std::cout<<" Reco cut1 = "<<mc1<<std::endl;
     std::cout<<" Reco cut2 = "<<mc2<<std::endl;
     std::cout<<" Reco cut3 = "<<mc3<<std::endl;
     std::cout<<" Reco cut4 = "<<mc4<<std::endl;
     std::cout<<" Reco cut5 = "<<mc5<<std::endl;
     std::cout<<" Reco cut6 = "<<mc6<<std::endl;
     std::cout<<" Reco cut7 = "<<mc7<<std::endl;
     std::cout<<" Reco cut8 = "<<mc8<<std::endl;
   std::cout<<"**********************************"<<std::endl;
     std::cout<<"Printing the Truth Summary "<<std::endl;
     std::cout<<" Truth cut0 = "<<mc_truth0<<std::endl;
     std::cout<<" Truth cut1 = "<<mc_truth1<<std::endl;
     std::cout<<" Truth cut2 = "<<mc_truth2<<std::endl;
     std::cout<<" Truth cut3 = "<<mc_truth3<<std::endl;
   //	return variables;
}
//============================================================================================================================
// Main
