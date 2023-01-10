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
#include "VariableEff.h" 
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

//#include "../../include/systematics/Systematics.h"

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

typedef std::map<std::string, std::vector<CVUniverse*> > SystMap;

std::map<std::string, std::vector<CVUniverse*> > GetErrorBands(PlotUtils::ChainWrapper* chain) {
  
  // return map
  SystMap error_bands;

  // CV 
  error_bands[std::string("CV")].push_back(new CVUniverse(chain));

  return error_bands;
}


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
   //int targetID = atoi(argv[2]);
   //nt targetZ = atoi(argv[3]);
  int targetID = 12345; int targetZ = 0; 
	
  // TString dir(argv[1]);
  // int targetID = 1;
  // int targetZ = 26;
  // const string playlist= argv[4];
   
  const std::string mc_file_list("../../include/playlists/minervame1E_mc_DualVertex_FullDetector.txt");
  const std::string reco_tree_name("MasterAnaDev");
  
  bool doDIS=false;
  const std::string plist_string("minervame1E");
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
  //PlotUtils::MinervaUniverse::SetReadoutVolume("Nuke");
  //Neutron CV reweight is on by default (recommended you keep this on)
  //PlotUtils::MinervaUniverse::SetMHRWeightNeutronCVReweight(true);
  //Elastics are on by default (recommended you keep this on)
  //PlotUtils::MinervaUniverse::SetMHRWeightElastics(true);


     
   NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
   NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
   NukeCC_Binning  *binsDef = new NukeCC_Binning();
  
    PlotUtils::ChainWrapper* chainTruth = util.m_truth;
    PlotUtils::ChainWrapper* chainMC = util.m_mc;
    HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);
    cout<<"Helicity"<<helicity<<endl;

    double MCPot=  util.m_mc_pot; 
 //  double total_pot_data,total_pot_mc;
  // utils->getPOT(total_pot_data,total_pot_mc);  
  // double  MCscale=1.0;
 
  std::cout<<" MC POT= "<<MCPot<<std::endl; 
   std::vector<Var*> variablesMC,variablesTruth; 
   std::vector<Var2D*> variables2DMC,variables2DTruth;

  TString histFileName;
  if(RunCodeWithSystematics){
    histFileName = utils->GetHistFileName( "Efficiency_ME1E_FullDet_sys", FileType::kAny, targetID, targetZ, helicity ); 
  }

  else{
    histFileName = utils->GetHistFileName( "Efficiency_ME1E_FullDet_nosys", FileType::kAny, targetID, targetZ, helicity ); 
  } 
   
     
   TFile fout(dir.Append(histFileName),"RECREATE");	
   
   // For 1D variables 
   std::cout<<"I am here"<<std::endl;
   FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true, targetID, targetZ, plist_string,doDIS);
       
   for (auto v : variablesMC){ 
    v->m_selected_mc_MLreco.SyncCVHistos();
    v->m_selected_mc_reco.SyncCVHistos();
   }
   for (auto v : variables2DMC) v->m_selected_mc_reco.SyncCVHistos();
   
   FillVariable(chainTruth, helicity, utils, cutter,binsDef,variablesTruth,variables2DTruth,false, targetID, targetZ, plist_string,doDIS);
   
   for (auto v : variablesTruth) v->m_selected_truth_reco.SyncCVHistos();
   for (auto v : variables2DTruth) v->m_selected_truth_reco.SyncCVHistos();
 
   for (auto v : variablesMC) {
     v->WriteAllHistogramsToFile(fout, true);
   }


   for (auto v : variablesTruth) {
     v->WriteAllHistogramsToFile(fout, false);
   }

   // Plotting If you want for 1D
   /*   for(int i=0; i< variablesMC.size();i++){
     PlotCVAndError(variablesData[i]->data_reco.hist,variablesMC[i]->mc_reco.hist,variablesMC[i]->GetName(),MCscale);
       
     PlotErrorSummary(variablesMC[i]->mc_reco.hist, variablesMC[i]->GetName());
     PlotStacked(variablesData[i]->data_reco_sb.hist,variablesMC[i]->mc_sb.GetHistArray(),MCscale, variablesMC[i]->mc_sb.GetName(), variablesMC[i]->mc_sb.GetName());
   }//End 1D plotting 
   
   for (auto v : variables2DMC) {
     v->WriteAllHistogramsToFile(fout,true);
   }

   for (auto v : variables2DTruth) {
     v->WriteAllHistogramsToFile(fout,false);
   }
   */

  fout.cd();
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
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

  std::vector<double> mcRunBin;
  
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
    mcRunBin = binsDef->GetEnergyBins("mc_run"); 
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

  Var* mcRun = new Var("mc_run", "mc_run", mcRunBin, &CVUniverse::GetMCRunN, &CVUniverse::GetMCRunN);

  variables = {mcRun}; //{enu,ehad}; 

  // 2D Variables 
  Var2D* pTmu_pZmu = new Var2D(*pTmu, *pZmu);
  Var2D* W_Q2 = new Var2D(*W, *Q2);
  Var2D* enu_ehad = new Var2D(*enu, *ehad);
  Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
  Var2D* x_y = new Var2D(*x, *y);  // y var
  Var2D* x_Q2 = new Var2D(*x, *Q2);  // y var
  
  variables2d = {};
   
  //smakefor (auto v : variables2d) v->InitializeAllHistograms(error_bands);
  for (auto v : variables2d) v->InitializeAllHistograms(error_bands);
  for (auto v : variables) v->InitializeAllHistograms(error_bands);

   int mc0=0;
   int mc1=0;
   int mc2=0;
   int mc3=0;
     int mc_truth0=0.0;
     int mc_truth1=0.0;

   
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

          if(!cutter->PassiveTargetRegionTrue(universe)) continue; // true region
          //std::cout<< "MC run " << universe->GetInt("mc_run") << std::endl;

          mc1++;
          for (auto v : variables){
            if(cutter->PassiveTargetRegion(universe)){
              v->m_selected_mc_MLreco.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetWeight());
              mc2++;
            }
            if(cutter->PassiveTargetRegionTBV(universe)){
              v->m_selected_mc_reco.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetWeight());
              mc3++;
            }
          }
        } 
        else if(!isMC ){ // DENOMINATOR
            mc_truth0++;

            if(!cutter->PassiveTargetRegionTrue(universe)) continue; // true tracker
            mc_truth1++;   

          for (auto v : variables){
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
     std::cout<<" All events = "<<mc0<<std::endl;
     std::cout<<" True passive target region = "<<mc1<<std::endl;
     std::cout<<" True + ML reco passive target region = "<<mc2<<std::endl;
     std::cout<<" True + TBV reco passive target region = "<<mc3<<std::endl;
   std::cout<<"**********************************"<<std::endl;
     std::cout<<"Printing the Truth Summary "<<std::endl;
     std::cout<<" All true events = "<<mc_truth0<<std::endl;
     std::cout<<" True passive target region = "<<mc_truth1<<std::endl;
   //	return variables;
}
//============================================================================================================================
// Main
