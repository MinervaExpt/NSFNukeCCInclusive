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
//#include "include/CVUniverse_faiza.h"
#include "../../../NUKECCSRC/include/CVUniverse.h"
#include "../../include/Variable_plasticSB.h"
#include "../../include/playlists/playlists.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../../NUKECCSRC/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "../../../NUKECCSRC/include/LateralSystematics.h"
#include <iostream>
#include <stdlib.h>
//#include "Cintex/Cintex.h"
#include "../../../NUKECCSRC/include/NukeCCUtilsNSF.h"
#include "../../../NUKECCSRC/include/NukeCC_Cuts.h"
#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
//#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "TParameter.h"

// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we shield it.
#ifndef __CINT__
#include "../../include/plotting_functions.h"
#endif
#include "PlotUtils/MacroUtil.h" 
//using namespace globalV;
using namespace NUKECC_ANA;
//======================================================================
typedef VarLoop::Variable Var;

void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables, bool isMC, int targetID, int targetZ, const string playlist, bool doDIS=false);
 
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
    PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);//NSFDefaults::beamThetaX_Err,NSFDefaults::beamThetaY_Err);
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
// Main
int main(int argc, char *argv[]){
   //ROOT::Cintex::Cintex::Enable();
   TH1::AddDirectory(false);
	
   if(argc==1){
     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
     std::cout<<"MACROS HELP:\n\n"<<
       "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n"<<
       "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
       "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" <<
       "\t-Playlist\t \t =\t eg. minervame1A"<< std::endl;
     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
     return 0;
   }

   TString dir(argv[1]);
   int targetID = atoi(argv[2]);
   int targetZ = atoi(argv[3]);
   const string playlist= argv[4];
   
   bool doDIS=false;
  //const std::string mc_file_list( get_mc_files(playlist, targetID) );
  //const std::string data_file_list( get_data_files(playlist) );
  
  const std::string mc_file_list("../../include/playlists/NukeCC_MC_minervame6A_MuonKludged.txt");
  const std::string data_file_list("../../include/playlists/NukeCC_Data_minervame6A_MuonKludged.txt");
 
  //const std::string plist_string("minervame1A");
   const std::string plist_string(playlist);
   const std::string reco_tree_name("NukeCC");
   const bool wants_truth = false;
   const bool is_grid = false;

   PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth, is_grid);

   util.PrintMacroConfiguration("main");
  //=========================================
  // Systematics
  //=========================================

   PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
   PlotUtils::MinervaUniverse::SetNuEConstraint(true);
   PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
   PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
   PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
   PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);


   NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
   NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
   NukeCC_Binning  *binsDef = new NukeCC_Binning();
      
   PlotUtils::ChainWrapper* chainData = util.m_data;
   PlotUtils::ChainWrapper* chainMC   = util.m_mc;
   HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);
   double DataPot=  util.m_data_pot; 
   double MCPot=  util.m_mc_pot; 
   double  MCscale=DataPot/MCPot;
   std::cout<<" MCScale= "<<MCscale<<std::endl; 

   std::vector<Var*> variablesMC,variablesData; 

   TString histFileName;
   if(RunCodeWithSystematics){
       histFileName = utils->GetHistFileName( "PlasticSideband", FileType::kAny, targetID, targetZ, helicity ); 
     }
   else{
       histFileName = utils->GetHistFileName( "PlasticBackgd_without", FileType::kAny, targetID, targetZ, helicity );
     } 
   TFile fout(dir.Append(histFileName),"RECREATE");	
   
   // For 1D variables 
   FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,true,targetID, targetZ, plist_string,doDIS);
       
   for (auto v : variablesMC){
      v->m_hists_US_Fe.SyncCVHistos();
      v->m_hists_US_Pb.SyncCVHistos();
      v->m_hists_US_C.SyncCVHistos();
      v->m_hists_US_other.SyncCVHistos();
      v->m_hists_US_regUS.SyncCVHistos();
      v->m_hists_US_regDS.SyncCVHistos();

      v->m_hists_DS_Fe.SyncCVHistos();
      v->m_hists_DS_Pb.SyncCVHistos();
      v->m_hists_DS_C.SyncCVHistos();
      v->m_hists_DS_other.SyncCVHistos();
      v->m_hists_DS_regUS.SyncCVHistos();
      v->m_hists_DS_regDS.SyncCVHistos();

      v->m_hists_tgt_Fe.SyncCVHistos();
      v->m_hists_tgt_Pb.SyncCVHistos();
      v->m_hists_tgt_C.SyncCVHistos();
      v->m_hists_tgt_other.SyncCVHistos();
      v->m_hists_tgt_regUS.SyncCVHistos();
      v->m_hists_tgt_regDS.SyncCVHistos();
   }

   FillVariable(chainData, helicity, utils, cutter,binsDef,variablesData,false,targetID, targetZ, plist_string,doDIS);
 
   for (auto v : variablesData){
     v->m_selected_data_reco_US.SyncCVHistos();
     v->m_selected_data_reco_DS.SyncCVHistos();
     v->m_selected_data_reco_tgt.SyncCVHistos();
  }
 
   for (auto v : variablesMC) {
     v->WriteAllHistogramsToFile(fout, true);
   }

   for (auto v : variablesData) {
     v->WriteAllHistogramsToFile(fout, false);
   }

  //Writing POT to the HistFile
  fout.cd();
  auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  auto Scale = new TParameter<double>("DataMCScale", MCscale);
  dataPOTOut->Write();
  Scale->Write();
  mcPOTOut->Write(); 


/*
   // Plotting If you want for 1D
      for(int i=0; i< variablesMC.size();i++){
       
      PlotStacked(variablesData[i]->m_selected_data_reco_US.hist, variablesMC[i]->m_selected_mc_US.GetHistArray(), MCscale, variablesMC[i]->m_selected_mc_US.GetName(), variablesMC[i]->m_selected_mc_US.GetName());
      PlotStacked(variablesData[i]->m_selected_data_reco_DS.hist, variablesMC[i]->m_selected_mc_DS.GetHistArray(), MCscale, variablesMC[i]->m_selected_mc_DS.GetName(), variablesMC[i]->m_selected_mc_DS.GetName());
      PlotStacked(variablesData[i]->m_selected_data_reco_tgt.hist, variablesMC[i]->m_selected_mc_tgt.GetHistArray(), MCscale, variablesMC[i]->m_selected_mc_tgt.GetName(), variablesMC[i]->m_selected_mc_tgt.GetName());
	
   }//End 1D plotting 
*/
cout<<"End of the macro"<<endl;
}//End Main

void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){

  std::map<std::string, std::vector<CVUniverse*> > error_bands =GetErrorBands(chain);

  std::vector<double> Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin;
  std::vector<double> planeDNNbin;

   if (doDIS){
     Enubin = binsDef->GetDISBins("Enu"); 
     Emubin = binsDef->GetDISBins("Emu"); 
     Ehadbin = binsDef->GetDISBins("Ehad");
     xbin = binsDef->GetDISBins("x");
     ybin = binsDef->GetDISBins("y");
     Q2bin = binsDef->GetDISBins("Q2");
     Wbin = binsDef->GetDISBins("W");
     planeDNNbin = binsDef->GetDISBins("planeDNN");
     }
   else{
     Enubin = binsDef->GetEnergyBins("Enu"); 
     Emubin = binsDef->GetEnergyBins("Emu"); 
     Ehadbin = binsDef->GetEnergyBins("Ehad");
     xbin = binsDef->GetEnergyBins("x");
     ybin = binsDef->GetEnergyBins("y");
     planeDNNbin = binsDef->GetEnergyBins("planeDNN");
   }
   Q2bin = binsDef->GetSidebandBins("Q2");
   Wbin = binsDef->GetSidebandBins("W");

   Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
   Var* ehad = new Var("Ehad", "Ehad (GeV)", Ehadbin, &CVUniverse::GetEhadGeV, &CVUniverse::GetEhadTrueGeV);
   Var* Q2 = new Var("Q2", "Q2 (GeV^2)", Q2bin, &CVUniverse::GetQ2RecoGeV, &CVUniverse::GetQ2TrueGeV);
   Var* W = new Var("W", "W (GeV)", Wbin, &CVUniverse::GetWRecoGeV, &CVUniverse::GetWTrueGeV);
   Var* emu = new Var("Emu", "Emu (GeV)", Emubin, &CVUniverse::GetMuonEGeV, &CVUniverse::GetMuonETrueGeV);
   Var* x = new Var("x", "x", xbin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
   Var* y = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);
   Var* planeDNN = new Var("planeDNN", "planeDNN", planeDNNbin, &CVUniverse::GetplaneDNNReco, &CVUniverse::GetplaneDNNTrue);
   
   variables = {planeDNN, enu, ehad, emu, x, y, W, Q2}; 
      
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


   for (auto v : variables) v->InitializeAllHistograms(error_bands,targetID,targetZ);


   int reco0=0;
   int reco00=0;
   int reco1=0;
   int reco2=0;
   int reco3=0;
   int reco4=0; 
   int reco5=0; 
   int allcuts=0;
   int us_1=0;
   int us_2=0;
   int us_3=0;
   int us_4=0;
   int us_5=0;
   int us_6=0;
   int us_7=0;
   int us_8=0;
   int us_9=0;
   int us_10=0;
   int ds_1=0;
   int ds_2=0;
   int ds_3=0;
   int ds_4=0;
   int ds_5=0;
   int ds_6=0;
   int ds_7=0;
   int ds_8=0;
   int tgt_1=0;
   int tgt_2=0;
   int tgt_3=0;
   int tgt_4=0;
   int tgt_5=0;
   int tgt_6=0;
   int end_1=0; 

   int us_fe=0;
   int us_pb=0;
   int us_c=0;
   int us_other=0;
   int us_regus=0;
   int us_regds=0;
   int ds_fe=0;
   int ds_pb=0;
   int ds_c=0;
   int ds_other=0;
   int ds_regus=0;
   int ds_regds=0;
   int data_US=0;
   int data_DS=0;
   int data_tgt=0;


   CVUniverse *dataverse = new CVUniverse(chain,0);
    
   int n = 0; //how far away do we want to cut the sideband from planes immediately US/DS the passive target
   //=========================================
   // Entry Loop
   //=========================================
   for(int target = 0;target<targetIDs.size();target++){    //target loop
   std::cout<<"# of entries = "<<chain->GetEntries()<<std::endl;
   for(int i=0; i<chain->GetEntries(); ++i){
     if(i%500000==0) std::cout << (i/1000) << "k " << std::endl;
     //=========================================
     // For every systematic, loop over the universes, and fill the
     // appropriate histogram in the MnvH1D
     //=========================================

       reco00++; //number of entries before any cuts

       if(isMC){     
       for (auto band : error_bands){
       std::vector<CVUniverse*> error_band_universes = band.second;
       for (auto universe : error_band_universes){
	 // Tell the Event which entry in the TChain it's looking at
	 //=========================================
	 // CUTS in each universe
	 //========================================
	 universe->SetEntry(i);
         reco0++; //number of entries before any cuts
	 if(!cutter->PassReco(universe,helicity)) continue;
	 reco1++; // no. of entries after reco cuts
	 if(!cutter->IsInMaterial(universe,targetIDs[target],targetZ, /*anyTrakerMod*/false)) continue;
	 reco2++; // no. of entries after Material cut
	 if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;	   
	 reco3++; // no. of entries after ML prob cut

           if( universe->Var("planeDNN" ) < universe->GetTargetMinBin(targetIDs[target])+n ) continue;
           if( universe->Var("planeDNN" ) > universe->GetTargetMaxBin(targetIDs[target])-(n+1) ) continue;
           if( targetIDs[target] == 1 && universe->Var("planeDNN") <= universe->GetTargetUSPlane(targetIDs[target]) ) continue;
           reco4++; //number of entries after Plastic SB cut

	   for (auto v : variables){  //1D variable loop
	     if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(universe)) continue;
           //  reco5++;
	     if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(universe))continue;
           //  reco6++;

	    //if(!cutter->PassDISCut(universe)) continue;
//------------------------------------------------------------PLASTIC SB STUFF START--------------------------------------------------------------------------------------------------------------------
            int planeVal = universe->Var("planeDNN");
            double VtxZ = universe->GetVecElem("mc_vtx",2); 
            int targetA = universe->GetAtomicNumber();
            const double wgt = universe->GetWeight();
            if( planeVal <= universe->GetTargetUSPlane(targetIDs[target]) && planeVal >= universe->GetTargetDSPlane(targetIDs[target]-1) && targetIDs[target] < 10 ){ //-----------1   
                     us_1++;
               if( planeVal != universe->GetTargetUSPlane(targetIDs[target]) && planeVal != universe->GetTargetDSPlane(targetIDs[target]-1) && planeVal != universe->GetTargetMinBin(targetIDs[target])+n && targetIDs[target] < 10 ){
                             us_2++;
                  if( universe->GetTargetZStart(targetIDs[target]) <= VtxZ && VtxZ <= universe->GetTargetZEnd(targetIDs[target]) ){ //if its upstream bg, consider what the actual target was //--------------1a
                            us_3++;
                   if(cutter->IsInTrueMaterial(universe, targetIDs[target], targetZ)){
                     if(targetA == 56 || targetA == 55){ us_fe++; 
                        //v->m_selected_mc_US.GetComponentHist("Fe")->Fill(v->GetRecoValue(*universe, 0), wgt);
                        v->m_hists_US_Fe.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt); 
                      }
                     else if(targetA == 207){   us_pb++; 
                        //v->m_selected_mc_US.GetComponentHist("Pb")->Fill(v->GetRecoValue(*universe, 0), wgt);
                        v->m_hists_US_Pb.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                      }
                     else if(targetA == 12 || targetA == 14 || targetA == 16){ us_c++; 
                        //v->m_selected_mc_US.GetComponentHist("C")->Fill(v->GetRecoValue(*universe, 0), wgt);
                        v->m_hists_US_C.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                      }
                     else std::cout << "couldn't find target material for target " << targetIDs[target] << "and Z " << targetZ <<  "; true atomic mass is = " << targetA << std::endl;
                   }
                 else{  us_other++; 
                       //v->m_selected_mc_US.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                       v->m_hists_US_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                     }
              } //----------------1a   

                else if( universe->GetTargetZStart(targetIDs[target]) > VtxZ ){ //--------------1b
                       us_4++;
                 if( targetA == 56 ){   us_other++;
                      //v->m_selected_mc_US.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                      v->m_hists_US_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                  }
                 else if(targetA == 207){ us_other++;
                      //v->m_selected_mc_US.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                      v->m_hists_US_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                  }
                 else if( cutter->IsInTrueMaterial(universe, 3, 6) || targetA == 16 || targetA == 14 ){ us_other++;
                      //v->m_selected_mc_US.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                      v->m_hists_US_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                  }
                else{  us_regus++; 
                      //v->m_selected_mc_US.GetComponentHist("regUS")->Fill(v->GetRecoValue(*universe, 0), wgt);
                      v->m_hists_US_regUS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                  }
              } //-----------------1b

               else if( universe->GetTargetZEnd(targetIDs[target]) < VtxZ ){ //--------------1c
                    us_5++;
                 if( targetA == 56 ){ us_other++; 
                     //v->m_selected_mc_US.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt); 
                     v->m_hists_US_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                  }
                 else if(targetA == 207){ us_other++; 
                     //v->m_selected_mc_US.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_US_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                  }
                 else if( cutter->IsInTrueMaterial(universe, 3, 6) || targetA == 16 || targetA == 14 ){ 
 us_other++;                     //v->m_selected_mc_US.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_US_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                 }
                else{ us_regus++; 
                     //v->m_selected_mc_US.GetComponentHist("regDS")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_US_regDS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                 }
              } //-----------------1c
              else{     us_6++; 
		    cout << "Missing Event in filling US sideband region!" << endl;
               }                    
            }
         }   //------------------1

      else if( planeVal <= universe->GetTargetUSPlane(targetIDs[target]+1) && planeVal >= universe->GetTargetDSPlane(targetIDs[target]) && targetIDs[target] < 10 ){  //-----------------2
                  ds_1++;
             if( planeVal != universe->GetTargetUSPlane(targetIDs[target]+1) && planeVal != universe->GetTargetDSPlane(targetIDs[target]) && planeVal != universe->GetTargetMaxBin(targetIDs[target])-(n+1) && targetIDs[target] < 10  ){
                    ds_2++;
               if( universe->GetTargetZStart(targetIDs[target]) <= VtxZ && VtxZ <= universe->GetTargetZEnd(targetIDs[target]) ){ //if its upstream bg, consider what the actual target was //--------------2a
                        //ds_3++;
                 if(cutter->IsInTrueMaterial(universe, targetIDs[target], targetZ)){
                   if(targetA == 56 || targetA == 55){   ds_fe++; 
                     //v->m_selected_mc_DS.GetComponentHist("Fe")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_DS_Fe.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                   }
                   else if(targetA == 207){   ds_pb++;
                     //v->m_selected_mc_DS.GetComponentHist("Pb")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_DS_Pb.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                   }
                   else if(targetA == 12 || targetA == 14 || targetA == 16){   ds_c++;
                     //v->m_selected_mc_DS.GetComponentHist("C")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_DS_C.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                   }
                   else std::cout << "couldn't find target material for target " << targetIDs[target] << "and Z " << targetZ <<  "; true atomic mass is = " << targetA << std::endl;
                 }
                 else{   ds_other++;
                    //v->m_selected_mc_DS.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                    v->m_hists_DS_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                 }
              } //----------------2a   
                else if( universe->GetTargetZStart(targetIDs[target]) > VtxZ ){ //--------------2b
                    ds_4++;
                 if( targetA == 56 ){   ds_other++; 
                     //v->m_selected_mc_DS.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_DS_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                   }
                 else if(targetA == 207){  ds_other++;
                     //v->m_selected_mc_DS.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_DS_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                   }
                 else if( cutter->IsInTrueMaterial(universe, 3, 6) || targetA == 16 || targetA == 14 ){    ds_other++;
                     //v->m_selected_mc_DS.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_DS_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                   }
                else{  ds_regus++;  
                     //v->m_selected_mc_DS.GetComponentHist("regUS")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_DS_regUS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                   }
              } //-----------------2b
               else if( universe->GetTargetZEnd(targetIDs[target]) < VtxZ ){ //--------------2c
                   ds_5++;
                 if( targetA == 56 ){    ds_other++;
                     //v->m_selected_mc_DS.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_DS_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                   }
                 else if(targetA == 207){     ds_other++;
                     //v->m_selected_mc_DS.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_DS_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                   }
                 else if( cutter->IsInTrueMaterial(universe, 3, 6) || targetA == 16 || targetA == 14 ){    ds_other++;
                     //v->m_selected_mc_DS.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_DS_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                   }
                else{  ds_regds++; 
                     //v->m_selected_mc_DS.GetComponentHist("regDS")->Fill(v->GetRecoValue(*universe, 0), wgt);
                     v->m_hists_DS_regDS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                }
            } //-----------------2c
               else{   ds_6++;  
                    cout << "Missing Event in filling DS sideband region!" << endl;
               }     
           }             
       }//---------------2
      else if(planeVal == universe->GetTargetPlane(targetIDs[target]) && targetIDs[target] < 10){ //--------------3                
             tgt_1++;
              if( universe->GetTargetZStart(targetIDs[target]) <= VtxZ && VtxZ <= universe->GetTargetZEnd(targetIDs[target]) ){ //if its upstream bg, consider what the actual target was  //--------------3a
                        tgt_2++;
                 if(cutter->IsInTrueMaterial(universe, targetIDs[target], targetZ)){
                     if(targetA == 56 || targetA == 55){ 
                          //v->m_selected_mc_tgt.GetComponentHist("Fe")->Fill(v->GetRecoValue(*universe, 0), wgt);
                          v->m_hists_tgt_Fe.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                       }
                     else if(targetA == 207){ 
                          //v->m_selected_mc_tgt.GetComponentHist("Pb")->Fill(v->GetRecoValue(*universe, 0), wgt);
                          v->m_hists_tgt_Pb.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                       }
                     else if(targetA == 12 || targetA == 14 || targetA == 16){ 
                          //v->m_selected_mc_tgt.GetComponentHist("C")->Fill(v->GetRecoValue(*universe, 0), wgt);
                          v->m_hists_tgt_C.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt); 
                       }
                     else std::cout << "couldn't find target material for target " << targetIDs[target] << "and Z " << targetZ <<  "; true atomic mass is = " << targetA << std::endl;
                 }
                     else{ 
                          //v->m_selected_mc_tgt.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                          v->m_hists_tgt_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                      }
              } //----------------3a   
                else if( universe->GetTargetZStart(targetIDs[target]) > VtxZ ){  //--------------3b
                        tgt_3++;
                 if( targetA == 56 ){ 
                          //v->m_selected_mc_tgt.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                          v->m_hists_tgt_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                      }
                 else if(targetA == 207){ 
                          //v->m_selected_mc_tgt.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                          v->m_hists_tgt_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                      }
                 else if( cutter->IsInTrueMaterial(universe, 3, 6) || targetA == 16 || targetA == 14 ){ 
                          //v->m_selected_mc_tgt.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                          v->m_hists_tgt_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                      }
                else{ 
                          //v->m_selected_mc_tgt.GetComponentHist("regUS")->Fill(v->GetRecoValue(*universe, 0), wgt);
                          v->m_hists_tgt_regUS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                  }
              } //-----------------3b
               else if( universe->GetTargetZEnd(targetIDs[target]) < VtxZ ){  //--------------3c
                      tgt_4++;
                 if( targetA == 56 ){ 
                       //v->m_selected_mc_tgt.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                       v->m_hists_tgt_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                    }
                 else if(targetA == 207){ 
                        //v->m_selected_mc_tgt.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                        v->m_hists_tgt_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                    }
                 else if( cutter->IsInTrueMaterial(universe, 3, 6) || targetA == 16 || targetA == 14 ){ 
                        //v->m_selected_mc_tgt.GetComponentHist("Other")->Fill(v->GetRecoValue(*universe, 0), wgt);
                        v->m_hists_tgt_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                    }
                else{ 
                        //v->m_selected_mc_tgt.GetComponentHist("regDS")->Fill(v->GetRecoValue(*universe, 0), wgt);
                        v->m_hists_tgt_regDS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                   }
              } //-----------------3c
               else{  tgt_5++; 
                      cout << "Missing Event in filling Target region!" << endl;
               }
       }//---------------3
      else{    end_1++;
             cout << "Missing Event!" << endl;
             cout << "z vertex, target A, plane = " << universe->GetVecElem("mc_vtx",2) << ", " << targetA << ", " << universe->Var("planeDNN", true) << endl;
       }
//------------------------------------------------------------PLASTIC SB STUFF END-------------------

         } //End 1D variables loop
       } // End band's universe loop
     }// End Band loop
   }//End isMC loop

else{
        dataverse->SetEntry(i);
         reco00++;
	 if(!cutter->PassReco(dataverse,helicity)) continue;
	 reco1++;
	 if(!cutter->IsInMaterial(dataverse,targetIDs[target],targetZ, false)) continue;
	 reco2++;
	 if( dataverse->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;	   
	   //trial for matching with NukeCC
           //if(targetIDs[target]<10 && dataverse->GetInt("NukeCC_targetID") != targetIDs[target]) continue;
         reco3++;
           if( dataverse->Var("planeDNN" ) < dataverse->GetTargetMinBin(targetIDs[target])+n ) continue;
           if( dataverse->Var("planeDNN" ) > dataverse->GetTargetMaxBin(targetIDs[target])-(n+1) ) continue;
           if( targetIDs[target] == 1 && dataverse-> Var("planeDNN") <= dataverse->GetTargetUSPlane(targetIDs[target]) ) continue; //NEW
         reco0++;
	   for (auto v : variables){
	     if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(dataverse)) continue;
             reco4++;
	     if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(dataverse))continue;
             reco5++;
 	     //if(!cutter->PassDISCut(dataverse)) continue;
//--------------------------------------------------------Plactic SB Data Histograms--------------------
      int planeVal = dataverse->Var("planeDNN");
 			
      if( planeVal <= dataverse->GetTargetUSPlane(targetIDs[target]) && planeVal >= dataverse->GetTargetDSPlane(targetIDs[target]-1) && targetIDs[target] < 10 ){ //-----------1
             if( planeVal != dataverse->GetTargetUSPlane(targetIDs[target]) && planeVal != dataverse->GetTargetDSPlane(targetIDs[target]-1) && planeVal != dataverse->GetTargetMinBin(targetIDs[target])+n && targetIDs[target] < 10 ){
              data_US++;
                  v->m_selected_data_reco_US.hist->Fill(v->GetRecoValue(*dataverse, 0));
              }
       }   //------------------1

      else if( planeVal <= dataverse->GetTargetUSPlane(targetIDs[target]+1) && planeVal >= dataverse->GetTargetDSPlane(targetIDs[target]) && targetIDs[target] < 10 ){  //-----------------2
             if( planeVal != dataverse->GetTargetUSPlane(targetIDs[target]+1) && planeVal != dataverse->GetTargetDSPlane(targetIDs[target]) && planeVal != dataverse->GetTargetMaxBin(targetIDs[target])-(n+1) && targetIDs[target] < 10  ){
              data_DS++;
                  v->m_selected_data_reco_DS.hist->Fill(v->GetRecoValue(*dataverse, 0));
              }             
       }//---------------2
      else if(planeVal == dataverse->GetTargetPlane(targetIDs[target]) && targetIDs[target] < 10){ //--------------3
                data_tgt++;
                  v->m_selected_data_reco_tgt.hist->Fill(v->GetRecoValue(*dataverse, 0));
       }//---------------3
      else{
             cout << "PlasticBG:: DATA Missing Event!" << endl;
       }
//--------------------------------------------------------Plactic SB Data Histograms end--------------------
	
        	}//End 1D var loop       
         } //End else 
   }//End entries loop
 }//End target loop 

for(auto band : error_bands){
    std::vector<CVUniverse*> band_universes = band.second;
    for(unsigned int i_universe = 0; i_universe < band_universes.size(); ++i_universe)
      delete band_universes[i_universe];
  } 
      delete dataverse;

   std::cout<<"**********************************"<<std::endl;
   std::cout<<"Printing the ";
      isMC? std::cout<<"MC": std::cout<<"Data";
   std::cout<<" Summary "<<std::endl;
   std::cout<<" No1 cuts = "<<reco00<<std::endl;
   std::cout<<" Reco Cut = "<<reco1<<std::endl;
   std::cout<<" Material Cut = "<<reco2<<std::endl;
   std::cout<<" Plane/targetID Cuts = "<<reco3<<std::endl;
   std::cout<<" SB only = "<<reco4<<std::endl;
   std::cout<<" Muon energy Cut = "<<reco4<<std::endl;
   std::cout<<" Muon theta Cut = "<<reco5<<std::endl;
   std::cout<<"**********************************"<<std::endl;
   std::cout<<" US first cut = "<<us_1<<std::endl;
   std::cout<<" US second cut = "<<us_2<<std::endl;
   std::cout<<" US Fe = "<<us_fe<<std::endl;
   std::cout<<" US Pb = "<<us_pb<<std::endl;
   std::cout<<" US C = "<<us_c<<std::endl;
   std::cout<<" US Other = "<<us_other<<std::endl;
   std::cout<<" US regUS = "<<us_regus<<std::endl;
   std::cout<<" US regDS = "<<us_regds<<std::endl;
   std::cout<<" DS first cut = "<<ds_1<<std::endl;
   std::cout<<" DS second cut = "<<ds_2<<std::endl;
   std::cout<<" DS Fe = "<<ds_fe<<std::endl;
   std::cout<<" DS Pb = "<<ds_pb<<std::endl;
   std::cout<<" DS C = "<<ds_c<<std::endl;
   std::cout<<" DS Other = "<<ds_other<<std::endl;
   std::cout<<" DS regUS = "<<ds_regus<<std::endl;
   std::cout<<" DS regDS = "<<ds_regds<<std::endl;
   std::cout<<"**********************************"<<std::endl;
   std::cout<<" Data US = "<<data_US<<std::endl;
   std::cout<<" Data DS = "<<data_DS<<std::endl;
   std::cout<<" Data Tgt = "<<data_tgt<<std::endl;
   std::cout<<"**********************************"<<std::endl;

}
//============================================================================================================================
// Main

