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

// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we
// shield it.
#ifndef __CINT__
#include "plotting_functions.h"
#endif
//using namespace globalV;
using namespace NUKECC_ANA;

//============================================================================================================================
// Main
int main(int argc, char * argv[]){
   ROOT::Cintex::Cintex::Enable();
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
   
   ROOT::Cintex::Cintex::Enable();
   TH1::AddDirectory(false);
   
     
   typedef VarLoop::Variable Var;
   typedef Var2DLoop::Variable2D Var2D;
     
  typedef PlotUtils::HistWrapper<NUKECC_ANA::CVUniverse> HW;
  typedef PlotUtils::MnvH1D MH1D;

  typedef PlotUtils::Hist2DWrapper<NUKECC_ANA::CVUniverse> HW2D;  
  typedef PlotUtils::MnvH2D MH2D;   


   NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(playlist);
  
     
    std::vector<Var*> variables,variablesData; 
     std::vector<Var2D*> variables2DMC;    

   	
    // We have to delete the events ourselves. Should use some sort of
    // smart pointer to fix that
TString hname1 = Form("Hists_EventSelection_t1_z26_Nu_v1_.root");
     double total_pot_data,total_pot_mc;
    utils->getPOT(total_pot_data,total_pot_mc);  
    double  MCscale=total_pot_data/total_pot_mc;
  
TFile fin("Hists_EventSelection_t1_z26_Nu_v1_.root", "READ");
// TFile *f1 = new TFile(hname1,"read" );
  fin.cd();
    MH1D *MCreco,*datareco,*datasb;  
    MH2D *MCreco2D,*datareco2D;    

    MCreco  = (MnvH1D*)fin.Get("selected_mc_reco_Enu"); 
    datareco  = (MnvH1D*)fin.Get("selected_data_reco_Enu"); 
     
   
    MCreco  = (MnvH1D*)fin.Get("selected_mc_reco_Enu"); 
    datareco  = (MnvH1D*)fin.Get("selected_data_reco_Enu"); 
    datasb  = (MnvH1D*)fin.Get("selected_data_reco_sb_Enu"); 
    
    
    MCreco2D  = (MnvH2D*)fin.Get("selected_mc_reco_Enu_Ehad"); 
    datareco2D  = (MnvH2D*)fin.Get("selected_data_reco_Enu_Ehad"); 
    

      PlotCVAndError(datareco,MCreco,"Enu",MCscale);       
       PlotErrorSummary(MCreco,"Enu");
  // Stacking
  PlotUtils::HistFolio<PlotUtils::MnvH1D> Enu(PlotUtils::LoadHistFolioFromFile<PlotUtils::MnvH1D>(fin, std::string("selected_mc_sb_Enu")), "Enu");
  PlotStacked(datasb,Enu.GetHistArray(),MCscale, Enu.GetName(), Enu.GetName());
   
       Plot2D(MCreco2D,"Enu_vs_Ehad", "Enu","Ehad");

  // PlotStacked(m_selected_data_reco_sb.hist,m_selected_mc_sb.GetHistArray(),MCscale,
    //            m_selected_mc_sb.GetName(), m_selected_mc_sb.GetName());
//	}

/*
   for(int i=0; i< variables2DMC.size();i++){
       Plot2D(variables2DMC[i]->m_selected_mc_reco2d.hist, variables2DMC[i]->GetName(), variables2DMC[i]->GetNameX(),
           variables2DMC[i]->GetNameY());

       Plot2D(variables2DMC[i]->m_selected_data_reco2d.hist, variables2DMC[i]->GetName(), variables2DMC[i]->GetNameX(),
           variables2DMC[i]->GetNameY());
*/
//	}

}//End Main
