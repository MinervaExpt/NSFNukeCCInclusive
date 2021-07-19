#ifndef VARIABLE_H
#define VARIABLE_H

#include <iterator>
#include "include/CVUniverse.h"
//#include "Histogram.h"
#include "PlotUtils/HistFolio.h"
#include "PlotUtils/HistWrapper.h"

#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MnvH2D.h"
#ifndef __CINT__  // CINT doesn't know about std::function
#include "PlotUtils/VariableBase.h"
#include <bitset>
#include "PlotUtils/Variable2DBase.h"
#include "MinervaUnfold/MnvResponse.h"

#include "PlotUtils/AnaBinning.h"

#endif  // __CINT__

namespace VarLoop {

class Variable : public PlotUtils::VariableBase<NUKECC_ANA::CVUniverse> {
 private:
  typedef PlotUtils::HistWrapper<NUKECC_ANA::CVUniverse> HW;
  typedef PlotUtils::MnvH1D MH1D;

 public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  Variable(ARGS... args) : PlotUtils::VariableBase<NUKECC_ANA::CVUniverse>(args...) {}

  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // HISTWRAPPER
  // selected mc reco histwrapper
  HW m_selected_mc_reco,m_selected_data_reco,m_selected_data_reco_sb, m_selected_mc_true;

  typedef PlotUtils::Hist2DWrapper<NUKECC_ANA::CVUniverse> HW2D;
  HW2D mresp1D;
  std::map<std::string,MinervaUnfold::MnvResponse*> Response1D;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr2; 

  MnvH2D *migrationH2D = NULL;
  MnvH1D *h_reco1D = NULL;
  MnvH1D *h_truth1D = NULL;
 // HISTFOLIO
  // selected mc reco - signal background histfolio
  PlotUtils::HistFolio<MH1D> m_selected_mc_sb;
  // PlotUtils::MH1D* m_selected_data_sb;
  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHistograms(T univs) {
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    const bool clear_bands = true;  // we want empty histograms

    // HISTWRAPPER
    // selected mc reco histwrapper
    MH1D* dummy_selected_mc_reco = new MH1D(Form("selected_mc_reco_%s", name), name,
                                      GetNBins(), bins.data());
    m_selected_mc_reco = HW(dummy_selected_mc_reco, univs, clear_bands);
    //for 1D analysis
    MH1D* dummy_selected_mc_true = new MH1D(Form("selected_mc_true_%s", name), name,
                                      GetNBins(), bins.data());
    m_selected_mc_true = HW(dummy_selected_mc_true, univs, clear_bands);
    
    MH1D* dummy_selected_data_reco = new MH1D(Form("selected_data_reco_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco = HW(dummy_selected_data_reco, univs, clear_bands);


    MH1D* selected_data_reco_sb = new MH1D(Form("selected_data_reco_sb_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco_sb = HW(selected_data_reco_sb, univs, clear_bands);
  
  // HISTFOLIO
    // selected mc reco - signal background histfolio
    
   m_selected_mc_sb = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("selected_mc_sb_%s", name), name, GetNBins(), bins.data());
    
 // PlotUtils::MnvH1D* data = new PlotUtils::MnvH1D(
   //   "dummy", "dummy", plotting::nbins, plotting::xmin, plotting::xmax);
  //  m_selected_data_sb = PlotUtils::HistFolio<PlotUtils::MnvH1D>(
    //    Form("selected_data_sb_%s", name), name, GetNBins(), bins.data());
   
    m_selected_mc_sb.AddComponentHist("DIS");
    m_selected_mc_sb.AddComponentHist("MC");
   // m_selected_data_sb.AddComponentHist("Data");
delete dummy_selected_mc_reco;
delete dummy_selected_mc_true;
delete dummy_selected_data_reco;
  }

void SetupResponse1D(std::map<const std::string, int> systematics){
//void SetupResponse(T univs){
	   std::vector<double> bins = GetBinVec();
	   const char* name = GetName().c_str();

	   axis_binning bin_x;
	   bin_x.uniform=false;
	   	
     	   vector<double> vx;
	   //for(int i=0; i<GetNBins()+1; i++){vx.push_back(GetBinVec().data()[i+1]-GetBinVec().data()[i]);}
	   for(int i=0; i<=GetNBins(); i++){vx.push_back(GetBinVec().data()[i]);}
     	   
	   bin_x.bin_edges = vx;
	   bin_x.nbins	   = GetNBins();
	   bin_x.min 	   = GetBinVec().data()[0];
	   bin_x.max       = GetBinVec().data()[GetNBins()];
           
           cout<<"bins min = "<<bin_x.min<<"bins max = "<<bin_x.max<<"Number of bins = "<<bin_x.nbins<<endl;
           
	   Response1D.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("selected_mc_response1d_%s", name), name, bin_x, bin_x, systematics)));
	   //Response1D.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("selected_mc_response1d_%s", name), name, GetNBins(), bins.data(), systematics)));
	  // Response =  new MinervaUnfold::MnvResponse(Form("selected_mc_response2d_%s", name), name, bin_x, bin_y, bin_x, bin_y, systematics);
}

//////////Additional Stuff from Huma///////////////

 	 //=======================================================================================
 	   // When plotting, don't make Cuts on this Variable's axes
 	   //=======================================================================================
 	   //TODO: Move into PlotUtils?
 	   #define cut_map_t std::bitset<64>
 	   std::bitset<64> ignoreTheseCuts;
 	   //Populate ignoreTheseCuts.
 	   template <class CUTS_T> //Kill two bird templates with one stone template
 	   cut_map_t MatchCutsToVars(const CUTS_T& cuts)
 	   {
 	      const char* name = GetName().c_str();
 	     //Stupid compiler!  Let me use generic lambdas!
 	     using cut_t = decltype(*cuts.begin());
 	     ignoreTheseCuts.reset(); //Make sure I start off ignoring NO Cuts.
 	                              //Only cuts with 1 will be ignored.
 	     const auto found = std::find_if(cuts.begin(), cuts.end(),
 	                                     [this](const cut_t cut)
 	                                     {
 	                                       return cut->getName() == this->GetName();
 	                                     });
 	   
 	     if(found != cuts.end()) //If we found a match
 	     {
 	       ignoreTheseCuts.flip(std::distance(cuts.begin(), found));
 	     }
 	     else std::cerr << "Failed to find ANY Cut that matches " << GetName() << ".  You should double-check your Cut names.\n";
 	     return ignoreTheseCuts; //Just here in case it helps with debugging
 	   }
 	   //Use ignoreTheseCuts
 	   inline cut_map_t IgnoreMyVars(const cut_map_t allCuts) const
 	   {
 	     return allCuts | ignoreTheseCuts;
 	   }
 	 



//===================================================================================
void FillResponse1D(double x_reco, double x_true, const std::string name, double w, int unv){
	
	//std::cout<<name<<std::endl;
 	for(mnv_itr = Response1D.begin(); mnv_itr != Response1D.end(); ++mnv_itr){
		(mnv_itr->second)->Fill(x_reco,x_true,name,unv, w);
	}		
}




//=====================================

template <typename T>
void getResponseObjects1D(T univs)
{
//  bool status = false;
  for(mnv_itr2 = Response1D.begin(); mnv_itr2 != Response1D.end(); ++mnv_itr2){
                (mnv_itr2->second)->GetMigrationObjects( migrationH2D, h_reco1D, h_truth1D );
        }
  const bool clear_bands = true;  
  mresp1D = HW2D(migrationH2D, univs, clear_bands);
}
	
  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
 void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC) { m_selected_mc_reco.hist->Write();
               m_selected_mc_true.hist->Write();
               mresp1D.hist->Write();
 //for (auto responses:Response1D)responses.second->GetMigrationMatrix()->Write();   
}
    else m_selected_data_reco.hist->Write();

    // selected mc  histfolio fir Hist Stacking
   if(isMC) m_selected_mc_sb.WriteToFile(f);
   else m_selected_data_reco_sb.hist->Write();
  }
};

}  // namespace VarLoop



namespace Var2DLoop {
class Variable2D : public PlotUtils::Variable2DBase<NUKECC_ANA::CVUniverse> {
 private:
  //=======================================================================================
  // TYPEDEFS CONVENIENCE
  //=======================================================================================
  typedef PlotUtils::Hist2DWrapper<NUKECC_ANA::CVUniverse> HW2D;
  //typedef Histogram<NUKECC_ANA::CVUniverse> HW2D;
  typedef PlotUtils::MnvH2D MH2D;

 public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  Variable2D(ARGS... args) : Variable2DBase<NUKECC_ANA::CVUniverse>(args...) {}

  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // HISTWRAPPER
  HW2D m_selected_mc_reco,m_selected_data_reco, m_selected_Migration, mresp;
  //MinervaUnfold::MnvResponse* Response;
  std::map<std::string,MinervaUnfold::MnvResponse*> Response;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr2;

  MnvH2D *migrationH = NULL;
  MnvH2D *h_reco = NULL;
  MnvH2D *h_truth = NULL;
  //// HISTFOLIO
  // PlotUtils::HistFolio<MH2D> m_selected_mc_sb;

  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHistograms(T univs) {
    const bool clear_bands = true;  // we want empty histograms
    const char* name = GetName().c_str();

    // HISTWRAPPER
    // selected mc reco histwrapper
    MH2D* dummy_selected_mc_reco =
        new MH2D(Form("selected_mc_reco2d_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco = HW2D(dummy_selected_mc_reco, univs, clear_bands);
    
   
    MH2D* dummy_selected_data_reco =
        new MH2D(Form("selected_data2d_reco_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_data_reco = HW2D(dummy_selected_data_reco, univs, clear_bands);
    
    /*MH2D* dummy_selected_Migration =
        new MH2D(Form("selected_Migration_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsX(), GetBinVecX().data());

    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);*/
    /////For 2D analysis
    MH2D* dummy_selected_Migration =
        new MH2D(Form("selected_Migration_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);///For 2D Migration
   ///For 1D analysis
    /*MH2D* dummy_selected_Migration =
        new MH2D(Form("selected_Migration_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsX(), GetBinVecX().data());

    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);*/

    //delete Response; 
    delete dummy_selected_mc_reco;
    delete dummy_selected_Migration;
    delete dummy_selected_data_reco;
  }
//=====



void SetupResponse(std::map<const std::string, int> systematics){
//void SetupResponse(T univs){

	   const char* name = GetName().c_str();
	   axis_binning bin_x, bin_y;
	   bin_x.uniform=false;
          	
     	vector<double> vx;
	   //for(int i=0; i<GetNBinsX()+1; i++){vx.push_back(GetBinVecX().data()[i+1]-GetBinVecX().data()[i]);}
	   for(int i=0; i<=GetNBinsX(); i++){vx.push_back(GetBinVecX().data()[i]);}
     	   
	   vector<double> vy;
	   //for(int j=0; j<GetNBinsY()+1; j++){vy.push_back(GetBinVecY().data()[j+1]-GetBinVecY().data()[j]);}
	   for(int j=0; j<=GetNBinsY(); j++){vy.push_back(GetBinVecY().data()[j]);}
	   bin_x.bin_edges = vx;
	   bin_x.nbins	    = GetNBinsX();
	   bin_x.min 	    = GetBinVecX().data()[0];
	   bin_x.max       = GetBinVecX().data()[GetNBinsX()];  
	   bin_y.bin_edges = vy;
	   bin_y.nbins	    = GetNBinsY();
	   bin_y.min 	    = GetBinVecY().data()[0];
	   bin_y.max       = GetBinVecY().data()[GetNBinsY()]; 
	   Response.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("selected_mc_response2d_%s", name), name, bin_x, bin_y, bin_x, bin_y, systematics))); 
	  // Response =  new MinervaUnfold::MnvResponse(Form("selected_mc_response2d_%s", name), name, bin_x, bin_y, bin_x, bin_y, systematics);
}

//===================================================================================
//
//===================================================================================
void FillResponse(double x_reco, double y_reco, double x_true, double y_true, const std::string name, double w, int unv){
	
	//std::cout<<name<<std::endl;
 	for(mnv_itr = Response.begin(); mnv_itr != Response.end(); ++mnv_itr){
		(mnv_itr->second)->Fill(x_reco,y_reco,x_true,y_true,name,unv, w);
	}		
	
	
}
//===================================================================================
//
//===================================================================================
template <typename T>
void getResponseObjects(T univs)
{
 // bool status = false;
  for(mnv_itr2 = Response.begin(); mnv_itr2 != Response.end(); ++mnv_itr2){
                (mnv_itr2->second)->GetMigrationObjects( migrationH, h_reco, h_truth );;
        }
  const bool clear_bands = true;  
  mresp = HW2D(migrationH, univs, clear_bands);
}

//=======================================================================================
// WRITE ALL HISTOGRAMS
//=======================================================================================
void WriteAllHistogramsToFile(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ m_selected_mc_reco.hist->Write();
        m_selected_Migration.hist->Write();
// for (auto responses:Response)responses.second->GetMigrationMatrix()->Write();   
mresp.hist->Write(); 
   }//  Response->GetMigrationMatrix()->Write();}   
    else m_selected_data_reco.hist->Write();
    // selected mc reco
  }
};
}  // namespace Var2DLoop


#endif  // VARInVecXGetBinVecXGetBinVecXABLE_H
