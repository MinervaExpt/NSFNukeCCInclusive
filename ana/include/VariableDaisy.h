#ifndef VARIABLE_H
#define VARIABLE_H

#include <iterator>
#include "../../NUKECCSRC/include/CVUniverse.h"
#include "PlotUtils/HistFolio.h"
#include "PlotUtils/HistWrapper.h"

#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MnvH2D.h"
#ifndef __CINT__  // CINT doesn't know about std::function
#include "PlotUtils/VariableBase.h"
#include "PlotUtils/Variable2DBase.h"
#include "MinervaUnfold/MnvResponse.h"
#include <bitset>
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
  HW m_selected_mc_reco, m_selected_data_reco, m_selected_truth_reco, m_selected_data_reco_sb, m_selected_mc_true;

  typedef PlotUtils::Hist2DWrapper<NUKECC_ANA::CVUniverse> HW2D;
  HW2D mresp1D;
  typedef PlotUtils::MnvH2D MH2D;
  HW2D m_selected_Migration;
  std::map<std::string,MinervaUnfold::MnvResponse*> Response1D;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr2; 

  MnvH2D *migrationH2D = NULL;
  MnvH1D *h_reco1D = NULL;
  MnvH1D *h_truth1D = NULL;

  HW daisy_petal_num_hists[12];
  HW daisy_petal_denom_hists[12];
  HW daisy_petal_denom_hists_QE[12];
  HW daisy_petal_denom_hists_RES[12];
  HW daisy_petal_denom_hists_DIS[12];
  HW daisy_petal_denom_hists_Other[12];
  HW daisy_petal_denom_hists_2p2h[12];
  HW2D daisy_petal_migration_hists[12];
 // HISTFOLIO
  // selected mc reco - signal background histfolio
  PlotUtils::HistFolio<MH1D> m_selected_mc_sb;
  PlotUtils::HistFolio<MH1D> daisy_petal_denom_hists_sb[12];
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
    MH1D* dummy_selected_mc_reco = new MH1D(Form("h_mc_%s", name), name,GetNBins(), bins.data());
    m_selected_mc_reco = HW(dummy_selected_mc_reco, univs, clear_bands);

    //for 1D analysis
    MH1D* dummy_selected_mc_true = new MH1D(Form("h_mc_true_%s", name), name,GetNBins(), bins.data());
    m_selected_mc_true = HW(dummy_selected_mc_true, univs, clear_bands);
    
    //////////For Efficiency 
    MH1D* dummy_selected_truth_reco = new MH1D(Form("h_truth_%s", name), name, GetNBins(), bins.data());
    m_selected_truth_reco = HW(dummy_selected_truth_reco, univs, clear_bands);

    // To extract in daisy bins
    for(int petal=0; petal<12; petal++){
      MH1D* dummy_selected_mc_reco_daisy = new MH1D(Form("h_mc_daisy_%d_%s", petal, name), name, GetNBins(), bins.data());
      daisy_petal_num_hists[petal] = HW(dummy_selected_mc_reco_daisy, univs, clear_bands);

      MH1D* dummy_selected_truth_reco_daisy = new MH1D(Form("h_truth_daisy_%d_%s", petal, name), name, GetNBins(), bins.data());
      daisy_petal_denom_hists[petal] = HW(dummy_selected_truth_reco_daisy, univs, clear_bands);

      MH1D* dummy_selected_truth_reco_daisy_QE = new MH1D(Form("h_truth_daisy_QE_%d_%s", petal, name), name, GetNBins(), bins.data());
      daisy_petal_denom_hists_QE[petal] = HW(dummy_selected_truth_reco_daisy_QE, univs, clear_bands);
      
      MH1D* dummy_selected_truth_reco_daisy_RES = new MH1D(Form("h_truth_daisy_RES_%d_%s", petal, name), name, GetNBins(), bins.data());
      daisy_petal_denom_hists_RES[petal] = HW(dummy_selected_truth_reco_daisy_RES, univs, clear_bands);

      MH1D* dummy_selected_truth_reco_daisy_DIS = new MH1D(Form("h_truth_daisy_DIS_%d_%s", petal, name), name, GetNBins(), bins.data());
      daisy_petal_denom_hists_DIS[petal] = HW(dummy_selected_truth_reco_daisy_DIS, univs, clear_bands);

      MH1D* dummy_selected_truth_reco_daisy_Other = new MH1D(Form("h_truth_daisy_Other_%d_%s", petal, name), name, GetNBins(), bins.data());
      daisy_petal_denom_hists_Other[petal] = HW(dummy_selected_truth_reco_daisy_Other, univs, clear_bands);

      MH1D* dummy_selected_truth_reco_daisy_2p2h = new MH1D(Form("h_truth_daisy_2p2h_%d_%s", petal, name), name, GetNBins(), bins.data());
      daisy_petal_denom_hists_2p2h[petal] = HW(dummy_selected_truth_reco_daisy_2p2h, univs, clear_bands);

      daisy_petal_denom_hists_sb[petal] = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("h_truth_daisy_sb_%d_%s", petal, name), name, GetNBins(), bins.data());

      //daisy_petal_denom_hists_sb[petal].AddComponentHist("QE");
      //daisy_petal_denom_hists_sb[petal].AddComponentHist("RES");
      //daisy_petal_denom_hists_sb[petal].AddComponentHist("DIS");
      //daisy_petal_denom_hists_sb[petal].AddComponentHist("Other");
      //daisy_petal_denom_hists_sb[petal].AddComponentHist("2p2h");

      delete dummy_selected_mc_reco_daisy;
      delete dummy_selected_truth_reco_daisy;
      delete dummy_selected_truth_reco_daisy_QE;
      delete dummy_selected_truth_reco_daisy_RES;
      delete dummy_selected_truth_reco_daisy_DIS;
      delete dummy_selected_truth_reco_daisy_Other;
      delete dummy_selected_truth_reco_daisy_2p2h;
    }

    //For Data
    MH1D* dummy_selected_data_reco = new MH1D(Form("h_data_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco = HW(dummy_selected_data_reco, univs, clear_bands);


    //For Data in sidebands
    MH1D* selected_data_reco_sb = new MH1D(Form("h_data_sb_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco_sb = HW(selected_data_reco_sb, univs, clear_bands);
  
    //HISTFOLIO
    //selected mc reco - signal background histfolio
    m_selected_mc_sb = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("selected_mc_sb_%s", name), name, GetNBins(), bins.data());
    
 //PlotUtils::MnvH1D* data = new PlotUtils::MnvH1D(
 //   "dummy", "dummy", plotting::nbins, plotting::xmin, plotting::xmax);
 //  m_selected_data_sb = PlotUtils::HistFolio<PlotUtils::MnvH1D>(
 //    Form("selected_data_sb_%s", name), name, GetNBins(), bins.data());
   
    m_selected_mc_sb.AddComponentHist("DIS");
    m_selected_mc_sb.AddComponentHist("MC");
   // m_selected_data_sb.AddComponentHist("Data");

    MH2D* dummy_selected_Migration = new MH2D(Form("selected_Migration_%s", name), name, GetNBins(),GetBinVec().data(), GetNBins(), GetBinVec().data());
    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);

    // To extract in daisy bins
    for(int petal=0; petal<12; petal++){
      MH2D* dummy_selected_Migration_daisy = new MH2D(Form("selected_Migration_daisy_%d_%s", petal, name), name, GetNBins(),GetBinVec().data(), GetNBins(), GetBinVec().data());
      daisy_petal_migration_hists[petal] = HW2D(dummy_selected_Migration_daisy, univs, clear_bands);

      delete dummy_selected_Migration_daisy;
    }


delete dummy_selected_mc_reco;
delete dummy_selected_truth_reco;
delete dummy_selected_Migration;
delete dummy_selected_mc_true;
delete dummy_selected_data_reco;
  }


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


void SetupResponse1D(std::map<const std::string, int> systematics){
//void SetupResponse(T univs){
	   std::vector<double> bins = GetBinVec();
	   const char* name = GetName().c_str();

	   axis_binning bin_x;
	   bin_x.uniform=false;
	
     vector<double> vx;
     	   
	   for(int i=0; i<=GetNBins(); i++){vx.push_back(GetBinVec().data()[i]);}
     	   
	   bin_x.bin_edges = vx;
	   bin_x.nbins	= GetNBins();
	   bin_x.min 	  = GetBinVec().data()[0];
	   bin_x.max    = GetBinVec().data()[GetNBins()];
           
      cout<<"bins min = "<<bin_x.min<<"bins max = "<<bin_x.max<<"Number of bins = "<<bin_x.nbins<<endl;
            
	   //Response1D.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("selected_mc_response1d_%s", name), name, bin_x, bin_x, systematics))); 
	   Response1D.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("response1d_%s", name), name, bin_x, bin_x, systematics))); 
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
  // WRITE ALL HISTOGRAMS FOR EVENTLOOP
  //=======================================================================================
 void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC) { m_selected_mc_reco.hist->Write();
}
    else {m_selected_data_reco.hist->Write();
}
    // selected mc  histfolio fir Hist Stacking
   //if(isMC) m_selected_mc_sb.WriteToFile(f);
   //else m_selected_data_reco_sb.hist->Write();
  }
  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
 void WriteAllHistogramsToFileEff(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC) { 
      m_selected_mc_reco.hist->Write();
      for(int petal=0; petal<12; petal++){
        daisy_petal_num_hists[petal].hist->Write();
      }
    }
    else {
      m_selected_truth_reco.hist->Write();
      for(int petal=0; petal<12; petal++){
        daisy_petal_denom_hists[petal].hist->Write();
        daisy_petal_denom_hists_QE[petal].hist->Write();
        daisy_petal_denom_hists_RES[petal].hist->Write();
        daisy_petal_denom_hists_DIS[petal].hist->Write();
        daisy_petal_denom_hists_Other[petal].hist->Write();
        daisy_petal_denom_hists_2p2h[petal].hist->Write();
        daisy_petal_denom_hists_sb[petal].WriteToFile(f);
      }
    }
}
  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
 void WriteAllHistogramsToFileMig(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC) { m_selected_mc_reco.hist->Write();
               m_selected_mc_true.hist->Write();
               m_selected_Migration.hist->Write();
               mresp1D.hist->Write();
    for(int petal=0; petal<12; petal++){
      daisy_petal_num_hists[petal].hist->Write();
      daisy_petal_migration_hists[petal].hist->Write();

    }
    
    }
    else {m_selected_data_reco.hist->Write();
}
    // selected mc  histfolio fir Hist Stacking
   //if(isMC) m_selected_mc_sb.WriteToFile(f);
   //else m_selected_data_reco_sb.hist->Write();
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
  HW2D m_selected_mc_reco,m_selected_data_reco,m_selected_truth_reco,m_selected_Migration, mresp;
  //MinervaUnfold::MnvResponse* Response;
  std::map<std::string,MinervaUnfold::MnvResponse*> Response;
  //std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr;
  //std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr2;

  MnvH2D *migrationH = NULL;
  MnvH2D *h_reco = NULL;
  MnvH2D *h_truth = NULL;
  //// HISTFOLIO
  // PlotUtils::HistFolio<MH2D> m_selected_mc_sb;

  HW2D daisy_petal_num2d_hists[12];
  HW2D daisy_petal_denom2d_hists[12];
  HW2D daisy_petal_denom2d_hists_QE[12];
  HW2D daisy_petal_denom2d_hists_RES[12];
  HW2D daisy_petal_denom2d_hists_DIS[12];
  HW2D daisy_petal_denom2d_hists_Other[12];
  HW2D daisy_petal_denom2d_hists_2p2h[12];

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
        new MH2D(Form("h_mc_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco = HW2D(dummy_selected_mc_reco, univs, clear_bands);
    
    
    MH2D* dummy_selected_truth_reco =
        new MH2D(Form("h_truth_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_truth_reco = HW2D(dummy_selected_truth_reco, univs, clear_bands);
  
  for(int petal=0; petal<12; petal++){
    MH2D* dummy_selected_mc2d_reco = new MH2D(Form("h_mc2d_%d_%s", petal, name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    daisy_petal_num2d_hists[petal] = HW2D(dummy_selected_mc2d_reco, univs, clear_bands);
    
    MH2D* dummy_selected_truth2d_reco = new MH2D(Form("h_truth2d_%d_%s", petal, name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    daisy_petal_denom2d_hists[petal] = HW2D(dummy_selected_truth2d_reco, univs, clear_bands);

    MH2D* dummy_selected_truth2d_reco_daisy_QE = new MH2D(Form("h_truth2d_daisy_QE_%d_%s", petal, name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    daisy_petal_denom2d_hists_QE[petal] = HW2D(dummy_selected_truth2d_reco_daisy_QE, univs, clear_bands);
    
    MH2D* dummy_selected_truth2d_reco_daisy_RES = new MH2D(Form("h_truth2d_daisy_RES_%d_%s", petal, name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    daisy_petal_denom2d_hists_RES[petal] = HW2D(dummy_selected_truth2d_reco_daisy_RES, univs, clear_bands);

    MH2D* dummy_selected_truth2d_reco_daisy_DIS = new MH2D(Form("h_truth2d_daisy_DIS_%d_%s", petal, name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    daisy_petal_denom2d_hists_DIS[petal] = HW2D(dummy_selected_truth2d_reco_daisy_DIS, univs, clear_bands);

    MH2D* dummy_selected_truth2d_reco_daisy_Other = new MH2D(Form("h_truth2d_daisy_Other_%d_%s", petal, name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    daisy_petal_denom2d_hists_Other[petal] = HW2D(dummy_selected_truth2d_reco_daisy_Other, univs, clear_bands);

    MH2D* dummy_selected_truth2d_reco_daisy_2p2h = new MH2D(Form("h_truth2d_daisy_2p2h_%d_%s", petal, name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    daisy_petal_denom2d_hists_2p2h[petal] = HW2D(dummy_selected_truth2d_reco_daisy_2p2h, univs, clear_bands);

    delete dummy_selected_mc2d_reco;
    delete dummy_selected_truth2d_reco;
    delete dummy_selected_truth2d_reco_daisy_QE;
    delete dummy_selected_truth2d_reco_daisy_RES;
    delete dummy_selected_truth2d_reco_daisy_DIS;
    delete dummy_selected_truth2d_reco_daisy_Other;
    delete dummy_selected_truth2d_reco_daisy_2p2h;

  }
    
   // MH2D* dummy_selected_truth_reco =
   //     new MH2D(Form("selected_truth2d_reco_%s", name), name, setBinLogX(), GetBinVecX().data(), setBinLogY(), GetBinVecY().data());
   //              GetBinVecX().data(), setBinLogY(), GetBinVecY().data());
   // m_selected_truth_reco = HW2D(dummy_selected_truth_reco, univs, clear_bands);

    MH2D* dummy_selected_data_reco =
        new MH2D(Form("h_data_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_data_reco = HW2D(dummy_selected_data_reco, univs, clear_bands);
    

    /*MH2D* dummy_selected_Migration =
        new MH2D(Form("selected_Migration_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsX(), GetBinVecX().data());
    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);*/
    

    /////For 2D analysis
    MH2D* dummy_selected_Migration = new MH2D(Form("Migration_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsX(), GetBinVecX().data());
    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);///For 2D Migration
   

    ///For 1D analysis
    /*MH2D* dummy_selected_Migration =
        new MH2D(Form("selected_Migration_%s", name), name, GetNBinsX(),GetBinVecX().data(), GetNBinsX(), GetBinVecX().data());
    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);*/

    //delete Response; 
    delete dummy_selected_mc_reco;
    delete dummy_selected_truth_reco;
    delete dummy_selected_Migration;
    delete dummy_selected_data_reco;
  }
//=====



MinervaUnfold::MnvResponse* SetupResponse(MinervaUnfold::MnvResponse* Response, int petal, std::map<const std::string, int> systematics){
//void SetupResponse(T univs){

  const char* name = GetName().c_str();
  axis_binning bin_x, bin_y;
  bin_x.uniform=false;

      vector<double> vx;
  
  for(int i=0; i<=GetNBinsX(); i++){vx.push_back(GetBinVecX().data()[i]);}
      
  vector<double> vy;
  for(int j=0; j<=GetNBinsY(); j++){vy.push_back(GetBinVecY().data()[j]);}
  bin_x.bin_edges = vx;
  bin_x.nbins	    = GetNBinsX();
  bin_x.min 	    = GetBinVecX().data()[0];
  bin_x.max       = GetBinVecX().data()[GetNBinsX()];  
  bin_y.bin_edges = vy;
  bin_y.nbins	    = GetNBinsY();
  bin_y.min 	    = GetBinVecY().data()[0];
  bin_y.max       = GetBinVecY().data()[GetNBinsY()]; 

  //Response.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("selected_mc_response2d_%s", name), name, bin_x, bin_y, bin_x, bin_y, systematics))); 
  Response = new MinervaUnfold::MnvResponse(Form("response2d_daisy_%d_%s", petal, name), name, bin_x, bin_y, bin_x, bin_y, systematics);

  return Response;
}


//===================================================================================
//
//===================================================================================
void FillResponse(MinervaUnfold::MnvResponse* Response, double x_reco, double y_reco, double x_true, double y_true, const std::string name, double w, int unv){
  Response->Fill(x_reco,y_reco,x_true,y_true,name,unv, w);
}		

//===================================================================================
//
//===================================================================================
template <typename T>
void getResponseObjects(MinervaUnfold::MnvResponse* Response, T univs)
{
   
  bool check = Response->GetMigrationObjects( migrationH, h_reco, h_truth );
  const bool clear_bands = false;  
  mresp = HW2D(migrationH, univs, clear_bands);
}

//=======================================================================================
// WRITE ALL HISTOGRAMS
//=======================================================================================
void WriteAllHistogramsToFile(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ m_selected_mc_reco.hist->Write();
   }//  Response->GetMigrationMatrix()->Write();}   
    else { m_selected_data_reco.hist->Write();
     }
    // selected mc reco
  }
//=======================================================================================
// WRITE ALL HISTOGRAMS
//=======================================================================================
void WriteAllHistogramsToFileEff(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ 
         m_selected_mc_reco.hist->Write();
         for(int petal=0; petal<12; petal++){
          daisy_petal_num2d_hists[petal].hist->Write();
          }
        }  
    else { 
      m_selected_truth_reco.hist->Write();
      for(int petal=0; petal<12; petal++){
        daisy_petal_denom2d_hists[petal].hist->Write();
        daisy_petal_denom2d_hists_QE[petal].hist->Write();
        daisy_petal_denom2d_hists_RES[petal].hist->Write();
        daisy_petal_denom2d_hists_DIS[petal].hist->Write();
        daisy_petal_denom2d_hists_Other[petal].hist->Write();
        daisy_petal_denom2d_hists_2p2h[petal].hist->Write();  
    }
     }
    // selected mc reco
  }
//=======================================================================================
// WRITE ALL HISTOGRAMS
//=======================================================================================
void WriteAllHistogramsToFileMig(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ m_selected_mc_reco.hist->Write();
                 m_selected_Migration.hist->Write();
      //for (auto responses:Response)responses.second->GetMigrationMatrix()->Write();   
        mresp.hist->Write(); 
   }//  Response->GetMigrationMatrix()->Write();}   
    else {m_selected_data_reco.hist->Write();
     }
  }
};
}  // namespace Var2DLoop


#endif  // VARIABLE_H