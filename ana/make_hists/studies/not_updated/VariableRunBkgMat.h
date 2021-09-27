#ifndef VARIABLE_H
#define VARIABLE_H

#include "include/CVUniverse.h"
#include "PlotUtils/HistFolio.h"
#include "PlotUtils/HistWrapper.h"

#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MnvH2D.h"
#ifndef __CINT__  // CINT doesn't know about std::function
#include "PlotUtils/VariableBase.h"

#include "PlotUtils/Variable2DBase.h"
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
  HW m_selected_mc_reco,m_selected_data_reco,m_selected_data_reco_sb;
  HW m_selected_mc_reco_signal, m_selected_mc_reco_bkg, m_selected_mc_reco_bkg_NC, m_selected_mc_reco_bkg_WrongSign, m_selected_mc_reco_bkg_muon;
  HW m_selected_mc_reco_bkg_fidcut, m_selected_mc_reco_bkg_water;

  // HISTFOLIO
  // selected mc reco - signal background histfolio
  PlotUtils::HistFolio<MH1D> m_selected_mc_sb;

  // histfolio for different plastic sources
  PlotUtils::HistFolio<MH1D> m_selected_mc_plastic;
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
    // all reco
    MH1D* dummy_selected_mc_reco = new MH1D(Form("selected_mc_reco_%s", name), name,
                                      GetNBins(), bins.data());
    m_selected_mc_reco = HW(dummy_selected_mc_reco, univs, clear_bands);

    // signal
    MH1D* dummy_selected_mc_reco_signal = new MH1D(Form("selected_mc_reco_signal_%s", name), name,
                                      GetNBins(), bins.data());
    m_selected_mc_reco_signal = HW(dummy_selected_mc_reco_signal, univs, clear_bands);

    // all bkg
    MH1D* dummy_selected_mc_reco_bkg = new MH1D(Form("selected_mc_reco_bkg_%s", name), name,
                                      GetNBins(), bins.data());
    m_selected_mc_reco_bkg = HW(dummy_selected_mc_reco_bkg, univs, clear_bands);

    // bkg NC
    MH1D* dummy_selected_mc_reco_bkg_NC = new MH1D(Form("selected_mc_reco_bkg_NC_%s", name), name,
                                      GetNBins(), bins.data());
    m_selected_mc_reco_bkg_NC = HW(dummy_selected_mc_reco_bkg_NC, univs, clear_bands);

    // bkg wrong sign
    MH1D* dummy_selected_mc_reco_bkg_WrongSign = new MH1D(Form("selected_mc_reco_bkg_WrongSign_%s", name), name,
                                      GetNBins(), bins.data());
    m_selected_mc_reco_bkg_WrongSign = HW(dummy_selected_mc_reco_bkg_WrongSign, univs, clear_bands);

    // bkg fid cut
    MH1D* dummy_selected_mc_reco_bkg_fidcut= new MH1D(Form("selected_mc_reco_bkg_fidcut_%s", name), name,
                                      GetNBins(), bins.data());
    m_selected_mc_reco_bkg_fidcut= HW(dummy_selected_mc_reco_bkg_fidcut, univs, clear_bands);

    // bkg muon kinematics
    MH1D* dummy_selected_mc_reco_bkg_muon = new MH1D(Form("selected_mc_reco_bkg_muon_%s", name), name,
                                      GetNBins(), bins.data());
    m_selected_mc_reco_bkg_muon = HW(dummy_selected_mc_reco_bkg_muon, univs, clear_bands);

    // bkg water
    MH1D* dummy_selected_mc_reco_bkg_water = new MH1D(Form("selected_mc_reco_bkg_water_%s", name), name,
                                      GetNBins(), bins.data());
    m_selected_mc_reco_bkg_water = HW(dummy_selected_mc_reco_bkg_water, univs, clear_bands);


    // DATA
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

    // material background type
    // all material bkg
    m_selected_mc_sb.AddComponentHist("MatBkg");
    m_selected_mc_sb.AddComponentHist("FeT3");
    m_selected_mc_sb.AddComponentHist("PbT3");
    m_selected_mc_sb.AddComponentHist("CT3");
   
    m_selected_mc_sb.AddComponentHist("otherFe");
    m_selected_mc_sb.AddComponentHist("otherPb");
    m_selected_mc_sb.AddComponentHist("otherC");
    m_selected_mc_sb.AddComponentHist("Plastic");
    // tops 8 subcategories

    // differrent sources of plastic
    m_selected_mc_plastic = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("selected_mc_plastic_%s", name), name, GetNBins(), bins.data());
    m_selected_mc_plastic.AddComponentHist("H");
    m_selected_mc_plastic.AddComponentHist("C");
    m_selected_mc_plastic.AddComponentHist("O");
    m_selected_mc_plastic.AddComponentHist("Al");
    m_selected_mc_plastic.AddComponentHist("Si");
    m_selected_mc_plastic.AddComponentHist("Ti");
    m_selected_mc_plastic.AddComponentHist("He");
    m_selected_mc_plastic.AddComponentHist("N");

delete dummy_selected_mc_reco;
delete dummy_selected_mc_reco_signal;
delete dummy_selected_mc_reco_bkg;
delete dummy_selected_mc_reco_bkg_NC;
delete dummy_selected_mc_reco_bkg_WrongSign;
delete dummy_selected_mc_reco_bkg_fidcut;
delete dummy_selected_mc_reco_bkg_muon;
delete dummy_selected_mc_reco_bkg_water;
delete dummy_selected_data_reco;
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
 void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC) {
      m_selected_mc_reco.hist->Write();
      m_selected_mc_reco_signal.hist->Write();
      m_selected_mc_reco_bkg.hist->Write();
      m_selected_mc_reco_bkg_NC.hist->Write();
      m_selected_mc_reco_bkg_WrongSign.hist->Write();
      m_selected_mc_reco_bkg_fidcut.hist->Write();
      m_selected_mc_reco_bkg_muon.hist->Write();
      m_selected_mc_reco_bkg_water.hist->Write();
    }

    else m_selected_data_reco.hist->Write();

    // selected mc  histfolio fir Hist Stacking
   if(isMC) {
    m_selected_mc_sb.WriteToFile(f);
    m_selected_mc_plastic.WriteToFile(f);
   }
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
  HW2D m_selected_mc_reco, m_selected_mc_reco_Sg, m_selected_mc_reco_Fe, m_selected_mc_reco_Pb, m_selected_mc_reco_C, m_selected_mc_reco_Plastic, m_selected_data_reco;
  HW2D m_selected_mc_reco_other, m_selected_mc_reco_otherFe, m_selected_mc_reco_otherPb, m_selected_mc_reco_otherC;

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
    // all reco
    MH2D* dummy_selected_mc_reco =
        new MH2D(Form("selected_mc_reco2d_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_mc_reco = HW2D(dummy_selected_mc_reco, univs, clear_bands);

    // signal
    MH2D* dummy_selected_mc_reco_Sg =
        new MH2D(Form("selected_mc_reco2d_Sg_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_mc_reco_Sg = HW2D(dummy_selected_mc_reco_Sg, univs, clear_bands);

    // target 3, different material
    // Fe in T3
    MH2D* dummy_selected_mc_reco_Fe =
        new MH2D(Form("selected_mc_reco2d_Fe_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_mc_reco_Fe = HW2D(dummy_selected_mc_reco_Fe, univs, clear_bands);

    // Pb in T3
    MH2D* dummy_selected_mc_reco_Pb =
        new MH2D(Form("selected_mc_reco2d_Pb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_mc_reco_Pb = HW2D(dummy_selected_mc_reco_Pb, univs, clear_bands);

    // C in T3
    MH2D* dummy_selected_mc_reco_C =
        new MH2D(Form("selected_mc_reco2d_C_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_mc_reco_C = HW2D(dummy_selected_mc_reco_C, univs, clear_bands);

    // Other targets + plastic
    // other ALL
    MH2D* dummy_selected_mc_reco_other =
        new MH2D(Form("selected_mc_reco2d_other_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_mc_reco_other = HW2D(dummy_selected_mc_reco_other, univs, clear_bands);

     // other target Fe
    MH2D* dummy_selected_mc_reco_otherFe =
        new MH2D(Form("selected_mc_reco2d_otherFe_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_mc_reco_otherFe = HW2D(dummy_selected_mc_reco_otherFe, univs, clear_bands);

    // other target Pb
    MH2D* dummy_selected_mc_reco_otherPb =
        new MH2D(Form("selected_mc_reco2d_otherPb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_mc_reco_otherPb = HW2D(dummy_selected_mc_reco_otherPb, univs, clear_bands);

    // otehr target C
    MH2D* dummy_selected_mc_reco_otherC =
        new MH2D(Form("selected_mc_reco2d_otherC_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_mc_reco_otherC = HW2D(dummy_selected_mc_reco_otherC, univs, clear_bands);

    // plastic
    MH2D* dummy_selected_mc_reco_Plastic =
        new MH2D(Form("selected_mc_reco2d_Plastic_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());

    m_selected_mc_reco_Plastic = HW2D(dummy_selected_mc_reco_Plastic, univs, clear_bands);
    
    
    MH2D* dummy_selected_data_reco =
        new MH2D(Form("selected_data2d_reco_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_data_reco = HW2D(dummy_selected_data_reco, univs, clear_bands);

    delete dummy_selected_mc_reco;
    delete dummy_selected_mc_reco_Sg;
    delete dummy_selected_mc_reco_Fe;
    delete dummy_selected_mc_reco_Pb;
    delete dummy_selected_mc_reco_C;
    delete dummy_selected_mc_reco_other;
    delete dummy_selected_mc_reco_otherFe;
    delete dummy_selected_mc_reco_otherPb;
    delete dummy_selected_mc_reco_otherC;
    delete dummy_selected_mc_reco_Plastic;

    delete dummy_selected_data_reco;
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f,bool isMC) const {
    f.cd();
       if(isMC) {m_selected_mc_reco.hist->Write();
       m_selected_mc_reco_Sg.hist->Write();
       m_selected_mc_reco_Fe.hist->Write();
       m_selected_mc_reco_Pb.hist->Write();
       m_selected_mc_reco_C.hist->Write();
       m_selected_mc_reco_other.hist->Write();
       m_selected_mc_reco_otherFe.hist->Write();
       m_selected_mc_reco_otherPb.hist->Write();
       m_selected_mc_reco_otherC.hist->Write();
       m_selected_mc_reco_Plastic.hist->Write();
       }
       else m_selected_data_reco.hist->Write();
    // selected mc reco
  }
};
}  // namespace Var2DLoop


#endif  // VARIABLE_H
