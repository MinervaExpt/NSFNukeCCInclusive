#ifndef MNV_GLOBALINCLUDES_h
#define MNV_GLOBALINCLUDES_h 1

//-----------------------------------------------------------
// Included files needed from ROOT and C++
//-----------------------------------------------------------
// C++ headers
#include<iomanip>
#include<iostream>
#include<stdlib.h>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<string>
#include<sstream>
#include<time.h>
#include<zlib.h>
#include<algorithm>
#include<map>
#include<vector>

using std::cout;
using std::cin;
using std::endl;
using std::flush;
using std::setw;
using std::streamsize;
using std::setprecision;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ios;
using std::sort;
using std::map;
using std::pair;
using std::vector;

// ROOT headers
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TMath.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TPostScript.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TProfile.h"
#include "TVector3.h"
//#include "TStyle.h"
#include "TText.h"
//#include "TStyle.h"
#include "TRandom.h"
#include "TArrow.h"
#include "TLine.h"
#include "TMarker.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TVector.h"
#include "TArray.h"
#include "TLegend.h"
#include "THStack.h"
#include "TChainElement.h"
#include "TCollection.h"
#include "TFractionFitter.h"
#include "TParameter.h" 
#endif

#ifndef ROOT6
//#include "Cintex/Cintex.h"
#endif

//-----------------------------------------------------------
// Global definitions 
//-----------------------------------------------------------
namespace NUKECC_ANA{
  static const unsigned int nHistos = 54;
  static const unsigned int nHistosMy = 7;
  static const unsigned int nHistosPID = 5;
  static const unsigned int nHistosInc = 11;
  //static const bool RunCodeWithSystematics = true;// This variable will turn off and on systematics. This includes making all the lat/vert error bands. You CANNOT mix with/without processings
  //static const bool RunCodeWithSysFluxUniverses = false; // This variable will apply a non-cv flux to systematic universes (if available). This is done by an extra correction factor wgt given by fluxreweighter.
  //static const bool neutrinoMode = true; // Use: True = Neutrino Mode /  False = Antineutrino Mode	 
  static const bool neutrinoMode = false; // Use: True = Neutrino Mode /  False = Antineutrino Mode	 
  static const bool meMode = true; // Use: True = Medium Energy Mode / False =  Low Energy Mode 	 
  static const double bindingE = (neutrinoMode)? 34.0:30.0;

  static const int fsi_pdg = (neutrinoMode)? 2212: 2112; //2212: proton, 2112:neutron
  static const int fsi_config = 2; 
  static const bool fsi_rwQE = true;

  enum kAnalysis{//Analyses with differing signal definitions
    k2DNuLepton, //LE or ME neutrino mode 2D lepton results
    k2DAntiNuLepton, // LE or ME anti-neutrino mode 2D lepton result
    k2DTransverse, //Future LE or ME neutrino 2D transverse result
    kInclusive,//Future 2D inclusive result
    k3D,//Future 3D result. pt/pl/Eavail qelike
    k3DInclusive,//Future 3D result inclusive
    kNeutron//Future results with neutrons
  };


  /* Split MC in different templates and different ways
  * The different combinations are:
  *  kMC = kQE + kQENot
  *  kMC = kQELike + kQELikeNot
  *  kMC = kQE + kRES + kDIS + kOTH
  *  kMC = kQELike_QE + kQELike_RES + kQELike_OTH + kQELikeNot
  *  kMC = kQELike_QE + kQELike_RES + kQELike_OTH + kQELikeNot_PositivePions + kQELikeNot_NegativePions + kQELikeNot_NeutralPions + kQELikeNot_NoPions
  *  kMC = kQELike_QE + kQELike_RES + kQELike_OTH + kQELikeNot_PositivePions_TrueMichel + kQELikeNot_PositivePions_TrueMichelNot + kQELikeNot_NegativePions + kQELikeNot_NeutralPions + kQELikeNot_NoPions
  *  kMC = kQELike + kQELikeNot_PositivePions + kQELikeNot_NegativePions + kQELikeNot_NeutralPions + kQELikeNot_NoPions
  *  kMC = kQELike + kQELikeNot_PositivePions_TrueMichel + kQELikeNot_PositivePions_TrueMichelNot + kQELikeNot_NegativePions + kQELikeNot_NeutralPions + kQELikeNot_NoPions
  *  kMC = kQELike_QE + kQELike_QENot + kQELikeNot_QE + kQELikeNot_QENot
  *  kMC = kQELike + kQELikeNot_ChargedPions + kQELikeNot_NoChargedPions
  *  kMC = kQELike + kQELikeNot_ChargedPions_TrueMichel + kQELikeNot_NoChargedPions_TrueMichelNot + kQELikeNot_NoChargedPions
  */
  //kHistos is the driver of the CCQELike analysis
  enum kHistos {
    kData,
    kMC,
    kQE,
    kQE_H,
    kQE_OTH,
    kQENot,
    k2p2h,
    kRES,
    kDIS,
    kOTH,
    kQENot_PositivePions,
    kQENot_PositivePions_TrueMichel,
    kQENot_PositivePions_TrueMichelNot,
    kQENot_NegativePions,
    kQENot_NeutralPions,
    kQENot_NoPions,
    kQELike_QE,
    kQELike_QE_H,
    kQELike_QE_OTH,
    kQELike_2p2h,
    kQELike_RES,
    kQELike_DIS,
    kQELike_OTH,
    kQELike,
    kQELikeNot,
    kQELikeNot_NoPions,
    kQELikeNot_PositivePions,
    kQELikeNot_PositivePions_TrueMichel,
    kQELikeNot_PositivePions_TrueMichelNot,
    kQELikeNot_NegativePions,
    kQELikeNot_NeutralPions,
    kQELikeNot_ChargedPions,
    kQELikeNot_ChargedPions_TrueMichel,
    kQELikeNot_ChargedPions_TrueMichelNot,
    kQELikeNot_ChargedPionsNot,
    kQELike_QENot,
    kQELikeNot_QE,
    kQELikeNot_QENot,
    kQELikeNot_2p2h,
    kQELikeNot_SingleChargedPion,
    kQELikeNot_SingleNeutralPion,
    kQELikeNot_MultiPion,
    kQELikeNot_Other,
    kQELikeNot_DIS,
    kQELikeNot_RES,
    kQELikeNot_Coh,
    kQELike_RES_NeutronFSI,
    kQELike_RES_ProtonFSI,
    kQELike_RES_NeutronIS,
    kQELike_RES_ProtonIS,
    kQELike_QE_ProtonFSI,
    kQELike_QE_NeutronFSI,
    kQELike_2p2h_np,
    kQELike_2p2h_nn,
  };

  enum kHistosInc {//Inclusive mode of CCQENu
    kIncData,
    kIncMC,
    kIncCC,
    kIncNC,
    kIncQE,
    kIncRES,
    kIncDIS,
    kInc2p2h,
    kIncDIS_SIS,
    kIncDIS_DIS,
    kIncOth
  };

  enum kHistosMy {//Inclusive mode of CCQENu
      kIncDataMy,
        
       kIncMCMy,
      kIncDISMy,
 //    kIncCCMy,
  //  kIncQEMy,
      kInclowWMy,
     kInclowQ2My,
     kInclowQ2transMy,
     kInctransMy
     
  };

  /* Split MC depending on the reconstructed 
  *primary proton track true PID information
  *  pMC = pProton + pPion + pOther
  */
  enum pHistos { //proton score histos
    pData,
    pMC,
    pProton,
    pPion,
    pOther
  };

  static string names[nHistos] = {
    "data",
    "mc",
    "qe",
    "qe_h",
    "qe_oth",
    "qenot",
    "2p2h",
    "res",
    "dis",
    "oth",
    "qenot_positivepions",
    "qenot_positivepions_truemichel",
    "qenot_positivepions_truemichelnot",
    "qenot_negativepions",
    "qenot_neutralpions",
    "qenot_nopions",
    "qelike_qe",
    "qelike_qe_h",
    "qelike_qe_oth",
    "qelike_2p2h",
    "qelike_res",
    "qelike_dis",
    "qelike_oth",
    "qelike",
    "qelikenot",
    "qelikenot_nopions",
    "qelikenot_positivepions",
    "qelikenot_positivepions_truemichel",
    "qelikenot_positivepions_truemichelnot",
    "qelikenot_negativepions",
    "qelikenot_neutralpions",
    "qelikenot_chargedpions",
    "qelikenot_chargedpions_truemichel",
    "qelikenot_chargedpions_truemichelnot",
    "qelikenot_chargedpionsnot",
    "qelike_qenot",
    "qelikenot_qe",
    "qelikenot_qenot",
    "qelikenot_2p2h",
    "qelikenot_singlechargedpion",
    "qelikenot_singleneutralpion",
    "qelikenot_multipion",
    "qelikenot_not_scp_snp_mp",
    "qelikenot_dis",
    "qelikenot_res",
    "qelikenot_coh",
    "qelike_res_neutron_fsi",
    "qelike_res_proton_fsi",
    "qelike_res_neutron_is",
    "qelike_res_proton_is",
    "qelike_qe_proton_fsi",
    "qelike_qe_neutron_fsi",
    "qelike_2p2h_np",
    "qelike_2p2h_nn",

  };

  static string namesPID[nHistosPID] = {
    "data",
    "mc",
    "proton",
    "pion",
    "oth"
  };

  static string namesMy[nHistosMy] = {
    "data",
     
    "mc",
    "dis",
 //  "CC",
  //  "qe",
    "lowW",
    "lowQ2",
   "lowQ2Trans",
    "trans" 
  };



  static string namesInc[nHistosInc] = {
    "data",
    "mc",
    "CC",
    "NC",
    "qe",
    "res",
    "dis",
    "2p2h",
    "sis",
    "true_dis",
    "oth"
  };
}

namespace Playlist{ 
  const std::string minerva1 = "minerva1"; 
  const std::string minerva5 = "minerva5"; 
  const std::string minerva7 = "minerva7"; 
  const std::string minerva9 = "minerva9"; 
  const std::string minerva10 = "minerva10"; 
  const std::string minerva13A = "minerva13A"; 
  const std::string minerva13B = "minerva13B"; 
  const std::string minerva13C = "minerva13C"; 
  const std::string minerva13D = "minerva13D"; 
  const std::string minerva13E = "minerva13E"; 
  const std::string minerva13  = "minerva13"; 
  const std::string minervame1A  = "minervame1A"; 
  const std::string minervame1B  = "minervame1B"; 
  const std::string minervame1C  = "minervame1C"; 
  const std::string minervame1D  = "minervame1D"; 
  const std::string minervame1E  = "minervame1E"; 
  const std::string minervame1F  = "minervame1F"; 
  const std::string minervame1G  = "minervame1G"; 
  const std::string minervame1H  = "minervame1H"; 
  const std::string minervame1L  = "minervame1L"; 
  const std::string minervame1M  = "minervame1M"; 
  const std::string minervame1N  = "minervame1N"; 
  const std::string minervame1O  = "minervame1O"; 
  const std::string minervame1P  = "minervame1P"; 
  const std::string minervame5A  = "minervame5A"; 
  const std::string minervame6A  = "minervame6A"; 
  const std::string minervame6B  = "minervame6B"; 
 
  const std::string all_nu = "Nu"; 
  const std::string all_antinu = "AntiNu"; 
}

namespace Analysis{ 
  const std::string NukeCC = "NukeCC"; 
}

const std::string DEFAULT_PLAYLIST = Playlist::all_nu; 


#endif
