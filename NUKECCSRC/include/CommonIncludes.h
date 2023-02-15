#ifndef _CommonIncludes_H
#define _CommonIncludes_H 

//From PlotUtils
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>
#include <PlotUtils/MnvH3D.h>
//#include <PlotUtils/MnvRecoShifter.h>
#include <PlotUtils/MnvPlotter.h>
#include <PlotUtils/HistogramUtils.h>
#include <PlotUtils/MnvApplication.h>

//note - don't include namespaces or headers with function definitions, or risk linker failures
//#include <PlotUtils/MnvNormalization.h>
//#include <PlotUtils/MnvAnaTuple.h>
//#include <PlotUtils/ArachneUtils.h>

//#include <Acceptance/TAcceptanceTable.h>
//#include <MinervaUnfold/MnvUnfold.h>

//yeah... I want to be in these namespaces for everything
using namespace std;
using namespace PlotUtils;

//From ROOT
#include <Rtypes.h>
#include <TSystem.h>
#include <TError.h>
#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TCut.h>
#include <TText.h>
#include <TLegendEntry.h>
#include <TFitResult.h>
#include <TF1.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TH2D.h>
#include <TH1F.h>
#include <TVector2.h>
#include <TNamed.h>
#include <THStack.h>
#include <TPie.h>
#include <TGraphAsymmErrors.h>
//#include <TStyle.h>
#include <TMath.h>
#include <TLatex.h>
//#include <TObject.h>
//c++
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//Use file or not...
//const bool USE_NUKEONLYMC =  false; //false; //


namespace NUKECC_ANA
{

 
  //static const bool RunCodeWithSystematics = false;// This variable will turn off and on systematics. This includes making all the lat/vert error bands. You CANNOT mix with/without processings
  static const bool RunCodeWithSystematics = true;// This variable will turn off and on systematics. This includes making all the lat/vert error bands. You CANNOT mix with/without processings
  static const bool RunCodeWithSysFluxUniverses = true; // This variable will apply a non-cv flux to systematic universes (if available). This is done by an extra correction factor wgt given by fluxreweighter.
  //static const bool neutrinoMode = false; // Use: True = Neutrino Mode /  False = Antineutrino Mode	 
  //vertexing
  //const bool USE_ML_VERTEX = true;

  //typedefs
  typedef std::pair<int,int> MatPair;
  typedef std::vector< MatPair > MatPairs;
  typedef std::pair< std::string, const double* > ErrorVar;
  typedef std::vector<ErrorVar> ErrorVars;

  //constants
  const double PI = atan(1.0)*4;
  const double deg_to_rad = PI / 180.0;
  const double rad_to_deg = 1. / deg_to_rad;
  const double mev_to_gev = .001;
  const double avagadro   = 6.0221415E23;
  const double EPSILON    = 1e-5;

  const bool NO_SPREAD_ERROR = true;
  const unsigned int N_VERT_UNIVERSES = 100;
  const unsigned int N_LAT_UNIVERSES  = 60;

  //reco range is wider because kinematic signal may be smeared out of the accepted range
  //const double MIN_E = 1.6*1000; //in MeV
  //const double MAX_E = 25.*1000; //in MeV
  // ---- no longer used, cut on Enu instead -----

  //for absolute cross sections of Enu, use 5-50 GeV
  const double MIN_RECO_E = 0.*1000; //in MeV // cut out events with neutrino energy smaller than this
  const double MAX_RECO_E = 120.*1000; //in MeV // cut out events with neutrino energy larger than this

  const double MIN_RECO_E_MU = 2.*1000; //in MeV // cut out events with muon energy smaller than this
  const double MAX_RECO_E_MU = 20.*1000; //in MeV // cut out events with muon energy larger than this

  // Inelastic Cut
  const bool USE_INEL_CUT = false;           // should I use the inelastic cut to get a dis-like sample?
  const bool TRIANGLE_INELCUT = true;        // switch between two diff 2D-exclusion cuts (triangle vs step)
  const double MIN_CCQE_RECOIL_E = 1000.0;         // MeV // default for flat cut on ccqe recoil
  const double LOW_Q2_CCQE_RECOIL_E_CUT  = 1000.0; // MeV // events with Q2<1 must have more ccqerecoil than this
  const double HIGH_Q2_CCQE_RECOIL_E_CUT = 500.0;  // MeV // events with 1<Q2 must have more ccqerecoil than this
  const double MIN_NUE_INEL_CUT = 2*1000.0;        // in MeV. for min Enu for the inelastic cut to get a dis-like sample
  const double MIN_PROB_PLANE_CUT = 0.2; //Using Dipak's most recent number
  // DIS Cut in GeV
  const double MIN_DIS_Q2 = 1.0; 
  const double MIN_DIS_W  = 2.0; 
  
  // Low W Cut in GeV
  const double MIN_LOW_W_W  = 1.5;
  const double MAX_LOW_W_W  = 1.9; 

  // Low Q2 Cut in GeV
  const double MIN_LOW_Q2_Q2  = 0.5;
  const double MAX_LOW_Q2_Q2  = 0.9; 

  // MINOS Cuts
  const double MIN_MINOS_CURVE = 5.0;
  const int MIN_MINOS_PLANE   = 100;
//  const double MAX_MINOS_CURVE = 0.0;
  const double MINOS_COIL_RADIUS = 210; //mm
  const double MAX_MINOS_RADIUS = 2500; //mm

  //same value because we assume negligable theta smearing
  const double MAX_THETA_MU = 17.0; //in deg // cut out events with muon theta wrt beam larger than this
  const double MAX_RECO_THETA_MU = 17.0; //in deg // cut out events with muon theta wrt beam larger than this

  const int FIRST_TRACKER_MOD = 27;
  const int LAST_TRACKER_MOD = 80;

  const double M_neutron = 939.56;
  const double M_proton  = 938.27;
  const double M_nucleon  = ( 1.5*M_neutron + M_proton ) / 2.5; //weighted average because xsec is bigger on n
  const double M_muonMeV = 105.6583755; //wiki
  const double M_muon    = .1057;
  const double M_muon2   = M_muon*M_muon;

  const int PROTON_PDG = 2212;
  const int NEUTRON_PDG = 2112;

  const double Z_scint = .5344;
  const double N_scint = .4656;
  const double Z_carbon = 6.;
  const double N_carbon = 6.;
  const double Z_iron = 26.;
  const double N_iron = 30.;
  const double Z_lead = 82.;
  const double N_lead = 125.;

  const int TRACKER_Z = 99;
  const int TRACKER_Z_EVEN = 98;
  const int TRACKER_Z_ODD = 97;
  const int TRACKER_Z_C  = 91;
  const int TRACKER_Z_FE = 92;
  const int TRACKER_Z_PB = 93;

  const int DEUTERON_Z = 2;

  //=========================
  // Plotting
  //=========================
  const bool WRITE_PRELIMINARY = true; //write preliminary on the plots?
  const bool ADD_HISTO_TITLES  = true; //true for talks and approved plots, false for latex
  const bool WRITE_POT         = true; //
  //const bool GRAYSCALE         = true;  //everything in grayscale
  const bool GRAYSCALE         = false;  //everything in grayscale
  const bool THINSTYLE         = true; //shave down some line widths for vector format plots

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
  const std::string minervame6C  = "minervame6C";  
  const std::string minervame6D  = "minervame6D"; 
  const std::string minervame6E  = "minervame6E"; 
  const std::string minervame6F  = "minervame6F"; 
  const std::string minervame6G  = "minervame6G"; 
  const std::string minervame6H  = "minervame6H"; 
  const std::string minervame6I  = "minervame6I"; 
  const std::string minervame6J  = "minervame6J"; 
 
  const std::string all_nu = "Nu"; 
  const std::string all_antinu = "AntiNu"; 
}

namespace Analysis{ 
  const std::string NukeCC = "NukeCC"; 
}

const std::string DEFAULT_PLAYLIST = Playlist::all_nu; 


  //==========================
  // Namespaces and enums
  //=========================
  /*
  namespace Playlist
  {
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
    
    const std::string minervame1A = "minervame1A";
    const std::string minervame1B = "minervame1B";
    const std::string minervame1C = "minervame1C";
    const std::string minervame1D = "minervame1D";
    const std::string minervame1E = "minervame1E";
    const std::string minervame1F = "minervame1F";
    const std::string minervame1G = "minervame1G";
    const std::string minervame1H = "minervame1H";
    const std::string minervame1L = "minervame1L";
    const std::string minervame1M = "minervame1M";
    const std::string minervame1N = "minervame1N";
    const std::string minervame1O = "minervame1O";
    const std::string minervame1P = "minervame1P";
    const std::string minervame2 = "minervame2";
    const std::string minervame3 = "minervame3";
    
    const std::string minervame5A = "minervame5A";
    const std::string minervame6A = "minervame6A";
    const std::string minervame6B = "minervame6B";
    
    const std::string all_nu = "Nu";
    const std::string all_antinu = "AntiNu";
    
    const std::string all_me_nu = "AllNuME";
    const std::string all_me_antinu = "AntiNuME";
    
    const std::string all_nu_A = "NuA"; //1,7,9,13B - no water
    const std::string all_nu_B = "NuB"; //13C,13D,13E - water
  }
  
  //default playlist
  //const std::string DEFAULT_PLAYLIST = Playlist::minerva13C;
  //const std::string DEFAULT_PLAYLIST = Playlist::all_nu;
  const std::string DEFAULT_PLAYLIST = Playlist::minervame6A;
//    const std::string DEFAULT_PLAYLIST = Playlist::all_me_antinu;
  */
  namespace BG
  {
    const int n_extrapolations = 9;
  }


  //! Rock muon event scan for tgt2 paylist 1
  namespace RockMuon
  {
    const double N_tot_scanned = 975.;
    const double N_rock_found  = 60.;

    //! contamination is the fraction of events that were rock muons
    const double contamination = N_rock_found / N_tot_scanned;

    //! error on contamination as gaussian approximation.  
    const double contamination_err = 
      TMath::Sqrt( 
          TMath::Power( TMath::Sqrt(N_rock_found)/N_tot_scanned, 2 ) +
          TMath::Power( TMath::Sqrt(N_tot_scanned)*N_rock_found/TMath::Power(N_tot_scanned,2), 2 ) 
          )
      * 3.; //scale up by 3 because events will be distributed among ~9 bins
    //! bayesian error up (calcualted with TEfficiency)
    const double contamination_errUp   = 0.008693 * 3.;
    //! bayesian error down (calcualted with TEfficiency)
    const double contamination_errDown = 0.007734 * 3.;
  }

  //=============================
  // Systematic shifted samples
  //=============================
  const int VTX_BLOB_SHIFT = 0; // 0 is CV. shifts 1-4

  namespace MCVariation
  {
    const std::string cv = "central_value";
    const std::string BirksUp    = "BirksConstantUp30pc";
    const std::string BirksDown  = "BirksConstantDown30pc";
    const std::string EFNUCRDown = "EFNUCR-1Sigma";
    const std::string EFNUCRUp   = "EFNUCR+1Sigma";
    const std::string FZONEDown  = "FZONE-1Sigma";
    const std::string FZONEUp    = "FZONE+1Sigma";
    const std::string FZONEAlt   = "FZONE_Alt1";
    const std::string Hadronization_Alt = "Hadronization_Alt1";
    const std::string NoFSI      = "NoFSI";
    const std::string XtalkDown  = "XtalkShift-1Sigma";
    const std::string XtalkUp    = "XtalkShift+1Sigma";

    const std::string DEFAULT = cv;
  }
  //==================================//
  // Helicity enum
  //==================================//
  namespace HelicityType
  {
    enum t_HelicityType
    {
      kAny = 0,
      kNeutrino,
      kAntiNeutrino
    };
  }

  //==================================//
  // File type enum
  //==================================//
  namespace FileType
  {
    enum t_FileType
    {
      kData = 0,    // Default: Assume NukeCC
      kMC,
      kTruth,
      kAny,
      kNukeOnlyMC,
      kNukeOnlyTruth,
      kSpecialMLVFData,
      kSpecialMLVFMC,
      kDNNData,
      kDNNMC,
      kDCNNData,
      kDCNNMC,
      kDANNData,
      kDANNMC


    };
  }


  //=============================
  // Define Channels
  //=============================
  namespace Channels
  {
    const unsigned int size = 10;

    namespace idx
    {
      enum t_ChannelIdx
      {
        CC      = 0, //NumuCC 
        QE      = 1, //NumuCC and GENIE Channel = 1 no charm
        lowW    = 2, //NumuCC  W < 1.3 GeV
	lowQ2Trans  = 3, //NumuCC, Q2 < 1, 1.3 < W < 2 GeV 
        trans   = 4, //NumuCC, Q2 > 1, 1.3 < W < 2 GeV
        softDIS = 5, //NumuCC, Q2 < 1 GeV2 and W > 2 GeV
        DIS     = 6, //NumuCC, Q2 > 1 GeV2 W > 2 GeV
        NC      = 7, //NC and wrong sign events
        WS      = 8, //wrong sign events
        otherCC = 9  //NumuCC Non-delta resonance, coherent and strange production
      };
    }

    namespace Color
    {
      enum t_ChannelColors
      {
        CC      = kBlack,
        QE      = kOrange-4,
        lowW    = kRed-3,
	lowQ2Trans   = kOrange +2,
        trans   = kRed + 2,
        softDIS = kBlue + 2,
        DIS     = kMagenta+3,
        NC      = kBlue + 2,
        WS      = kGreen + 2,
        otherCC = kGreen + 2
      };
    }
  }


  //=============================
  // Define Regions
  //=============================
  namespace Regions
  {
    const unsigned int size = 7; 
    namespace idx
    {
      enum t_RegionIdx
      {
        trueMat = 0,
        Fe    = 1,
        Pb    = 2,
        US    = 3,
        DS    = 4,
        C     = 5,
        Other = 6,
        CH    = 7
      };
    }

    namespace Color
    {
      enum t_RegionColors
      {
        trueMat   = kGreen+2,
        Fe  = kRed+1,
        Pb  = kAzure+4,
        US  = kMagenta+1,
        DS  = kOrange+1,
        C   = kMagenta+3,
        Other = kViolet+2,
        CH  = kPink+1
      };
    }
    namespace FillStyle
    {
      enum t_FillStyle
      {
        trueMat = 3001,
        Fe  = 3001,
        Pb  = 3001,
        US  = 3001,
        DS  = 3001,
        C   = 3001,
        Other = 3001,
        CH  = 3001
      };
    }
  } //end of regions


  
  //=============================
  // Define Multiplicities
  //=============================
  namespace Multiplicities
  {
    const unsigned int size = 4;
    namespace idx
    {
      enum t_Idx
      {
        BadVtx    = 0,
        OneTrk    = 1,
        ExactFit    = 2,
        CloseFit    = 3
      };
    }

    namespace Color
    {
      enum t_Colors
      {
        BadVtx  = kRed+1,
        OneTrk  = kBlue+3,
        ExactFit  = kGreen+2,
        CloseFit  = kMagenta+2,
      };
    }

    namespace FillStyle
    {
      enum t_FillStyle
      {
        BadVtx  = 3001,
        OneTrk  = 3004,
        ExactFit  = 3005,
        CloseFit  = 3006,
      };
    }

    namespace MarkerStyle
    {
      enum t_MarkerStyle
      {
        BadVtx = 24,
        OneTrk = 25,
        ExactFit = 27,
        CloseFit = 28
      };
    }

  } //end of multiplicities


  namespace FSParticles
  {
    const unsigned int size = 10;

    namespace idx
    {
      //index is order in stack
      enum t_Idx
      {
        kOther     = 0,
        kXT        = 1,
        kMuon      = 2,
        kEM        = 3,
        kMeson     = 4,
        kLowNeut   = 5,
        kMidNeut   = 6,
        kHighNeut  = 7,
        kProton    = 8,
        kTotal     = 9
      };
    }

    namespace Color
    {
      enum t_Colors
      {
        kOther     = kGreen+4,
        kXT        = kPink+1,
        kMuon      = kOrange+1,
        kEM        = kYellow+1,
        kMeson     = kBlue+1,
        kLowNeut   = kGray+1,
        kMidNeut   = kCyan+1,
        kHighNeut  = kViolet+1,
        kProton    = kRed+1,
        kTotal     = kGray+1
      };
    }
    const double Colors[size] = { Color::kOther, Color::kXT, Color::kMuon, Color::kEM, Color::kMeson, Color::kLowNeut, Color::kMidNeut, Color::kHighNeut, Color::kProton, Color::kTotal };

    namespace Error
    {
      const double kOther    = .20;
      const double kXT       = .20;
      const double kMuon     = .024;
      const double kEM       = .03;
      const double kMeson    = .05;
      const double kLowNeut  = .25;
      const double kMidNeut  = .10;
      const double kHighNeut = .20;
      const double kProton   = .035;
  }
  const double Errors[size] = { Error::kOther, Error::kXT, Error::kMuon, Error::kEM, Error::kMeson, Error::kLowNeut, Error::kMidNeut, Error::kHighNeut, Error::kProton, 0. };

  namespace FillStyle
  {
    enum t_FillStyle
    {
      kOther     = 3002,
      kXT        = 1001,
      kMuon      = 3144,
      kEM        = 3003,
      kMeson     = 3005,
      kLowNeut   = 1001,
      kMidNeut   = 3144,
      kHighNeut  = 3004,
      kProton    = 3002,
      kTotal     = 0
    };
  }
  const double FillStyles[size] = { FillStyle::kOther, FillStyle::kXT, FillStyle::kMuon, FillStyle::kEM, FillStyle::kMeson, FillStyle::kLowNeut, FillStyle::kMidNeut, FillStyle::kHighNeut, FillStyle::kProton, FillStyle::kTotal };
}

/// What type of information do the bins store?
namespace BinType
{

  enum t_BinType
  {
    kEnergy = 0,
    kPosition,
    kResidual,
    kResolution,
    kSideband,
    kTruth,
    kDeltaE
  };
}


//==============================================
// Define nuclei
//==============================================
namespace Nuclei
{
  const unsigned int size = 4;
  namespace Color
  {
    enum t_Colors
    {
      CH  = kOrange-4,
      Fe  = kRed-2,
      Pb  = kBlue+3,
      C   = kGreen
    };
  }
}


// Unfold Methods
namespace Unfold
{   
  enum t_UnfoldMethod
  {
    kBayes    = 0,
    kSVD      = 1,
    kTUnfold  = 2,
    kInvert   = 3,
    kBinByBin = 4
  };
};  


} //end of NukeCC_Ana namespace

#endif
