#ifndef MNV_CCQENUUTILSNSF_h
#define MNV_CCQENUUTILSNSF_h 1


#include "GlobalIncludes.h"
//#include "GlobalParameters.h" 
//#include "CCQENuCutsNSF.h"

#include "include/NukeCC_Cuts.h"
#include "CCQENuBinning.h"
//#include "include/ComputeUtils.h"
//#include "include/NeutronBlob.h"
//#include "include/NeutronBlob.h"
//#include "include/NeutronBlobBinning.h"
//#include "include/NeutronBlobCuts.h"
#include "CVUniverse.h"
#include "PlotUtils/DefaultCVUniverse.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"

#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvH3D.h"
#include "PlotUtils/AnaBinning.h"



using namespace PlotUtils;

// HadronReweight
namespace PlotUtils{
  class MnvHadronReweight;
  class FluxReweighter;
}

namespace NUKECC_ANA{

  class CCQENuUtilsNSF {

    typedef std::map<std::string,std::vector<double>*> StrVecMap;
    //! string-MnvND type for persistent map definition  
    //outside event loop
    typedef std::map<std::string, HistWrapper<CVUniverse>**> StrHist1DWrapperMap;
    typedef std::map<std::string, Hist2DWrapper<CVUniverse>**> StrHist2DWrapperMap;

private:
       //Michel Weights 
      std::vector<double> michel_weights; 

      //Target Mass Weights 
      std::vector<double> target_mass_weights; 

     //Proton Tracking Efficiency Weights 
      std::vector<double> track_eff_weights;  
      bool use_merged_files;


public:
    //! Default Constructor
    
    CCQENuUtilsNSF();
    
    CCQENuUtilsNSF(string playlist);

    //! Constructor with parameters
    CCQENuUtilsNSF( bool use_mergedFiles,string playlist );

    CCQENuUtilsNSF( bool use_mergedFiles, bool useFluxConstraint );

    //! Default Destructor
    ~CCQENuUtilsNSF();

    //! singleton gettor
    static CCQENuUtilsNSF& Get();

    public:
    // Some units
    StrHist1DWrapperMap cvhistos1D;
    StrHist2DWrapperMap cvhistos2D;
    
    //cutter
      NukeCC_Cuts *cutter;
    
    int incoming_pdg;
    //MinervaUnits::numi_beam_angle_rad
  /*  static constexpr double numi_beam_angle_rad =  -0.05887;

    static constexpr double M_mu    = 105.6583;
    //CLHEP/Units/PhysicalConstants.h
    static constexpr double M_p     = 938.272013;
    static constexpr double M_n     = 939.56536; 
    static constexpr double M_nucleon  = ( 1.5*M_n + M_p ) / 2.5;

    static constexpr int PROTON_PDG = 2212;
    static constexpr int NEUTRON_PDG = 2112;

    //Genie tune variables - default is on!//
    bool   doGenieTuning;
    static constexpr double genieMaRes              = 1.12;
    static constexpr double genieMaRes1sig          = 0.2 * 1.12;
    // GENIE central value MvRES from electroproduction data fit
    static constexpr double genieMvRes              = 0.84;
    static constexpr double genieMvRes1sig          = 0.1 * 0.84;
    // Reduced MvRES error from electroproduction data fit
    static constexpr double electroProdMvRes1sig    = 0.03 * 0.84;
    // Pion production parameters and errors from deuterium fit
    static constexpr double deuteriumMaRes          = 0.94;
    static constexpr double deuteriumMaRes1sig      = 0.05;
    static constexpr double deuteriumNonResNorm     = 0.43; //See docdb 11524 and 11567 for central value and error band modification values
    static constexpr double deuteriumNonResNorm1sig = 0.04;
    static constexpr double deuteriumResNorm        = 1.15;
    static constexpr double deuteriumResNorm1sig    = 0.07;
    */



    public:
    
    std::string SetPlaylist(std::string playlist);

HelicityType::t_HelicityType GetHelicityFromPlaylist( std::string playlist ) const;


std::vector<std::string> GetStdPlaylists( HelicityType::t_HelicityType helicity ) const;

 TString GetHistFileName( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ, HelicityType::t_HelicityType helicity, bool forceDisk = false, bool usebluearc = false ) const;

//TString GetHistFileName( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ, HelicityType::t_HelicityType helicity, bool forceDisk, bool usebluearc ) const;

TString HistDir( bool forceDisk  = false , bool ignorePlaylist = false, bool needsbluearc = false )const;

TString GetVariationTag() const;
std::string GetHelicityString( HelicityType::t_HelicityType helicity ) const;

TString  GetHistName( const std::string& histType, FileType::t_FileType fType, const std::string& var, int targetID, int targetZ ) const;
 
std::string GetFileTypeString( FileType::t_FileType fType ) const;
     double GetTotalMomentum( double px, double py, double pz ){ }
    //==============================================================
    // BookHistos (to be used int make_hists macros)
    //==============================================================
//==============================================================
    // BookHistos (to be used int make_hists macros)
    //==============================================================
      void bookHistosCV( HistWrapper<CVUniverse>** h, string var_name, string title, std::vector<double> nbin, double min, double max,std::map< std::string, std::vector<CVUniverse*> > error_bands );
/*      
      void bookHistosCV( HistWrapper<CVUniverse>** h, string var_name, string title, axis_binning xbins,std::map< std::string, std::vector<CVUniverse*> > error_bands );
      
      void bookHistosCVByPID( HistWrapper<CVUniverse>** h, string var_name, string title, int nbins, double min, double max,std::map< std::string, std::vector<CVUniverse*> > error_bands );
      
      void bookHistosCVByPID( HistWrapper<CVUniverse>** h, string var_name, string title, axis_binning xbins,std::map< std::string, std::vector<CVUniverse*> > error_bands );
      
      void bookHistosCV( Hist2DWrapper<CVUniverse>** h,string var_name,string title,axis_binning xbins, axis_binning ybins,std::map< std::string, std::vector<CVUniverse*> > error_bands );
      
      void bookHistosCV( Hist2DWrapper<CVUniverse>** h,string var_name,string title,int xnbins,double xmin,double xmax, int ynbins, double ymin, double ymax,std::map< std::string, std::vector<CVUniverse*> > error_bands );       
      
      void bookHistosCVByPID( Hist2DWrapper<CVUniverse>** h,string var_name,string title,axis_binning xbins, axis_binning ybins,std::map< std::string, std::vector<CVUniverse*> > error_bands );
      
      void bookHistosCVByPID( Hist2DWrapper<CVUniverse>** h,string var_name,string title,int xnbins,double xmin,double xmax, int ynbins, double ymin, double ymax,std::map< std::string, std::vector<CVUniverse*> > error_bands );      
 
 
*/

  //fill functions....
     //these ones to fill the mc
    // void fillHistos( Hist2DWrapper<CVUniverse>** h, double var_x, double var_y,CVUniverse* universe, double w );
   
     void fillHistos( HistWrapper<CVUniverse>** h, double var,CVUniverse* universe, double w);  
  
  /*
     //these ones to fill the data
     void fillHistos( Hist2DWrapper<CVUniverse>** h, double var_x, double var_y );
     
     void fillHistos( HistWrapper<CVUniverse>** h, double var);
    */ 

     
  //sync functions...
  void syncCVHistos(HistWrapper<CVUniverse>**h);   
  
  
 //getting the Chainwrapper pointers....
  
  ChainWrapper* GetChainWrapperMCPointer(string playslist,string tree_name = "NukeCC");
  ChainWrapper* GetChainWrapperDataPointer(string playlist,string tree_name = "NukeCC");
  
  
  //get the error bands.....
  std::map< std::string,  std::vector<CVUniverse*>> GetErrorBands(ChainWrapper* chain);
  
  //POT stuff...

  
  void getPOT(double& total_pot_data,double& total_pot_mc );
  double setPOTData(ChainWrapper* chain_data);
  double setPOTMC(ChainWrapper* chain_mc);
  void writePOT(TFile *f);
  

  };
}
#endif
