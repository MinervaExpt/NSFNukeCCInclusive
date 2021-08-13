#ifndef MNV_NUKECCUTILSNSF_h
#define MNV_NUKECCUTILSNSF_h 1


//#include "GlobalIncludes.h"
//#include "GlobalParameters.h" 
//#include "CCQENuCutsNSF.h"

#include "include/NukeCC_Cuts.h"
//#include "CCQENuBinning.h"
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
#include "PlotUtils/MacroUtil.h"
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

  class NukeCCUtilsNSF {

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
    
    NukeCCUtilsNSF();
    
    NukeCCUtilsNSF(string playlist);

    //! Constructor with parameters
    NukeCCUtilsNSF( bool use_mergedFiles,string playlist );

    NukeCCUtilsNSF( bool use_mergedFiles, bool useFluxConstraint );

    //! Default Destructor
    ~NukeCCUtilsNSF();

    //! singleton gettor
    static NukeCCUtilsNSF& Get();

    public:
    // Some units
    StrHist1DWrapperMap cvhistos1D;
    StrHist2DWrapperMap cvhistos2D;
    
    //cutter
    NukeCC_Cuts *cutter;
    
    int incoming_pdg;

    //set the pot profile....
    double C_global_used_pot_mc;
    double C_global_tot_pot_mc;



    public:
    
    std::string SetPlaylist(std::string playlist);

   void getPOT(double& total_pot_data,double& total_pot_mc );
   void writePOT(TFile *f, double DataPOT, double MCPOT);

HelicityType::t_HelicityType GetHelicityFromPlaylist( std::string playlist ) const;


std::vector<std::string> GetStdPlaylists( HelicityType::t_HelicityType helicity ) const;

TString GetHistFileName( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ, HelicityType::t_HelicityType helicity) const;

//NEW for grid
TString GetHistFileNameGrid( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ) const;

TString HistDir( bool forceDisk  = false , bool ignorePlaylist = false, bool needsbluearc = false )const;

TString GetVariationTag() const;
std::string GetHelicityString( HelicityType::t_HelicityType helicity ) const;

TString  GetHistName( const std::string& histType, FileType::t_FileType fType, const std::string& var, int targetID, int targetZ ) const;
 
std::string GetFileTypeString( FileType::t_FileType fType ) const;
  
 //getting the Chainwrapper pointers....
  
  //POT stuff...


 

  };
}
#endif
