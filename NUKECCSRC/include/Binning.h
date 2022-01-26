//#ifndef MNV_NUKECC_cxx
//#define MNV_NUKECC_cxx 1

#ifndef MNV_NUKECC_BINNING_h
#define MNV_NUKECC_BINNING_h 1

#include "GlobalIncludes.h"
#include "PlotUtils/MinervaUniverse.h"
#include "PlotUtils/TargetUtils.h"
#include "CVUniverse.h"
#include "CommonIncludes.h"
//#include "include/GlobalIncludes.h"


using namespace std;

namespace NUKECC_ANA{

class NukeCC_Binning{

public:

//constructor
NukeCC_Binning(); 

//Destructor

~NukeCC_Binning();

static NukeCC_Binning &Get();

       
      static double fDetBins[94];
      static const unsigned int fNDetBins = 93;
      //! What is the bin width normalized to "Event per X"
      double GetVarNormWidth( const std::string& var ) const;

      //! What is the lowest value we want to show
      double GetVarMinVal( const std::string& var ) const;

      //! What is the highest value we want to show
      double GetVarMaxVal( const std::string& var ) const;
      
      //! Min Val for DIS analysis
      double GetDISVarMinVal( const std::string& var ) const;

      //! Max Val for DIS analysis
      double GetDISVarMaxVal( const std::string& var ) const;

      //! Add bins to the vector in preparation of variable-sized bins array
      void AddBins( std::vector<double>& binsLowEdge, const int nBins, const double binWidth ) const;

      //! Add bins to the vector in preparation of variable-sized bins array
      void AddBins( std::vector<double>& binsLowEdge, const double binWidth, const double lowBin, const double highBin ) const;


      //!Given a value and a particular var, return the bin number of that value.
      unsigned int FindBin(std::string var, double Value) const;

      //! Generic get bins with specified type of thing to bin
      std::vector<double> GetBins( const std::string& var, BinType::t_BinType type = BinType::kEnergy ) const;

      //! Get bins for plots of energy variables or theta mu
      std::vector<double> GetEnergyBins( const std::string& var, bool DISbins = true ) const;

      //! Get bins for plots of residual of energy variables (reco - true)
      std::vector<double> GetResidualBins( const std::string& var ) const;

      //! Get bins for plots of resolution of energy variables (reco - true) / true
      std::vector<double> GetResolutionBins( const std::string& var ) const;

      //! Get bins for a position plot
      std::vector<double> GetPosBins( const std::string& var, const int targetID = -1 ) const;
      
      double GetBGZ( const int targetID ) const;
      //! Get bins for DIS
      std::vector<double> GetDISBins( const std::string& var ) const;
      
      //! Get bins for sidebands
      std::vector<double> GetSidebandBins( const std::string& var ) const;
      
      //! Get truth-level bins
      vector<double> GetTruthBins(const std::string& var) const;

      //! Get delta E bins -- because I'm lazy
      vector<double> GetDeltaEBins(const std::string& var) const;
      
      

      //! What is the X axis title for this variable
      std::string GetXaxisTitle( const std::string& var ) const;

      //! What is the Y axis title for this variable
      std::string GetYaxisTitle( const std::string& var ) const;
};


}
#endif //MNV_NUKECC_cxx

