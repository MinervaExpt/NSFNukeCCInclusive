#ifndef Systematics_h
#define Systematics_h

//==============================================================================
// Get Several standard MINERvA systematics
//==============================================================================

#include "../../../NUKECCSRC/include/CommonIncludes.h"
#include "../../../NUKECCSRC/include/CVUniverse.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/FluxReweighter.h"
using namespace NUKECC_ANA;

typedef std::map<std::string, std::vector<CVUniverse*> > SystMap;

// Error Bands
std::map<std::string, std::vector<CVUniverse*> > GetErrorBands(PlotUtils::ChainWrapper* chain) {
  
  // return map
  SystMap error_bands;

  // CV 
  error_bands[std::string("CV")].push_back(new CVUniverse(chain));

  if(RunCodeWithSystematics) {
  
    //========================================================================
    // FLUX
    //========================================================================
    SystMap flux_systematics = 
    PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain,CVUniverse::GetNFluxUniverses());
  
    error_bands.insert(flux_systematics.begin(), flux_systematics.end());

  }

  return error_bands;
}
#endif  // Systematics_h