#include "CVUniverse.h"
#include <iostream>

// An example of a lateral shift, where we have to change the value of
// one variable (in this case, muon energy). We need to give the
// number of sigma to the constructor
namespace NUKECC_ANA{
class MuonERangeCurvatureShiftUniverse: public CVUniverse
{
public:
  MuonERangeCurvatureShiftUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
    : CVUniverse(chw, nsigma)
    {}

  // MeV
  virtual double GetMuonE() const override { 
    //double muon_E_shift = GetDouble("CCNuPionInc_minosRangeCurveShift");
    double muon_E_shift = GetDouble("NukeCC_sys_muon_energy_shift");
    //double muon_E_shift = NukeCC_sys_muon_energy_shift[1];
    double shift_val    = m_nsigma*muon_E_shift;
    return shift_val+CVUniverse::GetMuonE();
  }

  virtual std::string ShortName() const { return "EmuRangeCurve"; }
  virtual std::string LatexName() const { return "MINOS Muon Energy - Range & Curvature"; }
};
}
