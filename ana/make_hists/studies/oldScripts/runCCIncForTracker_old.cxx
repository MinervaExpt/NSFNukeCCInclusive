// Copied by a.klustova20@imperial.ac.uk from code written by Andrew Olivier, Phil Rodriges, Rik Gran, Minerba Betancort, and maybe others.
// January 2022

#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "Math/AxisAngle.h"
#include "Math/Vector3D.h"

#include "TH1.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

#include <cstdlib>

typedef unsigned int uint;

template<class T> T sqr(T x) { return x*x; }

class CCIncXSec : public XSec
{
public:

  CCIncXSec(const char* name)
    : XSec(name)
    {}

    bool passesKinematics(ChainWrapper& chw, int entry)
    {
        const double PI = atan(1.0)*4;
        const double deg_to_rad = PI / 180.0;
        const double rad_to_deg = 1. / deg_to_rad;
        
        double muon_E     = (double)chw.GetValue("mc_primFSLepton",entry,3);
        double muon_theta = (double)chw.GetValue("truth_muon_theta",entry);
        ROOT::Math::AxisAngle toBeamFrame(ROOT::Math::XYZVector(1., 0., 0.), -0.05887); //NuMI beam angle in mrad from PlotUtils
        ROOT::Math::XYZVector muon(chw.GetValue("mc_primFSLepton", entry, 0), //In MeV
                               chw.GetValue("mc_primFSLepton", entry, 1),
                               chw.GetValue("mc_primFSLepton", entry, 2));
        
        const double thetamu = (toBeamFrame * muon).theta();
        if(thetamu*rad_to_deg>=17.0) return false; //This is a YAML parameter
        
        if( 2000.0 >= muon_E || muon_E >= 50000.0 ) return false;
        
        //if( muon_theta*rad_to_deg >= 17.0 ) return false;

        return true;
    }

    virtual bool passesCuts(ChainWrapper& chw, int entry)
    {
        if((int)chw.GetValue("mc_incoming", entry)!=-14) return false;
        if((int)chw.GetValue("mc_current", entry)!=1) return false;
        if(!passesKinematics(chw, entry)) return false;

        // Tracker x cut, energy cut, angle cut, CC, antineutrino

        return true;
    }

    protected:
        double m_max_E_avail;
};

int runXSecLooper(const bool antinu, const double Emin, const double Emax, const std::vector<const char*> fileNames)
{
  const char* fileName = "GENIEXSecExtract_CCInclusive_Tracker.root";
  auto outFile = TFile::Open(fileName, "CREATE");
  if(!outFile)
  {
    std::cerr << "GENIEXSecExtract failed to open a file for output named "
              << fileName << ".  Aborting the cross section loop...\n";
    return 1;
  }

  // Create the XSecLooper and tell it the input files
  if(fileNames.empty())
  {
    std::cerr << "No files for GENIEXSecExtract to read.  Aborting the cross section loop...\n";
    return 2;
  }

  XSecLooper loop(fileNames[0]);
  for(auto whichFile = fileNames.begin()+1; whichFile != fileNames.end(); ++whichFile)
  {
    loop.addFiles(*whichFile);
  }

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(antinu ? -14: 14);

  //Setting the number of Universes in the GENIE error band, default 100 universes put 0 if you do not want universes to be included
  loop.setNumUniv(0); 

  loop.setFiducial(5980, 8422);
  loop.setPlaylist(PlotUtils::FluxReweighter::minervame6A);

  const int nBins=14;
  double bins[nBins+1]={ 2.0, 3.0, 4.0, 5.0, 6.25, 7.5, 8.75, 10., 12.5, 15., 20., 25.0, 30., 40., 50. };

  CCIncXSec* total_xSec_E = new CCIncXSec("Enu");
  total_xSec_E->setVariable(XSec::kENu);
  total_xSec_E->setDimension(1);
  total_xSec_E->setBinEdges(nBins, bins);
  total_xSec_E->setIsFluxIntegrated(false); // total cross-section
  total_xSec_E->setFluxIntLimits(Emin, Emax);
  total_xSec_E->setNormalizationType(XSec::kPerNucleon);
  total_xSec_E->setUniverses(0);//default value, put 0 if you do not want universes to be included.

  loop.addXSec(total_xSec_E);

  // Once everything's set up, actually run the thing
  loop.runLoop();

  // Get the output histograms and save them to file
  outFile->cd();
  for(uint i=0; i<loop.getXSecs().size(); ++i){
    loop.getXSecs()[i]->getXSecHist()->Write();
    loop.getXSecs()[i]->getEvRateHist()->Write();
  }

  return 0;
}

int main(int argc, char** argv)
{

  TH1::AddDirectory(kFALSE); //Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).

  std::vector<const char*> fileNames(argv+1, argv+argc);

  bool antinu=true;

  int Emin=2; //GeV
  int Emax=50; //GeV 

  return runXSecLooper(antinu, Emin, Emax, fileNames);
}  