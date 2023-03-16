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

#include "NukeCC_bins.cc"

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
        //double muon_theta = (double)chw.GetValue("truth_muon_theta",entry);
        ROOT::Math::AxisAngle toBeamFrame(ROOT::Math::XYZVector(1., 0., 0.), -0.05887); //NuMI beam angle in mrad from PlotUtils
        ROOT::Math::XYZVector muon(chw.GetValue("mc_primFSLepton", entry, 0), //In MeV
                              chw.GetValue("mc_primFSLepton", entry, 1),
                              chw.GetValue("mc_primFSLepton", entry, 2));
        
        const double thetamu = (toBeamFrame * muon).theta();
        if(thetamu*rad_to_deg>=17.0) return false; //This is a YAML parameter
        
        if( 2000.0 >= muon_E || muon_E >= 20000.0 ) return false;
        
        //if( muon_theta*rad_to_deg >= 17.0 ) return false;

        return true;
    }

    virtual bool InHex( ChainWrapper& chw, int entry, double apothem )
    {
      if( apothem == 0. )
        return false;
        
      double x = (double)chw.GetValue("mc_vtx",entry,0);
      double y = (double)chw.GetValue("mc_vtx",entry,1);
        
      //Hexagon is symmetric about its x and y
      x = fabs(x);
      y = fabs(y);
        
      if( pow(x,2) + pow(y,2) < pow(apothem,2) )
        return true;
        
      double lenOfSide = apothem * ( 2 / sqrt(3) );
        
      if( x > apothem )
        return false;
        
      if( y < lenOfSide/2.0 )
        return true;
        
      double slope = (lenOfSide / 2.0) / apothem;
      if( y < lenOfSide - x*slope )
        return true;
        
      return false;
    }

    virtual bool passesCuts(ChainWrapper& chw, int entry)
    {
        if( (double)chw.GetValue("mc_vtx", entry, 2) < 5980. || (double)chw.GetValue("mc_vtx", entry, 2) > 8422. ) return false;
        if(!InHex(chw,entry,850.0)) return false;
        // fiducial
        if((int)chw.GetValue("mc_incoming", entry)!=-14) return false;
        if((int)chw.GetValue("mc_current", entry)!=1) return false;
        if(!passesKinematics(chw, entry)) return false;

        // Tracker x cut, energy cut, angle cut, CC, antineutrino

        return true;
    }

    protected:
        double m_max_E_avail;
};

int runXSecLooper(const bool antinu, const double Emin, const double Emax, const std::vector<const char*> fileNames, string outdir)
{
  const char* fileName = Form("%s/GENIEXSecExtract_CCInclusive_t99_z99.root", outdir.c_str());
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

  vector<XSec::EVariable> vars;
  vars.push_back( XSec::kENu );
  vars.push_back( XSec::kxExp );
  vars.push_back( XSec::kThetaLep );

  for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ){
                
    string varName = "";
    switch(*var)
    {
        case XSec::kxExp:
            varName = "x";
            break;
        case XSec::kENu:
            varName = "Enu";
            break;
        case XSec::kThetaLep:
            varName = "ThetamuDeg";
            break;
    }

    char* name = Form( "tracker_99_%s_std", varName.c_str());
    cout<<"making cross section for tracker" <<endl;
    CCIncXSec* xsec = new CCIncXSec(name);
    
    //! container for bins
    std::vector<double> bins;   
    GetBins( *var, bins, false, false );
    int num_bins = (int)bins.size()-1;
    //! set cross section data
    xsec->setBinEdges(num_bins,&bins.front());
    xsec->setVariable(*var);
    //if(*var !=XSec::kENu) {xsec->setIsFluxIntegrated(false);}
    //else 
    xsec->setIsFluxIntegrated(true);
    xsec ->setDimension(1); 
    xsec->setFluxIntLimits(0.,120.);
    xsec->setUniverses(0);
    xsec->setNormalizationType(XSec::kPerNucleon);
    //! add cross section
    loop.addXSec(xsec);
  }
  //const int nBinsEnu=14;
  //double binsEnu[nBinsEnu+1]={ 2.0, 3.0, 4.0, 5.0, 6.25, 7.5, 8.75, 10., 12.5, 15., 20., 25.0, 30., 40., 50. };

  //CCIncXSec* total_xSec_E = new CCIncXSec("total_Enu");
  //total_xSec_E->setVariable(XSec::kENu);
  //total_xSec_E->setDimension(1);
  //total_xSec_E->setBinEdges(nBinsEnu, binsEnu);
  //total_xSec_E->setIsFluxIntegrated(false); // total cross-section
  //total_xSec_E->setFluxIntLimits(Emin, Emax);
  //total_xSec_E->setNormalizationType(XSec::kPerNucleon);
  //total_xSec_E->setUniverses(0);//default value, put 0 if you do not want universes to be included.

  /*CCIncXSec* dif_xSec_E = new CCIncXSec("tracker_99_Enu_std");
  dif_xSec_E ->setVariable(XSec::kENu);
  dif_xSec_E ->setDimension(1);
  dif_xSec_E ->setBinEdges(nBinsEnu, binsEnu);
  dif_xSec_E ->setIsFluxIntegrated(true); // total cross-section
  dif_xSec_E ->setFluxIntLimits(Emin, Emax);
  dif_xSec_E ->setNormalizationType(XSec::kPerNucleon);
  dif_xSec_E ->setUniverses(0);//default value, put 0 if you do not want universes to be included.

  const int nBinsx=6;
  double binsx[nBinsx+1]={ 0.001, 0.05, 0.1, 0.2, 0.4, 1, 2.2 };

  CCIncXSec* dif_xSec_x = new CCIncXSec("tracker_99_x_std");
  dif_xSec_x ->setVariable(XSec::kxExp); 
  dif_xSec_x ->setDimension(1);  
  dif_xSec_x ->setBinEdges(nBinsx, binsx);
  dif_xSec_x ->setIsFluxIntegrated(true); // total cross-section
  dif_xSec_x ->setFluxIntLimits(Emin, Emax);
  dif_xSec_x ->setNormalizationType(XSec::kPerNucleon);
  dif_xSec_x ->setUniverses(0);//default value, put 0 if you do not want universes to be included.
  

  //loop.addXSec(total_xSec_E);
  loop.addXSec(dif_xSec_E);
  loop.addXSec(dif_xSec_x);
  */


  // Once everything's set up, actually run the thing
  loop.runLoop();

  // Get the output histograms and save them to file
  outFile->cd();
  for(uint i=0; i<loop.getXSecs().size(); ++i){
    loop.getXSecs()[i]->getXSecHist()->Write();
    loop.getXSecs()[i]->getEvRateHist()->Write();
  }
  loop.getFluxHist()->Write();

  return 0;
}

int main(int argc, char** argv)
{

  TH1::AddDirectory(kFALSE); //Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).

  std::vector<const char*> fileNames(argv+2, argv+argc);

  bool antinu=true;

  int Emin=0; //GeV
  int Emax=120; //GeV 

  string outdir = argv[1];

  return runXSecLooper(antinu, Emin, Emax, fileNames, outdir);
}  