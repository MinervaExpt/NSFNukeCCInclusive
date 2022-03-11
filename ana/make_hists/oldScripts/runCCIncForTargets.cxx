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
#include "PlotUtils/TargetUtils.h"

#include <cstdlib>


typedef unsigned int uint;

template<class T> T sqr(T x) { return x*x; }

double GetNormValue( int targetID, int targetZ, int targetNucleon = 0 )
{
  PlotUtils::TargetUtils targetInfo;
  //const double tracker_mass = PlotUtils::TargetUtils::Get().GetTrackerMass( 0, PlotUtils::TargetProp::Tracker::Back, true, 850.);
  //const double trackerAtomsC  = tracker_mass * PlotUtils::TargetUtils::Get().GetTrackerElementMassFraction( 6, true ) * PlotUtils::TargetProp::AtomsPerGram::C;
  double trackerAtomsC = targetInfo.Get().GetTrackerElementNAtoms( 6, PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, true, 850.0);
  //double trackerAtomsC = PlotUtils::TargetUtils::Get().GetTrackerElementNAtoms( 6, 5990., 8340., true, 850. );
  //targetInfo.GetTrackerElementNAtoms( 6, 5980, 8422, true );

  //double trackerAtomsC = targetInfo.GetTrackerElementNAtoms( 6, 5991.37, 8343.28, true );
  //5990., 8340., 1 ); //z 5990-8340 = [27,78]=52 modules
  std::cout <<  "ATTENTION" << endl;
  std::cout<<trackerAtomsC<<endl;  
  double passiveNucleons = -999.;    
    
  if ( targetNucleon == 0 ){     
        
    if(targetID < 10 ){
            
       passiveNucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons( targetID, targetZ, true, 850. );
       //targetInfo.GetPassiveTargetNNucleons(targetID,targetZ, true);
       std::cout <<  "ATTENTION" << endl;
       std::cout <<   passiveNucleons  << endl;
       
            
    }
    if(targetID > 10){
      passiveNucleons = targetInfo.GetTrackerNNucleons( 12., true, 850. );
    }
        
  }

  else
  {
    assert( false && "Target nucleons can be all, proton or neutron only" );
  }
    
  if( passiveNucleons < 0 )
    assert( false && "Normalizations only known for Z = 1,6,26,82 and 0(scint)" );
    
  cout<<"The normalization factor for this analysis is "<<trackerAtomsC / passiveNucleons<<endl;
    
  return trackerAtomsC / passiveNucleons;
}


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
        
        if( 2000.0 >= muon_E || muon_E >= 50000.0 ) return false;
        
        //if( muon_theta*rad_to_deg >= 17.0 ) return false;

        return true;
    }

    bool InHex( ChainWrapper& chw, int entry, double apothem )
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

    bool PassTrueDistToDivisionCut(ChainWrapper& chw, int entry, double xySep /* = 25. */ )
    {
      double true_target_dist_to_division = (double)chw.GetValue("truth_target_dist_to_division",entry);
      int true_targetID = (int)chw.GetValue("truth_targetID",entry);

      //only relevant for passive targets 1235
      if( 0 < true_targetID && true_targetID  < 10 && 4 != true_targetID)
        return ( xySep < true_target_dist_to_division);
    
      return true;
    }

    bool IsInTrueMaterial(ChainWrapper& chw, int entry, const int i_targetID, const int i_targetZ )
    {
      if(i_targetID < 10 )
      {
          // If targetID < 10, then you are looking for a passive target event.
          // Require that the event has the same targetID and targetZ.
          //if( truth_targetID == i_targetID )
          if( (int)chw.GetValue("truth_targetID",entry) == i_targetID )
          {
              if( i_targetZ > 0 )
                return (int)chw.GetValue("truth_targetZ",entry) == i_targetZ;
              else
                return true;
          }
      }
      return false;
    }


    virtual bool passesCuts(ChainWrapper& chw, int entry)
    {
      if((int)chw.GetValue("mc_incoming", entry)!=-14) return false;
      if((int)chw.GetValue("mc_current", entry)!=1) return false;
      if(!InHex(chw,entry,850.0)) return false;
      if(!PassTrueDistToDivisionCut( chw, entry, 25.0)) return false;
      if(!IsInTrueMaterial( chw, entry, 3, 26)) return false;
      if(!passesKinematics(chw, entry)) return false;

      // Target

      return true;
    }

    int m_target;
    protected:
        double m_max_E_avail;
};


int runXSecLooper(const bool antinu, const double Emin, const double Emax, const std::vector<const char*> fileNames)
{
  const char* fileName = "GENIEXSecExtract_CCInclusive_T3Iron.root";
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

  loop.setFiducial(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, 850);
  //loop.setFiducial( true, PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back);//z 5990-8340 = [27,78]=52 modules
  //loop.setFiducial(false, 5980, 8422);//z 5990-8340 = [27,78]=52 modules
  loop.setPlaylist(PlotUtils::FluxReweighter::minervame6A);

  const int nBinsEnu=14;
  double binsEnu[nBinsEnu+1]={ 2.0, 3.0, 4.0, 5.0, 6.25, 7.5, 8.75, 10., 12.5, 15., 20., 25.0, 30., 40., 50. };

  CCIncXSec* total_xSec_E = new CCIncXSec("total_Enu");
  total_xSec_E->setVariable(XSec::kENu);
  total_xSec_E->setDimension(1);
  total_xSec_E->setBinEdges(nBinsEnu, binsEnu);
  total_xSec_E->setNormalizationValue( GetNormValue( 3, 26 ) );
  total_xSec_E->setIsFluxIntegrated(false); // total cross-section
  //total_xSec_E->setFluxIntLimits(Emin, Emax);
  std::cout <<  "ATTENTION" << endl;
  std::cout <<  "ATTENTION" << endl;
  std::cout <<  GetNormValue( 3, 26 ) << endl;
  total_xSec_E->setNormalizationType(XSec::kSelfNorm);
  total_xSec_E->setUniverses(0);//default value, put 0 if you do not want universes to be included.

  CCIncXSec* dif_xSec_E = new CCIncXSec("dif_Enu");
  dif_xSec_E ->setVariable(XSec::kENu);
  dif_xSec_E ->setBinEdges(nBinsEnu, binsEnu);
  dif_xSec_E ->setIsFluxIntegrated(true); // total cross-section
  dif_xSec_E->setNormalizationValue( GetNormValue( 3, 26 ) );
  dif_xSec_E ->setDimension(1);
  dif_xSec_E ->setFluxIntLimits(Emin, Emax);
  dif_xSec_E ->setNormalizationType(XSec::kSelfNorm);
  dif_xSec_E ->setUniverses(0);//default value, put 0 if you do not want universes to be included.

  const int nBinsx=6;
  double binsx[nBinsx+1]={ 0.001, 0.05, 0.1, 0.2, 0.4, 1, 2.2 };

  CCIncXSec* dif_xSec_x = new CCIncXSec("dif_x");
  dif_xSec_x ->setVariable(XSec::kxExp); 
  dif_xSec_x ->setDimension(1);  
  dif_xSec_x ->setBinEdges(nBinsx, binsx);
  dif_xSec_x ->setIsFluxIntegrated(true); // total cross-section
  dif_xSec_x ->setFluxIntLimits(Emin, Emax);
  dif_xSec_x->setNormalizationValue( GetNormValue( 3, 26 ) );
  dif_xSec_x ->setNormalizationType(XSec::kSelfNorm);
  dif_xSec_x ->setUniverses(0);//default value, put 0 if you do not want universes to be included.
  

  loop.addXSec(total_xSec_E);
  loop.addXSec(dif_xSec_E);
  loop.addXSec(dif_xSec_x);


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

  std::vector<const char*> fileNames(argv+1, argv+argc);

  bool antinu=true;

  int Emin=0; //GeV
  int Emax=120; //GeV 

  return runXSecLooper(antinu, Emin, Emax, fileNames);
}  