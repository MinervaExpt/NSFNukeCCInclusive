#include "GENIEXSecExtract/XSecLooper.h"

#include <cstdlib>
#include "TFile.h"
#include "TH1.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/TargetUtils.h>
#include "Math/AxisAngle.h"
#include "Math/Vector3D.h"

#include "NukeCC_bins.cc"

using namespace std;
using namespace PlotUtils;

namespace {
    const double TRACKER_ZMIN = 6117;
    const double TRACKER_ZMAX = 8193;
    const int targetZ = 26;
}
bool debug=false;
bool anyTracker=false;//true;

double GetNormValue(int targetZ, int targetNucleon = 0 )
{
  cout<<"GetNormValue::l21"<<endl;
  double trackerAtomsC = TargetUtils::Get().GetTrackerElementNAtoms( 6, 5970, 8450, true, 850.0);

  //double trackerAtomsC = TargetUtils::Get().GetTrackerElementNAtoms( 6, 5990., 8340., true, 850. ); //z 5990-8340 = [27,78]=52 modules
  //double trackerAtomsC = TargetUtils::Get().GetTrackerElementNAtoms( 6, 106., 1 ); //[27,79] is 53 modules, not 52.
  //double trackerAtomsC = TargetUtils::Get().GetTrackerElementNAtoms( 6, 104., 0, 850. );//no idea why you would want ot set isMC to false but atm I'm just trying to reproduce Annes flux. 
  cout<<"trackerAtomsC: " <<trackerAtomsC<<endl;
    //cout<<trackerAtomsC<<"  "<<2.22311e+27 * 104.<<endl;
  //cout<<trackerAtomsC<<"  "<<2.22311e+27 * 108.<<endl;//not sure what this number is, its not the C12 number for a single plane. 
    
    
    //    const double trackerAtomsC = 2.22311e+27 * 92.;
    //     const double trackerAtomsC = 2.22311e+27 * 108.;
    double passiveNucleons = 0;
    
    if ( targetNucleon == 0 ){

      if(targetZ == 6 ){
        cout<<"GetNormValue::l34"<<endl;
        passiveNucleons = TargetUtils::Get().GetPassiveTargetNNucleons(3,targetZ, true, 850.);
      }
        
      if(targetZ == 26 ){
        cout<<"GetNormValue::l34"<<endl;
        passiveNucleons = TargetUtils::Get().GetPassiveTargetNNucleons(2, targetZ, true, 850.) +  
                          TargetUtils::Get().GetPassiveTargetNNucleons(3,targetZ, true, 850.) + 
                          TargetUtils::Get().GetPassiveTargetNNucleons(5,targetZ, true, 850.);
      }
      if(targetZ == 82 ){
        cout<<"GetNormValue::l34"<<endl;
        passiveNucleons = TargetUtils::Get().GetPassiveTargetNNucleons(2, targetZ, true, 850.) +  
                          TargetUtils::Get().GetPassiveTargetNNucleons(3,targetZ, true, 850.) + 
                          TargetUtils::Get().GetPassiveTargetNNucleons(4,targetZ, true, 850.) + 
                          TargetUtils::Get().GetPassiveTargetNNucleons(5,targetZ, true, 850.);
      }
      cout<<"target "<<targetZ<<" N nucleons: "<<passiveNucleons<<endl;
	    cout<<"GetNormValue::l36"<<endl;
      
        if(targetZ > 90){
            passiveNucleons = TargetUtils::Get().GetTrackerNNucleons( 12., true, 850.0 );
	    if(anyTracker) passiveNucleons*=9;
	    cout<<"N nucleons: "<<passiveNucleons<<endl;
        }
        
    }

    else
    {
        assert( false && "Target nucleons can be all, proton or neutron only" );
    }
    
    if( passiveNucleons < 0 )
        assert( false && "Normalizations only known for Z = 1,6,26,82 and 0(scint)" );
    
    cout<<"The normalization factor for this analysis is "<<trackerAtomsC / passiveNucleons<<endl;
    
    return  trackerAtomsC/passiveNucleons;
}






//======================================
//! NukeCC XSec
//======================================
class NukeCCXSec : public XSec
{
    
public:
    
    NukeCCXSec( const char* name, int targetZ = 0, bool passesCuts = 0 ) :
    XSec( name ),
    m_nucleus( targetZ ), 
    m_cuts( passesCuts){}; 

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
    
    
    bool PassTrueDistToDivisionCut( ChainWrapper& chw, int entry, int targetID, int targetZ, double xySep /* = 25. */ )
    {
        if(debug) cout<<"PassTrueDistToDivisionCut::l174"<<endl;
        double true_target_dist_to_division = (double)chw.GetValue("truth_target_dist_to_division",entry);
	
        int true_targetID = (int)chw.GetValue("truth_targetID",entry);
        int true_targetZ  = (int)chw.GetValue("truth_targetZ",entry);
        int true_module   = (int)chw.GetValue("truth_vtx_module",entry);
        
        //        cout<<true_targetID<< "  target ID "<<targetID<<endl;
        if( targetID < 10 && true_targetID != targetID) return false;
        //        cout<<true_targetZ<< "  target ID "<<num_targetZ<<endl;
        if( targetID < 10 && true_targetZ != targetZ) return false;
	if(targetID >10 && anyTracker) {
	  if( true_module < 27 || true_module > 80 ) return false;
	}//NukeCC really does cut off at 80 for all tracker, but 81 for 94 it looks like. But I compared the two and they were the same?? I think there might be a number of targets issue here. (maybe result of the plane code number change faiza did?). No I'm pretty sure that target 94 is just only supposed to go up to 80. Its 1 module bigger than every other one. 
	else{
	  if( targetID == 14 && ( true_module < 27 || true_module > 32)) return false;
	  if( targetID == 24 && ( true_module < 33 || true_module > 38)) return false;
	  if( targetID == 34 && ( true_module < 39 || true_module > 44)) return false;
	  if( targetID == 44 && ( true_module < 45 || true_module > 50)) return false;
	  if( targetID == 54 && ( true_module < 51 || true_module > 56)) return false;
	  if( targetID == 64 && ( true_module < 57 || true_module > 62)) return false;
	  if( targetID == 74 && ( true_module < 63 || true_module > 68)) return false;
	  if( targetID == 84 && ( true_module < 69 || true_module > 74)) return false;
	  if( targetID == 94 && ( true_module < 75 || true_module > 80)) return false;//true_module > 81)) return false;
	}
        if(debug) cout<<"PassTrueDistToDivisionCut::l194 \t target_module: "<<true_module<<endl;
        //            only relevant for passive targets 1235
	if( targetID<10 && targetID!=4) if(debug) cout<<"true_dist_to_division: "<<true_target_dist_to_division<<"\t true_targetID: "<<true_targetID<<"\t true_targetZ: "<<true_targetZ<<endl;
        if( 0 < true_targetID && true_targetID < 10 && 4 != true_targetID)
            return ( xySep < true_target_dist_to_division );
        //        cout<<"fiducial event "<<targetID<<endl;
        
	if(debug) cout<<"PassTrueDistToDivisionCut::l200"<<endl;
        return true;
    }
    
    

    bool passesCuts( ChainWrapper& chw, int entry ) {
        
        if((int)chw.GetValue("mc_incoming", entry)!=-14) return false;
        if((int)chw.GetValue("mc_current", entry)!=1) return false;
        if(!InHex(chw,entry,850.0)) return false;
        if( (int)chw.GetValue("truth_targetZ",entry) != targetZ) return false;

        if(targetZ == 6){ // carbon
          
          if((int)chw.GetValue("truth_targetID",entry) == 3) {
            if(!PassTrueDistToDivisionCut( chw, entry, 3, targetZ, 25.0)) return false;
          }
          
          else{return false;}
        }

        if(targetZ == 26){ // iron
        
          if ((int)chw.GetValue("truth_targetID",entry) == 2){
            if(!PassTrueDistToDivisionCut( chw, entry, 2, targetZ, 25.0)) return false;
          }
          
          else if((int)chw.GetValue("truth_targetID",entry) == 3) {
            if(!PassTrueDistToDivisionCut( chw, entry, 3, targetZ, 25.0)) return false;
          }
          
          else if (((int)chw.GetValue("truth_targetID",entry) == 5)){
            if(!PassTrueDistToDivisionCut( chw, entry, 5, targetZ, 25.0)) return false;
          }
          
          else{return false;}
        }

        if(targetZ == 82){ //lead
        
          if ((int)chw.GetValue("truth_targetID",entry) == 2){
            if(!PassTrueDistToDivisionCut( chw, entry, 2, targetZ, 25.0)) return false;
          }
          
          else if((int)chw.GetValue("truth_targetID",entry) == 3) {
            if(!PassTrueDistToDivisionCut( chw, entry, 3, targetZ, 25.0)) return false;
          }
          
          else if((int)chw.GetValue("truth_targetID",entry) == 4){
            if(!PassTrueDistToDivisionCut( chw, entry, 4, targetZ, 25.0)) return false;
          }
          else if (((int)chw.GetValue("truth_targetID",entry) == 5)){
            if(!PassTrueDistToDivisionCut( chw, entry, 5, targetZ, 25.0)) return false;
          }
          
          else{return false;}
        }
        

        //if ((int)chw.GetValue("truth_targetID",entry) != 3) return false;
        //if( (int)chw.GetValue("truth_targetZ",entry) != 6) return false;
        //if(!PassTrueDistToDivisionCut( chw, entry, 3, m_nucleus, 25.0)) return false;
        if(!passesKinematics(chw, entry)) return false;
        
        return true;
        
    } //! end function passesCuts
    
    //! input parameters
    int  m_nucleus;
    bool m_cuts;
    
};

//======================

int runXSecLooper_NukeCC(const bool antinu, const double Emin, const double Emax, const std::vector<const char*> fileNames)
{
  string nuc_name = "";
  if(targetZ==6) nuc_name = "carbon";
  if(targetZ==26) nuc_name = "iron";
  if(targetZ==82) nuc_name = "lead";

  const char* fileName = Form("GENIEXSecExtract_CCInclusive_multi%s.root", nuc_name.c_str()) ;
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

    loop.setNuPDG(antinu ? -14: 14);

    loop.setFiducial(5970, 8450);
    //loop.setFiducial(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, 850.0);

    loop.setPlaylist(PlotUtils::FluxReweighter::minervame6A);

    //! Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
    loop.setNumUniv(0); //Tammy used 100
    
    //! store cross section names
    const char* channel = "CC";
    
    string processes[1] = { "inclusive"};
    
    vector<XSec::EVariable> vars;
    vars.push_back( XSec::kENu );
    vars.push_back( XSec::kxExp );
            
      //! loop over the cross section extraction
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
        }
        
        
        char* name = Form( "%s_%s_std", nuc_name.c_str(), varName.c_str());
        cout<<"making cross section for "<< targetZ<<endl;
        NukeCCXSec* xsec = new NukeCCXSec( name, targetZ, true );
        
        //! container for bins
        std::vector<double> bins;   
        GetBins( *var, bins, false, false );
        int num_bins = (int)bins.size()-1;
        //! set cross section data
        xsec->setBinEdges(num_bins,&bins.front());
        xsec->setNormalizationValue( GetNormValue(targetZ) );
        xsec->setVariable(*var);
        //if(*var !=XSec::kENu) {xsec->setIsFluxIntegrated(false);}
        //else 
        xsec->setIsFluxIntegrated(true);
        xsec ->setDimension(1); 
        xsec->setFluxIntLimits(0.,120.);
        xsec->setUniverses(0);
        xsec->setNormalizationType(XSec::kSelfNorm);
        //! add cross section
        loop.addXSec(xsec);
              
        //! end loop over cross sections
      }
    
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

int main( int argc, char *argv[] )
{
    cout << "Enter running the GENIE Xsection Extraction for the NukeCC analysis" << endl;

    cout << "running genie xsection extraction" << endl;

    TH1::AddDirectory(kFALSE); //Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).

    std::vector<const char*> fileNames(argv+1, argv+argc);

    bool antinu=true;

    int Emin=0; //GeV
    int Emax=120; //GeV 

    return runXSecLooper_NukeCC(antinu, Emin, Emax, fileNames);
    
    cout << "Exit running the GENIE XSection Extraction for the NukeCC analysis" << endl;
}
