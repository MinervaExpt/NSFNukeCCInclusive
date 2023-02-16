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
}
bool debug=false;
bool anyTracker=false;//true;

double GetNormValue( int targetID, int targetZ, int targetNucleon = 0 )
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
        
        if(targetID < 10 ){
              cout<<"GetNormValue::l34"<<endl;
            passiveNucleons = TargetUtils::Get().GetPassiveTargetNNucleons(targetID,targetZ, true, 850.);
	    cout<<"target "<<targetID<<" N nucleons: "<<passiveNucleons<<endl;
	    cout<<"GetNormValue::l36"<<endl;
        }
        if(targetID > 10){
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
    
    NukeCCXSec( const char* name, int target = 0, int targetZ = 0, bool isDIS = 0, bool InHex = 0, bool PassTrueDistToDivisionCut = 0 ) :
    XSec( name ),
    m_target( target ),
    m_nucleus( targetZ ),
    m_isDIS( isDIS),
    m_isinHex( InHex),
    m_truedisttoDiv(PassTrueDistToDivisionCut ){};
   
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
    
    
    bool PassTrueDistToDivisionCut( ChainWrapper& chw, int entry, int targetID, int num_targetZ, double xySep /* = 25. */ )
    {
        if(debug) cout<<"PassTrueDistToDivisionCut::l174"<<endl;
        double true_target_dist_to_division = (double)chw.GetValue("truth_target_dist_to_division",entry);
	
        int true_targetID = (int)chw.GetValue("truth_targetID",entry);
        int true_targetZ  = (int)chw.GetValue("truth_targetZ",entry);
        int true_module   = (int)chw.GetValue("truth_vtx_module",entry);
        
        //        cout<<true_targetID<< "  target ID "<<targetID<<endl;
        if( targetID < 10 && true_targetID != targetID) return false;
        //        cout<<true_targetZ<< "  target ID "<<num_targetZ<<endl;
        if( targetID < 10 && true_targetZ != num_targetZ) return false;
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
    
    

    virtual bool passesCuts( ChainWrapper& chw, int entry ) {
        
        if((int)chw.GetValue("mc_incoming", entry)!=-14) return false;
        if((int)chw.GetValue("mc_current", entry)!=1) return false;
        if(!InHex(chw,entry,850.0)) return false;
        if(!PassTrueDistToDivisionCut( chw, entry, m_target, m_nucleus, 25.0)) return false;
        if(!passesKinematics(chw, entry)) return false;
        
        return true;
        
    } //! end function passesCuts
    
    //! input parameters
    int m_target;
    int  m_nucleus;
    bool m_isDIS;
    bool m_isinHex;
    bool m_truedisttoDiv;
    
};

//======================

int runXSecLooper_NukeCC(const bool antinu, const double Emin, const double Emax, const std::vector<const char*> fileNames, int targetID, int targetZ, string outdir)
{
  const char* fileName = Form("%s/GENIEXSecExtract_CCInclusive_t%d_z%02d.root", outdir.c_str(), targetID, targetZ);
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

    
    vector<int> targetZs;
    //targetZs.push_back(6);
    //targetZs.push_back(26);
    targetZs.push_back(targetZ);
    
    vector<int> targetIDs;
    //targetIDs.push_back(1);
    targetIDs.push_back(targetID);
    //targetIDs.push_back(3);
    //targetIDs.push_back(4);
    //targetIDs.push_back(5);
    //targetIDs.push_back(14);
    //targetIDs.push_back(24);
/*   targetIDs.push_back(34);
    targetIDs.push_back(44);
    targetIDs.push_back(54);
    targetIDs.push_back(64);
    targetIDs.push_back(74);
    targetIDs.push_back(84);
    targetIDs.push_back(94);*/
    for(int i=0; i <targetIDs.size(); ++i){
        
        for(int j=0; j <targetZs.size(); ++j){
            
            int num_target = targetIDs[i];
            int num_targetZ = targetZs[j];
            
            if(targetIDs[i]!=3 && targetZs[j]==6)
                continue;
            if(targetIDs[i]==4 && targetZs[j]!=82)
                continue;
            if(targetIDs[i]>10 && targetZs[j]!=82)
                continue;
            
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
                string nuc_name = "";
                if(num_targetZ==6) nuc_name = "carbon";
                if(num_targetZ==26) nuc_name = "iron";
                if(num_targetZ==82 && num_target< 10) nuc_name = "lead";
                if(num_targetZ==82 && num_target> 10) nuc_name = "tracker";
                
                
                char* name = Form( "%s_%d_%s_std", nuc_name.c_str(), num_target, varName.c_str());
                cout<<"making cross section for "<<num_target<<"  "<<num_targetZ<<endl;
                NukeCCXSec* xsec = new NukeCCXSec( name, num_target, num_targetZ, true, true, true );
                
                //! container for bins
                std::vector<double> bins;   
                GetBins( *var, bins, false, false );
                int num_bins = (int)bins.size()-1;
                //! set cross section data
                xsec->setBinEdges(num_bins,&bins.front());
                xsec->setNormalizationValue( GetNormValue( num_target, num_targetZ ) );
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
                
            } //! end loop over cross sections
        }
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

  std::vector<const char*> fileNames(argv+4, argv+argc);

  bool antinu=true;

  int Emin=0; //GeV
  int Emax=120; //GeV 

  string outdir = argv[1];
  int targetID = atoi(argv[2]);
  int targetZ = atoi(argv[3]);

  cout << "Outdir " << outdir << std::endl;
  cout << "Target ID " << targetID << std::endl;
  cout << "Target Z " << targetZ << std::endl;

  return runXSecLooper_NukeCC(antinu, Emin, Emax, fileNames, targetID, targetZ, outdir);
  
  cout << "Exit running the GENIE XSection Extraction for the NukeCC analysis" << endl;
}
