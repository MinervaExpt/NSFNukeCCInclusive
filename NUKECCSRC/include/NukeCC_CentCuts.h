#ifndef RECO_NUKECC_CENTCUTS_h
#define RECO_NUKECC_CENTCUTS_h 1

//#include "include/NukeCC_Cuts.h"
//#include "../include/NukeCCvars.h"
#include "PlotUtils/TargetUtils.h"
#include "include/CVUniverse.h"
//#include "PlotUtils/Cut.h"
//#include "include/GlobalIncludes.h" 
#include "PlotUtils/Cutter.h"
#include "include/CommonIncludes.h"
//#include "../include/NukeCCvars.h"
//#include "CCQENuUtilsNSF.h"
//#include "include/NukeUtils.h"
//#include "Acceptance/TAcceptanceTable.h"
#include <PlotUtils/MnvNormalization.h>
#include <PlotUtils/NuclModUtils.h>
#include <PlotUtils/FluxReweighter.h>
#include <PlotUtils/FluxReweighterWithWiggleFit.h>
//#include "include/CondorInput.h"
#include <PlotUtils/MnvNuclearModelWeight.h>
//#include <PlotUtils/MnvNormalizerME.h>
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
//#include "NukeCCsetbranchaddrs.h"

#include "TFileCollection.h"


#include "TVector3.h"

using namespace std;

using namespace NUKECC_ANA;

namespace reco
{	

bool IsInTargetSection( int targetID=1, int targetZ=1, double x=1.0, double y=1.0 );
int  GetTargetPlane( int targetID=1);
       double GetC(double x=1.0, double y=1.0); 
   double GetD(double x=1.0, double y=1.0); 
    double GetU(double x=1.0, double y=1.0);
 template <class UNIVERSE, class EVENT = detail::empty>    
 class PassMuCurveCut: public Cut<UNIVERSE, EVENT>
 {
	public:
		PassMuCurveCut():Cut<UNIVERSE, EVENT>(" Muon Curve"){}
	private:
		bool checkCut(const UNIVERSE& univ,  EVENT& /*evt*/) const override
		{
		    return univ.GetMuonCurve() <= 0 ;
		}
 };


 template <class UNIVERSE, class EVENT = detail::empty>    
 class PassMuCoilCut:public Cut<UNIVERSE, EVENT>
{
    public:
          PassMuCoilCut():Cut<UNIVERSE, EVENT>("Mu CoilCut"){}
    private:
          bool checkCut(const UNIVERSE& univ,  EVENT& /*evt*/) const override   
    {
    const double coilXPos = 1219.0;
    const double coilYPos = 393.0;
    const double minos_x = univ.GetDouble("NukeCC_minos_trk_end_x") + coilXPos;
    const double minos_y = univ.GetDouble("NukeCC_minos_trk_end_y") + coilYPos;
    double minosR = sqrt(pow(minos_x,2) + pow(minos_y,2) );
    
    return (minosR > MINOS_COIL_RADIUS && minosR < MAX_MINOS_RADIUS );
     }
};


 template <class UNIVERSE, class EVENT = detail::empty>    
 class PassMuQualityCut:public Cut<UNIVERSE, EVENT> 
{
   public:
     PassMuQualityCut(int qual=2):Cut<UNIVERSE, EVENT>("Mu QualityCut"),m_qual(qual){}
   private: 
          bool checkCut(const UNIVERSE& univ,  EVENT& /*evt*/) const override   
    {   
   //! If the minos track quality is unknown or not set, return false
    if( univ.GetDouble("NukeCC_minos_trk_quality") <= 0 )
        return false;
    //! If the minos track quality is equal to or better than requested, return true
    return ( univ.GetDouble("NukeCC_minos_trk_quality") <= m_qual );
    }
    int m_qual;
};


 template <class UNIVERSE, class EVENT = detail::empty>    
 class PassMuEnergyCut:public Cut<UNIVERSE, EVENT> 
{
   public:
   PassMuEnergyCut():Cut<UNIVERSE, EVENT>("Emu"){}
   private:
  
          bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override   
      {
    return ( MIN_RECO_E_MU < univ.GetEmu() && univ.GetEmu() < MAX_RECO_E_MU );
      }
};


 template <class UNIVERSE, class EVENT = detail::empty>    
 class PassThetaCut:public Cut<UNIVERSE, EVENT> 
{
   public:
   PassThetaCut():Cut<UNIVERSE, EVENT>("Thetamu"){}
   private:
  
          bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override   
      {
  return ( 0. <= univ.GetThetamu()*rad_to_deg && univ.GetThetamu()*rad_to_deg < MAX_RECO_THETA_MU );
      }
};


  template <class UNIVERSE, class EVENT = detail::empty>    
  class PassZDistCut:public Cut<UNIVERSE, EVENT>    
{  
  public:
    PassZDistCut():Cut<UNIVERSE, EVENT>("Pass ZDistCut"){}
    private:
          bool checkCut(const UNIVERSE& univ,  EVENT& /*evt*/) const override  
     { 
    double minZ = 4290.; // z center of mod -5 plane 1 is 4293.04mm, while plane 2 is 4313.68mm
    double maxZ = 6000.;
    if ( univ.GetInt("NukeCC_targetID") == 0 && ( minZ <= univ.GetVecElem("NukeCC_vtx",2) && univ.GetVecElem("NukeCC_vtx",2) <= maxZ ) ){
       return true;
    }
    //passive target ZDist cut is made with nPlanes cut in framework
    if( 1 <= univ.GetInt("NukeCC_targetID") && univ.GetInt("NukeCC_targetID") <= 5 )
    {
        // if we only want fitted vtx events with vertex really in target
        //if( IsMultiTrack() && 1.0E-4 < fabs(CVUniverse_target_zDist) )
        //  return false;
        return true;
    }
   
    //all events in the tracker pass for faux targets
    if( FIRST_TRACKER_MOD <= univ.GetInt("NukeCC_vtx_module") && univ.GetInt("NukeCC_vtx_module") <= LAST_TRACKER_MOD )
        return true;
    
    return false;
}
};

  template <class UNIVERSE, class EVENT = detail::empty>
  class PassDistToDivisionCut:public Cut<UNIVERSE, EVENT>
{
  public:
       PassDistToDivisionCut(double xySep  = 25. ):Cut<UNIVERSE, EVENT>("Pass DistToDivision"),m_xySep(xySep){}
  private:
     bool checkCut(const UNIVERSE& univ,  EVENT& /*evt*/) const override  
  {
    //only relevant for passive targets 1235
    if( 0 < univ.GetInt("NukeCC_targetID") && univ.GetInt("NukeCC_targetID") < 10 && 4 != univ.GetInt("NukeCC_targetID") )
        return ( m_xySep < univ.GetDouble("NukeCC_target_dist_to_division" ));
    
    return true;
   }
const double m_xySep;

};
  
 
  template <class UNIVERSE, class EVENT = detail::empty>
  class IsInMaterial:public Cut<UNIVERSE, EVENT>
{
   public:
       IsInMaterial(int i_targetID,  int i_targetZ, bool anyTrackerMod  = false):Cut<UNIVERSE, EVENT>("IsIn Material"), m_targetID(i_targetID),m_targetZ(i_targetZ),m_anyTrackerMod(anyTrackerMod){}
   private:
        bool checkCut( const UNIVERSE& univ,   EVENT& /*evt*/)const override
       {
 	return IsInMaterial_b(univ, m_targetID, m_targetZ, m_anyTrackerMod);
}


bool IsInMaterial_b(const UNIVERSE& univ,int m_targetID,  int m_targetZ, bool m_anyTrackerMod  = false ) const {
    if( m_targetID < 0)
    {
        // if targetID < 0, then we want any event in the nuclear target region mods -5-26
         double z = univ.GetVecElem("NukeCC_vtx",2);
        double minZ = 4290.; // z center of mod -5 plane 1 is 4293.04mm, while plane 2 is 4313.68mm
        double maxZ = 6000.;
        if ( !( minZ <= z && z <= maxZ ) )
            return false;
        
        // -targetID is the reference target
        // THIS DOESN'T WORK FOR LOCAL PLASTIC SIDEBAND
        // Instead check if the x,y position is in target material
        /*int refTarg = -i_targetID;
        if( i_targetZ > 0 && CVUniverse_ref_targZ[refTarg-1] != i_targetZ )
          return false;*/

        if( m_targetZ > 0 && ! IsInTargetSection(m_targetID, m_targetZ, univ.GetVecElem("NukeCC_vtx",0),univ.GetVecElem("NukeCC_vtx",1) ) )
          return false;

        return true;
    }
    
    if(m_targetID < 10 )
    {
       //THIS IS A SPECIAL CASE FOR PLASTIC BACKGROUND
      if(univ.GetInt("NukeCC_targetID") == 0 && univ.GetDouble("NukeCC_vtx_module") < 27 )
      {
	if( m_targetZ > 0 && IsInTargetSection( m_targetID, m_targetZ, univ.GetVecElem("NukeCC_vtx",0),univ.GetVecElem("NukeCC_vtx",1) ) )
	  return true;
      }
      //cout << "targetID < 10" << endl;
      // If targetID < 10, then you are looking for a passive target event.
      // Require that the event has the same targetID and targetZ.
      if( univ.GetInt("NukeCC_targetID") == m_targetID )
      {
	if( m_targetZ > 0 )
	  return univ.GetInt("NukeCC_targetZ") == m_targetZ;
	else
	  return true;
      }//targetID< 10 looking as passive target. Event doesn't have right targetID
    }
    else if( m_targetID < 100 )
    {
        // If 10 < targetID < 100, then we are looking for an event in a plastic reference target.
        // Say targetID = AT, then the event must be in the Ath active target group and the reference target is T.
        int refTarg = m_targetID % 10;
        
        // The starting module of the 4-module reference target is 6*A + 21
        int refModStart = 6*((m_targetID - refTarg )/10) + 21;
        
        int refTargID = refTarg*10000 + 6*1000 + refModStart;
	//return true;
        return IsInMaterial_b(univ,refTargID, m_targetZ, m_anyTrackerMod );

    }
    else
    {
        int refTarg = (m_targetID - m_targetID % 10000 ) / 10000;
        if( m_targetZ > 0 && univ.GetVecElem("NukeCC_ref_targZ",(refTarg-1)) != m_targetZ )
            return false;
        
        int refModStart = m_targetID % 1000;
        int refNMod = ( ( m_targetID - refModStart ) / 1000 % 10 );
        int firstMod = refModStart;
        int lastMod  = refModStart + refNMod - 1;
        
        if( m_anyTrackerMod )
        {
            firstMod = FIRST_TRACKER_MOD;
            lastMod  = LAST_TRACKER_MOD;
        }
        
        // OK if the vertex module is within the range specified
        if( firstMod <= univ.GetDouble("NukeCC_vtx_module") && univ.GetDouble("NukeCC_vtx_module") <= lastMod )
            return true;
    }
    return false;

}

int m_targetID,m_targetZ;
bool m_anyTrackerMod;
};
//bool  reco::IsInTargetSection( int targetID=1, int targetZ=1, double x=1.0, double y=1.0 );

 //template <class UNIVERSE, class EVENT = detail::empty>
bool  IsInTargetSection( int targetID, int targetZ, double x, double y )
{
    // refTarg is targetID for passives
    //     ... targetID%10 for shorthand faux
    int refTarg = targetID % 10;
    
    //     ... and this for scint chunks
    if( targetID > 100 )
        refTarg = (targetID - targetID % 10000 ) / 10000;
    
    //everyone's a winner in target 4!
    if( 4 == refTarg ) return true;
    
    if( 1 == refTarg || 5 == refTarg )
    {
        double u = GetU( x, y );
        //iron is u < 205mm
        if( 26 == targetZ ) return ( u < 205 );
        //lead is 205mm < u
        if( 82 == targetZ ) return ( 205 < u );
    }
    if( 2 == refTarg )
    {
        double d = GetD( x, y );
        //iron is d < 205mm
        if( 26 == targetZ ) return ( d < 205 );
        //lead is 205mm < d
        if( 82 == targetZ ) return ( 205 < d );
    }
    if( 3 == refTarg )
    {
        double c = GetC(x,y);
        //carbon is 0mm < c
        if( 6 == targetZ ) return ( 0 < c );
        //iron is c < 0mm and x < 0
        if( 26 == targetZ ) return ( c < 0 && x < 0 );
        //lead is c < 0mm and 0 < x
        if( 82 == targetZ ) return ( c < 0 && 0 < x );
    }
    
    Error( "IsInTargetSection", Form("Given invalid targetID %d and targetZ %d", targetID, targetZ ) );
    return false;

}

    double GetU(double x, double y) 
{
    return -x*cos(PI/6) + y*sin(PI/6); // perp to divide in target 1/5.  up and to the left.
}


   double GetD(double x, double y)  
{
    return x*cos(PI/6) + y*sin(PI/6); //points to MINERvA racks (negative x), perp to divide in target 2. down and to the right.
}


       double GetC(double x, double y)  
{
    return x*sin(PI/6) + y*cos(PI/6); // points to MINERvA racks (negative x), perp to carbon of target 3
}


        int  GetTargetPlane( int targetID )
{ 
    if( targetID == 1 ) return 9;
    else if( targetID == 2 ) return 19;
    else if( targetID == 3 ) return 30;
    else if( targetID == 4 ) return 49;
    else if( targetID == 5 ) return 55;
    else{
        std::cerr << "Bad targetID choice in NukeCC_Cuts::GetTargetMaxBin! Try again. Hint: we only have 5 passive target" << std::endl;
    }
    return -999;
}

 template <class UNIVERSE, class EVENT = detail::empty>
class PassHelicityCut: public Cut<UNIVERSE, EVENT>
{
   public:
     PassHelicityCut(HelicityType::t_HelicityType h ):Cut<UNIVERSE, EVENT>("PassHelicity"),m_h(h){ }
   private:
        bool checkCut(const UNIVERSE& univ,  EVENT& /*evt*/) const override 
{   
   //int helicity = univ.GetHelicity();
   if(m_h==0) return true;
   else 
     return HelicityType::t_HelicityType(univ.GetHelicity()) == m_h; 
}
	HelicityType::t_HelicityType m_h;
};


  template <class UNIVERSE, class EVENT = detail::empty> 
  class Deadtime: public Cut<UNIVERSE, EVENT> 
  {
    public:
       Deadtime():Cut<UNIVERSE, EVENT>("tDead Cut") {}
     private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override       {
       if(univ.GetTdead()<=1)
       //if(1>univ.GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj"))
       return true;
       }
};
       


  template <class UNIVERSE, class EVENT = detail::empty> 
  class FiducialCut: public Cut<UNIVERSE, EVENT> 
  {
    public:
       FiducialCut():Cut<UNIVERSE, EVENT>("Fiducial Cut") {}
     private:
      bool checkCut(const UNIVERSE& univ,  EVENT& /*evt*/) const override       {
      // if(univ.GetFiducial())
       return univ.GetFiducial();
       }
};
       

 template <class UNIVERSE, class EVENT = detail::empty>
 class IsDIS: public Cut<UNIVERSE, EVENT>
 {
    public:
      // Constructor
      IsDIS(const double Q2Min = 1, const double WMin = 2): Cut<UNIVERSE, EVENT>("Is DIS"), m_Q2Min(Q2Min), m_WMin(WMin){ }
    private:
      // THE cut function
      bool checkCut(const UNIVERSE& univ,  EVENT& /*evt*/) const override
      {
        // Call a CVUniverse member function to make the cut
        return  univ.GetQ2RecoGeV() <= m_Q2Min && univ.GetWRecoGeV() >= m_WMin;
      }
//    private:
      const double m_Q2Min; //GeV^2
      const double m_WMin; //GeV/c
  };


  template <class UNIVERSE, class EVENT = detail::empty> 
  class TargetIDCut: public Cut<UNIVERSE, EVENT> 
  {
    public:
       TargetIDCut(int targetID):Cut<UNIVERSE, EVENT>("TargetID Cut"),m_targetID(targetID) {}
     private:
      bool checkCut(const UNIVERSE& univ,  EVENT& /*evt*/) const override       {
      // if(univ.GetFiducial())
       if (m_targetID<10 && univ.GetTargetID() != m_targetID)
       return true;
       }
int m_targetID;
};



  template <class UNIVERSE, class EVENT = detail::empty> 
  class ANNProbCut: public Cut<UNIVERSE, EVENT> 
  {
    public:
       ANNProbCut(const double MIN_PROB_PLANE_CUT = 0.2):Cut<UNIVERSE, EVENT>("ANNProb Cut"),MIN_PROB(MIN_PROB_PLANE_CUT) {}
     private:
      bool checkCut(const UNIVERSE& univ,  EVENT& /*evt*/) const override       {
      // if(univ.GetFiducial())
       if (univ.GetVecElem("ANN_plane_probs",0) > MIN_PROB)
       return true;
       }
const double MIN_PROB;
};


                  //   if(targetID<10 && universe->GetInt("NukeCC_targetID") != targetID) continue;
                    // if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;	   

template <class UNIVERSE, class EVENT = detail::empty>
cuts_t<UNIVERSE, EVENT> getCCInclusiveCuts( int i_targetID, int i_targetZ, HelicityType::t_HelicityType h,bool anyTrackerMod)
  {
    cuts_t<UNIVERSE, EVENT> inclusive_cuts;
      inclusive_cuts.emplace_back(new PassHelicityCut<UNIVERSE, EVENT>(h));
      inclusive_cuts.emplace_back(new FiducialCut<UNIVERSE, EVENT>());
      inclusive_cuts.emplace_back(new PassZDistCut<UNIVERSE, EVENT>());
       inclusive_cuts.emplace_back(new PassDistToDivisionCut<UNIVERSE, EVENT>());
       inclusive_cuts.emplace_back(new Deadtime<UNIVERSE, EVENT>());
       inclusive_cuts.emplace_back(new PassMuCurveCut<UNIVERSE, EVENT>());
       inclusive_cuts.emplace_back(new PassMuCoilCut<UNIVERSE, EVENT>());
      
     inclusive_cuts.emplace_back(new IsInMaterial<UNIVERSE, EVENT>(i_targetID, i_targetZ, false));
      inclusive_cuts.emplace_back(new TargetIDCut<UNIVERSE, EVENT>(i_targetID));
      inclusive_cuts.emplace_back(new ANNProbCut<UNIVERSE, EVENT>());
      inclusive_cuts.emplace_back(new PassMuEnergyCut<UNIVERSE, EVENT>());
      inclusive_cuts.emplace_back(new PassThetaCut<UNIVERSE, EVENT>());
   // inclusive_cuts.emplace_back(new IsDIS<UNIVERSE, EVENT>());

  return inclusive_cuts;
  }

}
#endif 
