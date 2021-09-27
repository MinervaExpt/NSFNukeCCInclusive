#ifndef RECO_NUKECC_TRUTHCUTS_h
#define RECO_NUKECC_TRUTHCUTS_h 1

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

namespace truth
{	

template <class UNIVERSE>
  class IsNeutrino: public PlotUtils::SignalConstraint<UNIVERSE>	
  {
    public:
      IsNeutrino(): PlotUtils::SignalConstraint<UNIVERSE>("IsNeutrino")
      {
      }
    private:	
      bool checkConstraint(const UNIVERSE& univ) const override
      {	
        return univ.GetTruthNuPDG() == 14;
      }	
  };


template <class UNIVERSE>
  class IsCC: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      IsCC(): PlotUtils::SignalConstraint<UNIVERSE>("IsCC")
      {
      }
    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetCurrent() == 1;
      }
  };


 template <class UNIVERSE>
 class PassTrueMuEnergyCut: public PlotUtils::SignalConstraint<UNIVERSE>
 {
    public:
 PassTrueMuEnergyCut():PlotUtils::SignalConstraint<UNIVERSE>("TrueMuEnergyCut") {}
    private:

  bool checkConstraint(const UNIVERSE& univ) const override  
   //data always passes
{
 //   if( ! isMC )
   //     return true;
    
    //same as reco limits
    return ( MIN_RECO_E_MU <  univ.GetElepTrue() &&  univ.GetElepTrue() < MAX_RECO_E_MU );
}

};


 template <class UNIVERSE>
 class PassTrueThetaCut: public PlotUtils::SignalConstraint<UNIVERSE>
 {
    public:
 PassTrueThetaCut():PlotUtils::SignalConstraint<UNIVERSE>("True Theta Cut") {}
    private:

  bool checkConstraint(const UNIVERSE& univ) const override  
{
  return ( 0. <= univ.GetThetalepTrue()*rad_to_deg && univ.GetThetalepTrue()*rad_to_deg < MAX_RECO_THETA_MU );
}
};


 template <class UNIVERSE>
 class IsInTrueMaterial: public PlotUtils::SignalConstraint<UNIVERSE>
{
  public:
       IsInTrueMaterial(int i_targetID,  int i_targetZ, bool anyTrackerMod  = false):PlotUtils::SignalConstraint<UNIVERSE>("IsInTrue Material"), m_targetID(i_targetID),m_targetZ(i_targetZ),m_anyTrackerMod(anyTrackerMod){}
   private:
 
        bool checkConstraint( const UNIVERSE& univ)const override
{ 

 	return IsInTrueMaterial_b(univ, m_targetID, m_targetZ, m_anyTrackerMod);
} 

bool IsInTrueMaterial_b(const UNIVERSE& univ,int m_targetID,  int m_targetZ, bool m_anyTrackerMod  = false ) const {
//int truth_targetID=cv->GetInt("truth_targetID");
//  int truth_vtx_module=cv->GetInt("truth_vtx_module");
//  int truth_targetZ=cv->GetInt("truth_targetZ");
    if(m_targetID < 10 )
    {
        // If targetID < 10, then you are looking for a passive target event.
        // Require that the event has the same targetID and targetZ.
        if(univ.GetInt("truth_targetID") == m_targetID )
        {
            if( m_targetZ > 0 )
                return univ.GetInt("truth_targetZ") == m_targetZ;
            else
                return true;
        }
    }
    else if( m_targetID < 100 )
    {
        // If 10 < targetID < 100, then we are looking for an event in a plastic reference target.
        // Say targetID = AT, then the event must be in the Ath active target group and the reference target is T.
        int refTarg = m_targetID % 10;
        
        // The starting module of the 4-module reference target is 6*A + 21
        int refModStart = 6*((m_targetID - refTarg )/10) + 21;
        
        int refTargID = refTarg*10000 + 6*1000 + refModStart;
        return IsInTrueMaterial_b( univ,refTargID, m_targetZ, m_anyTrackerMod );
    }
    else
    {
        int refTarg = (m_targetID - m_targetID % 10000 ) / 10000;
        //if( m_targetZ > 0 && univ.GetVecElem("truth_ref_targZ",[refTarg-1]) != m_targetZ )
          
        if( m_targetZ > 0 && univ.GetVecElem("truth_ref_targZ",(refTarg-1)) != m_targetZ )
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
        if( firstMod <= univ.GetInt("truth_vtx_module") && univ.GetInt("truth_vtx_module") <= lastMod )
            return true;
        
    }
    return false;
}


int m_targetID,m_targetZ;
bool m_anyTrackerMod;
};

 template <class UNIVERSE>
 class IsTrueDIS: public PlotUtils::SignalConstraint<UNIVERSE> 
 {
    public:
      // Constructor
      IsTrueDIS(const double Q2Min = 1, const double WMin = 2): PlotUtils::SignalConstraint<UNIVERSE>("IsTrue DIS"), m_Q2Min(Q2Min), m_WMin(WMin){ }
    private:
      // THE cut function
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        // Call a CVUniverse member function to make the cut
        return  univ.GetQ2TrueGeV() <= m_Q2Min && univ.GetWTrueGeV()>= m_WMin;
      }
//    private:
      const double m_Q2Min; //GeV^2
      const double m_WMin; //GeV/c
  };



 template <class UNIVERSE>
 class PassLowQ2CutTrue: public PlotUtils::SignalConstraint<UNIVERSE> 
 {
    public:
      // Constructor
      PassLowQ2CutTrue(const double Q2Min = 1, const double WMin = 2): PlotUtils::SignalConstraint<UNIVERSE>("LowQ2CutTrue"), m_Q2Min(Q2Min), m_WMin(WMin){ }
    private:
      // THE cut function
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        // Call a CVUniverse member function to make the cut
        return  univ.GetQ2TrueGeV() < m_Q2Min && univ.GetWTrueGeV()>= m_WMin;
      }
//    private:
      const double m_Q2Min; //GeV^2
      const double m_WMin; //GeV/c
  };



// lowQ2Trans True

 template <class UNIVERSE>
 class PassLowQ2TransTrue: public PlotUtils::SignalConstraint<UNIVERSE> 
 {
    public:
      // Constructor
      PassLowQ2TransTrue():PlotUtils::SignalConstraint<UNIVERSE>("LowQ2Trans True"){ }
    private:
      // THE cut function
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        // Call a CVUniverse member function to make the cut
        return  univ.GetQ2TrueGeV() < 1.0 && 1.5 <= univ.GetWTrueGeV() && univ.GetWTrueGeV() < 2.0;
      }
//    private:
     // const double m_Q2Min; //GeV^2
     // const double m_WMin; //GeV/c
  };


 template <class UNIVERSE>
 class PasstransTrue: public PlotUtils::SignalConstraint<UNIVERSE> 
 {
    public:
      // Constructor
      PasstransTrue():PlotUtils::SignalConstraint<UNIVERSE>("trans True"){ }
    private:
      // THE cut function
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        // Call a CVUniverse member function to make the cut
        return  univ.GetQ2TrueGeV() >= 1.0 && 1.5 <= univ.GetWTrueGeV() && univ.GetWTrueGeV() < 2.0;
      }
//    private:
     // const double m_Q2Min; //GeV^2
     // const double m_WMin; //GeV/c
  };


/*
bool NukeCC_Cuts::PassTruth(CVUniverse* cv, HelicityType::t_HelicityType h)
{
    if( ! PassTrueCC(cv,h) ) return false; //true nu_mu (bar) CC event?
    if( ! PassTrueFiducial(cv ) ) return false; //true fiducial volume cut
    
    return true;
}


bool NukeCC_Cuts::PassTrueCC(CVUniverse* cv, HelicityType::t_HelicityType h )
{
    if( 1 != cv->GetInt("mc_current") ) return false;
    if( ! PassTrueHelicityCut(cv, h ) ) return false;
    
    return true;
}

bool NukeCC_Cuts::PassTrueFiducial(CVUniverse* cv )
{
    if( ! InHexagonTrue(cv, 850.) ) return false;
    if( ! PassTrueDistToDivisionCut(cv ) ) return false;
    
    return true;
}

bool NukeCC_Cuts::PassTrueInelasticCut(CVUniverse* cv)
{
    if( 1 == cv->GetInt("mc_intType") ) return false;
    if( 2 == cv->GetInt("mc_intType") ) return false;
	 
     return true;
}
*/
/*
template <class UNIVERSE>
class PassTrueHelicityCut: public PlotUtils::SignalConstraint<UNIVERSE>
{
   public:
     PassTrueHelicityCut(HelicityType::t_HelicityType h ):PlotUtils::SignalConstraint<UNIVERSE>("TrueHelicity"),m_h(h){ }
   private:
        bool checkCut(const UNIVERSE& univ) const override 
{   
   int helicity = univ.GetTrueHelicity();
   if(m_h==0) return true;
    else if (HelicityType::kNeutrino == m_h)
     return helicity ==14;
  else if (HelicityType::kAntiNeutrino == m_h)
     return helicity == -14;
  else{
     Warning( "PassTrueHelicityCut", "Helicity type unknown, should not have compiled.  Returning false..." );
     return false;
}
	HelicityType::t_HelicityType m_h;
};

 template <class UNIVERSE>
 class InHexagonTrue: public PlotUtils::SignalConstraint<UNIVERSE>
 {
   public:
      // Constructor
      InHexagonTrue(const double apothem): PlotUtils::SignalConstraint<UNIVERSE>(std::string("Apothem")+std::to_string(apothem)),fapothem { }
    private:
    {
      // THE cut function
      bool checkCut(const UNIVERSE& univ) const override
   
    //vector<double> newVertex = GetMLVertex();
   // double x = cv->GetVecElem("mc_vtx",0);
   // double y = cv->GetVecElem("mc_vtx",1);
    
    if( pow(univ.GetVecElem("mc_vtx",0),2) + pow(univ.GetVecElem("mc_vtx",1),2) < pow(fapothem,2) )
        return true;
    
    //Hexagon is symmetric about its x and y
    x = fabs(univ.GetVecElem("mc_vtx",0));
    y = fabs(univ.GetVecElem("mc_vtx",1));
    
    double lenOfSide = fapothem * ( 2 / sqrt(3) );
    
    if( x > fapothem )
        return false;
    
    if( y < lenOfSide/2.0 )
        return true;
    
    double slope = (lenOfSide / 2.0) / fapothem;
    if( y < lenOfSide - x*slope )
        return true;
    return false;
}
const double fApothem;

};
*/
/*
bool NukeCC_Cuts::passTrueCCQE( CVUniverse* cv ){
  return passTrueCCQE( cv->GetInt("mc_incoming"), cv->GetInt("mc_current"), cv->GetInt("mc_intType"), cv->GetInt("mc_charm") );
}

bool NukeCC_Cuts::passTrueCCQE( int pdg, int current, int type, bool charm ){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 1 && !charm );
  else return ( pdg == -14 && current == 1 && type == 1 && !charm );
}

bool NukeCC_Cuts::passTrueCCRES(CVUniverse* cv){
  return passTrueCCRES(cv->GetInt("mc_incoming"),cv->GetInt("mc_current"),cv->GetInt("mc_intType"));

}

bool NukeCC_Cuts::passTrueCCRES( int pdg, int current, int type ){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 2 );
  else  return ( pdg == -14 && current == 1 && type == 2 );
}


bool NukeCC_Cuts::passTrueMEC(CVUniverse* cv){
  return passTrueMEC(cv->GetInt("mc_incoming"),cv->GetInt("mc_current"),cv->GetInt("mc_intType"),cv->GetBool("mc_charm"));
}


bool NukeCC_Cuts::passTrueMEC( int pdg, int current, int type, bool charm ){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 8 && !charm );
  else  return ( pdg == -14 && current == 1 && type == 8 && !charm );
}

bool NukeCC_Cuts::passTrueCoh(CVUniverse* cv){
  return passTrueCoh(cv->GetInt("mc_incoming"),cv->GetInt("mc_current"),cv->GetBool("mc_intType"));

}

bool NukeCC_Cuts::passTrueCoh( int pdg, int current, int type ){ 
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 4 ); 
  else return ( pdg == -14 && current == 1 && type == 4 ); //antinu 
}   

bool NukeCC_Cuts::passTrueCCDIS(CVUniverse *cv){
  return passTrueCCDIS(cv->GetInt("mc_incoming"),cv->GetInt("mc_current"),cv->GetInt("mc_intType"));

}

bool NukeCC_Cuts::passTrueCCDIS( int pdg, int current, int type ){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 3 );
  else  return ( pdg == -14 && current == 1 && type == 3 );
}


bool NukeCC_Cuts::passTrueCCTrueDIS( CVUniverse* cv ){
  return passTrueCCTrueDIS( cv->GetInt("mc_incoming"), cv->GetInt("mc_current"), cv->GetInt("mc_intType"), cv->GetQ2True(), cv->GetWTrue() );
}

bool NukeCC_Cuts::passTrueCCTrueDIS( int pdg, int current, int type, double Q2, double W){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 3 && Q2>=1000 && W >=2000);
  else  return ( pdg == -14 && current == 1 && type==3 && Q2>=1000 && W >=2000);
}

bool NukeCC_Cuts::passTrueCCTrueSIS( CVUniverse* cv ){
  return passTrueCCTrueSIS( cv->GetInt("mc_incoming"), cv->GetInt("mc_current"), cv->GetInt("mc_intType"), cv->GetQ2True(), cv->GetWTrue() );
}

bool NukeCC_Cuts::passTrueCCTrueSIS( int pdg, int current, int type, double Q2, double W){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 3 && !(Q2>=1000 && W >=2000));
  else  return ( pdg == -14 && current == 1 && type == 3 && !(Q2>=1000 && W >=2000));
}
*/
template <class UNIVERSE>
  PlotUtils::constraints_t<UNIVERSE> GetCCInclusive2DPhaseSpace(int i_targetID, int i_targetZ, HelicityType::t_HelicityType h,bool anyTrackerMod)
  {
    PlotUtils::constraints_t<UNIVERSE> signalDef;
    
     signalDef.emplace_back(new IsInTrueMaterial<UNIVERSE>(i_targetID, i_targetZ, anyTrackerMod));
//     signalDef.emplace_back(new InHexagonTrue<UNIVERSE>(850.));
     signalDef.emplace_back(new PassTrueThetaCut<UNIVERSE>());
     signalDef.emplace_back(new PassTrueMuEnergyCut<UNIVERSE>());

   // signalDef.emplace_back(new ZRange<UNIVERSE>("Tracker", 5980, 8422));
    return signalDef;
  }

  template <class UNIVERSE>
  PlotUtils::constraints_t<UNIVERSE> GetCCInclusive2DSignal()
  {
    PlotUtils::constraints_t<UNIVERSE> signalDef;
    signalDef.emplace_back(new IsNeutrino<UNIVERSE>());
    signalDef.emplace_back(new IsCC<UNIVERSE>());
    return signalDef;
  }


}
#endif 
