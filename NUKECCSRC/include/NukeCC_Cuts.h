//#ifndef MNV_NUKECC_cxx
//#define MNV_NUKECC_cxx 1

#ifndef MNV_NUKECC_CUTS_h
#define MNV_NUKECC_CUTS_h 1

#include "GlobalIncludes.h"
#include "PlotUtils/MinervaUniverse.h"
#include "PlotUtils/TargetUtils.h"
#include "CVUniverse.h"
#include "CommonIncludes.h"
//#include "include/GlobalIncludes.h"


using namespace std;

namespace NUKECC_ANA{

class NukeCC_Cuts{

public:

//constructor
NukeCC_Cuts(); 

//Destructor

~NukeCC_Cuts();

static NukeCC_Cuts &Get();


    //!Set up a vector of bad events
      std::vector< pair< pair<int, int>, pair<int, int> > > InitBadEventVector();

      //===== Public member variables =====//
      //! Is this tree made from MC?
      bool isMC;

      
      //===== BEGIN ANALYSIS CUTS =====//


       
      bool PassMuCurveCut( CVUniverse *cv, double minCut, HelicityType::t_HelicityType h);
      bool PassMuCurveCut( CVUniverse *cv, HelicityType::t_HelicityType h);
      
       bool PassMuCoilCut( CVUniverse *cv );

      bool PassMuQualityCut(CVUniverse *cv, int qual = 2 );

      
        bool PassHelicityCut(CVUniverse *cv, HelicityType::t_HelicityType h );
      //! Do the true nu and true theta mu cuts
      bool PassTrueMuEnergyCut(CVUniverse *cv );


   //   virtual bool PassMuEnergyCut( );
       bool PassMuEnergyCut(CVUniverse *cv );

       bool PassThetaCut(CVUniverse *cv);
       bool PassTrueThetaCut(CVUniverse *cv);

       //void Show(CVUniverse *cv,Long64_t entry);
       bool PassGoodTrackingCut(CVUniverse *cv);
       bool PassZDistCut( CVUniverse *cv,double lowZ = 1001., double highZ = 1001. );

       bool PassDistToDivisionCut( CVUniverse *cv,double xySep = 25. );

      //! Is the true vertex far enough away form a division of target sections?
       bool PassTrueDistToDivisionCut( CVUniverse *cv,double xySep  = 25. );

       bool PassReco(CVUniverse *cv,HelicityType::t_HelicityType h);
       //Truth cuts
       bool PassTruth(CVUniverse *cv,HelicityType::t_HelicityType h);
       bool PassTrueCC(CVUniverse *cv,HelicityType::t_HelicityType h= HelicityType::kNeutrino);
       bool PassTrueFiducial(CVUniverse *cv);
       bool PassTrueInelasticCut(CVUniverse *cv);
       bool PassTrueHelicityCut(CVUniverse *cv, HelicityType::t_HelicityType h);
       bool InHexagonTrue(CVUniverse *cv, double apothem = 850. ); 


      //! Check that the ccqe recoil energy is greater than this value in MeV
       bool PassCCQERecoilCut(CVUniverse *cv, const double minE = MIN_CCQE_RECOIL_E );
      //! Check if passes the inelastic cuts using branch variables (when USE_INEL_CUT is true) 
       bool PassInelasticCut(CVUniverse* cv);
      //! Check if passes the inelastic cuts using given variables in MeV (when USE_INEL_CUT is true) 
       bool PassInelasticCut( double q2, double ccqeRecoil );

      //! Check if passes the DIS cuts using reco branch variables 
      bool PassDISCut(CVUniverse *cv);
      
       bool PassDISCut( double q2, double W );
      //! Check if passes the DIS cuts using given variables in MeV 
    //  bool PassDISCut( double q2, double W );

      

      //! Check if passes the DIS cuts using true variables in MeV 
       bool PassTrueDISCut(CVUniverse *cv );


       bool PassTrueDISCut( double q2, double W );


	bool PassLowWCut(CVUniverse* cv);

	bool PassLowWCut(double W );
       
	bool PassLowQ2Cut(CVUniverse* cv);
	bool PassLowQ2Cut( double q2, double W );

       
        bool PassLowQ2Trans(CVUniverse* cv);
         
        bool PassLowQ2Trans( double q2, double W );


        
         bool Passtrans(CVUniverse* cv);         
         bool Passtrans( double q2, double W );


	bool PassLowWCutTrue(CVUniverse* cv);

	bool PassLowWCutTrue(double W );
       
	bool PassLowQ2CutTrue(CVUniverse* cv);
	bool PassLowQ2CutTrue( double q2, double W );

       
        bool PassLowQ2TransTrue(CVUniverse* cv);
         
        bool PassLowQ2TransTrue( double q2, double W );


        
         bool PasstransTrue(CVUniverse* cv);         
         bool PasstransTrue( double q2, double W );


	 //True interaction type category cuts
	 //bool passTrueCC(CVUniverse *cv);
	 //bool passTrueCC(int pdg, int current);
	 bool passTrueCCQE(CVUniverse *cv);
	 bool passTrueCCQE(int pdg, int current, int type,bool charm);
	 bool passTrueCCRES(CVUniverse *cv);
	 bool passTrueCCRES(int pdg,int current,int type);
	 bool passTrueCCDIS(CVUniverse *cv);
	 bool passTrueCCDIS(int pdg,int current,int type);
	 bool passTrueMEC(CVUniverse *cv);
	 bool passTrueMEC(int pdg,int current,int type, bool charm);
	 bool passTrueCoh(CVUniverse *cv);
	 bool passTrueCoh(int pdg,int current,int type);
	 bool passTrueCCTrueDIS( CVUniverse* cv );
	 bool passTrueCCTrueDIS( int pdg, int current, int type, double Q2, double W);
	 bool passTrueCCTrueSIS( CVUniverse* cv );
	 bool passTrueCCTrueSIS( int pdg, int current, int type, double Q2, double W);

      //@}
      //===== END ANALYSIS CUTS =====//


      bool IsInMaterial(CVUniverse *cv,  int targID, int targZ, bool anyTrackerMod = false);
      //virtual bool IsInMaterial( const MatPair& mat, bool anyTrackerMod = true );


       bool IsInTrueMaterial(CVUniverse *cv, const int targID, const int targZ, bool anyTrackerMod = false );


      //@}
      //===== END ANALYSIS FUNCTIONS =====//

      //! coord to U (perdendicular to section divide in tgt 1,5)
       double GetU(double x, double y) const;
      //! coord transform to D (perdendicular to section divide in tgt 2)
      double GetD(double x, double y) const;
      //! coord transform to C (perdendicular to carbon divide in target 3)
       double GetC(double x, double y) const;

      
      //! Machine Learning Stuffs
       bool InHexagon(CVUniverse *cv, double apothem = 850. ); 
       
       bool IsInTargetSection( int targetID, int targetZ, double x, double y ) const;
       int GetTargetPlane( int targetID ) const;
};


}
#endif //MNV_NUKECC_cxx

