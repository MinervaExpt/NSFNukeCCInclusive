// ========================================================================
// Base class for an un-systematically shifted (i.e. CV) universe.  Implement
// your own base class with the functions you need. I've implemented GetEnu(),
// GetMuonE() and GetHadronE() as examples: you'll have other variables you
// want.
//
// To add a systematic, inherit from this class and override whichever
// functions you need to. For a "vertical" error, this will mean overriding the
// GetWeight() function to modify the event weight. For a "lateral" error, this
// will mean overriding the function that calculates the quantity that is being
// shifted (muon energy, or hadronic energy or whatever).
//
// For examples of each of those two cases, see ./LateralSystematics.h and
// PlotUtils/GenieSystematics.h. For an example of how to put the whole thing
// together and actually *use* the classes, see the runEventLoop.C macro in
// this directory. `root -l -b load.C+ runEventLoop.C+`
// ========================================================================
#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

//#include "NukeCCsetbranchaddrs.h"
#include "CommonIncludes.h"
#include "PlotUtils/MinervaUniverse.h"
#include "PlotUtils/ChainWrapper.h"
#include <iostream>
//#include "CCQENuUtilsNSF.h"
//using namespace globalV;

//#include "NukeCCvars.h"
namespace NUKECC_ANA
{
class CVUniverse : public PlotUtils::MinervaUniverse {
  public:
    #include "PlotUtils/MuonFunctions.h" // GetMinosEfficiencyWeight
    #include "PlotUtils/TruthFunctions.h" //Getq3True
    #include "PlotUtils/WeightFunctions.h"
    // Constructor
    CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma=0);
    CVUniverse();

    // Destructor
    virtual ~CVUniverse();
    // All functions we write here should be 'virtual', so that the universes
    // that inherit from CVUniverse can override them.

    // ========================================================================
    // Get Weight
    //
    // We override the various weight getting functions herein in different
    // vertical systematic universe classes.
    // ========================================================================
  

    virtual double GetWeight() const;
    virtual double GetWeightEmu() const;
    virtual double GetWeightQ2() const;
    virtual double GetWeighty() const;
    virtual double GetTruthWeight() const;
    virtual double GetTruthWeightFlux() const;
    virtual Long64_t GetMaxEntries(); 
    virtual double GetMuonE() const { return GetVecElem("NukeCC_leptonE",3); }

    //      virtual void  Show(Long64_t entry = -1);


    ////// Dipak's study equivalent: ML vertex plane prob cut /////
    double virtual GetANNPlaneProb() const;

    /////////////////////////////////////////////////////////////////////////////////
    // MC initial momentum of the struck nucleons 
    // mc_initNucVec (px,py,pz,E) in MeV/c

    virtual double GetInNuclMomentumX() const;
    virtual double GetInNuclMomentumY() const;
    virtual double GetInNuclMomentumZ() const;
    virtual double GetInNuclMomentum() const;

    /////////////////////////////////////////////////////////////////////////////////
    // NEW for vertex position
    double virtual GetVertexXMy() const;
    double virtual GetVertexXTrueMy() const;
  
    double virtual GetVertexYMy() const;
    double virtual GetVertexYTrueMy() const;

    double virtual GetVertexZMy() const;
    double virtual GetVertexZTrueMy() const;
    ////////////////////////////////////////////////////////////////////////////////////


    virtual double calcq3(const double Q2,const double Enu,const double Emu	)const;
    virtual double calcq0(const double Enu,const double Emu	)const;
    //! Turn off branches that I don't use
    //virtual bool SetBranches( bool useWeights = true );
    virtual double GetMuonCurve()  const;
    virtual double GetFiducial()  const;
    virtual double GetTdead()  const;
    virtual double GetTargetID()  const;
    virtual double GetHelicity()  const;
    virtual double GetTrueHelicity()  const;
    virtual double GetRecoilEnergy() const;
    virtual double GetEhadGeV() const;
    virtual double GetMuonEGeV() const;
    virtual double GetMuonETrueGeV() const;
    virtual double GetThetamuTrueDeg()  const ;     
    //virtual double GetThetamu()  const ;
    virtual double GetMuonPt()  const ;
    virtual double GetMuonPz() const;
    virtual double GetlepPtTrue()  const ;
    virtual double GetlepPzTrue()  const ;


    virtual double GetQ2Reco()     const ;
    virtual double GetWReco()     const;
    virtual double GetxReco()     const;
    virtual double GetyReco()     const;
    virtual double GetplaneDNNReco() const;

    virtual double GetQ2RecoGeV()     const ;
    virtual double GetWRecoGeV()     const;
    virtual double GetEnu() const;  
    virtual double GetEnuGeV() const;  
    virtual double GetEnuTrueGeV() const;  
    // for true variables
    virtual double GetTruthNuPDG() const;
    virtual double GetCurrent() const;
    virtual double GetEhadTrue() const;
    virtual double GetEhadTrueGeV() const;
    virtual double GetThetamuTrue( )const;

    virtual double GetQ2IncTrue()     const ;
    virtual double GetWTrue()     const;
    virtual double GetyTrue()     const;
    virtual double GetxTrue()     const;
    virtual double GetplaneDNNTrue()   const; 
 
    virtual double GetQ2TrueGeV()     const ;
    virtual double GetWTrueGeV()     const;
    virtual double GetThetamuDeg()  const;
    //virtual double Getq3True()const; 
    //virtual double Getq0True()const; 
        
    virtual double calcRecoQ2(const double Enu,const double Emu,const double Thetamu)const;	
    virtual double calcWReco(const double Q2,const double Ehad) const;

    virtual double calcTrueQ2(const double EnuTrue,const double EmuTrue,const double ThetamuTrue)const;	
    virtual double calcWTrue(const double Q2True,const double EhadTrue) const;


    virtual double calcXTrue(const double Q2True,const double EnuTrue, const double EmuTrue)const;
    virtual double readXTrue() const;

    virtual double calcXReco(const double Q2,const double Enu, const double Emu)
    const ;


    virtual double calcYTrue(const double EnuTrue,const double EhadTrue)const;	

    virtual double calcYReco(const double Enu,const double Ehad)const;	

    virtual double Var( const std::string& var, bool useTrue = false );       
       
    virtual  int GetTargetFromSegment( int segment, int& vtx_module, int& vtx_plane ) const;
    
    virtual double GetVtxEnergy(  unsigned int shift = VTX_BLOB_SHIFT );  
	
    virtual double GetCCQERecoil( unsigned int shift = VTX_BLOB_SHIFT );
	
    virtual double GetMuonP() const;
		
    virtual std::vector<std::string> GetStdEnergyVars( ) const;
    virtual std::vector< std::pair<std::string,std::string> > GetStd2DVars( )const;
 
    virtual int GetAtomicNumber();
    //! Map Atomic Number to Its Atomic Mass
    virtual int AtomicNumberToMass( int targetZ ) const;
    //! Map Target to the Minimum Plane Number Bin
    virtual int GetTargetMinBin( int targetID ) const;
    //! Map Target to the Maximum Plane Number Bin
    virtual int GetTargetMaxBin( int targetID ) const;
    //! Map Target to its Plane Number
    virtual int GetTargetPlane( int targetID ) const;
    //! Map Plane to its respective Target ID
    virtual int GetPlaneTargetID( int plane ) const;
    //! Map Target to the Plane immediately upstream of it
    virtual int GetTargetUSPlane( int targetID ) const;
    //! Map Target to the Plane immediately downstream of it
    virtual int GetTargetDSPlane( int targetID ) const;
    //! Map Target to the Plane immediately downstream of it
    virtual double GetTargetZStart( int targetID ) const;
    //! Map Target to the Plane immediately downstream of it
    virtual double GetTargetZEnd( int targetID ) const;

    // ARACHNE EVENT DISPLAY
    // data
    virtual int GetRunN() const;
    virtual int GetSubRunN() const;
    virtual int GetGateN() const;
    virtual int GetSliceN() const;

    // MC
    virtual int GetMCRunN() const;
    virtual int GetMCSubRunN() const;
    virtual int GetMCGateN() const;
    virtual int GetMCSliceN() const;



  };
}

#endif
