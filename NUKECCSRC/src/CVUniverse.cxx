#ifndef CVUNIVERSE_cxx
#define CVUNIVERSE_cxx 1

#include "../include/CVUniverse.h" 
#include "../include/Cuts.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
//#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
#include "PlotUtils/PlotUtilsPhysicalConstants.h"

#include "PlotUtils/MinervaUniverse.h"
#include "PlotUtils/ChainWrapper.h"
//using namespace globalV;
using namespace NUKECC_ANA;
//again the constructor.....
CVUniverse::CVUniverse(PlotUtils::ChainWrapper *chw,double nsigma):PlotUtils::MinervaUniverse(chw,nsigma){
  //just go with 2 universe world
 // recoShifter = new MnvRecoShifter(neutrinoMode,2);
  

}

//destructor....
CVUniverse::~CVUniverse(){
 // delete recoShifter;
}

Long64_t CVUniverse::GetMaxEntries()
{
   
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    //! If MAX_ENTRIES environmental is set to something positive, use it.  Otherwise use them all.
    if( getenv("MAX_ENTRIES") )
    {
        int maxEntries = atoi(getenv( "MAX_ENTRIES") );
        if( maxEntries > 0  )
        return maxEntries;
    }
    return fChain->GetEntries();
}

////// Dipak's study equivalent: ML vertex plane prob cut /////
double CVUniverse::GetANNPlaneProb() const {
    double prob = GetVecElem("ANN_plane_probs",0); // y position in cm...
    return prob;
}
/////////////////////////////////////////////////////////////////////////////////
// MC initial momentum of the struck nucleons 
// mc_initNucVec (px,py,pz,E) in MeV/c

double CVUniverse::GetInNuclMomentumX() const {
    double initNuclpX = GetVecElem("mc_initNucVec",0); // in MeV ..
    return initNuclpX;
}

double CVUniverse::GetInNuclMomentumY() const {
    double initNuclpY = GetVecElem("mc_initNucVec",1); // in MeV ..
    return initNuclpY;
}

double CVUniverse::GetInNuclMomentumZ() const {
    double initNuclpZ = GetVecElem("mc_initNucVec",2); // in MeV ..
    return initNuclpZ;
}

double CVUniverse::GetInNuclMomentum() const {
    double initNuclpX = GetVecElem("mc_initNucVec",0); // in MeV ..
    double initNuclpY = GetVecElem("mc_initNucVec",1); // in MeV ..
    double initNuclpZ = GetVecElem("mc_initNucVec",2); // in MeV ..

    double InNuclMomentum = sqrt(initNuclpX*initNuclpX + initNuclpY*initNuclpY + initNuclpZ*initNuclpZ);

    return InNuclMomentum;
}


/////////////////////////////////////////////////////////////////////////////////
// NEW for vertex position
// Reco

double CVUniverse::GetVertexXMy() const {
    double vertex_y = GetVecElem("ANN_vtx",0)/10; // y position in cm...
    return vertex_y;
}

double CVUniverse::GetVertexXTrueMy() const { 
    return GetVecElem("mc_vtx",0)/10; 
}


double CVUniverse::GetVertexYMy() const {
    double vertex_y = GetVecElem("ANN_vtx",1)/10; // y position
    return vertex_y;
}

double CVUniverse::GetVertexYTrueMy() const { 
    return GetVecElem("mc_vtx",1)/10; 
}


double CVUniverse::GetVertexZMy() const {
    double vertex = GetVecElem("ANN_vtx", 2)/10; // in this current ntuple, Z position in mm not cm, hence divide
    return vertex;
}

double CVUniverse::GetVertexZTrueMy() const { 
    return GetVecElem("mc_vtx", 2)/10; 
}


// ORIGINAL
//calling branches and variables
double CVUniverse::GetMuonCurve()  const { return 1/GetDouble((GetAnaToolName() + "_minos_trk_eqp_qp").c_str()); }
double CVUniverse::GetHelicity()  const { return GetInt((GetAnaToolName() + "_nuHelicity").c_str()); }
double CVUniverse::GetTrueHelicity()  const { return GetInt("mc_incomig"); }
double CVUniverse::GetFiducial()  const { return GetInt((GetAnaToolName() + "_ANN_in_fiducial_area").c_str()); }
double CVUniverse::GetTdead()  const { return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj"); }
double CVUniverse::GetTargetID()   const{return GetInt((GetAnaToolName() + "_ANN_targetID").c_str());}

//double CVUniverse::GetRecoilEnergy()  const { return GetDouble((GetAnaToolName() + "_ANN_recoil_E").c_str()); }
// ======================================================================================================================================
//  ENERGY FUNCTIONS
// ======================================================================================================================================

// Recoil energy [GeV]
// ===================
/* Calorimetric recoil energy directly from anatuples (i.e. not corrected by multiplicative factors nor splines).
* This method of defining recoil energy is a bit ambiguous, but basically can be explained this in a way that
* we can satisfy the functions that are in [MAT-MINERvA/calculators/RecoilEnergyFunctions.h]:
* 
* 'GetCalRecoilEnergy()' returns the calorimetric recoil of the event, in my case it's all calorimetry energy since
* I don't use any non-calorimetric recoil reconstruction. For this, I use the "Recoil_Ecalo" branch of my anatool.
* 
* 'GetNonCalRecoilEnergy()' returns the non calorimetric recoil. In my case it will return zero
* since I don't perform any hadron energy reconstruction via non-calorimetric methods (like proton/pion tracks).
* 
* GetRecoilEnergy() comes from RecoilEnergyFunctions.h in MAT-MINERvA/calculators by adding the above up
*/


double CVUniverse::GetCalRecoilEnergy() const {
    return GetDouble((GetAnaToolName() + "_ANN_recoil_E").c_str());  // [MeV]
    // usual inclusive recoil energy should just be calorimetric energy
    // instead of Gonzalo's return GetDouble("Recoil_Ecalo") / 1000.0;  // [GeV]
}

double CVUniverse::GetNonCalRecoilEnergy() const {
    return 0.0;
}


//double CVUniverse::GetThetamu()       const {return GetDouble("muon_theta");}
double CVUniverse::GetEnu()           const {return GetEmu()+ GetRecoilEnergy();} 
double CVUniverse::GetQ2Reco()        const { return calcRecoQ2(GetEnu(),  GetEmu(), GetThetamu()); }
double CVUniverse::GetWReco()         const{return calcWReco(GetQ2Reco(), GetRecoilEnergy());}
double CVUniverse::GetxReco()         const{return calcXReco(GetQ2Reco(), GetEnu(),GetEmu());}
double CVUniverse::GetyReco()         const{return calcYReco(GetEnu(), GetRecoilEnergy());}

// Plastic bkg
double CVUniverse::GetplaneDNNReco()  const{
           int ANN_vtx_module, ANN_vtx_plane;
           int Segment = GetVecElem("ANN_segments",0);
           int targetID = GetTargetFromSegment( Segment, ANN_vtx_module, ANN_vtx_plane );
           double ANN_VTX_module = (double)ANN_vtx_module; 
           double ANN_VTX_plane = (double)ANN_vtx_plane;
           double test=ANN_VTX_module*2.+ANN_VTX_plane+10.; 
//           return  ANN_VTX_module*2.+ANN_VTX_plane+10.;}
           return test;}

// Function to calculate the difference between the original reconstructed plane and the plane that one gets using the segment with the second highest probability
double CVUniverse::GetDiffPlaneReco() const{
    // reconstructed plane
    double reco_plane_0 = GetplaneDNNReco();

    // get reconstructed plane using the second segment in the vector (with the second highest probability)
    int ANN_vtx_module, ANN_vtx_plane;
    int Segment = GetVecElem("ANN_segments",1);
    int targetID = GetTargetFromSegment( Segment, ANN_vtx_module, ANN_vtx_plane );
    double ANN_VTX_module = (double)ANN_vtx_module; 
    double ANN_VTX_plane = (double)ANN_vtx_plane;
    double reco_plane_1 = ANN_VTX_module*2.+ANN_VTX_plane+10.; 

    double dif = reco_plane_0 - reco_plane_1;

    return dif;
}

double CVUniverse::GetDiffPlaneRecoTrue() const{
    // reconstructed plane (most probable)
    double reco_plane_0 = GetplaneDNNReco();
    // get true plane
    double true_plane = GetplaneDNNTrue();
    

    double dif = reco_plane_0 - true_plane;

    return dif;
}

// for true variables
double CVUniverse::GetTruthNuPDG() const{return GetInt("mc_incoming");}
double CVUniverse::GetCurrent()    const{return GetInt("mc_current");} 
double CVUniverse::GetEhadTrue()    const{return GetEnuTrue() - GetElepTrue();}
double CVUniverse::GetThetamuTrue( )const{return GetThetalepTrue();}

double CVUniverse::GetQ2IncTrue()      const{return calcTrueQ2(GetEnuTrue(),  GetElepTrue(), GetThetamuTrue()); }
double CVUniverse::GetWTrue()       const{return calcWTrue(GetQ2IncTrue(), GetEhadTrue());}

double CVUniverse::GetxTrue()       const{return calcXTrue(GetQ2IncTrue(), GetEnuTrue(),GetElepTrue());}
double CVUniverse::GetyTrue()       const{return calcYTrue(GetEnuTrue(), GetEhadTrue());}
//double CVUniverse::Getq3True()      const{return calcq3(GetQ2True(),  GetEnuTrue(), GetElepTrue());}
//double CVUniverse::Getq0True()      const{return calcq0(GetEnuTrue(), GetElepTrue());}
double CVUniverse::GetThetamuDeg()  const{return GetThetamu()*rad_to_deg;}
double CVUniverse::GetThetamuTrueDeg()  const{return GetThetamuTrue()*rad_to_deg;}
double CVUniverse::GetMuonP()       const{return GetPmu()/1000;} // in GeV
double CVUniverse::GetlepPTrue()       const{return GetPlepTrue()/1000;} // in GeV
double CVUniverse::GetplaneDNNTrue()         const{return  static_cast<double>(GetInt("truth_vtx_module")*2+GetInt("truth_vtx_plane")+10);}

// in GeV
double CVUniverse::GetMuonEGeV()           const {return GetEmu()*mev_to_gev;}
double CVUniverse::GetMuonETrueGeV()           const {return GetElepTrue()*mev_to_gev;}
double CVUniverse::GetEnuGeV()           const {return (GetEmu()+ GetRecoilEnergy())*mev_to_gev;}

double CVUniverse::GetEnuTrueGeV()           const {return (GetEnuTrue())*mev_to_gev;}
double CVUniverse::GetEhadGeV()  const { return GetRecoilEnergy()*mev_to_gev; }

double CVUniverse::GetEhadTrueGeV()    const{return (GetEnuTrue() - GetElepTrue())*mev_to_gev;}

double CVUniverse::GetQ2RecoGeV()        const { return (calcRecoQ2(GetEnu(),  GetEmu(), GetThetamu()))*mev_to_gev*mev_to_gev; }
double CVUniverse::GetWRecoGeV()         const{return (calcWReco(GetQ2Reco(), GetRecoilEnergy()))*mev_to_gev;}

double CVUniverse::GetQ2TrueGeV()      const{return (calcTrueQ2(GetEnuTrue(),  GetElepTrue(), GetThetamuTrue()))*mev_to_gev*mev_to_gev; }
double CVUniverse::GetWTrueGeV()       const{return (calcWTrue(GetQ2IncTrue(), GetEhadTrue()))*mev_to_gev;}

double CVUniverse::GetMuonPt()const{ // in GeV
   return GetPmu()/1000. * sin(GetThetamu());
}

double CVUniverse::GetMuonPz()const{ // in GeV
  return  GetPmu()/1000. * cos(GetThetamu());
}

double CVUniverse::GetlepPtTrue()const{ // in GeV
   return GetPlepTrue()/1000*sin(GetThetalepTrue());
}

double CVUniverse::GetlepPzTrue()  const{ // in GeV
    return GetPlepTrue()/1000*cos(GetThetalepTrue()); 
} // Gev/c 

double CVUniverse::calcTrueQ2(const double EnuTrue,const double EmuTrue,const double ThetamuTrue)	
const {
double Q2 = 4*EnuTrue*EmuTrue*pow( sin( ThetamuTrue/2), 2.);// - pow(M_muonMeV,2); // EnuTruue in MeV
//std::cout<<" True Enu = "<< EnuTrue <<"\t"<<" True Emu = "<< EmuTrue <<"\t"<<" True ThetaMu = "<< ThetamuTrue <<"\t"<<" True Q2 = "<<Q2<<std::endl;
return Q2; 
}
  
double CVUniverse::calcWTrue(const double Q2True,const double EhadTrue)
const {
double nuclMass = 1.0;
if( NEUTRON_PDG  == GetInt("mc_targetNucleon") )
    nuclMass = M_neutron;
else if( PROTON_PDG == GetInt("mc_targetNucleon") )
    nuclMass = M_proton;
else
    nuclMass = M_nucleon;

double W2 = ( pow(nuclMass, 2) +  2. * ( EhadTrue ) * nuclMass - Q2True);
//if(W2>0){std::cout<<Q2True<<"\t"<<sqrt(W2)<<std::endl;}
//std::cout<<" True Ehad = "<< EhadTrue <<"\t"<<" Nucleon Mass =  "<<nuclMass<<"\t"<<" True Q2 = "<<Q2True<<std::endl;
W2 = W2 > 0 ? sqrt(W2) : 0.0;
return W2;
}  

double CVUniverse::calcRecoQ2(const double Enu,const double Emu,const double Thetamu)	
const {
    double Q2Reco = 4*Enu*Emu*pow( sin( Thetamu/2), 2.);// - pow(M_muonMeV,2);
    return Q2Reco; 
}
  
double CVUniverse::calcWReco(const double Q2,const double Ehad)
const {
double nuclMass = M_nucleon;
 //if( NEUTRON_PDG  == GetInt("mc_targetNucleon") )
                //nuclMass = M_neutron;
           // else if( PROTON_PDG == GetInt("mc_targetNucleon") )
                //nuclMass = M_proton;

double W2Reco = ( pow(nuclMass, 2) +  2. * ( Ehad ) * nuclMass - Q2);
W2Reco = W2Reco > 0 ? sqrt(W2Reco) : 0.0;
return W2Reco;
//return GetDouble("NukeCC_W");
}  


double CVUniverse::calcXTrue(const double Q2True,const double EnuTrue, const double EmuTrue)
const {
double nuclMass = 1.0;
if( NEUTRON_PDG  == GetInt("mc_targetNucleon") )
    nuclMass = M_neutron;
else if( PROTON_PDG == GetInt("mc_targetNucleon") )
    nuclMass = M_proton;
else
    nuclMass = M_nucleon;
double x = Q2True / (2. * ( EnuTrue -  EmuTrue) * nuclMass);
return x;
}

double CVUniverse::readXTrue() const{ // read true GENIE value from Monte Carlo
    double x = GetDouble("mc_Bjorkenx");
    return x;
}

double CVUniverse::calcXReco(const double Q2,const double Enu, const double Emu)
const {
double nuclMass = M_nucleon;
 //if( NEUTRON_PDG  == GetInt("mc_targetNucleon") )
                //nuclMass = M_neutron;
            //else if( PROTON_PDG == GetInt("mc_targetNucleon") )
                //nuclMass = M_proton;
double x = Q2 / (2. * ( Enu -  Emu) * nuclMass);
return x;
}  


double CVUniverse::calcYTrue(const double EnuTrue,const double EhadTrue)	
const {
return EhadTrue/EnuTrue; 
}
	
double CVUniverse::calcYReco(const double Enu,const double Ehad)	
const {
return Ehad/Enu; 
}

double CVUniverse::calcq3(const double 	Q2,const double Enu,const double Emu	 )const
{
	return sqrt(Q2 + pow(Enu - Emu,2.0));
}

double CVUniverse::calcq0(const double Enu,const double Emu)const
{
	return Enu-Emu;
}

double CVUniverse::GetWeightEmu() const {
   //const bool do_warping = true;
   double wgt_flux=1., wgt_2p2h=1.;
   double wgt_rpa=1.,   wgt_nrp=1.,  wgt_lowq2=1.;
   double wgt_genie=1., wgt_mueff=1.;
   //double wgt_anisodd=1.;
   double wgt_emu=1.;
 
   // genie
   wgt_genie = GetGenieWeight();
   double Enu  = GetDouble("mc_incomingE")*1e-3;
   int nu_type = GetInt("mc_incoming");
   // flux
   wgt_flux = GetFluxAndCVWeight(Enu, nu_type);
   ///For Emu Iron
   //wgt_emu = 1.06991e-07 * GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 1.02357e-05* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()+ 0.000378062* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 0.00713933*GetMuonETrueGeV()*GetMuonETrueGeV() + 0.0670028*GetMuonETrueGeV() + 0.919948;
   //For Iron 2D
    //wgt_emu = 1.34*(5.4546e-05 * GetMuonETrueGeV()* GetMuonETrueGeV() + 9.65024e-06 * GetEhadTrueGeV() * GetEhadTrueGeV() + 0.760717);
    //For Tracker 2D
    wgt_emu = 1.34*(4.21486e-05 * GetMuonETrueGeV()* GetMuonETrueGeV() + 7.51935e-05 * GetEhadTrueGeV() * GetEhadTrueGeV() + 0.734119);
   //For Emu Lead
   //wgt_emu = 2.22052e-07 * GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 2.39719e-05* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()+ 0.000979981* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 0.0188648*GetMuonETrueGeV()*GetMuonETrueGeV() + 0.163442*GetMuonETrueGeV() + 0.668094;
   //For Emu Tracker
   //wgt_emu = 1.06991e-07 * GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 1.02357e-05* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()+ 0.000378062* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 0.00713933*GetMuonETrueGeV()*GetMuonETrueGeV() + 0.0670028*GetMuonETrueGeV() + 0.919948;
   // 2p2h
   wgt_2p2h = GetLowRecoil2p2hWeight();
   // rpa
   wgt_rpa = GetRPAWeight();
   //MINOS Efficiency
   wgt_mueff = GetMinosEfficiencyWeight(); 
   
   return wgt_flux*wgt_genie*wgt_rpa*wgt_nrp*wgt_lowq2*wgt_mueff*wgt_2p2h*wgt_emu;
   //cout<< "\t NonResWeight: "<<GetGenieWeight()<<"\t 2p2hWeight: "<<GetLowRecoil2p2hWeight()<<"\t RPAWeight: "<<GetRPAWeight()<<"\t MINOS Eff: "<<GetMinosEfficiencyWeight()<<"\t\t  wgt_emu"<<wgt_emu<<endl; 

}

double CVUniverse::GetWeightQ2() const {
   //const bool do_warping = true;
   double wgt_flux=1., wgt_2p2h=1.;
   double wgt_rpa=1.,   wgt_nrp=1.  /*, wgt_lowq2=1.*/;
   double wgt_genie=1., wgt_mueff=1.;
   double wgt_q2=1.;
 
   // genie
   wgt_genie = GetGenieWeight();
   double Enu  = GetDouble("mc_incomingE")*1e-3;
   int nu_type = GetInt("mc_incoming");
   // flux
   wgt_flux = GetFluxAndCVWeight(Enu, nu_type);
   //For Iron
   wgt_q2 = 0.00143089*GetQ2TrueGeV()+1.16667;
   //For Lead
   //wgt_q2 = 0.0013294*GetQ2TrueGeV()+1.11518;
   //For Tracker
   //wgt_q2 = -0.000423594*GetQ2TrueGeV()+1.12566;
   // 2p2h
   wgt_2p2h = GetLowRecoil2p2hWeight();
   // rpa
   wgt_rpa = GetRPAWeight();
   //MINOS Efficiency
   wgt_mueff = GetMinosEfficiencyWeight(); 
   
   return wgt_flux*wgt_genie*wgt_rpa*wgt_nrp*wgt_q2*wgt_mueff*wgt_2p2h;

}

double CVUniverse::GetWeighty() const {
   //const bool do_warping = true;
   double wgt_flux=1., wgt_2p2h=1.;
   double wgt_rpa=1.,   wgt_nrp=1.,  wgt_lowq2=1.;
   double wgt_genie=1., wgt_mueff=1.;
   //double wgt_anisodd=1.;
   double wgt_y=1.;
 
   // genie
   wgt_genie = GetGenieWeight();
   double Enu  = GetDouble("mc_incomingE")*1e-3;
   int nu_type = GetInt("mc_incoming");
   // flux
   wgt_flux = GetFluxAndCVWeight(Enu, nu_type);
   //For 2D Iron with xy weight function
   //wgt_y = 1.3*(0.00117*GetxTrue()*GetxTrue() + 0.00724* GetyTrue()*GetyTrue() + 0.00146);
   //For 2D Tracker with xy weight function
   wgt_y = 1.3*(-9.76e-05*GetxTrue()*GetxTrue() + 0.0180* GetyTrue()*GetyTrue() + 0.00362);
   //For Iron
   //wgt_y = 1.24836*GetyTrue()*GetyTrue()* GetyTrue() - 1.64309* GetyTrue()*GetyTrue() + 0.49131*GetyTrue() + 1.16857;
   //For Lead
   //wgt_y = 1.95476*GetyTrue()*GetyTrue()* GetyTrue() - 2.80867* GetyTrue()*GetyTrue() + 1.19048*GetyTrue() + 0.979472;
   //For Lead
   //wgt_y = 3.02641*GetyTrue()*GetyTrue()* GetyTrue() - 4.12468* GetyTrue()*GetyTrue() + 1.5317*GetyTrue() + 0.988857;
   // 2p2h
   wgt_2p2h = GetLowRecoil2p2hWeight();
   // rpa
   wgt_rpa = GetRPAWeight();
   //MINOS Efficiency
   wgt_mueff = GetMinosEfficiencyWeight(); 
   
   return wgt_flux*wgt_genie*wgt_rpa*wgt_nrp*wgt_lowq2*wgt_mueff*wgt_2p2h*wgt_y;

}

// Functions for coherent pion reweight
// Anezka 24/04/2023 using Mehreen's def

double CVUniverse::GetTrueHighEpi() const {
    int nFSpi = GetInt("mc_nFSPart");
    double pionE = -1.0;
    double pionKE = -1.0;
    for (int i = 0; i < nFSpi; i++){
        int pdg = GetVecElem("mc_FSPartPDG",i);
        if(pdg != -211) continue;
        double energy = GetVecElem("mc_FSPartE", i);
        double mass = 139.569;
        double tpi = energy - mass;
        if (tpi >= pionKE){
            pionKE = tpi;
            pionE = energy;
        }
    } 
    //std::cout << "Printing energy of pion " << pionE << std::endl;
    return pionE; // MeV
}

double CVUniverse::thetaWRTBeam(double x, double y, double z) const{
    double pyp = -1.0 *sin( MinervaUnits::numi_beam_angle_rad)*z + cos( MinervaUnits::numi_beam_angle_rad )*y;
    double pzp = cos( MinervaUnits::numi_beam_angle_rad )*z + sin( MinervaUnits::numi_beam_angle_rad )*y;
    double denom2 = pow(x,2) + pow(pyp,2) + pow(pzp,2);
    if( 0. == denom2 ) return -9999.;
    else return acos(pzp / sqrt(denom2) );
}

double CVUniverse::GetTrueAngleHighTpi() const {
    int nFSpi = GetInt("mc_nFSPart");
    double angle = -9999.; //WRTbeam and in degrees
    double pionKE = 0.0;
    int idk = -9999;
    for (int i = 0; i < nFSpi; i++){
        int pdg = GetVecElem("mc_FSPartPDG",i);
        if(pdg != -211) continue;
        double energy = GetVecElem("mc_FSPartE", i);
        double mass = 139.569;
        double tpi = energy - mass;
        if (tpi >= pionKE) {
            pionKE = tpi;
            TVector3 pimomentumvec(GetVecElem("mc_FSPartPx", i), GetVecElem("mc_FSPartPz", i),GetVecElem("mc_FSPartPz", i));
            double deg_wrtb = thetaWRTBeam(pimomentumvec.X(), pimomentumvec.Y(), pimomentumvec.Z()); //rad
        
            angle = deg_wrtb; //*180./M_PI;
        }
    }
    //Making sure angle is only between 0 and pi
    if (angle < 0.0) angle = -1.0*angle;
    if (angle > PI) angle = 2.0*PI - angle; 
    return angle*180./PI; // Degrees
}

double CVUniverse::GetWeight() const {
    //const bool do_warping = true;
    double wgt_flux=1., wgt_2p2h=1.;
    double wgt_rpa=1.,   wgt_nrp=1.,  wgt_lowq2=1.;
    double wgt_genie=1., wgt_mueff=1.;
    double wgt_target_mass = 1.;
    double wgt_geant = 1.;
    double wgt_coh = 1.;
    double wgt_fsi = 1.;

    double Enu  = GetDouble("mc_incomingE")*1e-3;
    int nu_type = GetInt("mc_incoming");
 
    // genie
    wgt_genie = GetGenieWeight();
    // flux
    wgt_flux = GetFluxAndCVWeight(Enu, nu_type);
    // 2p2h
    wgt_2p2h = GetLowRecoil2p2hWeight();
    // rpa
    wgt_rpa = GetRPAWeight();
    // MINOS efficiency
    wgt_mueff = GetMinosEfficiencyWeight(); 
    // target mass systermatics
    wgt_target_mass = GetTargetMassWeight();
    // GEANT hadron weight
    wgt_geant = GetGeantHadronWeight();
    //Aaron's LowQ2 weights
    wgt_lowq2 = GetLowQ2PiWeight("MENU1PI");
    // coherent pion reweight
    if(GetInt("mc_intType") == 4){
       double angle = GetTrueAngleHighTpi();//*180./M_PI; //this is now in degrees
       double KE = GetTrueHighEpi()/1000.; // This is supposed to be the energy of the pion!!!!
       if (KE < 0){
        wgt_coh = 1.;
       }
       else {
	    wgt_coh = GetCoherentPiWeight(angle, KE); //Inputs are in Degrees and GeV
        //std::cout << "Printing COHerent weight for COH event for angle: " << angle << " degrees and KE : " << KE << " GeV. And weight: " << wgt_coh << std::endl;
       }
    }
    //GENIE FSI bug fix
    //wgt_fsi = GetFSIWeight(0);

    return wgt_genie*wgt_flux*wgt_2p2h*wgt_rpa*wgt_mueff*wgt_target_mass*wgt_geant*wgt_lowq2*wgt_coh*wgt_fsi;

}

double CVUniverse::GetTruthWeight()const{
     
    double wgt_flux=1., wgt_2p2h=1.;
    double wgt_rpa=1.,   wgt_nrp=1.,  wgt_lowq2=1.;
    double wgt_genie=1., wgt_mueff=1.;
    double wgt_target_mass = 1;
    double wgt_coh = 1.;
    double wgt_fsi = 1.;

 
    //There need to be flags added to this to turn on and off different tunes. Same goes for GetWeight().  -- ANF 2020-3-18
    double Enu  = GetDouble("mc_incomingE")*1e-3;
    int nu_type = GetInt("mc_incoming");
    // genie
    wgt_genie = GetGenieWeight();
    // flux
    wgt_flux = GetFluxAndCVWeight(Enu, nu_type);
    // 2p2h
    wgt_2p2h = GetLowRecoil2p2hWeight();
    // rpa
    wgt_rpa = GetRPAWeight();
    // target mass systematics
    wgt_target_mass = GetTargetMassWeight();
    // Aaron's lowQ2 weight
    wgt_lowq2 = GetLowQ2PiWeight("MENU1PI");
    // coherent pion reweight
    if(GetInt("mc_intType") == 4){
       double angle = GetTrueAngleHighTpi();//*180./M_PI; //this is now in degrees
       double KE = GetTrueHighEpi()/1000.; // This is supposed to be the energy of the pion!!!!
       if (KE < 0){
        wgt_coh = 1.;
       }
       else {
	    wgt_coh = GetCoherentPiWeight(angle, KE); //Inputs are in Degrees and GeV
        //std::cout << "Printing COHerent weight for COH event for angle: " << angle << " degrees and KE : " << KE << " GeV. And weight: " << weight << std::endl;
       }
    }
    //GENIE FSI bug fix
    //wgt_fsi = GetFSIWeight(0);

   // Note: truth weight has no GEANT hadron and MINOS weight      
   return wgt_genie*wgt_flux*wgt_2p2h*wgt_rpa*wgt_target_mass*wgt_lowq2*wgt_coh*wgt_fsi;
   }

double CVUniverse::GetTruthWeightFlux()const{
     
    double wgt_flux=1.;
    double wgt_genie=1.;

    double Enu  = GetDouble("mc_incomingE")*1e-3;
    int nu_type = GetInt("mc_incoming");
    
    // flux
    wgt_flux = GetFluxAndCVWeight(Enu, nu_type);
    // genie
    wgt_genie = GetGenieWeight();

   return wgt_flux*wgt_genie;
}

bool CVUniverse::IsInHexagon(double apothem /*= 850. */ ) const
{
   double x = GetVecElem((GetAnaToolName()+"_vtx").c_str(),0);
   double y = GetVecElem((GetAnaToolName()+"_vtx").c_str(),1);
   double lenOfSide = apothem*(2/sqrt(3)); 
   double slope     = (lenOfSide/2.0)/apothem;
   double xp        = fabs(x);
   double yp        = fabs(y);
   
   if( (xp*xp + yp*yp) < apothem*apothem )             return true;
   else if( xp <= apothem && yp*yp < lenOfSide/2.0 )   return true; 
   else if( xp <= apothem && yp < lenOfSide-xp*slope ) return true;

   return false;
}

bool CVUniverse::IsInHexagonTrue(double apothem /*= 850. */ ) const
{
   double x = GetVecElem("mc_vtx",0);
   double y = GetVecElem("mc_vtx",1);
   double lenOfSide = apothem*(2/sqrt(3)); 
   double slope     = (lenOfSide/2.0)/apothem;
   double xp        = fabs(x);
   double yp        = fabs(y);
   
   if( (xp*xp + yp*yp) < apothem*apothem )             return true;
   else if( xp <= apothem && yp*yp < lenOfSide/2.0 )   return true; 
   else if( xp <= apothem && yp < lenOfSide-xp*slope ) return true;

   return false;
}

double CVUniverse::GetDaisyPetal() const
{
  double x = GetVecElem((GetAnaToolName()+"_vtx").c_str(),0);
  double y = GetVecElem((GetAnaToolName()+"_vtx").c_str(),1);
  double apothem = 850.;
  if( !IsInHexagon( apothem ) ) return -1;
  double angle = atan2(y,x) / ( 2 * TMath::Pi() );
  if ( angle < 0 ) angle += 1;
  int out = floor(12*angle);
  return out;
} 


double CVUniverse::GetDaisyPetalTrue( ) const
{ 
  double x = GetVecElem("mc_vtx",0);
  double y = GetVecElem("mc_vtx",1);
  double apothem = 850.;
  if( !IsInHexagonTrue(apothem ) ) return -1;
  double angle = atan2(y,x) / ( 2 * TMath::Pi() );
  if ( angle < 0 ) angle += 1;
  int out = floor(12*angle);
  return out;
} 

// Plastic background

double CVUniverse::Var( const std::string& name, bool useTrue )
{
    if( name == "Emu" )
    {
        if( useTrue )      return GetElepTrue() * mev_to_gev;
        else               return GetEmu() * mev_to_gev ;
    }
    if ( name == "Enu" || name == "E" )
    {
        if( useTrue )      return GetEnuTrue() * mev_to_gev;
        else               return GetEnu() * mev_to_gev;
    }
    if ( name == "Ehad" )
    {
        if( useTrue )      return ( GetEhadTrue() ) * mev_to_gev;
        else               return GetRecoilEnergy()* mev_to_gev;
    }
    if ( name == "CCQE-Recoil")
    {
        //not sure what true would mean for this
        if( useTrue)   return ( GetEhadTrue() ) * mev_to_gev;
        else           return GetCCQERecoil() * mev_to_gev;
    }
    if( name == "VtxEnergy")
    {
        //true means nothing here
        if( useTrue ) return 0.;
        else          return GetVtxEnergy() * mev_to_gev;
    }
    if( name == "ETheta" )
    {
     if( useTrue ) return ( (GetElepTrue() * mev_to_gev ) * ( 1. - cos(GetThetamuTrue()) ));
        else  return ((GetEmu()  * mev_to_gev ) * ( 1. - cos( GetThetamu())) );
    }
    
    if( name == "x" ){
       
            if( useTrue ) return GetxTrue();
            else return GetxReco();
    }
    
    if( name == "y" )
    {
        if( useTrue )    return GetyTrue(); 
        else             return GetyReco();
    }
    if( name == "Q2" || name == "q2" )
    {
        if( useTrue ){
     
            return (GetQ2IncTrue());
        }
        
        else             return ( GetQ2Reco()  );
    }
    
    
    if( name == "W" )
    {
        if( useTrue ) return GetWTrue();
        else return GetWReco() ;
    }
    
    
    if( name == "PhiMu" )
    {
        if( useTrue )     return GetDouble("truth_muon_phi")*rad_to_deg;
        else              return       GetDouble("muon_phi")*rad_to_deg;
    }
    
    if( name == "ThetaMu" )
    {
        if( useTrue )      return GetThetamuTrue()*rad_to_deg;
        else               return     GetThetamu()*rad_to_deg;
    }
    
    if( name == "muonPt" )
    {
        if( useTrue )      return GetlepPtTrue()*rad_to_deg;
                       return     GetMuonPt()*rad_to_deg;
    }
    
    if( name == "CosThetaMu" )
    {
        if( useTrue )      return cos(GetThetamuTrue());
        else               return cos(GetThetamu());
    }
    
    
    if( name == "moduleDNN" )
    {
        if( useTrue ) return GetInt("truth_vtx_module");
        else          return GetInt("ANN_vtx_module");
    }
    if( name == "planeDNN" )
    {
        if( useTrue) return GetInt("truth_vtx_module")*2+GetInt("truth_vtx_plane")+10;
        else{
           //return ANN_vtx_modules[0]*2+ANN_vtx_planes[0]+10;
           int ANN_vtx_module, ANN_vtx_plane;
           int targetID = GetTargetFromSegment( GetVecElem("ANN_segments",0), ANN_vtx_module, ANN_vtx_plane );
    return ANN_vtx_module*2+ANN_vtx_plane+10;
        }
    }
    
    
    cout << "Error [CVUniverse::Var] : I do not know how to interpret that variable \"" << name << "\", so I am throwing an exception." << endl;
    throw 1;
    
}


int CVUniverse::GetAtomicNumber(){
   return GetInt("mc_targetA");
}

//========================================================================================================
//    Adding Bunch of Mapping Function To Split The Nuclear Target Region Based on The Target Number
//========================================================================================================
int CVUniverse::AtomicNumberToMass(int targetZ ) const{
    if( targetZ == 6 ) return 12;
    if( targetZ == 26 ) return 56;
    if( targetZ == 207 ) return 207;
    return -999;
}

int CVUniverse::GetTargetMinBin( int targetID ) const{ 
    if( targetID == 1 ) return 1;
    else if( targetID == 2 ) return 11;
    else if( targetID == 3 ) return 21;
    else if( targetID == 4 ) return 41;
    else if( targetID == 5 ) return 51;
    //else if( targetID == 6 ) return 57;  //faiza
    else{
        std::cerr << "Bad targetID choice in CVUniverse::GetTargetMinBin! Try again. Hint: we only have 5 passive target" << std::endl;
    //}
        return -999;
   }
}

int CVUniverse::GetTargetMaxBin( int targetID ) const{ 
    if( targetID == 1 ) return 19;
    else if( targetID == 2 ) return 29;
    else if( targetID == 3 ) return 41;
    else if( targetID == 4 ) return 55;
    else if( targetID == 5 ) return 65;
    //else if( targetID == 0 ) return 9;   //faiza
    else{
        std::cerr << "Bad targetID choice in CVUniverse::GetTargetMaxBin! Try again. Hint: we only have 5 passive target" << std::endl;
   // }
        return -999;
    }
}

int CVUniverse::GetTargetPlane( int targetID ) const{
    if( targetID == 1 ) return 9;
    else if( targetID == 2 ) return 19;
    else if( targetID == 3 ) return 30;
    else if( targetID == 4 ) return 49;
    else if( targetID == 5 ) return 55;
    else{
        std::cerr << "Bad targetID choice in CVUniverse::GetTargetMaxBin! Try again. Hint: we only have 5 passive target" << std::endl;
     }
    return -999;
    //}
}

int CVUniverse::GetPlaneTargetID( int plane ) const{
    if( plane == 9 ) return 1;
    else if( plane == 19 ) return 2;
    else if( plane == 30 ) return 3;
    else if( plane == 49 ) return 4;
    else if( plane == 55 ) return 5;
    //else{
    //std::cerr << "Bad targetID choice in CVUniverse::GetPlaneTargetID! Try again." << std::endl;
    //}
    return -999;
}

int CVUniverse::GetTargetUSPlane( int targetID ) const{
    if( targetID == 1 ) return 8;
    else if( targetID == 2 ) return 18;
    else if( targetID == 3 ) return 28;
    else if( targetID == 4 ) return 48;
    else if( targetID == 5 ) return 54;
    else if( targetID == 6 ) return 64;
    //  else{
    //    std::cerr << "Bad targetID choice in CVUniverse::GetTargetUSPlane! Try again. Hint: we only have 5 passive target" << std::endl;
    //  }
    return -999;
}

int CVUniverse::GetTargetDSPlane( int targetID ) const{
    if( targetID == 1 ) return 11;
    else if( targetID == 2 ) return 21;
    else if( targetID == 3 ) return 33;
    else if( targetID == 4 ) return 51;
    else if( targetID == 5 ) return 57;
    else if( targetID == 0 ) return 1;
    //  else{
    //    std::cerr << "Bad targetID choice in CVUniverse::GetTargetUSPlane! Try again. Hint: we only have 5 passive target" << std::endl;
    //  }
    return -999;
}

double CVUniverse::GetTargetZStart( int targetID ) const
{
    if( targetID == 1 ) return 4469.46;
    else if( targetID == 2 ) return 4689.73;
    else if( targetID == 3 ) return 4907.82;
    else if( targetID == 4 ) return 5637.73;
    else if( targetID == 5 ) return 5769.88;
    //  else{
    //    std::cerr << "Bad targetID choice in CVUniverse::GetTargetUSPlane! Try again. Hint: we only have 5 passive target" << std::endl;
    //  }
    return -999.9;
}

double CVUniverse::GetTargetZEnd( int targetID ) const
{
    if( targetID == 1 ) return 4495.2;
    else if( targetID == 2 ) return 4715.47;
    else if( targetID == 3 ) return 4984.96;
    else if( targetID == 4 ) return 5645.73;
    else if( targetID == 5 ) return 5782.87;
    //  else{
    //    std::cerr << "Bad targetID choice in CVUniverse::GetTargetUSPlane! Try again. Hint: we only have 5 passive target" << std::endl;
    //  }
    return -999.9;
}

int CVUniverse::GetTargetFromSegment( int segment, int& vtx_module, int& vtx_plane ) const
{
    //int targetID = -999;
    int targetID = -1;
    vtx_module = -999;
    vtx_plane = -999;
    
    if( segment == 0 ){
        vtx_module = -5;
        vtx_plane = 0;
        targetID = -1;
    }
    else if( segment == 1 ){
        vtx_module = -5;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 2 ){
        vtx_module = -5;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 3 ){
        vtx_module = -4;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 4 ){
        vtx_module = -4;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 5 ){
        vtx_module = -3;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 6 ){
        vtx_module = -3;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 7 ){
        vtx_module = -2;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 8 ){
        vtx_module = -2;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 1;
    }
    else if( segment == 9 ){
        vtx_module = -1;
        vtx_plane = 1;
        //targetID = -1;
        targetID = 1;
    }
    else if( segment == 10 ){
        vtx_module = 0;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 1;
    }
    else if( segment == 11 ){
        vtx_module = 0;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 1;
    }
    else if( segment == 12 ){
        vtx_module = 1;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 13 ){
        vtx_module = 1;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 14 ){
        vtx_module = 2;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 15 ){
        vtx_module = 2;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 16 ){
        vtx_module = 3;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 17 ){
        vtx_module = 3;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 2;
    }
    else if( segment == 18 ){
        vtx_module = 4;
        vtx_plane = 1;
        //targetID = -1;
        targetID = 2;
    }
    else if( segment == 19 ){
        vtx_module = 5;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 2;
    }
    else if( segment == 20 ){
        vtx_module = 5;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 2;
    }
    else if( segment == 21 ){
        vtx_module = 6;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 22 ){
        vtx_module = 6;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 23 ){
        vtx_module = 7;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 24 ){
        vtx_module = 7;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 25 ){
        vtx_module = 8;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 26 ){
        vtx_module = 8;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 3;
    }
    else if( segment == 27 ){
        vtx_module = 9;
        vtx_plane = 2;
        //targetID = -1;
        targetID = 3;
    }
    else if( segment == 28 ){
        vtx_module = 11;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 3;
    }
    else if( segment == 29 ){
        vtx_module = 11;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 3;
    }
    else if( segment == 30 ){
        vtx_module = 12;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 31 ){
        vtx_module = 12;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 32 ){
        vtx_module = 13;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 33 ){
        vtx_module = 13;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 34 ){
        vtx_module = 14;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 35 ){
        vtx_module = 14;
        vtx_plane = 2;
        targetID = -1;
    }
   //segment 36 is water!!
    else if( segment == 36 ){
        vtx_module = -999;
        vtx_plane = -999;
        targetID = 6;
    }
    else if( segment == 37 ){
        vtx_module = 15;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 38 ){
        vtx_module = 15;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 39 ){
        vtx_module = 16;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 40 ){
        vtx_module = 16;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 41 ){
        vtx_module = 17;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 42 ){
        vtx_module = 17;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 43 ){
        vtx_module = 18;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 4;
    }
    else if( segment == 44 ){
        vtx_module = 18;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 4;
    }
    else if( segment == 45 ){
        vtx_module = 19;
        vtx_plane = 1;
        //targetID = -1;
        targetID = 4;
    }
    else if( segment == 46 ){
        vtx_module = 20;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 4;
    }
    else if( segment == 47 ){
        vtx_module = 20;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 48 ){
        vtx_module = 21;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 5;
    }
    else if( segment == 49 ){
        vtx_module = 21;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 5;
    }
    else if( segment == 50 ){
        vtx_module = 22;
        vtx_plane = 1;
        //targetID = -1;
        targetID = 5;
    }
    else if( segment == 51 ){
        vtx_module = 23;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 5;
    }
    else if( segment == 52 ){
        vtx_module = 23;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 53 ){
        vtx_module = 24;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 54 ){
        vtx_module = 24;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 55 ){
        vtx_module = 25;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 56 ){
        vtx_module = 25;
        vtx_plane = 2; 
        targetID = -1;
    }
    else if( segment == 57 ){
        vtx_module = 26;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 58 ){
        vtx_module = 26;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 59 ){
        vtx_module = 27;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 60 ){
        vtx_module = 27;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 61 ){
        vtx_module = 28;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 62 ){
        vtx_module = 28;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 63 ){
        vtx_module = 29;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 64 ){
        vtx_module = 29;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 65 ){
        vtx_module = 30;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 66 ){
        vtx_module = 30;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 67 ){
        vtx_module = 31;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 68 ){
        vtx_module = 31;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 69 ){
        vtx_module = 32;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 70 ){
        vtx_module = 32;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 71 ){
        vtx_module = 33;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 72 ){
        vtx_module = 33;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 73 ){
        vtx_module = 34;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 74 ){
        vtx_module = 34;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 75 ){
        vtx_module = 35;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 76 ){
        vtx_module = 35;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 77 ){ 
        vtx_module = 36;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 78 ){
        vtx_module = 36;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 79 ){
        vtx_module = 37;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 80 ){
        vtx_module = 37;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 81 ){
        vtx_module = 38;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 82 ){
        vtx_module = 38;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 83 ){
        vtx_module = 39;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 84 ){
        vtx_module = 39;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 85 ){
        vtx_module = 40;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 86 ){
        vtx_module = 40;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 87 ){
        vtx_module = 41;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 88 ){
        vtx_module = 41;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 89 ){
        vtx_module = 42;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 90 ){
        vtx_module = 42;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 91 ){
        vtx_module = 43;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 92 ){
        vtx_module = 43;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 93 ){
        vtx_module = 44;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 94 ){
        vtx_module = 44;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 95 ){
        vtx_module = 45;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 96 ){
        vtx_module = 45;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 97 ){
        vtx_module = 46;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 98 ){
        vtx_module = 46;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 99 ){
        vtx_module = 47;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 100 ){
        vtx_module = 47;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 101 ){
        vtx_module = 48;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 102 ){
        vtx_module = 48;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 103 ){
        vtx_module = 49;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 104 ){
        vtx_module = 49;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 105 ){
        vtx_module = 50;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 106 ){
        vtx_module = 50;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 107 ){
        vtx_module = 51;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 108 ){
        vtx_module = 51;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 109 ){
        vtx_module = 52;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 110 ){
        vtx_module = 52;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 111 ){
        vtx_module = 53;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 112 ){
        vtx_module = 53;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 113 ){
        vtx_module = 54;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 114 ){
        vtx_module = 54;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 115 ){
        vtx_module = 55;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 116 ){
        vtx_module = 55;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 117 ){
        vtx_module = 56;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 118 ){
        vtx_module = 56;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 119 ){
        vtx_module = 57;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 120 ){
        vtx_module = 57;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 121 ){
        vtx_module = 58;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 122 ){
        vtx_module = 58;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 123 ){
        vtx_module = 59;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 124 ){
        vtx_module = 59;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 125 ){
        vtx_module = 60;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 126 ){
        vtx_module = 60;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 127 ){
        vtx_module = 61;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 128 ){
        vtx_module = 61;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 129 ){
        vtx_module = 62;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 130 ){
        vtx_module = 62;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 131 ){
        vtx_module = 63;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 132 ){
        vtx_module = 63;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 133 ){
        vtx_module = 64;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 134 ){
        vtx_module = 64;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 135 ){
        vtx_module = 65;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 136 ){
        vtx_module = 65;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 137 ){
        vtx_module = 66;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 138 ){
        vtx_module = 66;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 139 ){
        vtx_module = 67;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 140 ){
        vtx_module = 67;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 141 ){
        vtx_module = 68;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 142 ){
        vtx_module = 68;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 143 ){
        vtx_module = 69;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 144 ){
        vtx_module = 69;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 145 ){
        vtx_module = 70;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 146 ){
        vtx_module = 70;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 147 ){
        vtx_module = 71;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 148 ){
        vtx_module = 71;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 149 ){
        vtx_module = 72;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 150 ){
        vtx_module = 72;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 151 ){
        vtx_module = 73;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 152 ){
        vtx_module = 73;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 153 ){
        vtx_module = 74;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 154 ){
        vtx_module = 74;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 155 ){
        vtx_module = 75;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 156 ){
        vtx_module = 75;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 157 ){
        vtx_module = 76;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 158 ){
        vtx_module = 76;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 159 ){
        vtx_module = 77;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 160 ){
        vtx_module = 77;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 161 ){
        vtx_module = 78;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 162 ){
        vtx_module = 78;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 163 ){
        vtx_module = 79;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 164 ){
        vtx_module = 79;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 165 ){
        vtx_module = 80;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 166 ){
        vtx_module = 80;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 167 ){
        vtx_module = 81;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 168 ){
        vtx_module = 81;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 169 ){
        vtx_module = 82;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 170 ){
        vtx_module = 82;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 171 ){
        vtx_module = 83;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 172 ){
        vtx_module = 83;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 173 ){
        vtx_module = 84;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 174 ){
        vtx_module = 84;
        vtx_plane = 2;
        targetID = -1;
    }
  return targetID;
 
}

//=====================
// Systematic shift gets
//=====================
//double CVUniverse::GetCCQERecoil( unsigned int shift )
double CVUniverse::GetCCQERecoil( unsigned int shift )
{
    // if( 0 == shift ) return blob_iso_E_nucl + blob_iso_E_tracker + blob_iso_E_ecal + blob_disp_E_nucl + blob_disp_E_tracker + blob_disp_E_ecal;
    // if( 1 == shift ) return blob_iso_alt01_E_nucl + blob_iso_alt01_E_tracker + blob_iso_alt01_E_ecal + blob_disp_alt01_E_nucl + blob_disp_alt01_E_tracker + blob_disp_alt01_E_ecal;
    // if( 2 == shift ) return blob_iso_alt02_E_nucl + blob_iso_alt02_E_tracker + blob_iso_alt02_E_ecal + blob_disp_alt02_E_nucl + blob_disp_alt02_E_tracker + blob_disp_alt02_E_ecal;
    // if( 3 == shift ) return blob_iso_alt03_E_nucl + blob_iso_alt03_E_tracker + blob_iso_alt03_E_ecal + blob_disp_alt03_E_nucl + blob_disp_alt03_E_tracker + blob_disp_alt03_E_ecal;
    // if( 4 == shift ) return blob_iso_alt04_E_nucl + blob_iso_alt04_E_tracker + blob_iso_alt04_E_ecal + blob_disp_alt04_E_nucl + blob_disp_alt04_E_tracker + blob_disp_alt04_E_ecal;
    
    Error("CVUniverse::GetCCQERecoil", "Valid values for the shift are 0-4");
    throw 1;
}

double CVUniverse::GetVtxEnergy( unsigned int shift )
{
    if( 0 == shift ) return GetDouble("blob_vtx_E");
    // if( 1 == shift ) return blob_vtx_alt01_E;
    // if( 2 == shift ) return blob_vtx_alt02_E;
    // if( 3 == shift ) return blob_vtx_alt03_E;
    // if( 4 == shift ) return blob_vtx_alt04_E;
    
    Error("CVUniverse::GetVtxEnergy", "Valid values for the shift are 0-4");
    throw 1;
}


std::vector<std::string> CVUniverse::GetStdEnergyVars( ) const
{
    std::vector<std::string> vars;
    vars.push_back( "Enu" );
    vars.push_back( "planeDNN" );
    vars.push_back( "Ehad" );
    vars.push_back( "Emu" );
   // vars.push_back( "Emum" );
   //   vars.push_back( "W" );
   //  vars.push_back( "Q2" );
   //  vars.push_back( "x" );
  //  vars.push_back( "y" );
  //  vars.push_back( "ThetaMu" );
    
    return vars;
}

std::vector< std::pair<std::string,std::string> > CVUniverse::GetStd2DVars( ) const
{
    std::vector< std::pair<std::string,std::string> > vars;
//    vars.push_back( std::pair<std::string,std::string>( "Emu", "ThetaMu" ) );
 //   vars.push_back( std::pair<std::string,std::string>( "ThetaMu", "PhiMu" ) );
    //vars.push_back( std::pair<std::string,std::string>( "Enu", "ThetaMu" ) );
    //vars.push_back( std::pair<std::string,std::string>( "Enu", "Ehad" ) );
    vars.push_back( std::pair<std::string,std::string>( "planeDNN", "Ehad" ) );
    //vars.push_back( std::pair<std::string,std::string>( "Ehad", "ThetaMu" ) );
    //vars.push_back( std::pair<std::string,std::string>( "muonPt", "ThetaMu" ) );
    //vars.push_back( std::pair<std::string,std::string>( "W", "Q2" ) );
    //vars.push_back( std::pair<std::string,std::string>( "x", "y" ) );
    return vars;
}

// ARACHNE EVENT DISPLAY
// data
int CVUniverse::GetRunN() const{return GetInt("ev_run");}
int CVUniverse::GetSubRunN() const{return GetInt("ev_subrun");}
int CVUniverse::GetGateN() const{return GetInt("ev_gate");}
int CVUniverse::GetSliceN() const{return GetVecElem("slice_numbers",0);}

// MC
int CVUniverse::GetMCRunN() const{return GetInt("mc_run");}
int CVUniverse::GetMCSubRunN() const{return GetInt("mc_subrun");}
int CVUniverse::GetMCGateN() const{return GetInt("mc_nthEvtInFile");}
int CVUniverse::GetMCSliceN() const{return GetVecElem("slice_numbers",0);}


#endif
