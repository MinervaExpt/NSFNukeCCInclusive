#ifndef MNV_COMPUTEUTILS_cxx
#define MNV_COMPUTEUTILS_cxx 1

#include "include/ComputeUtils.h"

#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"

using namespace Neutron_ANA;

//=============================================================================
// Cone
//=============================================================================

Cone::Cone( XYZVector &vtx, XYZVector &axis, double angle )
{
  SetVtx( vtx );
  SetAxis( axis );
  SetAngle( angle );
  m_xz.SetXYZ( 0, 1 , 0);  
  m_uz.SetXYZ( Sqrt(3)/2, .5 , 0); 
  m_xz.SetXYZ( -Sqrt(3)/2, .5 , 0);
}

XYZVector Cone::GetNorm( int view )
{
  try{
    if ( view > 3 || view < 1 ) throw "View not defined!";
  } catch ( const char* msg )
  {
    cout<<msg<<endl;
  }
  XYZVector norm;
  switch(view)
  {
    case 1: norm = m_xz;
    case 2: norm = m_uz;
    case 3: norm = m_vz;
  }
  return norm;
}

double Cone::GetTPosAngle( int view, double tpos, double zpos )
{
  XYZVector norm = this->GetNorm( view );
  return GetTPosAngle( norm, tpos, zpos );
}


double Cone::GetTPosAngle( XYZVector norm, double tpos, double zpos )
{
  XYZVector tpos_dir = norm.Cross( XYZVector(0,0,1) ).Unit();

  XYZVector vtx_2DT( m_vtx.Dot(tpos_dir),0,m_vtx.Z() );
  XYZVector projectedAxis_2DT = XYZVector( m_axis.Dot( tpos_dir ), 0, m_axis.Z() ).Unit();

  XYZVector blobPos_2DT( tpos, 0,  zpos );
  XYZVector dirWRTvtx = (blobPos_2DT-vtx_2DT).Unit();

  double cosTheta = dirWRTvtx.Dot( projectedAxis_2DT );
  return TMath::ACos( cosTheta );
}


bool Cone::InsideConeAbsPos( int view, double tpos, double zpos ) // view = 1,2,3 for x,u,v; vec = absolute position
{
  XYZVector norm = this->GetNorm( view );
  XYZVector projectedAxis = (this->m_axis - this->m_axis.Dot( norm )*norm).Unit();
  //// see Jaewon's nice figure in Doc 7757 (p. 22) for an explanation of this angle projection
  double projectedAngle = atan2( tan(m_angle), projectedAxis.Dot( m_axis ) );

  return this->GetTPosAngle( norm, tpos, zpos ) <  projectedAngle;
}

bool Cone::InsideCone( XYZVector & posWRTvtx )
{
  double cosTheta = posWRTvtx.Unit().Dot( this->m_axis.Unit() );
  double theta = TMath::ACos( cosTheta );
  return abs(theta) < this->m_angle;
}

bool Cone::InsideConeAbsPos( XYZVector & pos )
{
  XYZVector posWRTvtx = pos - this->m_vtx;
  return InsideCone( posWRTvtx );
}

//=============================================================================
// Detector utils
//=============================================================================

DetectorUtils::DetectorUtils()
:m_beamAngle(-0.05887)
{
    m_beamAxis.SetXYZ( 0, Sin( m_beamAngle ), Cos( m_beamAngle ) );
    m_beamX.SetXYZ(1,0,0);
    m_beamZ = m_beamAxis;
    m_beamY = m_beamZ.Cross( m_beamX );
    m_xz.SetXYZ( 0, 1 , 0);  m_uz.SetXYZ( Sqrt(3)/2, .5 , 0); m_xz.SetXYZ( -Sqrt(3)/2, .5 , 0);
}

XYZVector DetectorUtils::GetBeamCoord( XYZVector &pos )
{
  return XYZVector( pos.Dot(this->m_beamX), pos.Dot(this->m_beamY), pos.Dot(this->m_beamZ) );
}

bool DetectorUtils::IsInHexagon( double x, double y, double apothem, double rotation )
{
  if( x*x + y*y < apothem *apothem ) return true;
  
  double xp = Abs( x*Cos( -rotation ) + y * Sin( -rotation ) );
  double yp = Abs( -x*Sin(-rotation) +y*Cos(-rotation ) );
  double lenOfSide = apothem * 2 / Sqrt(3);
  if ( xp > apothem ) return false;
  if ( yp < lenOfSide / 2 ) return true;
  double slope = lenOfSide/2/apothem;
  if (yp<lenOfSide - xp*slope ) return true;
  return false;
}

XYZVector DetectorUtils::ViewProjection( XYZVector pos, int view )
{
  try{
    if ( view > 3 || view < 1 ) throw "View not defined!";
  } catch ( const char* msg )
  {
    cout<<msg<<endl;
  }
  XYZVector norm;
  switch(view)
  {
    case 1: norm = m_xz;
    case 2: norm = m_uz;
    case 3: norm = m_vz;
  }
  XYZVector projectedAxis = pos - norm * norm.Dot(pos);
  return projectedAxis;
}

//=============================================================================
// GeoUtils
//=============================================================================







std::vector<double> GeoUtils::ComputeNeutronAngularVars( XYZVector &nu, XYZVector &expVec, XYZVector &targetVec )
{	
	XYZVector tgt = targetVec;
  XYZVector tgt_expected = tgt.Unit()*expVec.R(); // this is the inferred neutron momentum, it has same magnitude as the expected momentum, but with the measured direction as the targetVec
	XYZVector coordz = expVec.Unit();
	//XYZVector coordy = coordz.Cross( nu ).Unit();
	//XYZVector coordx = coordy.Cross( coordz );
	XYZVector coordx = coordz.Cross( nu ).Unit();
	XYZVector coordy = coordz.Cross( coordx );
	double dx = coordx.Dot( tgt ), dy = coordy.Dot(tgt), dz = coordz.Dot(tgt);
	double dxi = coordx.Dot( tgt_expected ), dyi = coordy.Dot(tgt_expected), dzi = coordz.Dot(tgt_expected);
  double dPPerp = dx, dPReact = dy;
  double dPPerpI = dxi, dPReactI = dyi;
	
	double ReactPlaneAngle = TMath::ATan2( dy, dz );
	//double PerpPlaneAngle = TMath::ATan2( dy, TMath::Power( (dz*dz+dx*dx), 0.5 ) );
	double PerpPlaneAngle = TMath::ATan2( dx, dz );
  double Theta = ACos( expVec.Unit().Dot( tgt ) );
	vector<double> ret({PerpPlaneAngle, ReactPlaneAngle, Theta, dPPerp, dPReact, dPPerpI, dPReactI});
	//ret.push_back(PerpPlaneAngle);
	//ret.push_back(ReactPlaneAngle);
	return ret;
}

XYZTVector GeoUtils::ComputeExpectedNucleon( XYZVector &nu, XYZTVector &muon, double ISMass,double FSMass, double BindingE ) // vbar p --> mu n, Initial Mass, Final state mass
{
	//rotate to nu = (0,0,1)
	XYZVector nuAxis(0,0,1);
	AxisAngle R1( nu.Cross( nuAxis) , ACos( nuAxis.Dot( nu.Unit() )) );
	
	XYZTVector muonR1 = (R1*muon);
	
	//(v00v) + (M000) = (Em pxpy pz ) + (EM -px-py pmz)
	// v-EM = Em-M, v- pmz = pz
	//EM - pmz = pz - Em + M0
	// A = pz - Em + M0
	// M1^2 + px^2 + py^2 + pmz^2 = A^2 + pmz^2 + 2A*pmz, 
	// A = pz - Em + Mi, Mi = M0 - Eb
	// B^2 = M1^2 + px^2 + py^2
	// B^2 = A^2 + 2A*pmz ===> pmz = (B^2 - A^2 )/(2A)
	double B2 = FSMass*FSMass + muonR1.X()*muonR1.X() + muonR1.Y()*muonR1.Y();
	double A = muonR1.Z() - muonR1.E() + (ISMass-BindingE);
	double pMz = (B2/A-A)/2.;
	double pMx = -muonR1.X(), pMy = -muonR1.Y();
	XYZTVector nucleonR1( pMx, pMy, pMz, Sqrt( pMx*pMx + pMy*pMy + pMz*pMz + FSMass*FSMass ) );
	
	//rotate back
	return ( R1.Inverse()*nucleonR1 );	
}

XYZVector GeoUtils::BeamAxis( double bias )
{
  double angle = this->beam_angle() + bias;
  return XYZVector( 0, Sin( angle ), Cos( angle ) );
}

double GeoUtils::DistanceFromPointToLine(XYZVector &pt, XYZVector &line )
{
	return pt.Cross( line.Unit() ).R();
}

vector<double> GeoUtils::GetTransverseVariables( XYZVector &nu3P, XYZTVector &muon, XYZTVector &nucleon , double NucleusMass, double bindingE)
{
  double RemnantNucleusMass = NucleusMass - this->MnGeV() + bindingE;

  XYZVector  muon3P = muon.Vect(), nucleon3P = nucleon.Vect();
  double Enu = -999., Emuon = muon.E(), Enucleon = nucleon.E();

  XYZVector coordz = nu3P.Unit();
  XYZVector coordy = -(muon3P-muon3P.Dot(coordz)*coordz ).Unit();
  XYZVector coordx = coordz.Cross(coordy);


  //prime coordinate of coord(x,y,z)
  XYZVector pNucleon3P( nucleon3P.Dot(coordx),nucleon3P.Dot(coordy),nucleon3P.Dot(coordz));
  XYZVector pMuon3P( muon3P.Dot(coordx),muon3P.Dot(coordy),muon3P.Dot(coordz));

  XYZVector pNucleonPT( pNucleon3P.X(), pNucleon3P.Y(), 0 );
  XYZVector pMuonPT( pMuon3P.X(), pMuon3P.Y(), 0 );

  XYZVector dpTvec = pNucleonPT+pMuonPT;

  double R = NucleusMass + pMuon3P.Z() + pNucleon3P.Z() - Emuon - Enucleon;
  double dPL = 0.5*( R- ( NucleusMass*NucleusMass + dpTvec.Perp2() )/R  );
  double Pn = TMath::Power(dpTvec.Perp2() + dPL*dPL, 0.5 );
  Enu = pMuon3P.Z() + pNucleon3P.Z()-dPL;

  double dalphaT = TMath::ACos( -pMuonPT.Unit().Dot( dpTvec.Unit() )) * 180/TMath::Pi();
  double dpT = dpTvec.R();
  double dpTx = dpTvec.X();
  double dpTy = dpTvec.Y();
  double dphiT = TMath::ACos( -pMuonPT.Unit().Dot( pNucleonPT.Unit() ) ) *180/TMath::Pi();
  double sign = (dpTx>0)? +1 : -1;

  return vector<double>({ Enu, Pn, dpT, dpTx, dpTy, dalphaT, dphiT, sign });
}



PhysicsUtils::PhysicsUtils():
  M_p(938.272013/1000),
  M_n(939.56536/1000),
  M_mu(105.6583/1000)
{

}


double PhysicsUtils::nuEnergyCCQE( XYZVector&nu3P, XYZTVector& lepton, int charge, double bindingE) 
{
  double lep_energy = lepton.E();
  double lep_p = lepton.P();
  double lep_theta = TMath::ACos( nu3P.Unit().Dot( lepton.Vect().Unit() ) );
  double lep_mass = lepton.M();
  double nu_energy = -1000;
  if( charge > 0 ) {
    double nu_energy_num = pow(this->M_n,2) - pow(this->M_p - bindingE,2)
      - pow(lep_mass,2) + 2.0*(this->M_p - bindingE)*lep_energy;
    double nu_energy_den = 2.0*(this->M_p - bindingE - lep_energy + lep_p*cos(lep_theta));
    if( nu_energy_den ) nu_energy = nu_energy_num / nu_energy_den;
  } else if ( charge < 0 ) {
    double nu_energy_num = pow(this->M_p,2) - pow(this->M_n - bindingE,2)
      - pow(lep_mass,2) + 2.0*(this->M_n - bindingE)*lep_energy;
    double nu_energy_den = 2.0*(this->M_n - bindingE - lep_energy + lep_p*cos(lep_theta));
    if( nu_energy_den ) nu_energy = nu_energy_num / nu_energy_den;
  }
  return nu_energy;
}



double PhysicsUtils::qsqCCQE( XYZVector&nu3P, XYZTVector& lepton, int charge, double bindingE) 
{
  double lep_energy = lepton.E();
  double lep_p = lepton.P();
  double lep_theta = TMath::ACos( nu3P.Unit().Dot( lepton.Vect().Unit() ) );
  double lep_mass = lepton.M();

  double nu_energy = nuEnergyCCQE( nu3P, lepton, charge, bindingE );

  if( nu_energy < 0 ) return -1000000;
  double Qsquared  = 2.0 * nu_energy * (lep_energy - lep_p * TMath::Cos(lep_theta)) - TMath::Power(lep_mass,2);
  return Qsquared;
}











#endif
