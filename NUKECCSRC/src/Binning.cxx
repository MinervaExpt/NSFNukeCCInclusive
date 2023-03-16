//#ifndef MNV_NUKECC_cxx
//#define MNV_NUKECC_cxx 1

#ifndef MNV_NUKECC_BINNING_CXX
#define MNV_NUKECC_BINNING_CXX 1

#include "../include/Binning.h"
//#include "../include/MasterAnaDevvars.h"
#include "PlotUtils/TargetUtils.h"
#include "../include/CVUniverse.h"
#include "../include/GlobalIncludes.h" 
#include "../include/CommonIncludes.h"
#include <PlotUtils/MnvNormalization.h>
#include <PlotUtils/NuclModUtils.h>
#include <PlotUtils/FluxReweighter.h>
//#include <PlotUtils/FluxReweighterWithWiggleFit.h>
#include <PlotUtils/MnvNuclearModelWeight.h>
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"

#include "TFileCollection.h"


#include "TVector3.h"

using namespace std;

using namespace NUKECC_ANA;

NukeCC_Binning::NukeCC_Binning()
{
}

NukeCC_Binning::~NukeCC_Binning()
{
}


//======================
// get singleton
//======================
NukeCC_Binning& NukeCC_Binning::Get()
{
    static NukeCC_Binning singleton;
    return singleton;
}

/*
axis_binning NukeCC_Binning::GetQ2BinsGeV() const{
  axis_binning tmp;

  std::vector<double> bins;
  tmp.default_width = 0.05;
  tmp.uniform = false;
  
  bins.push_back(0.0);
  bins.push_back(0.05);
  bins.push_back(0.1);
  bins.push_back(0.2);
  bins.push_back(0.3);
  bins.push_back(0.5);
  bins.push_back(0.75);
  bins.push_back(1.0);
  bins.push_back(1.5);
  bins.push_back(3.0);
  bins.push_back(5.0);
  bins.push_back(8.0);
  bins.push_back(50.0);

  tmp.bin_edges = bins;
  tmp.nbins = bins.size()-1;
  tmp.min = bins.front();
  tmp.max = bins.back();

  return tmp;

}

axis_binning NukeCC_Binning::GetWBinsGeV() const{
  axis_binning tmp;

  std::vector<double> bins;
  tmp.default_width = 0.05;
  tmp.uniform = false;

  bins.push_back(0.00);
  bins.push_back(1.1);
  bins.push_back(1.3);
  bins.push_back(1.7);
  bins.push_back(2.0);
  bins.push_back(2.5);
  bins.push_back(3.5);
  bins.push_back(5.0);
  bins.push_back(10.0);
  tmp.bin_edges = bins;
  tmp.nbins = bins.size()-1;
  tmp.min = bins.front();
  tmp.max = bins.back();

  return tmp;

}


axis_binning NukeCC_Binning::GetEnuBinsGeV() const{
  axis_binning tmp;

  std::vector<double> bins;
  tmp.default_width = 0.05;
  tmp.uniform = false;

  bins.push_back(2.00);
  bins.push_back(3.0);
  bins.push_back(4.0);
  bins.push_back(5.0);
  bins.push_back(6.25);
  bins.push_back(7.5);
  bins.push_back(10.0);
  bins.push_back(25.0);
  bins.push_back(50.0);

  tmp.bin_edges = bins;
  tmp.nbins = bins.size()-1;
  tmp.min = bins.front();
  tmp.max = bins.back();

  return tmp;

}


axis_binning NukeCC_Binning::GetEhadBinsGeV() const{
  axis_binning tmp;

  std::vector<double> bins;
  tmp.default_width = 0.05;
  tmp.uniform = false;

  bins.push_back(0.00);
  bins.push_back(0.25);
  bins.push_back(0.75);
  bins.push_back(1.50);
  bins.push_back(2.50);
  bins.push_back(3.50);
  bins.push_back(5.00);
  bins.push_back(8.00);
  bins.push_back(18.00);
  bins.push_back(25.0);

  tmp.bin_edges = bins;
  tmp.nbins = bins.size()-1;
  tmp.min = bins.front();
  tmp.max = bins.back();

  return tmp;

}
*/
void NukeCC_Binning::AddBins( std::vector<double>& binsLowEdge, const int nBins, const double binWidth ) const
{
    //! If the bins have no defined low entry, assume it should be 0
    double x = 0.;
    if( binsLowEdge.size() > 0 )
        x = binsLowEdge[ binsLowEdge.size() - 1 ];
    else
        binsLowEdge.push_back( x );
    
    //! add nBins spaced by this width
    for( int iBin = 0; iBin < nBins; ++iBin )
    {
        x += binWidth;
        binsLowEdge.push_back( x );
    }
}

void NukeCC_Binning::AddBins( std::vector<double>& binsLowEdge, const double binWidth, const double lowBin, const double highBin ) const
{
    double x = lowBin;
    if( binsLowEdge.size() && lowBin == binsLowEdge.back() )
        x += binWidth;
    while( x <= highBin )
    {
        binsLowEdge.push_back( x );
        x += binWidth;
    }
}


std::vector<double>NukeCC_Binning::GetBins( const std::string& var, BinType::t_BinType type ) const
{
    vector<double> binsLowEdge;
    if(      type == BinType::kEnergy )
        binsLowEdge = GetEnergyBins( var );
    else if( type == BinType::kResolution )
        binsLowEdge = GetResolutionBins( var );
    else if( type == BinType::kResidual )
        binsLowEdge = GetResidualBins( var );
    else if( type == BinType::kPosition )
        binsLowEdge = GetPosBins( var );
    else if( type == BinType::kSideband )
        binsLowEdge = GetSidebandBins( var );
    else if( type == BinType::kTruth )
        binsLowEdge = GetTruthBins( var );
    else if( type == BinType::kDeltaE )
        binsLowEdge = GetDeltaEBins( var );
    
    if( binsLowEdge.empty() )
        Warning( "NukeCC_Binning::GetBins", Form("No bins found for variable %s, bintype %d", var.c_str(), int(type) ) );
    
    return binsLowEdge;
}

std::vector<double> NukeCC_Binning::GetEnergyBins( const std::string& var, bool DISbins /*=true*/ ) const
{
    vector<double> binsLowEdge;
    //if( DISbins ) return binsLowEdge = GetDISBins( var );

    if( var == "Q2"){
      double tmpBins[] = { 0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75,1.0, 1.5, 3.0, 5.0, 8.0, 50. };
      //double tmpBins[] = { 0, 1000.0 };
      //double tmpBins[] = { 0, 0.5, 1.0, 1.5, 5.0, 8.0};
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "W" ){
      double tmpBins[]  = {0, 1.1, 1.3, 1.7, 2., 2.5, 3.5, 5., 8., 10.};
      //double tmpBins[]  = {0, 2., 5., 8., 10.};
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "Emu" || var == "Emum" ){
      //double tmpBins[]  = { 0. ,1. ,2. ,3. ,4. ,5. ,7. ,9. ,12. ,15. ,18. ,22. ,36. ,50. ,75. ,100. ,120. };
      //double tmpBins[]  = {2.0, 3.5, 4.75, 5.5, 7.5, 10., 13., 16., 20., 25., 35., 50 }; //original
      double tmpBins[]  = {2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0}; // changed 20/07/2022
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "Enu" ){
      double tmpBins[]  = { 2.0, 3.0, 4.0, 5.0, 6.25, 7.5, 8.75, 10., 12.5, 15., 20.}; //, 25.0, 30., 40., 50. }; // changed on 22/02/2023
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "EnuFlux" ){
      // 0.1 GeV bins from 0 to 100 GeV
      AddBins( binsLowEdge, 0.1, 0., 100. );
    }

    else if( var == "vtxz_all" ){
      // 0.1 GeV bins from 0 to 100 GeV
      //AddBins( binsLowEdge, 0.05, 0., 1000. );
      double tmpBins[]  = 
        {4293.04/10,
        4313.68/10,
        4337.25/10,
        4357.9/10,
        4381.47/10,
        4402.11/10,
        4425.68/10,
        4446.33/10,
        4514.11/10,
        4534.76/10,
        4558.33/10,
        4578.97/10,
        4602.54/10,
        4623.19/10,
        4646.76/10,
        4667.4/10,
        4735.19/10,
        4755.83/10,
        4779.4/10,
        4800.05/10,
        4823.62/10,
        4844.26/10,
        4867.83/10,
        4888.48/10,
        5000.48/10,
        5021.12/10,
        5044.69/10,
        5065.34/10,
        5088.91/10,
        5109.55/10,
        5133.12/10,
        5153.77/10,
        5456.74/10,
        5477.38/10,
        5500.95/10,
        5521.6/10,
        5545.17/10,
        5565.81/10,
        5589.38/10,
        5610.02/10,
        5677.81/10,
        5698.45/10,
        5722.03/10,
        5742.67/10,
        5810.45/10,
        5831.1/10,
        5855.68/10,
        5876.33/10,
        5900.91/10,
        };
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }


    else if( var == "Ehad" ){
      double tmpBins[]  = { 0., 0.25, 0.75, 1.5, 2.5, 3.5, 5.0, 8.0, 12., 18., 25. };
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "x" ){

      binsLowEdge.push_back(0.001);
      binsLowEdge.push_back(0.05);
      binsLowEdge.push_back(0.1);
      binsLowEdge.push_back(0.2);
      binsLowEdge.push_back(0.4);
      binsLowEdge.push_back(1); 
      binsLowEdge.push_back(2.2);
    }
    else if( var == "x09" ){

      binsLowEdge.push_back(0.);
      binsLowEdge.push_back(0.05);
      binsLowEdge.push_back(0.1);
      binsLowEdge.push_back(0.2);
      binsLowEdge.push_back(0.4);
      binsLowEdge.push_back(0.9);
      binsLowEdge.push_back(2.2);
    }
    else if( var == "xBrian" ){

      binsLowEdge.push_back(0.001);
      binsLowEdge.push_back(0.1);
      binsLowEdge.push_back(0.3);
      binsLowEdge.push_back(0.7);
      binsLowEdge.push_back(0.9);
      binsLowEdge.push_back(1.1);
      binsLowEdge.push_back(1.5);

    }
     else if( var == "xfine" ){

      binsLowEdge.push_back(0.001);
      binsLowEdge.push_back(0.01);
      binsLowEdge.push_back(0.02);
      binsLowEdge.push_back(0.03);
      binsLowEdge.push_back(0.04);
      binsLowEdge.push_back(0.05);
      binsLowEdge.push_back(0.1);
      binsLowEdge.push_back(0.3);
      binsLowEdge.push_back(0.7);
      binsLowEdge.push_back(0.9);
      binsLowEdge.push_back(1.1);
      binsLowEdge.push_back(1.5);

    }
    else if( var == "y" )
    {
       double tmpBins[] = { 0, 0.03, 0.15, 0.3, 0.45, 0.6, 1.};
       binsLowEdge.assign(tmpBins,tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }


    else if( var == "CCQE-Recoil" || var == "VtxEnergy" )
    {
        //.5 * recoil for now
        binsLowEdge = GetBins( "Ehad" );
        for( vector<double>::iterator i = binsLowEdge.begin(); i != binsLowEdge.end(); ++i )
            *i = .5*(*i);
    }
    else if( var == "ETheta" )
    {
        AddBins( binsLowEdge, 4, .02 );
        AddBins( binsLowEdge, 2, .04 );
        AddBins( binsLowEdge, 2, .08 );
    }
    else if( var == "PhiMu" )
    {
        //4 deg bins from -180 to 180
        AddBins( binsLowEdge, 12., -180, 180 );
    }
    else if( var == "ThetaMu")
    {
     double tmpBins[]  = {0, 1, 2, 3, 5, 7, 9, 11, 14, 17}; // changed 24/07/2022
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "ThetaMuX" || var == "ThetaMuY" )
    {
        //3 deg bins from -29 to -20
        AddBins( binsLowEdge, 3., -29., -20. );
        //---underflow---//
        //2 deg bins from -20 to -10
        AddBins( binsLowEdge, 2., -20., -10. );
        //1 deg bins from -10 to 10
        AddBins( binsLowEdge, 1., -10., 10. );
        //2 deg bins from 10 to 20
        AddBins( binsLowEdge, 2., 10., 20. );
        //--- overflow ---//
        //3 deg bins from 20 to 29
        AddBins( binsLowEdge, 3., 20., 29. );
    }
    else if( var == "CosThetaMu" )
    {
        binsLowEdge = GetBins( "ThetaMu" );
        for( vector<double>::iterator i = binsLowEdge.begin(); i != binsLowEdge.end(); ++i )
            *i = cos( deg_to_rad * (*i) );
        reverse( binsLowEdge.begin(), binsLowEdge.end() );
    }
    
    else if( var == "NTracks" )
    {
        // 1-10
        AddBins( binsLowEdge, 10, 1. );
    }
    
    else if( var == "plane" || var == "planeDNN"  )
    {
        // 0-66
        double tmpBins[] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 11, 12, 13., 14., 15., 16., 17., 18., 19., 21., 22., 23., 24., 25., 26., 27., 28., 29., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 51., 52., 53., 54., 55., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66.,
        67., 68., 69., 70., 71., 72.,73.,74.,75.,76.,77.,78.,79.,80.,81.,82.,83.,84.,85.,86.,87.,88.,89.,90.,91.,92.,93.,94.,95.,96.,97.,98.,99.,100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,110.,111.,112.,113.,114.,115.,116.,117.,118.,119.,120., 121.,122.,123.,124.,125.,126.,127.,128.,129.,130.,131.,132.,133.,134.,
        135.,136.,137.,138.,139.,140.,141.,142.,143.,144.,145.,146.,147.,148.,149.,150.,151.,152.,153.,154.,155.,156.,157.,158.,159.,160.,161.,162.,163.,164.,165.,166.,167.,168.,169.,170.,171.,172.,173. };
        //double tmpBins[] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 11., 12., 16., 17., 18., 19., 21., 22., 28., 29., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 51., 52., 53., 54., 55., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66. };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "vtxz" || var == "vtxy" || var == "vtxx"){
        binsLowEdge = GetPosBins( var ); 
    }

    else if(var == "muonPt"){
       double tmpBins[] = { 0.0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 1.25, 1.5, 2.5, 4.5};
       binsLowEdge.assign(tmpBins,tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    } 

    else if(var == "muonPt1D"){
       double tmpBins[] = { 0.0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 1.25, 1.5, 2.5}; // 16/03/2023 pT1D
       binsLowEdge.assign(tmpBins,tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    } 

    else if(var == "muonPz"){
       double tmpBins[] = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60};
       binsLowEdge.assign(tmpBins, tmpBins+sizeof(tmpBins) / sizeof(tmpBins[0]));
   }

   else if(var == "muonPz1D"){
       double tmpBins[] = {2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20}; //, 40, 60}; // Changed on 22/02/2023
       binsLowEdge.assign(tmpBins, tmpBins+sizeof(tmpBins) / sizeof(tmpBins[0]));
   }

   else if(var == "muonP"){
       double tmpBins[] = {0, 2, 4, 6, 8, 10, 12, 14, 15, 16, 18, 20, 30, 40, 60};
       binsLowEdge.assign(tmpBins, tmpBins+sizeof(tmpBins) / sizeof(tmpBins[0]));
   }

    else if(var=="ANNPlaneProb"){ //finer binning to see the differences
        binsLowEdge.push_back( 0 ); // NEW
        AddBins(binsLowEdge, 20, 0.05); //NEw

    }

    else if(var=="petal"){ //finer binning to see the differences
        double tmpBins[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        binsLowEdge.assign(tmpBins, tmpBins+sizeof(tmpBins) / sizeof(tmpBins[0]));
    }

    else if(var=="planeDiff"){ //finer binning to see the differences
        double tmpBins[] = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5};
        //double tmpBins[] = {-5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
        binsLowEdge.assign(tmpBins, tmpBins+sizeof(tmpBins) / sizeof(tmpBins[0]));
    }

    else if (var == "mc_run"){
    
        AddBins(binsLowEdge, 1., 111325, 111375); 
        // take binning from here, depending on the playlist
        //https://cdcvs.fnal.gov/redmine/projects/minerva-sw/wiki/Inextinguishable_Production 
    }

    return binsLowEdge;
}


std::vector<double> NukeCC_Binning::GetDISBins(const std::string& var ) const
{


    vector<double> binsLowEdge;
    if( var == "Q2"){
        double tmpBins[] = { 1.0, 1.25, 3.0, 5.0, 8.0, 30. };
        //For cross check with Maya
        //double tmpBins[] = { 1.0, 1.25, 3.0, 5.0, 8.0 };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "W" )
    {
        //for DIS: one underflow bin (0 - 2 GeV) same otherwise
        double tmpBins[]  = {2., 2.5, 3.5, 5., 8., 16.};
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "Emu" || var == "Emum" )
    {
        //for DIS: one underflow bin (0 - 2 GeV) same otherwise
        //double tmpBins[]  = { 2., 3., 4., 5., 6., 8., 10., 13., 16., 20., 25., 35., 50., 120. };
        double tmpBins[]  = {2.0, 3.5, 4.75, 5.5, 7.5, 10., 13., 16., 20., 25., 35., 50 };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "Enu" )
    {
        double tmpBins[]  = { 5., 7.5, 10., 12.5, 15., 20., 30., 40., 50. };
       //For Data Quality 
       //double tmpBins[]  = { 0., 1.0, 2.0, 3.0, 4.0,5.0, 6.0, 7.0, 8.0, 9.0, 10., 11.0, 12.0, 13.0, 14.0, 15.0, 17.5, 20.0, 22.5, 25., 27.5, 30. };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "Ehad" )
    {
        //for DIS: one underflow bin (0 - 2 GeV) same otherwise
        double tmpBins[]  = { 1.5, 3., 5., 8., 11., 15., 20., 25. };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "x" )
    {
        
        //from Jorge
        //I guess this should be overflow since we don't have anything generated below 0.005
        //binsLowEdge.push_back(0.);
        
        // 0-.1 (shadowing)
        binsLowEdge.push_back(0.005);
        binsLowEdge.push_back(0.05);
        binsLowEdge.push_back(0.1);
        
        //.1-.3 (antishadowing)
        binsLowEdge.push_back(0.2);
        binsLowEdge.push_back(0.3);
        
        //.3-.8 (EMC)
        binsLowEdge.push_back(0.5);
        binsLowEdge.push_back(0.8);
        //binsLowEdge.push_back(1.5);
    }
    else if( var == "y" )
    {
       double tmpBins[] = { 0, 0.03, 0.15, 0.3, 0.45, 0.6, 1.};
       //For check with Maya studies
       //double tmpBins[] = { 0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
       binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }


    else if( var == "ThetaMu")
    {
        // 1 deg bins up to 8
        //AddBins( binsLowEdge, 10, 1.5 );
        AddBins( binsLowEdge, 8, 1.0);
        // 2 deg bins up to 14,
        //AddBins( binsLowEdge, 3, 2.0 );
        // 3 deg bins up to 26
        AddBins( binsLowEdge, 4, 3. );
        //--- overflow --//
        // 5 deg bins up to 51
        AddBins( binsLowEdge, 5, 5. );
        // 2 20 deg bins up to 91
        AddBins( binsLowEdge, 2, 20. );
    }
    return binsLowEdge;
    
}


std::vector<double> NukeCC_Binning::GetSidebandBins(const std::string& var ) const
{
    vector<double> binsLowEdge;
    
    if( var == "W" )
    {
        AddBins( binsLowEdge, .05, 0., 25. );
        //for DIS: one underflow bin (0 - 2 GeV) same otherwise
        //double tmpBins[]  = {0, 1.3, 1.5, 1.7, 2.0, 2.5, 3., 3.5, 4., 4.5, 5., 6., 10., 25, 100.};
        //binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "Q2" )
    {
        AddBins( binsLowEdge, .05, 0., 25. );
        //for DIS: one underflow bin (0 - 2 GeV) same otherwise
        //double tmpBins[]  = {0, 1.3, 1.5, 1.7, 2.0, 2.5, 3., 3.5, 4., 4.5, 5., 6., 10., 25, 100.};
        //binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else
        binsLowEdge =  GetEnergyBins( var );
    
    return binsLowEdge;
    
}

std::vector<double> NukeCC_Binning::GetResidualBins( const std::string& var ) const
{
    vector<double> binsLowEdge;
    if( var == "Emu" )
    {
        // .2 GeV bins from -2 to 2
        AddBins( binsLowEdge, 0.2, -6.1, 6.1 );
    }
    else if( var == "Enu" )
    {
        // .2 GeV bins from -5 to 5
        AddBins( binsLowEdge, 0.2, -4.1, 4.1 );
    }
    else if( var == "Emum")
    {
        // .2 GeV bins from -2 to 2
        AddBins( binsLowEdge, .2, -2.1, 2.1 );
    }
    else if( var == "Ehad" || var == "CCQE-Recoil" || var == "VtxEnergy" )
    {
        AddBins( binsLowEdge, 0.4, -10, 10 );
    }
    else if( var == "PhiMu" )
    {
        // -10 to 10 deg in bins of .25 deg
        AddBins( binsLowEdge, .5, -15.25, 15.15 );
    }
    else if( var == "ThetaMu" || var == "ThetaMuX" || var == "ThetaMuY" )
    {
        // -5 deg to 5 deg in bins of .2 deg
        AddBins( binsLowEdge, .2, -5.1, 5.1 );
    }
    else if( var == "x" )
    {
        AddBins( binsLowEdge, .05, -1.025, 1.025 );
    }
    else if( var == "y" )
    {
        AddBins( binsLowEdge, 0.02, -.31, .31 );
   }
    else if( var == "Q2" )
    {
        AddBins( binsLowEdge, .1, -1.05, 1.05 );
    }
    else if( var == "W" )
    {
        AddBins( binsLowEdge, .2, -2.1, 2.1 );
    }
    else if( var == "vtxt" )
    {
        // 10ns bins from -200 to 200
        AddBins( binsLowEdge, 10., -205., 205.);
    }
    else if( //all transverse directions get the same bins
            string::npos != var.find("vtxx") ||
            string::npos != var.find("vtxy") ||
            string::npos != var.find("vtxd") ||
            string::npos != var.find("vtxu") ||
            string::npos != var.find("vtxc")
            )
    {
        AddBins( binsLowEdge, .05, -2.525,2.525);
    }
    else if( string::npos != var.find("vtxz" ) )
    {
        //Try to make these ~1 strip wide
        AddBins( binsLowEdge, 17, -150, 420 );
    }
    else if ( string::npos != var.find("vtxr" ) )
    {
        AddBins( binsLowEdge, .05, -.075, 2.425);
    }
    
    
    return binsLowEdge;
}

std::vector<double> NukeCC_Binning::GetResolutionBins( const std::string& var ) const
{
    vector<double> binsLowEdge;
    
    AddBins( binsLowEdge, .05, -1.25, 2.25 );
    
    return binsLowEdge;
}

std::vector<double> NukeCC_Binning::GetDeltaEBins( const std::string& var ) const
{
    vector<double> binsLowEdge;
    
    AddBins( binsLowEdge, .01, -0.35, 1.25 );
    
    return binsLowEdge;
}


//==========================================
// start of Det bins
//==========================================
const unsigned int NukeCC_Binning::fNDetBins;

double NukeCC_Binning::fDetBins[94] = { 
    424.11224, //-5
    428.22447, //-4
    432.64595, //-3
    437.06743, //-2
    441.48891, //tgt1
    445.91039, //0
    450.33188,
    454.75336,
    459.17484, //3
    463.59632, //tgt2
    468.01780, //5
    472.43928,
    476.86077,
    481.28225, //8
    485.70373, //tgt3a
    490.12521, //tgt3b
    494.54669, //11
    498.96817,
    503.38965,
    507.81114,
    512.23262,
    516.44336,
    544.59410,
    549.01558,
    553.43706,
    557.85854,
    562.28002,
    566.70151,
    571.12299,
    575.54447,
    579.91524,
    584.43814,
    588.96104,
    593.48394,
    598.00684,
    602.52974,
    607.05264,
    611.57554,
    616.09844,
    620.62134,
    625.14424,
    629.66714,
    634.19004,
    638.71294,
    643.23584,
    647.75874,
    652.28164,
    656.80454,
    661.32744,
    665.85034,
    670.37324,
    674.89614,
    679.41904,
    683.94194,
    688.46484,
    692.98774,
    697.51064,
    702.03354,
    706.55644,
    711.07934,
    715.60224,
    720.12514,
    724.64804,
    729.17094,
    733.69384,
    738.21674,
    742.73964,
    747.26254,
    751.78544,
    756.30834,
    760.83124,
    765.35414,
    769.87704,
    774.39994,
    778.92284,
    783.44574,
    787.96864,
    792.49154,
    797.01444,
    801.53734,
    806.06024,
    810.58314,
    815.10604,
    819.62894,
    824.15184,
    828.67474,
    833.19764,
    837.72054,
    842.24344,
    846.76634,
    851.28924,
    855.81214,
    860.35635,
    864.83663 //end of mod 85
};
//==========================================
// end of Det bins
//==========================================


//distances here in cm, ns
std::vector<double> NukeCC_Binning::GetPosBins( const std::string& var, const int targetID /* = -1 */ ) const
{
    vector<double> binsLowEdge;
    
    const int xyspace = 2.5; //cm
    const double widthXUDC = 200.;
    const double widthY = 220.;
    if( string::npos != var.find("vtxt") )
    {
        // from 0 to 10000 in 100ns intervals
        static const double deltaT = 100;
        AddBins( binsLowEdge, 10000/deltaT, deltaT );
    }
    if(
       string::npos != var.find("vtxx" ) ||
       string::npos != var.find("vtxu" ) ||
       string::npos != var.find("vtxd" ) ||
       string::npos != var.find("vtxc" )
       )
    {
        binsLowEdge.push_back( -widthXUDC/2. );
        AddBins( binsLowEdge, widthXUDC/xyspace, xyspace );
    }
    else if( string::npos != var.find("vtxy" ) )
    {
        binsLowEdge.push_back( -widthY/2. );
        AddBins( binsLowEdge, widthY/xyspace, xyspace );
    }
    else if (string::npos != var.find("vtxr") )
    {
        binsLowEdge.push_back( 0 );
        AddBins( binsLowEdge, 1.17E2/xyspace, xyspace );
    }
    else if( string::npos != var.find("vtxz") )
    {
        int zspace = 2.;//cm
        if( targetID == -1 )
        {
            binsLowEdge.assign(fDetBins, fDetBins + fNDetBins); //return the detector bins /CHANGE ABBEY
            
            //binsLowEdge.push_back( 0 ); // NEW
            //AddBins(binsLowEdge, (900-0)/5, 5); //NEw
        }
        else if( targetID == 1 )
        {
            // from 440 to 460
            binsLowEdge.push_back( 440. );
            AddBins( binsLowEdge, (460-440)/zspace, zspace );
        }
        else if( targetID == 2 )
        {
            // from 462 to 482
            binsLowEdge.push_back( 462. );
            AddBins( binsLowEdge, (482-462)/zspace, zspace );
        }
        else if( targetID == 3 )
        {   
            // from 485 to 505
            binsLowEdge.push_back( 485 );
            AddBins( binsLowEdge, (505-485)/zspace, zspace );

        }
        else if( targetID == 4 )
        {
            // from 555 to 575
            binsLowEdge.push_back( 555 );
            AddBins( binsLowEdge, (575-555)/zspace, zspace );
        }
        else if( targetID == 5 )
        {
            // from 575 to 585
            binsLowEdge.push_back( 575 );
            AddBins( binsLowEdge, (585-575)/zspace, zspace );
        }
        else
        {
            const int totalLowestMod = 27;
            const double totalLowestZ = 595.;
            const double deltaZPerMod = 4.5; //4.5cm per module
            
            int thisLowestMod = targetID % 1000;
            if( targetID < 100 )
            {
                int fauxSet = ( targetID - targetID%10 ) / 10;
                thisLowestMod = totalLowestMod + 6*(fauxSet-1);
            }
            double lowZ = totalLowestZ + deltaZPerMod * ( thisLowestMod - totalLowestMod );
            double highZ = lowZ + deltaZPerMod * 6;
            
            binsLowEdge.push_back( lowZ );
            AddBins( binsLowEdge, (int(highZ)-int(lowZ))/zspace, zspace );
        }
    }
    else if( "DistToTarg" == var )
    {
        const double targPos = GetBGZ( targetID );
        
        //skip this many modules before getting to the first (last) active faux module
        //mod 85 is the most DS module included in module z bins
        //mod 80 is the most DS module of the most DS faux target
        const int firstPlaneCount = 6;
        
        int counter = 0;
        int planesSinceLast = 0; //trick loop into making an active bin at the first chance
        
        std::vector<double> planeBins = GetPosBins( "vtxz", -1 ); //get bins by plane
        for( std::vector<double>::reverse_iterator i = planeBins.rbegin(); i != planeBins.rend(); ++i, ++counter )
        {
            if( counter < firstPlaneCount )
                continue;
            
            //faux targets every 6 modules
            if( 0 == planesSinceLast % 6 )
                binsLowEdge.push_back( targPos - *i );
            
            ++planesSinceLast;
        }
    }
    else if( var == "NTracks" )
    {
        // 1-10
        AddBins( binsLowEdge, 10, 1. );
    }
    else if( var == "plane" || var == "planeDNN"  )
    {
        // 0-66
        double tmpBins[] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 11, 12, 13., 14., 15., 16., 17., 18., 19., 21., 22., 23., 24., 25., 26., 27., 28., 29., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 51., 52., 53., 54., 55., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66. };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    
    return binsLowEdge;
}

double NukeCC_Binning::GetBGZ( const int targetID ) const
{
    //all in cm
    if( targetID == 1 )  return 452.4;
    if( targetID == 2 )  return 474.5 ;
    if( targetID == 3 )  return 501.05;
    if( targetID == 4 )  return 568.75;
    if( targetID == 5 )  return 582.05;
    
    Error( "GetBGZ", Form("Cannot get background Z for targetID = %d", targetID ) );
    throw 1;
}

vector<double> NukeCC_Binning::GetTruthBins(const std::string& var) const{
    
    vector<double> binsLowEdge;
    
    if(var == "Enu" || var == "Emu" || var == "Emum" ){
        AddBins( binsLowEdge, 5, 1. );
        //2 GeV bin up to 7 GeV
        /*AddBins( binsLowEdge, 1, 2. );
         //3 GeV bin up to 10 GeV
         AddBins( binsLowEdge, 1, 3. );*/
        //5GeV bins up to 30GeV
        AddBins( binsLowEdge, 5, 5. );
        //10GeV Bins up to 50 GeV
        AddBins( binsLowEdge, 2, 10 );
        AddBins( binsLowEdge, 1, 50. ); //50-100
        //Overflow
        AddBins( binsLowEdge, 1, 100);
    }
    
    else if(var == "Ehad")
        NukeCC_Binning::Get().AddBins(binsLowEdge, 100, 0.5);
    
    else if(var == "Q2"){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 200, 0.1);
        
    }
    
    else if(var == "x" ){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 110, 0.1);
        
    }
    
    else if(var == "xGen" ){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 110, 0.1);
        
    }
    
    else if(var == "y"){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 110, 0.01);
        
    }
    
    else if(var == "W"){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 200, 0.1);
        
    }
    
    else if(var == "ThetaMu"){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 180, 1.0);
        
    }
    
    else if (var == "NTracks"){
        
        AddBins( binsLowEdge, 10, 1. );
    }
    
    else if (var == "plane" || var == "planeDNN" ){
        // 0-66
        double tmpBins[] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 11, 12, 13., 14., 15., 16., 17., 18., 19., 21., 22., 23., 24., 25., 26., 27., 28., 29., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 51., 52., 53., 54., 55., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66. };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    
    else{
        Error( "NukeCC_Binning::GetTruthBins", Form("No bins found for variable %s", var.c_str() ) );
    }
    
    return binsLowEdge;
    
    
}//End of GetTruthBins

std::string NukeCC_Binning::GetXaxisTitle( const std::string& var ) const
{
    if( var == "Emu" )     return "Emu(GeV)";
    else if( var == "Ehad" )    return "Ehad (GeV)";
    else if( var == "CCQE-Recoil" )   return "CCQE-style Recoil (GeV)";
    else if( var == "VtxEnergy" )     return "Vertex Energy (GeV)";
    else if( var == "Enu" )     return "Enu (GeV)";
    else if( var == "Q2" )      return "Q^{2} (GeV/c)^{2}";
    else if( var == "W" )      return "W (GeV/c^{2})";
    else if( var == "y" )       return "y";
    else if( var == "x" )       return " x";
    else if( var == "PhiMu" )   return "Muon Phi wrt Beam (deg)";
    else if( var == "CosThetaMu" ) return "Cos(Muon Theta wrt Beam)";
    else if( var == "ThetaMu" ) return "Muon Theta wrt Beam (deg)";
    else if( var == "muonPt" ) return "Muon Pt(GeV)";
    else if( var == "muonPz" ) return "Muon Pz(GeV)";
    else if( var == "ThetaMuX" )return "Muon Theta_{X} wrt Beam (deg)";
    else if( var == "ThetaMuY" )return "Muon Theta_{Y} wrt Beam (deg)";
    else if( var == "Emum" )    return "Muon Energy at MINOS (GeV)";
    else if( var == "ETheta" )  return "E_{\\mu}*(1-cos(\\theta_{\\mu}))";
    
    else if( string::npos != var.find("vtxt") )   return "Vertex Time (ns)";
    else if( string::npos != var.find("vtxx") )   return "Vertex X (cm)";
    else if( string::npos != var.find("vtxy") )   return "Vertex Y (cm)";
    else if( string::npos != var.find("vtxz") )   return "Vertex Z (cm)";
    else if( string::npos != var.find("vtxu") )   return "Vertex U (cm)";
    else if( string::npos != var.find("vtxd") )   return "Vertex D (cm)";
    else if( string::npos != var.find("vtxc") )   return "Vertex C (cm)";
    
    else if( var == "plane" ) return "Plane Number Track-Based";
    else if( var == "planeDNN" )   return "Plane Number DNN";
    else if( var == "NTracks" )   return "Number of Tracks";
    else if( var == "segment" )   return "Segment Number DNN";
    
    cout << "ERROR [CVUniverse::GetXaxisTitle] : Could not get title for variable: " << var << endl;
    throw 1;
}

std::string NukeCC_Binning::GetYaxisTitle( const std::string& var ) const
{
    double normWidth = GetVarNormWidth( var );
    char words[256];
    if( var == "Emu" )          sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "Ehad" )    sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "CCQE-Recoil" )    sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "VtxEnergy" )      sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "Enu" )     sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "muonPt" )     sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "muonPz" )     sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "Emum" )    sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "ETheta" )  sprintf(words, "N Events / %.01f MeV", normWidth*1000. );
    else if( var == "Q2" )      sprintf(words, "N Events / %.1f (GeV/c)^{2}", normWidth );
    else if( var == "W" )       sprintf(words, "N Events / %.1f (GeV/c^{2})", normWidth );
    //else if( var == "x" )       sprintf(words, "N Events / %.0f ", normWidth ); //norm 1?
    else if( var == "x" )       sprintf(words, "N Events / %.1f ", normWidth );
    else if( var == "y" )       sprintf(words, "N Events / %.2f ", normWidth );
    
    else if( var == "PhiMu" )   sprintf(words, "N Event / %.0f deg", normWidth );
    
    else if( var == "ThetaMu"
            ||   var == "CosThetaMu"
            ||   var == "ThetaMuX"
            ||   var == "ThetaMuY" ) sprintf(words, "N Events / %.1f deg", normWidth );
    
    else if( string::npos != var.find("vtxt") ) sprintf(words, "N Events / %.0f ns", normWidth );
    else if( string::npos != var.find("vtxx") ) sprintf(words, "N Events / %.0f cm", normWidth );
    else if( string::npos != var.find("vtxy") ) sprintf(words, "N Events / %.0f cm", normWidth );
    else if( string::npos != var.find("vtxz") ) sprintf(words, "N Events / %.1f cm", normWidth );
    else if( string::npos != var.find("vtxu") ) sprintf(words, "N Events / %.0f cm", normWidth );
    else if( string::npos != var.find("vtxd") ) sprintf(words, "N Events / %.0f cm", normWidth );
    else if( string::npos != var.find("vtxc") ) sprintf(words, "N Events / %.0f cm", normWidth );
    
    else if( var == "plane" || var == "planeDNN" ) sprintf(words, "N Events" );
    else if( var == "NTracks") sprintf(words, "N Events" );
    else if( var == "segment") sprintf(words, "N Events" );
    
    else
    {
        cout << "ERROR [CVUnivers::GetYaxisTitle] : Could not get title for variable: " << var << endl;
        throw 1;
    }
   
    return string(words);
}


double NukeCC_Binning::GetVarNormWidth( const std::string& var ) const
{
    if( var == "Emu" )    return 1.0;
    if( var == "Ehad" )   return .4;
    if( var == "CCQE-Recoil" )   return .25;
    if( var == "VtxEnergy" )     return .25;
    if( var == "Enu" )    return 1.0;
    if( var == "muonPt" )    return 1.0;
    if( var == "muonPz" ) return 1.0;
    if( var == "Emum" )   return 1.0;
    if( var == "ETheta" ) return .02;
    if( var == "Q2" )     return .1;
    if( var == "W" )     return 1.0;
    //if( var == "y" )      return .1;//?
    if( var == "y" )      return 1;//?
    //if( var == "x" )      return 1.;//?//if norm 1
    if( var == "x" )      return .1;//?
    
    if( var == "PhiMu" )  return 6.;
    if( var == "ThetaMu" || var == "ThetaMuX" || var == "ThetaMuY" ) return  1.;
    if( var == "CosThetaMu" )  return 1.;
    
    if( string::npos != var.find("vtxt") )      return 100.;
    if( string::npos != var.find("vtxx") )      return 10.;
    if( string::npos != var.find("vtxy") )      return 10.;
    if( string::npos != var.find("vtxz") )      return 1.7; //if it uses module/plane bins then you need to set the norm width to -1 yourself
    if( string::npos != var.find("vtxd") )      return 10.;
    if( string::npos != var.find("vtxu") )      return 10.;
    if( string::npos != var.find("vtxc") )      return 10.;
    
    if( var == "plane" || var == "planeDNN"  ) return 1.;
    if( var == "NTracks" || var == "segment" )   return 1.;
    
    
    Error("GetVarNormWidth", Form("Could not get the normalization width for variable: %s", var.c_str() ) );
    throw 1;
}

double NukeCC_Binning::GetDISVarMinVal( const std::string& var ) const
{
    if( var == "Emu" )    return 2.0;
    if( var == "Ehad" )   return 1.0;
    if( var == "CCQE-Recoil" )   return 0.;
    if( var == "VtxEnergy" )     return 0.;
    if( var == "Enu" )    return 5.0;
    if( var == "Emum" )   return 2.0;
    if( var == "ETheta" ) return 0.;
    if( var == "Q2" )     return 1.0;
    if( var == "W" )      return 2.0;
    if( var == "y" )      return 0.;
    if( var == "x" )      return 0.005;
    
    if( var == "PhiMu" ) return -180.;
    
    if( var == "ThetaMu" ) return 0.;
    if( var == "ThetaMuX" || var == "ThetaMuY" ) return  -20.;
    
    if( var == "plane" || var == "planeDNN"  ) return 0;
    if( var == "NTracks" ) return 0;
    if( var == "segment" ) return 0.;
    
    if( var == "vtxx" || var == "vtxy" ) return -100.;
    if( var == "vtxz" ) return 400.0;
    
    Error("GetVarMinVal", Form("Could not get the min value to plot for variable: %s", var.c_str() ) );
    throw 1;
}
double NukeCC_Binning::GetDISVarMaxVal( const std::string& var ) const
{
    const double e = EPSILON;
    if( var == "Emu" )    return 25.-e;
    if( var == "Ehad" )   return 20.-e;
    if( var == "CCQE-Recoil" )   return 3.-e;
    if( var == "VtxEnergy" )     return 3.-e;
    if( var == "Enu" )    return 50.-e;
    if( var == "Emum" )   return 50.-e;
    if( var == "ETheta" ) return .5-e;
    if( var == "Q2" )     return 8.-e;
    if( var == "W" )      return 8.-e;
    if( var == "y" )      return 1.-e;
    if( var == "x" )      return 0.8 - e;
    
    if( var == "PhiMu" ) return 180.-e;
    if( var == "plane" || var == "planeDNN"  ) return 68-e;
    if( var == "NTracks" ) return 10-e;
    if( var == "segment" ) return 174.;
    
    if( var == "ThetaMu" ) return 17.-e; //was 24
    if( var == "ThetaMuX" || var == "ThetaMuY" ) return  20.-e;
    
    if( var == "vtxx" || var == "vtxy" ) return 100.-e;
    if( var == "vtxz" ) return 850.-e;
    
    Error("GetVarMaxVal", Form("Could not get the max value to plot for variable: %s", var.c_str() ) );
    throw 1;
}
#endif 

