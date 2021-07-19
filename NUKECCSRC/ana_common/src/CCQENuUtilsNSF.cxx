#ifndef MNV_CCQENUUTILSNSF_cxx
#define MNV_CCQENUUTILSNSF_cxx 1

//#include "CCQENuUtilsNSF.h"
#include "include/CCQENuUtilsNSF.h"

//#include "../include/NukeCC_Cuts.h"
//#include "CCQENuUtilsNSF.h"
#include "include/CVUniverse.h"
#include "include/GlobalIncludes.h" 
#include "include/LateralSystematics.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"

using namespace NUKECC_ANA; 

//--Default constructor
/*
CCQENuUtilsNSF::CCQENuUtilsNSF(  bool useFluxConstraint,string playlist )  
{ 
   cvhistos1D.clear();
   cvhistos2D.clear();
   
   //cutter = new CVUniverse();
  // GlobalParameters::Get().m_useFluxConstraint = useFluxConstraint;
  DefaultCVUniverse::SetNuEConstraint(useFluxConstraint);
  DefaultCVUniverse::SetPlaylist(playlist);
  // GlobalParameters::Get().m_usePPFX1Flux = false;
  // use_merged_files = use_mergedFiles;
   
   if(neutrinoMode){

   incoming_pdg = 14;
   cout<<"Initializing CCQENuUtils in Neutrino Mode"<<endl;
   }
   else{
   incoming_pdg = -14;
   cout<<"Initializing CCQENuUtils in Antineutrino Mode"<<endl;
   
   }
  DefaultCVUniverse::SetAnalysisNuPDG(incoming_pdg);
  
  if( neutrinoMode ){
    //Use neutrino weights for general constructor
    //Vector declared in header file
    michel_weights.push_back(0.986);
    michel_weights.push_back(1.013);
  } else{
    //For antineutrino, set weights to 0. at present 
    michel_weights.push_back(1.0);
    michel_weights.push_back(1.0);
  }   
   
 }
*/
 
CCQENuUtilsNSF::CCQENuUtilsNSF(string playlist){
    cvhistos1D.clear();
    cvhistos2D.clear();
     cutter = new NukeCC_Cuts();
//   GlobalParameters::Get().m_useFluxConstraint = false;
   DefaultCVUniverse::SetPlaylist(playlist);
  // GlobalParameters::Get().m_usePPFX1Flux = false;
   use_merged_files = false;
   
   if(neutrinoMode){
   incoming_pdg = 14;
   cout<<"Initializing CCQENuUtils in Neutrino Mode"<<endl;
   }
   else{
   incoming_pdg = -14;
   cout<<"Initializing CCQENuUtils in Antineutrino Mode"<<endl;
   
   }
   DefaultCVUniverse::SetAnalysisNuPDG(incoming_pdg);
  if( neutrinoMode ){
    //Use neutrino weights for general constructor
    //Vector declared in header file
    michel_weights.push_back(0.986);
    michel_weights.push_back(1.013);
  } else{
    //For antineutrino, set weights to 0. at present 
    michel_weights.push_back(1.0);
    michel_weights.push_back(1.0);
  }    

}


CCQENuUtilsNSF::CCQENuUtilsNSF()
{
    cvhistos1D.clear();
    cvhistos2D.clear();
     
     cutter = new NukeCC_Cuts();
 //   cutter = new CVUniverse();
  // GlobalParameters::Get().m_useFluxConstraint = false;
  // DefaultCVUniverse::SetPlaylist("minervame1A");
  // GlobalParameters::Get().m_usePPFX1Flux = false;
   use_merged_files = false;
    
   if(neutrinoMode){
   incoming_pdg = 14;
   cout<<"Initializing CCQENuUtils in Neutrino Mode"<<endl;
   }
   else{
   incoming_pdg = -14;
   cout<<"Initializing CCQENuUtils in Antineutrino Mode"<<endl;
   
   }
   DefaultCVUniverse::SetAnalysisNuPDG(incoming_pdg);
  if( neutrinoMode ){
    //Use neutrino weights for general constructor
    //Vector declared in header file
    michel_weights.push_back(0.986);
    michel_weights.push_back(1.013);
  } else{
    //For antineutrino, set weights to 0. at present 
    michel_weights.push_back(1.0);
    michel_weights.push_back(1.0);
  }    

}

 //--Destructor
 CCQENuUtilsNSF::~CCQENuUtilsNSF()
 {
    delete cutter;
     cvhistos1D.clear();
     cvhistos2D.clear();
 }



std::vector<std::string> CCQENuUtilsNSF::GetStdPlaylists( HelicityType::t_HelicityType helicity ) const
{
    std::vector<std::string> playlists;
    if( HelicityType::kAntiNeutrino != helicity )
    {
        playlists.push_back( "minerva1" );
        playlists.push_back( "minerva7" );
        playlists.push_back( "minerva9" );
                   playlists.push_back( "minerva13C" );
                                        }
        //                                    
          if( HelicityType::kNeutrino != helicity )
          {
        playlists.push_back( "minervame1A" );

}
    
    return playlists;
}

HelicityType::t_HelicityType CCQENuUtilsNSF::GetHelicityFromPlaylist( std::string playlist ) const
{
 //   if ("" == playlist)
 //       playlist = fPlaylist;
    
        if(Playlist::minervame1A == playlist ||
           Playlist::minervame1B == playlist ) {
        return HelicityType::kNeutrino;
    }
       if (Playlist::minervame5A == playlist ||
        Playlist::minervame6A == playlist
        ) {
        return HelicityType::kAntiNeutrino;
    }
    
    const std::vector<std::string> nus = GetStdPlaylists( HelicityType::kNeutrino );
    if ( find( nus.begin(), nus.end(), playlist ) != nus.end() )
        return HelicityType::kNeutrino;
    
    const std::vector<std::string> antinus = GetStdPlaylists( HelicityType::kAntiNeutrino );
    if ( find( antinus.begin(), antinus.end(), playlist ) != antinus.end() )
        return HelicityType::kAntiNeutrino;
    
    Warning( "CCQENuUtilsNSF::GetHelicityFromPlaylist", Form( "Playlist '%s' is not a standard nu or antinu playlist. Using any neutrino.", playlist.c_str() ) );
    return HelicityType::kAny;
}

TString CCQENuUtilsNSF::GetHistFileName( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ, HelicityType::t_HelicityType helicity, bool forceDisk, bool usebluearc ) const
{
    TString histFileName = HistDir(forceDisk,false,usebluearc);
    
    TString variation = GetVariationTag();
    if( FileType::kAny == fType )
    {
        histFileName += Form("/Hists_%s_t%d_z%02d_%s_%s_%s.root",
                             histType.c_str(),
                             targetID, targetZ,
                             GetHelicityString( helicity ).c_str(),
                             getenv("NUKECC_TAG"),variation.Data());
    }
    else
    {
        histFileName += Form("/Hists_%s_%s_t%d_z%02d_%s_%s_%s.root",
                             histType.c_str(),
                             GetFileTypeString(fType).c_str(),
                             targetID, targetZ,
                             GetHelicityString( helicity ).c_str(),
                             getenv("NUKECC_TAG"),variation.Data());
                             
    }
    return histFileName; 
}


TString CCQENuUtilsNSF::HistDir( bool forceDisk /* = false */, bool ignorePlaylist /*= false*/, bool needsbluearc /*= false*/ ) const
{
    TString CONDOR_DIR_HISTS = getenv("CONDOR_DIR_HISTS");
    if( ! ( CONDOR_DIR_HISTS.IsNull() || forceDisk ) )
    {
        return TString( CONDOR_DIR_HISTS + "/");
    }
    else
    {
        
        TString HISTS_ROOT;
        
        if(!needsbluearc) {
            HISTS_ROOT = getenv("HISTS");
        }
        if(needsbluearc){
            HISTS_ROOT = getenv("HISTS_BLUEARC");
        }
        cout<<" The hists directory is "<<HISTS_ROOT<<endl;
        if( HISTS_ROOT.IsNull() )
        {
            Error( "CCQENuUtilsNSF::HistDir", "Environmental variable HISTS not defined." );
            throw 1;
        }
        
        //replace minervagli with GRID_USER
        const TString minervagli = "minervagli";
        int idx = HISTS_ROOT.Index(minervagli);
        if( idx != kNPOS )
        {
            TString GRID_USER = getenv("GRID_USER");
            HISTS_ROOT.Replace(idx, idx + minervagli.Length(), GRID_USER );
            cout << "  REPLACE minervagli with GRID_USER = " << GRID_USER << ", HISTS_ROOT = " << HISTS_ROOT << endl;
        }
        
        TString NUKECC_TAG  = getenv("NUKECC_TAG");
        if( NUKECC_TAG.IsNull() )
        {
            Error( "CCQENuUtilsNSF::HistDir", "Environmental variable NUKECC_TAG not defined." );
            throw 1;
        }
        
        TString histdir(  HISTS_ROOT + "/" + NUKECC_TAG + "/"  );
       // if( fPlaylist.empty() || ignorePlaylist )
         //   histdir = TString( HISTS_ROOT + "/" + NUKECC_TAG + "/" );
        
        //cout << system( Form( "test -d %s", histdir.Data() ) ) <<endl;
        if( 0 != system( Form( "test -d %s", histdir.Data() ) ) )
        {
            int madedir = system( Form( "mkdir -m 755 -p %s", histdir.Data() ) );
            
            if( 0 != madedir )
                Error( "CCQENuUtilsNSF::HistDir", Form("Could not make hist directory '%s'", histdir.Data() ) );
        }
        
        return histdir;
    }
}


TString CCQENuUtilsNSF::GetVariationTag() const
{
    TString var = "";
  /* if( fSysShiftTag.empty() )
    {
        if( USE_INEL_CUT )
            var += "_inel";
        
        if( 0 != VTX_BLOB_SHIFT )
            var += Form("_VtxBlobShift%d", VTX_BLOB_SHIFT );
     */   return var;
  //  }
   // else 
  //  {
    //    var += "_" + fSysShiftTag;
  //  }
    
   // return TString( fSysShiftTag.c_str() );
}

std::string CCQENuUtilsNSF::GetHelicityString( HelicityType::t_HelicityType helicity ) const
{
    if(      HelicityType::kNeutrino     == helicity ) return "Nu";
    else if( HelicityType::kAntiNeutrino == helicity ) return "AntiNu";
    else return "AnyHelicity";
}
TString CCQENuUtilsNSF::GetHistName( const std::string& histType, FileType::t_FileType fType, const std::string& var, int targetID, int targetZ ) const
{
    TString histName;
    
    if( FileType::kAny == fType )
    {
        histName += Form( "%s_%s_t%d_z%02d",
                         histType.c_str(),
                         var.c_str(),
                         targetID,
                         targetZ
                         );
    }
    else
    {
        histName += Form( "%s_%s_%s_t%d_z%02d",
                         histType.c_str(),
                         GetFileTypeString(fType).c_str(),
                         var.c_str(),
                         targetID,
                         targetZ
                         );
    }
    return histName;
}

std::string CCQENuUtilsNSF::GetFileTypeString( FileType::t_FileType fType ) const
{
    if(      FileType::kData  == fType ) return "Data";
    else if( FileType::kMC    == fType ) return "MC";
    else if( FileType::kTruth == fType ) return "Truth";
    else if( FileType::kAny   == fType ) return "Any";
    else if( FileType::kDNNData == fType ) return "DNNData";
    else if( FileType::kDNNMC == fType ) return "DNNMC";
    else if( FileType::kDCNNData == fType ) return "DCNNData";
    else if( FileType::kDCNNMC == fType ) return "DCNNMC";
    else if( FileType::kDANNData == fType ) return "DANNData";
    else if( FileType::kDANNMC == fType ) return "DANNMC";
    else if( FileType::kNukeOnlyMC == fType ) return "NukeOnlyMC";
    else if( FileType::kNukeOnlyTruth == fType ) return "NukeOnlyTruth";
    else throw( "Unknown FileType found in GetFileTypeString" );
}

 //============================================================================
 // bookHistosCV
 //===========================================================================
 void CCQENuUtilsNSF::bookHistosCV( HistWrapper<CVUniverse>** h, string var_name, string title, std::vector<double> nbin, double min, double max, std::map< std::string, std::vector<CVUniverse*> > error_bands ){
   MnvH1D *hist[nHistosMy];
   //if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
     for( unsigned int i = 0; i < nHistosMy; ++i ){
       hist[i] = new MnvH1D( ( var_name + '_' + namesMy[i] ). c_str(),title.c_str(),
			  nbin.size(), min, max ); 
       h[i] = new HistWrapper<CVUniverse>( hist[i],error_bands );
       
     }
  // }
  // else{
    // for( unsigned int i = 0; i < nHistosInc; ++i )
      // h[i] = new HistWrapper<CVUniverse>( hist[i],error_bands );
  // }
    cvhistos1D[var_name] = h;
 }
/*
 void CCQENuUtilsNSF::bookHistosCV( HistWrapper<CVUniverse>** h, string var_name, string title, axis_binning xbins, std::map< std::string, std::vector<CVUniverse*> > error_bands ){
  MnvH1D *hist[nHistos];
  // if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
     for( unsigned int i = 0; i < nHistos; ++i ){
       
       hist[i] = new  MnvH1D( ( var_name + '_' + names[i] ).c_str(), title.c_str(),
			   xbins.nbins, &( xbins.bin_edges[0] ), xbins.default_width );
			   
       h[i] = new HistWrapper<CVUniverse>(hist[i],error_bands);
	   
	}		   
  // }
   else{
     for( unsigned int i = 0; i < nHistosInc; ++i ){
       hist[i] = new  MnvH1D( ( var_name + '_' + namesInc[i] ).c_str(), title.c_str(),
			   xbins.nbins, &( xbins.bin_edges[0] ), xbins.default_width);
       h[i] = new HistWrapper<CVUniverse>(hist[i],error_bands)	;   
			   
      }
   }
    cvhistos1D[var_name] = h;
 }*/
/*
  void CCQENuUtilsNSF::bookHistosCV( Hist2DWrapper<CVUniverse>** h, string var_name, string title, int xnbins, double xmin, double xmax, int ynbins, double ymin, double ymax, std::map< std::string, std::vector<CVUniverse*> > error_bands ){
   MnvH2D *hist[nHistos]; 
  // if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
     for( unsigned int i = 0; i < nHistos; ++i ){
       hist[i] = new MnvH2D( ( var_name + '_' + names[i] ). c_str(),title.c_str			  xnbins, xmin, xmax, ynbins, ymin, ymax);
			  xnbins, xmin, xmax, ynbins, ymin, ymax);
      h[i] = new Hist2DWrapper<CVUniverse>(hist[i],error_bands);	
      
      }		 
  // }
   else{
     for( unsigned int i = 0; i < nHistosInc; ++i ){
       hist[i] = new MnvH2D( ( var_name + '_' + namesInc[i] ). c_str(),title.c_str(),
			  xnbins, xmin, xmax, ynbins, ymin, ymax);
	
       h[i] = new Hist2DWrapper<CVUniverse>(hist[i],error_bands);
   }

    // }
    cvhistos2D[var_name] = h;
 }*/
 /*
 void CCQENuUtilsNSF::bookHistosCV( Hist2DWrapper<CVUniverse>** h,string var_name,string title,axis_binning xbins, axis_binning ybins, std::map< std::string, std::vector<CVUniverse*> > error_bands ){
   Double_t xbinsarray[xbins.nbins];
   MnvH2D *hist[nHistos];
   for( unsigned int i = 0; i <= xbins.nbins; i++ )
     xbinsarray[i] = xbins.bin_edges[i];
   Double_t ybinsarray[ybins.nbins];
   for( unsigned int i = 0; i <= ybins.nbins; i++ )
     ybinsarray[i] = ybins.bin_edges[i];
//   if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
     for( unsigned int i = 0; i < nHistos; i++ ){
       hist[i] = new MnvH2D((var_name+Form("_%s",names[i].c_str())).c_str(),title.c_str(),xbins.nbins,xbinsarray,ybins.nbins,ybinsarray);
       
       h[i] = new Hist2DWrapper<CVUniverse>(hist[i],error_bands); 
       
     }  
  // }
 
  else{
     for( unsigned int i = 0; i < nHistosInc; i++ ){
       hist[i] = new MnvH2D((var_name+Form("_%s",namesInc[i].c_str())).c_str(),title.c_str(),xbins.nbins,xbinsarray,ybins.nbins,ybinsarray);
       h[i] = new Hist2DWrapper<CVUniverse>(hist[i],error_bands);
   }
     }
    cvhistos2D[var_name] = h;
 }
*/

//I will only have the fill for 1D and 2D histograms for now....

 void CCQENuUtilsNSF::fillHistos( HistWrapper<CVUniverse>**h, double var, CVUniverse* universe, double w ){  // fill mc
  // if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
 //   NukeCC_Cuts* cutter = new NukeCC_Cuts();
//cout<<"Huma I ams here"<<endl;
        h[kIncCCMy]->univHist(universe)->Fill( var, w );
//cout<<"Huma I ams here 1"<<endl;   
      if(universe->ShortName()=="cv"){ //only fill the CV histogram
      // if(cutter->PassLowWCut(universe))h[kIncCCMy]->univHist(universe)->Fill( var, w );

//cout<<"Huma I ams here 2"<<endl;   
//     if( cutter->PassMuEnergyCut( universe ) )    h[kIncQEMy]->univHist(universe)->Fill(var, w );

//cout<<"Huma I ams here 3"<<endl;   
       if( cutter->PassDISCut( universe ) )  h[kIncDISMy]->univHist(universe)->Fill( var, w );

//cout<<"Huma I ams here 3a"<<endl;   
     if( cutter->PassLowWCut( universe ) )  h[kInclowW]->univHist(universe)->Fill( var, w );

//cout<<"Huma I ams here 4"<<endl;   
   //  if( !cutter->passTrueCCQE( universe ) && !cutter->passTrueCCRES( universe ) && !cutter->passTrueCCDIS( universe ) && !cutter->passTrueMEC( universe ) )
     //  h[kOTH]->univHist(universe)->Fill( var, w);
     
//cout<<"Huma I ams here 5"<<endl;   
      if( cutter->PassLowQ2Cut( universe ) ) h[kInclowQ2]->univHist(universe)->Fill( var, w );
     
//cout<<"Huma I ams here 6"<<endl;   
   }

 }

/*
void CCQENuUtilsNSF::fillHistos( Hist2DWrapper<CVUniverse>**h, double var_x, double var_y,CVUniverse* universe, double w ){  // fill mc
  //if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
    //std::cout<<"Filling histograms "<<var_x<<" "<<var_y<<" "<<w<<std::endl;
    h[kMC]->univHist(universe)->Fill( var_x, var_y, w );
    
    if( universe->PassMuEnergyCut()  )  h[kQE]->univHist(universe)->Fill( var_x, var_y, w );

    else                                  h[kQENot]->univHist(universe)->Fill( var_x, var_y, w );

    if( cutter->passTrueMEC( universe ) )  h[k2p2h]->univHist(universe)->Fill( var_x, var_y, w );
    if( cutter->passTrueCCRES( universe ) )  h[kRES]->univHist(universe)->Fill( var_x, var_y, w );
    if( cutter->passTrueCCDIS( universe ) )  h[kDIS]->univHist(universe)->Fill( var_x, var_y, w );
    if( !cutter->passTrueCCQE( universe ) && !cutter->passTrueCCRES( universe ) && !cutter->passTrueCCDIS( universe ) && !cutter->passTrueMEC( universe ) )
      h[kOTH]->univHist(universe)->Fill( var_x, var_y, w);
    if( cutter->passTrueCCQELike( universe ) ) h[kQELike]->univHist(universe)->Fill( var_x, var_y, w );
    else                                    h[kQELikeNot]->univHist(universe)->Fill( var_x, var_y, w );

    if( cutter->passTrueCCQELike( universe ) && cutter->passTrueCCQE( universe )  ){
      h[kQELike_QE]->univHist(universe)->Fill( var_x, var_y, w );
      if( cutter->passTrueFSProton( universe ) ) h[kQELike_QE_ProtonFSI]->univHist(universe)->Fill( var_x, var_y, w );
      else h[kQELike_QE_NeutronFSI]->univHist(universe)->Fill( var_x, var_y, w );
    }
    if( cutter->passTrueCCQELike( universe ) && cutter->passTrueMEC( universe )  ){
      h[kQELike_2p2h]->univHist(universe)->Fill( var_x, var_y, w );
      if(cutter->passMECnp( universe ) )h[kQELike_2p2h_np]->univHist(universe)->Fill( var_x, var_y, w );
      if(cutter->passMECnn( universe ) ) h[kQELike_2p2h_nn]->univHist(universe)->Fill( var_x, var_y, w );
    }
    if( cutter->passTrueCCQELike( universe ) && cutter->passTrueCCRES( universe ) ){
      h[kQELike_RES]->univHist(universe)->Fill( var_x, var_y, w );
      
      if( cutter->passTrueFSProton( universe ) ) h[kQELike_RES_ProtonFSI]->univHist(universe)->Fill( var_x, var_y, w );
      else h[kQELike_RES_NeutronFSI]->univHist(universe)->Fill( var_x, var_y, w );
      if( cutter->passTrueISProton( universe ) ) h[kQELike_RES_ProtonIS]->univHist(universe)->Fill( var_x, var_y, w );
      else h[kQELike_RES_NeutronIS]->univHist(universe)->Fill( var_x, var_y, w );
    }
    if( cutter->passTrueCCQELike( universe ) && cutter->passTrueCCDIS( universe ) )  h[kQELike_DIS]->univHist(universe)->Fill( var_x, var_y, w );
    if( cutter->passTrueCCQELike( universe ) && (!cutter->passTrueCCQE( universe ) && !cutter->passTrueCCRES( universe ) && !cutter->passTrueCCDIS( universe ) &&!cutter->passTrueMEC( universe ) ) )
      h[kQELike_OTH]->univHist(universe)->Fill( var_x, var_y, w );
  
    if( !cutter->passTrueCCQELike( universe ) && cutter->passTruePositivePion( universe ) ){
      h[kQELikeNot_PositivePions]->univHist(universe)->Fill( var_x, var_y, w );
      if( universe->GetBool("truth_has_michel_electron")==true )
	h[kQELikeNot_PositivePions_TrueMichel]->univHist(universe)->Fill( var_x, var_y, w );
      else 
	h[kQELikeNot_PositivePions_TrueMichelNot]->univHist(universe)->Fill( var_x, var_y, w );
    }
    if( !cutter->passTrueCCQELike( universe ) && cutter->passTrueNegativePion( universe ) )  h[kQELikeNot_NegativePions]->univHist(universe)->Fill( var_x, var_y, w );
    if( !cutter->passTrueCCQELike( universe ) && cutter->passTrueNeutralPion( universe ) )   h[kQELikeNot_NeutralPions]->univHist(universe)->Fill( var_x, var_y, w );
    if( !cutter->passTrueCCQELike( universe ) && !cutter->passTruePion( universe ) )         h[kQELikeNot_NoPions]->univHist(universe)->Fill( var_x, var_y, w );

    if( !cutter->passTrueCCQELike( universe ) && (cutter->passTruePositivePion( universe ) || cutter->passTrueNegativePion( universe )) ){
      h[kQELikeNot_ChargedPions]->univHist(universe)->Fill( var_x, var_y, w );
      if( universe->GetBool("truth_has_michel_electron")==true )
	h[kQELikeNot_ChargedPions_TrueMichel]->univHist(universe)->Fill( var_x, var_y, w );
      else 
	h[kQELikeNot_ChargedPions_TrueMichelNot]->univHist(universe)->Fill( var_x, var_y, w );
    }
    if( !cutter->passTrueCCQELike( universe ) &&  !cutter->passTruePositivePion( universe ) && !cutter->passTrueNegativePion( universe ) )
      h[kQELikeNot_ChargedPionsNot]->univHist(universe)->Fill( var_x, var_y, w );
    
    if( !cutter->passTrueCCQE( universe ) && cutter->passTruePositivePion( universe ) ){
      h[kQENot_PositivePions]->univHist(universe)->Fill( var_x, var_y, w );
      if( universe->GetBool("truth_has_michel_electron")==true )
	h[kQENot_PositivePions_TrueMichel]->univHist(universe)->Fill( var_x, var_y, w );
      else 
	h[kQENot_PositivePions_TrueMichelNot]->univHist(universe)->Fill( var_x, var_y, w );
    }
    
    if( !cutter->passTrueCCQE( universe ) && cutter->passTrueNegativePion( universe ) )  h[kQENot_NegativePions]->univHist(universe)->Fill( var_x, var_y, w );
    if( !cutter->passTrueCCQE( universe ) && cutter->passTrueNeutralPion( universe ) )   h[kQENot_NeutralPions]->univHist(universe)->Fill( var_x, var_y, w );
    if( !cutter->passTrueCCQE( universe ) && !cutter->passTruePion( universe ) )         h[kQENot_NoPions]->univHist(universe)->Fill( var_x, var_y, w );
  
    if( cutter->passTrueCCQELike( universe ) && !cutter->passTrueCCQE( universe ) )  h[kQELike_QENot]->univHist(universe)->Fill( var_x, var_y, w );
    if( !cutter->passTrueCCQELike( universe ) && cutter->passTrueCCQE( universe ) )  h[kQELikeNot_QE]->univHist(universe)->Fill( var_x, var_y, w );
    if( !cutter->passTrueCCQELike( universe ) && !cutter->passTrueCCQE( universe ) )  h[kQELikeNot_QENot]->univHist(universe)->Fill( var_x, var_y, w );

    if( !cutter->passTrueCCQELike( universe ) && cutter->passTrueSingleChargedPion( universe ) )  h[kQELikeNot_SingleChargedPion]->univHist(universe)->Fill( var_x, var_y, w );
    if( !cutter->passTrueCCQELike( universe ) && cutter->passTrueSingleNeutralPion( universe ) )  h[kQELikeNot_SingleNeutralPion]->univHist(universe)->Fill( var_x, var_y, w );
    if( !cutter->passTrueCCQELike( universe ) && cutter->passTrueMultiPion( universe ) )  h[kQELikeNot_MultiPion]->univHist(universe)->Fill( var_x, var_y, w );
    if( !cutter->passTrueCCQELike( universe ) && !cutter->passTrueSingleChargedPion( universe ) && !cutter->passTrueSingleNeutralPion( universe ) && !cutter->passTrueMultiPion( universe ) )  h[kQELikeNot_Other]->univHist(universe)->Fill( var_x, var_y, w );

    if( !cutter->passTrueCCQELike( universe ) && cutter->passTrueMEC( universe ) )  h[kQELikeNot_2p2h]->univHist(universe)->Fill( var_x, var_y, w );	 
    if( !cutter->passTrueCCQELike( universe ) && cutter->passTrueCCDIS( universe ) )  h[kQELikeNot_DIS]->univHist(universe)->Fill( var_x, var_y, w );	 
    if( !cutter->passTrueCCQELike( universe ) && cutter->passTrueCCRES( universe ) )  h[kQELikeNot_RES]->univHist(universe)->Fill( var_x, var_y, w );	 
    if( !cutter->passTrueCCQELike( universe ) && cutter->passTrueCoh( universe ) )  h[kQELikeNot_Coh]->univHist(universe)->Fill( var_x, var_y, w );
  }
  else{
    h[kIncMC]->univHist(universe)->Fill( var_x, var_y, w );
    if( cutter->passTrueCCQE( universe )  )  h[kIncQE]->univHist(universe)->Fill( var_x, var_y, w );
    if( cutter->passTrueMEC( universe ) )    h[kInc2p2h]->univHist(universe)->Fill(var_x, var_y, w );
    if( cutter->passTrueCCRES( universe ) )  h[kIncRES]->univHist(universe)->Fill( var_x, var_y, w );
    if( cutter->passTrueCCDIS( universe ) ){
      h[kIncDIS]->univHist(universe)->Fill( var_x, var_y, w );
      if(cutter->passTrueCCTrueDIS( universe ) ) h[kIncDIS_DIS]->univHist(universe)->Fill( var_x, var_y, w);
      if(cutter->passTrueCCTrueSIS( universe ) ) h[kIncDIS_SIS]->univHist(universe)->Fill( var_x, var_y, w);
    }

    if( !cutter->passTrueCCQE( universe ) && !cutter->passTrueCCRES( universe ) && !cutter->passTrueCCDIS( universe ) && !cutter->passTrueMEC( universe ) && cutter->passTrueCC( universe ) ){
      h[kIncOth]->univHist(universe)->Fill( var_x, var_y, w);
    }
    if( cutter->passTrueCC( universe )){
      h[kIncCC]->univHist(universe)->Fill( var_x, var_y, w);
    }
    else{
      h[kIncNC]->univHist(universe)->Fill( var_x, var_y, w);
    }
 // }
  
}*/

void CCQENuUtilsNSF::syncCVHistos(HistWrapper<CVUniverse>**h){
    //the sync hist function....
    for(unsigned int i=kIncCCMy;i<nHistosMy;++i)h[i]->SyncCVHistos();
}

/*
void CCQENuUtilsNSF::fillHistos( HistWrapper<CVUniverse>** h, double var ){  // fill data
 // if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
    h[kData]->hist->Fill( var );
  }
 // else{
   // h[kIncData]->hist->Fill( var );
 // }
// }

void CCQENuUtilsNSF::fillHistos( Hist2DWrapper<CVUniverse>** h, double var_x, double var_y ){  // fill data
  //if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
    h[kData]->hist->Fill( var_x,var_y );
  }
 // else{
   // h[kIncData]->hist->Fill( var_x,var_y );
 // }
// } 
 */
ChainWrapper* CCQENuUtilsNSF::GetChainWrapperMCPointer(string playlist,string tree_name){
  std::string mybase = getenv("MY_NSFNUKECC");
 // std::string mc_playlist_path = mybase + "/whole_"+playlist+"Anne_me1A.txt";
  std::string mc_playlist_path = mybase+"/include/playlists/NukeCC_"+playlist+"_MC_Inextinguishable_merged.txt";
// ChainWrapper* chain = makeChainWrapperPtr("../../../playlists/whole_Anne_me1A.txt","NukeCC");
  ChainWrapper* chain = makeChainWrapperPtr(mc_playlist_path,tree_name); 
 return chain;
} 
 
 
ChainWrapper* CCQENuUtilsNSF::GetChainWrapperDataPointer(string playlist,string tree_name){
  std::string mybase = getenv("MY_CCQENU");
  std::string data_playlist_path = mybase+"/include/playlists/CCQENu_"+playlist+"_DATA_Inextinguishable_merged.txt";
  ChainWrapper* chain = makeChainWrapperPtr(data_playlist_path, tree_name);
  return chain;
} 


void CCQENuUtilsNSF::getPOT(double& total_pot_data,double& total_pot_mc){
     ChainWrapper *chain_mc =  GetChainWrapperMCPointer(DefaultCVUniverse::GetPlaylist(),"Meta");
     ChainWrapper *chain_data = GetChainWrapperDataPointer(DefaultCVUniverse::GetPlaylist(),"Meta");
     total_pot_data = setPOTData(chain_data);
     total_pot_mc = setPOTMC(chain_mc);
     //delete chain_data;
     //delete chain_mc;
}

double CCQENuUtilsNSF::setPOTData(ChainWrapper* _chain){
   TChain *chain = _chain->GetChain();
   double tot_pot_data = 0.0;
   double pot_data = 0.0;
   TBranch *b_pot_data;
   chain->SetBranchAddress("POT_Used",&pot_data,&b_pot_data);
   Int_t tot_entries = chain->GetEntries();
   double global_pot_data = 0.0;
   double global_pot=0.0;
   TBranch *b_global_pot;
   chain->SetBranchAddress("POT_Total",&global_pot,&b_global_pot);
   for(int i =0;i<tot_entries;i++){
     chain->GetEntry(i);
     tot_pot_data+= pot_data;
     global_pot_data+= global_pot;
   
   }
   std::cout<<"Total Data POT: "<<global_pot_data<<std::endl;
   std::cout<<"Used Data POT: "<<tot_pot_data<<std::endl;
   delete chain;
   //delete b_global_pot;
   return tot_pot_data;
   
}

double CCQENuUtilsNSF::setPOTMC(ChainWrapper* _chain){
   TChain *chain = _chain->GetChain();
   double tot_pot_data = 0.0;
   double pot_data = 0.0;
   TBranch *b_pot_data;
   chain->SetBranchAddress("POT_Used",&pot_data,&b_pot_data);
   Int_t tot_entries = chain->GetEntries();
   double global_pot_data = 0.0;
   double global_pot=0.0;
   TBranch *b_global_pot;
   chain->SetBranchAddress("POT_Total",&global_pot,&b_global_pot);
   for(int i =0;i<tot_entries;i++){
     chain->GetEntry(i);
     tot_pot_data+= pot_data;
     global_pot_data+= global_pot;
   
   }
   std::cout<<"Total MC POT: "<<tot_pot_data<<std::endl;
   std::cout<<"Used MC POT: "<<global_pot_data<<std::endl;
   delete chain;
   //delete b_global_pot;
   return tot_pot_data;
   
}


void CCQENuUtilsNSF::writePOT( TFile *f ){
  
  double data = 0.0;
  double mc   = 0.0;
  this->getPOT(data,mc);
  TVector2 *pot = new TVector2( data, mc );
   f->WriteTObject( pot, "pot" );
   
 }

std::map<std::string,std::vector<CVUniverse*>>CCQENuUtilsNSF::GetErrorBands(ChainWrapper*chain){
  typedef std::map< std::string, std::vector<CVUniverse*> > SystMap;

  SystMap error_bands;

  // CV
  error_bands[std::string("CV")].push_back( new CVUniverse(chain,0) );

  //if RunCodeWithSystematics is not set....dont put these systematics...
  if(RunCodeWithSystematics){

  //Detector systematics, lateral shifts
  error_bands[std::string("Muon_Energy")].push_back(
      new MuonERangeCurvatureShiftUniverse(chain, -1));

  error_bands[std::string("Muon_Energy")].push_back(
      new MuonERangeCurvatureShiftUniverse(chain, +1));

//  error_bands[std::string("Pion_Response")].push_back(
  //    new RecoilShiftPionResponse(chain, -1));

 // error_bands[std::string("Pion_Response")].push_back(
   //   new RecoilShiftPionResponse(chain, +1));

 // error_bands[std::string("Proton_Response")].push_back(
   //   new RecoilShiftProtonResponse(chain, -1));

 // error_bands[std::string("Proton_Response")].push_back(
/*   //   new RecoilShiftProtonResponse(chain, +1));
      
  error_bands[std::string("Other_Response")].push_back(
      new RecoilShiftOtherResponse(chain, -1));

  error_bands[std::string("Other_Response")].push_back(
      new RecoilShiftOtherResponse(chain, +1));      
  
  error_bands[std::string("Muon_Energy_Resolution")].push_back(
      new MuonEnergyResolution(chain, -1));

  error_bands[std::string("Muon_Energy_Resolution")].push_back(
      new MuonEnergyResolution(chain, +1));
      
  */    
      
  
  //Flux
  int n_flux_universes = 50;
  SystMap flux_systematics = 
      PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain,n_flux_universes);
  error_bands.insert(flux_systematics.begin(), flux_systematics.end());

  //GENIE
  SystMap genie_systematics = 
      PlotUtils::GetGenieSystematicsMap<CVUniverse>(chain,true);// change that true to a switch on do_nonrespi_tune
  error_bands.insert(genie_systematics.begin(), genie_systematics.end());
  
  
  
  }

  return error_bands;

}
 
#endif
