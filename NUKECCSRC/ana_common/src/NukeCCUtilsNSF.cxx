#ifndef MNV_NUKECCUTILSNSF_cxx
#define MNV_NUKECCUTILSNSF_cxx 1

//#include "NukeCCUtilsNSF.h"
#include "include/NukeCCUtilsNSF.h"

//#include "../include/NukeCC_Cuts.h"
//#include "CCQENuUtilsNSF.h"
#include "include/CVUniverse.h"
//#include "include/GlobalIncludes.h" 
#include "include/LateralSystematics.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/MinosEfficiencySystematics.h"


#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "MinervaUnfold/MnvResponse.h"
#include "MinervaUnfold/MnvUnfold.h"

using namespace NUKECC_ANA; 

NukeCCUtilsNSF::NukeCCUtilsNSF(string playlist){
    cvhistos1D.clear();
    cvhistos2D.clear();
     cutter = new NukeCC_Cuts();
//   GlobalParameters::Get().m_useFluxConstraint = false;
   DefaultCVUniverse::SetPlaylist(playlist);
  // GlobalParameters::Get().m_usePPFX1Flux = false;
   use_merged_files = false;
   
   if(neutrinoMode){
   incoming_pdg = 14;
   cout<<"Initializing NukeCCUtils in Neutrino Mode"<<endl;
   }
   else{
   incoming_pdg = -14;
   cout<<"Initializing NukeCCUtils in Antineutrino Mode"<<endl;
   
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
 


NukeCCUtilsNSF::NukeCCUtilsNSF()
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
   cout<<"Initializing NukeCCUtils in Neutrino Mode"<<endl;
   }
   else{
   incoming_pdg = -14;
   cout<<"Initializing NukeCCUtils in Antineutrino Mode"<<endl;
   
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
 NukeCCUtilsNSF::~NukeCCUtilsNSF()
 {
    delete cutter;
     cvhistos1D.clear();
     cvhistos2D.clear();
 }



std::vector<std::string> NukeCCUtilsNSF::GetStdPlaylists( HelicityType::t_HelicityType helicity ) const
{
    std::vector<std::string> playlists;
    if( HelicityType::kAntiNeutrino != helicity )
    {
        playlists.push_back( "minerva1" );
        playlists.push_back( "minerva7" );
        playlists.push_back( "minerva9" );
        playlists.push_back( "minervame6A" );
        playlists.push_back( "minervame6B" );
        playlists.push_back( "minervame6H" );
        playlists.push_back( "minervame6I" );
        playlists.push_back( "minervame6J" );
        playlists.push_back( "minerva13C" );
                                        }
        //                                    
          if( HelicityType::kNeutrino != helicity )
          {
        playlists.push_back( "minervame1A" );

}
    
    return playlists;
}

HelicityType::t_HelicityType NukeCCUtilsNSF::GetHelicityFromPlaylist( std::string playlist ) const
{
 //   if ("" == playlist)
 //       playlist = fPlaylist;
    
        if(Playlist::minervame1A == playlist ||
           Playlist::minervame1B == playlist || 
           Playlist::minervame1L == playlist ) {
        return HelicityType::kNeutrino;
    }
       if (Playlist::minervame5A == playlist ||
        Playlist::minervame6A == playlist    ||
        Playlist::minervame6B == playlist    ||
        Playlist::minervame6H == playlist    
        ) {
        return HelicityType::kAntiNeutrino;
    }
    
    const std::vector<std::string> nus = GetStdPlaylists( HelicityType::kNeutrino );
    if ( find( nus.begin(), nus.end(), playlist ) != nus.end() )
        return HelicityType::kNeutrino;
    
    const std::vector<std::string> antinus = GetStdPlaylists( HelicityType::kAntiNeutrino );
    if ( find( antinus.begin(), antinus.end(), playlist ) != antinus.end() )
        return HelicityType::kAntiNeutrino;
    
    Warning( "NukeCCUtilsNSF::GetHelicityFromPlaylist", Form( "Playlist '%s' is not a standard nu or antinu playlist. Using any neutrino.", playlist.c_str() ) );
    return HelicityType::kAny;
}

TString NukeCCUtilsNSF::GetHistFileName( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ, HelicityType::t_HelicityType helicity) const
{
    TString histFileName;
    
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

// NEW for grid!!
TString NukeCCUtilsNSF::GetHistFileNameGrid( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ) const
{
    TString histFileName;
    
    TString variation = GetVariationTag();
    if( FileType::kAny == fType )
    {
        histFileName += Form("/Hists_%s_t%d_z%02d.root",
                             histType.c_str(),
                             targetID, targetZ);
    }
    else
    {
        histFileName += Form("/Hists_%s_%s_t%d_z%02d.root",
                             histType.c_str(),
                             GetFileTypeString(fType).c_str(),
                             targetID, targetZ);
                             
    }
    return histFileName; 
}


TString NukeCCUtilsNSF::HistDir( bool forceDisk /* = false */, bool ignorePlaylist /*= false*/, bool needsbluearc /*= false*/ ) const
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
            Error( "NukeCCUtilsNSF::HistDir", "Environmental variable HISTS not defined." );
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
            Error( "NukeCCUtilsNSF::HistDir", "Environmental variable NUKECC_TAG not defined." );
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
                Error( "NukeCCUtilsNSF::HistDir", Form("Could not make hist directory '%s'", histdir.Data() ) );
        }
        
        return histdir;
    }
}

TString NukeCCUtilsNSF::GetVariationTag() const
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

std::string NukeCCUtilsNSF::GetHelicityString( HelicityType::t_HelicityType helicity ) const
{
    if(      HelicityType::kNeutrino     == helicity ) return "Nu";
    else if( HelicityType::kAntiNeutrino == helicity ) return "AntiNu";
    else return "AnyHelicity";
}
TString NukeCCUtilsNSF::GetHistName( const std::string& histType, FileType::t_FileType fType, const std::string& var, int targetID, int targetZ ) const
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

std::string NukeCCUtilsNSF::GetFileTypeString( FileType::t_FileType fType ) const
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

 
  void writePOT(TFile *f , double DataPOT, double MCPOT){
    //PlotUtils::ChainWrapper* chainData = util.m_data;
    //PlotUtils::ChainWrapper* chainMC = util.m_mc;
  double Data = 0;
  double MC   = 0;
  //this->getPOT(data,mc);
  TVector2 *pot = new TVector2( Data, MC );
   f->WriteTObject( pot, "pot" );
   }

 
#endif
