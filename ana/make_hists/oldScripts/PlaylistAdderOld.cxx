#include <iostream>
#include <stdlib.h>
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"	
#include "PlotUtils/HistogramUtils.h"	
#include "TParameter.h" 
#include "../../NUKECCSRC/include/UtilsNSF.h"
#ifndef __CINT__
#endif

using namespace NUKECC_ANA;

double GetPOTScale( string basefilename, std::string playlist )
{
    double data_pot = 0.0;
    double mc_pot = 0.0;

    const std::string HISTS( getenv("HISTS") );
    const std::string NUKECC_TAG( getenv("NUKECC_TAG") );
    std::string histogram_path = HISTS + "/" + NUKECC_TAG;

    string filename = histogram_path + "/" + playlist + "/" + basefilename + ".root";
    TFile *histFile = new TFile( filename.c_str(), "read" );
    TVector2 *read_pot = (TVector2*)histFile->Get("POT_Used");
    data_pot = read_pot->X();
    mc_pot = read_pot->Y();
    delete histFile;

    const double dataMCScale = data_pot/mc_pot;
    if( dataMCScale > 1 ) cout << "I expected more POT in MC than data. Make sure you ran your jobs correctly or know this is true" << endl;
    cout << " The Data POT for " << playlist << " is " << data_pot << endl;
    cout << " The MC POT for " << playlist << " is " << mc_pot << endl;
    cout << " Scale " << playlist << "  " << dataMCScale << endl;

    return dataMCScale;
}

double GetPOTTotalScale( string basefilename, std::vector<string> playlists )
{
    double Data_sum = 0.0;
    double MC_sum = 0.0;

    const std::string HISTS( getenv("HISTS") );
    const std::string NUKECC_TAG( getenv("NUKECC_TAG") );
    std::string histogram_path = HISTS + "/" + NUKECC_TAG;

    for( unsigned i = 0; i < playlists.size(); i++ )
    {
         string filename = histogram_path + "/" + playlists[i] + "/" + basefilename + ".root";
         TFile *histFile = new TFile( filename.c_str(), "read" );
         TVector2 *read_pot = (TVector2*)histFile->Get("POT_Used");
         Data_sum += read_pot->X();
         MC_sum += read_pot->Y();
         delete histFile;
    }

    cout << " The Total Data POT for your samples is " << Data_sum << endl;
    cout << " The Total MC POT for your samples is " << MC_sum << endl;
    cout << " The Ratio is " << Data_sum/MC_sum << endl;
    
    return Data_sum/MC_sum;
}


int PlaylistAdder( bool antinu, string basefilename, std::vector<string> playlists )
{
    //Need to open files
    const int num_playlists = playlists.size(); //Full par input vector less 3 header values

    cout << "I'm gonna add your POT and get a scale factor!" << endl;
    double sum_scale = GetPOTTotalScale(basefilename, playlists);
    cout << " sum_scale    " << sum_scale << endl;

    const std::string HISTS( getenv("HISTS") );
    cout << "HISTS: " << HISTS << endl;

    const std::string NUKECC_TAG( getenv("NUKECC_TAG") );
    cout << "NUKECC_TAG: " << NUKECC_TAG << endl;

    TFile *f_input[num_playlists];
    std::string histogram_path = HISTS + "/" + NUKECC_TAG;
    cout << "The path for your histograms is " << histogram_path << endl;
    cout << "You are combining " << num_playlists << " playlists" << endl;

    //Read files to get histogram
    for( int i = 0; i < num_playlists; i++ )
    {
         cout << basefilename << "  " << playlists[i] << endl;
         string filename = histogram_path + "/" + playlists[i] + "/" + basefilename + ".root";
         cout << "Opening " << filename << endl;
         f_input[i] = new TFile(filename.c_str());
         
         if( f_input[i]->IsZombie() || f_input[i]->GetListOfKeys()->IsEmpty() )
         {
             Error("Playlistadder", "Could not get histogram ROOT file or it was empty.");
             return 1;
         }
    }
  
    string outname = "";
    if( antinu )
        outname = histogram_path + "/AllNubarME/" + basefilename + ".root";
    else   
        outname = histogram_path + "/AllNuME/" + basefilename + ".root";

    TFile *f_output = new TFile(outname.c_str(), "RECREATE");
   
    //Add hist
    cout << "file name " << f_input[0]->GetName() << endl;
    TList *file0List = f_input[0]->GetListOfKeys();
  
    TIter next(file0List);
    TObject *ob = 0;
 
    while( (ob = next()) )
    {
        bool dopotnorm = true;
        string name = ob->GetName();
        string classname = ((TKey*)ob)->GetClassName();

        cout << name << endl;
        if( name.find("Data") != std::string::npos || name.find("data") != std::string::npos || name.find("POT_Used") != std::string::npos ) 
        {
            cout << "The object name contains 'Data' or 'data' or 'POT_Used', I'm not going to scale it by POT " << endl;
            dopotnorm = false;
        }

        cout << "for " << name << " in " << basefilename << "  am I going to normalize by pot? " << dopotnorm << endl;
        
        if( classname == "TH1D" )
        {
            cout << "I'm summing a set of TH1D for " << name << endl;
            TH1D *hist1D = (TH1D*)f_input[0]->Get(name.c_str());
            //if( dopotnorm ) hist1D->Scale( GetPOTScale(basefilename, playlists[0]) );

            for( int i = 1; i < num_playlists; i++ )
            {
                 TH1D *tmp1D = (TH1D*)f_input[i]->Get(name.c_str());
                 //if( dopotnorm ) tmp1D->Scale( GetPOTScale(basefilename, playlists[i]) );
                 hist1D->Add(tmp1D);
                 delete tmp1D;
            }

            f_output->cd();
            hist1D->Write(name.c_str());
            delete hist1D;
        }

        else if( classname == "PlotUtils::MnvH1D" )
        {
            cout << "I'm summing a set of PlotUtils::MnvH1D for " << name << endl;
            PlotUtils::MnvH1D *hist1D = (PlotUtils::MnvH1D*)f_input[0]->Get(name.c_str());
            //if( dopotnorm ) hist1D->Scale( GetPOTScale(basefilename, playlists[0]) );

            for( int i = 1; i < num_playlists; i++ )
            {
                 PlotUtils::MnvH1D *tmp1D = (PlotUtils::MnvH1D*)f_input[i]->Get(name.c_str());
                 //if( dopotnorm ) tmp1D->Scale( GetPOTScale(basefilename, playlists[i]) );
                 hist1D->Add(tmp1D);
                 delete tmp1D;
            }

            f_output->cd();
            hist1D->Write(name.c_str());
            delete hist1D;
        }

        else if( classname == "TH2D" )
        {
            cout << "I'm summing a set of TH2D for " << name << endl;
            TH2D *hist2D = (TH2D*)f_input[0]->Get(name.c_str());
            //if( dopotnorm ) hist2D->Scale( GetPOTScale(basefilename, playlists[0]) );

            for( int i = 1; i < num_playlists; i++ )
            {
                 TH2D *tmp2D = (TH2D*)f_input[i]->Get(name.c_str());
                 //if( dopotnorm ) tmp2D->Scale( GetPOTScale(basefilename, playlists[i]) );
                 hist2D->Add(tmp2D);
                 delete tmp2D;
            }

            f_output->cd();
            hist2D->Write(name.c_str());
            delete hist2D;
        }

        else if( classname == "PlotUtils::MnvH2D" )
        {
            cout << "I'm summing a set of PlotUtils::MnvH2D for " << name << endl;
            PlotUtils::MnvH2D *hist2D = (PlotUtils::MnvH2D*)f_input[0]->Get(name.c_str());
            //if( dopotnorm ) hist2D->Scale( GetPOTScale(basefilename, playlists[0]) );

            for( int i = 1; i < num_playlists; i++ )
            {
                 PlotUtils::MnvH2D *tmp2D = (PlotUtils::MnvH2D*)f_input[i]->Get(name.c_str());
                 //if( dopotnorm ) tmp2D->Scale( GetPOTScale(basefilename, playlists[i]) );
                 hist2D->Add(tmp2D);
                 delete tmp2D;
            }

            f_output->cd();
            hist2D->Write(name.c_str());
            delete hist2D;
        }

        else if( name == "POT_Used" )
        {
            cout << "I'm summing a set of TVector2 for " << name << endl;

            TVector2 *read_pot = (TVector2*)f_input[0]->Get("POT_Used");
            double Data_sum = read_pot->X();
            double MC_sum = read_pot->Y();
            delete read_pot;

            for( int i = 1; i < num_playlists; i++ )
            {
                 TVector2 *tmp_read_pot = (TVector2*)f_input[i]->Get("POT_Used");
                 Data_sum += tmp_read_pot->X();
                 MC_sum += tmp_read_pot->Y();
                 delete tmp_read_pot;
            }

            TVector2 *pot = new TVector2( Data_sum, MC_sum );
            f_output->cd();
            pot->Write( name.c_str() );
            delete pot;
        }

        else
        {
            cout << "This is a type of class I don't know how to handle: " << classname << "\tskipping. Get your life right. " << endl;
        }

    }//obj next

    ob = NULL;
    file0List = NULL;
    f_output->Close();
 
    return 0;
};

int main( int argc, char *argv[] )
{
    TH1::AddDirectory(false);
    if( argc == 1 )
    {
        std::cout << "-----------------------------------------------------------------------------------------------" << std::endl;
        std::cout << "MACROS HELP:\n\n" <<
        "PlaylistAdder outfile_root  file_genie_variation_xsection\n\n"<<
        "\t-Base string of files. Assumes files are of the form XXXXXX_minervaY.root where Y is any playlist(s).\n" <<
        "\t-Type of file we are adding. Valid choices: MuonEventSelection, EffPurity, BackgroundTemples, Migration\n" <<
        "\t-MCOnly sample. Valid choices: 0 = no special mconly sample (2p2h for instance), 1 = yes\n" <<
        "\t-List of playlists. Example would be minerva1 minerva7 minerva 9. This will add 1,7,9 together.\n" << std::endl;
        std::cout << "---------------------------------------------------------------------------------" << std::endl;
        return 0;
    
    }

    //! Default parameters
    std::vector<std::string> par;
    std::vector<std::string> playlists;
    par.push_back("PlaylistAdder");
    par.push_back(argv[1]); //0 nu, 1 nubar 
    par.push_back(argv[2]);

    //! Set user parameters
    for( int i = 3; i < argc; ++i )
    {
         par.push_back(argv[i]);
         playlists.push_back(argv[i]);
    }

    for( unsigned int i = 0; i < par.size(); ++i )
         std::cout << "Parameter " << i << ": " << par[i] << std::endl;
 
    return PlaylistAdder( bool(atoi(argv[1])), par[2], playlists );

}
