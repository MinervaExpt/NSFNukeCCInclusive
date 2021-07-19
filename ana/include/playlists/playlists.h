

namespace NUKECC_ANA{
  
  bool useMADtuples=false;
  bool useNukeOnly=false;

  string get_mc_files( string playlist, int targetID){
    //cout<<"Playlist: "<<playlist<<endl;
    if (!(playlist=="minervame1A" || playlist=="minervame1B" || playlist=="minervame1C" ||playlist=="minervame1D" || playlist=="minervame1E" || playlist=="minervame1F" || playlist=="minervame1G" || playlist=="minervame1L" || playlist=="minervame1M" || playlist=="minervame1N" || playlist=="minervame1O" || playlist=="minervame1P" || playlist=="minervame5A" || playlist=="minervame6A" || playlist=="minervame6B") ){ 
cout << "The specified playlist is not in the list for playlists.h" << endl;
      return "ERROR: playlist does not exist";
    }
    string dir = getenv("MY_NSFNUKECC");
    string base= useMADtuples? "MasterAnaDev":"NukeCC";
    string pname= dir+"/include/playlists/"+base+"_";
    (useNukeOnly && targetID<10)? pname+="MCNuke":pname+="MC";
    pname+="_"+playlist+"_MuonKludged.txt";
    cout<<"Playlist: " <<pname<<endl;
    return pname;

  }

 string get_data_files( string playlist){
    //cout<<"Playlist: "<<playlist<<endl;
    if (!(playlist=="minervame1A" || playlist=="minervame1B" || playlist=="minervame1C" ||playlist=="minervame1D" || playlist=="minervame1E" || playlist=="minervame1F" || playlist=="minervame1G" || playlist=="minervame1L" || playlist=="minervame1M" || playlist=="minervame1N" || playlist=="minervame1O" || playlist=="minervame1P" || playlist=="minervame5A" || playlist=="minervame6A" || playlist=="minervame6B") ){ 
cout << "The specified playlist is not in the list for playlists.h" << endl;
      return "ERROR: playlist does not exist";
    }
    string dir = getenv("MY_NSFNUKECC");
    string base= useMADtuples? "MasterAnaDev":"NukeCC";
    string pname= dir+"/include/playlists/"+base+"_Data_"+playlist+"_MuonKludged.txt";
    cout<<"Playlist: " <<pname<<endl;
    return pname;

  }
}
