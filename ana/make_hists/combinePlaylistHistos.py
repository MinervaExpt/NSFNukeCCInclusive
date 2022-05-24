# File:         combinePlaylistHistos.py
# Brief:        Program to combine per-playlist (=MINERvA run period) root files filled with histograms.
#               This is applicable for eventSelection or background subtracted event selection, efficiency,
#               migration and plastic sideband file. 
#               Note on systematics: If the error bands are not consistent across the files to be merged, 
#               the macro will result in an error.
# Procedure:    Data histograms are simply added. MC histograms are scaled using the following: 
#               h_mc_me6A*data_pot_me6A/mc_pot_me6A + h_mc_me6B*data_pot_me6B/mc_pot_me6B+....
#               Total MC POT and total Data POT is stored in the final root file.
# Usage:        python combinePlaylistHistos.py /path/to/dir/with/files/to/merge outfile.root
# Author:       Anezka Klustova a.klustova20@imperial.ac.uk

import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import TParameter
import os,sys,math
import time

ROOT.TH1.AddDirectory(False) # to get rid of weird python errors

path = str(sys.argv[1])
allfiles = [f for f in os.listdir(path) if f.endswith('.root')]
#os.listdir(path)
# make sure that this path only has the files that needs to be merged

# Names of the combined files
print("Files in the directory to be merged:")
print("-------------------------------------------------------------")
for filename in allfiles:
    if filename.endswith(".root"): 
        print(str(filename))
print("-------------------------------------------------------------")

_file0 = ROOT.TFile(path+"/"+allfiles[0],"r")
print("Adding: " + str(_file0.GetName()))

if _file0.IsZombie() or _file0.GetListOfKeys().IsEmpty():
    raise TypeError("Could not get histogram ROOT file or it was empty.")

# initialize variables to store total pot
total_mc_pot = 0
total_data_pot = 0

# first file pot
mcPOT = _file0.Get("MCPOT").GetVal()
dataPOT = _file0.Get("DataPOT").GetVal()
print mcPOT
print dataPOT

# add to total pot
total_mc_pot += mcPOT
total_data_pot += dataPOT

# define outfile
_outfile = ROOT.TFile(path+"/"+str(sys.argv[2]),"RECREATE") # write to the directory of files that were merged
_outfile.cd()

outkeys={}
add_type = ["PlotUtils::MnvH1D", "PlotUtils::MnvH2D", "TParameter<double>" ]

# first file
for key in _file0.GetListOfKeys():
    _name = key.GetName()
    temp_key = (_file0.Get(_name)).Clone()
    if temp_key.ClassName() in add_type and 'mc' in _name:
        #print temp_key.ClassName()
        #print _name
        temp_key.Scale(dataPOT/mcPOT) # scale using MC scale
        outkeys[_name] = temp_key
    elif temp_key.ClassName() in add_type and 'data' in _name:
    # data files are left unchanged
        outkeys[_name] = temp_key
    elif temp_key.ClassName() in add_type:
        outkeys[_name] = temp_key
    else: 
        raise TypeError("Type not defined in the add list.")

# rest of the files
for i in range(1,len(allfiles)):
    _file = ROOT.TFile(path+"/"+allfiles[i],"r")
    print("Adding: "+ str(allfiles[i]))
    if _file.IsZombie() or _file.GetListOfKeys().IsEmpty():
        raise TypeError("Could not get histogram ROOT file or it was empty.")

    # adding to total mc and data pot
    mcPOT = _file.Get("MCPOT").GetVal()
    dataPOT = _file.Get("DataPOT").GetVal()
    print mcPOT
    print dataPOT

    total_mc_pot += mcPOT
    total_data_pot += dataPOT
    
    for key in _file.GetListOfKeys():
        _name = key.GetName()
        temp_key = _file.Get(_name)
        # adding data files (just sum)
        if temp_key.ClassName() in add_type and 'data' in _name: # data histos are left unchanged just added
            #print temp_key.ClassName()
            #print _name
            key_mem = outkeys.get(_name)
            key_mem.Add(temp_key)
            _outfile.WriteTObject(key_mem,_name)

        if temp_key.ClassName() in add_type and 'mc' in _name: 
            # mc histos are combined using 
            # h_mc_me6A*data_pot_me6A/mc_pot_me6A + h_mc_me6B*data_pot_me6B/mc_pot_me6B procedure
            #print temp_key.ClassName()
            #print _name
            temp_key.Scale(dataPOT/mcPOT)
            key_mem = outkeys.get(_name)
            key_mem.Add(temp_key)
            _outfile.WriteTObject(key_mem,_name)

# print total POT
print("-------------------------------------------------------------")
print("Total MC POT: " + str(total_mc_pot))
print("Total Data POT: " + str(total_data_pot))

# write total pot to the out file
_outfile.cd()
mcPOTout = TParameter(float)("MCPOT", total_mc_pot)
mcPOTout.Write()
dataPOTout = TParameter(float)("DataPOT", total_data_pot)
dataPOTout.Write()
# close file
_outfile.Close()
