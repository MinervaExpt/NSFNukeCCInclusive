import os,sys,time
from optparse import OptionParser
import datetime
# made to sit just one up of the NSFNukeCCInclusive

###########################################################
#  Script: JobSubmission python script for submitting
#          analysis jobs to the grid
###########################################################

# Write the command you used to run your analysis
def writeEventLoop(mywrapper,outdir):
    #mywrapper.write("ldd runEventLoopALL_trueCCwo2p2h_tracker\n")
    #mywrapper.write("./runEventLoopALL_trueCCwo2p2h " + str(output_dir) + " 3 " + "6" + "\n")
    #mywrapper.write("./runEventLoopALL_trueCCwo2p2h_tracker " + str(output_dir) + " minervame1A\n")
    #mywrapper.write("./runEventLoop_tracker " + str(output_dir) + "\n")
    mywrapper.write("source runOver6GHIJPlaylists.sh\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)


def writeSetups(mywrapper,basedir):
    topdir = "$CONDOR_DIR_INPUT" + basedir
    mywrapper.write("cd "+topdir+"\n")
    mywrapper.write("cd "+topdir+"\n")
    mywrapper.write("source setup_utils.sh\n")
    mywrapper.write("cd "+topdir+"\n")
    mywrapper.write("source /cvmfs/minerva.opensciencegrid.org/minerva/setup/setup_minerva_products.sh\n")
    mywrapper.write("cd "+topdir+"\n")
    # if topdir MAT_Github
    mywrapper.write("cd NSFNukeCCInclusive\n")
    mywrapper.write("export NUKECCSRCROOT=${CONDOR_DIR_INPUT}/minerva/app/users/$USER/MAT_GitHub/NSFNukeCCInclusive/NUKECCSRC\n")
    mywrapper.write("export NUKECCSRC_ANA=${NUKECCSRCROOT}\n")
    mywrapper.write("export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CONDOR_DIR_INPUT}/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/NUKECCSRC/src\n")
    mywrapper.write("echo ${LD_LIBRARY_PATH}\n")

    mywrapper.write("cd "+topdir+"\n")
    mywrapper.write("cd NSFNukeCCInclusive/ana/make_hists\n")
    mywrapper.write("cd fluxConstraint\n")
    

# Define function to create tarball
def createTarball(outdir,tag,basedir):
    print("I'm inside createTarball()")
    found = os.path.isfile("%s/myareatar_%s.tar.gz"%(outdir,tag))
    if(not found):
        cmd = "tar -czf /minerva/app/users/$USER/myareatar_%s.tar.gz %s"%(tag,basedir)
        print "Making tar",cmd
        os.system(cmd)
        cmd2 = "cp /minerva/app/users/$USER/myareatar_%s.tar.gz %s/"%(tag,outdir)
        print "Copying tar",cmd2
        os.system(cmd2)


# Define function to unpack tarball
def writeTarballProceedure(mywrapper,tag,basedir):
    print "I will be making the tarball upacking with this version" 
    print "Path is",basedir
    mywrapper.write("cd $CONDOR_DIR_INPUT\n")
    mywrapper.write("tar -xvzf myareatar_%s.tar.gz\n"%tag)
    mywrapper.write("cd $CONDOR_DIR_INPUT\n")
    mywrapper.write("export MINERVA_PREFIX=${CONDOR_DIR_INPUT}/minerva/app/users/anezkak/MAT_GitHub/opt\n")
    writeSetups(mywrapper,basedir)

def writeOptions(parser):
    print "Now write options to the parser" 
    # Directory to write output
    parser.add_option('--outdir', dest='outdir', help='Directory to write output to', default = "/pnfs/minerva/persistent/users/"+_user_+"/default_analysis_loc/")
    parser.add_option('--basedir', dest='basedir', help='Base directory for making tarball', default = "NONE")
    parser.add_option('--stage', dest='stage', help='Process type', default="NONE")
    parser.add_option('--sample', dest='sample', help='Sample type', default="NONE")
    parser.add_option('--playlist', dest='playlist', help='Playlist type', default="NONE")
    ##########################################################################
    #  Options for making tar files....Basically you can make tarfiles 
    #######################################################################
    parser.add_option('--tag', dest='tag', help="Tag your release",default="tag_")
    parser.add_option('--mail',dest='mail',help="Want mail after job completion or failure",default=False,action="store_true")
    parser.add_option('--sametar',dest='sametar',help="Recycle the same tar file for jobsubmission",default=False,action="store_true")
    parser.add_option('--tarfilename',dest='tarfilename',help='Name of the tarfile you want to use',default="NA")
    parser.add_option('--notimestamp',dest='notimestamp',help='Flag to TURN OFF time stamp in the tag',default=False,action="store_true")

# Valid stages for the neutrinos (you can add more for your own analysis)
valid_stages=["eventLoop"]
avail_playlists=["minervame1A","minervame1D", "minervame1E", "minervame1F", "minervame1G","minervame1L","minervame1M","minervame1N", "minervame5A", "minervame6A","minervame6B", "minervame6D", "minervame6F", "minervame6G"]
# Get user's name
_user_ = os.getenv("USER")
usage = "usage: %prog[opts]" 
parser = OptionParser(usage=usage)
writeOptions(parser)
(opts,args) = parser.parse_args()

# Now print some information for user's selected options
if opts.mail:
    print "############################################" 
    print "You have opted to Get mails once the jobs are finished of failed." 
    print "#############################################" 

##############################################
#  TODO: Here I want to put the tag scheme...
##############################################
tag_name = str(opts.tag)
time_stamp = int(time.time())
if tag_name=="tag_":
    print "YOU DIDNT SPECIFY ANY ANY TAG...SO I WILL USE MY OWN TAGGING SCHEME" 
    # You can add in your own tagging scheme
    tag_name += "default_" 

print "Is everything fine till here*********" 
# Add the time stamp to tag
if not opts.notimestamp:
    tag_name += str(time_stamp)

print "Tag for this version is ",tag_name
print "******************************************************" 

# This is the output directory after the job is finished
output_dir = "$CONDOR_DIR_HISTS/" 

# Make outdir if not exist
if(not os.path.isdir(opts.outdir)):
    os.makedirs(opts.outdir)

memory = 2000

# Check the stage is valid or not
if opts.stage not in valid_stages:
    print opts.stage,"Selected stage is not valid. Here are the valid options",valid_stages
    sys.exit()

##############################################
#  Create wrapper
##############################################

wrapper_name = "%s_%s_wrapper_%s.sh"%(opts.stage,opts.playlist,tag_name)

mywrapper = open(wrapper_name,"w")
mywrapper.write("#!/bin/sh\n")
# Now create tarball
if (opts.sametar==False):
    createTarball(opts.outdir,tag_name,opts.basedir)
else:
    tarname = str(opts.tarfilename)
    if not os.path.exists(opts.outdir+opts.tarfilename):
        print "Tar File "+opts.outdir+opts.tarfilename+" doesn't Exist!" 
        sys.exit()
    #change the tag to the current one...
    cmd="cp "+opts.outdir+opts.tarfilename+" "+opts.outdir+"myareatar_"+tag_name+".tar.gz" 
    os.system(cmd)

# This will unpack the tarball we just made above
writeTarballProceedure(mywrapper,tag_name,opts.basedir)

# Now the add the command to run event loop 
if(opts.stage=="eventLoop"):
    writeEventLoop(mywrapper,opts.outdir)

mywrapper.close()

# Making the wrapper we just created readable, writable and executable by everyone
os.system("chmod 777 %s"%(wrapper_name))

# TODO: not sure if this is needed
configstring = "" 
# Get gcc release
#gccstring = "x86_64-slc7-gcc49-opt" 
# Now add the execute command
cmd = "" 
cmd += "jobsub_submit --group minerva " #Group of experiment
#cmd += "--cmtconfig "+gccstring+" " #Setup minerva soft release built with minerva configuration
#cmd += "--OS sl7 " #Operating system #Not needed in SL7
cmd += " --lines '+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\\\"' "


if opts.mail:
    cmd += "-M " #this option to make decide if you want the mail or not
#cmd += "--subgroup=Nightly " #This is only for high priority jobs
cmd += "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC " 
cmd += "--role=Analysis " 
cmd += "--expected-lifetime 36h " 
cmd += "--memory "+str(memory)+"MB " 
#cmd += configstring+" " #the environments for the tunes to bee applied
cmd += "-f "+opts.outdir+"/myareatar_"+tag_name+".tar.gz " 
#cmd += "-i /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1"+" " 
cmd += "file://"+os.environ["PWD"]+"/"+wrapper_name
print cmd

os.system(cmd)

print "Deleting the app area tar..... " 
os.system("rm -rf /minerva/app/users/$USER/myareatar_"+tag_name+".tar.gz")
print "Sleeping" 

time.sleep(1)
