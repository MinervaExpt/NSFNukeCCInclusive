lcgcmt_installarea=without_installarea
CMTPATH=/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/MINERVA/MINERVA_v22r1p1:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/LHCB/LHCB_v33r0p1b_lcgcmake:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/GEANT4/GEANT4_v94r2p2b_lcgcmake:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/MDBASE/MDBASE_v22r1p1:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/GAUDI/GAUDI_v22r4b_lcgcmake:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61
CMT_tag=$(tag)
CMTROOT=/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/lcgcmake/lcg_61a_forSL7/external/cmt/v1r20p20090520/x86_64-slc7-gcc49-opt/CMT/v1r20p20090520
CMT_root=/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/lcgcmake/lcg_61a_forSL7/external/cmt/v1r20p20090520/x86_64-slc7-gcc49-opt/CMT/v1r20p20090520
CMTVERSION=v1r20p20090520
cmt_hardware_query_command=uname -m
cmt_hardware=`$(cmt_hardware_query_command)`
cmt_system_version_query_command=${CMTROOT}/mgr/cmt_linux_version.sh | ${CMTROOT}/mgr/cmt_filter_version.sh
cmt_system_version=`$(cmt_system_version_query_command)`
cmt_compiler_version_query_command=${CMTROOT}/mgr/cmt_gcc_version.sh | ${CMTROOT}/mgr/cmt_filter3_version.sh
cmt_compiler_version=`$(cmt_compiler_version_query_command)`
PATH=/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew/InstallArea/${CMTCONFIG}/bin:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/MINERVA/MINERVA_v22r1p1/InstallArea/${CMTCONFIG}/bin:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/LHCB/LHCB_v33r0p1b_lcgcmake/InstallArea/${CMTCONFIG}/bin:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/GEANT4/GEANT4_v94r2p2b_lcgcmake/InstallArea/${CMTCONFIG}/bin:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/GAUDI/GAUDI_v22r4b_lcgcmake/InstallArea/${CMTCONFIG}/bin:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/lcgcmake/lcg_61a_forSL7/external/cmt/v1r20p20090520/x86_64-slc7-gcc49-opt/CMT/v1r20p20090520/${CMTBIN}:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/GAUDI/GAUDI_v22r4b_lcgcmake/InstallArea/x86_64-slc7-gcc49-opt/scripts:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1/postgres_client/x86_64-slc7-gcc49-opt/bin:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1/GENIE/x86_64-slc7-gcc49-opt/bin:/cvmfs/minerva.opensciencegrid.org/product/prd/gcc/v4_9_3/Linux64bit-3-10-2-17/bin:/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/cigetcert/v1_16_1/Linux64bit-3-10-2-17/bin:/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/jobsub_client/v1_2_10/NULL:/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/cpn/v1_7/NULL/bin:/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/awscli/v1_7_15/Linux64bit-2-6-2-12/bin:/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/ifdhc/v2_2_3/Linux64bit-3-10-2-17-python27/bin:/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/sam_web_client/v2_0/NULL/bin:/cvmfs/minerva.opensciencegrid.org/product/prd/sam/v8_8_5/Linux-2/bin:/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/ups/v6_0_7/Linux64bit-3-10-2-17/bin:/cvmfs/minerva.opensciencegrid.org/minerva/condor//condor_wrappers:/cvmfs/minerva.opensciencegrid.org/minerva/condor/:/cvmfs/minerva.opensciencegrid.org/minerva/valgrind-3.5.0/bin:/nashome/a/anezkak/.vscode-server/bin/c3f126316369cd610563c75b1b1725e0679adfb3/bin:/opt/puppetlabs/bin:/usr/local/sbin:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/MinervaScripts:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/cmake_sl7/bin:/cvmfs/minerva.opensciencegrid.org/minerva/ack:/usr/krb5/bin:/usr/local/bin:/usr/bin:/usr/sbin:/etc:/usr/etc:/bin:/sbin:.:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/MINERVA/MINERVA_v22r1p1/Tools/MinervaUtils/python:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/MINERVA/MINERVA_v22r1p1/Sim/MNVDetectormc/x86_64-slc7-gcc49-opt:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1/MINEXTERNAL/scripts
CLASSPATH=/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/lcgcmake/lcg_61a_forSL7/external/cmt/v1r20p20090520/x86_64-slc7-gcc49-opt/CMT/v1r20p20090520/java
debug_option=-g
cc=gcc
ccomp=$(cc) -c $(includes) $(cdebugflags) $(cflags) $(pp_cflags)
clink=$(cc) $(clinkflags) $(cdebugflags)
ppcmd=-I
preproc=c++ -MD -c 
cpp=c++
cppflags=-pipe -ansi -pedantic -W -Wall -Wwrite-strings -Wpointer-arith -Woverloaded-virtual 
pp_cppflags=-D_GNU_SOURCE
cppcomp=$(cpp) -c $(includes) $(cppdebugflags) $(cppflags) $(pp_cppflags)
cpplink=$(cpp) $(cpplinkflags) $(cppdebugflags)
for=g77
fflags=$(debug_option)
fcomp=$(for) -c $(fincludes) $(fflags) $(pp_fflags)
flink=$(for) $(flinkflags)
javacomp=javac -classpath $(src):$(CLASSPATH) 
javacopy=cp
jar=jar
X11_cflags=-I/usr/include
Xm_cflags=-I/usr/include
X_linkopts=-L/usr/X11R6/lib -lXm -lXt -lXext -lX11 -lm
lex=lex $(lexflags)
yaccflags= -l -d 
yacc=yacc $(yaccflags)
ar=ar cr
ranlib=ranlib
make_shlib=${CMTROOT}/mgr/cmt_make_shlib_common.sh extract
shlibsuffix=so
shlibbuilder=g++ $(cmt_installarea_linkopts) 
shlibflags=-shared
symlink=/bin/ln -fs 
symunlink=/bin/rm -f 
library_install_command=${symlink}
build_library_links=$(cmtexe) build library_links -tag=$(tags)
remove_library_links=$(cmtexe) remove library_links -tag=$(tags)
cmtexe=${CMTROOT}/${CMTBIN}/cmt.exe
build_prototype=$(cmtexe) build prototype
build_dependencies=$(cmtexe) -tag=$(tags) build dependencies
build_triggers=$(cmtexe) build triggers
format_dependencies=${CMTROOT}/mgr/cmt_format_deps.sh
implied_library_prefix=-l
SHELL=/bin/sh
q="
src=../src/
doc=../doc/
inc=../src/
mgr=../cmt/
application_suffix=.exe
library_prefix=lib
lock_command=lockfile 
unlock_command=rm -rf 
lock_name=cmt
lock_suffix=.lock
lock_file=${lock_name}${lock_suffix}
svn_checkout_command=/usr/bin/python ${CMTROOT}/mgr/cmt_svn_checkout.py 
MAKEFLAGS= --no-print-directory 
gmake_hosts=lx1 rsplus lxtest as7 dxplus ax7 hp2 aleph hp1 hpplus papou1-fe atlas
make_hosts=virgo-control1 rio0a vmpc38a
everywhere=hosts
install_command=cp 
uninstall_command=/bin/rm -f 
cmt_installarea_command=ln -fs 
cmt_uninstallarea_command=/bin/rm -f 
cmt_install_area_command=$(cmt_installarea_command)
cmt_uninstall_area_command=$(cmt_uninstallarea_command)
cmt_install_action=$(CMTROOT)/mgr/cmt_install_action.sh
cmt_installdir_action=$(CMTROOT)/mgr/cmt_installdir_action.sh
cmt_uninstall_action=$(CMTROOT)/mgr/cmt_uninstall_action.sh
cmt_uninstalldir_action=$(CMTROOT)/mgr/cmt_uninstalldir_action.sh
mkdir=mkdir
cmt_cvs_protocol_level=v1r1
cmt_installarea_prefix=InstallArea
CMT_PATH_remove_regexp=/[^/]*/
CMT_PATH_remove_share_regexp=/share/
NEWCMTCONFIG=x86_64-sl79-gcc493
NSFNukeCCInclusive_tag=$(tag)
NSFNUKECCINCLUSIVEROOT=/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive
NSFNukeCCInclusive_root=/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive
NSFNUKECCINCLUSIVEVERSION=${MINERVA_RELEASE}
NSFNukeCCInclusive_cmtpath=/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew
NSFNukeCCInclusive_offset=Ana
NSFNukeCCInclusive_project=MINERVA_v22r1p1_MADNew
NSFNukeCCInclusive_project_release=Minerva_v22r1p1_MADNew
NUKECCSRC_ANA=${NUKECCSRCROOT}/ana_common
MY_NSFNUKECC=/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive/ana
NUKECCSRC_INCLUDE_PATH= -I"../src"  
tag=x86_64-slc7-gcc49-opt
package=NSFNukeCCInclusive
version=${MINERVA_RELEASE}
PACKAGE_ROOT=$(NSFNUKECCINCLUSIVEROOT)
srcdir=../src
bin=../$(NSFNukeCCInclusive_tag)/
BIN=/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive/$(NSFNukeCCInclusive_tag)/
javabin=../classes/
mgrdir=cmt
project=MINERVA_v22r1p1_MADNew
cmt_installarea_paths= $(cmt_installarea_prefix)/$(CMTCONFIG)/bin $(GAUDI_installarea_prefix)/$(CMTCONFIG)/lib $(GAUDI_installarea_prefix)/share/lib $(GAUDI_installarea_prefix)/share/bin $(GEANT4_installarea_prefix)/$(CMTCONFIG)/lib $(GEANT4_installarea_prefix)/share/lib $(GEANT4_installarea_prefix)/share/bin $(LHCB_installarea_prefix)/$(CMTCONFIG)/lib $(LHCB_installarea_prefix)/share/lib $(LHCB_installarea_prefix)/share/bin $(MINERVA_installarea_prefix)/$(CMTCONFIG)/lib $(MINERVA_installarea_prefix)/share/lib $(MINERVA_installarea_prefix)/share/bin $(MINERVA_v22r1p1_MADNew_installarea_prefix)/$(CMTCONFIG)/lib $(MINERVA_v22r1p1_MADNew_installarea_prefix)/share/lib $(MINERVA_v22r1p1_MADNew_installarea_prefix)/share/bin
use_linkopts= $(cmt_installarea_linkopts)   $(NSFNukeCCInclusive_linkopts) 
LCGCMT_installarea_prefix=$(cmt_installarea_prefix)
LCGCMT_installarea_prefix_remove=$(LCGCMT_installarea_prefix)
LD_LIBRARY_PATH=/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew/InstallArea/${CMTCONFIG}/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/MINERVA/MINERVA_v22r1p1/InstallArea/${CMTCONFIG}/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/LHCB/LHCB_v33r0p1b_lcgcmake/InstallArea/${CMTCONFIG}/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/GEANT4/GEANT4_v94r2p2b_lcgcmake/InstallArea/${CMTCONFIG}/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/GAUDI/GAUDI_v22r4b_lcgcmake/InstallArea/${CMTCONFIG}/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../clhep/1.9.4.7/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../tcmalloc/2.6.2/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../libunwind/5c2cade/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../RELAX/RELAX_1_3_0p/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../ROOT/5.34.36/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../xrootd/3.3.6/x86_64-slc7-gcc49-opt/lib64:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../Boost/1.57.0_python2.7/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1/ppfx/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1/MinosGeom/x86_64-slc7-gcc49-opt:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1/NumiBeamDB/x86_64-slc7-gcc49-opt:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1/psycopg2/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1/postgres_client/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1/GDML/x86_64-slc7-gcc49-opt:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../sqlite/3070900/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../COOL/COOL_2_8_20/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../HepMC/2.06.08/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1/GENIE/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../MCGenerators_hepmc2.06.08/lhapdf/5.9.1/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../MCGenerators_hepmc2.06.08/pythia6/427.2/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/EXTERNAL/EXTERNAL_v22r1p1/LOG4CPP/x86_64-slc7-gcc49-opt:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../pytools/1.6_python2.7/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../mysql/5.5.14/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../HepPDT/2.06.01/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../POOL/POOL_2_9_19/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../CORAL/CORAL_2_3_27a/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../GSL/1.10/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../XercesC/3.1.1p1/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../uuid/1.38p1/x86_64-slc7-gcc49-opt/lib:/cvmfs/minerva.opensciencegrid.org/product/prd/gcc/v4_9_3/Linux64bit-3-10-2-17/lib64:/cvmfs/minerva.opensciencegrid.org/product/prd/gcc/v4_9_3/Linux64bit-3-10-2-17/lib:/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/ifdhc/v2_2_3/Linux64bit-3-10-2-17-python27/lib:/usr/lib64:/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew/Ana/UnfoldUtils/x86_64-slc7-gcc49-opt:/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/LCGCMT/LCGCMT_61/LCG_Settings/../../../Python/2.7.11/x86_64-slc7-gcc49-opt/lib
GAUDI_installarea_prefix=$(cmt_installarea_prefix)
GAUDI_installarea_prefix_remove=$(GAUDI_installarea_prefix)
MDBASE_installarea_prefix=$(cmt_installarea_prefix)
MDBASE_installarea_prefix_remove=$(MDBASE_installarea_prefix)
GEANT4_installarea_prefix=$(cmt_installarea_prefix)
GEANT4_installarea_prefix_remove=$(GEANT4_installarea_prefix)
EXTERNAL_installarea_prefix=$(cmt_installarea_prefix)
EXTERNAL_installarea_prefix_remove=$(EXTERNAL_installarea_prefix)
LHCB_installarea_prefix=$(cmt_installarea_prefix)
LHCB_installarea_prefix_remove=$(LHCB_installarea_prefix)
MINERVA_installarea_prefix=$(cmt_installarea_prefix)
MINERVA_installarea_prefix_remove=$(MINERVA_installarea_prefix)
MINERVA_v22r1p1_MADNew_installarea_prefix=$(cmt_installarea_prefix)
MINERVA_v22r1p1_MADNew_installarea_prefix_remove=$(MINERVA_v22r1p1_MADNew_installarea_prefix)
cmt_installarea_linkopts= -L/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew/$(MINERVA_v22r1p1_MADNew_installarea_prefix)/$(CMTCONFIG)/lib  -L/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/minerva/MINERVA/MINERVA_v22r1p1/$(MINERVA_installarea_prefix)/$(CMTCONFIG)/lib  -L/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/LHCB/LHCB_v33r0p1b_lcgcmake/$(LHCB_installarea_prefix)/$(CMTCONFIG)/lib  -L/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/GEANT4/GEANT4_v94r2p2b_lcgcmake/$(GEANT4_installarea_prefix)/$(CMTCONFIG)/lib  -L/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lhcb/GAUDI/GAUDI_v22r4b_lcgcmake/$(GAUDI_installarea_prefix)/$(CMTCONFIG)/lib 
CMTINSTALLAREA=/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew/$(cmt_installarea_prefix)
use_requirements=requirements $(CMT_root)/mgr/requirements 
use_includes= 
use_fincludes= $(use_includes)
use_stamps=  $(NSFNukeCCInclusive_stamps) 
use_cflags=  $(NSFNukeCCInclusive_cflags) 
use_pp_cflags=  $(NSFNukeCCInclusive_pp_cflags) 
use_cppflags=  $(NSFNukeCCInclusive_cppflags) 
use_pp_cppflags=  $(NSFNukeCCInclusive_pp_cppflags) 
use_fflags=  $(NSFNukeCCInclusive_fflags) 
use_pp_fflags=  $(NSFNukeCCInclusive_pp_fflags) 
includes= $(ppcmd)"$(srcdir)" $(use_includes)
fincludes= $(includes)
make_GUID={88BF15AB-5A2D-4bea-B64F-02752C2A1F4F}
constituents= 
all_constituents= $(constituents)
constituentsclean= 
all_constituentsclean= $(constituentsclean)
cmt_actions_constituents= make 
cmt_actions_constituentsclean= makeclean 
