# echo "Setting NUKECCSRC ${MINERVA_RELEASE} in /minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/lcgcmake/lcg_61a_forSL7/external/cmt/v1r20p20090520/x86_64-slc7-gcc49-opt/CMT/v1r20p20090520
endif
source ${CMTROOT}/mgr/setup.csh

set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=NUKECCSRC -version=${MINERVA_RELEASE} -path=/minerva/app/users/anezkak/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive  -no_cleanup $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

