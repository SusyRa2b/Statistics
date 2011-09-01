#!/bin/tcsh

#
# Args: mg mlsp Xsec iter
#

  set jobdir = `pwd`

  echo jobdir is $jobdir

  cd /afs/cern.ch/user/o/owen/rel-dirs/CMSSW_4_2_5

  echo doing cmsenv

  cmsenv

  which root

  echo setting up root 5.28

  echo first
# source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.csh
##----------------------------------------
set gcc_config_version = 4.3.2
set mpfr_config_version = 2.3.1
set gmp_config_version=4.2.2
set LCGPLAT = x86_64-slc5-gcc34-opt
set LCG_lib_name = lib64

set LCG_contdir = /afs/cern.ch/sw/lcg/contrib
set LCG_gcc_home = ${LCG_contdir}/gcc/${gcc_config_version}/${LCGPLAT}
set LCG_mpfr_home = ${LCG_contdir}/mpfr/${mpfr_config_version}/${LCGPLAT}
set LCG_gmp_home=${LCG_contdir}/gmp/${gmp_config_version}/${LCGPLAT}

setenv PATH ${LCG_gcc_home}/bin:${PATH}

setenv LD_LIBRARY_PATH ${LCG_gcc_home}/${LCG_lib_name}:${LCG_mpfr_home}/lib:${LCG_gmp_home}/lib:${LD_LIBRARY_PATH}
#setenv LD_LIBRARY_PATH ${LCG_gcc_home}/${LCG_lib_name}:${LCG_mpfr_home}/lib:${LCG_gmp_home}/lib
##----------------------------------------

  echo second
#  source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.28.00c/x86_64-slc5-gcc43-opt/root/bin/thisroot.csh
##----------------------------------------
if ($?ROOTSYS) then
   setenv OLD_ROOTSYS "$ROOTSYS"
endif

setenv ROOTSYS /afs/cern.ch/sw/lcg/app/releases/ROOT/5.28.00c/x86_64-slc5-gcc43-opt/root

if ($?OLD_ROOTSYS) then
   if ( ! -e $ROOTSYS/bin/drop_from_path ) then
      echo "ERROR: the utility drop_from_path has not been build yet. Do:"
      echo "make bin/drop_from_path"
      exit 1
   endif   
   setenv PATH `$ROOTSYS/bin/drop_from_path -e "$OLD_ROOTSYS/bin"`
   if ($?LD_LIBRARY_PATH) then
      setenv LD_LIBRARY_PATH `$ROOTSYS/bin/drop_from_path -D -e -p "$LD_LIBRARY_PATH" "$OLD_ROOTSYS/lib"`
   endif
   if ($?DYLD_LIBRARY_PATH) then
      setenv DYLD_LIBRARY_PATH `$ROOTSYS/bin/drop_from_path -D -e -p "$DYLD_LIBRARY_PATH" "$OLD_ROOTSYS/lib"`
   endif
   if ($?SHLIB_PATH) then
      setenv SHLIB_PATH `$ROOTSYS/bin/drop_from_path -D -e -p "$SHLIB_PATH" "$OLD_ROOTSYS/lib"`
   endif
   if ($?LIBPATH) then
      setenv LIBPATH `$ROOTSYS/bin/drop_from_path -D -e -p "$LIBPATH" "$OLD_ROOTSYS/lib"`
   endif
   if ($?PYTHONPATH) then
      setenv PYTHONPATH `$ROOTSYS/bin/drop_from_path -D -e -p "$PYTHONPATH" "$OLD_ROOTSYS/lib"`
   endif
   if ($?MANPATH) then
      setenv MANPATH `$ROOTSYS/bin/drop_from_path -D -e -p "$MANPATH" "$OLD_ROOTSYS/man"`
   endif
endif


if ($?MANPATH) then
# Nothing to do
else
   # Grab the default man path before setting the path to avoid duplicates 
   if ( -X manpath ) then
      set default_manpath = `manpath`
   else
      set default_manpath = `man -w`
   endif
endif

set path = ($ROOTSYS/bin $path)

if ($?LD_LIBRARY_PATH) then
   setenv LD_LIBRARY_PATH $ROOTSYS/lib:$LD_LIBRARY_PATH      # Linux, ELF HP-UX
else
   setenv LD_LIBRARY_PATH $ROOTSYS/lib
endif

if ($?DYLD_LIBRARY_PATH) then
   setenv DYLD_LIBRARY_PATH $ROOTSYS/lib:$DYLD_LIBRARY_PATH  # Mac OS X
else
   setenv DYLD_LIBRARY_PATH $ROOTSYS/lib
endif

if ($?SHLIB_PATH) then
   setenv SHLIB_PATH $ROOTSYS/lib:$SHLIB_PATH                # legacy HP-UX
else
   setenv SHLIB_PATH $ROOTSYS/lib
endif

if ($?LIBPATH) then
   setenv LIBPATH $ROOTSYS/lib:$LIBPATH                      # AIX
else
   setenv LIBPATH $ROOTSYS/lib
endif

if ($?PYTHONPATH) then
   setenv PYTHONPATH $ROOTSYS/lib:$PYTHONPATH
else
   setenv PYTHONPATH $ROOTSYS/lib
endif

if ($?MANPATH) then
   setenv MANPATH `dirname $ROOTSYS/man/man1`:$MANPATH
else
   setenv MANPATH `dirname $ROOTSYS/man/man1`:$default_manpath
endif
##----------------------------------------

  echo going to working dir

  cd /afs/cern.ch/user/o/owen/ra2b-likelihood/UserCode/SusyAnalysis/RA2b/Statistics

  which root

  echo args $1 $2


  echo ensuring output directories exist

  mkdir -p log-files-ge2btight-t1bbbb-iter$4
# mkdir -p output-files-ge2btight-t1bbbb-iter$4



#---------
  root -b -q cls-scripts-v2/runcls_ln_ge2btight_t1bbbb.c\($1,$2,$3\) >& $jobdir/cls-ge2btight-t1bbbb-$1-$2-iter$4.log

  echo Done with running root.

  echo listing of jobdir $jobdir is
  ls -l $jobdir


  echo Log file is cls-ge2btight-t1bbbb-$1-$2-iter$4.log   Checking it...
  ls -l $jobdir/cls-ge2btight-t1bbbb-$1-$2-iter$4.log

  echo doing the greps
  grep "test stat"    $jobdir/cls-ge2btight-t1bbbb-$1-$2-iter$4.log >  log-files-ge2btight-t1bbbb-iter$4/short-cls-ge2btight-t1bbbb-$1-$2-iter$4.log
  grep "p-value of"   $jobdir/cls-ge2btight-t1bbbb-$1-$2-iter$4.log >> log-files-ge2btight-t1bbbb-iter$4/short-cls-ge2btight-t1bbbb-$1-$2-iter$4.log
  grep "final result" $jobdir/cls-ge2btight-t1bbbb-$1-$2-iter$4.log >> log-files-ge2btight-t1bbbb-iter$4/short-cls-ge2btight-t1bbbb-$1-$2-iter$4.log

  gzip -9v log-files-ge2btight-t1bbbb-iter$4/short-cls-ge2btight-t1bbbb-$1-$2-iter$4.log
#---------
# root -b -q cls-scripts-v2/runcls_ln_ge2btight_t1bbbb.c\($1,$2,$3\) >& log-files-ge2btight-t1bbbb-iter$4/cls-ge2btight-t1bbbb-$1-$2-$4.log

# grep "test stat"    log-files-ge2btight-t1bbbb-iter$4/cls-ge2btight-t1bbbb-$1-$2-$4.log >  log-files-ge2btight-t1bbbb-iter$4/short-cls-ge2btight-t1bbbb-$1-$2-$4.log
# grep "p-value of"   log-files-ge2btight-t1bbbb-iter$4/cls-ge2btight-t1bbbb-$1-$2-$4.log >> log-files-ge2btight-t1bbbb-iter$4/short-cls-ge2btight-t1bbbb-$1-$2-$4.log
# grep "final result" log-files-ge2btight-t1bbbb-iter$4/cls-ge2btight-t1bbbb-$1-$2-$4.log >> log-files-ge2btight-t1bbbb-iter$4/short-cls-ge2btight-t1bbbb-$1-$2-$4.log

# gzip -9v log-files-ge2btight-t1bbbb-iter$4/cls-ge2btight-t1bbbb-$1-$2-$4.log
#---------



  echo all done.
