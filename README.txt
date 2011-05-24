
   This directory is for statistical software for the RA2b analysis.

   The global likelihood is in ra2bRoostatsClass.c/h and is described
   in ra2b-roostats-class.tex

   If you are running at CERN, do this to get the right version of root:

   cmsenv
   source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.csh
   source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.28.00c/x86_64-slc5-gcc43-opt/root/bin/thisroot.csh

   If you are running on your own computer, you must have a release of root
   with roofit support built in.  On my Macs, I built from source with

     ./configure macosx --enable-roofit --disable-xrootd



