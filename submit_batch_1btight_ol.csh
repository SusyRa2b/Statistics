#!/bin/tcsh

#
# Args: mg mlsp Xsec iter
#

  which root

  echo args $1 $2 $3

  mkdir -p output-files

  $ROOTSYS/bin/root -b -q 'cls-scripts-v2/runcls_ln_ge1btight_t1bbbb_7_ol.c('`echo $1`','`echo $2`','`echo $3`')'

  echo all done.
