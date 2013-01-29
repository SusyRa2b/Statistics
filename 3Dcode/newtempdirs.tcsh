#!/bin/tcsh

  set savename = $1


  if ( $#argv == 0 ) then
     printf "\n\n Usage:  newtmpdirs.tcsh  <savename>\n\n"
     exit
  endif

  printf "\n\n savename is %s\n\n" $savename

  mv datfiles    datfiles-$savename
  mv logfiles    logfiles-$savename
  mv outputfiles outputfiles-$savename
  mv rootfiles   rootfiles-$savename

  mkdir datfiles
  mkdir logfiles
  mkdir outputfiles
  mkdir rootfiles


