#!/bin/tcsh -f

  # the script takes as arguments the gluino and the LSP masses,
  # and the number of signal events to start the scan from

  mkdir -p results
  mkdir -p workspaces


  # first make the ws

  setenv wsFname workspaces/ws-T1bbbb-$1-$2.root

  $ROOTSYS/bin/root -b -q make_ws_rootfile3b.c\(\"datfiles_19fb/data-vals-unblind.dat\",\"datfiles_18fb/sigcounts.T1bbbbMG.txt\",$1,$2,\"datfiles_18fb/btagEff_T1bbbb_madgraph_SIG_SL_LDP.txt\",\"datfiles_18fb/lightMistag_T1bbbb_madgraph_SIG_SL_LDP.txt\",4,\"$wsFname\",\"null\",\"datfiles_18fb/sigsystematics.T1bbbbMG.JES.txt\"\,\"datfiles_18fb/sigsystematics.T1bbbbMG.PDF.txt\",\"datfiles_18fb/sigsystematics.T1bbbbMG.ISR.txt\",\"datfiles_18fb/wjets-xsec-shapesyst.txt\",\"datfiles_18fb/singletop-xsec-shapesyst.txt\"\)

  setenv fname results/T1bbbb_CLsAsy_$1_$2.txt

  
  # then submit the asymptotic CLs scan

  set nPoints  = $3
  set startVal = $4  
  set stopVal  = $5

  $ROOTSYS/bin/root -l -b -q RunSimpleAsyScan.C\(\"$wsFname\",$nPoints,$startVal,$stopVal,\"$fname\"\)


