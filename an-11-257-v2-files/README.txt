

  August 6, 2011

  Owen Long, UCR

  This directory and its subdirectories contains all of the input files and
  analysis macros needed to run the RA2b likelihood analysis and produce
  susy exclusion plots (though currently only for the old susy scan MC).

  In everything below, I'll assume you are on a unix/linux machine, that
  you have the path to root set up, and that you are using tcsh (not bash).
  If your shell is not tcsh, just type "tcsh" to start one.

  First, set up the output directories by doing this command from this
  directory

    source sourceme-setup-directories


  All other things should be run from the directory below this one,
  which is UserCode/SusyAnalysis/RA2b/Statistics.


  To run the likelihood analysis using the data inputs

    source an-11-257-v2-files/analysis-macros/sourceme-data-results


  To run the likelihood analysis using the MC inputs

    source an-11-257-v2-files/analysis-macros/sourceme-mc-results


  To run the susy scans

    source an-11-257-v2-files/analysis-macros/sourceme-scans


