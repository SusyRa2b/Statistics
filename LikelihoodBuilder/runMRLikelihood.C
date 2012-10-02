{

  gROOT->ProcessLine(".x setup.C");
  
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  // THIS GREATLY SIMPLIFIES THE TAU->HAD PORTION

  gROOT->ProcessLine(".L metReweightingBuilderSIMPLETAU.C++");
  // inputs are: output root file, input list of bins and their content
  buildMRLikelihood( "testoutput.root", "MRSimpleSetupFile.txt" );
 
}
