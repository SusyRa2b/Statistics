#ifndef METREWEIGHTINGBUILDERSIMPLETAU_H
#define METREWEIGHTINGBUILDERSIMPLETAU_H


void makePrediction( RooWorkspace& wspace, TString thisBin, bool standalone );

void makePolarizationConstraintsPredictions( RooWorkspace& wspace, TString binname, TString trigeffname );

void makeDileptonConstraintsPredictions( RooWorkspace& wspace, TString binname,  TString binname_outside, TString trigeffname );

void makeTauHadBinPrediction( RooWorkspace& wspace, TString binname, TString binname_outside, TString trigeffname );

void buildMRLikelihood( RooWorkspace& wspace, TString outputFile, TString setupFileName, bool standalone ); 


#endif
