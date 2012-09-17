#ifndef MR
#define MR


void makeBigBin( RooWorkspace& wspace, TString binname, ifstream &scalefactorsFile, ifstream &sigfracFile );


void makeDileptonConstraintsPrediction( RooWorkspace& wspace, TString binname, ifstream &scaleFactorsFile, ifstream &sigfracFile );



void makeOutsideUnbinnedConstraints( RooWorkspace& wspace, int nbs, double oneTightMusigfrac, double twoTightMusigfrac, double twoTightMuLooseMusigfrac, double twoTightMuLooseEsigfrac );


void makeTauHadBinPrediction( RooWorkspace& wspace, TString binname, vector<TString> allbinnames, double sfactor1Tau, double sfactor2Tau, double sfactorMuTau, double sfactorETau );


void buildMRLikelihood( TString outputFile, TString bincontentFileName, TString outsidebincontentFileName, TString scalefactorsFileName, TString tauhadscalefactorsFileName, TString signalfractionsFileName, TString outsidesignalfractionsFileName, bool standaloneMode=false, RooWorkspace* ws = 0 );

#endif
