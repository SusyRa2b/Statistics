#include <iostream>

void GenerateDummySusyFile() {

  gROOT->Reset();

  int nBinsMET   = 3 ;
  int nBinsHT    = 3 ;
  int nBinsBjets = 3 ;   // this must always be 3

  // dummy masses
  int minGlMass = 200 ;
  int maxGlMass = 400 ;
  int mLsp = 0 ;

  double dummyYield = 9.9 ;
  double dummyCorr = 1. ;
  double dummyErr = 10. ; // error is in %
  long dummyEvts = 10000 ;
  
  ofstream inFile;
  inFile.open("dummy_Susy.dat");

  // loop over gluino masses

  for ( int mGl = minGlMass ; mGl < maxGlMass ; mGl = mGl + 25 ) {
    for ( int mLsp = 50 ; mLsp < ( mGl - 25 ) ; mLsp = mLsp + 25 ) {

      inFile << mGl << " " << mLsp << " " << dummyEvts << " " ;

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBjets ; k++) {
	    inFile << dummyYield << " " << dummyYield << " " << dummyYield << " " ;
	  }
	}
      }

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBjets ; k++) {
	    inFile << dummyCorr << " " << dummyCorr << " " << dummyCorr << " " ;
	  }
	}
      }

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBjets ; k++) {
	    inFile << dummyErr << " " << dummyErr << " " << dummyErr << " " ;
	  }
	}
      }

      inFile << endl ;

    }
  }



  return;

}
