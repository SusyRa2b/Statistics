#include <iostream>

void GenerateDummyDeffDbtag() {

  gROOT->Reset();

  int nBinsMET   = 2 ;
  int nBinsHT    = 2 ;
  int nBinsBjets = 3 ;   // this must always be 3

  // dummy masses
  int minGlMass = 200 ;
  int maxGlMass = 400 ;
  int mLsp = 0 ;

  // derivatives for the 1b, 2b, 3b cases:
  double dummyDer_0l[3] = {-0.1,0.0,0.1} ;
  double dummyDer_1l[3] = {-0.05,0.0,0.05} ;
  double dummyDer_ldp[3] = {-0.01,0.0,0.01} ;
  
  ofstream inFile;
  inFile.open("dummy_DeffDbtag.dat");

  // loop over gluino masses

  for ( int mGl = minGlMass ; mGl < maxGlMass ; mGl = mGl + 25 ) {
    for ( int mLsp = 50 ; mLsp < ( mGl - 25 ) ; mLsp = mLsp + 25 ) {

      inFile << mGl << " " << mLsp << " " ;

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBjets ; k++) {
	    inFile << dummyDer_0l[k] << " " ;
	  }
	}
      }

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBjets ; k++) {
	    inFile << dummyDer_1l[k] << " " ;
	  }
	}
      }

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBjets ; k++) {
	    inFile << dummyDer_ldp[k] << " " ;
	  }
	}
      }

      inFile << endl ;

    }
  }



  return;

}
