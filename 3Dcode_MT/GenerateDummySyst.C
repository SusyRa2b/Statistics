#include <iostream>
#include <fstream>
#include <sstream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TROOT.h"

void GenerateDummySyst() {

  gROOT->Reset();

  // get binning from "Binning.txt"

  int in_version, in_nBinsMET, in_nBinsHT, in_nBinsBjets ;
  TString label ;

  ifstream inBinning ;
  inBinning.open("Binning.txt") ;

  inBinning >> label >> in_nBinsMET ;
  inBinning >> label >> in_nBinsHT ;
  inBinning >> label >> in_nBinsBjets ;

  const int nBinsMET   = in_nBinsMET ;
  const int nBinsHT    = in_nBinsHT ;
  const int nBinsBjets = in_nBinsBjets ;

  inBinning.close() ;


  int nBins = nBinsMET * nBinsHT * nBinsBjets ;

  // just generate one line, gluino mass, lsp mass (for the point 1000, 0),
  // followed by 4*Nbins dummy numbers 

  double dummySig   = 0.01 ;
  double dummySLSig = 0.015 ;
  double dummySL    = 0.005 ;
  double dummyLdp   = 0.02 ;

  char outfile[10000] ;
  sprintf( outfile, "datfiles/DummySyst-met%d-ht%d-nB%d.dat", nBinsMET, nBinsHT, nBinsBjets ) ;

  ofstream outDummy;
  outDummy.open(outfile);


  // generate dummy systematics for a bunch of points 

  int mGls[5] = {-1,350,500,600,700} ;
  int mLsps[5] = {-1,0,0,0,0} ;

  for ( int iGl = 0 ; iGl < 5 ; iGl++ ) {

    outDummy << mGls[iGl] << " " << mLsps[iGl] << " " ;

    for ( int i = 0 ; i < nBins ; i++ ) {
      outDummy << dummySig << " " ;
    }
    
    for ( int i = 0 ; i < nBins ; i++ ) {
      outDummy << dummySLSig << " " ;
    }
    
    for ( int i = 0 ; i < nBins ; i++ ) {
      outDummy << dummySL << " " ;
    }
    
    for ( int i = 0 ; i < nBins ; i++ ) {
      outDummy << dummyLdp << " " ;
    }

    outDummy << endl ;
    
  }

  return;

}
