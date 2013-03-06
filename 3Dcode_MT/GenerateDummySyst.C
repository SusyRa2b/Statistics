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

  int nBinsMET   = 4 ;
  int nBinsHT    = 4 ;
  int nBinsBjets = 3 ;

  int nBins = nBinsMET * nBinsHT * nBinsBjets ;

  // just generate one line, gluino mass, lsp mass (for the point 1225, 225),
  // followed by 4*Nbins dummy numbers 

  double dummySig   = 0.01 ;
  double dummySLSig = 0.015 ;
  double dummySL    = 0.005 ;
  double dummyLdp   = 0.02 ;

  ofstream outDummy;
  outDummy.open("DummySyst.txt");

  outDummy << "1225 " << "225 " ;

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


  return;

}
