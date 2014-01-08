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

  int in_version, in_nBinsVar1, in_nBinsVar2, in_nBinsBjets ;
  TString label, in_Var1, in_Var2, in_Var3 ;

  ifstream inBinning ;
  inBinning.open("Binning.txt") ;

  inBinning >> label >> in_Var1 ;
  inBinning >> label >> in_Var2 ;
  inBinning >> label >> in_Var3 ;

  inBinning >> label >> in_nBinsVar1 ;
  inBinning >> label >> in_nBinsVar2 ;
  inBinning >> label >> in_nBinsBjets ;

  const int nBinsVar1  = in_nBinsVar1 ;
  const int nBinsVar2  = in_nBinsVar2 ;
  const int nBinsBjets = in_nBinsBjets ;

  inBinning.close() ;

  TString aVar1 = in_Var1 ;
  TString aVar2 = in_Var2 ;
  if ( aVar1 == "MET/sqrt(HT)" ) aVar1 = "MET_div_sqrtHT" ;
  if ( aVar2 == "MET/sqrt(HT)" ) aVar2 = "MET_div_sqrtHT" ;

  int nBins = nBinsVar1 * nBinsVar2 * nBinsBjets ;

  double dummySig   = 0.01 ;
  double dummySLSig = 0.015 ;
  double dummySL    = 0.005 ;
  double dummyLdp   = 0.02 ;

  char outfile[10000] ;
  sprintf( outfile, "datfiles/DummySyst-%s-%d-%s-%d-nB%d.dat", aVar1.Data(), nBinsVar1, aVar2.Data(), nBinsVar2, nBinsBjets ) ;

  ofstream outDummy;
  outDummy.open(outfile);


  // generate dummy systematics for a bunch of points 

  int mGls[9] = {150,200,250,300,350,400,450,500,550};
  int mLsps[9] = {0,0,0,0,0,0,0,0,0};


  outDummy << -1 << " " << -1 << " " ;

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


  for ( int iGl = 0 ; iGl < 9 ; iGl++ ) {

    int mGl  = mGls[iGl];
    int mLsp = mLsps[iGl];

    outDummy << mGl << " " << mLsp << " " ;

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
