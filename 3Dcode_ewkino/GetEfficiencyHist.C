#include <iostream>
#include <fstream>

void GetEfficiencyHist(TString inFile, TString outFile) {

  gROOT->Reset() ;

  ifstream infp ;
  infp.open(inFile) ;

  TH1D *heff = new TH1D("heff","efficiency",15,137.5,512.5);

  while ( infp.good() ) {

    //const int ArraySize = 219 ;
    const int ArraySize = 291 ;
    double ArrayContent[ArraySize] ;

    for (int i = 0; infp && i < ArraySize; ++ i) {
      infp >> ArrayContent[i];
    }

    int mgl = ArrayContent[0] ;
    int mlsp = ArrayContent[1] ;
    int ngen = ArrayContent[2] ;

    cout << "ngen = " << ngen << endl ;

    double sig_yield = 0. ;

    for (int i = 3; i < 39 ; i++ ) {
      sig_yield += ArrayContent[i] ;
    }

    double efficiency = sig_yield / ngen ;

    int binX = int(mgl/25) - 5 ;

    heff->SetBinContent(binX,efficiency);

  }

  heff->Draw("text");
  
  TFile *epsF = new TFile(outFile,"recreate");
  epsF->cd();
  heff->Write();
  epsF->Close();

  return ;

}
