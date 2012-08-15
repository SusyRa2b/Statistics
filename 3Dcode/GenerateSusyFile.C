#include <iostream>
#include <fstream>
#include <sstream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TROOT.h"

void GenerateSusyFile() {

//gluino cross sections in pb.

//  TFile prospino("referenceXSecs.root");
//  TH1F *gluinoxsec = gluino;

  TFile prospino("referenceXSecs.root");
  ///// TH1F *gluinoxsec = (TH1F*) prospino.Get("gluino") ;
  TH1F *gluinoxsec = (TH1F*) prospino.Get("gluino_NLONLL") ;
  if ( gluinoxsec == 0 ) {
     printf("\n\n *** Can't find histogram with name gluino in referenceXSecs.root.\n\n") ;
     return ;
  } else {
     gluinoxsec->Print("all") ;
  }

// bin 1 = gluino mass = 100 GeV, 2 = 125, 3 = 150, 4 = 175, ...
// so gluino mass = 75+nbin*25; or nbin = (gluinomass-75)/25.

  TChain chainT1tttt("tree");
  chainT1tttt.Add("files5fb_MT/T1tttt.root");

  gROOT->Reset();

  const int nBinsBjets = 3 ;   // this must always be 3
  const int nJetsCut = 3 ;     // #jets >= nJetsCut

  //-- met2-ht1-v1
//const int nBinsMET   = 2 ;
//const int nBinsHT    = 1 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,250.,99999.};
//float Hbins[nBinsHT+1] = {400.,99999.};

////-- met2-ht2-v1
//const int nBinsMET   = 2 ;
//const int nBinsHT    = 2 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,250.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,99999.};

////-- met2-ht8-v1
//const int nBinsMET   = 2 ;
//const int nBinsHT    = 8 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,250.,99999.};
//float Hbins[nBinsHT+1] = {400.,500.,600.,700.,800.,900.,1000.,1200.,99999.};

  //-- met3-ht2-v1
//const int nBinsMET   = 3 ;
//const int nBinsHT    = 2 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,250.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,800.,99999.};

//-- met3-ht3-v1
const int nBinsMET   = 3 ;
const int nBinsHT    = 3 ;
    const int version = 1;
float Mbins[nBinsMET+1] = {150.,250.,350.,99999.};
float Hbins[nBinsHT+1] = {400.,600.,1000.,99999.};

//-- met3-ht3-v2
//const int nBinsMET   = 3 ;
//const int nBinsHT    = 3 ;
//    const int version = 2;
//float Mbins[nBinsMET+1] = {150.,250.,350.,99999.};
//float Hbins[nBinsHT+1] = {300.,500.,1000.,99999.};

////-- met3-ht4-v1
//const int nBinsMET   = 3 ;
//const int nBinsHT    = 4 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,250.,350.,99999.};
//float Hbins[nBinsHT+1] = {200, 300.,500.,1000.,99999.};

////-- met3-ht5-v1
//const int nBinsMET   = 3 ;
//const int nBinsHT    = 5 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,250.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,500.,600.,800.,1000.,99999.};

////-- met4-ht3-v1
//const int nBinsMET   = 4 ;
//const int nBinsHT    = 3 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,250.,350.,500.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,1000.,99999.};

////-- met4-ht3-v2
//const int nBinsMET   = 4 ;
//const int nBinsHT    = 3 ;
//    const int version = 2;
//float Mbins[nBinsMET+1] = {125, 150.,250.,350.,99999.};
//float Hbins[nBinsHT+1] = {300.,500.,1000.,99999.};

  //-- met5-ht4-v1
//const int nBinsMET   = 5 ;
//const int nBinsHT    = 4 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,800.,1000.,99999.};

  //-- met4-ht4-v1
  //   const int nBinsMET   = 4 ;
  //  const int nBinsHT	 = 4 ;
  //  const int version = 1;
  //  float Mbins[nBinsMET+1] = {150.,200.,250.,300.,99999.};
  //  float Hbins[nBinsHT+1] = {400.,500.,600.,800.,99999.};

  //-- met4-ht4-v2
//    const int nBinsMET   = 4 ;
//    const int nBinsHT    = 4 ;
//    const int version = 2;
//    float Mbins[nBinsMET+1] = {150.,250.,350.,450.,99999.};
//    float Hbins[nBinsHT+1] = {400.,500.,600.,800.,99999.};

  //-- met4-ht4-v3
//    const int nBinsMET   = 4 ;
//    const int nBinsHT    = 4 ;
//    const int version = 3;
//    float Mbins[nBinsMET+1] = {150.,200.,250.,300.,99999.};
//    float Hbins[nBinsHT+1] = {400.,600.,800.,1000.,99999.};

  //-- met4-ht4-v4
//    const int nBinsMET   = 4 ;
//    const int nBinsHT    = 4 ;
//    const int version = 4;
//    float Mbins[nBinsMET+1] = {150.,250.,350.,450.,99999.};
//    float Hbins[nBinsHT+1] = {400.,600.,800.,1000.,99999.};

  //-- met4-ht4-v5
//    const int nBinsMET   = 4 ;
//    const int nBinsHT	 = 4 ;
//    const int version = 5;
//    float Mbins[nBinsMET+1] = {150.,175.,200.,400.,99999.};
//    float Hbins[nBinsHT+1] = {400.,450.,550.,850.,99999.};

  //-- met4-ht4-v6
//    const int nBinsMET   = 4 ;
//    const int nBinsHT    = 4 ;
//    const int version = 6;
//    float Mbins[nBinsMET+1] = {150.,200.,350.,450.,99999.};
//    float Hbins[nBinsHT+1] = {400.,550.,800.,950.,99999.};

////-- met4-ht4-v7
//const int nBinsMET   = 4 ;
//const int nBinsHT    = 4 ;
//    const int version = 7;
//float Mbins[nBinsMET+1] = {125, 150.,250.,350.,99999.};
//float Hbins[nBinsHT+1] = {200, 300.,500.,1000.,99999.};

////-- met4-ht5-v1
//const int nBinsMET   = 4 ;
//const int nBinsHT    = 5 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,200.,300.,400.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,800.,1000.,1200.,99999.};

////-- met5-ht5-v1
//const int nBinsMET   = 5 ;
//const int nBinsHT    = 5 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,500.,600.,800.,1000.,99999.};

////-- met5-ht3-v1
//const int nBinsMET   = 5 ;
//const int nBinsHT    = 3 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,1000.,99999.};

////-- met5-ht3-v2
//const int nBinsMET   = 5 ;
//const int nBinsHT    = 3 ;
//    const int version = 2;
//float Mbins[nBinsMET+1] = {150.,175.,200.,350.,450.,99999.};
//float Hbins[nBinsHT+1] = {400.,550.,800.,99999.};

  //-- met6-ht6-v1
//const int nBinsMET   = 6 ;
//const int nBinsHT    = 6 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,400.,99999.};
//float Hbins[nBinsHT+1] = {400.,500.,600.,700.,800.,900.,99999.};

  //-- met7-ht7-v1
//const int nBinsMET   = 7 ;
//const int nBinsHT    = 7 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,400.,500.,99999.};
//float Hbins[nBinsHT+1] = {400.,500.,600.,700.,800.,900.,1000.,99999.};

////-- met8-ht8-v1
//const int nBinsMET   = 8 ;
//const int nBinsHT    = 8 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,400.,450.,600.,99999.};
//float Hbins[nBinsHT+1] = {400.,500.,600.,700.,800.,900.,1000.,1200.,99999.};

////-- met8-ht2-v1
//const int nBinsMET   = 8 ;
//const int nBinsHT    = 2 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,400.,450.,600.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,99999.};



//TString sMbins[nBinsMET];
//TString sHbins[nBinsHT];
//TString sBbins[3] = {"_1b","_2b","_3b"};

//TString cMbins[nBinsMET];
//TString cHbins[nBinsHT];
  ////// TString cBbins[3] = {"&&nB==1","&&nB==2","&&nB>=3"};
  TString cBbins[3] = {"nB==1","nB==2","nB>=3"};


  // dummy masses
  int minGlMass = 950 ;
  int maxGlMass = 960 ;


//double dummyYield = 9.9 ;
  double dummyCorr = 1. ;
  double dummyErr = 10. ; // error is in % 
  long dummyEvts = 10000 ;

  ofstream inFile;
  char outfile[10000] ;
  sprintf( outfile, "datfiles/T1tttt-met%d-ht%d-v%d.dat", nBinsMET, nBinsHT, version ) ;
  inFile.open( outfile );

  // loop over gluino masses

  TH1F* ht = new TH1F("ht","ht",10,0,10000);
  TString cutsSig = "minDelPhiN>4&&nMu==0&&nEl==0&&";
  TString cutsSL =  "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nEl==1&&nMu==0) )&&";
  TString cutsLDP = "minDelPhiN<4&&nMu==0&&nEl==0&&";


  TH2F* h_susy_sig[10] ;
  TH2F* h_susy_sl[10] ;
  TH2F* h_susy_ldp[10] ;
  for ( int bi=0; bi<nBinsBjets; bi++ ) {
     char hname[1000] ;
     sprintf( hname, "h_susy_sig_%db", bi+1 ) ;
     h_susy_sig[bi] = new TH2F( hname, hname, nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_susy_sig[bi] -> Sumw2() ;
     sprintf( hname, "h_susy_sl_%db", bi+1 ) ;
     h_susy_sl[bi] = new TH2F( hname, hname, nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_susy_sl[bi] -> Sumw2() ;
     sprintf( hname, "h_susy_ldp_%db", bi+1 ) ;
     h_susy_ldp[bi] = new TH2F( hname, hname, nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_susy_ldp[bi] -> Sumw2() ;
  }


  stringstream njcut ; njcut << nJetsCut;
  TString cutsNjets = "&&nJets>=";
  cutsNjets += njcut.str();

  float xsec = -1.;
  for ( int mGl = minGlMass ; mGl < maxGlMass ; mGl = mGl + 25 ) {

    int theBin = gluinoxsec->FindBin( mGl ) ;					      
    if ( theBin <=0 || theBin > gluinoxsec->GetNbinsX() ) {			      
       printf("\n\n *** can't find bin for mgl=%d.  Returned %d\n\n", mGl, theBin ) ; 
       return ; 								      
    }										      
    xsec = gluinoxsec->GetBinContent( theBin ) ;				      

    ////// for ( int mLsp = 50 ; mLsp < ( mGl - 25 ) ; mLsp = mLsp + 25 ) {
    for ( int mLsp = 550 ; mLsp < 560 ; mLsp = mLsp + 25 ) {

      inFile << mGl << " " << mLsp << " " << dummyEvts << " " ;
      printf(" mGl=%4d, mLsp=%4d\n", mGl, mLsp ) ; cout << flush ;


      printf("\n\n") ;
      for (int k = 0 ; k < nBinsBjets ; k++) {

         TString cutSMS = "mgluino>";
         cutSMS += mGl-1;
         cutSMS += "&&mgluino<";
         cutSMS += mGl+1;
         cutSMS += "&&mlsp>";
         cutSMS += mLsp-1;
         cutSMS += "&&mlsp<";
         cutSMS += mLsp+1;
         cutSMS += "&&";

         TString cut = cBbins[k] ;

         TString allSigCuts = cutSMS+cutsSig+cut+cutsNjets ;
         TString allSLCuts  = cutSMS+cutsSL+cut+cutsNjets ;
         TString allLDPCuts = cutSMS+cutsLDP+cut+cutsNjets ;

         char hname[1000] ;

         sprintf( hname, "h_susy_sig_%db", k+1 ) ;
         chainT1tttt.Project( hname,"HT:MET",allSigCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  SIG selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_sig[k]->Integral() ) ; cout << flush ;

         sprintf( hname, "h_susy_sl_%db", k+1 ) ;
         chainT1tttt.Project( hname,"HT:MET",allSLCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  SL  selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_sl[k]->Integral() ) ; cout << flush ;

         sprintf( hname, "h_susy_ldp_%db", k+1 ) ;
         chainT1tttt.Project( hname,"HT:MET",allLDPCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  LDP selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_ldp[k]->Integral() ) ; cout << flush ;

      } // k (nBjets)
      printf("\n\n") ;

      float totalSUSYyield = 0;
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

         // each point has 10k events generated. The sigma is in pb and I want to normalized to 5 fb-1. 
         // so multiple cross section by 0.5 to get events in 5 fb-1

             inFile << 0.5*xsec*(h_susy_sig[k] -> GetBinContent( i+1, j+1 )) << " " ;
             inFile << 0.5*xsec*(h_susy_sl[k]  -> GetBinContent( i+1, j+1 )) << " " ;
             inFile << 0.5*xsec*(h_susy_ldp[k] -> GetBinContent( i+1, j+1 )) << " " ;

             printf ( " mGl=%d, mLsp=%d: MET,HT (%d,%d) SIG = %9.1f, SL=%9.1f, LDP=%9.1f\n",
                 mGl, mLsp, i+1, j+1, 
                 h_susy_sig[k] -> GetBinContent( i+1, j+1 ),
                 h_susy_sl[k]  -> GetBinContent( i+1, j+1 ),
                 h_susy_ldp[k] -> GetBinContent( i+1, j+1 )  ) ;
             totalSUSYyield += (h_susy_sig[k] -> GetBinContent( i+1, j+1 )*0.5*xsec);
          } // k
        } // j
      } // i
      printf("Total SUSY yield within current binning = %9.1f", totalSUSYyield);
      printf("\n\n") ;
      


  //----------------------------------------------------------------------------

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

    } // mLsp
  } // mGl


  // print out header line:
  inFile << "Using HT bins:  " ;
  for (int j = 0 ; j <= nBinsHT ; j++ ) {
    inFile << Hbins[j] ;
    if ( j < nBinsHT ) inFile << "-" ;
  }

  inFile << "\t Using MET bins: " ;
  for (int i = 0 ; i <= nBinsMET ; i++ ) {
    inFile << Mbins[i] ;
    if ( i < nBinsMET ) inFile << "-" ;
  }
  inFile << endl ;


  return;

}
