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
  TH1F *gluinoxsec = (TH1F*) prospino.Get("gluino") ;
  if ( gluinoxsec == 0 ) {
     printf("\n\n *** Can't find histogram with name gluino in referenceXSecs.root.\n\n") ;
     return ;
  } else {
     gluinoxsec->Print("all") ;
  }

// bin 1 = gluino mass = 100 GeV, 2 = 125, 3 = 150, 4 = 175, ...
// so gluino mass = 75+nbin*25; or nbin = (gluinomass-75)/25.

  TChain chainT1bbbb("tree");
  chainT1bbbb.Add("files5fb_lowHT/T1bbbb.root");

  gROOT->Reset();

  const int nBinsBjets = 3 ;   // this must always be 3
  const int nJetsCut = 3 ;     // #jets >= nJetsCut


  //-- met2-ht1-v1
//const int nBinsMET   = 2 ;
//const int nBinsHT    = 1 ;
//float Mbins[nBinsMET+1] = {150.,250.,99999.};
//float Hbins[nBinsHT+1] = {400.,99999.};

////-- met2-ht2-v1
//const int nBinsMET   = 2 ;
//const int nBinsHT    = 2 ;
//float Mbins[nBinsMET+1] = {150.,250.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,99999.};


  //-- met3-ht3-v1
//const int nBinsMET   = 3 ;
//const int nBinsHT    = 3 ;
//float Mbins[nBinsMET+1] = {150.,250.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,1000.,99999.};

    //-- met4-ht4-v1
    const int nBinsMET   = 4 ;
    const int nBinsHT    = 4 ;
    float Mbins[nBinsMET+1] = {150.,200.,250.,300.,99999.};
    float Hbins[nBinsHT+1] = {400.,500.,600.,800.,99999.};

////-- met5-ht5-v1
//const int nBinsMET   = 5 ;
//const int nBinsHT    = 5 ;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,500.,600.,800.,1000.,99999.};



//TString sMbins[nBinsMET];
//TString sHbins[nBinsHT];
//TString sBbins[3] = {"_1b","_2b","_3b"};

//TString cMbins[nBinsMET];
//TString cHbins[nBinsHT];
  ////// TString cBbins[3] = {"&&nB==1","&&nB==2","&&nB>=3"};
  TString cBbins[3] = {"nB==1","nB==2","nB>=3"};


  // dummy masses
  int minGlMass = 850 ;
  int maxGlMass = 910 ;


//double dummyYield = 9.9 ;
  double dummyCorr = 1. ;
  double dummyErr = 10. ; // error is in % 
  long dummyEvts = 10000 ;

  ofstream inFile;
  inFile.open("Susy.dat");

  // loop over gluino masses

  TH1F* ht = new TH1F("ht","ht",10,0,10000);
  TString cutsSig = "minDelPhiN>4&&nMu==0&&nEl==0&&";
  TString cutsSL = "minDelPhiN>4&&(nMu==1||nEl==1)&&";
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
  xsec = gluinoxsec->GetBinContent((mGl-75)/25);
    ////// for ( int mLsp = 50 ; mLsp < ( mGl - 25 ) ; mLsp = mLsp + 25 ) {
    for ( int mLsp = 200 ; mLsp < 310 ; mLsp = mLsp + 25 ) {

      inFile << mGl << " " << mLsp << " " << dummyEvts << " " ;
      printf(" mGl=%4d, mLsp=%4d\n", mGl, mLsp ) ; cout << flush ;

  //----------------------------------------------------------------------------
  /// for (int i = 0 ; i < nBinsMET ; i++) {
  ///   for (int j = 0 ; j < nBinsHT ; j++) {
  ///     for (int k = 0 ; k < nBinsBjets ; k++) {
  ///    // the three entries are the signal region, 1L region and very low met sideband

  ///       TString cutSMS = "mgluino>";
  ///       cutSMS += mGl-1;
  ///       cutSMS += "&&mgluino<";
  ///       cutSMS += mGl+1;
  ///       cutSMS += "&&mlsp>";
  ///       cutSMS += mLsp-1;
  ///       cutSMS += "&&mlsp<";
  ///       cutSMS += mLsp+1;
  ///       cutSMS += "&&";

  ///       TString cut = "HT>";
  ///       cut += Hbins[j];
  ///       cut += "&&HT<";
  ///       cut += Hbins[j+1];
  ///       cut += "&&MET>";
  ///       cut += Mbins[i];
  ///       cut += "&&MET<";
  ///       cut += Mbins[i+1];
  ///       cut += cBbins[k];

  ///       // cross section in pb = gluinoxsec->GetBinContent((mGl-75)/25);
  ///       //cout << " bin = " << (mGl-75)/25 << " for mGl = " << mGl << " and xsec = " << gluinoxsec->GetBinContent((mGl-75)/25) << endl;

  ///       // each point has 10k events generated. The sigma is in pb and I want to normalized to 5 fb-1. 
  ///       // so multiple cross section by 0.5 to get events in 5 fb-1

  ///       TString allSigCuts = cutSMS+cutsSig+cut+cutsNjets ;
  ///       TString allSLCuts  = cutSMS+cutsSL+cut+cutsNjets ;
  ///       TString allLDPCuts = cutSMS+cutsLDP+cut+cutsNjets ;

  ///       chainT1bbbb.Project("ht","HT",allSigCuts);
  ///       double nselSig = 0.5*xsec*ht->GetSumOfWeights() ;
  ///       inFile << nselSig << " ";
  ///       ht->Reset() ;
  ///       printf("  mGl=%4d, mLsp=%4d : bins in met,HT,nB (%d,%d,%d) : nSig = %7.1f : cuts=%s\n",
  ///            mGl, mLsp, i, j, k, nselSig, allSigCuts.Data() ) ; cout << flush ;

  ///       chainT1bbbb.Project("ht","HT",allSLCuts);
  ///       double nselSL = 0.5*xsec*ht->GetSumOfWeights() ;
  ///       inFile << nselSL << " ";
  ///       ht->Reset() ;
  ///       printf("  mGl=%4d, mLsp=%4d : bins in met,HT,nB (%d,%d,%d) : nSL  = %7.1f : cuts=%s\n",
  ///            mGl, mLsp, i, j, k, nselSL, allSLCuts.Data() ) ; cout << flush ;

  ///       chainT1bbbb.Project("ht","HT",allLDPCuts);
  ///       double nselLDP = 0.5*xsec*ht->GetSumOfWeights() ;
  ///       inFile << nselLDP << " ";
  ///       ht->Reset() ;
  ///       printf("  mGl=%4d, mLsp=%4d : bins in met,HT,nB (%d,%d,%d) : nLDP = %7.1f : cuts=%s\n",
  ///            mGl, mLsp, i, j, k, nselLDP, allLDPCuts.Data() ) ; cout << flush ;

  ///     }
  ///   }
  /// }
  //-------- slow way above, faster way below.  --------------------------------

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
         chainT1bbbb.Project( hname,"HT:MET",allSigCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  SIG selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_sig[k]->Integral() ) ; cout << flush ;

         sprintf( hname, "h_susy_sl_%db", k+1 ) ;
         chainT1bbbb.Project( hname,"HT:MET",allSLCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  SL  selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_sl[k]->Integral() ) ; cout << flush ;

         sprintf( hname, "h_susy_ldp_%db", k+1 ) ;
         chainT1bbbb.Project( hname,"HT:MET",allLDPCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  LDP selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_ldp[k]->Integral() ) ; cout << flush ;

      } // k (nBjets)
      printf("\n\n") ;


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

          } // k
        } // j
      } // i
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
