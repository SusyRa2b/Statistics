#include <iostream>
#include <fstream>
#include <sstream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TROOT.h"

void GenerateSusyFile( double flatDummyErr = 10. ) {  //-- flat error in %.  If negative, use MC stat err.

//gluino cross sections in pb.


  TFile prospino("referenceXSecs.root");
  TH1F *gluinoxsec = (TH1F*) prospino.Get("gluino_NLONLL") ;
  if ( gluinoxsec == 0 ) {
     printf("\n\n *** Can't find histogram with name gluino_NLONLL in referenceXSecs.root.\n\n") ;
     return ;
  }
  TH1F *gluinoxsec8TeV = (TH1F*) prospino.Get("gluino8TeV_NLONLL") ;
  if ( gluinoxsec8TeV == 0 ) {
     printf("\n\n *** Can't find histogram with name gluino8TeV_NLONLL in referenceXSecs.root.\n\n") ;
     return ;
  }

// bin 1 = gluino mass = 100 GeV, 2 = 125, 3 = 150, 4 = 175, ...
// so gluino mass = 75+nbin*25; or nbin = (gluinomass-75)/25.

  TChain chainT1bbbb("tree");
  chainT1bbbb.Add("filesHCP_53_v3_QCDweights/T1bbbb.root");

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
//const int nBinsMET   = 3 ;
//const int nBinsHT    = 3 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,250.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,1000.,99999.};

////-- met3-ht3-v5
//    const int nBinsMET = 3 ;
//    const int nBinsHT  = 3 ;
//    const int version = 5;
//    float Mbins[nBinsMET+1] = { 125, 200,  350, 99999. } ;
//    float Hbins[nBinsHT+1]  = { 400, 600, 1000, 99999. } ;

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

  //-- met4-ht4-v15
      const int nBinsMET   = 4 ;
      const int nBinsHT    = 4 ;
      const int version = 15;
      float Mbins[nBinsMET+1] = {125.,150.,250.,350.,99999.};
      float Hbins[nBinsHT+1] = {400.,500.,800.,1000.,99999.};

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


  // jet pt thresholds
  double minLeadJetPt = 70. ;
  double min3rdJetPt = 50. ;


  // dummy masses
  int minGlMass = 800 ;
  int maxGlMass = 1210 ;


//double dummyYield = 9.9 ;
  double dummyCorr = 1. ;
  long dummyEvts = 10000 ;

  ofstream inFile;
  char outfile[10000] ;
  sprintf( outfile, "datfiles/T1bbbb-met%d-ht%d-v%d.dat", nBinsMET, nBinsHT, version ) ;
  inFile.open( outfile );

  // loop over gluino masses

  // TH1F* ht = new TH1F("ht","ht",10,0,10000);
  TString cutsSig = "minDelPhiN>4&&nMu==0&&nEl==0&&";
  TString cutsSL =  "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nEl==1&&nMu==0) )&&";
  TString cutsLDP = "minDelPhiN<4&&nMu==0&&nEl==0&&";
  TString jetpt = "pt_1st_leadJet>";
  jetpt+= minLeadJetPt;
  jetpt+= "&&pt_2nd_leadJet>";
  jetpt+= minLeadJetPt;
  jetpt+= "&&pt_3rd_leadJet>";
  jetpt+= min3rdJetPt;
  jetpt+= "&&";
  cutsSig += jetpt;
  cutsSL += jetpt;
  cutsLDP += jetpt;


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
  float xsec8TeV = -1. ;
  for ( int mGl = minGlMass ; mGl < maxGlMass ; mGl = mGl + 100 ) {

    int theBin = gluinoxsec->FindBin( mGl ) ;					      
    if ( theBin <=0 || theBin > gluinoxsec->GetNbinsX() ) {			      
       printf("\n\n *** can't find bin for mgl=%d.  Returned %d\n\n", mGl, theBin ) ; 
       return ; 								      
    }										      
    xsec = gluinoxsec->GetBinContent( theBin ) ;				      

    int theBin8TeV = gluinoxsec8TeV->FindBin( mGl ) ;					      
    if ( theBin8TeV <=0 || theBin8TeV > gluinoxsec8TeV->GetNbinsX() ) {			      
       printf("\n\n *** can't find bin for mgl=%d.  Returned %d\n\n", mGl, theBin ) ; 
       return ; 								      
    }										      
    xsec8TeV = gluinoxsec8TeV->GetBinContent( theBin8TeV ) ;				      

    printf("\n\n  SUSY Xsecs:  7 TeV = %f,  8 TeV = %f\n\n", xsec, xsec8TeV ) ;

    ////// for ( int mLsp = 50 ; mLsp < ( mGl - 25 ) ; mLsp = mLsp + 25 ) {
    for ( int mLsp = 300 ; mLsp < 710 ; mLsp = mLsp + 400 ) {

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

      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

             printf ( " Raw MC counts: mGl=%d, mLsp=%d: MET,HT (%d,%d) nb=%d   SIG = %9.0f, SL=%9.0f, LDP=%9.0f\n",
                 mGl, mLsp, i+1, j+1, k+1,
                 h_susy_sig[k] -> GetBinContent( i+1, j+1 ),
                 h_susy_sl[k]  -> GetBinContent( i+1, j+1 ),
                 h_susy_ldp[k] -> GetBinContent( i+1, j+1 )  ) ;

          } // k
          printf("----------------\n") ;
        } // j
      } // i
      printf("\n\n") ;


      float totalSUSYyield = 0;
      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

         //-- Aug 31, 2012: This is for using the 7TeV T1bbbb to get 8TeV 15 fb normalization.
         //
         // each point has 10k events generated. The sigma is in pb and I want to normalized to 15 fb-1. 
         // so multiple cross section by 1.5 to get events in 15 fb-1

             inFile << 1.5*xsec8TeV*(h_susy_sig[k] -> GetBinContent( i+1, j+1 )) << " " ;
             inFile << 1.5*xsec8TeV*(h_susy_sl[k]  -> GetBinContent( i+1, j+1 )) << " " ;
             inFile << 1.5*xsec8TeV*(h_susy_ldp[k] -> GetBinContent( i+1, j+1 )) << " " ;

             totalSUSYyield += (h_susy_sig[k] -> GetBinContent( i+1, j+1 )*1.5*xsec8TeV);

                double nsel_sig = h_susy_sig[k] -> GetBinContent( i+1, j+1 ) ;
                double nsel_sl  = h_susy_sl[k]  -> GetBinContent( i+1, j+1 ) ;
                double nsel_ldp = h_susy_ldp[k] -> GetBinContent( i+1, j+1 ) ;
                double nevt_err_sig = 1 ;
                double nevt_err_sl  = 1 ;
                double nevt_err_ldp = 1 ;
                if ( nsel_sig > 0. ) { nevt_err_sig = 1.5*xsec8TeV*sqrt(nsel_sig) ; }
                if ( nsel_sl  > 0. ) { nevt_err_sl  = 1.5*xsec8TeV*sqrt(nsel_sl ) ; }
                if ( nsel_ldp > 0. ) { nevt_err_ldp = 1.5*xsec8TeV*sqrt(nsel_ldp) ; }

             printf ( " xsec8TeV weighted events: mGl=%d, mLsp=%d: MET,HT (%d,%d) nb=%d   SIG = %6.1f +/- %4.1f,   SL=%6.1f +/- %4.1f,   LDP=%6.1f +/- %4.1f\n",
                 mGl, mLsp, i+1, j+1, k+1,
                 1.5*xsec8TeV*(h_susy_sig[k] -> GetBinContent( i+1, j+1 )), nevt_err_sig,
                 1.5*xsec8TeV*(h_susy_sl[k]  -> GetBinContent( i+1, j+1 )), nevt_err_sl,
                 1.5*xsec8TeV*(h_susy_ldp[k] -> GetBinContent( i+1, j+1 )), nevt_err_ldp  ) ;

         //-- Aug 31, 2012: This is for 7TeV T1bbbb with 5 fb normalization.
         //
         // each point has 10k events generated. The sigma is in pb and I want to normalized to 5 fb-1. 
         // so multiple cross section by 0.5 to get events in 5 fb-1

      ////   inFile << 0.5*xsec*(h_susy_sig[k] -> GetBinContent( i+1, j+1 )) << " " ;
      ////   inFile << 0.5*xsec*(h_susy_sl[k]  -> GetBinContent( i+1, j+1 )) << " " ;
      ////   inFile << 0.5*xsec*(h_susy_ldp[k] -> GetBinContent( i+1, j+1 )) << " " ;

      ////   totalSUSYyield += (h_susy_sig[k] -> GetBinContent( i+1, j+1 )*0.5*xsec);

      ////      double nsel_sig = h_susy_sig[k] -> GetBinContent( i+1, j+1 ) ;
      ////      double nsel_sl  = h_susy_sl[k]  -> GetBinContent( i+1, j+1 ) ;
      ////      double nsel_ldp = h_susy_ldp[k] -> GetBinContent( i+1, j+1 ) ;
      ////      double nevt_err_sig = 1 ;
      ////      double nevt_err_sl  = 1 ;
      ////      double nevt_err_ldp = 1 ;
      ////      if ( nsel_sig > 0. ) { nevt_err_sig = 0.5*xsec*sqrt(nsel_sig) ; }
      ////      if ( nsel_sl  > 0. ) { nevt_err_sl  = 0.5*xsec*sqrt(nsel_sl ) ; }
      ////      if ( nsel_ldp > 0. ) { nevt_err_ldp = 0.5*xsec*sqrt(nsel_ldp) ; }

      ////   printf ( " Xsec weighted events: mGl=%d, mLsp=%d: MET,HT (%d,%d) nb=%d   SIG = %6.1f +/- %4.1f,   SL=%6.1f +/- %4.1f,   LDP=%6.1f +/- %4.1f\n",
      ////       mGl, mLsp, i+1, j+1, k+1,
      ////       0.5*xsec*(h_susy_sig[k] -> GetBinContent( i+1, j+1 )), nevt_err_sig,
      ////       0.5*xsec*(h_susy_sl[k]  -> GetBinContent( i+1, j+1 )), nevt_err_sl,
      ////       0.5*xsec*(h_susy_ldp[k] -> GetBinContent( i+1, j+1 )), nevt_err_ldp  ) ;

          } // k
          printf("----------------\n") ;
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

      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {
             if ( flatDummyErr >= 0 ) {

                inFile << flatDummyErr << " " << flatDummyErr << " " << flatDummyErr << " " ;

             } else {

                //-- compute approximate stat err.
                //-- This is 100* sig_eff / eff = 100 * [sqrt(sel)/N]/[sel/N] = 100/sqrt(sel).
                double nsel_sig = h_susy_sig[k] -> GetBinContent( i+1, j+1 ) ;
                double nsel_sl  = h_susy_sl[k]  -> GetBinContent( i+1, j+1 ) ;
                double nsel_ldp = h_susy_ldp[k] -> GetBinContent( i+1, j+1 ) ;
                double frerr_sig = 100 ;
                double frerr_sl  = 100 ;
                double frerr_ldp = 100 ;
                if ( nsel_sig > 0. ) { frerr_sig = 100./sqrt(nsel_sig) ; }
                if ( nsel_sl  > 0. ) { frerr_sl  = 100./sqrt(nsel_sl ) ; }
                if ( nsel_ldp > 0. ) { frerr_ldp = 100./sqrt(nsel_ldp) ; }

                inFile << frerr_sig << " " << frerr_sl << " " << frerr_ldp << " " ;

                printf ( " MC sig_eff/eff (%%): mGl=%d, mLsp=%d: MET,HT (%d,%d) nb=%d   SIG = %5.1f,   SL=%5.1f,   LDP=%5.1f\n",
                    mGl, mLsp, i+1, j+1, k+1,
                    frerr_sig, frerr_sl, frerr_ldp ) ;
             }
          }
          printf("----------------\n") ;
        }
      }
      printf("\n\n") ;

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
