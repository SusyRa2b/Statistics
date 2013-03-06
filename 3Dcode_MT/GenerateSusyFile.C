#include <iostream>
#include <fstream>
#include <sstream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TROOT.h"

void GenerateSusyFile( double flatDummyErr = 0.00001 ) {  //-- flat error in %.  If negative, use MC stat err.

//gluino cross sections in pb.


  TFile prospino("referenceXSecs.root");

  TH1F *gluinoxsec8TeV = (TH1F*) prospino.Get("gluino8TeV_NLONLL") ;
  if ( gluinoxsec8TeV == 0 ) {
     printf("\n\n *** Can't find histogram with name gluino8TeV_NLONLL in referenceXSecs.root.\n\n") ;
     return ;
  }

  // bin 1 = gluino mass = 100 GeV, 2 = 125, 3 = 150, 4 = 175, ...
  // so gluino mass = 75+nbin*25; or nbin = (gluinomass-75)/25.

  TChain chainT2tt("tree");
  chainT2tt.Add("filesMoriond_v3/T2tt.root");

  gROOT->Reset();

  const int nJetsCut = 3 ;     // #jets >= nJetsCut
  const int MTbCut   = 0 ;   // MTb cut

  //-- met4-ht4-v15
  const int nBinsMET   = 4 ;
  const int nBinsHT    = 4 ;
  const int nBinsBjets = 3 ;   
  const int version = 15;
  float Mbins[nBinsMET+1] = {125.,150.,250.,350.,99999.};
  float Hbins[nBinsHT+1] = {400.,500.,800.,1000.,99999.};

  TString cBbins[nBinsBjets];
  
  cBbins[0] = "nB==1" ;

  if ( nBinsBjets == 2 ) {
    cBbins[1] = "nB>=2" ;
  }
  else if ( nBinsBjets == 3 ) {
    cBbins[1] = "nB==2" ;    
    cBbins[2] = "nB>=3" ;
  }
  else {
    cout << "\n\n This number of #b-jets bins is not implemented! Exiting ... " << endl ;
    return ;
  }

  // jet pt thresholds
  double minLeadJetPt = 70. ;
  double min3rdJetPt = 50. ;


  // dummy masses
  int minGlMass = 200 ;
  int maxGlMass = 1400 ;


//double dummyYield = 9.9 ;
  double dummyCorr = 1. ;
  long dummyEvts = 10000 ;

  ofstream inFile;
  char outfile[10000] ;
  sprintf( outfile, "datfiles/T2tt-met%d-ht%d-v%d.dat", nBinsMET, nBinsHT, version ) ;
  inFile.open( outfile );

  // loop over gluino masses

  // TH1F* ht = new TH1F("ht","ht",10,0,10000);
  TString cutsSig   = "minDelPhiN>4&&nMu==0&&nEl==0&&";
  TString cutsSLSig = "MT>100&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nEl==1&&nMu==0) )&&";
  TString cutsSL    = "MT<100&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nEl==1&&nMu==0) )&&";
  TString cutsLDP   = "minDelPhiN<4&&nMu==0&&nEl==0&&";

  char commoncuts[10000] ;
  sprintf( commoncuts, "maxChNMultDiff<40&&pfOcaloMET<2.0&&nJets>=%d&&MT_bestCSV>%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f&&pt_3rd_leadJet>%.0f)&&",
           nJetsCut, MTbCut, minLeadJetPt, minLeadJetPt, min3rdJetPt ) ;

  cutsSig   += commoncuts ;
  cutsSLSig += commoncuts ;
  cutsSL    += commoncuts ;
  cutsLDP   += commoncuts ;


  TH2F* h_susy_sig[10] ;
  TH2F* h_susy_slsig[10] ;
  TH2F* h_susy_sl[10] ;
  TH2F* h_susy_ldp[10] ;
  for ( int bi=0; bi<nBinsBjets; bi++ ) {
     char hname[1000] ;
     sprintf( hname, "h_susy_sig_%db", bi+1 ) ;
     h_susy_sig[bi] = new TH2F( hname, hname, nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_susy_sig[bi] -> Sumw2() ;
     sprintf( hname, "h_susy_slsig_%db", bi+1 ) ;
     h_susy_slsig[bi] = new TH2F( hname, hname, nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_susy_slsig[bi] -> Sumw2() ;
     sprintf( hname, "h_susy_sl_%db", bi+1 ) ;
     h_susy_sl[bi] = new TH2F( hname, hname, nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_susy_sl[bi] -> Sumw2() ;
     sprintf( hname, "h_susy_ldp_%db", bi+1 ) ;
     h_susy_ldp[bi] = new TH2F( hname, hname, nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_susy_ldp[bi] -> Sumw2() ;
  }

  float xsec8TeV = -1. ;

  int mGls[4] = {350,450,550,650} ;
  int mLsps[4] = {0,0,0,0} ;

  for ( int iGl = 0 ; iGl < 1/*4*/ ; iGl++ ) {

  int mGl = mGls[iGl] ;
  int mLsp = mLsps[iGl] ;

    int theBin8TeV = gluinoxsec8TeV->FindBin( mGl ) ;					      
    if ( theBin8TeV <=0 || theBin8TeV > gluinoxsec8TeV->GetNbinsX() ) {			      
       printf("\n\n *** can't find bin for mgl=%d.  Returned %d\n\n", mGl, theBin8TeV ) ; 
       return ; 								      
    }										      
    xsec8TeV = gluinoxsec8TeV->GetBinContent( theBin8TeV ) ;				      

    printf("\n\n  SUSY Xsecs:  8 TeV = %f\n\n", xsec8TeV ) ;

    ////// for ( int mLsp = 50 ; mLsp < ( mGl - 25 ) ; mLsp = mLsp + 25 ) {
    //    for ( int mLsp = 300 ; mLsp < 710 ; mLsp = mLsp + 400 ) {



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

         TString allSigCuts   = cutSMS+cutsSig+cut ;
         TString allSLSigCuts = cutSMS+cutsSLSig+cut ;
         TString allSLCuts    = cutSMS+cutsSL+cut ;
         TString allLDPCuts   = cutSMS+cutsLDP+cut ;

         char hname[1000] ;

         sprintf( hname, "h_susy_sig_%db", k+1 ) ;
         chainT2tt.Project( hname,"HT:MET",allSigCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  SIG selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_sig[k]->Integral() ) ; cout << flush ;

         sprintf( hname, "h_susy_slsig_%db", k+1 ) ;
         chainT2tt.Project( hname,"HT:MET",allSLSigCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  SL Sig selection %6.1f events.\n", mGl, mLsp, k+1, h_susy_slsig[k]->Integral() ) ; cout << flush ;

         sprintf( hname, "h_susy_sl_%db", k+1 ) ;
         chainT2tt.Project( hname,"HT:MET",allSLCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  SL  selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_sl[k]->Integral() ) ; cout << flush ;

         sprintf( hname, "h_susy_ldp_%db", k+1 ) ;
         chainT2tt.Project( hname,"HT:MET",allLDPCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  LDP selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_ldp[k]->Integral() ) ; cout << flush ;

      } // k (nBjets)
      printf("\n\n") ;

      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

             printf ( " Raw MC counts: mGl=%d, mLsp=%d: MET,HT (%d,%d) nb=%d   SIG = %9.0f, SLSIG =%9.0f, SL=%9.0f, LDP=%9.0f\n",
                 mGl, mLsp, i+1, j+1, k+1,
                 h_susy_sig[k]  -> GetBinContent( i+1, j+1 ),
                 h_susy_slsig[k]-> GetBinContent( i+1, j+1 ),
                 h_susy_sl[k]   -> GetBinContent( i+1, j+1 ),
                 h_susy_ldp[k]  -> GetBinContent( i+1, j+1 )  ) ;

          } // k
          printf("----------------\n") ;
        } // j
      } // i
      printf("\n\n") ;


      // 0-lep yields

      float totalSUSYyield = 0;
      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

             inFile << 1.5*xsec8TeV*(h_susy_sig[k]  -> GetBinContent( i+1, j+1 )) << " " ;

             totalSUSYyield += (h_susy_sig[k] -> GetBinContent( i+1, j+1 )*1.5*xsec8TeV);
	    	     
          } // k
        } // j
      } // i
      printf("Total SUSY yield within current binning = %9.1f",totalSUSYyield);
      printf("\n\n") ;


      // SLSig yields:
      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

             inFile << 1.5*xsec8TeV*(h_susy_slsig[k]-> GetBinContent( i+1, j+1 )) << " " ;

	  }
	}
      }


      // SL yields:
      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

             inFile << 1.5*xsec8TeV*(h_susy_sl[k]   -> GetBinContent( i+1, j+1 )) << " " ;


	  }
	}
      }


      // Ldp yields:
      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

             inFile << 1.5*xsec8TeV*(h_susy_ldp[k]  -> GetBinContent( i+1, j+1 )) << " " ;

	     double nsel_sig   = h_susy_sig[k]  -> GetBinContent( i+1, j+1 ) ;
	     double nsel_slsig = h_susy_slsig[k]-> GetBinContent( i+1, j+1 ) ;
	     double nsel_sl    = h_susy_sl[k]   -> GetBinContent( i+1, j+1 ) ;
	     double nsel_ldp   = h_susy_ldp[k]  -> GetBinContent( i+1, j+1 ) ;
	     
	     double nevt_err_sig    = 1 ;
	     double nevt_err_slsig  = 1 ;
	     double nevt_err_sl     = 1 ;
	     double nevt_err_ldp    = 1 ;
	     
	     if ( nsel_sig    > 0. ) { nevt_err_sig   = 1.5*xsec8TeV*sqrt(nsel_sig)   ; }
	     if ( nsel_slsig  > 0. ) { nevt_err_slsig = 1.5*xsec8TeV*sqrt(nsel_slsig) ; }
	     if ( nsel_sl     > 0. ) { nevt_err_sl    = 1.5*xsec8TeV*sqrt(nsel_sl )   ; }
	     if ( nsel_ldp    > 0. ) { nevt_err_ldp   = 1.5*xsec8TeV*sqrt(nsel_ldp)   ; }
 
             printf ( " xsec8TeV weighted events: mGl=%d, mLsp=%d: MET,HT (%d,%d) nb=%d   SIG = %6.1f +/- %4.1f,   SLSIG=%6.1f +/- %4.1f,   SL=%6.1f +/- %4.1f,   LDP=%6.1f +/- %4.1f\n",
		      mGl, mLsp, i+1, j+1, k+1,
		      1.5*xsec8TeV*(h_susy_sig[k]  -> GetBinContent( i+1, j+1 )), nevt_err_sig,
		      1.5*xsec8TeV*(h_susy_slsig[k]-> GetBinContent( i+1, j+1 )), nevt_err_slsig,
		      1.5*xsec8TeV*(h_susy_sl[k]   -> GetBinContent( i+1, j+1 )), nevt_err_sl,
		      1.5*xsec8TeV*(h_susy_ldp[k]  -> GetBinContent( i+1, j+1 )), nevt_err_ldp  ) ;
	     
	     
	  }
	}
      }




  //----------------------------------------------------------------------------

      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {
             if ( flatDummyErr < 0 ) {

                //-- compute approximate stat err.
                //-- This is 100* sig_eff / eff = 100 * [sqrt(sel)/N]/[sel/N] = 100/sqrt(sel).
                double nsel_sig   = h_susy_sig[k]  -> GetBinContent( i+1, j+1 ) ;
                double nsel_slsig = h_susy_slsig[k]-> GetBinContent( i+1, j+1 ) ;
                double nsel_sl    = h_susy_sl[k]   -> GetBinContent( i+1, j+1 ) ;
                double nsel_ldp   = h_susy_ldp[k]  -> GetBinContent( i+1, j+1 ) ;
                double frerr_sig   = 100 ;
                double frerr_slsig = 100 ;
                double frerr_sl    = 100 ;
                double frerr_ldp   = 100 ;
                if ( nsel_sig   > 0. ) { frerr_sig   = 100./sqrt(nsel_sig)    ; }
                if ( nsel_slsig > 0. ) { frerr_slsig = 100./sqrt(nsel_slsig ) ; }
                if ( nsel_sl    > 0. ) { frerr_sl    = 100./sqrt(nsel_sl )    ; }
                if ( nsel_ldp   > 0. ) { frerr_ldp   = 100./sqrt(nsel_ldp)    ; }

                inFile << frerr_sig << " " << frerr_slsig << " " << frerr_sl << " " << frerr_ldp << " " ;

                printf ( " MC sig_eff/eff (%%): mGl=%d, mLsp=%d: MET,HT (%d,%d) nb=%d   SIG = %5.1f,   SLSIG=%5.1f,   SL=%5.1f,   LDP=%5.1f\n",
			 mGl, mLsp, i+1, j+1, k+1,
			 frerr_sig, frerr_slsig, frerr_sl, frerr_ldp ) ;
             }
          }
          printf("----------------\n") ;
        }
      }
      printf("\n\n") ;


      // dump the errors to file


      // Sig uncertainty
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

	    if ( flatDummyErr >= 0 ) {
	      inFile << flatDummyErr << " " ;
	    }
	    else {
	      inFile << frerr_sig << " " ;
	    }

	  }
	}
      }


      // SLSig uncertainty
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

	    if ( flatDummyErr >= 0 ) {
	      inFile << flatDummyErr << " " ;
	    }
	    else {
	      inFile << frerr_slsig << " " ;
	    }

	  }
	}
      }


      // SL uncertainty
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

	    if ( flatDummyErr >= 0 ) {
	      inFile << flatDummyErr << " " ;
	    }
	    else {
	      inFile << frerr_sl << " " ;
	    }

	  }
	}
      }


      // Ldp uncertainty
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

	    if ( flatDummyErr >= 0 ) {
	      inFile << flatDummyErr << " " ;
	    }
	    else {
	      inFile << frerr_ldp << " " ;
	    }

	  }
	}
      }


      inFile << endl ;

      //    } // mLsp
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
