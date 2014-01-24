#include <iostream>
#include <fstream>
#include <sstream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TROOT.h"

void GenerateSusyFile( int Sample = 0, double flatDummyErr = 0.00001 ) {  //-- flat error in %.  If negative, use MC stat err.

//gluino cross sections in pb.


  TFile prospino("referenceXSecs.root");

  TH1F *gluinoxsec8TeV = (TH1F*) prospino.Get("gluino8TeV_NLONLL") ;
  if ( gluinoxsec8TeV == 0 ) {
     printf("\n\n *** Can't find histogram with name gluino8TeV_NLONLL in referenceXSecs.root.\n\n") ;
     return ;
  }

  // bin 1 = gluino mass = 100 GeV, 2 = 125, 3 = 150, 4 = 175, ...
  // so gluino mass = 75+nbin*25; or nbin = (gluinomass-75)/25.

  TChain chainTChiHH("tree");
  if ( Sample == 0 ) {
    chainTChiHH.Add("files_20fb_v71_wip/TChiHH_1.root");
    chainTChiHH.Add("files_20fb_v71_wip/TChiHH_2.root");
  }
  else if ( Sample == 1 ) {
    chainTChiHH.Add("files_20fb_v71_wip/TChiHH_WWbb.root");
  }
  else if ( Sample == 2 ) {
    chainTChiHH.Add("files_20fb_v71_wip/TChiHH_2taubb.root");
  }
  else { cout << "Unknown Sample ... exiting " << endl ; return ; }


  gROOT->Reset();

  const int nJetsCut = 3 ;     // #jets >= nJetsCut
  const int MTbCut   = 0 ;   // MTb cut


  // get binning definition from "Binning.txt"

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

  float Mbins[nBinsVar1+1] ;
  float Hbins[nBinsVar2+1] ;

  for ( int i = 0 ; i < nBinsVar1 ; i++ ) {
    inBinning >> label >> Mbins[i] ;
    if ( i > 0 && Mbins[i] < Mbins[i-1] ) { cout << "\n\n Mismatch in Var1 binning, check Binning.txt! \n\n" ; return ; }
  }

  for ( int i = 0 ; i < nBinsVar2 ; i++ ) {
    inBinning >> label >> Hbins[i] ;
    if ( i > 0 && Hbins[i] < Hbins[i-1] ) { cout << "\n\n Mismatch in Var2 binning, check Binning.txt! \n\n" ; return ; }
  }

  const TString sVar1 = in_Var1 ;
  const TString sVar2 = in_Var2 ;

  TString aVar1 = in_Var1 ;
  TString aVar2 = in_Var2 ;
  if ( aVar1 == "MET/sqrt(HT)" ) aVar1 = "MET_div_sqrtHT" ;
  if ( aVar2 == "MET/sqrt(HT)" ) aVar2 = "MET_div_sqrtHT" ;

  TString s2DVars = sVar2;
  s2DVars += ":";
  s2DVars += sVar1;

  Mbins[nBinsVar1] = 999999. ;
  Hbins[nBinsVar2] = 999999. ;

  inBinning >> label >> in_version ;
  const int version = in_version ;

  if ( !label.Contains("version") ) {
    cout << "\n\n Found inconsistency in Binning.txt, check number of bins and Var1 and Var2 lower bounds\n\n" << endl ;
    return ;
  }

  inBinning.close();



  TString cBbins[nBinsBjets];
  
  cBbins[0] = "nB==1" ;

  if ( nBinsBjets == 2 ) {
    cBbins[1] = "nB>=2" ;
  }
  else if ( nBinsBjets == 3 ) {
    cBbins[1] = "nB==2" ;    
    cBbins[2] = "nB>=3" ;
  }
  else if ( nBinsBjets == 4 ) {
    cBbins[1] = "nB==2" ;    
    cBbins[2] = "nB==3" ;
    cBbins[3] = "nB>=4" ;
  }
  else {
    cout << "\n\n This number of #b-jets bins is not implemented! Exiting ... " << endl ;
    return ;
  }

  // jet pt thresholds
  double minLeadJetPt = 70. ;
  double min3rdJetPt = 50. ;

  bool ExcludeHiggs = false ;
  //TString sLooseHiggsCuts = "";
  //TString sLooseHiggsCuts = "njets20<=7&&deltaRmax_hh<2.4&&((higgsMbb1MassDiff>95&&higgsMbb1MassDiff<145&&higgsMbb2MassDiff==-1)||(higgsMbb1MassDiff>95&&higgsMbb1MassDiff<145&&higgsMbb2MassDiff>95&&higgsMbb2MassDiff<145))&&METsig>30.&&deltaPhiStar>0.1&&pt_1st_leadJet<500&&MET>125&&HT>300&&" ;
  TString sLooseHiggsCuts = "njets20<=7&&deltaRmax_hh<2.4&&(higgsMbb1MassDiff>90&&higgsMbb1MassDiff<150&&higgsMbb2MassDiff>90&&higgsMbb2MassDiff<150)&&METsig>30.&&deltaPhiStar>0.1&&pt_1st_leadJet<500&&MET>125&&HT>300&&" ;


  // dummy masses
  int minGlMass = 150 ;
  int maxGlMass = 500 ;


//double dummyYield = 9.9 ;
  double dummyCorr = 1. ;
  long dummyEvts = 10000 ;

  ofstream inFile;
  char outfile[10000] ;

  if ( Sample == 0 ) sprintf( outfile, "datfiles/TChiHH_bbbb-v%d.dat", version ) ;
  if ( Sample == 1 ) sprintf( outfile, "datfiles/TChiHH_WWbb-v%d.dat", version ) ;
  if ( Sample == 2 ) sprintf( outfile, "datfiles/TChiHH_2taubb-v%d.dat", version ) ;

  inFile.open( outfile );

  // loop over gluino masses

  // TH1F* ht = new TH1F("ht","ht",10,0,10000);
  TString cutsSig   = "minDelPhiN>4&&nMu==0&&nEl==0&&";
  TString cutsSLSig = "MT>100&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nEl==1&&nMu==0) )&&";
  TString cutsSL    = "MT<100&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nEl==1&&nMu==0) )&&";
  TString cutsLDP   = "minDelPhiN<4&&nMu==0&&nEl==0&&";

  char commoncuts[10000] ;

  if ( ExcludeHiggs ) {
    sprintf( commoncuts, "!passHiggsSel&&maxChNMultDiff<40&&pfOcaloMET<2.0&&nJets>=%d&&MT_bestCSV>%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f&&pt_3rd_leadJet>%.0f)&&",
	     nJetsCut, MTbCut, minLeadJetPt, minLeadJetPt, min3rdJetPt ) ;
  }
  else {
    sprintf( commoncuts, "maxChNMultDiff<40&&pfOcaloMET<2.0&&nJets>=%d&&MT_bestCSV>%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f&&pt_3rd_leadJet>%.0f)&&",
	     nJetsCut, MTbCut, minLeadJetPt, minLeadJetPt, min3rdJetPt ) ;
  }


  cutsSig   += commoncuts ;
  cutsSLSig += commoncuts ;
  cutsSL    += commoncuts ;
  cutsLDP   += commoncuts ;

  cutsSig   += sLooseHiggsCuts ;
  cutsSLSig += sLooseHiggsCuts ;
  cutsSL    += sLooseHiggsCuts ;
  cutsLDP   += sLooseHiggsCuts ;


  TH2F* h_susy_sig[10] ;
  TH2F* h_susy_slsig[10] ;
  TH2F* h_susy_sl[10] ;
  TH2F* h_susy_ldp[10] ;
  for ( int bi=0; bi<nBinsBjets; bi++ ) {
     char hname[1000] ;
     sprintf( hname, "h_susy_sig_%db", bi+1 ) ;
     h_susy_sig[bi] = new TH2F( hname, hname, nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_susy_sig[bi] -> Sumw2() ;
     sprintf( hname, "h_susy_slsig_%db", bi+1 ) ;
     h_susy_slsig[bi] = new TH2F( hname, hname, nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_susy_slsig[bi] -> Sumw2() ;
     sprintf( hname, "h_susy_sl_%db", bi+1 ) ;
     h_susy_sl[bi] = new TH2F( hname, hname, nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_susy_sl[bi] -> Sumw2() ;
     sprintf( hname, "h_susy_ldp_%db", bi+1 ) ;
     h_susy_ldp[bi] = new TH2F( hname, hname, nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_susy_ldp[bi] -> Sumw2() ;
  }

  float xsec8TeV = -1. ;

  
  for ( int mGl = minGlMass ; mGl <= maxGlMass ; mGl += 25 ) {

    // keep only the 200, 300, and 500 points
    //if ( mGl != 200 && mGl != 300 && mGl != 500 ) continue ;

    mLsp = 1 ;
    xsec8TeV = 1. ;

    printf("\n\n  SUSY Xsecs:  8 TeV = %f\n\n", xsec8TeV ) ;

    ////// for ( int mLsp = 50 ; mLsp < ( mGl - 25 ) ; mLsp = mLsp + 25 ) {
    //    for ( int mLsp = 300 ; mLsp < 710 ; mLsp = mLsp + 400 ) {

    int EvtsGen(0);

    TString cutSMS = "mgluino>";
    cutSMS += mGl-1;
    cutSMS += "&&mgluino<";
    cutSMS += mGl+1;
    cutSMS += "&&mlsp>";
    cutSMS += mLsp-1;
    cutSMS += "&&mlsp<";
    cutSMS += mLsp+1;
    EvtsGen = chainTChiHH.Draw("MET",cutSMS);
    cutSMS += "&&";


    cout << "EvtsGen = " << EvtsGen << endl ;

      printf("\n\n") ;
      for (int k = 0 ; k < nBinsBjets ; k++) {

         TString cut = cBbins[k] ;

         TString allSigCuts   = cutSMS+cutsSig+cut ;
         TString allSLSigCuts = cutSMS+cutsSLSig+cut ;
         TString allSLCuts    = cutSMS+cutsSL+cut ;
         TString allLDPCuts   = cutSMS+cutsLDP+cut ;

         char hname[1000] ;

         sprintf( hname, "h_susy_sig_%db", k+1 ) ;
         chainTChiHH.Project( hname,s2DVars,allSigCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  SIG selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_sig[k]->Integral()/19399. ) ; cout << flush ;

         sprintf( hname, "h_susy_slsig_%db", k+1 ) ;
         chainTChiHH.Project( hname,s2DVars,allSLSigCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  SL Sig selection %6.1f events.\n", mGl, mLsp, k+1, h_susy_slsig[k]->Integral()/19399. ) ; cout << flush ;

         sprintf( hname, "h_susy_sl_%db", k+1 ) ;
         chainTChiHH.Project( hname,s2DVars,allSLCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  SL  selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_sl[k]->Integral()/19399. ) ; cout << flush ;

         sprintf( hname, "h_susy_ldp_%db", k+1 ) ;
         chainTChiHH.Project( hname,s2DVars,allLDPCuts);
         printf("   mGl=%d, mLsp=%d, nBjets = %d,  LDP selection %9.1f events.\n", mGl, mLsp, k+1, h_susy_ldp[k]->Integral()/19399. ) ; cout << flush ;

      } // k (nBjets)
      printf("\n\n") ;



      inFile << mGl << " " << mLsp-1 << " " << EvtsGen << " " ;
      printf(" mGl=%4d, mLsp=%4d\n", mGl, mLsp ) ; cout << flush ;


      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsVar1 ; i++) {
        for (int j = 0 ; j < nBinsVar2 ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

             printf ( " Raw MC counts: mGl=%d, mLsp=%d: Var1,Var2 (%d,%d) nb=%d   SIG = %9.0f, SLSIG =%9.0f, SL=%9.0f, LDP=%9.0f\n",
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
      for (int i = 0 ; i < nBinsVar1 ; i++) {
        for (int j = 0 ; j < nBinsVar2 ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

             inFile << xsec8TeV*(h_susy_sig[k]  -> GetBinContent( i+1, j+1 ))/19399. << " " ;

             totalSUSYyield += (h_susy_sig[k] -> GetBinContent( i+1, j+1 )*xsec8TeV/19399);
	    	     
          } // k
        } // j
      } // i
      printf("Total SUSY yield within current binning = %9.1f",totalSUSYyield);
      printf("\n\n") ;


      // SLSig yields:
      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsVar1 ; i++) {
        for (int j = 0 ; j < nBinsVar2 ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

             inFile << xsec8TeV*(h_susy_slsig[k]-> GetBinContent( i+1, j+1 )/19399.) << " " ;

	  }
	}
      }


      // SL yields:
      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsVar1 ; i++) {
        for (int j = 0 ; j < nBinsVar2 ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

             inFile << xsec8TeV*(h_susy_sl[k]   -> GetBinContent( i+1, j+1 )/19399.) << " " ;


	  }
	}
      }


      // Ldp yields:
      printf("----------------\n") ;
      for (int i = 0 ; i < nBinsVar1 ; i++) {
        for (int j = 0 ; j < nBinsVar2 ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

             inFile << xsec8TeV*(h_susy_ldp[k]  -> GetBinContent( i+1, j+1 )/19399.) << " " ;

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
 
             printf ( " xsec8TeV weighted events: mGl=%d, mLsp=%d: Var1,Var2 (%d,%d) nb=%d   SIG = %6.1f +/- %4.1f,   SLSIG=%6.1f +/- %4.1f,   SL=%6.1f +/- %4.1f,   LDP=%6.1f +/- %4.1f\n",
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
      for (int i = 0 ; i < nBinsVar1 ; i++) {
        for (int j = 0 ; j < nBinsVar2 ; j++) {
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

                printf ( " MC sig_eff/eff (%%): mGl=%d, mLsp=%d: Var1,Var2 (%d,%d) nb=%d   SIG = %5.1f,   SLSIG=%5.1f,   SL=%5.1f,   LDP=%5.1f\n",
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
      for (int i = 0 ; i < nBinsVar1 ; i++) {
        for (int j = 0 ; j < nBinsVar2 ; j++) {
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
      for (int i = 0 ; i < nBinsVar1 ; i++) {
        for (int j = 0 ; j < nBinsVar2 ; j++) {
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
      for (int i = 0 ; i < nBinsVar1 ; i++) {
        for (int j = 0 ; j < nBinsVar2 ; j++) {
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
      for (int i = 0 ; i < nBinsVar1 ; i++) {
        for (int j = 0 ; j < nBinsVar2 ; j++) {
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
  inFile << "Using Var2 bins:  " ;
  for (int j = 0 ; j <= nBinsVar2 ; j++ ) {
    inFile << Hbins[j] ;
    if ( j < nBinsVar2 ) inFile << "-" ;
  }

  inFile << "\t Using Var1 bins: " ;
  for (int i = 0 ; i <= nBinsVar1 ; i++ ) {
    inFile << Mbins[i] ;
    if ( i < nBinsVar1 ) inFile << "-" ;
  }
  inFile << endl ;


  return;

}
