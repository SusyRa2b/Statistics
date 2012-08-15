#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TRegexp.h"
#include "TSystem.h"
#include "TH2F.h"

  using std::stringstream ;
  using std::ofstream ;
  using std::endl ;

  void saveHist(const char* filename, const char* pat) ;
  TH1F* bookHist(const char* hname, const char* htitle, const char* selstring, int nbjet, int nBinsMET, int nBinsHT ) ;

// to add in: nMu, nEl, minDelPhi

void GenerateInputFile( double mgl=-1., double mlsp=-1., double target_susy_all0lep=-1. ) {

  TChain* dyTree = new TChain("treeZ") ;
  int nAdded = dyTree->Add("files5fb_MT/DY.root") ;
  if ( nAdded <= 0 ) {
     printf("\n\n\n *** No treeZ in files5fb_MT/DY.root\n\n\n") ;
     return ;
  }

  double t1bbbbWeight(0.) ;
  TChain chainT1bbbb("tree") ;
  char susycutstring[1000] ;
  sprintf( susycutstring, "&&mgluino==%.0f&&mlsp==%.0f", mgl, mlsp ) ;
  TString susycut( susycutstring ) ;
  if ( mgl>0. && mlsp>0. ) {
     nAdded = chainT1bbbb.Add("files5fb_MT/T1bbbb.root") ;
     if ( nAdded <= 0 ) {
        printf("\n\n\n *** No tree in files5fb_MT/T1bbbb.root\n\n\n") ;
        return ;
     }
     TFile f("referenceXSecs.root") ;
     TH1F* xsechist = (TH1F*) f.Get("gluino_NLONLL") ;
     if ( xsechist==0x0 ) { printf("\n\n *** can't find reference Xsec histogram in referenceXSecs.root.\n\n") ; return ; }
     int theBin = xsechist->FindBin( mgl ) ;
     if ( theBin <=0 || theBin > xsechist->GetNbinsX() ) {
        printf("\n\n *** can't find bin for mgl=%g.  Returned %d\n\n", mgl, theBin ) ;
        return ;
     }
     double xsec = xsechist->GetBinContent( theBin ) ;
     printf("\n\n T1bbbb xsec for mgl=%g is %g\n\n", mgl, xsec ) ;
     t1bbbbWeight = 0.5*xsec ;  //in pb. scan has 10k events, so nScan*0.5*xsec = events in 5fb-1
     //////  t1bbbbWeight = 0.1*xsec ;  //in pb. T1tttt scan has 50k events, so nScan*0.1*xsec = events in 5fb-1
     printf("\n\n Susy ttree cut: %s\n\n", susycutstring ) ;
  }


  TChain chainQCD("tree") ;
   //--- these have high weight
//chainQCD.Add("files5fb_MT/QCD-50to80.root");
//chainQCD.Add("files5fb_MT/QCD-80to120.root");
//chainQCD.Add("files5fb_MT/QCD-120to170.root");
//chainQCD.Add("files5fb_MT/QCD-170to300.root");
   //--- below here, these have weight less than one.
  chainQCD.Add("files5fb_MT/QCD-300to470.root");
  chainQCD.Add("files5fb_MT/QCD-470to600.root");
  chainQCD.Add("files5fb_MT/QCD-600to800.root");
  chainQCD.Add("files5fb_MT/QCD-800to1000.root");
  chainQCD.Add("files5fb_MT/QCD-1000to1400.root");
  chainQCD.Add("files5fb_MT/QCD-1400to1800.root");
  chainQCD.Add("files5fb_MT/QCD-1800.root");

  TChain chainTT("tree") ;
  chainTT.Add("files5fb_MT/TT.root") ;

  TChain chainZnn("tree") ;
  chainZnn.Add("files5fb_MT/Zinv.root") ;

  TChain chainWJets("tree") ;
  chainWJets.Add("files5fb_MT/WJets-300.root") ;

  TChain chainAll("tree");
//chainAll.Add("files5fb_MT/QCD-50to80.root");
//chainAll.Add("files5fb_MT/QCD-80to120.root");
//chainAll.Add("files5fb_MT/QCD-120to170.root");
//chainAll.Add("files5fb_MT/QCD-170to300.root");
  chainAll.Add("files5fb_MT/QCD-300to470.root");
  chainAll.Add("files5fb_MT/QCD-470to600.root");
  chainAll.Add("files5fb_MT/QCD-600to800.root");
  chainAll.Add("files5fb_MT/QCD-800to1000.root");
  chainAll.Add("files5fb_MT/QCD-1000to1400.root");
  chainAll.Add("files5fb_MT/QCD-1400to1800.root");
  chainAll.Add("files5fb_MT/QCD-1800.root");
  chainAll.Add("files5fb_MT/TT.root");
  chainAll.Add("files5fb_MT/Zinv.root");

  TChain chainTZ("tree");
  chainTZ.Add("files5fb_MT/TT.root");
  chainTZ.Add("files5fb_MT/Zinv.root");

  gROOT->Reset();

  const int nBinsBjets = 3 ;   // this must always be 3
  const int nJetsCut = 3 ;     // #jets >= nJetsCut

  double minLeadJetPt = 50. ;

  //-- met2-ht1-v1
//const int nBinsMET   = 2 ;
//const int nBinsHT    = 1 ;
//    const int version = 1;
//float Mbins[nBinsMET+1] = {150.,250.,99999.};
//float Hbins[nBinsHT+1] = {400.,99999.};

//-- met2-ht2-v1
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

////-- met3-ht3-v2
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

////-- met4-ht4-v1
//    const int nBinsMET   = 4 ;
//    const int nBinsHT    = 4 ;
//    const int version = 1;
//    float Mbins[nBinsMET+1] = {150.,200.,250.,300.,99999.};
//    float Hbins[nBinsHT+1] = {400.,500.,600.,800.,99999.};

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
//    const int nBinsHT    = 4 ;
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



  TString sMbins[nBinsMET];
  TString sHbins[nBinsHT];
  TString sBbins[3] = {"_1b","_2b","_3b"};

  TString cMbins[nBinsMET];
  TString cHbins[nBinsHT];

  for (int i = 0 ; i < nBinsMET ; i++) {
    TString base = "_M";
    stringstream sbin;
    sbin << i+1;
    base += sbin.str();
    sMbins[i] = base;
    base = "&&MET>";
    stringstream cbin;
    cbin << Mbins[i] << "&&MET<" << Mbins[i+1];
    base += cbin.str();
    cMbins[i] = base;
  }

  for (int j = 0 ; j < nBinsHT ; j++) {
    TString base = "_H";
    stringstream sbin;
    sbin << j+1;
    base += sbin.str();
    sHbins[j] = base;
    base = "&&HT>";
    stringstream cbin;
    cbin << Hbins[j] << "&&HT<" << Hbins[j+1];
    base += cbin.str();
    cHbins[j] = base;
  }

  TString leadJetPtCutString ;
  {
     leadJetPtCutString = "(pt_1st_leadJet>" ;
     stringstream number ;
     number << minLeadJetPt ;
     leadJetPtCutString += number.str() ;
     leadJetPtCutString += "&&pt_2nd_leadJet>" ;
     leadJetPtCutString += number.str() ;
     leadJetPtCutString += ")" ;
  }

//int dummyInt = 99;
//float dummyFloat = 9.999;
  float dummyZero = 0.;
  float dummyOne = 1.0;
  float dummyPoint999 = 0.999 ;
//float dummyErr = 0.1;
//float dummyErr = 0.001;
  float dummyErr = 0.0;

  float sl_frac2b_val[nBinsMET][nBinsHT];
  float sl_frac2b_err[nBinsMET][nBinsHT];
  float sl_frac3b_val[nBinsMET][nBinsHT];
  float sl_frac3b_err[nBinsMET][nBinsHT];


  ofstream inFile;
  char outfile[10000] ;
  if ( mgl > 0. && mlsp > 0. ) {
     if ( target_susy_all0lep > 0. ) {
        sprintf( outfile, "InputWT1bbbb-mgl%.0f-mlsp%.0f-%.0fevts-met%d-ht%d-v%d.dat", mgl, mlsp, target_susy_all0lep, nBinsMET, nBinsHT, version  ) ;
     } else {
        sprintf( outfile, "InputWT1bbbb-mgl%.0f-mlsp%.0f-met%d-ht%d-v%d.dat", mgl, mlsp, nBinsMET, nBinsHT, version  ) ;
     }
  } else {
     sprintf( outfile, "Input-met%d-ht%d-v%d.dat", nBinsMET, nBinsHT, version  ) ;
  }
  inFile.open( outfile );

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


   int nSel(3) ;
   char selname[3][100] = { "0lep", "1lep", "ldp" } ;

   char selcuts[3][10000] = {
        "minDelPhiN>4&&nMu==0&&nEl==0&&",
        "TypeI_MT<100&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&",
        "minDelPhiN<4&&nMu==0&&nEl==0&&" } ;


  //--- Output histograms.

   TH1F* hmctruth_susy[3][3] ;
   TH1F* hmctruth_ttwj[3][3] ;
   TH1F* hmctruth_ttbar[3][3] ;
   TH1F* hmctruth_wjets[3][3] ;
   TH1F* hmctruth_qcd[3][3] ;
   TH1F* hmctruth_znn[3][3] ;
   TH1F* hmctruth_allsm[3][3] ;
   TH1F* hmctruth_all[3][3] ;

   for ( int si=0; si<nSel; si++ ) {
      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

         char hname[1000] ;
         char htitle[1000] ;
         sprintf( htitle, "%s, %d btag", selname[si], bbi+1 ) ;

         sprintf( hname, "hmctruth_susy_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_susy[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_ttwj_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_ttwj[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_ttbar_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_ttbar[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_wjets_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_wjets[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_qcd_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_qcd[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_znn_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_znn[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_allsm_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_allsm[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_all_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_all[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;


      } // bbi.
   } // si.




     TH1F* hmctruth_qcd_lsb_pass = new TH1F("hmctruth_qcd_lsb_pass","LSB, QCD pass", 1+nBinsBjets*(nBinsHT+1), 0.5, 1+nBinsBjets*(nBinsHT+1)+0.5 )  ;
     TH1F* hmctruth_qcd_lsb_fail = new TH1F("hmctruth_qcd_lsb_fail","LSB, QCD fail", 1+nBinsBjets*(nBinsHT+1), 0.5, 1+nBinsBjets*(nBinsHT+1)+0.5 )  ;



     TH1F* hmctruth_fit_zee_1b  = bookHist("hmctruth_fit_zee_1b" , "Zee" , "Zee", 1, nBinsMET, nBinsHT ) ;
     TH1F* hmctruth_fit_zmm_1b  = bookHist("hmctruth_fit_zmm_1b" , "Zmm" , "Zmm", 1, nBinsMET, nBinsHT ) ;






  //--- histograms used in getting the observables.

    TH2F* h_tt[10] ;
    TH2F* h_wjets[10] ;
    TH2F* h_qcd[10] ;
    TH2F* h_znn[10] ;
    TH2F* h_susy[10] ;
    TH2F* h_mc[10] ;

    for ( int bi=0; bi<nBinsBjets; bi++ ) {

       char hname[100] ;

       sprintf( hname, "h_tt_%db", bi+1 ) ;
       h_tt[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
       h_tt[bi] -> Sumw2() ;

       sprintf( hname, "h_wjets_%db", bi+1 ) ;
       h_wjets[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
       h_wjets[bi] -> Sumw2() ;

       sprintf( hname, "h_qcd_%db", bi+1 ) ;
       h_qcd[bi]  = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
       h_qcd[bi] -> Sumw2() ;

       sprintf( hname, "h_znn_%db", bi+1 ) ;
       h_znn[bi]  = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
       h_znn[bi] -> Sumw2() ;

       sprintf( hname, "h_susy_%db", bi+1 ) ;
       h_susy[bi] = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
       h_susy[bi] -> Sumw2() ;

       sprintf( hname, "h_mc_%db", bi+1 ) ;
       h_mc[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
       h_mc[bi] -> Sumw2() ;

    }

  float nSusyTotal = 0;

  for ( int si=0 ; si<nSel ; si++ ) {


  printf("\n\n-----------------------------------------------------------------\n\n") ;

      for (int k = 0 ; k < nBinsBjets ; k++) {

        char allcuts[10000] ;
        char allsusycuts[10000] ;

        if ( k < (nBinsBjets-1) ) {
	  sprintf( allcuts, "%snB==%d&&nJets>=%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f)", selcuts[si], k+1, nJetsCut, minLeadJetPt, minLeadJetPt ) ;
	  sprintf( allsusycuts, "%snB==%d&&nJets>=%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f)%s", selcuts[si], k+1, nJetsCut, minLeadJetPt, minLeadJetPt, susycut.Data() ) ;
        } else {
	  sprintf( allcuts, "%snB>=%d&&nJets>=%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f)", selcuts[si], k+1, nJetsCut, minLeadJetPt, minLeadJetPt ) ;
	  sprintf( allsusycuts, "%snB>=%d&&nJets>=%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f)%s", selcuts[si], k+1, nJetsCut, minLeadJetPt, minLeadJetPt, susycut.Data() ) ;
        }


        printf("\n\n N_%s -- nbjet bin (%d): cuts=%s\n\n", selname[si], k, allcuts) ; cout << flush ;

        char hname[100] ;

        sprintf( hname, "h_tt_%db", k+1 ) ;
        chainTT.Project(hname,"HT:MET",allcuts);
        printf("    %12s %7.1f events\n", hname, h_tt[k]->Integral() ) ; cout << flush ;

        sprintf( hname, "h_wjets_%db", k+1 ) ;
        chainWJets.Project(hname,"HT:MET",allcuts);
        printf("    %12s %7.1f events\n", hname, h_wjets[k]->Integral() ) ; cout << flush ;

        sprintf( hname, "h_qcd_%db", k+1 ) ;
        chainQCD.Project(hname,"HT:MET",allcuts);
        printf("    %12s %7.1f events\n", hname, h_qcd[k]->Integral() ) ; cout << flush ;

        sprintf( hname, "h_znn_%db", k+1 ) ;
        chainZnn.Project(hname,"HT:MET",allcuts);
        printf("    %12s %7.1f events\n", hname, h_znn[k]->Integral() ) ; cout << flush ;
        if ( mgl > 0. ) {
           sprintf( hname, "h_susy_%db", k+1 ) ;
           chainT1bbbb.Project(hname,"HT:MET",allsusycuts);
           h_susy[k]->Scale( t1bbbbWeight ) ;
           printf("    %12s %7.1f events\n", hname, h_susy[k]->Integral() ) ; cout << flush ;
	   if (si==0) nSusyTotal += h_susy[k]->Integral();
        }
      } // k
      if (si==0) printf("N_Susy_Total = %7.1f events", nSusyTotal); cout << flush;
      printf("\n\n") ;


      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {


             char obsname[1000] ;
             sprintf( obsname, "N_%s_M%d_H%d_%db", selname[si], i+1, j+1, k+1 ) ;


             double ttval = h_tt[k] -> GetBinContent( i+1, j+1 ) ;
             double tterr = h_tt[k] -> GetBinError(   i+1, j+1 ) ;

             double wjetsval = h_wjets[k] -> GetBinContent( i+1, j+1 ) ;
             double wjetserr = h_wjets[k] -> GetBinError(   i+1, j+1 ) ;

             double qcdval = h_qcd[k] -> GetBinContent( i+1, j+1 ) ;
             double qcderr = h_qcd[k] -> GetBinError(   i+1, j+1 ) ;

             double znnval = h_znn[k] -> GetBinContent( i+1, j+1 ) ;
             double znnerr = h_znn[k] -> GetBinError(   i+1, j+1 ) ;

             double susyval = h_susy[k] -> GetBinContent( i+1, j+1 ) ;
             double susyerr = h_susy[k] -> GetBinError(   i+1, j+1 ) ;

             printf(" N_%s, tt     met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, ttval,tterr) ; cout << flush ;
             printf(" N_%s, wjets  met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, wjetsval,wjetserr) ; cout << flush ;
             printf(" N_%s, qcd    met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, qcdval,qcderr) ; cout << flush ;
             printf(" N_%s, znn    met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, znnval,znnerr) ; cout << flush ;
             if ( mgl>0. ) {
                if ( target_susy_all0lep > 0. ) {
                   susyval = susyval * (target_susy_all0lep/nSusyTotal);
                   susyerr = susyerr * (target_susy_all0lep/nSusyTotal);
                }
                printf(" N_%s, susy   met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, susyval,susyerr) ; cout << flush ;
             }
             printf("\n") ;

             double allval = ttval + wjetsval + qcdval + znnval + susyval ;

             //// inFile << obsname << "  \t" << (int)allval << endl;
             inFile << obsname << "  \t" << allval << endl;

             int histbin = 1 + (nBinsHT+1)*i + j + 1 ;

             hmctruth_ttwj[si][k]  -> SetBinContent( histbin, ttval + wjetsval  ) ;
             hmctruth_ttwj[si][k]  -> SetBinError(   histbin, sqrt( pow(tterr,2) + pow(wjetserr,2) )  ) ;
             hmctruth_ttbar[si][k] -> SetBinContent( histbin, ttval ) ;
             hmctruth_ttbar[si][k] -> SetBinError(   histbin, tterr ) ;
             hmctruth_wjets[si][k] -> SetBinContent( histbin, wjetsval  ) ;
             hmctruth_wjets[si][k] -> SetBinError(   histbin, wjetserr  ) ;
             hmctruth_qcd[si][k]   -> SetBinContent( histbin, qcdval ) ;
             hmctruth_qcd[si][k]   -> SetBinError(   histbin, qcderr ) ;
             hmctruth_znn[si][k]   -> SetBinContent( histbin, znnval ) ;
             hmctruth_znn[si][k]   -> SetBinError(   histbin, znnerr ) ;
             hmctruth_susy[si][k]  -> SetBinContent( histbin, susyval  ) ;
             hmctruth_susy[si][k]  -> SetBinError(   histbin, susyerr  ) ;
             hmctruth_allsm[si][k] -> SetBinContent( histbin, ttval+wjetsval+qcdval+znnval ) ;
             hmctruth_allsm[si][k] -> SetBinError(   histbin, sqrt( pow(tterr,2) + pow(wjetserr,2) + pow(qcderr,2) + pow(znnerr,2) ) ) ;
             hmctruth_all[si][k]   -> SetBinContent( histbin, ttval+wjetsval+qcdval+znnval+susyval ) ;
             hmctruth_all[si][k]   -> SetBinError(   histbin, sqrt( pow(tterr,2) + pow(wjetserr,2) + pow(qcderr,2) + pow(znnerr,2) + pow(susyerr,2) ) ) ;


	     // compute fractions of SL 2b/1b and 3b/1b

	     if ( si == 1 && k > 0 ) {

	       double ttval_1b = h_tt[0] -> GetBinContent( i+1, j+1 ) ;
	       double wjetsval_1b = h_wjets[0] -> GetBinContent( i+1, j+1 ) ;

	       double ttwjval = ttval + wjetsval ;
	       double ttwjval_1b = ttval_1b + wjetsval_1b ;

	       if ( k == 1 ) {
		 sl_frac2b_val[i][j] = ( ttwjval ) / ( ttwjval_1b ) ;
		 sl_frac2b_err[i][j] = 0.01 ;   // start with arbitrary errors, first
	       }

	       if ( k == 2 ) {
		 sl_frac3b_val[i][j] = ( ttwjval ) / ( ttwjval_1b ) ;
		 sl_frac3b_err[i][j] = 0.01 ;   // start with arbitrary errors, first
	       }

	     }

            } // k
          } // j
        } // i



          for (int k = 0 ; k < nBinsBjets ; k++) {
             h_tt[k] -> Reset() ;
             h_wjets[k] -> Reset() ;
             h_qcd[k] -> Reset() ;
             h_znn[k] -> Reset() ;
             h_susy[k] -> Reset() ;
          }


     } // si.







    printf("\n\n-----------------------------------------------------------------\n\n") ;
  
    TH1F* ht_pass = new TH1F("ht_pass","ht_pass",10,0,10000);
    TH1F* ht_fail = new TH1F("ht_fail","ht_fail",10,0,10000);
  
    ht_pass -> Sumw2() ;
    ht_fail -> Sumw2() ;

    //inFile << "Note: I've removed the b jet dependance of this result since it's always with zero b jets" << endl;
    // R_lsb  very low met sideband (50-100) ratio of mdp>4/mdp<4 (with zero b ratio)

    stringstream njcut ; njcut << nJetsCut ;
    TString cutslsb = "MET>50&&MET<100&&nMu==0&&nEl==0&&nB==0&&nJets>=";
    cutslsb += njcut.str();
    cutslsb += "&&";

    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {
  
        TString Rlsb = "R_lsb" ;
        Rlsb = Rlsb+sHbins[j]+sBbins[k] ;
        
        TString cut = "HT>";  
        cut += Hbins[j];
        cut += "&&HT<";
        cut += Hbins[j+1];
        
        TString pass = "&&minDelPhiN>4";
        TString fail = "&&minDelPhiN<4";
        TString allcutspass = cutslsb+cut+pass ;
     // if ( k==0 ) { chainAll.Project("ht_pass","HT",allcutspass); } // only do it once in bjet loop, since using nB==0.
        double npasserr(0.) ;
        float npass = 100. ;
     // npass = ht_pass->IntegralAndError(1,10,npasserr) ;
        printf(" R_lsb -- HT,MET bins (%d,%d): npass=%10.1f, cuts=%s\n", j,k,npass, allcutspass.Data()) ; cout << flush ;
        hmctruth_qcd_lsb_pass->SetBinContent( 1+k*(nBinsHT+1)+j+1, npass ) ;
        hmctruth_qcd_lsb_pass->SetBinError(   1+k*(nBinsHT+1)+j+1, npasserr ) ;
        char passbinlabel[1000] ;
        sprintf( passbinlabel, "ldp_H%d_%db_pass", j+1, k+1 ) ;
        hmctruth_qcd_lsb_pass->GetXaxis()->SetBinLabel( 1+k*(nBinsHT+1)+j+1, passbinlabel ) ;
  
  
        TString allcutsfail = cutslsb+cut+fail ;
     // if ( k==0 ) { chainAll.Project("ht_fail","HT",allcutsfail); } // only do it once in the bjet loop, since using nB==0.
        double nfailerr(0.) ;
        float nfail = 1000. ;
     // nfail = ht_fail->IntegralAndError(1,10,nfailerr) ;
        printf(" R_lsb -- HT,MET bins (%d,%d): nfail=%10.1f, cuts=%s\n", j,k,nfail, allcutsfail.Data()) ; cout << flush ;
        hmctruth_qcd_lsb_fail->SetBinContent( 1+k*(nBinsHT+1)+j+1, nfail ) ;
        hmctruth_qcd_lsb_fail->SetBinError(   1+k*(nBinsHT+1)+j+1, nfailerr ) ;
        char failbinlabel[1000] ;
        sprintf( failbinlabel, "ldp_H%d_%db_fail", j+1, k+1 ) ;
        hmctruth_qcd_lsb_fail->GetXaxis()->SetBinLabel( 1+k*(nBinsHT+1)+j+1, failbinlabel ) ;
        
        inFile << Rlsb << "      \t" << npass/nfail << endl;
        
        float error = TMath::Sqrt( (1/npass) + (1/nfail) )*(npass/nfail);
        Rlsb = Rlsb+"_err" ;
        inFile << Rlsb << "  \t" << error << endl;
  
      }  
      ht_pass->Reset() ;
      ht_fail->Reset() ;
    }
  
    
    printf("\n\n-----------------------------------------------------------------\n\n") ; cout << flush ;
  
    // Z -> ee observables 
  
    TH1F* ht = new TH1F("ht","ht",10,0,10000);
  
    ht -> Sumw2() ;

    TString cutszee = "cat==2&&minDelPhiNee>4&&nVLB>=1&&nJets>=";
    cutszee += njcut.str();
    cutszee += "&&";

    for (int i = 0 ; i < nBinsMET ; i++) {
      for (int j = 0 ; j < nBinsHT ; j++) {
  
        TString obs_Zee = "N_Zee" ;
        obs_Zee = obs_Zee+sMbins[i]+sHbins[j] ;
        
        TString cut = "HT>";  
        cut += Hbins[j];
        cut += "&&HT<";
        cut += Hbins[j+1];
        cut += "&&METee>";
        cut += Mbins[i];
        cut += "&&METee<";
        cut += Mbins[i+1];
  
        TString allcutsZ = cutszee+cut ;
        allcutsZ += "&&" ;
        allcutsZ += leadJetPtCutString ;
  
        dyTree->Project("ht","HT",allcutsZ);
        double allerr(0.) ;
        double allval = ht->IntegralAndError(1,10,allerr) ;
        printf(" N_Zee -- HT,MET bins (%d,%d): events=%7.1f +/- %6.1f, cuts=%s\n", j,i,allval,allerr,allcutsZ.Data() ) ; cout << flush ;
          ht->Reset() ;
  
        ////// inFile << obs_Zee << "  \t" << (int)allval << endl;
        inFile << obs_Zee << "  \t" << allval << endl;

        //Z->ee counts, with 1 VLb and sig selection, so so MET>250, HT>400, mindelphi>4, 2e, 0mu, nJets >= 3
  
        int histbin = 1 + (nBinsHT+1)*i + j + 1 ;
  
        hmctruth_fit_zee_1b -> SetBinContent( histbin, allval ) ;
        hmctruth_fit_zee_1b -> SetBinError(   histbin, allerr ) ;
  
      }
    }
  
    
    printf("\n\n-----------------------------------------------------------------\n\n") ; cout << flush ;
  
    // Z -> mm observables
  
    TString cutszmm = "cat==1&&minDelPhiNmm>4&&nVLB>=1&&nJets>=";
    cutszmm += njcut.str();
    cutszmm += "&&";

    for (int i = 0 ; i < nBinsMET ; i++) {
      for (int j = 0 ; j < nBinsHT ; j++) {
  
        TString obs_Zmm = "N_Zmm" ;
        obs_Zmm = obs_Zmm+sMbins[i]+sHbins[j] ;
        
        TString cut = "HT>";  
        cut += Hbins[j];
        cut += "&&HT<";
        cut += Hbins[j+1];
        cut += "&&METmm>";
        cut += Mbins[i];
        cut += "&&METmm<";
        cut += Mbins[i+1];
        
        TString allcutsZ = cutszmm+cut ;
        allcutsZ += "&&" ;
        allcutsZ += leadJetPtCutString ;
  
        dyTree->Project("ht","HT",allcutsZ);
        double allerr(0.) ;
        double allval = ht->IntegralAndError(1,10,allerr) ;
        printf(" N_Zmm -- HT,MET bins (%d,%d): events=%7.1f +/- %6.1f, cuts=%s\n", j,i,allval,allerr,allcutsZ.Data() ) ; cout << flush ;
          ht->Reset() ;
  
        ////// inFile << obs_Zmm << "  \t" << (int)allval << endl;
        inFile << obs_Zmm << "  \t" << allval << endl;

        //Z->mm counts, with 1 VLb and sig selection, so so MET>250, HT>400, mindelphi>4, 2mu, 0e, nJets >= 3
  
        int histbin = 1 + (nBinsHT+1)*i + j + 1 ;
  
        hmctruth_fit_zmm_1b -> SetBinContent( histbin, allval ) ;
        hmctruth_fit_zmm_1b -> SetBinError(   histbin, allerr ) ;
  
  
      }
    }
  
    //inFile << "Why are all three of these MC categories separate? I've just put them together (only QCD, Zinv, tt)" << endl;
    // Nttbarsingletopzjetsmc_ldp
  
  
    printf("\n\n-----------------------------------------------------------------\n\n") ; cout << flush ;



    for (int k = 0 ; k < nBinsBjets ; k++) {

       char allcuts[100000] ;
       if ( k < (nBinsBjets-1) ) {
          sprintf( allcuts, "nJets>=%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f)&&minDelPhiN<4&&nMu==0&&nEl==0&&nB==%d", nJetsCut, minLeadJetPt, minLeadJetPt, k+1 ) ;
       } else {
          sprintf( allcuts, "nJets>=%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f)&&minDelPhiN<4&&nMu==0&&nEl==0&&nB>=%d", nJetsCut, minLeadJetPt, minLeadJetPt, k+1 ) ;
       }

       printf("\n cuts : %s\n", allcuts ) ; cout << flush ;

       char hname[1000] ;
       sprintf( hname, "h_mc_%db", k+1 ) ;
       chainTZ.Project( hname, "HT:MET", allcuts ) ;
       printf(" ttbar+singletop+zjets  %12s %7.1f events\n", hname, h_mc[k]->Integral() ) ; cout << flush ;

    } // k.
    printf("\n\n") ;

    for (int i = 0 ; i < nBinsMET ; i++) {
      for (int j = 0 ; j < nBinsHT ; j++) {
        for (int k = 0 ; k < nBinsBjets ; k++) {

           char obsname[1000] ;
           sprintf( obsname, "N_ttbarsingletopzjetsmc_ldp_M%d_H%d_%db", i+1, j+1, k+1 ) ;

           double val, err ;
           val = h_mc[k] -> GetBinContent( i+1, j+1 ) ;
           err = h_mc[k] -> GetBinError(   i+1, j+1 ) ;

           printf(" %s : %7.1f +/- %7.1f\n", obsname, val, err ) ;

           inFile << obsname << "  \t" << val << endl ;

        }
        printf("\n") ;
      }
      printf("\n") ;
    }








    // NWJmc_ldp
  
    for (int i = 0 ; i < nBinsMET ; i++) {
      for (int j = 0 ; j < nBinsHT ; j++) {
        for (int k = 0 ; k < nBinsBjets ; k++) {
  
  	TString obs_WJmc_ldp = "N_WJmc_ldp" ;
  	obs_WJmc_ldp = obs_WJmc_ldp+sMbins[i]+sHbins[j]+sBbins[k] ;
  	
  	inFile << obs_WJmc_ldp << "  \t" << dummyZero << endl;
  // signal selection, but ldp, so MET>250, HT>400, >=1 b, mindelphi<4, 0L, nJets >= 3
        }
      }
    }
  
  
    // NZnnmc_ldp
  
    for (int i = 0 ; i < nBinsMET ; i++) {
      for (int j = 0 ; j < nBinsHT ; j++) {
        for (int k = 0 ; k < nBinsBjets ; k++) {
  
  	TString obs_Znnmc_ldp = "N_Znnmc_ldp" ;
  	obs_Znnmc_ldp = obs_Znnmc_ldp+sMbins[i]+sHbins[j]+sBbins[k] ;
  	
  	inFile << obs_Znnmc_ldp << "  \t" << dummyZero << endl;
  // signal selection, but ldp, so MET>250, HT>400, >=1 b, mindelphi<4, 0L, nJets >= 3
        }
      }
    }
  
  
    // various parameters needed for Z -> invis.
  
    // Z -> ee acceptance
  
    for (int i = 0 ; i < nBinsMET ; i++) {
  
      TString Zee_acc = "acc_Zee";
      Zee_acc = Zee_acc+sMbins[i];
  
      inFile << Zee_acc << "  \t" << dummyPoint999 << endl;
  
      Zee_acc = Zee_acc+"_err";
      inFile << Zee_acc << "  \t" << dummyErr << endl;
  
    }
  
  
    // Z -> mm acceptance
  
    for (int i = 0 ; i < nBinsMET ; i++) {
  
      TString Zmm_acc = "acc_Zmm";
      Zmm_acc = Zmm_acc+sMbins[i];
  
      inFile << Zmm_acc << "  \t" << dummyPoint999 << endl;
  
      Zmm_acc = Zmm_acc+"_err";
      inFile << Zmm_acc << "  \t" << dummyErr << endl;
  
    }
  
  
    // Z -> ll efficiencies
  
  // use 2011 values for now.
    inFile << "Z_ee_eff  \t" << 0.6774 << endl;
    inFile << "Z_ee_eff_err  \t" << 0.0580 << endl;
    inFile << "Z_mm_eff  \t" << 0.7217 << endl;
    inFile << "Z_mm_eff_err  \t" << 0.0506 << endl;
  
  
    // Z -> ee VL to nominal scale factors
  // would it make more sense to count events in only the loosest HT and/or MET bin and
  // use the scale factors to translate between the different HT and MET bins?

   //    for (int k = 0 ; k < nBinsBjets ; k++) {
   //    TString knn_ee = "knn_ee" ;
   //    knn_ee = knn_ee+sBbins[k] ;
   //    inFile << knn_ee << "  \t" << dummyOne << endl;
   //    knn_ee = knn_ee+"_err" ;
   //    inFile << knn_ee << "  \t" << dummyErr << endl;
   //  }
   //
   //
   //  // use 2011 values for now.
   //  inFile << "knn_1b     \t" << 0.401 << endl;
   //  inFile << "knn_1b_err \t" << 0.018 << endl;
   //  inFile << "knn_2b     \t" << 0.067 << endl;
   //  inFile << "knn_2b_err \t" << 0.009 << endl;
   //  inFile << "knn_3b     \t" << 0.009 << endl;
   //  inFile << "knn_3b_err \t" << 0.003 << endl;
   //

    // updated SF's to improve MC closure (the errors are the same as before)
    inFile << "knn_1b     \t" << 0.394  << endl;
    inFile << "knn_1b_err \t" << 0.018  << endl;
    inFile << "knn_2b     \t" << 0.0626 << endl;
    inFile << "knn_2b_err \t" << 0.009  << endl;
    inFile << "knn_3b     \t" << 0.0036 << endl;
    inFile << "knn_3b_err \t" << 0.003  << endl;

    // Z -> ll purity
  
    // use 2011 values for now.
    inFile << "Z_ee_pur  \t" << 0.911 << endl;
    inFile << "Z_ee_pur_err  \t" << 0.104 << endl;
    inFile << "Z_mm_pur  \t" << 0.866 << endl;
    inFile << "Z_mm_pur_err  \t" << 0.079 << endl;
  
  
    // scale factors:
  
    inFile << "sf_mc  \t" << dummyOne << endl ; 
    inFile << "sf_mc_err  \t" << dummyErr << endl; 
  
  
    // sf_qcd
  
    for (int i = 0 ; i < nBinsMET ; i++) {
      for (int j = 0 ; j < nBinsHT ; j++) {
        for (int k = 0 ; k < nBinsBjets ; k++) {
  
  	TString sf_qcd = "sf_qcd" ;
  	sf_qcd = sf_qcd+sMbins[i]+sHbins[j]+sBbins[k] ;
  	
  	inFile << sf_qcd << "  \t" << dummyOne << endl;	
  
  	sf_qcd = sf_qcd+"_err" ;
  	inFile << sf_qcd << "  \t" << dummyErr << endl;	
  
        }
      }
    }
  
  
    // sf_ttwj
  
    for (int i = 0 ; i < nBinsMET ; i++) {
      for (int j = 0 ; j < nBinsHT ; j++) {
        for (int k = 0 ; k < nBinsBjets ; k++) {
  
  	TString sf_ttwj = "sf_ttwj" ;
  	sf_ttwj = sf_ttwj+sMbins[i]+sHbins[j]+sBbins[k] ;
  	
  	inFile << sf_ttwj << "  \t" << dummyOne << endl;	
  
  	sf_ttwj = sf_ttwj+"_err" ;
  	inFile << sf_ttwj << "  \t" << dummyErr << endl;	
  
        }
      }
    }
  
  
    // sf_ee
  
    for (int k = 0 ; k < nBinsBjets ; k++) {
  
      TString sf_ee = "sf_ee" ;
      sf_ee = sf_ee+sBbins[k] ;
  
      inFile << sf_ee << "  \t" << dummyOne << endl;
      
      sf_ee = sf_ee+"_err" ;
      inFile << sf_ee << "  \t" << dummyErr << endl;
  
    }
  
  
    // sf_mm
  
    for (int k = 0 ; k < nBinsBjets ; k++) {
  
      TString sf_mm = "sf_mm" ;
      sf_mm = sf_mm+sBbins[k] ;
  
      inFile << sf_mm << "  \t" << dummyOne << endl;
      
      sf_mm = sf_mm+"_err" ;
      inFile << sf_mm << "  \t" << dummyErr << endl;
  
    }


    // btag eff err (Note: this was missing until Aug 3, 2012, but it's apparently not used.)
    inFile << "btageff_err" << " \t" << dummyErr << endl ;



    //--- Addding ttwj and znn LDP/ZL MC values

    //--- ttwj MC LDP/ZL

    for (int mbi = 0 ; mbi < nBinsMET ; mbi++) {
      for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
        int hbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
        for (int bbi = 0 ; bbi < nBinsBjets ; bbi++) {

            float ldpval = hmctruth_ttwj[2][bbi] -> GetBinContent( hbin ) ;
            float ldperr = hmctruth_ttwj[2][bbi] -> GetBinError(   hbin ) ;
            float zlval  = hmctruth_ttwj[0][bbi] -> GetBinContent( hbin ) ;
            float zlerr  = hmctruth_ttwj[0][bbi] -> GetBinError(   hbin ) ;

            float ldpoverzl = 0. ;
            float ldpoverzlerr = 0. ;

            if ( zlval > 0. && ldpval > 0. ) {
               ldpoverzl = ldpval / zlval ;
               ldpoverzlerr = ldpoverzl * sqrt( pow((zlerr/zlval),2) + pow((ldperr/ldpval),2) ) ;
            }

            char parname[1000] ;

            sprintf( parname, "ttwj_mc_ldpover0lep_ratio_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
            printf(" %s  :  %6.3f +/- %5.3f\n", parname, ldpoverzl, ldpoverzlerr ) ;
            inFile << parname << "  \t" << ldpoverzl << endl;

            sprintf( parname, "ttwj_mc_ldpover0lep_ratio_M%d_H%d_%db_err", mbi+1, hbi+1, bbi+1 ) ;
            inFile << parname << "  \t" << ldpoverzlerr << endl;

        } // bbi
      } // hbi
    } // mbi


    //--- Znn MC LDP/ZL
    //--- Note: the 2b and >=3b MC stats are too low.  Use 1b values for all 3.

    for (int mbi = 0 ; mbi < nBinsMET ; mbi++) {
      for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
        int hbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

        float ldpval = hmctruth_znn[2][0] -> GetBinContent( hbin ) ;
        float ldperr = hmctruth_znn[2][0] -> GetBinError(   hbin ) ;
        float zlval  = hmctruth_znn[0][0] -> GetBinContent( hbin ) ;
        float zlerr  = hmctruth_znn[0][0] -> GetBinError(   hbin ) ;

        float ldpoverzl = 0. ;
        float ldpoverzlerr = 0. ;

        if ( zlval > 0. && ldpval > 0. ) {
           ldpoverzl = ldpval / zlval ;
           ldpoverzlerr = ldpoverzl * sqrt( pow((zlerr/zlval),2) + pow((ldperr/ldpval),2) ) ;
        }

        char parname[1000] ;

        for (int bbi = 0 ; bbi < nBinsBjets ; bbi++) {
            sprintf( parname, "znn_mc_ldpover0lep_ratio_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
            printf(" %s  :  %6.3f +/- %5.3f\n", parname, ldpoverzl, ldpoverzlerr ) ;
            inFile << parname << "  \t" << ldpoverzl << endl;

            sprintf( parname, "znn_mc_ldpover0lep_ratio_M%d_H%d_%db_err", mbi+1, hbi+1, bbi+1 ) ;
            inFile << parname << "  \t" << ldpoverzlerr << endl;
        } // bbi

      } // hbi
    } // mbi


    // add here the fractions (bin by bin) of 2b/1b and 3b/1b SL events (just for the MC for now)

    for (int mbi = 0 ; mbi < nBinsMET ; mbi++) {
      for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {

	TString Sl2bstring     = "sl_frac_2b_val";
	TString Sl2bstring_err = "sl_frac_2b_err";
	TString Sl3bstring     = "sl_frac_3b_val";
	TString Sl3bstring_err = "sl_frac_3b_err";

	Sl2bstring     += sMbins[mbi]+sHbins[hbi] ;
	Sl2bstring_err += sMbins[mbi]+sHbins[hbi] ;
	Sl3bstring     += sMbins[mbi]+sHbins[hbi] ;
	Sl3bstring_err += sMbins[mbi]+sHbins[hbi] ;

	inFile << Sl2bstring     << "   \t" << sl_frac2b_val[mbi][hbi] << endl ;
	inFile << Sl2bstring_err << "   \t" << sl_frac2b_err[mbi][hbi] << endl ;
	inFile << Sl3bstring     << "   \t" << sl_frac3b_val[mbi][hbi] << endl ;
	inFile << Sl3bstring_err << "   \t" << sl_frac3b_err[mbi][hbi] << endl ;

      }
    }




    gSystem->Exec("mkdir -p rootfiles") ;
    char outHistName[1000] ;
    if ( mgl>0. && mlsp>0. ) {
       if ( target_susy_all0lep > 0 ) {
          sprintf( outHistName, "rootfiles/gi-plots-wsusy-mgl%.0f-mlsp%.0f-%.0fevts-met%d-ht%d-v%d.root", mgl, mlsp, target_susy_all0lep, nBinsMET, nBinsHT, version ) ;
       } else {
          sprintf( outHistName, "rootfiles/gi-plots-wsusy-mgl%.0f-mlsp%.0f-met%d-ht%d-v%d.root", mgl, mlsp, nBinsMET, nBinsHT, version ) ;
       }
    } else {
       sprintf( outHistName, "rootfiles/gi-plots-met%d-ht%d-v%d.root", nBinsMET, nBinsHT, version ) ;
    }
    saveHist( outHistName, "hmc*" ) ;


    inFile.close();
    return;
  
  }
  //==========================================================================================
  
  void saveHist(const char* filename, const char* pat)
  {
  
    cout << "\n\n Saving histograms matching " << pat << " in file " << filename << "\n\n" << flush ;
  
    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();
  
    TRegexp re(pat,kTRUE) ;
  
    TFile outf(filename,"RECREATE") ;
    TObject* obj ;
    while((obj=iter->Next())) {
      if (TString(obj->GetName()).Index(re)>=0) {
        obj->Write() ;
        std::cout << "." ;
      }
    }
    std::cout << std::endl ;
    outf.Close() ;
  
    delete iter ;
  }
  
  //==========================================================================================
  
  
    TH1F* bookHist(const char* hname, const char* htitle, const char* selstring, int nbjet, int nBinsMET, int nBinsHT ) {
  
       int nbins = nBinsMET*(nBinsHT+1) + 1 ;
  
       TH1F* retVal = new TH1F( hname, htitle, nbins, 0.5 + 0.1*(nbjet-2), nbins+0.5 + 0.1*(nbjet-2) ) ;
       TAxis* xaxis = retVal->GetXaxis() ;
  
       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
             int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
             char binlabel[1000] ;
             sprintf( binlabel, "%s_M%d_H%d_%db", selstring, mbi+1, hbi+1, nbjet ) ;
             xaxis->SetBinLabel( histbin, binlabel ) ;
          } // hbi.
       } // mbi.
  
       retVal->SetLabelSize(0.055,"x") ;
       xaxis->LabelsOption("v") ;
  
       return retVal ;
  
    }
  
  //==========================================================================================
 
