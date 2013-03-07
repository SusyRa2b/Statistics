#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TString.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TRegexp.h"
#include "TSystem.h"
#include "TH2F.h"
#include "TCanvas.h"
//#include "SmallTree.C"

  using std::stringstream ;
  using std::ofstream ;
  using std::endl ;
  using std::cout;
  using std::endl;
  using std::flush;

  void saveHist(const char* filename, const char* pat) ;
  TH1F* bookHist(const char* hname, const char* htitle, const char* selstring, int nbjet, int nBinsMET, int nBinsHT ) ;
////void FillHTMET(TChain *chain, TH2F *histo, int si, int k) ;
  
// to add in: nMu, nEl, minDelPhi

void GenerateInputFile( double mgl=-1., double mlsp=-1., double target_susy_all0lep=-1. ) {

  TChain* dyTree = new TChain("treeZ") ;
  int nAdded = dyTree->Add("filesMoriond_v3/DY-400.root") ;
  if ( nAdded <= 0 ) {
     printf("\n\n\n *** No treeZ in filesMoriond_v3/DY.root\n\n\n") ;
     return ;
  }

  double t2ttWeight(0.) ;
  TChain chainT2tt("tree") ;
  char susycutstring[1000] ;
  sprintf( susycutstring, "&&mgluino==%.0f&&mlsp==%.0f", mgl, mlsp ) ;
  TString susycut( susycutstring ) ;
  if ( mgl>0. && mlsp>0. ) {
     nAdded = chainT2tt.Add("filesMoriond_v3/T2tt.root") ;
     if ( nAdded <= 0 ) {
        printf("\n\n\n *** No tree in filesMoriond_v3/T2tt.root\n\n\n") ;
        return ;
     }
     TFile f("referenceXSecs.root") ;
     TH1F* xsechist = (TH1F*) f.Get("gluino8TeV_NLONLL") ;
     if ( xsechist==0x0 ) { printf("\n\n *** can't find reference Xsec histogram in referenceXSecs.root.\n\n") ; return ; }
     int theBin = xsechist->FindBin( mgl ) ;
     if ( theBin <=0 || theBin > xsechist->GetNbinsX() ) {
        printf("\n\n *** can't find bin for mgl=%g.  Returned %d\n\n", mgl, theBin ) ;
        return ;
     }
     double xsec = xsechist->GetBinContent( theBin ) ;
     printf("\n\n T2tt xsec for mgl=%g is %g\n\n", mgl, xsec ) ;
     t2ttWeight = 1.5*xsec ;  //in pb. scan has 10k events, so nScan*1.5*xsec = events in 15fb-1
     //////  t2ttWeight = 0.1*xsec ;  //in pb. T2tt scan has 50k events, so nScan*0.1*xsec = events in 5fb-1
     printf("\n\n Susy ttree cut: %s\n\n", susycutstring ) ;
  }


  TChain chainQCD("tree") ;
  chainQCD.Add("filesMoriond_v4/QCD-120to170.root");
  chainQCD.Add("filesMoriond_v4/QCD-170to300.root");
  chainQCD.Add("filesMoriond_v4/QCD-300to470.root");
  chainQCD.Add("filesMoriond_v4/QCD-470to600.root");
  chainQCD.Add("filesMoriond_v4/QCD-600to800.root");
  chainQCD.Add("filesMoriond_v4/QCD-800to1000.root");
  chainQCD.Add("filesMoriond_v4/QCD-1000to1400.root");
  chainQCD.Add("filesMoriond_v4/QCD-1400to1800.root");
  chainQCD.Add("filesMoriond_v4/QCD-1800.root");
  double kfactor_qcd = 1.8 ;
  printf("\n\n Rescaling QCD by %5.3f\n\n", kfactor_qcd ) ;

  TChain chainZnn("tree") ;
  chainZnn.Add("filesMoriond_v4/Zinv-100to200.root") ;
  chainZnn.Add("filesMoriond_v4/Zinv-200to400.root") ;
  chainZnn.Add("filesMoriond_v4/Zinv-400.root") ;

  TChain chainTT("tree") ;
  TChain chainTTPowheg("tree") ;
  TChain chainTTMCaNLO("tree") ;

  chainTT.Add("filesMoriond_v4/TT_FullLept.root") ;
  chainTT.Add("filesMoriond_v4/TT_SemiLept.root") ;
  chainTT.Add("filesMoriond_v4/TT_FullHad.root") ;

  //-------
  chainTTPowheg.Add("filesMoriond_v4/TT-powheg.root");
  //-------
  chainTTMCaNLO.Add("filesMoriond_v4/TT-MCatNLO.root");
  //-------


  double kfactor_tt = 0.90 ;
  printf("\n\n Rescaling ttbar by %5.3f\n\n", kfactor_tt ) ;

  TChain chainWJets("tree") ;
  chainWJets.Add("filesMoriond_v4/WJets-250to300.root") ;
  chainWJets.Add("filesMoriond_v4/WJets-300to400.root") ;
  chainWJets.Add("filesMoriond_v4/WJets-400.root") ;
  chainWJets.Add("filesMoriond_v4/T-s.root") ;
  chainWJets.Add("filesMoriond_v4/T-t.root") ;
  chainWJets.Add("filesMoriond_v4/T-tW.root") ;
  chainWJets.Add("filesMoriond_v4/Tbar-s.root") ;
  chainWJets.Add("filesMoriond_v4/Tbar-t.root") ;
  chainWJets.Add("filesMoriond_v4/Tbar-tW.root") ;
  double kfactor_wjets = 0.90 ;
  printf("\n\n Rescaling wjets by %5.3f\n\n", kfactor_wjets ) ;


  //-- make chains of W+jets and single top separately.

  TChain chainWJetsOnly("tree") ;
  chainWJetsOnly.Add("filesMoriond_v4/WJets-250to300.root") ;
  chainWJetsOnly.Add("filesMoriond_v4/WJets-300to400.root") ;
  chainWJetsOnly.Add("filesMoriond_v4/WJets-400.root") ;
  double kfactor_wjetsonly = 0.90 ;

  TChain chainSingletop("tree") ;
  chainSingletop.Add("filesMoriond_v4/T-s.root") ;
  chainSingletop.Add("filesMoriond_v4/T-t.root") ;
  chainSingletop.Add("filesMoriond_v4/T-tW.root") ;
  chainSingletop.Add("filesMoriond_v4/Tbar-s.root") ;
  chainSingletop.Add("filesMoriond_v4/Tbar-t.root") ;
  chainSingletop.Add("filesMoriond_v4/Tbar-tW.root") ;
  double kfactor_singletop = 0.90 ;

  TChain chainSingletop_s("tree") ;
  chainSingletop_s.Add("filesMoriond_v4/T-s.root") ;
  chainSingletop_s.Add("filesMoriond_v4/Tbar-s.root") ;
  double kfactor_singletops = 0.90 ;

  TChain chainSingletop_t("tree") ;
  chainSingletop_t.Add("filesMoriond_v4/T-t.root") ;
  chainSingletop_t.Add("filesMoriond_v4/Tbar-t.root") ;
  double kfactor_singletopt = 0.90 ;

  TChain chainSingletop_tw("tree") ;
  chainSingletop_tw.Add("filesMoriond_v4/T-tW.root") ;
  chainSingletop_tw.Add("filesMoriond_v4/Tbar-tW.root") ;
  double kfactor_singletoptw = 0.90 ;





//include Z->ll in VV contribution
  TChain chainVV("tree");
  chainVV.Add("filesMoriond_v4/WW.root"); 
  chainVV.Add("filesMoriond_v4/WZ.root");
  chainVV.Add("filesMoriond_v4/ZZ.root");
  chainVV.Add("filesMoriond_v4/DY-200to400.root");
  chainVV.Add("filesMoriond_v4/DY-400.root");





  char qcdinputfile[9][1000] = {
    "filesMoriond_v4/QCD-120to170.root"
    ,"filesMoriond_v4/QCD-170to300.root"
    ,"filesMoriond_v4/QCD-300to470.root"
    ,"filesMoriond_v4/QCD-470to600.root"
    ,"filesMoriond_v4/QCD-600to800.root"
    ,"filesMoriond_v4/QCD-800to1000.root"
    ,"filesMoriond_v4/QCD-1000to1400.root"
    ,"filesMoriond_v4/QCD-1400to1800.root"
    ,"filesMoriond_v4/QCD-1800.root"
  } ;

  char qcdsamplename[9][100] = {
    "qcd_0120_to_0170"
    ,"qcd_0170_to_0300"
    ,"qcd_0300_to_0470"
    ,"qcd_0470_to_0600"
    ,"qcd_0600_to_0800"
    ,"qcd_0800_to_1000"
    ,"qcd_1000_to_1400"
    ,"qcd_1400_to_1800"
    ,"qcd_1800_to_9999"
  } ;



  gROOT->Reset();

  const int nJetsCut = 3 ;     // #jets >= nJetsCut
  const int MTbCut = 0 ;       // cut on MTb

  double minLeadJetPt = 70. ;
  double min3rdJetPt = 50. ;


  //-- met4-ht4-v15
  const int nBinsMET   = 4 ;
  const int nBinsHT    = 4 ;
  const int nBinsBjets = 3 ;
  const int version = 15;
  float Mbins[nBinsMET+1] = {125.,150.,250.,350.,99999.};
  float Hbins[nBinsHT+1] = {400.,500.,800.,1000.,99999.};
  

  TString sMbins[nBinsMET];
  TString sHbins[nBinsHT];
  TString sBbins[4] = {"_1b","_2b","_3b","_4b"};

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

//int dummyInt = 99;
//float dummyFloat = 9.999;
  float dummyZero = 0.;
  float dummyOne = 1.0;
  ///  float dummyPoint999 = 0.999 ;
  float dummyErr = 0.1;

  float sl_frac2b_val[nBinsMET][nBinsHT];
  float sl_frac2b_err[nBinsMET][nBinsHT];
  float sl_frac3b_val[nBinsMET][nBinsHT];
  float sl_frac3b_err[nBinsMET][nBinsHT];


  ofstream inFile;
  char outfile[10000] ;
  if ( mgl > 0. && mlsp > 0. ) {
     if ( target_susy_all0lep > 0. ) {
        sprintf( outfile, "InputWT2tt-mgl%.0f-mlsp%.0f-%.0fevts-met%d-ht%d-v%d.dat", mgl, mlsp, target_susy_all0lep, nBinsMET, nBinsHT, version  ) ;
     } else {
        sprintf( outfile, "InputWT2tt-mgl%.0f-mlsp%.0f-met%d-ht%d-v%d.dat", mgl, mlsp, nBinsMET, nBinsHT, version  ) ;
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


  bool useBtagSF = true ;
  char bcut[nBinsBjets][100] ;
  char bcutSF[nBinsBjets][100] ;

  sprintf(bcut[0],"nB==1");
  sprintf(bcutSF[0],"prob1");

  // this will have to be modified once we get the probge4 in the tiny trees

  if ( nBinsBjets == 2 ) {
    sprintf(bcut[1],"nB>=2");
    sprintf(bcutSF[1],"prob2+probge3");
  }
  if ( nBinsBjets == 3 ) {
    sprintf(bcut[1],"nB==2");
    sprintf(bcut[2],"nB>=3");
    sprintf(bcutSF[1],"prob2");
    sprintf(bcutSF[2],"probge3");
  }
  if ( nBinsBjets == 4 ) {
    sprintf(bcut[1],"nB==2");
    sprintf(bcut[2],"nB==3");
    sprintf(bcut[3],"nB>=4");
    sprintf(bcutSF[1],"prob2");
    sprintf(bcutSF[2],"prob3");
    sprintf(bcutSF[3],"probge4");
  }



  char commoncuts[10000] ;
  sprintf( commoncuts, "maxChNMultDiff<40&&pfOcaloMET<2.0&&nJets>=%d&&MT_bestCSV>%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f&&pt_3rd_leadJet>%.0f)",
	   nJetsCut, MTbCut, minLeadJetPt, minLeadJetPt, min3rdJetPt ) ;

   int nSel(4) ;
   char selname[4][100] = { "0lep", "1lepSig", "1lep", "ldp" } ;

   char selcuts[4][10000] ;
   sprintf( selcuts[0], "minDelPhiN>4&&nMu==0&&nEl==0&&nIsoTrk==0" ) ; //--- 0lep
   sprintf( selcuts[1], "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&MT>100" ) ; //--- 1lepSig
   sprintf( selcuts[2], "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&MT<100" ) ; //--- 1lep
   sprintf( selcuts[3], "minDelPhiN<4&&nMu==0&&nEl==0&&nIsoTrk==0" ) ; //--- ldp



   //--- Output histograms.

   TH1F* hmctruth_susy[nSel][nBinsBjets] ;
   TH1F* hmctruth_ttwj[nSel][nBinsBjets] ;
   TH1F* hmctruth_ttwjpowheg[nSel][nBinsBjets] ;
   TH1F* hmctruth_ttwjmcanlo[nSel][nBinsBjets] ;
   TH1F* hmctruth_ttbar[nSel][nBinsBjets] ;
   TH1F* hmctruth_wjets[nSel][nBinsBjets] ;
   TH1F* hmctruth_wjetsonly[nSel][nBinsBjets] ;
   TH1F* hmctruth_singletop[nSel][nBinsBjets] ;
   TH1F* hmctruth_singletops[nSel][nBinsBjets] ;
   TH1F* hmctruth_singletopt[nSel][nBinsBjets] ;
   TH1F* hmctruth_singletoptw[nSel][nBinsBjets] ;
   TH1F* hmctruth_qcd[nSel][nBinsBjets] ;
   TH1F* hmctruth_znn[nSel][nBinsBjets] ;
   TH1F* hmctruth_vv[nSel][nBinsBjets] ;
   TH1F* hmctruth_allsm[nSel][nBinsBjets] ;
   TH1F* hmctruth_all[nSel][nBinsBjets] ;

   for ( int si=0; si<nSel; si++ ) {
      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

         char hname[1000] ;
         char htitle[1000] ;
         sprintf( htitle, "%s, %d btag", selname[si], bbi+1 ) ;

         sprintf( hname, "hmctruth_susy_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_susy[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_ttwj_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_ttwj[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_ttwjpowheg_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_ttwjpowheg[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_ttwjmcanlo_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_ttwjmcanlo[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_ttbar_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_ttbar[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_wjets_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_wjets[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_wjetsonly_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_wjetsonly[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_singletop_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_singletop[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_singletops_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_singletops[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_singletopt_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_singletopt[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_singletoptw_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_singletoptw[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_qcd_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_qcd[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_znn_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_znn[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_vv_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_vv[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_allsm_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_allsm[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_all_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_all[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;


      } // bbi.
   } // si.


   //--- histograms used in getting the observables.

   TH2F* h_tt[10] ;
   TH2F* h_ttpowheg[10] ;
   TH2F* h_ttmcanlo[10] ;
   TH2F* h_wjets[10] ;
   TH2F* h_wjetsonly[10] ;
   TH2F* h_singletop[10] ;
   TH2F* h_singletops[10] ;
   TH2F* h_singletopt[10] ;
   TH2F* h_singletoptw[10] ;
   TH2F* h_qcd[10] ;
   TH2F* h_znn[10] ;
   TH2F* h_vv[10];
   TH2F* h_susy[10] ;
   TH2F* h_mc[10] ;

   TH2F *h_znn_0lep_1b;
   h_znn_0lep_1b  = new TH2F( "h_znn_0lep_1b", "h_znn_0lep_1b", nBinsMET, Mbins, nBinsHT, Hbins ) ;
   h_znn_0lep_1b->Sumw2() ;

   for ( int bi=0; bi<nBinsBjets; bi++ ) {

     char hname[100] ;
     
     sprintf( hname, "h_tt_%db", bi+1 ) ;
     h_tt[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_tt[bi] -> Sumw2() ;
     
     sprintf( hname, "h_ttpowheg_%db", bi+1 ) ;
     h_ttpowheg[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_ttpowheg[bi] -> Sumw2() ;
     
     sprintf( hname, "h_ttmcanlo_%db", bi+1 ) ;
     h_ttmcanlo[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_ttmcanlo[bi] -> Sumw2() ;
     
     sprintf( hname, "h_wjets_%db", bi+1 ) ;
     h_wjets[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_wjets[bi] -> Sumw2() ;
     
     sprintf( hname, "h_wjetsonly_%db", bi+1 ) ;
     h_wjetsonly[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_wjetsonly[bi] -> Sumw2() ;
     
     sprintf( hname, "h_singletop_%db", bi+1 ) ;
     h_singletop[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_singletop[bi] -> Sumw2() ;
     
     sprintf( hname, "h_singletops_%db", bi+1 ) ;
     h_singletops[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_singletops[bi] -> Sumw2() ;
     
     sprintf( hname, "h_singletopt_%db", bi+1 ) ;
     h_singletopt[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_singletopt[bi] -> Sumw2() ;
     
     sprintf( hname, "h_singletoptw_%db", bi+1 ) ;
     h_singletoptw[bi]   = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_singletoptw[bi] -> Sumw2() ;

     sprintf( hname, "h_qcd_%db", bi+1 ) ;
     h_qcd[bi]  = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_qcd[bi] -> Sumw2() ;

     sprintf( hname, "h_znn_%db", bi+1 ) ;
     h_znn[bi]  = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_znn[bi] -> Sumw2() ;
     
     sprintf( hname, "h_vv_%db", bi+1 ) ;
     h_vv[bi]  = new TH2F( hname, hname , nBinsMET, Mbins, nBinsHT, Hbins ) ;
     h_vv[bi] -> Sumw2() ;
     
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

       if(useBtagSF) {
	 sprintf( allcuts,     "%s*weightPU*(%s&&%s)"    , bcutSF[k], commoncuts, selcuts[si] ) ;
	 sprintf( allsusycuts, "%s*weightPU*(%s&&%s&&%s)", bcutSF[k], commoncuts, selcuts[si], susycut.Data() ) ;
       }
       else {
	 sprintf( allcuts,     "weightPU*(%s&&%s&&%s)"    , commoncuts, selcuts[si], bcut[k] ) ;
	 sprintf( allsusycuts, "weightPU*(%s&&%s&&%s&&%s)", commoncuts, selcuts[si], bcut[k], susycut.Data() ) ;
       }


       printf("\n\n N_%s -- nbjet bin (%d): cuts=%s\n\n", selname[si], k, allcuts) ; cout << flush ;
       if (mgl > 0 ) printf("\n\n N_%s -- nbjet bin (%d): cuts=%s\n\n", selname[si], k, allsusycuts) ; cout << flush ;
       
       char hname[100] ;
       sprintf( hname, "h_tt_%db", k+1 ) ;
       chainTT.Project (hname,"HT:MET",allcuts);
       
       h_tt[k] -> Scale( kfactor_tt ) ;
       printf("    %12s %7.1f events\n", hname, h_tt[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_ttpowheg_%db", k+1 ) ;
       chainTTPowheg.Project (hname,"HT:MET",allcuts);
       h_ttpowheg[k] -> Scale( kfactor_tt ) ;
       printf("    %12s %7.1f events\n", hname, h_ttpowheg[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_ttmcanlo_%db", k+1 ) ;
       chainTTMCaNLO.Project (hname,"HT:MET",allcuts);
       h_ttmcanlo[k] -> Scale( kfactor_tt ) ;
       printf("    %12s %7.1f events\n", hname, h_ttmcanlo[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_wjets_%db", k+1 ) ;
       chainWJets.Project(hname,"HT:MET",allcuts);
       h_wjets[k] -> Scale( kfactor_wjets ) ;
       printf("    %12s %7.1f events\n", hname, h_wjets[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_wjetsonly_%db", k+1 ) ;
       chainWJetsOnly.Project(hname,"HT:MET",allcuts);
       h_wjetsonly[k] -> Scale( kfactor_wjetsonly ) ;
       printf("    %12s %7.1f events\n", hname, h_wjetsonly[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_singletop_%db", k+1 ) ;
       chainSingletop.Project(hname,"HT:MET",allcuts);
       h_singletop[k] -> Scale( kfactor_singletop ) ;
       printf("    %12s %7.1f events\n", hname, h_singletop[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_singletops_%db", k+1 ) ;
       chainSingletop_s.Project(hname,"HT:MET",allcuts);
       h_singletops[k] -> Scale( kfactor_singletops ) ;
       printf("    %12s %7.1f events\n", hname, h_singletops[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_singletopt_%db", k+1 ) ;
       chainSingletop_t.Project(hname,"HT:MET",allcuts);
       h_singletopt[k] -> Scale( kfactor_singletopt ) ;
       printf("    %12s %7.1f events\n", hname, h_singletopt[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_singletoptw_%db", k+1 ) ;
       chainSingletop_tw.Project(hname,"HT:MET",allcuts);
       h_singletoptw[k] -> Scale( kfactor_singletoptw ) ;
       printf("    %12s %7.1f events\n", hname, h_singletoptw[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_qcd_%db", k+1 ) ;
       chainQCD.Project(hname,"HT:MET",allcuts);
       h_qcd[k] -> Scale( kfactor_qcd ) ;
       printf("    %12s %7.1f events\n", hname, h_qcd[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_znn_%db", k+1 ) ;
       chainZnn.Project(hname,"HT:MET",allcuts);
       printf("    %12s %7.1f events\n", hname, h_znn[k]->Integral() ) ; cout << flush ;
       
       if ( si == 0 && k == 0 ) {
	 chainZnn.Project("h_znn_0lep_1b","HT:MET",allcuts);
       }


       sprintf( hname, "h_vv_%db", k+1 ) ;
       chainVV.Project(hname,"HT:MET",allcuts);
       printf("    %12s %7.1f events\n", hname, h_vv[k]->Integral() ) ; cout << flush ;

       if ( mgl > 0. ) {
	 sprintf( hname, "h_susy_%db", k+1 ) ;
	 chainT2tt.Project(hname,"HT:MET",allsusycuts);
	 h_susy[k]->Scale( t2ttWeight ) ;
	 printf("    %12s %7.1f events\n", hname, h_susy[k]->Integral() ) ; cout << flush ;
	 if (si==0) nSusyTotal += h_susy[k]->Integral();
       }
      } // k
      if (si==0) printf("N_Susy_Total = %7.1f events", nSusyTotal); cout << flush;
      printf("\n\n") ;


  /// //--- Rescale W+jets component as a test variation -------------------------------
  ///
  /// for ( int bi=0; bi<nBinsBjets; bi++ ) {
  ///   /// h_wjets[bi] -> Scale(1.5) ;
  ///   h_wjets[bi] -> Scale(0.5) ;
  /// } // bi
  ///
  /// //--------------------------------------------------------------------------------

      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

	    char obsname[1000] ;
	    sprintf( obsname, "N_%s_M%d_H%d_%db", selname[si], i+1, j+1, k+1 ) ;
	    
	    double ttval = h_tt[k] -> GetBinContent( i+1, j+1 ) ;
	    double tterr = h_tt[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double ttpowhegval = h_ttpowheg[k] -> GetBinContent( i+1, j+1 ) ;
	    double ttpowhegerr = h_ttpowheg[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double ttmcanloval = h_ttmcanlo[k] -> GetBinContent( i+1, j+1 ) ;
	    double ttmcanloerr = h_ttmcanlo[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double wjetsval = h_wjets[k] -> GetBinContent( i+1, j+1 ) ;
	    double wjetserr = h_wjets[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double wjetsonlyval = h_wjetsonly[k] -> GetBinContent( i+1, j+1 ) ;
	    double wjetsonlyerr = h_wjetsonly[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double singletopval = h_singletop[k] -> GetBinContent( i+1, j+1 ) ;
	    double singletoperr = h_singletop[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double singletopsval = h_singletops[k] -> GetBinContent( i+1, j+1 ) ;
	    double singletopserr = h_singletops[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double singletoptval = h_singletopt[k] -> GetBinContent( i+1, j+1 ) ;
	    double singletopterr = h_singletopt[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double singletoptwval = h_singletoptw[k] -> GetBinContent( i+1, j+1 ) ;
	    double singletoptwerr = h_singletoptw[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double qcdval = h_qcd[k] -> GetBinContent( i+1, j+1 ) ;
	    double qcderr = h_qcd[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double znnval = h_znn[k] -> GetBinContent( i+1, j+1 ) ;
	    double znnerr = h_znn[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double vvval = h_vv[k] -> GetBinContent( i+1, j+1 ) ;
	    double vverr = h_vv[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double susyval = h_susy[k] -> GetBinContent( i+1, j+1 ) ;
	    double susyerr = h_susy[k] -> GetBinError(   i+1, j+1 ) ;


	    printf(" N_%s, tt          met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, ttval,tterr) ; cout << flush ;
	    printf(" N_%s, wjets       met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, wjetsval,wjetserr) ; cout << flush ;
	    printf(" N_%s,(wjets only) met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, wjetsonlyval,wjetsonlyerr) ; cout << flush ;
	    printf(" N_%s,(singletop ) met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, singletopval,singletoperr) ; cout << flush ;
	    printf(" N_%s,((sng.t,s))  met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, singletopsval,singletopserr) ; cout << flush ;
	    printf(" N_%s,((sng.t,t))  met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, singletoptval,singletopterr) ; cout << flush ;
	    printf(" N_%s,((sng.t,tW)) met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, singletoptwval,singletoptwerr) ; cout << flush ;
	    printf(" N_%s, qcd         met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, qcdval,qcderr) ; cout << flush ;
	    printf(" N_%s, znn         met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, znnval,znnerr) ; cout << flush ;
	    printf(" N_%s, vv          met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, vvval,vverr) ; cout << flush ;
	    if ( mgl>0. ) {
	      if ( target_susy_all0lep > 0. ) {
		susyval = susyval * (target_susy_all0lep/nSusyTotal);
		susyerr = susyerr * (target_susy_all0lep/nSusyTotal);
	      }
	      printf(" N_%s, susy   met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, susyval,susyerr) ; cout << flush ;
	    }
	    printf("\n") ;

	    double allval = ttval + wjetsval + qcdval + znnval + vvval + susyval ;

	    //// inFile << obsname << "  \t" << (int)allval << endl;
	    inFile << obsname << "  \t" << allval << endl;
	    
	    int histbin = 1 + (nBinsHT+1)*i + j + 1 ;
	    
	    hmctruth_ttwj[si][k]  -> SetBinContent( histbin, ttval + wjetsval  ) ;
	    hmctruth_ttwj[si][k]  -> SetBinError(   histbin, sqrt( pow(tterr,2) + pow(wjetserr,2) )  ) ;
	    hmctruth_ttwjpowheg[si][k]  -> SetBinContent( histbin, ttpowhegval + wjetsval  ) ;
	    hmctruth_ttwjpowheg[si][k]  -> SetBinError(   histbin, sqrt( pow(ttpowhegerr,2) + pow(wjetserr,2) )  ) ;
	    hmctruth_ttwjmcanlo[si][k]  -> SetBinContent( histbin, ttmcanloval + wjetsval  ) ;
	    hmctruth_ttwjmcanlo[si][k]  -> SetBinError(   histbin, sqrt( pow(ttmcanloerr,2) + pow(wjetserr,2) )  ) ;
	    hmctruth_ttbar[si][k] -> SetBinContent( histbin, ttval ) ;
	    hmctruth_ttbar[si][k] -> SetBinError(   histbin, tterr ) ;
	    hmctruth_wjets[si][k] -> SetBinContent( histbin, wjetsval  ) ;
	    hmctruth_wjets[si][k] -> SetBinError(   histbin, wjetserr  ) ;
	    hmctruth_wjetsonly[si][k] -> SetBinContent( histbin, wjetsonlyval  ) ;
	    hmctruth_wjetsonly[si][k] -> SetBinError(   histbin, wjetsonlyerr  ) ;
	    hmctruth_singletop[si][k] -> SetBinContent( histbin, singletopval  ) ;
	    hmctruth_singletop[si][k] -> SetBinError(   histbin, singletoperr  ) ;
	    hmctruth_singletops[si][k] -> SetBinContent( histbin, singletopsval  ) ;
	    hmctruth_singletops[si][k] -> SetBinError(   histbin, singletopserr  ) ;
	    hmctruth_singletopt[si][k] -> SetBinContent( histbin, singletoptval  ) ;
	    hmctruth_singletopt[si][k] -> SetBinError(   histbin, singletopterr  ) ;
	    hmctruth_singletoptw[si][k] -> SetBinContent( histbin, singletoptwval  ) ;
	    hmctruth_singletoptw[si][k] -> SetBinError(   histbin, singletoptwerr  ) ;
	    hmctruth_qcd[si][k]   -> SetBinContent( histbin, qcdval ) ;
	    hmctruth_qcd[si][k]   -> SetBinError(   histbin, qcderr ) ;
	    hmctruth_znn[si][k]   -> SetBinContent( histbin, znnval ) ;
	    hmctruth_znn[si][k]   -> SetBinError(   histbin, znnerr ) ;
	    hmctruth_vv[si][k]    -> SetBinContent( histbin, vvval ) ;
	    hmctruth_vv[si][k]    -> SetBinError(   histbin, vverr ) ;
	    hmctruth_susy[si][k]  -> SetBinContent( histbin, susyval  ) ;
	    hmctruth_susy[si][k]  -> SetBinError(   histbin, susyerr  ) ;
	    hmctruth_allsm[si][k] -> SetBinContent( histbin, ttval+wjetsval+qcdval+znnval+vvval ) ;
	    hmctruth_allsm[si][k] -> SetBinError(   histbin, sqrt( pow(tterr,2) + pow(wjetserr,2) + pow(qcderr,2) + pow(znnerr,2) + pow(vverr,2) ) ) ;
	    hmctruth_all[si][k]   -> SetBinContent( histbin, ttval+wjetsval+qcdval+znnval+vvval+susyval ) ;
	    hmctruth_all[si][k]   -> SetBinError(   histbin, sqrt( pow(tterr,2) + pow(wjetserr,2) + pow(qcderr,2) + pow(znnerr,2) + pow(vverr,2) + pow(susyerr,2) ) ) ;


	     // compute fractions of SL 2b/1b and 3b/1b
	     // the fractions are computed using only the MT < 100 sample

	     if ( si == 2 && k > 0 ) {

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
	h_wjetsonly[k] -> Reset() ;
	h_singletop[k] -> Reset() ;
	h_singletops[k] -> Reset() ;
	h_singletopt[k] -> Reset() ;
	h_singletoptw[k] -> Reset() ;
	h_qcd[k] -> Reset() ;
	h_znn[k] -> Reset() ;
	h_vv[k] -> Reset() ;
	h_susy[k] -> Reset() ;
	
      }


   } // si.


   // for systematics evaluation, print out fraction of events in 1L from non-ttwj sources					     
   float totalttwj = 0.0;													     
   float totalsm = 0.0;													     
   for (int i = 0 ; i < nBinsMET ; i++) {											     
     for (int j = 0 ; j < nBinsHT ; j++) {											     
       for (int k = 0 ; k < nBinsBjets ; k++) {										     
	 int histbin = 1 + (nBinsHT+1)*i + j + 1 ;										     
	 float ttwjcount = hmctruth_ttwj[2][k] -> GetBinContent( histbin );
	 float smcount = hmctruth_allsm[2][k] -> GetBinContent( histbin );							     
	 cout << "For bin M" << i+1 << "H" << j+1 << "b" << k+1 << " 1L non-ttwj fraction = " << (1-(ttwjcount/smcount)) << endl;  
	 totalttwj += ttwjcount;												     
	 totalsm += smcount;
       }
     }
   }
   cout << "Total 1L non-ttwj fraction = " << (1-(totalttwj/totalsm)) << endl;						     
   
   // for systematics evaluation, print out fraction of events in LDP from non-QCD sources					     
   float totalqcd = 0.0;													     
   totalsm = 0.0;														     
   for (int i = 0 ; i < nBinsMET ; i++) {											     
     for (int j = 0 ; j < nBinsHT ; j++) {											     
       for (int k = 0 ; k < nBinsBjets ; k++) {										     
	 int histbin = 1 + (nBinsHT+1)*i + j + 1 ;										     
	 float qcdcount = hmctruth_qcd[3][k] -> GetBinContent( histbin );
	 float smcount = hmctruth_allsm[3][k] -> GetBinContent( histbin );							     
	 cout << "For bin M" << i+1 << "H" << j+1 << "b" << k+1 << " LDP non-qcd fraction = " << (1-(qcdcount/smcount)) << endl;   
	 totalqcd += qcdcount;												     
	 totalsm += smcount;
       }
     }
   }
   cout << "Total LDP non-qcd fraction = " << (1-(totalqcd/totalsm)) << endl;						     



   //--- Insert dummy R_lsb lines for backwards compatibility.
   for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
     for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
       char dummyline[1000] ;
       sprintf( dummyline, "R_lsb_H%d_%db       0.1", hbi+1, bbi+1 ) ;
       inFile << dummyline << endl ;
       sprintf( dummyline, "R_lsb_H%d_%db_err   0.01", hbi+1, bbi+1 ) ;
       inFile << dummyline << endl ;
     } // bbi.
   } // hbi.
   


   printf("\n\n-----------------------------------------------------------------\n\n") ; cout << flush ;
   
   { //--- scoping bracket for QCD chunk.
     
     //--- Fill histograms to be used in QCD analysis (done in mcclosure4.c).
     
     const int nQcdSamples(9) ;
     
     TCanvas* cqcd = new TCanvas("cqcd","QCD") ;
     
     TH2F*   h0lep[nQcdSamples][nBinsBjets] ;
     TH2F*   hldp [nQcdSamples][nBinsBjets] ;
     
     TChain* qcdch[nQcdSamples] ;
     
     char hname[1000] ;
     char htitle[1000] ;
     
     TH2F* hdummy = new TH2F("hdummy","",2, Mbins[0], 1500., 2, Hbins[0], 1500. ) ;
     
     printf("\n\n") ;
     for ( int si=0; si<nQcdSamples; si++ ) {
       
       qcdch[si] = new TChain("tree") ;
       printf(" %2d : connecting to %s\n", si, qcdinputfile[si] ) ;
       qcdch[si] -> Add( qcdinputfile[si] ) ;
       
       for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
	 
	 sprintf( hname, "h_0lep_%db_%s", bbi+1, qcdsamplename[si] ) ;
	 sprintf( htitle, "QCD 0lep yield, nb=%d, %s", bbi+1, qcdsamplename[si] ) ;
	 printf("         booking hist %s : %s\n", hname, htitle ) ;
	 h0lep[si][bbi] = new TH2F( hname, htitle, nBinsMET, Mbins, nBinsHT, Hbins ) ;
	 h0lep[si][bbi] -> Sumw2() ;
	 sprintf( hname, "h_ldp_%db_%s", bbi+1, qcdsamplename[si] ) ;
	 sprintf( htitle, "QCD  LDP yield, nb=%d, %s", bbi+1, qcdsamplename[si] ) ;
	 printf("         booking hist %s  : %s\n", hname, htitle ) ;
	 hldp [si][bbi] = new TH2F( hname, htitle, nBinsMET, Mbins, nBinsHT, Hbins ) ;
	 hldp [si][bbi] -> Sumw2() ;
	 
       } // bbi.
       
     } // si.
     printf("\n\n") ;
     
     
     //--- NOTE: Small PU weights can screw up the ZL/LDP ratio for the same reason
     //          that the sample weights can, since the number of selected events
     //          is small in each ratio.
     //          Turning PU weighting off for QCD closure.
     
     
     for ( int si=0; si<nQcdSamples; si++ ) {
       
       printf(" %2d : %s : 0lep\n", si, qcdsamplename[si] ) ; cout << flush ;
       for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
	 
	 char arg1[1000] ;
	 
	 char cuts0lep[10000] ;
	 sprintf( cuts0lep, "(%s&&%s&&%s)"    , commoncuts, selcuts[0], bcut[bbi] ) ;
	 printf("     %db, 0lep cuts : %s\n", bbi+1, cuts0lep ) ;
	 sprintf( arg1, "HT:MET>>h_0lep_%db_%s", bbi+1, qcdsamplename[si] ) ;
	 qcdch[si] -> Draw( arg1, cuts0lep ) ;
	 hdummy->Draw() ;
	 h0lep[si][bbi]->Draw("samecolz") ;
	 cqcd->Update() ; cqcd->Draw() ;


	 char cutsldp[10000] ;
	 sprintf( cutsldp, "(%s&&%s&&%s)"    , commoncuts, selcuts[3], bcut[bbi] ) ;
	 printf("     %db, ldp  cuts : %s\n", bbi+1, cutsldp  ) ;
	 sprintf( arg1, "HT:MET>>h_ldp_%db_%s", bbi+1, qcdsamplename[si] ) ;
	 qcdch[si] -> Draw( arg1, cutsldp, "colz" ) ;
	 hdummy->Draw() ;
	 hldp[si][bbi]->Draw("samecolz") ;
	 cqcd->Update() ; cqcd->Draw() ;
	 
       } // bbi.
       
     } // si.
     
   } //--- scoping bracket for QCD chunk.

   printf("\n\n-----------------------------------------------------------------\n\n") ; cout << flush ;
  

   // derive Z -> ee and Z -> mm observables from Z -> inv. 1b counts
   // using the following dummy values for the other parameters:

   float BfZnn = 5.95 ;

   float dummy_AccEE = 0.75 ;
   float dummy_AccMM = 0.80 ;

   float dummy_EffEE = 0.50 ;
   float dummy_EffMM = 0.55 ;

   float dummy_PurEE = 0.80 ;
   float dummy_PurMM = 0.80 ;

   float dummy_Knn1b = 0.40 ;


   for (int i = 0 ; i < nBinsMET ; i++) {
     for (int j = 0 ; j < nBinsHT ; j++) {
       
       TString obs_Zee = "N_Zee" ;
       obs_Zee = obs_Zee+sMbins[i]+sHbins[j] ;

       double allval = h_znn_0lep_1b->GetBinContent( i+1, j+1 );
       allval = allval * ( dummy_AccEE * dummy_EffEE )/( dummy_PurEE * BfZnn * dummy_Knn1b ) ;

       inFile << obs_Zee << "  \t" << allval << endl;
       
     }
   }


   for (int i = 0 ; i < nBinsMET ; i++) {
     for (int j = 0 ; j < nBinsHT ; j++) {
       
       TString obs_Zmm = "N_Zmm" ;
       obs_Zmm = obs_Zmm+sMbins[i]+sHbins[j] ;

       double allval = h_znn_0lep_1b->GetBinContent( i+1, j+1 );

       allval = allval * ( dummy_AccMM * dummy_EffMM )/( dummy_PurMM * BfZnn * dummy_Knn1b ) ;

       inFile << obs_Zmm << "  \t" << allval << endl;
       
     }
   }

    
   printf("\n\n-----------------------------------------------------------------\n\n") ; cout << flush ;
  
  
   //inFile << "Why are all three of these MC categories separate? I've just put them together (only QCD, Zinv, tt)" << endl;
   // Nttbarsingletopzjetsmc_ldp
   
  
   printf("\n\n-----------------------------------------------------------------\n\n") ; cout << flush ;
    
    //--- Owen : these MC inputs are no longer used.  Insert dummy values for backwards compatibility in format.

   for (int i = 0 ; i < nBinsMET ; i++) {
     for (int j = 0 ; j < nBinsHT ; j++) {
       for (int k = 0 ; k < nBinsBjets ; k++) {
	 
	 char obsname[1000] ;
	 sprintf( obsname, "N_ttbarsingletopzjetsmc_ldp_M%d_H%d_%db", i+1, j+1, k+1 ) ;
	 
	 double val, err ;
	 val = 0. ;
	 err = 0. ;
	 
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
   
    printf("\n\n-----------------------------------------------------------------\n\n") ; cout << flush ;
  
    // various parameters needed for Z -> invis.

    inFile << "acc_Zee_M1       \t" << 0.75<< endl;
    inFile << "acc_Zee_M1_err  \t" << 0.01<< endl;
    inFile << "acc_Zee_M2       \t" << 0.75<< endl;
    inFile << "acc_Zee_M2_err  \t" << 0.01<< endl;
    inFile << "acc_Zee_M3       \t" << 0.75<< endl;
    inFile << "acc_Zee_M3_err  \t" << 0.01<< endl;
    inFile << "acc_Zee_M4       \t" << 0.75<< endl;
    inFile << "acc_Zee_M4_err \t" <<  0.01<< endl;
    inFile << "acc_Zmm_M1       \t" << 0.80<< endl;
    inFile << "acc_Zmm_M1_err  \t" << 0.01<< endl;
    inFile << "acc_Zmm_M2       \t" << 0.80<< endl;
    inFile << "acc_Zmm_M2_err  \t" << 0.01<< endl;
    inFile << "acc_Zmm_M3       \t" << 0.80<< endl;
    inFile << "acc_Zmm_M3_err \t" <<  0.01<< endl;
    inFile << "acc_Zmm_M4       \t" << 0.80<< endl;
    inFile << "acc_Zmm_M4_err  \t" << 0.01<< endl;

  
    // Z -> ll efficiencies (these are eff_reco**2 * eff_sel**2 * eff_trig
    //                       or Z_ee_eff*Z_ee_eff*Z_ee_trg*Z_ee_rec*Z_ee_rec from Zinv_inputs.dat)
    
    inFile << "Z_ee_eff          \t" << 0.50 << endl; // 
    inFile << "Z_ee_eff_err      \t" << 0.05 << endl; // 
    inFile << "Z_mm_eff          \t" << 0.55 << endl; // 
    inFile << "Z_mm_eff_err      \t" << 0.05 << endl; // 


    inFile << "knn_1b_M1        \t" << 0.40 << endl;
    inFile << "knn_1b_M1_err    \t" << 0.05 << endl;
    inFile << "knn_1b_M2        \t" << 0.40 << endl;
    inFile << "knn_1b_M2_err    \t" << 0.05 << endl;
    inFile << "knn_1b_M3       \t" <<  0.40 << endl;
    inFile << "knn_1b_M3_err    \t" << 0.05 << endl;
    inFile << "knn_1b_M4        \t" << 0.40 << endl;
    inFile << "knn_1b_M4_err    \t" << 0.05 << endl;
    inFile << "knn_2b           \t" << 0.05 << endl;
    inFile << "knn_2b_err       \t" << 0.02 << endl;

    if ( nBinsBjets > 2 ) {
      inFile << "knn_3b           \t" << 0.002<< endl;
      inFile << "knn_3b_err       \t" << 0.001<< endl;
    }

    if ( nBinsBjets > 3 ) {
      inFile << "knn_4b           \t" << 0.0002<< endl;
      inFile << "knn_4b_err       \t" << 0.0001<< endl;
    }

    // Z -> ll purity
  
    inFile << "Z_ee_pur  \t" << 0.80 << endl;
    inFile << "Z_ee_pur_err  \t" << 0.10 << endl;
    inFile << "Z_mm_pur  \t" << 0.80 << endl;
    inFile << "Z_mm_pur_err  \t" << 0.10 << endl;
  
    
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
    

    // sf_ttwj_slSig
  
    for (int i = 0 ; i < nBinsMET ; i++) {
      for (int j = 0 ; j < nBinsHT ; j++) {
        for (int k = 0 ; k < nBinsBjets ; k++) {
  
	  TString sf_ttwj_slSig = "sf_ttwj_slSig" ;
	  sf_ttwj_slSig = sf_ttwj_slSig+sMbins[i]+sHbins[j]+sBbins[k] ;
	  
	  inFile << sf_ttwj_slSig << "  \t" << dummyOne << endl;	
	  
	  sf_ttwj_slSig = sf_ttwj_slSig+"_err" ;
	  inFile << sf_ttwj_slSig << "  \t" << dummyErr << endl;	
	  
        }
      }
    }
    
  
    // sf_ee
    
    inFile << "sf_ee \t" << 1.0 << endl;
    inFile << "sf_ee_err \t" << 0.12 << endl;
    
    
    // sf_mm
  
    inFile << "sf_mm \t" << 1.0 << endl;
    inFile << "sf_mm_err \t" << 0.15 << endl;

    
    // btag eff err (Note: this was missing until Aug 3, 2012, but it's apparently not used.)
    inFile << "btageff_err" << " \t" << dummyErr << endl ;



    //--- Addding ttwj and znn LDP/ZL MC values

    //--- ttwj MC LDP/ZL

    for (int mbi = 0 ; mbi < nBinsMET ; mbi++) {
      for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
        int hbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
        for (int bbi = 0 ; bbi < nBinsBjets ; bbi++) {
	  
	  float ldpval = hmctruth_ttwj[3][bbi] -> GetBinContent( hbin ) ;
	  float ldperr = hmctruth_ttwj[3][bbi] -> GetBinError(   hbin ) ;
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
	
        float ldpval = hmctruth_znn[3][0] -> GetBinContent( hbin ) ;
        float ldperr = hmctruth_znn[3][0] -> GetBinError(   hbin ) ;
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

    // print out combined mc truth values for VV and DY (included in fit from MC)
    for ( int si=0; si<nSel; si++ ) {
      for (int mbi = 0 ; mbi < nBinsMET ; mbi++) {
    	for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
    	  int hbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
    	  for (int bbi = 0 ; bbi < nBinsBjets ; bbi++) {

	    float vvvalue = hmctruth_vv[si][bbi] -> GetBinContent( hbin ) ;
	    char parname[1000] ;
	    
	    sprintf( parname, "N_VVmc_%s_M%d_H%d_%db", selname[si], mbi+1, hbi+1, bbi+1 ) ;
	    printf(" %s  :  %6.3f \n", parname, vvvalue ) ;
	    inFile << parname << "  \t" << vvvalue << endl;
	    

    	  } // bbi
    	} // hbi
      } // mbi
    } // si

    // next write out measured trigger efficiencies. Values from plotEMuFrac.C from averaging e/mu bin by bin.
    float trigeff1LVal[nBinsHT][nBinsMET] = {{0.90418,0.983407,0.999797,0.999801},{0.952519,0.996088,0.999794,0.999807},{1,1,1,1},{1,1,1,1}};
    float trigeff1LErr[nBinsHT][nBinsMET] = {{0.0254576,0.0115725,0.0112111,0.0110318},{0.0153225,0.0101077,0.0101169,0.0102007},{0.0151531,0.0110547,0.013888,0.014495},{0.0108621,0.0108621,0.0108621,0.0108621}};
    float trigeff0LVal[nBinsHT][nBinsMET] = {{0.8,0.833333,1,1},{0.666667,1,1,1},{1,1,1,1},{1,1,1,1}};
    float trigeff0LErr[nBinsHT][nBinsMET] = {{0.136256,0.169997,0.0560853,0.0560853},{0.124213,0.0737418,0.0560853,0.0560853},{0.0207913,0.0207913,0.0207913,0.0207913},{0.0207913,0.0207913,0.0207913,0.0207913}};

    for (int mbi = 0 ; mbi < nBinsMET ; mbi++) {
      for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
    	char parname[1000] ;
    	sprintf( parname, "trigeff_val_0L_M%d_H%d", mbi+1, hbi+1 ) ;
    	printf(" %s  :  %6.3f \n", parname, trigeff0LVal[hbi][mbi] ) ;
    	inFile << parname << "  \t" << trigeff0LVal[hbi][mbi] << endl;
	sprintf( parname, "trigeff_err_0L_M%d_H%d", mbi+1, hbi+1 ) ;
    	printf(" %s  :  %6.3f \n", parname, trigeff0LErr[hbi][mbi] ) ;
    	inFile << parname << "  \t" << trigeff0LErr[hbi][mbi] << endl;
      }
    }
    for (int mbi = 0 ; mbi < nBinsMET ; mbi++) {
      for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
    	char parname[1000] ;
    	sprintf( parname, "trigeff_val_1L_M%d_H%d", mbi+1, hbi+1 ) ;
    	printf(" %s  :  %6.3f \n", parname, trigeff1LVal[hbi][mbi] ) ;
    	inFile << parname << "  \t" << trigeff1LVal[hbi][mbi] << endl;
    	sprintf( parname, "trigeff_err_1L_M%d_H%d", mbi+1, hbi+1 ) ;
    	printf(" %s  :  %6.3f \n", parname, trigeff1LErr[hbi][mbi] ) ;
    	inFile << parname << "  \t" << trigeff1LErr[hbi][mbi] << endl;
      }
    }






   //-- Oct 31, 2012 : add global uncertainties here (e.g. lumi, met cleaning, ...)

    inFile << "GU_luminosity   0.044" << endl ;
    inFile << "GU_metcleaning  0.031" << endl ;
    inFile << "GU_JER          0.020" << endl ;
    inFile << "GU_unclMET      0.010" << endl ;


   //-- Nov 14, 2012 : add QCD model parameters to avoid hardwiring things in ra2bRoostatsClass3D_3b
   //    see logfile of mcclosure4 for values.
   //    Uncertainties are (SFmet3-SFmet2)/2
   //                      (SFmet4-SFmet2)/2
   //                      (1-SFnb3)/2 added in quad with stat err from chi2 fit.

    inFile << "SFqcd_met3       1.49" << endl ;
    inFile << "SFqcd_met3_err   0.15" << endl ;
    inFile << "SFqcd_met4       2.14" << endl ;
    inFile << "SFqcd_met4_err   0.48" << endl ;
    inFile << "SFqcd_nb3        0.79" << endl ;
    inFile << "SFqcd_nb3_err    0.18" << endl ;

    inFile << "SFqcd_nb4        1.00" << endl ;     // these are dummy values for now!!!
    inFile << "SFqcd_nb4_err    0.10" << endl ;




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
    saveHist( outHistName, "h*" ) ;


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
  
//  //==========================================================================================
//
//  void FillHTMET(TChain *chain, TH2F *histo, int si, int k) {
//
//     TObjArray *fileElements=chain->GetListOfFiles();
//     TIter next(fileElements);
//     TChainElement *chEl=0;
//     while (( chEl=(TChainElement*)next() )) {
//        TFile f(chEl->GetTitle());
//        TTree *tree = (TTree*)f.Get("tree");
//	SmallTree *t = new SmallTree(tree);
//	t->Loop(histo, si, k);
//     }
//
//     SmallTree *t = new SmallTree(chain);
//     t->Loop(histo, si, k);
//
//    return;    
//  }
//
//  //==========================================================================================
