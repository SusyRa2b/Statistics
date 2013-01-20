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
///////#include "SmallTree.C"

  using std::stringstream ;
  using std::ofstream ;
  using std::endl ;
  using std::cout;
  using std::endl;
  using std::flush;

  void saveHist(const char* filename, const char* pat) ;
  TH1F* bookHist(const char* hname, const char* htitle, const char* selstring, int nbjet, int nBinsMET, int nBinsHT ) ;
  /////// void FillHTMET(TChain *chain, TH2F *histo, int si, int k) ;
  
// to add in: nMu, nEl, minDelPhi

void GenerateInputFile( double mgl=-1., double mlsp=-1., double target_susy_all0lep=-1. ) {

  TChain* dyTree = new TChain("treeZ") ;
  int nAdded = dyTree->Add("files15fb_8TeV_old1/DY.root") ;
  if ( nAdded <= 0 ) {
     printf("\n\n\n *** No treeZ in files15fb_8TeV_old1/DY.root\n\n\n") ;
     return ;
  }

  double t1bbbbWeight(0.) ;
  TChain chainT1bbbb("tree") ;
  char susycutstring[1000] ;
  sprintf( susycutstring, "mgluino==%.0f&&mlsp==%.0f", mgl, mlsp ) ;
  TString susycut( susycutstring ) ;
  if ( mgl>0. && mlsp>0. ) {
     nAdded = chainT1bbbb.Add("filesHCP_53_v6/T1bbbb.root") ;
     if ( nAdded <= 0 ) {
        printf("\n\n\n *** No tree in filesHCP_53_v6/T1bbbb.root\n\n\n") ;
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
     printf("\n\n T1bbbb xsec for mgl=%g is %g\n\n", mgl, xsec ) ;
     t1bbbbWeight = 1.5*xsec ;  //in pb. scan has 10k events, so nScan*1.5*xsec = events in 15fb-1
     //////  t1bbbbWeight = 0.1*xsec ;  //in pb. T1tttt scan has 50k events, so nScan*0.1*xsec = events in 5fb-1
     printf("\n\n Susy ttree cut: %s\n\n", susycutstring ) ;
  }


  TChain chainQCD("tree") ;
  chainQCD.Add("filesHCP_53_v6/QCD-120to170.root");
  chainQCD.Add("filesHCP_53_v6/QCD-170to300.root");
  chainQCD.Add("filesHCP_53_v6/QCD-300to470.root");
  chainQCD.Add("filesHCP_53_v6/QCD-470to600.root");
  chainQCD.Add("filesHCP_53_v6/QCD-600to800.root");
  chainQCD.Add("filesHCP_53_v6/QCD-800to1000.root");
  chainQCD.Add("filesHCP_53_v6/QCD-1000to1400.root");
  chainQCD.Add("filesHCP_53_v6/QCD-1400to1800.root");
  chainQCD.Add("filesHCP_53_v6/QCD-1800.root");
  double kfactor_qcd = 1.8 ;
  printf("\n\n Rescaling QCD by %5.3f\n\n", kfactor_qcd ) ;

  TChain chainZnn("tree") ;
  chainZnn.Add("filesHCP_53_v6/Zinv-100to200.root") ;
  chainZnn.Add("filesHCP_53_v6/Zinv-200to400.root") ;
  chainZnn.Add("filesHCP_53_v6/Zinv-400.root") ;

  TChain chainTT("tree") ;
  //-------
  //chainTT.Add("filesHCP_53_v6/TT.root") ;
  //-------
  //chainTT.Add("filesHCP_53_v6/TT-powheg.root");
  //-------
    chainTT.Add("filesHCP_53_v6/TT_FullLept.root") ;
    chainTT.Add("filesHCP_53_v6/TT_SemiLept.root") ;
    chainTT.Add("filesHCP_53_v6/TT_FullHad.root") ;
  //-------
  double kfactor_tt = 0.90 ;
  printf("\n\n Rescaling ttbar by %5.3f\n\n", kfactor_tt ) ;

  TChain chainWJets("tree") ;
  chainWJets.Add("filesHCP_53_v6/WJets-250to300.root") ;
  chainWJets.Add("filesHCP_53_v6/WJets-300to400.root") ;
  chainWJets.Add("filesHCP_53_v6/WJets-400.root") ;
  chainWJets.Add("filesHCP_53_v6/T-s.root") ;
  chainWJets.Add("filesHCP_53_v6/T-t.root") ;
  chainWJets.Add("filesHCP_53_v6/T-tW.root") ;
  chainWJets.Add("filesHCP_53_v6/Tbar-s.root") ;
  chainWJets.Add("filesHCP_53_v6/Tbar-t.root") ;
  chainWJets.Add("filesHCP_53_v6/Tbar-tW.root") ;
  double kfactor_wjets = 0.90 ;
  printf("\n\n Rescaling wjets by %5.3f\n\n", kfactor_wjets ) ;

//include Z->ll in VV contribution
  TChain chainVV("tree");
  chainVV.Add("filesHCP_53_v6/WW.root"); 
  chainVV.Add("filesHCP_53_v6/WZ.root");
  chainVV.Add("filesHCP_53_v6/ZZ.root");
  chainVV.Add("filesHCP_53_v6/DY-200to400.root");
  chainVV.Add("filesHCP_53_v6/DY-400.root");





      char qcdinputfile[9][1000] = {
        "filesHCP_53_v6/QCD-120to170.root"
       ,"filesHCP_53_v6/QCD-170to300.root"
       ,"filesHCP_53_v6/QCD-300to470.root"
       ,"filesHCP_53_v6/QCD-470to600.root"
       ,"filesHCP_53_v6/QCD-600to800.root"
       ,"filesHCP_53_v6/QCD-800to1000.root"
       ,"filesHCP_53_v6/QCD-1000to1400.root"
       ,"filesHCP_53_v6/QCD-1400to1800.root"
       ,"filesHCP_53_v6/QCD-1800.root"
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

  const int nBinsBjets = 3 ;   // this must always be 3
  const int nJetsCut = 3 ;     // #jets >= nJetsCut
  //////// const int nJetsCut = 5 ;     // #jets >= nJetsCut

  double minLeadJetPt = 70. ;
  double min3rdJetPt = 50. ;
  

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

////-- met3-ht3-v1
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

   bool useBtagSF = false;
   char bcut[3][100] = { "nB==1", "nB==2", "nB>=3" } ;
   char bcutSF[3][100] = { "prob1", "prob2", "probge3" } ;

   char commoncuts[10000] ;
   sprintf( commoncuts, "maxChNMultDiff<40&&pfOcaloMET<2.0&&nJets>=%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f&&pt_3rd_leadJet>%.0f)",
             nJetsCut, minLeadJetPt, minLeadJetPt, min3rdJetPt ) ;

   int nSel(3) ;
   char selname[3][100] = { "0lep", "1lep", "ldp" } ;

   char selcuts[3][10000] ;
   sprintf( selcuts[0], "minDelPhiN>4&&nMu==0&&nEl==0&&nIsoTrk==0" ) ; //--- 0lep
   sprintf( selcuts[1], "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&MT<100" ) ; //--- 1lep
   sprintf( selcuts[2], "minDelPhiN<4&&nMu==0&&nEl==0&&nIsoTrk==0" ) ; //--- ldp


  //--- Output histograms.

   TH1F* hmctruth_susy[3][3] ;
   TH1F* hmctruth_ttwj[3][3] ;
   TH1F* hmctruth_ttbar[3][3] ;
   TH1F* hmctruth_wjets[3][3] ;
   TH1F* hmctruth_qcd[3][3] ;
   TH1F* hmctruth_znn[3][3] ;
   TH1F* hmctruth_vv[3][3] ;
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

         sprintf( hname, "hmctruth_vv_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_vv[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_allsm_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_allsm[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;

         sprintf( hname, "hmctruth_all_%s_%db", selname[si], bbi+1 ) ;
         hmctruth_all[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsMET, nBinsHT ) ;


      } // bbi.
   } // si.


     bookHist("hmctruth_fit_zee_1b" , "Zee" , "Zee", 1, nBinsMET, nBinsHT ) ;
     bookHist("hmctruth_fit_zmm_1b" , "Zmm" , "Zmm", 1, nBinsMET, nBinsHT ) ;






  //--- histograms used in getting the observables.

    TH2F* h_tt[10] ;
    TH2F* h_wjets[10] ;
    TH2F* h_qcd[10] ;
    TH2F* h_znn[10] ;
    TH2F* h_vv[10];
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

        sprintf( hname, "h_wjets_%db", k+1 ) ;
        chainWJets.Project(hname,"HT:MET",allcuts);
        h_wjets[k] -> Scale( kfactor_wjets ) ;
        printf("    %12s %7.1f events\n", hname, h_wjets[k]->Integral() ) ; cout << flush ;

        sprintf( hname, "h_qcd_%db", k+1 ) ;
        chainQCD.Project(hname,"HT:MET",allcuts);
        h_qcd[k] -> Scale( kfactor_qcd ) ;
        printf("    %12s %7.1f events\n", hname, h_qcd[k]->Integral() ) ; cout << flush ;

        sprintf( hname, "h_znn_%db", k+1 ) ;
        chainZnn.Project(hname,"HT:MET",allcuts);
        printf("    %12s %7.1f events\n", hname, h_znn[k]->Integral() ) ; cout << flush ;

        sprintf( hname, "h_vv_%db", k+1 ) ;
        chainVV.Project(hname,"HT:MET",allcuts);
        printf("    %12s %7.1f events\n", hname, h_vv[k]->Integral() ) ; cout << flush ;

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

             double wjetsval = h_wjets[k] -> GetBinContent( i+1, j+1 ) ;
             double wjetserr = h_wjets[k] -> GetBinError(   i+1, j+1 ) ;

             double qcdval = h_qcd[k] -> GetBinContent( i+1, j+1 ) ;
             double qcderr = h_qcd[k] -> GetBinError(   i+1, j+1 ) ;

             double znnval = h_znn[k] -> GetBinContent( i+1, j+1 ) ;
             double znnerr = h_znn[k] -> GetBinError(   i+1, j+1 ) ;

             double vvval = h_vv[k] -> GetBinContent( i+1, j+1 ) ;
             double vverr = h_vv[k] -> GetBinError(   i+1, j+1 ) ;

             double susyval = h_susy[k] -> GetBinContent( i+1, j+1 ) ;
             double susyerr = h_susy[k] -> GetBinError(   i+1, j+1 ) ;

             printf(" N_%s, tt     met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, ttval,tterr) ; cout << flush ;
             printf(" N_%s, wjets  met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, wjetsval,wjetserr) ; cout << flush ;
             printf(" N_%s, qcd    met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, qcdval,qcderr) ; cout << flush ;
             printf(" N_%s, znn    met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, znnval,znnerr) ; cout << flush ;
             printf(" N_%s, vv     met,ht,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, vvval,vverr) ; cout << flush ;
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
             hmctruth_ttbar[si][k] -> SetBinContent( histbin, ttval ) ;
             hmctruth_ttbar[si][k] -> SetBinError(   histbin, tterr ) ;
             hmctruth_wjets[si][k] -> SetBinContent( histbin, wjetsval  ) ;
             hmctruth_wjets[si][k] -> SetBinError(   histbin, wjetserr  ) ;
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
           float ttwjcount = hmctruth_ttwj[1][k] -> GetBinContent( histbin );
           float smcount = hmctruth_allsm[1][k] -> GetBinContent( histbin );
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
           float qcdcount = hmctruth_qcd[2][k] -> GetBinContent( histbin );
           float smcount = hmctruth_allsm[2][k] -> GetBinContent( histbin );
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
               //--- Don't do this.
               //// if(useBtagSF) {
                 //// sprintf( cuts0lep, "%s*(%s&&%s)"    , bcutSF[bbi], commoncuts, selcuts[0] ) ;
               //// }
               //// else {
                 sprintf( cuts0lep, "(%s&&%s&&%s)"    , commoncuts, selcuts[0], bcut[bbi] ) ;
	       //// }
               printf("     %db, 0lep cuts : %s\n", bbi+1, cuts0lep ) ;
               sprintf( arg1, "HT:MET>>h_0lep_%db_%s", bbi+1, qcdsamplename[si] ) ;
               qcdch[si] -> Draw( arg1, cuts0lep ) ;
               hdummy->Draw() ;
               h0lep[si][bbi]->Draw("samecolz") ;
               cqcd->Update() ; cqcd->Draw() ;


	       char cutsldp[10000] ;
               //---- Don't do this.
               //// if(useBtagSF) {
                 //// sprintf( cutsldp, "%s*(%s&&%s)"    , bcutSF[bbi], commoncuts, selcuts[2] ) ;
               //// }
               else {
                 sprintf( cutsldp, "(%s&&%s&&%s)"    , commoncuts, selcuts[2], bcut[bbi] ) ;
               }
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
  
    // Z -> ee observables 
  
    TH1F* ht = new TH1F("ht","ht",10,0,10000);
  
    ht -> Sumw2() ;

    TString cutszee = "cat==2&&minDelPhiNee>4&&nVLB>=1&&nJets>=";
    stringstream njcut ; njcut << nJetsCut ;
    cutszee += njcut.str();
    cutszee += "&&";
/*
  TString leadJetPtCutString ;
  {
     leadJetPtCutString = "(pt_1st_leadJet>" ;
     stringstream number ;
     number << minLeadJetPt ;
     stringstream number2 ;
     number2 << min3rdJetPt ;
     leadJetPtCutString += number.str() ;
     leadJetPtCutString += "&&pt_2nd_leadJet>" ;
     leadJetPtCutString += number.str() ;
     leadJetPtCutString += "&&pt_3rd_leadJet>" ;
     leadJetPtCutString += number2.str() ;
     leadJetPtCutString += ")" ;
  }

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
  */
    
    printf("\n\n-----------------------------------------------------------------\n\n") ; cout << flush ;
  
    // Z -> mm observables
 /* 
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
  */

//these are the values from data
  
inFile << "N_Zee_M1_H1      \t" << 44<< endl;
inFile << "N_Zee_M1_H2      \t" << 19<< endl;
inFile << "N_Zee_M1_H3      \t" << 5<< endl;
inFile << "N_Zee_M1_H4      \t" << 3<< endl;
inFile << "N_Zee_M2_H1      \t" << 92<< endl;
inFile << "N_Zee_M2_H2      \t" << 80<< endl;
inFile << "N_Zee_M2_H3      \t" << 12<< endl;
inFile << "N_Zee_M2_H4      \t" << 2<< endl;
inFile << "N_Zee_M3_H1      \t" << 40<< endl;
inFile << "N_Zee_M3_H2      \t" << 43<< endl;
inFile << "N_Zee_M3_H3      \t" << 7<< endl;
inFile << "N_Zee_M3_H4      \t" << 5<< endl;
inFile << "N_Zee_M4_H1      \t" << 16<< endl;
inFile << "N_Zee_M4_H2      \t" << 28<< endl;
inFile << "N_Zee_M4_H3      \t" << 4<< endl;
inFile << "N_Zee_M4_H4      \t" << 3<< endl;
inFile << "N_Zmm_M1_H1      \t" << 55<< endl;
inFile << "N_Zmm_M1_H2      \t" << 47<< endl;
inFile << "N_Zmm_M1_H3      \t" << 1<< endl;
inFile << "N_Zmm_M1_H4      \t" << 0<< endl;
inFile << "N_Zmm_M2_H1      \t" << 134<< endl;
inFile << "N_Zmm_M2_H2      \t" << 124<< endl;
inFile << "N_Zmm_M2_H3      \t" << 12<< endl;
inFile << "N_Zmm_M2_H4      \t" << 3<< endl;
inFile << "N_Zmm_M3_H1      \t" << 49<< endl;
inFile << "N_Zmm_M3_H2      \t" << 51<< endl;
inFile << "N_Zmm_M3_H3      \t" << 8<< endl;
inFile << "N_Zmm_M3_H4      \t" << 3<< endl;
inFile << "N_Zmm_M4_H1      \t" << 15<< endl;
inFile << "N_Zmm_M4_H2      \t" << 37<< endl;
inFile << "N_Zmm_M4_H3      \t" << 7<< endl;
inFile << "N_Zmm_M4_H4      \t" << 4<< endl;

  
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
  
    // Z -> ee acceptance
/*  
    for (int i = 0 ; i < nBinsMET ; i++) {
  
      TString Zee_acc = "acc_Zee";
      Zee_acc = Zee_acc+sMbins[i];
  
      inFile << Zee_acc << "  \t" << dummyPoint999 << endl;
  
      Zee_acc = Zee_acc+"_err";
      inFile << Zee_acc << "  \t" << dummyErr << endl;
  
    }
*/
  
    // Z -> mm acceptance
  
  /*
    for (int i = 0 ; i < nBinsMET ; i++) {
  
      TString Zmm_acc = "acc_Zmm";
      Zmm_acc = Zmm_acc+sMbins[i];
  
      inFile << Zmm_acc << "  \t" << dummyPoint999 << endl;
  
      Zmm_acc = Zmm_acc+"_err";
      inFile << Zmm_acc << "  \t" << dummyErr << endl;
  
    }
*/
      



inFile << "acc_Zee_M1       \t" << 0.701<< endl;
inFile << "acc_Zee_M1_err   \t" << 0.014<< endl;
inFile << "acc_Zee_M2       \t" << 0.762<< endl;
inFile << "acc_Zee_M2_err   \t" << 0.016<< endl;
inFile << "acc_Zee_M3       \t" << 0.833<< endl;
inFile << "acc_Zee_M3_err   \t" << 0.017<< endl;
inFile << "acc_Zee_M4       \t" << 0.873<< endl;
inFile << "acc_Zee_M4_err   \t" << 0.018<< endl;
inFile << "acc_Zmm_M1       \t" << 0.674<< endl;
inFile << "acc_Zmm_M1_err   \t" << 0.014<< endl;
inFile << "acc_Zmm_M2       \t" << 0.733<< endl;
inFile << "acc_Zmm_M2_err   \t" << 0.014<< endl;
inFile << "acc_Zmm_M3       \t" << 0.810<< endl;
inFile << "acc_Zmm_M3_err   \t" << 0.015<< endl;
inFile << "acc_Zmm_M4       \t" << 0.857<< endl;
inFile << "acc_Zmm_M4_err   \t" << 0.018<< endl;

  
    // Z -> ll efficiencies (these are eff_reco**2 * eff_sel**2 * eff_trig
    //                       or Z_ee_eff*Z_ee_eff*Z_ee_trg*Z_ee_rec*Z_ee_rec from Zinv_inputs.dat)

///inFile << "Z_ee_eff          \t" << 0.516 << endl; // 
///inFile << "Z_ee_eff_err      \t" << 0.056 << endl; // 
///inFile << "Z_mm_eff          \t" << 0.557 << endl; // 
///inFile << "Z_mm_eff_err      \t" << 0.036 << endl; // 
inFile << "Z_ee_eff          \t" << 0.500 << endl; // 
inFile << "Z_ee_eff_err      \t" << 0.069 << endl; // 
inFile << "Z_mm_eff          \t" << 0.542 << endl; // 
inFile << "Z_mm_eff_err      \t" << 0.045 << endl; // 


  
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
  /*
    for (int i = 0 ; i < nBinsMET ; i++) {
  
      TString knn_1bS = "knn_1b";
      knn_1bS = knn_1bS+sMbins[i];
  
      inFile << knn_1bS << "  \t" << 0.394 << endl;
  
      knn_1bS = knn_1bS+"_err";
      inFile << knn_1bS << "  \t" << 0.018 << endl;
  
    }

    inFile << "knn_2b     \t" << 0.0626 << endl;
    inFile << "knn_2b_err \t" << 0.009  << endl;
    inFile << "knn_3b     \t" << 0.0036 << endl;
    inFile << "knn_3b_err \t" << 0.003  << endl;
*/

inFile << "knn_1b_M1        \t" << 0.454<< endl;
inFile << "knn_1b_M1_err    \t" << 0.067<< endl;
inFile << "knn_1b_M2        \t" << 0.438<< endl;
inFile << "knn_1b_M2_err    \t" << 0.064<< endl;
inFile << "knn_1b_M3       \t" <<  0.404<< endl;
inFile << "knn_1b_M3_err    \t" << 0.060<< endl;
inFile << "knn_1b_M4        \t" << 0.349<< endl;
inFile << "knn_1b_M4_err    \t" << 0.054<< endl;
inFile << "knn_2b           \t" << 0.112<< endl;
inFile << "knn_2b_err       \t" << 0.027<< endl;
inFile << "knn_3b           \t" << 0.0052<< endl;
inFile << "knn_3b_err       \t" << 0.0027<< endl;



    // Z -> ll purity
  
    inFile << "Z_ee_pur  \t" << 0.848 << endl;
    inFile << "Z_ee_pur_err  \t" << 0.109 << endl;
    inFile << "Z_mm_pur  \t" << 0.802 << endl;
    inFile << "Z_mm_pur_err  \t" << 0.088 << endl;
  
  
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
    float trigeff1LVal[nBinsMET][nBinsHT] = {{0.916547,0.987935,0.999278,0.999781},{0.952359,0.993993,0.998778,0.999788},{1,1,1,1},{1,1,1,1}};
    float trigeff1LErr[nBinsMET][nBinsHT] = {{0.0216782,0.0101562,0.0101094,0.0110318},{0.101613,0.0101077,0.0101169,0.0102007},{0.0151531,0.0110547,0.013888,0.014495},{0.0108621,0.0108621,0.0108621,0.0108621}};
    float trigeff0LVal[nBinsMET][nBinsHT] = {{0.859365,0.901058,1,1},{0.739399,1,1,1},{1,1,1,1},{1,1,1,1}};
    float trigeff0LErr[nBinsMET][nBinsHT] = {{0.107915,0.110075,0.0582974,0.0582974},{0.0871261,0.0829576,0.0582974,0.0582974},{0.0212098,0.0212098,0.0212098,0.0212098},{0.0212098,0.0212098,0.0212098,0.0212098}};
    for (int mbi = 0 ; mbi < nBinsMET ; mbi++) {
      for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
    	char parname[1000] ;
    	sprintf( parname, "trigeff_val_0L_M%d_H%d", mbi+1, hbi+1 ) ;
    	printf(" %s  :  %6.3f \n", parname, trigeff0LVal[mbi][hbi] ) ;
    	inFile << parname << "  \t" << trigeff0LVal[mbi][hbi] << endl;
	sprintf( parname, "trigeff_err_0L_M%d_H%d", mbi+1, hbi+1 ) ;
    	printf(" %s  :  %6.3f \n", parname, trigeff0LErr[mbi][hbi] ) ;
    	inFile << parname << "  \t" << trigeff0LErr[mbi][hbi] << endl;
      }
    }
    for (int mbi = 0 ; mbi < nBinsMET ; mbi++) {
      for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
    	char parname[1000] ;
    	sprintf( parname, "trigeff_val_1L_M%d_H%d", mbi+1, hbi+1 ) ;
    	printf(" %s  :  %6.3f \n", parname, trigeff1LVal[mbi][hbi] ) ;
    	inFile << parname << "  \t" << trigeff1LVal[mbi][hbi] << endl;
    	sprintf( parname, "trigeff_err_1L_M%d_H%d", mbi+1, hbi+1 ) ;
    	printf(" %s  :  %6.3f \n", parname, trigeff1LErr[mbi][hbi] ) ;
    	inFile << parname << "  \t" << trigeff1LErr[mbi][hbi] << endl;
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
  
////==========================================================================================
//
//void FillHTMET(TChain *chain, TH2F *histo, int si, int k) {

//   TObjArray *fileElements=chain->GetListOfFiles();
//   TIter next(fileElements);
//   TChainElement *chEl=0;
//   while (( chEl=(TChainElement*)next() )) {
//      TFile f(chEl->GetTitle());
//      TTree *tree = (TTree*)f.Get("tree");
//      SmallTree *t = new SmallTree(tree);
//      t->Loop(histo, si, k);
//   }

////     SmallTree *t = new SmallTree(chain);
////     t->Loop(histo, si, k);

//  return;    
//}
////==========================================================================================
