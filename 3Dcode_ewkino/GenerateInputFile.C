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

  using std::stringstream ;
  using std::ofstream ;
  using std::endl ;
  using std::cout;
  using std::endl;
  using std::flush;

  void saveHist(const char* filename, const char* pat) ;
  TH1F* bookHist(const char* hname, const char* htitle, const char* selstring, int nbjet, int nBinsVar1, int nBinsVar2 ) ;


void GenerateInputFile( double mgl=-1., double mlsp=-1., double target_susy_all0lep=-1. ) {

  TChain* dyTree = new TChain("treeZ") ;
  int nAdded = dyTree->Add("files_20fb_v71_wip/DY-400.root") ;
  if ( nAdded <= 0 ) {
     printf("\n\n\n *** No treeZ in files_20fb_v71_wip/DY.root\n\n\n") ;
     return ;
  }

  double t2ttWeight(0.) ;
  TChain chainT2tt("tree") ;
  char susycutstring[1000] ;
  sprintf( susycutstring, "&&mgluino==%.0f&&mlsp==%.0f", mgl, mlsp ) ;
  TString susycut( susycutstring ) ;
  if ( mgl>0. && mlsp>0. ) {
     nAdded = chainT2tt.Add("files_20fb_v71_wip/T2tt.root") ;
     if ( nAdded <= 0 ) {
        printf("\n\n\n *** No tree in files_20fb_v71_wip/T2tt.root\n\n\n") ;
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
  chainQCD.Add("files_20fb_v71_wip/QCD-120to170.root");
  chainQCD.Add("files_20fb_v71_wip/QCD-170to300.root");
  chainQCD.Add("files_20fb_v71_wip/QCD-300to470.root");
  chainQCD.Add("files_20fb_v71_wip/QCD-470to600.root");
  chainQCD.Add("files_20fb_v71_wip/QCD-600to800.root");
  chainQCD.Add("files_20fb_v71_wip/QCD-800to1000.root");
  chainQCD.Add("files_20fb_v71_wip/QCD-1000to1400.root");
  chainQCD.Add("files_20fb_v71_wip/QCD-1400to1800.root");
  chainQCD.Add("files_20fb_v71_wip/QCD-1800.root");
  double kfactor_qcd = 1.8 ;
  printf("\n\n Rescaling QCD by %5.3f\n\n", kfactor_qcd ) ;

  TChain chainZnn("tree") ;
  chainZnn.Add("files_20fb_v71_wip/Zinv-100to200.root") ;
  chainZnn.Add("files_20fb_v71_wip/Zinv-200to400.root") ;
  chainZnn.Add("files_20fb_v71_wip/Zinv-400.root") ;

  TChain chainTT("tree") ;
  chainTT.Add("files_20fb_v71_wip/TT_FullLept.root") ;
  chainTT.Add("files_20fb_v71_wip/TT_SemiLept.root") ;
  chainTT.Add("files_20fb_v71_wip/TT_FullHad.root") ;
  double kfactor_tt = 0.90 ;
  printf("\n\n Rescaling ttbar by %5.3f\n\n", kfactor_tt ) ;

  TChain chainWJets("tree") ;
  chainWJets.Add("files_20fb_v71_wip/Wjets-250to300.root") ;
  chainWJets.Add("files_20fb_v71_wip/Wjets-300to400.root") ;
  chainWJets.Add("files_20fb_v71_wip/Wjets-400.root") ;
  chainWJets.Add("files_20fb_v71_wip/T-s.root") ;
  chainWJets.Add("files_20fb_v71_wip/T-t.root") ;
  chainWJets.Add("files_20fb_v71_wip/T-tW.root") ;
  chainWJets.Add("files_20fb_v71_wip/Tbar-s.root") ;
  chainWJets.Add("files_20fb_v71_wip/Tbar-t.root") ;
  chainWJets.Add("files_20fb_v71_wip/Tbar-tW.root") ;
  double kfactor_wjets = 0.90 ;
  printf("\n\n Rescaling wjets by %5.3f\n\n", kfactor_wjets ) ;


  //-- make chains of W+jets and single top separately.

  TChain chainWJetsOnly("tree") ;
  chainWJetsOnly.Add("files_20fb_v71_wip/Wjets-250to300.root") ;
  chainWJetsOnly.Add("files_20fb_v71_wip/Wjets-300to400.root") ;
  chainWJetsOnly.Add("files_20fb_v71_wip/Wjets-400.root") ;
  double kfactor_wjetsonly = 0.90 ;

  TChain chainSingletop("tree") ;
  chainSingletop.Add("files_20fb_v71_wip/T-s.root") ;
  chainSingletop.Add("files_20fb_v71_wip/T-t.root") ;
  chainSingletop.Add("files_20fb_v71_wip/T-tW.root") ;
  chainSingletop.Add("files_20fb_v71_wip/Tbar-s.root") ;
  chainSingletop.Add("files_20fb_v71_wip/Tbar-t.root") ;
  chainSingletop.Add("files_20fb_v71_wip/Tbar-tW.root") ;
  double kfactor_singletop = 0.90 ;

  TChain chainSingletop_s("tree") ;
  chainSingletop_s.Add("files_20fb_v71_wip/T-s.root") ;
  chainSingletop_s.Add("files_20fb_v71_wip/Tbar-s.root") ;
  double kfactor_singletops = 0.90 ;

  TChain chainSingletop_t("tree") ;
  chainSingletop_t.Add("files_20fb_v71_wip/T-t.root") ;
  chainSingletop_t.Add("files_20fb_v71_wip/Tbar-t.root") ;
  double kfactor_singletopt = 0.90 ;

  TChain chainSingletop_tw("tree") ;
  chainSingletop_tw.Add("files_20fb_v71_wip/T-tW.root") ;
  chainSingletop_tw.Add("files_20fb_v71_wip/Tbar-tW.root") ;
  double kfactor_singletoptw = 0.90 ;


  //include Z->ll in VV contribution
  TChain chainVV("tree");
  chainVV.Add("files_20fb_v71_wip/WW.root"); 
  chainVV.Add("files_20fb_v71_wip/WZ.root");
  chainVV.Add("files_20fb_v71_wip/ZZ.root");
  chainVV.Add("files_20fb_v71_wip/DY-200to400.root");
  chainVV.Add("files_20fb_v71_wip/DY-400.root");


  /*
  char qcdinputfile[9][1000] = {
    "files_20fb_v71_wip/QCD-120to170.root"
    ,"files_20fb_v71_wip/QCD-170to300.root"
    ,"files_20fb_v71_wip/QCD-300to470.root"
    ,"files_20fb_v71_wip/QCD-470to600.root"
    ,"files_20fb_v71_wip/QCD-600to800.root"
    ,"files_20fb_v71_wip/QCD-800to1000.root"
    ,"files_20fb_v71_wip/QCD-1000to1400.root"
    ,"files_20fb_v71_wip/QCD-1400to1800.root"
    ,"files_20fb_v71_wip/QCD-1800.root"
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
  */


  gROOT->Reset();

  const int nJetsCut = 3 ;     // #jets >= nJetsCut
  const int MTbCut = 0 ;       // cut on MTb

  bool ExcludeHiggs = true ;
  //TString sLooseHiggsCuts = "";
  TString sLooseHiggsCuts = "njets20<=7&&deltaRmax_hh<2.4&&((higgsMbb1MassDiff>95&&higgsMbb1MassDiff<145&&higgsMbb2MassDiff==-1)||(higgsMbb1MassDiff>95&&higgsMbb1MassDiff<145&&higgsMbb2MassDiff>95&&higgsMbb2MassDiff<145))&&METsig>30.&&deltaPhiStar>0.1&&pt_1st_leadJet<500&&" ;

  double minLeadJetPt = 70. ;
  double min3rdJetPt = 50. ;


  // read in the binning definitions from Binning.txt

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


  if ( in_Var3 != "nB" ) {
    cout << "\nUsing as 3rd variable something that is not nB is not supported! Exiting...\n\n" ;
    return ;
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

  Mbins[nBinsVar1] = 999999. ;
  Hbins[nBinsVar2] = 999999. ;
  
  inBinning >> label >> in_version ;
  const int version = in_version ;

  if ( !label.Contains("version") ) {
    cout << "\n\n Found inconsistency in Binning.txt, check number of bins and Var1 and Var2 lower bounds\n\n" << endl ;
    return ;
  }

  inBinning.close();


  TString sMbins[nBinsVar1];
  TString sHbins[nBinsVar2];
  TString sBbins[4] = {"_1b","_2b","_3b","_4b"};

  TString cMbins[nBinsVar1];
  TString cHbins[nBinsVar2];

  for (int i = 0 ; i < nBinsVar1 ; i++) {
    TString base = "_M";
    stringstream sbin;
    sbin << i+1;
    base += sbin.str();
    sMbins[i] = base;
    base = "&&"+sVar1+">";
    stringstream cbin;
    cbin << Mbins[i] << "&&" << sVar1 << "<" << Mbins[i+1];
    base += cbin.str();
    cMbins[i] = base;
  }

  for (int j = 0 ; j < nBinsVar2 ; j++) {
    TString base = "_H";
    stringstream sbin;
    sbin << j+1;
    base += sbin.str();
    sHbins[j] = base;
    base = "&&"+sVar2+">";
    stringstream cbin;
    cbin << Hbins[j] << "&&" << sVar2 << "<" << Hbins[j+1];
    base += cbin.str();
    cHbins[j] = base;
  }

  float dummyOne = 1.0;
  float dummyErr = 0.1;

  ofstream inFile;
  char outfile[10000] ;
  if ( mgl > 0. && mlsp > 0. ) {
     if ( target_susy_all0lep > 0. ) {
       sprintf( outfile, "InputWT2tt-mgl%.0f-mlsp%.0f-%.0fevts-%s-%d-%s-%d-nB%d-v%d.dat", mgl, mlsp, target_susy_all0lep, aVar1.Data(), nBinsVar1, aVar2.Data(), nBinsVar2, nBinsBjets, version  ) ;
     } else {
       sprintf( outfile, "InputWT2tt-mgl%.0f-mlsp%.0f-%s-%d-%s-%d-nB%d-v%d.dat", mgl, mlsp, aVar1.Data(), nBinsVar1, aVar2.Data(), nBinsVar2, nBinsBjets, version  ) ;
     }
  } else {
    sprintf( outfile, "Input-%s-%d-%s-%d-nB%d-v%d.dat", aVar1.Data(), nBinsVar1, aVar2.Data(), nBinsVar2, nBinsBjets, version  ) ;
  }
  inFile.open( outfile );

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


  bool useBtagSF = true ;
  char bcut[nBinsBjets][100] ;
  char bcutSF[nBinsBjets][100] ;

  sprintf(bcut[0],"nB==1");
  sprintf(bcutSF[0],"prob1");


  if ( nBinsBjets == 2 ) {
    sprintf(bcut[1],"nB>=2");
    sprintf(bcutSF[1],"(prob2+probge3)");
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

  if ( ExcludeHiggs ) {
    sprintf( commoncuts, "%s!passHiggsSel&&maxChNMultDiff<40&&pfOcaloMET<2.0&&nJets>=%d&&MT_bestCSV>%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f&&pt_3rd_leadJet>%.0f)",
	     sLooseHiggsCuts.Data(),nJetsCut, MTbCut, minLeadJetPt, minLeadJetPt, min3rdJetPt ) ;
  }
  else {
    sprintf( commoncuts, "%smaxChNMultDiff<40&&pfOcaloMET<2.0&&nJets>=%d&&MT_bestCSV>%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f&&pt_3rd_leadJet>%.0f)",
	     sLooseHiggsCuts.Data(),nJetsCut, MTbCut, minLeadJetPt, minLeadJetPt, min3rdJetPt ) ;
  }


  const int nSel(3) ;
  char selname[3][100] = { "0lep", "1lepSig", "1lep" } ;
  
  char selcuts[3][10000] ;
  sprintf( selcuts[0], "minDelPhiN>4&&nMu==0&&nEl==0&&nIsoTrk==0" ) ; //--- 0lep
  sprintf( selcuts[1], "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&MT>100" ) ; //--- 1lepSig
  sprintf( selcuts[2], "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&MT<100" ) ; //--- 1lep


  //--- Output histograms.
  
  TH1F* hmctruth_susy[nSel][nBinsBjets] ;
  TH1F* hmctruth_ttwj[nSel][nBinsBjets] ;
  TH1F* hmctruth_ttbar[nSel][nBinsBjets] ;
  TH1F* hmctruth_wjets[nSel][nBinsBjets] ;
  TH1F* hmctruth_wjetsonly[nSel][nBinsBjets] ;
  TH1F* hmctruth_singletop[nSel][nBinsBjets] ;
  TH1F* hmctruth_singletops[nSel][nBinsBjets] ;
  TH1F* hmctruth_singletopt[nSel][nBinsBjets] ;
  TH1F* hmctruth_singletoptw[nSel][nBinsBjets] ;
  TH1F* hmctruth_qcd[nSel][nBinsBjets] ;
  TH1F* hmctruth_znn[nSel][nBinsBjets] ;
  TH1F* hmctruth_qzo[nSel][nBinsBjets] ;
  TH1F* hmctruth_allsm[nSel][nBinsBjets] ;
  TH1F* hmctruth_all[nSel][nBinsBjets] ;


  for ( int si=0; si<nSel; si++ ) {
    for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
      
      char hname[1000] ;
      char htitle[1000] ;
      sprintf( htitle, "%s, %d btag", selname[si], bbi+1 ) ;
      
      sprintf( hname, "hmctruth_susy_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_susy[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_ttwj_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_ttwj[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_ttbar_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_ttbar[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_wjets_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_wjets[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_wjetsonly_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_wjetsonly[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_singletop_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_singletop[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_singletops_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_singletops[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_singletopt_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_singletopt[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_singletoptw_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_singletoptw[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_qcd_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_qcd[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_znn_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_znn[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_qzo_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_qzo[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_allsm_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_allsm[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      sprintf( hname, "hmctruth_all_%s_%db", selname[si], bbi+1 ) ;
      hmctruth_all[si][bbi] = bookHist( hname, htitle, selname[si], bbi+1, nBinsVar1, nBinsVar2 ) ;
      
      
    } // bbi.
  } // si.
  

   //--- histograms used in getting the observables.

   TH2F* h_tt[10] ;
   TH2F* h_wjets[10] ;
   TH2F* h_wjetsonly[10] ;
   TH2F* h_singletop[10] ;
   TH2F* h_singletops[10] ;
   TH2F* h_singletopt[10] ;
   TH2F* h_singletoptw[10] ;
   TH2F* h_qcd[10] ;
   TH2F* h_znn[10] ;
   TH2F* h_qzo[10];
   TH2F* h_susy[10] ;
   TH2F* h_mc[10] ;

   TH2F *h_znn_0lep_1b;
   h_znn_0lep_1b  = new TH2F( "h_znn_0lep_1b", "h_znn_0lep_1b", nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
   h_znn_0lep_1b->Sumw2() ;

   for ( int bi=0; bi<nBinsBjets; bi++ ) {

     char hname[100] ;
     
     sprintf( hname, "h_tt_%db", bi+1 ) ;
     h_tt[bi]   = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_tt[bi] -> Sumw2() ;
     
     sprintf( hname, "h_wjets_%db", bi+1 ) ;
     h_wjets[bi]   = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_wjets[bi] -> Sumw2() ;
     
     sprintf( hname, "h_wjetsonly_%db", bi+1 ) ;
     h_wjetsonly[bi]   = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_wjetsonly[bi] -> Sumw2() ;
     
     sprintf( hname, "h_singletop_%db", bi+1 ) ;
     h_singletop[bi]   = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_singletop[bi] -> Sumw2() ;
     
     sprintf( hname, "h_singletops_%db", bi+1 ) ;
     h_singletops[bi]   = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_singletops[bi] -> Sumw2() ;
     
     sprintf( hname, "h_singletopt_%db", bi+1 ) ;
     h_singletopt[bi]   = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_singletopt[bi] -> Sumw2() ;
     
     sprintf( hname, "h_singletoptw_%db", bi+1 ) ;
     h_singletoptw[bi]   = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_singletoptw[bi] -> Sumw2() ;

     sprintf( hname, "h_qcd_%db", bi+1 ) ;
     h_qcd[bi]  = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_qcd[bi] -> Sumw2() ;

     sprintf( hname, "h_znn_%db", bi+1 ) ;
     h_znn[bi]  = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_znn[bi] -> Sumw2() ;
     
     sprintf( hname, "h_qzo_%db", bi+1 ) ;
     h_qzo[bi]  = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_qzo[bi] -> Sumw2() ;
     
     sprintf( hname, "h_susy_%db", bi+1 ) ;
     h_susy[bi] = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
     h_susy[bi] -> Sumw2() ;

     sprintf( hname, "h_mc_%db", bi+1 ) ;
     h_mc[bi]   = new TH2F( hname, hname , nBinsVar1, Mbins, nBinsVar2, Hbins ) ;
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
       chainTT.Project (hname,s2DVars,allcuts);
       
       h_tt[k] -> Scale( kfactor_tt ) ;
       printf("    %12s %7.1f events\n", hname, h_tt[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_wjets_%db", k+1 ) ;
       chainWJets.Project(hname,s2DVars,allcuts);
       h_wjets[k] -> Scale( kfactor_wjets ) ;
       printf("    %12s %7.1f events\n", hname, h_wjets[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_wjetsonly_%db", k+1 ) ;
       chainWJetsOnly.Project(hname,s2DVars,allcuts);
       h_wjetsonly[k] -> Scale( kfactor_wjetsonly ) ;
       printf("    %12s %7.1f events\n", hname, h_wjetsonly[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_singletop_%db", k+1 ) ;
       chainSingletop.Project(hname,s2DVars,allcuts);
       h_singletop[k] -> Scale( kfactor_singletop ) ;
       printf("    %12s %7.1f events\n", hname, h_singletop[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_singletops_%db", k+1 ) ;
       chainSingletop_s.Project(hname,s2DVars,allcuts);
       h_singletops[k] -> Scale( kfactor_singletops ) ;
       printf("    %12s %7.1f events\n", hname, h_singletops[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_singletopt_%db", k+1 ) ;
       chainSingletop_t.Project(hname,s2DVars,allcuts);
       h_singletopt[k] -> Scale( kfactor_singletopt ) ;
       printf("    %12s %7.1f events\n", hname, h_singletopt[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_singletoptw_%db", k+1 ) ;
       chainSingletop_tw.Project(hname,s2DVars,allcuts);
       h_singletoptw[k] -> Scale( kfactor_singletoptw ) ;
       printf("    %12s %7.1f events\n", hname, h_singletoptw[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_qcd_%db", k+1 ) ;
       chainQCD.Project(hname,s2DVars,allcuts);
       h_qcd[k] -> Scale( kfactor_qcd ) ;
       printf("    %12s %7.1f events\n", hname, h_qcd[k]->Integral() ) ; cout << flush ;
       
       sprintf( hname, "h_znn_%db", k+1 ) ;
       chainZnn.Project(hname,s2DVars,allcuts);
       printf("    %12s %7.1f events\n", hname, h_znn[k]->Integral() ) ; cout << flush ;
       
       if ( si == 0 && k == 0 ) {
	 chainZnn.Project("h_znn_0lep_1b",s2DVars,allcuts);
       }


       // this component includes QCD, Z -> inv. and all the other minor components
       sprintf( hname, "h_qzo_%db", k+1 ) ;
       chainVV.Project(hname,s2DVars,allcuts);
       h_qzo[k]->Add(h_qcd[k],1.);
       h_qzo[k]->Add(h_znn[k],1.);
       printf("    %12s %7.1f events\n", hname, h_qzo[k]->Integral() ) ; cout << flush ;

       if ( mgl > 0. ) {
	 sprintf( hname, "h_susy_%db", k+1 ) ;
	 chainT2tt.Project(hname,s2DVars,allsusycuts);
	 h_susy[k]->Scale( t2ttWeight ) ;
	 printf("    %12s %7.1f events\n", hname, h_susy[k]->Integral() ) ; cout << flush ;
	 if (si==0) nSusyTotal += h_susy[k]->Integral();
       }
      } // k
      if (si==0) printf("N_Susy_Total = %7.1f events", nSusyTotal); cout << flush;
      printf("\n\n") ;


      for (int i = 0 ; i < nBinsVar1 ; i++) {
        for (int j = 0 ; j < nBinsVar2 ; j++) {
          for (int k = 0 ; k < nBinsBjets ; k++) {

	    char obsname[1000] ;
	    sprintf( obsname, "N_%s_M%d_H%d_%db", selname[si], i+1, j+1, k+1 ) ;
	    
	    double ttval = h_tt[k] -> GetBinContent( i+1, j+1 ) ;
	    double tterr = h_tt[k] -> GetBinError(   i+1, j+1 ) ;
	    
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
	    
	    double qzoval = h_qzo[k] -> GetBinContent( i+1, j+1 ) ;
	    double qzoerr = h_qzo[k] -> GetBinError(   i+1, j+1 ) ;
	    
	    double susyval = h_susy[k] -> GetBinContent( i+1, j+1 ) ;
	    double susyerr = h_susy[k] -> GetBinError(   i+1, j+1 ) ;


	    printf(" N_%s, tt          var1,var2,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, ttval,tterr) ; cout << flush ;
	    printf(" N_%s, wjets       var1,var2,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, wjetsval,wjetserr) ; cout << flush ;
	    printf(" N_%s,(wjets only) var1,var2,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, wjetsonlyval,wjetsonlyerr) ; cout << flush ;
	    printf(" N_%s,(singletop ) var1,var2,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, singletopval,singletoperr) ; cout << flush ;
	    printf(" N_%s,((sng.t,s))  var1,var2,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, singletopsval,singletopserr) ; cout << flush ;
	    printf(" N_%s,((sng.t,t))  var1,var2,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, singletoptval,singletopterr) ; cout << flush ;
	    printf(" N_%s,((sng.t,tW)) var1,var2,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, singletoptwval,singletoptwerr) ; cout << flush ;
	    printf(" N_%s, qcd         var1,var2,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, qcdval,qcderr) ; cout << flush ;
	    printf(" N_%s, znn         var1,var2,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, znnval,znnerr) ; cout << flush ;
	    printf(" N_%s, qzo         var1,var2,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, qzoval,qzoerr) ; cout << flush ;
	    if ( mgl>0. ) {
	      if ( target_susy_all0lep > 0. ) {
		susyval = susyval * (target_susy_all0lep/nSusyTotal);
		susyerr = susyerr * (target_susy_all0lep/nSusyTotal);
	      }
	      printf(" N_%s, susy   var1,var2,nbjet bin (%d,%d,%d)  --  npass=%7.1f +/- %6.1f\n", selname[si], i,j,k, susyval,susyerr) ; cout << flush ;
	    }
	    printf("\n") ;

	    double allval = ttval + wjetsval + qzoval + susyval ;

	    inFile << obsname << "  \t" << floor(allval+0.5) << endl;
	    
	    int histbin = 1 + (nBinsVar2+1)*i + j + 1 ;
	    
	    hmctruth_ttwj[si][k]  -> SetBinContent( histbin, ttval + wjetsval  ) ;
	    hmctruth_ttwj[si][k]  -> SetBinError(   histbin, sqrt( pow(tterr,2) + pow(wjetserr,2) )  ) ;
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
	    hmctruth_qzo[si][k]    -> SetBinContent( histbin, qzoval ) ;
	    hmctruth_qzo[si][k]    -> SetBinError(   histbin, qzoerr ) ;
	    hmctruth_susy[si][k]  -> SetBinContent( histbin, susyval  ) ;
	    hmctruth_susy[si][k]  -> SetBinError(   histbin, susyerr  ) ;
	    hmctruth_allsm[si][k] -> SetBinContent( histbin, ttval+wjetsval+qzoval ) ;
	    hmctruth_allsm[si][k] -> SetBinError(   histbin, sqrt( pow(tterr,2) + pow(wjetserr,2) + pow(qzoerr,2) ) ) ;
	    hmctruth_all[si][k]   -> SetBinContent( histbin, ttval+wjetsval+qzoval+susyval ) ;
	    hmctruth_all[si][k]   -> SetBinError(   histbin, sqrt( pow(tterr,2) + pow(wjetserr,2) + pow(qzoerr,2) + pow(susyerr,2) ) ) ;


	  }
	}
      }

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
	h_qzo[k] -> Reset() ;
	h_susy[k] -> Reset() ;
	
      }


   } // si.


   printf("\n\n-----------------------------------------------------------------\n\n") ; cout << flush ;

   
  
   // sf_ttwj
   
   for (int i = 0 ; i < nBinsVar1 ; i++) {
     for (int j = 0 ; j < nBinsVar2 ; j++) {
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
   
   for (int i = 0 ; i < nBinsVar1 ; i++) {
     for (int j = 0 ; j < nBinsVar2 ; j++) {
       for (int k = 0 ; k < nBinsBjets ; k++) {
	 
	 TString sf_ttwj_slSig = "sf_ttwj_slSig" ;
	 sf_ttwj_slSig = sf_ttwj_slSig+sMbins[i]+sHbins[j]+sBbins[k] ;
	 
	 inFile << sf_ttwj_slSig << "  \t" << dummyOne << endl;	
	  
	 sf_ttwj_slSig = sf_ttwj_slSig+"_err" ;
	 inFile << sf_ttwj_slSig << "  \t" << dummyErr << endl;	
	 
       }
     }
   }
    


   //--- ttwj MC ZL 3b/2b ratios
   
   // 0-lep
   
   if ( nBinsBjets > 2 ) {
     
     for (int mbi = 0 ; mbi < nBinsVar1 ; mbi++) {
       for (int hbi = 0 ; hbi < nBinsVar2 ; hbi++) {
	 int hbin = 1 + (nBinsVar2+1)*mbi + hbi + 1 ;
	 
	 float zl2bval = hmctruth_ttwj[0][1] -> GetBinContent( hbin ) ;
	 float zl2berr = hmctruth_ttwj[0][1] -> GetBinError(   hbin ) ;
	 float zl3bval  = hmctruth_ttwj[0][2] -> GetBinContent( hbin ) ;
	 float zl3berr  = hmctruth_ttwj[0][2] -> GetBinError(   hbin ) ;
	 
	 float bbboverbb = 0. ;
	 float bbboverbberr = 0. ;
	 
	 if ( zl2bval > 0. && zl3bval > 0. ) {
	   bbboverbb = zl3bval / zl2bval ;
	   bbboverbberr = bbboverbb * sqrt( pow((zl2berr/zl2bval),2) + pow((zl3berr/zl3bval),2) ) ;
	 }
	 
	 char parname[1000] ;
	 
	 sprintf( parname, "ttwj_mc_3bover2b_ratio_M%d_H%d", mbi+1, hbi+1 ) ;
	 printf(" %s  :  %6.3f +/- %5.3f\n", parname, bbboverbb, bbboverbberr ) ;
	 inFile << parname << "  \t" << bbboverbb << endl;
	 
	 sprintf( parname, "ttwj_mc_3bover2b_ratio_M%d_H%d_err", mbi+1, hbi+1 ) ;
	 inFile << parname << "  \t" << bbboverbberr << endl;
	 
       } // hbi
     } // mbi
     

     if ( nBinsBjets > 3 ) {
       
       //--- ttwj MC ZL 4b/2b ratios
       
       for (int mbi = 0 ; mbi < nBinsVar1 ; mbi++) {
	 for (int hbi = 0 ; hbi < nBinsVar2 ; hbi++) {
	   int hbin = 1 + (nBinsVar2+1)*mbi + hbi + 1 ;
	   
	   float zl2bval = hmctruth_ttwj[0][1] -> GetBinContent( hbin ) ;
	   float zl2berr = hmctruth_ttwj[0][1] -> GetBinError(   hbin ) ;
	   float zl4bval  = hmctruth_ttwj[0][3] -> GetBinContent( hbin ) ;
	   float zl4berr  = hmctruth_ttwj[0][3] -> GetBinError(   hbin ) ;
	   
	   float bbbboverbb = 0. ;
	   float bbbboverbberr = 0. ;
	   
	   if ( zl2bval > 0. && zl4bval > 0. ) {
	     bbbboverbb = zl4bval / zl2bval ;
	     bbbboverbberr = bbbboverbb * sqrt( pow((zl2berr/zl2bval),2) + pow((zl4berr/zl4bval),2) ) ;
	   }
	   
	   char parname[1000] ;
	   
	   sprintf( parname, "ttwj_mc_4bover2b_ratio_M%d_H%d", mbi+1, hbi+1 ) ;
	   printf(" %s  :  %6.3f +/- %5.3f\n", parname, bbbboverbb, bbbboverbberr ) ;
	   inFile << parname << "  \t" << bbbboverbb << endl;
	   
	   sprintf( parname, "ttwj_mc_4bover2b_ratio_M%d_H%d_err", mbi+1, hbi+1 ) ;
	   inFile << parname << "  \t" << bbbboverbberr << endl;
	   
	 } // hbi
       } // mbi
       
     }
   }


   // 1-lepSig
   
   if ( nBinsBjets > 2 ) {
     
     for (int mbi = 0 ; mbi < nBinsVar1 ; mbi++) {
       for (int hbi = 0 ; hbi < nBinsVar2 ; hbi++) {
	 int hbin = 1 + (nBinsVar2+1)*mbi + hbi + 1 ;
	 
	 float zl2bval = hmctruth_ttwj[1][1] -> GetBinContent( hbin ) ;
	 float zl2berr = hmctruth_ttwj[1][1] -> GetBinError(   hbin ) ;
	 float zl3bval  = hmctruth_ttwj[1][2] -> GetBinContent( hbin ) ;
	 float zl3berr  = hmctruth_ttwj[1][2] -> GetBinError(   hbin ) ;
	 
	 float bbboverbb = 0. ;
	 float bbboverbberr = 0. ;
	 
	 if ( zl2bval > 0. && zl3bval > 0. ) {
	   bbboverbb = zl3bval / zl2bval ;
	   bbboverbberr = bbboverbb * sqrt( pow((zl2berr/zl2bval),2) + pow((zl3berr/zl3bval),2) ) ;
	 }
	 
	 char parname[1000] ;
	 
	 sprintf( parname, "ttwjSlSig_mc_3bover2b_ratio_M%d_H%d", mbi+1, hbi+1 ) ;
	 printf(" %s  :  %6.3f +/- %5.3f\n", parname, bbboverbb, bbboverbberr ) ;
	 inFile << parname << "  \t" << bbboverbb << endl;
	 
	 sprintf( parname, "ttwjSlSig_mc_3bover2b_ratio_M%d_H%d_err", mbi+1, hbi+1 ) ;
	 inFile << parname << "  \t" << bbboverbberr << endl;
	 
       } // hbi
     } // mbi
      

     if ( nBinsBjets > 3 ) {
       
       //--- ttwj MC ZL 4b/2b ratios
	
       for (int mbi = 0 ; mbi < nBinsVar1 ; mbi++) {
	 for (int hbi = 0 ; hbi < nBinsVar2 ; hbi++) {
	   int hbin = 1 + (nBinsVar2+1)*mbi + hbi + 1 ;
	   
	   float zl2bval = hmctruth_ttwj[1][1] -> GetBinContent( hbin ) ;
	   float zl2berr = hmctruth_ttwj[1][1] -> GetBinError(   hbin ) ;
	   float zl4bval  = hmctruth_ttwj[1][3] -> GetBinContent( hbin ) ;
	   float zl4berr  = hmctruth_ttwj[1][3] -> GetBinError(   hbin ) ;
	   
	   float bbbboverbb = 0. ;
	   float bbbboverbberr = 0. ;
	   
	   if ( zl2bval > 0. && zl4bval > 0. ) {
	     bbbboverbb = zl4bval / zl2bval ;
	     bbbboverbberr = bbbboverbb * sqrt( pow((zl2berr/zl2bval),2) + pow((zl4berr/zl4bval),2) ) ;
	   }
	   
	   char parname[1000] ;
	   
	   sprintf( parname, "ttwjSlSig_mc_4bover2b_ratio_M%d_H%d", mbi+1, hbi+1 ) ;
	   printf(" %s  :  %6.3f +/- %5.3f\n", parname, bbbboverbb, bbbboverbberr ) ;
	   inFile << parname << "  \t" << bbbboverbb << endl;
	   
	   sprintf( parname, "ttwjSlSig_mc_4bover2b_ratio_M%d_H%d_err", mbi+1, hbi+1 ) ;
	   inFile << parname << "  \t" << bbbboverbberr << endl;
	   
	 } // hbi
       } // mbi
       
     }
   }

    
   // print out combined mc truth values for QCD + Z -> inv. + others
   for ( int si=0; si<nSel; si++ ) {
     for (int mbi = 0 ; mbi < nBinsVar1 ; mbi++) {
       for (int hbi = 0 ; hbi < nBinsVar2 ; hbi++) {
	 int hbin = 1 + (nBinsVar2+1)*mbi + hbi + 1 ;
	 for (int bbi = 0 ; bbi < nBinsBjets ; bbi++) {
	   
	   float qzovalue = hmctruth_qzo[si][bbi] -> GetBinContent( hbin ) ;
	   char parname[1000] ;
	   
	   sprintf( parname, "N_QZOmc_%s_M%d_H%d_%db", selname[si], mbi+1, hbi+1, bbi+1 ) ;
	   printf(" %s  :  %6.3f \n", parname, qzovalue ) ;
	   inFile << parname << "  \t" << qzovalue << endl;
	    

	 } // bbi
       } // hbi
     } // mbi
   } // si


   float trigeff1LVal[nBinsVar2][nBinsVar1] ;
   float trigeff1LErr[nBinsVar2][nBinsVar1] ;
   float trigeff0LVal[nBinsVar2][nBinsVar1] ;
   float trigeff0LErr[nBinsVar2][nBinsVar1] ;
   
   for (int mbi = 0 ; mbi < nBinsVar1 ; mbi++) {
     for (int hbi = 0 ; hbi < nBinsVar2 ; hbi++) {
       
       trigeff1LVal[hbi][mbi] = 0.98 ;
       trigeff1LErr[hbi][mbi] = 0.015 ;
       trigeff0LVal[hbi][mbi] = 0.95 ;
       trigeff0LErr[hbi][mbi] = 0.05 ;
       
     }
   }
   

   for (int mbi = 0 ; mbi < nBinsVar1 ; mbi++) {
     for (int hbi = 0 ; hbi < nBinsVar2 ; hbi++) {
       char parname[1000] ;
       sprintf( parname, "trigeff_val_0L_M%d_H%d", mbi+1, hbi+1 ) ;
       printf(" %s  :  %6.3f \n", parname, trigeff0LVal[hbi][mbi] ) ;
       inFile << parname << "  \t" << trigeff0LVal[hbi][mbi] << endl;
       sprintf( parname, "trigeff_err_0L_M%d_H%d", mbi+1, hbi+1 ) ;
       printf(" %s  :  %6.3f \n", parname, trigeff0LErr[hbi][mbi] ) ;
       inFile << parname << "  \t" << trigeff0LErr[hbi][mbi] << endl;
     }
   }
   for (int mbi = 0 ; mbi < nBinsVar1 ; mbi++) {
     for (int hbi = 0 ; hbi < nBinsVar2 ; hbi++) {
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



    gSystem->Exec("mkdir -p rootfiles") ;
    char outHistName[1000] ;
    if ( mgl>0. && mlsp>0. ) {
       if ( target_susy_all0lep > 0 ) {
	 sprintf( outHistName, "rootfiles/gi-plots-wsusy-mgl%.0f-mlsp%.0f-%.0fevts-%s-%d-%s-%d-nB%d-v%d.root", mgl, mlsp, target_susy_all0lep, aVar1.Data(), nBinsVar1, aVar2.Data(), nBinsVar2, nBinsBjets, version ) ;
       } else {
	 sprintf( outHistName, "rootfiles/gi-plots-wsusy-mgl%.0f-mlsp%.0f-%s-%d-%s-%d-nB%d-v%d.root", mgl, mlsp, aVar1.Data(), nBinsVar1, aVar2.Data(), nBinsVar2, nBinsBjets, version ) ;
       }
    } else {
      sprintf( outHistName, "rootfiles/gi-plots-%s-%d-%s-%d-nB%d-v%d.root", aVar1.Data(), nBinsVar1, aVar2.Data(), nBinsVar2, nBinsBjets, version ) ;
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
  

    TH1F* bookHist(const char* hname, const char* htitle, const char* selstring, int nbjet, int nBinsVar1, int nBinsVar2 ) {
  
       int nbins = nBinsVar1*(nBinsVar2+1) + 1 ;
  
       TH1F* retVal = new TH1F( hname, htitle, nbins, 0.5 + 0.1*(nbjet-2), nbins+0.5 + 0.1*(nbjet-2) ) ;
       TAxis* xaxis = retVal->GetXaxis() ;
  
       for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
          for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
             int histbin = 1 + (nBinsVar2+1)*mbi + hbi + 1 ;
             char binlabel[1000] ;
             sprintf( binlabel, "%s_M%d_H%d_%db", selstring, mbi+1, hbi+1, nbjet ) ;
             xaxis->SetBinLabel( histbin, binlabel ) ;
          } // hbi.
       } // mbi.
  
       retVal->SetLabelSize(0.055,"x") ;
       xaxis->LabelsOption("v") ;
  
       return retVal ;
  
    }
