
#include "TStyle.h"
#include "TH1F.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>


    int nComps(3) ;
    char compname[3][100] ;
    int compcolor[3] ;
    double complumi[3] ;
    TChain* compchain[3] ;
    double compscale[3] ;

    TCanvas* dcan ;

   //----------------------------
   // prototypes

    void bookSet( const char* hname_base, const char* htitle_base, int nbins, double xmin, double xmax ) ;
    void fillSet( const char* hname_base, const char* varname, const char* cuts ) ;
    void saveHist(const char* filename, const char* pat) ;

   //----------------------------

    void dmc_fillplots() {

       gDirectory -> Delete("h*") ;

       gStyle -> SetOptStat(0) ;

       for ( int ci=0; ci<nComps; ci++ ) { compchain[ci] = new TChain("tree") ; }

       int compind ;

       compind = 0 ;
       sprintf( compname[compind], "data" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/MET_2012A_BLIND.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/MET_2012B_BLIND.root" ) ;
       complumi[compind] = 5.295 ;
       compcolor[compind] = 0 ;
       compscale[compind] = 1. ;

       compind = 1 ;
       sprintf( compname[compind], "ttbar" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/TT.root" ) ;
       complumi[compind] = 11.9 ;
       compcolor[compind] = 600-7 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;


       compind = 2 ;
       sprintf( compname[compind], "wjets" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WJets-250to300.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WJets-300to400.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WJets-400.root" ) ;
       complumi[compind] = 11.9 ;
       compcolor[compind] = 600-4 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;


       char basecuts_1lep[10000] ;
       sprintf( basecuts_1lep, "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nonb[10000] ;
       sprintf( basecuts_1lep_nonb, "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;



       dcan = (TCanvas*) gDirectory->FindObject("dcan") ;
       if ( dcan == 0x0 ) {
          dcan = new TCanvas("dcan","") ;
       }

       bookSet( "h_met_sl_all", "MET, SL", 30, 125., 500. ) ;
       bookSet( "h_met_sl_nb1", "MET, SL, nB=1", 30, 125., 500. ) ;
       bookSet( "h_met_sl_nb2", "MET, SL, nB=2", 30, 125., 500. ) ;
       bookSet( "h_met_sl_nb3", "MET, SL, nB>=3", 30, 125., 500. ) ;

       bookSet( "h_ht_sl_all", "HT, SL", 32, 400., 2000. ) ;
       bookSet( "h_ht_sl_nb1", "HT, SL, nB=1", 32, 400., 2000. ) ;
       bookSet( "h_ht_sl_nb2", "HT, SL, nB=2", 32, 400., 2000. ) ;
       bookSet( "h_ht_sl_nb3", "HT, SL, nB>=3", 32, 400., 2000. ) ;

       bookSet( "h_nb_sl_all", "nbjets, SL", 6, -0.5, 5.5 ) ;


       char cuts[10000] ;

       sprintf( cuts, "%s", basecuts_1lep ) ;
       fillSet( "h_met_sl_all", "MET", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_1lep ) ;
       fillSet( "h_met_sl_nb1", "MET", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_1lep ) ;
       fillSet( "h_met_sl_nb2", "MET", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_1lep ) ;
       fillSet( "h_met_sl_nb3", "MET", cuts ) ;

       sprintf( cuts, "%s", basecuts_1lep ) ;
       fillSet( "h_ht_sl_all", "HT", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_1lep ) ;
       fillSet( "h_ht_sl_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_1lep ) ;
       fillSet( "h_ht_sl_nb2", "HT", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_1lep ) ;
       fillSet( "h_ht_sl_nb3", "HT", cuts ) ;

       sprintf( cuts, "%s", basecuts_1lep_nonb ) ;
       fillSet( "h_nb_sl_all", "nB", cuts ) ;

       saveHist( "rootfiles/dmc_plots1.root", "h*" ) ;

    } // dmc_fillplots.

   //----------------------------


    void bookSet( const char* hname_base, const char* htitle_base, int nbins, double xmin, double xmax ) {

       printf("\n\n") ;
       for ( int ci=0; ci<nComps; ci++ ) {

          char hname[1000] ;
          char htitle[1000] ;

          sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
          sprintf( htitle, "%s, %s", htitle_base, compname[ci] ) ;

          printf(" bookSet: %s : %s\n", hname, htitle ) ;
          TH1F* hp = new TH1F( hname, htitle, nbins, xmin, xmax ) ;
          if ( compcolor[ci] > 0 ) { hp -> SetFillColor( compcolor[ci] ) ; }
          hp -> Sumw2() ;

       } // ci
       printf("\n\n") ;


    } // bookSet

   //----------------------------



    void fillSet( const char* hname_base, const char* varname, const char* cuts ) {

       printf("\n\n") ;
       for ( int ci=0; ci<nComps; ci++ ) {

          char hname[1000] ;
          char arg1[10000] ;

          sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
          sprintf( arg1, "%s>>%s", varname, hname ) ;
          printf(" %s : %s\n", arg1, cuts ) ;
          compchain[ci] -> Draw( arg1, cuts ) ;

          TH1F* hp = (TH1F*) gDirectory->FindObject( hname ) ;
          if ( hp == 0x0 ) { printf("\n *** missing hist! %s \n\n", hname ) ; return ; }

          hp -> Scale( compscale[ci] ) ;

          int nbinsx = hp->GetNbinsX() ;
          double last_val     = hp->GetBinContent( nbinsx  ) ;
          double last_err     = hp->GetBinError( nbinsx  ) ;
          double overflow_val = hp->GetBinContent( nbinsx + 1 ) ;
          double overflow_err = hp->GetBinError( nbinsx + 1 ) ;

          hp -> SetBinContent( nbinsx, last_val + overflow_val ) ;
          hp -> SetBinError( nbinsx, sqrt( pow(last_err,2) + pow(overflow_err,2) ) ) ;


          hp -> Draw("hist") ;
          hp -> Draw("esame") ;

          dcan->Update() ;

       } // ci
       printf("\n\n") ;


    } // fillSet

   //----------------------------

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

