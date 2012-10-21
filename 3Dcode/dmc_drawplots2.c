//
//  dmc_drawplots2( infile, savePdf, recycleCanvas, histnamepat )
//
//    If you run it with the default arguments, it will produce pdf files
//    for all of the plots, recycling one canvas.
//
//    If you set recycleCanvas = false, a new canvas will open for each plot.
//    this is not recommended unles you specify a subset of plots with
//    the histnamepat argument.  The total is too big now.
//
//    If histnamepat is set, only plots that contain the histnamepat
//    substring will be drawn.
//
//    The kfactor array contains fudge factors for the samples.
//    Set them all to 1 if you want the raw MC normalization.
//

#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TText.h"
#include "TString.h"

#include <iostream>


  //-- met4-ht4-v15
      const int nBinsMET   = 4 ;
      const int nBinsHT    = 4 ;
      float Mbins[nBinsMET+1] = {125.,150.,250.,350.,99999.};
      float Hbins[nBinsHT+1] = {400.,500.,800.,1000.,99999.};


    bool recycleCanvas(false) ;

    int nComps(7) ;
    char compname[7][100] ;
    double kfactor[7] ;

    bool kfactorAlreadyApplied ;

    bool savePdf ;
    char inrootfile[10000] ;

  //----------
  // prototypes

   void drawSet( const char* hname_base, const char* xtitle ) ;
   void drawSetLHB( const char* hname_base ) ;
   void flattenLHB( const char* hname_base ) ;
   void drawMetBinEdges( TH1F* hp ) ;
   void drawHtBinEdges( TH1F* hp ) ;

   void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;

  //----------

    void dmc_drawplots2( const char* infile = "rootfiles/dmc_plots2_all.root", bool arg_savePdf = true, bool arg_recycleCanvas=true, const char* histnamepat = "" ) {

       recycleCanvas = arg_recycleCanvas ;
       savePdf = arg_savePdf ;
       sprintf( inrootfile, "%s", infile ) ;

       gDirectory->Delete("h*") ;

       gStyle -> SetOptStat(0) ;
       gStyle -> SetTitleH(0.035 ) ;

       sprintf( compname[0], "data" ) ;      kfactor[0] = 1.00 ;
       sprintf( compname[1], "diboson" ) ;   kfactor[1] = 1.00 ;
       sprintf( compname[2], "znn" ) ;       kfactor[2] = 1.00 ;
       sprintf( compname[3], "qcd" ) ;       kfactor[3] = 1.80 ; //**********************
       sprintf( compname[4], "singlet" ) ;   kfactor[4] = 1.00 ;
       sprintf( compname[5], "wjets" ) ;     kfactor[5] = 0.90 ; //**********************
       sprintf( compname[6], "ttbar" ) ;     kfactor[6] = 0.90 ; //**********************


     //---------

       char histnamelist[500][100] ;
       char histxtitlelist[500][100] ;
       int hind(0) ;

       sprintf( histnamelist[hind], "h_met_sl_all" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nb2" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nb3" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_sl_htlt800_all" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_htlt800_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_htlt800_nb2" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_htlt800_nb3" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_sl_htgt800_all" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_htgt800_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_htgt800_nb2" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_htgt800_nb3" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_all" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_nb2" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_nb3" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht1_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht2_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht3_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht4_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht1_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht2_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht3_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht4_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_sl_ht1_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht2_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht3_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht4_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_sl_ht1_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht2_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht3_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht4_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht1_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht2_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht3_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht4_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht1_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht2_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht3_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_nomctrigcorr_ht4_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_sl_ht1_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht2_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht3_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht4_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_sl_ht1_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht2_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht3_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_sl_ht4_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_all" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nb2" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nb3" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_htlt800_all" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_htlt800_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_htlt800_nb2" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_htlt800_nb3" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_htgt800_all" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_htgt800_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_htgt800_nb2" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_htgt800_nb3" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_all" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_nb2" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_nb3" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht1_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht2_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht3_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht4_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht1_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht2_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht3_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht4_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_ht1_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht2_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht3_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht4_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_ht1_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht2_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht3_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht4_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht1_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht2_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht3_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht4_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht1_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht2_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht3_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_nomctrigcorr_ht4_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_ht1_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht2_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht3_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht4_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_ldp_ht1_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht2_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht3_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_ldp_ht4_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;


       sprintf( histnamelist[hind], "h_met_zl_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nb2" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nb3" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_nb2" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_nb3" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht1_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht2_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht3_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht4_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht1_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht2_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht3_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht4_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_zl_ht1_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht2_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht3_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht4_nb1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_zl_ht1_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht2_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht3_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht4_nb1_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht1_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht2_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht3_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht4_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht1_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht2_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht3_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_nomctrigcorr_ht4_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_zl_ht1_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht2_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht3_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht4_nb0" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_met_zl_ht1_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht2_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht3_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_met_zl_ht4_nb0_wide1" ) ; sprintf( histxtitlelist[hind], "MET (GeV)" ) ; hind++ ;



       sprintf( histnamelist[hind], "h_ht_sl_all" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nb2" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nb3" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_met1_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_met2_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_met3_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_met4_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_sl_met1_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_met2_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_met3_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_met4_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_met1_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_met2_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_met3_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_met4_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_sl_met1_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_met2_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_met3_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_met4_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_sl_metlt250_all" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_metlt250_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_metlt250_nb2" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_metlt250_nb3" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_sl_metgt250_all" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_metgt250_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_metgt250_nb2" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_metgt250_nb3" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_all" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_nb2" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_sl_nomctrigcorr_nb3" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;



       sprintf( histnamelist[hind], "h_ht_ldp_all" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nb2" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nb3" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_met1_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_met2_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_met3_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_met4_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_ldp_met1_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_met2_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_met3_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_met4_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_met1_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_met2_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_met3_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_met4_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_ldp_met1_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_met2_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_met3_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_met4_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_ldp_metlt250_all" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_metlt250_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_metlt250_nb2" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_metlt250_nb3" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_ldp_metgt250_all" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_metgt250_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_metgt250_nb2" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_metgt250_nb3" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_all" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_nb2" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_ldp_nomctrigcorr_nb3" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;


       sprintf( histnamelist[hind], "h_ht_zl_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nb2" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nb3" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_nb2" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_nb3" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_met1_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_met2_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_met3_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_met4_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_zl_met1_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_met2_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_met3_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_met4_nb1" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_met1_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_met2_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_met3_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_nomctrigcorr_met4_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_ht_zl_met1_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_met2_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_met3_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_ht_zl_met4_nb0" ) ; sprintf( histxtitlelist[hind], "HT (GeV)" ) ; hind++ ;



       sprintf( histnamelist[hind], "h_nb_sl_all" ) ; sprintf( histxtitlelist[hind], "N btags" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_nb_ldp_all" ) ; sprintf( histxtitlelist[hind], "N btags" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_nb_zl_all" ) ; sprintf( histxtitlelist[hind], "N btags" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_njets_sl_all" ) ; sprintf( histxtitlelist[hind], "N jets" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_njets_sl_nb0" ) ; sprintf( histxtitlelist[hind], "N jets" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_njets_sl_nb1" ) ; sprintf( histxtitlelist[hind], "N jets" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_njets_sl_nb2" ) ; sprintf( histxtitlelist[hind], "N jets" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_njets_sl_nb3" ) ; sprintf( histxtitlelist[hind], "N jets" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_njets_ldp_all" ) ; sprintf( histxtitlelist[hind], "N jets" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_njets_ldp_nb0" ) ; sprintf( histxtitlelist[hind], "N jets" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_njets_ldp_nb1" ) ; sprintf( histxtitlelist[hind], "N jets" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_njets_ldp_nb2" ) ; sprintf( histxtitlelist[hind], "N jets" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_njets_ldp_nb3" ) ; sprintf( histxtitlelist[hind], "N jets" ) ; hind++ ;





       sprintf( histnamelist[hind], "h_mt_sl_all" ) ; sprintf( histxtitlelist[hind], "MT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mt_sl_nb0" ) ; sprintf( histxtitlelist[hind], "MT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mt_sl_nb1" ) ; sprintf( histxtitlelist[hind], "MT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mt_sl_nb2" ) ; sprintf( histxtitlelist[hind], "MT (GeV)" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mt_sl_nb3" ) ; sprintf( histxtitlelist[hind], "MT (GeV)" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_npv_sl_nb0"  ) ; sprintf( histxtitlelist[hind], "N PV" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_npv_ldp_nb0" ) ; sprintf( histxtitlelist[hind], "N PV" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_npv_zl_nb0"  ) ; sprintf( histxtitlelist[hind], "N PV" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_npv_sl_nb1"  ) ; sprintf( histxtitlelist[hind], "N PV" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_npv_ldp_nb1" ) ; sprintf( histxtitlelist[hind], "N PV" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_npv_zl_nb1"  ) ; sprintf( histxtitlelist[hind], "N PV" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_mindphi_0lep_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_met1_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_met2_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_met3_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_met4_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_ht1_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_ht2_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_ht3_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_ht4_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_mindphi_0lep_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_met1_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_met2_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_met3_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_met4_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_ht1_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_ht2_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_ht3_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphi_0lep_ht4_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhi" ) ; hind++ ;


       sprintf( histnamelist[hind], "h_mindphin_0lep_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_met1_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_met2_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_met3_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_met4_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_ht1_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_ht2_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_ht3_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_ht4_nb0" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_mindphin_0lep_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_met1_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_met2_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_met3_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_met4_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_ht1_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_ht2_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_ht3_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_0lep_ht4_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_mindphin_1lep_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_1lep_met1_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_1lep_met2_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_1lep_met3_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_1lep_met4_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_1lep_ht1_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_1lep_ht2_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_1lep_ht3_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_mindphin_1lep_ht4_nb1" ) ; sprintf( histxtitlelist[hind], "minDeltaPhiN" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_pfOcaloMET_ldp" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_ldp_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_ldp_nb0" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_ldp_nb0_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_ldp_nb1" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_ldp_nb1_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_ldp_nb2" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_ldp_nb2_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_ldp_nb3" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_ldp_nb3_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_pfOcaloMET_sl" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_sl_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_sl_nb0" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_sl_nb0_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_sl_nb1" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_sl_nb1_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_sl_nb2" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_sl_nb2_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_sl_nb3" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_sl_nb3_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;

       sprintf( histnamelist[hind], "h_pfOcaloMET_zl_nb0" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_zl_nb0_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_zl_nb1" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;
       sprintf( histnamelist[hind], "h_pfOcaloMET_zl_nb1_wide" ) ; sprintf( histxtitlelist[hind], "pfMET/caloMET" ) ; hind++ ;



       int nhistlist( hind ) ;



     //---------

       char histnamelist_lhb[200][100] ;
       int hind_lhb(0) ;


       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_sl_nb0" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_sl_nb1" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_sl_nb2" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_sl_nb3" ) ; hind_lhb++ ;

       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_ldp_nb0" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_ldp_nb1" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_ldp_nb2" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_ldp_nb3" ) ; hind_lhb++ ;

       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_sl_nomctrigcorr_nb0" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_sl_nomctrigcorr_nb1" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_sl_nomctrigcorr_nb2" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_sl_nomctrigcorr_nb3" ) ; hind_lhb++ ;

       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_ldp_nomctrigcorr_nb0" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_ldp_nomctrigcorr_nb1" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_ldp_nomctrigcorr_nb2" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_ldp_nomctrigcorr_nb3" ) ; hind_lhb++ ;

       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_zl_nb0" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_zl_nb1" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_zl_nb2" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_zl_nb3" ) ; hind_lhb++ ;

       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_zl_nomctrigcorr_nb0" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_zl_nomctrigcorr_nb1" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_zl_nomctrigcorr_nb2" ) ; hind_lhb++ ;
       sprintf( histnamelist_lhb[hind_lhb], "h_lhb_zl_nomctrigcorr_nb3" ) ; hind_lhb++ ;

       int nhistlist_lhb( hind_lhb ) ;


     //---------


       loadHist( infile ) ;

       kfactorAlreadyApplied = false ;

     //--- log y scale.
       printf("\n\n Log scale plots.\n\n") ;

       gStyle->SetOptLogy(1) ;

       for ( int hi=0; hi<nhistlist; hi++ ) {
          if ( strlen(histnamepat) > 0 ) {
             TString str( histnamelist[hi] ) ;
             if ( str.Contains( histnamepat ) ) {
                drawSet( histnamelist[hi], histxtitlelist[hi] ) ;
             }
          } else {
             drawSet( histnamelist[hi], histxtitlelist[hi] ) ;
          }
       }

       for ( int hi=0; hi<nhistlist_lhb; hi++ ) {
          if ( strlen(histnamepat) > 0 ) {
             TString str( histnamelist_lhb[hi] ) ;
             if ( str.Contains( histnamepat ) ) {
                drawSetLHB( histnamelist_lhb[hi] ) ;
             }
          } else {
             drawSetLHB( histnamelist_lhb[hi] ) ;
          }
       }





       kfactorAlreadyApplied = true ;

     //--- linear y scale.
       printf("\n\n Linear scale plots.\n\n") ;

       gStyle->SetOptLogy(0) ;

       for ( int hi=0; hi<nhistlist; hi++ ) {
          if ( strlen(histnamepat) > 0 ) {
             TString str( histnamelist[hi] ) ;
             if ( str.Contains( histnamepat ) ) {
                drawSet( histnamelist[hi], histxtitlelist[hi] ) ;
             }
          } else {
             drawSet( histnamelist[hi], histxtitlelist[hi] ) ;
          }
       }

       for ( int hi=0; hi<nhistlist_lhb; hi++ ) {
          if ( strlen(histnamepat) > 0 ) {
             TString str( histnamelist_lhb[hi] ) ;
             if ( str.Contains( histnamepat ) ) {
                drawSetLHB( histnamelist_lhb[hi] ) ;
             }
          } else {
             drawSetLHB( histnamelist_lhb[hi] ) ;
          }
       }


       printf("\n\n Done.\n\n") ;


    } // dmc_drawplots

  //--------------------------------------------------------


   void drawSetLHB( const char* hname_base ) {

      printf(" drawSetLHB : %s\n", hname_base ) ;

      flattenLHB( hname_base ) ;

      bool islogy = gStyle->GetOptLogy() ;

      char cname[1000] ;
      if ( islogy ) {
         if ( recycleCanvas ) {
            sprintf( cname, "can_logy" ) ;
         } else {
            sprintf( cname, "can_logy_%s", hname_base ) ;
         }
      } else {
         if ( recycleCanvas ) {
            sprintf( cname, "can" ) ;
         } else {
            sprintf( cname, "can_%s", hname_base ) ;
         }
      }
      TCanvas* dmccan = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( dmccan == 0x0 ) {
         dmccan = new TCanvas( cname, hname_base, 600, 750 ) ;
      }
      dmccan->Clear() ;

      char hname[1000] ;
      sprintf( hname, "%s_%s_flat", hname_base, "data" ) ;
      TH1F* hdata = (TH1F*) gDirectory->FindObject( hname ) ;
      if ( hdata == 0x0 ) {
         printf("\n\n *** drawSetLHB: can't find data hist with name %s\n\n", hname ) ; return ;
      }

      sprintf( hname, "%s_mcstack_flat", hname_base ) ;
      THStack* hmcstack = new THStack() ;

      sprintf( hname, "%s_mcsum_flat", hname_base ) ;
      TH1F* hmcsum = (TH1F*) hdata->Clone( hname ) ;
      hmcsum -> Reset() ;

      hdata -> SetLineWidth(2) ;
      hdata -> SetMarkerStyle(20) ;

      TLegend* legend = new TLegend( 0.80, 0.67, 0.95, 0.92 ) ;
      legend->SetFillColor(kWhite) ;

      for ( int ci=1; ci<nComps; ci++ ) {

         sprintf( hname, "%s_%s_flat", hname_base, compname[ci] ) ;
         TH1F* hmc = (TH1F*) gDirectory->FindObject( hname ) ;
         if ( !kfactorAlreadyApplied ) { hmc->Scale( kfactor[ci] ) ; }
         if ( hmc == 0x0 ) { printf("\n\n *** drawSet: missing MC hist %s\n", hname ) ; return ; }
         hmcsum -> Add( hmc ) ;
         hmcstack -> Add( hmc ) ;
         legend -> AddEntry( hmc, compname[ci] ) ;

      }

      sprintf( hname, "%s_ratio_flat", hname_base ) ;
      TH1F* hratio = (TH1F*) hdata->Clone( hname ) ;
      hratio->Reset() ;

      for ( int bi=1; bi<=hdata->GetNbinsX(); bi++ ) {
         double data = hdata->GetBinContent(bi) ;
         double data_err = hdata->GetBinError(bi) ;
         double mc = hmcsum->GetBinContent(bi) ;
         double mc_err = hmcsum->GetBinError(bi) ;
         double val = 0. ;
         double err = 0. ;
         if ( mc > 0 ) {
            val = data / mc ;
            double errsq(0.) ;
            if ( mc > 0 ) { errsq += pow(mc_err/mc,2) ; }
            if ( data > 0 ) { errsq += pow(data_err/data,2) ; }
            err = sqrt(errsq) ;
         }
         hratio->SetBinContent(bi,val) ;
         hratio->SetBinError(bi,err) ;
      } // bi.

      hratio->SetMinimum(0.0) ;
      hratio->SetMaximum(2.0) ;




      double hmax = hdata->GetBinContent( hdata->GetMaximumBin() ) ;
      if ( hmcsum->GetBinContent( hmcsum->GetMaximumBin() ) > hmax ) { hmax = hmcsum->GetBinContent( hmcsum->GetMaximumBin() ) ; }
      if ( islogy ) {
         hmax = 3*hmax ;
      } else {
         hmax = 1.2*hmax ;
      }
      hdata->SetMaximum(hmax) ;

      char padname[1000] ;
      sprintf( padname, "tp_%s_flat", hname_base ) ;
      TPad* toppad = new TPad( padname, padname, 0.02, 0.3, 0.98, 0.98 ) ;
      toppad->Draw() ;
      sprintf( padname, "bp_%s_flat", hname_base ) ;
      TPad* bottompad = new TPad( padname, padname, 0.02, 0.02, 0.98, 0.28 ) ;
      bottompad->Draw() ;

      hmcsum->SetMarkerStyle(0) ;

      toppad->cd() ;
      gStyle->SetOptTitle(0) ;
      hdata->DrawCopy() ;
      hmcstack->Draw("samehist") ;
      hmcsum->Draw("esame") ;
      hdata->DrawCopy("same") ;
      legend->Draw() ;
      toppad->Update() ;
      toppad->Draw() ;
      TText* title = new TText() ;
      title->SetTextSize(0.040) ;
      title->DrawTextNDC( 0.05, 0.95, hdata->GetTitle() ) ;

      dmccan->Update() ;
      dmccan->Draw() ;


      bottompad->cd() ;
      hratio->UseCurrentStyle() ;
      gPad->SetLogy(0) ;
      gStyle->SetNdivisions(404,"y") ;
      gStyle->SetOptTitle(0) ;
      gStyle->SetLabelSize(0.10,"x") ;
      gStyle->SetLabelSize(0.10,"y") ;
      gStyle->SetTitleSize(0.11,"y") ;
      gStyle->SetTitleOffset(0.4,"y") ;
      hratio->SetLineWidth(2) ;
      hratio->SetMarkerStyle(20) ;
      hratio->Draw() ;
      hratio->SetYTitle("Data/MC") ;
      TLine* line = new TLine() ;
      line->SetLineStyle(2) ;
      double xl = hdata->GetBinLowEdge(1) ;
      double xh = hdata->GetBinLowEdge( hdata->GetNbinsX() ) + hdata->GetBinWidth(1) ;
      line->DrawLine( xl, 1., xh, 1. ) ;

      dmccan->Update() ;
      dmccan->Draw() ;


      if ( savePdf ) {
         TString dataset( inrootfile ) ;
         dataset.ReplaceAll("rootfiles/dmc_plots2_","") ;
         dataset.ReplaceAll(".root","") ;
         char filename[10000] ;
         if ( islogy ) {
            sprintf( filename, "outputfiles/%s_flat_%s_logy.pdf", hname_base, dataset.Data() ) ;
         } else {
            sprintf( filename, "outputfiles/%s_flat_%s.pdf", hname_base, dataset.Data() ) ;
         }
         dmccan->SaveAs( filename ) ;
      }

   } // drawSetLHB

   //=======================================================================================

   void drawSet( const char* hname_base, const char* xtitle ) {

      bool isMetPlot(false) ;
      bool isHtPlot(false) ;
      TString hnstr( hname_base ) ;

      if ( hnstr.Contains("h_met") ) {
         isMetPlot = true ;
      }
      if ( hnstr.Contains("h_ht") ) {
         isHtPlot = true ;
      }

      printf(" drawSet : %s\n", hname_base ) ;

      bool islogy = gStyle->GetOptLogy() ;

      char cname[1000] ;
      if ( islogy ) {
         if ( recycleCanvas ) {
            sprintf( cname, "can_logy" ) ;
         } else {
            sprintf( cname, "can_logy_%s", hname_base ) ;
         }
      } else {
         if ( recycleCanvas ) {
            sprintf( cname, "can" ) ;
         } else {
            sprintf( cname, "can_%s", hname_base ) ;
         }
      }
      TCanvas* dmccan = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( dmccan == 0x0 ) {
         dmccan = new TCanvas( cname, hname_base, 600, 750 ) ;
      }
      dmccan->Clear() ;

      char hname[1000] ;
      sprintf( hname, "%s_%s", hname_base, "data" ) ;
      TH1F* hdata = (TH1F*) gDirectory->FindObject( hname ) ;
      if ( hdata == 0x0 ) {
         printf("\n\n *** drawSet: can't find data hist with name %s\n\n", hname ) ; return ;
      }

      sprintf( hname, "%s_mcstack", hname_base ) ;
      THStack* hmcstack = new THStack() ;

      sprintf( hname, "%s_mcsum", hname_base ) ;
      TH1F* hmcsum = (TH1F*) hdata->Clone( hname ) ;
      hmcsum -> Reset() ;

      hdata -> SetLineWidth(2) ;
      hdata -> SetMarkerStyle(20) ;

      TLegend* legend = new TLegend( 0.80, 0.67, 0.95, 0.92 ) ;
      legend->SetFillColor(kWhite) ;

      for ( int ci=1; ci<nComps; ci++ ) {

         sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
         TH1F* hmc = (TH1F*) gDirectory->FindObject( hname ) ;
         if ( !kfactorAlreadyApplied ) { hmc->Scale( kfactor[ci] ) ; }
         if ( hmc == 0x0 ) { printf("\n\n *** drawSet: missing MC hist %s\n", hname ) ; return ; }
         hmcsum -> Add( hmc ) ;
         hmcstack -> Add( hmc ) ;
         legend -> AddEntry( hmc, compname[ci] ) ;

      }

      sprintf( hname, "%s_diff", hname_base ) ;
      TH1F* hratio = (TH1F*) hdata->Clone( hname ) ;
      hratio->Reset() ;

      for ( int bi=1; bi<=hdata->GetNbinsX(); bi++ ) {
         double data = hdata->GetBinContent(bi) ;
         double data_err = hdata->GetBinError(bi) ;
         double mc = hmcsum->GetBinContent(bi) ;
         double mc_err = hmcsum->GetBinError(bi) ;
         double val = 0. ;
         double err = 0. ;
         if ( mc > 0 ) {
            val = data / mc ;
            double errsq(0.) ;
            if ( mc > 0 ) { errsq += pow(mc_err/mc,2) ; }
            if ( data > 0 ) { errsq += pow(data_err/data,2) ; }
            err = sqrt(errsq) ;
         }
         hratio->SetBinContent(bi,val) ;
         hratio->SetBinError(bi,err) ;
      } // bi.

      hratio->SetMinimum(0.0) ;
      hratio->SetMaximum(2.0) ;




      double hmax = hdata->GetBinContent( hdata->GetMaximumBin() ) ;
      if ( hmcsum->GetBinContent( hmcsum->GetMaximumBin() ) > hmax ) { hmax = hmcsum->GetBinContent( hmcsum->GetMaximumBin() ) ; }
      if ( islogy ) {
         hmax = 3*hmax ;
      } else {
         hmax = 1.2*hmax ;
      }
      hdata->SetMaximum(hmax) ;
      hdata->SetXTitle( xtitle ) ;

      char padname[1000] ;
      sprintf( padname, "tp_%s", hname_base ) ;
      TPad* toppad = new TPad( padname, padname, 0.02, 0.3, 0.98, 0.98 ) ;
      toppad->Draw() ;
      sprintf( padname, "bp_%s", hname_base ) ;
      TPad* bottompad = new TPad( padname, padname, 0.02, 0.02, 0.98, 0.28 ) ;
      bottompad->Draw() ;

      hmcsum->SetMarkerStyle(0) ;

      toppad->cd() ;
      gStyle->SetOptTitle(0) ;
      hdata->DrawCopy() ;
      hmcstack->Draw("samehist") ;
      hmcsum->Draw("esame") ;
      hdata->DrawCopy("same") ;
      legend->Draw() ;
      toppad->Update() ;
      toppad->Draw() ;
      TText* title = new TText() ;
      title->SetTextSize(0.040) ;
      title->DrawTextNDC( 0.05, 0.95, hdata->GetTitle() ) ;
      if ( isMetPlot ) { drawMetBinEdges( hdata ) ; }
      if ( isHtPlot )  { drawHtBinEdges( hdata ) ; }


      dmccan->Update() ;
      dmccan->Draw() ;


      bottompad->cd() ;
      hratio->UseCurrentStyle() ;
      gPad->SetLogy(0) ;
      gStyle->SetNdivisions(404,"y") ;
      gStyle->SetOptTitle(0) ;
      gStyle->SetLabelSize(0.10,"x") ;
      gStyle->SetLabelSize(0.10,"y") ;
      gStyle->SetTitleSize(0.11,"y") ;
      gStyle->SetTitleOffset(0.4,"y") ;
      hratio->SetLineWidth(2) ;
      hratio->SetMarkerStyle(20) ;
      hratio->Draw() ;
      hratio->SetYTitle("Data/MC") ;
      TLine* line = new TLine() ;
      line->SetLineStyle(2) ;
      double xl = hdata->GetBinLowEdge(1) ;
      double xh = hdata->GetBinLowEdge( hdata->GetNbinsX() ) + hdata->GetBinWidth(1) ;
      line->DrawLine( xl, 1., xh, 1. ) ;
      if ( isMetPlot ) { drawMetBinEdges( hratio ) ; }
      if ( isHtPlot )  { drawHtBinEdges( hratio ) ; }

      dmccan->Update() ;
      dmccan->Draw() ;


      if ( savePdf ) {
         TString dataset( inrootfile ) ;
         dataset.ReplaceAll("rootfiles/dmc_plots2_","") ;
         dataset.ReplaceAll(".root","") ;
         char filename[10000] ;
         if ( islogy ) {
            sprintf( filename, "outputfiles/%s_%s_logy.pdf", hname_base, dataset.Data() ) ;
         } else {
            sprintf( filename, "outputfiles/%s_%s.pdf", hname_base, dataset.Data() ) ;
         }
         dmccan->SaveAs( filename ) ;
      }

   } // drawSet

   //=======================================================================================


   void flattenLHB( const char* hname_base ) {

      ///// int nBinsMET(0), nBinsHT(0), nflatbins(0) ;
      int nflatbins(0) ;

      for ( int ci=0; ci<nComps; ci++ ) {

         TH2F* hp(0x0) ;
         char hname[1000] ;

         sprintf( hname, "%s_%s_flat", hname_base, compname[ci] ) ;
         TH1F* checkhp = (TH1F*) gDirectory->FindObject( hname ) ;
         if ( checkhp != 0x0 ) { continue ; }


         sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
         hp = (TH2F*) gDirectory->FindObject( hname ) ;
         if ( hp == 0x0 ) { printf("\n\n *** drawSet: missing lhb hist %s\n", hname ) ; return ; }

         if ( ci==0 ) {
            //nBinsMET = hp -> GetNbinsX() ;
            //nBinsHT  = hp -> GetNbinsY() ;
            nflatbins = 1 + (nBinsHT+1)*nBinsMET ;
         }

         char flatname[1000] ;
         sprintf( flatname, "%s_flat", hname ) ;

         gStyle->SetLabelSize(0.05,"x") ;
         gStyle->SetLabelSize(0.05,"y") ;
         TH1F* hpf = new TH1F( flatname, hp->GetTitle(), nflatbins, 0.5, nflatbins+0.5 ) ;
         hpf->Sumw2() ;
         hpf->SetFillColor( hp->GetFillColor() ) ;

         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
               int fbi = 1 + mbi*(nBinsHT+1) + hbi + 1 ;
               hpf -> SetBinContent( fbi, hp->GetBinContent( mbi+1, hbi+1 ) ) ;
               hpf -> SetBinError(   fbi, hp->GetBinError  ( mbi+1, hbi+1 ) ) ;
               char binlabel[100] ;
               sprintf(binlabel, "M%d_H%d", mbi+1, hbi+1) ;
               hpf->GetXaxis()->SetBinLabel( fbi, binlabel ) ;
            } // hbi
         } // mbi

         hpf->SetLabelSize(0.04,"x") ;
         hpf->GetXaxis()->LabelsOption("v") ;

      } // ci.


   } // flattenLHB


   //=======================================================================================

   void drawMetBinEdges( TH1F* hp ) {

      double hmax = hp->GetMaximum() ;
      TLine* line = new TLine() ;
      line->SetLineStyle(2) ;
      for ( int mbi=1; mbi<nBinsMET; mbi++ ) {
         line->DrawLine( Mbins[mbi], 0., Mbins[mbi], hmax ) ;
      }

   } // drawMetBinEdges

   //=======================================================================================

   void drawHtBinEdges( TH1F* hp ) {

      double hmax = hp->GetMaximum() ;
      TLine* line = new TLine() ;
      line->SetLineStyle(2) ;
      for ( int mbi=1; mbi<nBinsHT; mbi++ ) {
         line->DrawLine( Hbins[mbi], 0., Hbins[mbi], hmax ) ;
      }

   } // drawHtBinEdges



   //=======================================================================================
void loadHist(const char* filename, const char* pfx, const char* pat, Bool_t doAdd, Double_t scaleFactor)
{
  cout << " Reading histograms from file: " << filename << endl << flush ;
  TFile inf(filename) ;
  //inf.ReadAll() ;
  TList* list = inf.GetListOfKeys() ;
  TIterator* iter = list->MakeIterator();

  TRegexp re(pat,kTRUE) ;
  std::cout << "pat = " << pat << std::endl ;

  gDirectory->cd("Rint:") ;

  TObject* obj ;
  TKey* key ;
  std::cout << "doAdd = " << (doAdd?"T":"F") << std::endl ;
  std::cout << "loadHist: reading." ;
  while((key=(TKey*)iter->Next())) {
   
    Int_t ridx = TString(key->GetName()).Index(re) ;    
    if (ridx==-1) {
      continue ;
    }

    obj = inf.Get(key->GetName()) ;
    TObject* clone ;
    if (pfx) {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	oldObj = gDirectory->Get(Form("%s_%s",pfx,obj->GetName())) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }
      if (oldObj) {
	clone = oldObj ;
        if ( scaleFactor > 0 ) {
           ((TH1*)clone)->Sumw2() ;
           ((TH1*)clone)->Add((TH1*)obj, scaleFactor) ;
        } else {
           ((TH1*)clone)->Add((TH1*)obj) ;
        }
      } else {
	clone = obj->Clone(Form("%s_%s",pfx,obj->GetName())) ;
      }


    } else {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	oldObj = gDirectory->Get(key->GetName()) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }

      if (oldObj) {
	clone = oldObj ;
        if ( scaleFactor > 0 ) {
           ((TH1*)clone)->Sumw2() ;
           ((TH1*)clone)->Add((TH1*)obj, scaleFactor) ;
        } else {
           ((TH1*)clone)->Add((TH1*)obj) ;
        }
      } else {
	clone = obj->Clone() ;
      }
    }
    if ( scaleFactor > 0 && !doAdd ) {
       ((TH1*) clone)->Sumw2() ;
       ((TH1*) clone)->Scale(scaleFactor) ;
    }
    if (!gDirectory->GetList()->FindObject(clone)) {
      gDirectory->Append(clone) ;
    }
    std::cout << "." ;
    std::cout.flush() ;
  }
  std::cout << std::endl;
  inf.Close() ;
  delete iter ;
}

//==========================================================================================
