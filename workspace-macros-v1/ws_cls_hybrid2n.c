
#include "TRandom2.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooStats/ModelConfig.h"

#include <string.h>

  using namespace RooFit;
  using namespace RooStats;

  //==============================================================================================

    //--- prototypes

    void printObservables() ;

    void saveObservables() ;

    void resetObservables() ;

    void generateObservables() ;

    void initializeFitpars() ;

  //==============================================================================================

    //--- Global variables.


       //--- pointers to the observables.

       RooRealVar* rrv_Nsig_1b(0x0) ;
       RooRealVar* rrv_Nsb_1b(0x0) ;
       RooRealVar* rrv_Nsig_sl_1b(0x0) ;
       RooRealVar* rrv_Nsb_sl_1b(0x0) ;
       RooRealVar* rrv_Nsig_ldp_1b(0x0) ;
       RooRealVar* rrv_Nsb_ldp_1b(0x0) ;

       RooRealVar* rrv_Nsig_2b(0x0) ;
       RooRealVar* rrv_Nsb_2b(0x0) ;
       RooRealVar* rrv_Nsig_sl_2b(0x0) ;
       RooRealVar* rrv_Nsb_sl_2b(0x0) ;
       RooRealVar* rrv_Nsig_ldp_2b(0x0) ;
       RooRealVar* rrv_Nsb_ldp_2b(0x0) ;

       RooRealVar* rrv_Nsig_3b(0x0) ;
       RooRealVar* rrv_Nsb_3b(0x0) ;
       RooRealVar* rrv_Nsig_sl_3b(0x0) ;
       RooRealVar* rrv_Nsb_sl_3b(0x0) ;
       RooRealVar* rrv_Nsig_ldp_3b(0x0) ;
       RooRealVar* rrv_Nsb_ldp_3b(0x0) ;


       RooRealVar* rrv_Nsig_ee(0x0) ;
       RooRealVar* rrv_Nsb_ee(0x0) ;
       RooRealVar* rrv_Nsig_mm(0x0) ;
       RooRealVar* rrv_Nsb_mm(0x0) ;


       //--- pointers to likelihood model predictions for the observables.

       RooAbsReal* rfv_n_sig_1b(0x0)       ;
       RooAbsReal* rfv_n_sb_1b(0x0)        ;
       RooAbsReal* rfv_n_sig_sl_1b(0x0)    ;
       RooAbsReal* rfv_n_sb_sl_1b(0x0)     ;
       RooAbsReal* rfv_n_sig_ldp_1b(0x0)   ;
       RooAbsReal* rfv_n_sb_ldp_1b(0x0)    ;

       RooAbsReal* rfv_n_sig_2b(0x0)       ;
       RooAbsReal* rfv_n_sb_2b(0x0)        ;
       RooAbsReal* rfv_n_sig_sl_2b(0x0)    ;
       RooAbsReal* rfv_n_sb_sl_2b(0x0)     ;
       RooAbsReal* rfv_n_sig_ldp_2b(0x0)   ;
       RooAbsReal* rfv_n_sb_ldp_2b(0x0)    ;

       RooAbsReal* rfv_n_sig_3b(0x0)       ;
       RooAbsReal* rfv_n_sb_3b(0x0)        ;
       RooAbsReal* rfv_n_sig_sl_3b(0x0)    ;
       RooAbsReal* rfv_n_sb_sl_3b(0x0)     ;
       RooAbsReal* rfv_n_sig_ldp_3b(0x0)   ;
       RooAbsReal* rfv_n_sb_ldp_3b(0x0)    ;


       RooAbsReal* rfv_n_sig_ee(0x0)    ;
       RooAbsReal* rfv_n_sb_ee(0x0)     ;
       RooAbsReal* rfv_n_sig_mm(0x0)    ;
       RooAbsReal* rfv_n_sb_mm(0x0)     ;


       //--- values of the observables.

       int actual_Nsig_1b(0) ;
       int actual_Nsb_1b(0) ;
       int actual_Nsig_sl_1b(0) ;
       int actual_Nsb_sl_1b(0) ;
       int actual_Nsig_ldp_1b(0) ;
       int actual_Nsb_ldp_1b(0) ;

       int actual_Nsig_2b(0) ;
       int actual_Nsb_2b(0) ;
       int actual_Nsig_sl_2b(0) ;
       int actual_Nsb_sl_2b(0) ;
       int actual_Nsig_ldp_2b(0) ;
       int actual_Nsb_ldp_2b(0) ;

       int actual_Nsig_3b(0) ;
       int actual_Nsb_3b(0) ;
       int actual_Nsig_sl_3b(0) ;
       int actual_Nsb_sl_3b(0) ;
       int actual_Nsig_ldp_3b(0) ;
       int actual_Nsb_ldp_3b(0) ;

       int actual_Nsig_ee(0) ;
       int actual_Nsb_ee(0) ;
       int actual_Nsig_mm(0) ;
       int actual_Nsb_mm(0) ;




       //--- pointers to fit parameters (so that I can initialize them).

         RooRealVar* fv_mu_susy_sig_1b(0x0) ;
         RooRealVar* fv_mu_ttwj_sb_1b(0x0) ;

         RooRealVar* fv_mu_qcd_sb_ldp_1b(0x0) ;
         RooRealVar* fv_mu_qcd_sig_ldp_1b(0x0) ;
         RooRealVar* fv_mu_ttwj_sb_sl_1b(0x0) ;
         RooRealVar* fv_mu_ttwj_sig_sl_1b(0x0) ;

         RooRealVar* fv_mu_qcd_sb_ldp_2b(0x0) ;
         RooRealVar* fv_mu_qcd_sig_ldp_2b(0x0) ;
         RooRealVar* fv_mu_ttwj_sb_sl_2b(0x0) ;
         RooRealVar* fv_mu_ttwj_sig_sl_2b(0x0) ;

         RooRealVar* fv_mu_qcd_sb_ldp_3b(0x0) ;
         RooRealVar* fv_mu_qcd_sig_ldp_3b(0x0) ;
         RooRealVar* fv_mu_ttwj_sb_sl_3b(0x0) ;
         RooRealVar* fv_mu_ttwj_sig_sl_3b(0x0) ;

         RooRealVar* fv_mu_znn_sig_1b(0x0) ;
         RooRealVar* fv_mu_znn_sb_1b(0x0) ;





       //--- pointers to parameters (needed for computing initialization values based on given observables).

         RooRealVar* znnoverll_bfratio(0x0) ;
         RooRealVar* dataoverll_lumiratio(0x0) ;



        //--- all below are functions.

         RooAbsReal* Rlsb_passfail(0x0) ;



         //-- Zll, common
         RooAbsReal* knn_sig_1b(0x0) ;
         RooAbsReal* knn_sig_2b(0x0) ;
         RooAbsReal* knn_sig_3b(0x0) ;
         RooAbsReal* knn_sb_1b(0x0) ;
         RooAbsReal* knn_sb_2b(0x0) ;
         RooAbsReal* knn_sb_3b(0x0) ;

         //-- Zee
         RooAbsReal* acc_ee_sig(0x0) ;
         RooAbsReal* acc_ee_sb(0x0) ;
         RooAbsReal* fsig_ee(0x0) ;
         RooAbsReal* eff_ee(0x0) ;
         RooAbsReal* sf_ee(0x0) ;

         //-- Zmm
         RooAbsReal* acc_mm_sig(0x0) ;
         RooAbsReal* acc_mm_sb(0x0) ;
         RooAbsReal* fsig_mm(0x0) ;
         RooAbsReal* eff_mm(0x0) ;
         RooAbsReal* sf_mm(0x0) ;



         //-- LDP MC inputs
         RooAbsReal* mu_ttwj_sig_ldp_1b(0x0) ;
         RooAbsReal* mu_ttwj_sb_ldp_1b(0x0) ;
         RooAbsReal* mu_znn_sig_ldp_1b(0x0) ;
         RooAbsReal* mu_znn_sb_ldp_1b(0x0) ;

         RooAbsReal* mu_ttwj_sig_ldp_2b(0x0) ;
         RooAbsReal* mu_ttwj_sb_ldp_2b(0x0) ;
         RooAbsReal* mu_znn_sig_ldp_2b(0x0) ;
         RooAbsReal* mu_znn_sb_ldp_2b(0x0) ;

         RooAbsReal* mu_ttwj_sig_ldp_3b(0x0) ;
         RooAbsReal* mu_ttwj_sb_ldp_3b(0x0) ;
         RooAbsReal* mu_znn_sig_ldp_3b(0x0) ;
         RooAbsReal* mu_znn_sb_ldp_3b(0x0) ;


         //-- Efficiency parameters
         RooAbsReal* eff_sf_sig_1b(0x0) ;
         RooAbsReal* eff_sf_sb_1b(0x0) ;
         RooAbsReal* eff_sf_sig_sl_1b(0x0) ;
         RooAbsReal* eff_sf_sb_sl_1b(0x0) ;
         RooAbsReal* eff_sf_sig_ldp_1b(0x0) ;
         RooAbsReal* eff_sf_sb_ldp_1b(0x0) ;

         RooAbsReal* eff_sf_sig_2b(0x0) ;
         RooAbsReal* eff_sf_sb_2b(0x0) ;
         RooAbsReal* eff_sf_sig_sl_2b(0x0) ;
         RooAbsReal* eff_sf_sb_sl_2b(0x0) ;
         RooAbsReal* eff_sf_sig_ldp_2b(0x0) ;
         RooAbsReal* eff_sf_sb_ldp_2b(0x0) ;

         RooAbsReal* eff_sf_sig_3b(0x0) ;
         RooAbsReal* eff_sf_sb_3b(0x0) ;
         RooAbsReal* eff_sf_sig_sl_3b(0x0) ;
         RooAbsReal* eff_sf_sb_sl_3b(0x0) ;
         RooAbsReal* eff_sf_sig_ldp_3b(0x0) ;
         RooAbsReal* eff_sf_sb_ldp_3b(0x0) ;


         RooAbsReal* btageff_sf_sig_1b(0x0) ;
         RooAbsReal* btageff_sf_sb_1b(0x0) ;
         RooAbsReal* btageff_sf_sig_sl_1b(0x0) ;
         RooAbsReal* btageff_sf_sb_sl_1b(0x0) ;
         RooAbsReal* btageff_sf_sig_ldp_1b(0x0) ;
         RooAbsReal* btageff_sf_sb_ldp_1b(0x0) ;

         RooAbsReal* btageff_sf_sig_2b(0x0) ;
         RooAbsReal* btageff_sf_sb_2b(0x0) ;
         RooAbsReal* btageff_sf_sig_sl_2b(0x0) ;
         RooAbsReal* btageff_sf_sb_sl_2b(0x0) ;
         RooAbsReal* btageff_sf_sig_ldp_2b(0x0) ;
         RooAbsReal* btageff_sf_sb_ldp_2b(0x0) ;

         RooAbsReal* btageff_sf_sig_3b(0x0) ;
         RooAbsReal* btageff_sf_sb_3b(0x0) ;
         RooAbsReal* btageff_sf_sig_sl_3b(0x0) ;
         RooAbsReal* btageff_sf_sb_sl_3b(0x0) ;
         RooAbsReal* btageff_sf_sig_ldp_3b(0x0) ;
         RooAbsReal* btageff_sf_sb_ldp_3b(0x0) ;


         //-- SUSY yields
         RooAbsReal* mu_susy_sb_1b(0x0) ;
         RooAbsReal* mu_susy_sig_sl_1b(0x0) ;
         RooAbsReal* mu_susy_sb_sl_1b(0x0) ;
         RooAbsReal* mu_susy_sig_ldp_1b(0x0) ;
         RooAbsReal* mu_susy_sb_ldp_1b(0x0) ;

         RooAbsReal* mu_susy_sig_2b(0x0) ;
         RooAbsReal* mu_susy_sb_2b(0x0) ;
         RooAbsReal* mu_susy_sig_sl_2b(0x0) ;
         RooAbsReal* mu_susy_sb_sl_2b(0x0) ;
         RooAbsReal* mu_susy_sig_ldp_2b(0x0) ;
         RooAbsReal* mu_susy_sb_ldp_2b(0x0) ;

         RooAbsReal* mu_susy_sig_3b(0x0) ;
         RooAbsReal* mu_susy_sb_3b(0x0) ;
         RooAbsReal* mu_susy_sig_sl_3b(0x0) ;
         RooAbsReal* mu_susy_sb_sl_3b(0x0) ;
         RooAbsReal* mu_susy_sig_ldp_3b(0x0) ;
         RooAbsReal* mu_susy_sb_ldp_3b(0x0) ;



         //-- Systematics
         RooAbsReal* sf_mc(0x0) ;
         RooAbsReal* sf_qcd_sig_1b(0x0) ;
         RooAbsReal* sf_qcd_sig_2b(0x0) ;
         RooAbsReal* sf_qcd_sig_3b(0x0) ;
         RooAbsReal* sf_qcd_sb_1b(0x0) ;
         RooAbsReal* sf_qcd_sb_2b(0x0) ;
         RooAbsReal* sf_qcd_sb_3b(0x0) ;
         RooAbsReal* sf_ttwj_sig_1b(0x0) ;
         RooAbsReal* sf_ttwj_sig_2b(0x0) ;
         RooAbsReal* sf_ttwj_sig_3b(0x0) ;
         RooAbsReal* sf_ttwj_sb_2b(0x0) ;
         RooAbsReal* sf_ttwj_sb_3b(0x0) ;


         RooAbsReal* mu_qcd_sb_1b(0x0) ;

       //--- nuisance parameter initial values.

       int np_count(0) ;
       double np_initial_val[1000] ;




       TRandom2* random_ng ;




  //==============================================================================================

   void ws_cls_hybrid2n( const char* wsfile = "output-files/expected-ws-lm9-2BL.root", bool isBgonlyStudy=false, double poiVal = 150.0, int nToys=100, bool makeTtree=true, int verbLevel=0 ) {



       TTree* toytt(0x0) ;
       TFile* ttfile(0x0) ;

       int    tt_gen_Nsig_1b ;
       int    tt_gen_Nsb_1b ;
       int    tt_gen_Nsig_sl_1b ;
       int    tt_gen_Nsb_sl_1b ;
       int    tt_gen_Nsig_ldp_1b ;
       int    tt_gen_Nsb_ldp_1b ;

       int    tt_gen_Nsig_2b ;
       int    tt_gen_Nsb_2b ;
       int    tt_gen_Nsig_sl_2b ;
       int    tt_gen_Nsb_sl_2b ;
       int    tt_gen_Nsig_ldp_2b ;
       int    tt_gen_Nsb_ldp_2b ;

       int    tt_gen_Nsig_3b ;
       int    tt_gen_Nsb_3b ;
       int    tt_gen_Nsig_sl_3b ;
       int    tt_gen_Nsb_sl_3b ;
       int    tt_gen_Nsig_ldp_3b ;
       int    tt_gen_Nsb_ldp_3b ;

       int    tt_gen_Nsig_ee ;
       int    tt_gen_Nsb_ee ;
       int    tt_gen_Nsig_mm ;
       int    tt_gen_Nsb_mm ;

       double tt_testStat ;
       double tt_dataTestStat ;
       double tt_hypo_mu_susy_sig_1b ;

       char ttname[1000] ;
       char tttitle[1000] ;

       if ( makeTtree ) {

          ttfile = gDirectory->GetFile() ;
          if ( ttfile == 0x0 ) { printf("\n\n\n *** asked for a ttree but no open file???\n\n") ; return ; }


          if ( isBgonlyStudy ) {
             sprintf( ttname, "toytt_%.0f_bgo", poiVal ) ;
             sprintf( tttitle, "Toy study for background only, mu_susy_sig_1b = %.0f", poiVal ) ;
          } else {
             sprintf( ttname, "toytt_%.0f_spb", poiVal ) ;
             sprintf( tttitle, "Toy study for signal+background, mu_susy_sig_1b = %.0f", poiVal ) ;
          }

          printf("\n\n Creating TTree : %s : %s\n\n", ttname, tttitle ) ;

          gDirectory->pwd() ;
          gDirectory->ls() ;

          toytt = new TTree( ttname, tttitle ) ;

          gDirectory->ls() ;

          toytt -> Branch(  "gen_Nsig_1b"         ,       &tt_gen_Nsig_1b         ,      "gen_Nsig_1b/I"         ) ;
          toytt -> Branch(  "gen_Nsb_1b"          ,       &tt_gen_Nsb_1b          ,      "gen_Nsb_1b/I"          ) ;
          toytt -> Branch(  "gen_Nsig_sl_1b"      ,       &tt_gen_Nsig_sl_1b      ,      "gen_Nsig_sl_1b/I"      ) ;
          toytt -> Branch(  "gen_Nsb_sl_1b"       ,       &tt_gen_Nsb_sl_1b       ,      "gen_Nsb_sl_1b/I"       ) ;
          toytt -> Branch(  "gen_Nsig_ldp_1b"     ,       &tt_gen_Nsig_ldp_1b     ,      "gen_Nsig_ldp_1b/I"     ) ;
          toytt -> Branch(  "gen_Nsb_ldp_1b"      ,       &tt_gen_Nsb_ldp_1b      ,      "gen_Nsb_ldp_1b/I"      ) ;

          toytt -> Branch(  "gen_Nsig_2b"         ,       &tt_gen_Nsig_2b         ,      "gen_Nsig_2b/I"         ) ;
          toytt -> Branch(  "gen_Nsb_2b"          ,       &tt_gen_Nsb_2b          ,      "gen_Nsb_2b/I"          ) ;
          toytt -> Branch(  "gen_Nsig_sl_2b"      ,       &tt_gen_Nsig_sl_2b      ,      "gen_Nsig_sl_2b/I"      ) ;
          toytt -> Branch(  "gen_Nsb_sl_2b"       ,       &tt_gen_Nsb_sl_2b       ,      "gen_Nsb_sl_2b/I"       ) ;
          toytt -> Branch(  "gen_Nsig_ldp_2b"     ,       &tt_gen_Nsig_ldp_2b     ,      "gen_Nsig_ldp_2b/I"     ) ;
          toytt -> Branch(  "gen_Nsb_ldp_2b"      ,       &tt_gen_Nsb_ldp_2b      ,      "gen_Nsb_ldp_2b/I"      ) ;

          toytt -> Branch(  "gen_Nsig_3b"         ,       &tt_gen_Nsig_3b         ,      "gen_Nsig_3b/I"         ) ;
          toytt -> Branch(  "gen_Nsb_3b"          ,       &tt_gen_Nsb_3b          ,      "gen_Nsb_3b/I"          ) ;
          toytt -> Branch(  "gen_Nsig_sl_3b"      ,       &tt_gen_Nsig_sl_3b      ,      "gen_Nsig_sl_3b/I"      ) ;
          toytt -> Branch(  "gen_Nsb_sl_3b"       ,       &tt_gen_Nsb_sl_3b       ,      "gen_Nsb_sl_3b/I"       ) ;
          toytt -> Branch(  "gen_Nsig_ldp_3b"     ,       &tt_gen_Nsig_ldp_3b     ,      "gen_Nsig_ldp_3b/I"     ) ;
          toytt -> Branch(  "gen_Nsb_ldp_3b"      ,       &tt_gen_Nsb_ldp_3b      ,      "gen_Nsb_ldp_3b/I"      ) ;

          toytt -> Branch(  "gen_Nsig_ee"      ,       &tt_gen_Nsig_ee      ,      "gen_Nsig_ee/I"      ) ;
          toytt -> Branch(  "gen_Nsb_ee"       ,       &tt_gen_Nsb_ee       ,      "gen_Nsb_ee/I"       ) ;
          toytt -> Branch(  "gen_Nsig_mm"      ,       &tt_gen_Nsig_mm      ,      "gen_Nsig_mm/I"      ) ;
          toytt -> Branch(  "gen_Nsb_mm"       ,       &tt_gen_Nsb_mm       ,      "gen_Nsb_mm/I"       ) ;

          toytt -> Branch(  "testStat"         ,       &tt_testStat         ,      "testStat/D"         ) ;
          toytt -> Branch(  "dataTestStat"     ,       &tt_dataTestStat     ,      "dataTestStat/D"     ) ;
          toytt -> Branch(  "hypo_mu_susy_sig_1b" ,       &tt_hypo_mu_susy_sig_1b ,      "hypo_mu_susy_sig_1b/D" ) ;

       }


     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;


       random_ng = new TRandom2(12345) ;





       TFile* wstf = new TFile( wsfile ) ;

       RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );

       ws->Print() ;






       RooDataSet* rds = (RooDataSet*) ws->obj( "ra2b_observed_rds" ) ;
       printf("\n\n\n  ===== RooDataSet ====================\n\n") ;

       rds->Print() ;
       rds->printMultiline(cout, 1, kTRUE, "") ;





       ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;
       RooAbsPdf* likelihood = modelConfig->GetPdf() ;

       const RooArgSet* nuisanceParameters = modelConfig->GetNuisanceParameters() ;

       RooRealVar* rrv_mu_susy_sig_1b = ws->var("mu_susy_sig_1b") ;
       if ( rrv_mu_susy_sig_1b == 0x0 ) {
          printf("\n\n\n *** can't find mu_susy_sig_1b in workspace.  Quitting.\n\n\n") ;
          return ;
       }












       //--- Get pointers to the model predictions of the observables.

       rfv_n_sig_1b       = ws->function("n_sig_1b") ;
       rfv_n_sb_1b        = ws->function("n_sb_1b") ;
       rfv_n_sig_sl_1b    = ws->function("n_sig_sl_1b") ;
       rfv_n_sb_sl_1b     = ws->function("n_sb_sl_1b") ;
       rfv_n_sig_ldp_1b   = ws->function("n_sig_ldp_1b") ;
       rfv_n_sb_ldp_1b    = ws->function("n_sb_ldp_1b") ;

       rfv_n_sig_2b       = ws->function("n_sig_2b") ;
       rfv_n_sb_2b        = ws->function("n_sb_2b") ;
       rfv_n_sig_sl_2b    = ws->function("n_sig_sl_2b") ;
       rfv_n_sb_sl_2b     = ws->function("n_sb_sl_2b") ;
       rfv_n_sig_ldp_2b   = ws->function("n_sig_ldp_2b") ;
       rfv_n_sb_ldp_2b    = ws->function("n_sb_ldp_2b") ;

       rfv_n_sig_3b       = ws->function("n_sig_3b") ;
       rfv_n_sb_3b        = ws->function("n_sb_3b") ;
       rfv_n_sig_sl_3b    = ws->function("n_sig_sl_3b") ;
       rfv_n_sb_sl_3b     = ws->function("n_sb_sl_3b") ;
       rfv_n_sig_ldp_3b   = ws->function("n_sig_ldp_3b") ;
       rfv_n_sb_ldp_3b    = ws->function("n_sb_ldp_3b") ;

       rfv_n_sig_ee    = ws->function("n_sig_ee") ;
       rfv_n_sb_ee     = ws->function("n_sb_ee") ;
       rfv_n_sig_mm    = ws->function("n_sig_mm") ;
       rfv_n_sb_mm     = ws->function("n_sb_mm") ;

       if (  rfv_n_sig_1b        == 0x0 ) { printf("\n\n\n *** can't find n_sig_1b       in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sb_1b         == 0x0 ) { printf("\n\n\n *** can't find n_sb_1b        in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sig_sl_1b     == 0x0 ) { printf("\n\n\n *** can't find n_sig_sl_1b    in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sb_sl_1b      == 0x0 ) { printf("\n\n\n *** can't find n_sb_sl_1b     in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sig_ldp_1b    == 0x0 ) { printf("\n\n\n *** can't find n_sig_ldp_1b   in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sb_ldp_1b     == 0x0 ) { printf("\n\n\n *** can't find n_sb_ldp_1b    in workspace.  Quitting.\n\n\n") ; return ; }

       if (  rfv_n_sig_2b        == 0x0 ) { printf("\n\n\n *** can't find n_sig_2b       in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sb_2b         == 0x0 ) { printf("\n\n\n *** can't find n_sb_2b        in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sig_sl_2b     == 0x0 ) { printf("\n\n\n *** can't find n_sig_sl_2b    in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sb_sl_2b      == 0x0 ) { printf("\n\n\n *** can't find n_sb_sl_2b     in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sig_ldp_2b    == 0x0 ) { printf("\n\n\n *** can't find n_sig_ldp_2b   in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sb_ldp_2b     == 0x0 ) { printf("\n\n\n *** can't find n_sb_ldp_2b    in workspace.  Quitting.\n\n\n") ; return ; }

       if (  rfv_n_sig_3b        == 0x0 ) { printf("\n\n\n *** can't find n_sig_3b       in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sb_3b         == 0x0 ) { printf("\n\n\n *** can't find n_sb_3b        in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sig_sl_3b     == 0x0 ) { printf("\n\n\n *** can't find n_sig_sl_3b    in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sb_sl_3b      == 0x0 ) { printf("\n\n\n *** can't find n_sb_sl_3b     in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sig_ldp_3b    == 0x0 ) { printf("\n\n\n *** can't find n_sig_ldp_3b   in workspace.  Quitting.\n\n\n") ; return ; }
       if (  rfv_n_sb_ldp_3b     == 0x0 ) { printf("\n\n\n *** can't find n_sb_ldp_3b    in workspace.  Quitting.\n\n\n") ; return ; }

       if ( rfv_n_sig_ee      == 0x0 ) { printf("\n\n\n *** can't find n_sig_ee    in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sb_ee       == 0x0 ) { printf("\n\n\n *** can't find n_sb_ee     in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sig_mm      == 0x0 ) { printf("\n\n\n *** can't find n_sig_mm    in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sb_mm       == 0x0 ) { printf("\n\n\n *** can't find n_sb_mm     in workspace.  Quitting.\n\n\n") ; return ; }


       //--- Get pointers to the fit parameters.

       fv_mu_susy_sig_1b     = ws->var("mu_susy_sig_1b"    ) ;
       fv_mu_ttwj_sb_1b      = ws->var("mu_ttwj_sb_1b"     ) ;

       fv_mu_qcd_sb_ldp_1b   = ws->var("mu_qcd_sb_ldp_1b"  ) ;
       fv_mu_qcd_sig_ldp_1b  = ws->var("mu_qcd_sig_ldp_1b" ) ;
       fv_mu_ttwj_sb_sl_1b   = ws->var("mu_ttwj_sb_sl_1b"     ) ;
       fv_mu_ttwj_sig_sl_1b  = ws->var("mu_ttwj_sig_sl_1b"    ) ;

       fv_mu_qcd_sb_ldp_2b   = ws->var("mu_qcd_sb_ldp_2b"  ) ;
       fv_mu_qcd_sig_ldp_2b  = ws->var("mu_qcd_sig_ldp_2b" ) ;
       fv_mu_ttwj_sb_sl_2b   = ws->var("mu_ttwj_sb_sl_2b"     ) ;
       fv_mu_ttwj_sig_sl_2b  = ws->var("mu_ttwj_sig_sl_2b"    ) ;

       fv_mu_qcd_sb_ldp_3b   = ws->var("mu_qcd_sb_ldp_3b"  ) ;
       fv_mu_qcd_sig_ldp_3b  = ws->var("mu_qcd_sig_ldp_3b" ) ;
       fv_mu_ttwj_sb_sl_3b   = ws->var("mu_ttwj_sb_sl_3b"     ) ;
       fv_mu_ttwj_sig_sl_3b  = ws->var("mu_ttwj_sig_sl_3b"    ) ;

       fv_mu_znn_sig_1b      = ws->var("mu_znn_sig_1b"        ) ;
       fv_mu_znn_sb_1b       = ws->var("mu_znn_sb_1b"         ) ;



       if ( fv_mu_susy_sig_1b     == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sig_1b    in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fv_mu_ttwj_sb_1b      == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sb_1b     in workspace.  Quitting.\n\n\n") ; return ; }

       if ( fv_mu_qcd_sb_ldp_1b   == 0x0 ) { printf("\n\n\n *** can't find mu_qcd_sb_ldp_1b  in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fv_mu_qcd_sig_ldp_1b  == 0x0 ) { printf("\n\n\n *** can't find mu_qcd_sig_ldp_1b in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fv_mu_ttwj_sb_sl_1b   == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sb_sl_1b  in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fv_mu_ttwj_sig_sl_1b  == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sig_sl_1b in workspace.  Quitting.\n\n\n") ; return ; }

       if ( fv_mu_qcd_sb_ldp_2b   == 0x0 ) { printf("\n\n\n *** can't find mu_qcd_sb_ldp_2b  in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fv_mu_qcd_sig_ldp_2b  == 0x0 ) { printf("\n\n\n *** can't find mu_qcd_sig_ldp_2b in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fv_mu_ttwj_sb_sl_2b   == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sb_sl_2b  in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fv_mu_ttwj_sig_sl_2b  == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sig_sl_2b in workspace.  Quitting.\n\n\n") ; return ; }

       if ( fv_mu_qcd_sb_ldp_3b   == 0x0 ) { printf("\n\n\n *** can't find mu_qcd_sb_ldp_3b  in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fv_mu_qcd_sig_ldp_3b  == 0x0 ) { printf("\n\n\n *** can't find mu_qcd_sig_ldp_3b in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fv_mu_ttwj_sb_sl_3b   == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sb_sl_3b  in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fv_mu_ttwj_sig_sl_3b  == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sig_sl_3b in workspace.  Quitting.\n\n\n") ; return ; }

       if ( fv_mu_znn_sig_1b      == 0x0 ) { printf("\n\n\n *** can't find mu_znn_sig_1b     in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fv_mu_znn_sb_1b       == 0x0 ) { printf("\n\n\n *** can't find mu_znn_sb_1b      in workspace.  Quitting.\n\n\n") ; return ; }






       //--- Get pointers to parameters.

       znnoverll_bfratio    = ws->var("znnoverll_bfratio") ;
       dataoverll_lumiratio = ws->var("dataoverll_lumiratio") ;

       if ( znnoverll_bfratio     == 0x0 ) { printf("\n\n\n *** can't find znnoverll_bfratio       in workspace.  Quitting.\n\n\n") ; return ; }
       if ( dataoverll_lumiratio  == 0x0 ) { printf("\n\n\n *** can't find dataoverll_lumiratio    in workspace.  Quitting.\n\n\n") ; return ; }

       knn_sig_1b               = ws->function("knn_sig_1b"              ) ;
       knn_sig_2b               = ws->function("knn_sig_2b"              ) ;
       knn_sig_3b               = ws->function("knn_sig_3b"              ) ;
       knn_sb_1b                = ws->function("knn_sb_1b"               ) ;
       knn_sb_2b                = ws->function("knn_sb_2b"               ) ;
       knn_sb_3b                = ws->function("knn_sb_3b"               ) ;
       acc_ee_sig               = ws->function("acc_ee_sig"              ) ;
       acc_ee_sb                = ws->function("acc_ee_sb"               ) ;
       fsig_ee                  = ws->function("fsig_ee"                 ) ;
       eff_ee                   = ws->function("eff_ee"                  ) ;
       sf_ee                    = ws->function("sf_ee"                   ) ;
       acc_mm_sig               = ws->function("acc_mm_sig"              ) ;
       acc_mm_sb                = ws->function("acc_mm_sb"               ) ;
       fsig_mm                  = ws->function("fsig_mm"                 ) ;
       eff_mm                   = ws->function("eff_mm"                  ) ;
       sf_mm                    = ws->function("sf_mm"                   ) ;
       mu_ttwj_sig_ldp_1b       = ws->function("mu_ttwj_sig_ldp_1b"      ) ;
       mu_ttwj_sb_ldp_1b        = ws->function("mu_ttwj_sb_ldp_1b"       ) ;
       mu_znn_sig_ldp_1b        = ws->function("mu_znn_sig_ldp_1b"       ) ;
       mu_znn_sb_ldp_1b         = ws->function("mu_znn_sb_ldp_1b"        ) ;
       mu_ttwj_sig_ldp_2b       = ws->function("mu_ttwj_sig_ldp_2b"      ) ;
       mu_ttwj_sb_ldp_2b        = ws->function("mu_ttwj_sb_ldp_2b"       ) ;
       mu_znn_sig_ldp_2b        = ws->function("mu_znn_sig_ldp_2b"       ) ;
       mu_znn_sb_ldp_2b         = ws->function("mu_znn_sb_ldp_2b"        ) ;
       mu_ttwj_sig_ldp_3b       = ws->function("mu_ttwj_sig_ldp_3b"      ) ;
       mu_ttwj_sb_ldp_3b        = ws->function("mu_ttwj_sb_ldp_3b"       ) ;
       mu_znn_sig_ldp_3b        = ws->function("mu_znn_sig_ldp_3b"       ) ;
       mu_znn_sb_ldp_3b         = ws->function("mu_znn_sb_ldp_3b"        ) ;
       eff_sf_sig_1b            = ws->function("eff_sf_sig_1b"           ) ;
       eff_sf_sb_1b             = ws->function("eff_sf_sb_1b"            ) ;
       eff_sf_sig_sl_1b         = ws->function("eff_sf_sig_sl_1b"        ) ;
       eff_sf_sb_sl_1b          = ws->function("eff_sf_sb_sl_1b"         ) ;
       eff_sf_sig_ldp_1b        = ws->function("eff_sf_sig_ldp_1b"       ) ;
       eff_sf_sb_ldp_1b         = ws->function("eff_sf_sb_ldp_1b"        ) ;
       eff_sf_sig_2b            = ws->function("eff_sf_sig_2b"           ) ;
       eff_sf_sb_2b             = ws->function("eff_sf_sb_2b"            ) ;
       eff_sf_sig_sl_2b         = ws->function("eff_sf_sig_sl_2b"        ) ;
       eff_sf_sb_sl_2b          = ws->function("eff_sf_sb_sl_2b"         ) ;
       eff_sf_sig_ldp_2b        = ws->function("eff_sf_sig_ldp_2b"       ) ;
       eff_sf_sb_ldp_2b         = ws->function("eff_sf_sb_ldp_2b"        ) ;
       eff_sf_sig_3b            = ws->function("eff_sf_sig_3b"           ) ;
       eff_sf_sb_3b             = ws->function("eff_sf_sb_3b"            ) ;
       eff_sf_sig_sl_3b         = ws->function("eff_sf_sig_sl_3b"        ) ;
       eff_sf_sb_sl_3b          = ws->function("eff_sf_sb_sl_3b"         ) ;
       eff_sf_sig_ldp_3b        = ws->function("eff_sf_sig_ldp_3b"       ) ;
       eff_sf_sb_ldp_3b         = ws->function("eff_sf_sb_ldp_3b"        ) ;
       btageff_sf_sig_1b        = ws->function("btageff_sf_sig_1b"       ) ;
       btageff_sf_sb_1b         = ws->function("btageff_sf_sb_1b"        ) ;
       btageff_sf_sig_sl_1b     = ws->function("btageff_sf_sig_sl_1b"    ) ;
       btageff_sf_sb_sl_1b      = ws->function("btageff_sf_sb_sl_1b"     ) ;
       btageff_sf_sig_ldp_1b    = ws->function("btageff_sf_sig_ldp_1b"   ) ;
       btageff_sf_sb_ldp_1b     = ws->function("btageff_sf_sb_ldp_1b"    ) ;
       btageff_sf_sig_2b        = ws->function("btageff_sf_sig_2b"       ) ;
       btageff_sf_sb_2b         = ws->function("btageff_sf_sb_2b"        ) ;
       btageff_sf_sig_sl_2b     = ws->function("btageff_sf_sig_sl_2b"    ) ;
       btageff_sf_sb_sl_2b      = ws->function("btageff_sf_sb_sl_2b"     ) ;
       btageff_sf_sig_ldp_2b    = ws->function("btageff_sf_sig_ldp_2b"   ) ;
       btageff_sf_sb_ldp_2b     = ws->function("btageff_sf_sb_ldp_2b"    ) ;
       btageff_sf_sig_3b        = ws->function("btageff_sf_sig_3b"       ) ;
       btageff_sf_sb_3b         = ws->function("btageff_sf_sb_3b"        ) ;
       btageff_sf_sig_sl_3b     = ws->function("btageff_sf_sig_sl_3b"    ) ;
       btageff_sf_sb_sl_3b      = ws->function("btageff_sf_sb_sl_3b"     ) ;
       btageff_sf_sig_ldp_3b    = ws->function("btageff_sf_sig_ldp_3b"   ) ;
       btageff_sf_sb_ldp_3b     = ws->function("btageff_sf_sb_ldp_3b"    ) ;
       mu_susy_sb_1b            = ws->function("mu_susy_sb_1b"           ) ;
       mu_susy_sig_sl_1b        = ws->function("mu_susy_sig_sl_1b"       ) ;
       mu_susy_sb_sl_1b         = ws->function("mu_susy_sb_sl_1b"        ) ;
       mu_susy_sig_ldp_1b       = ws->function("mu_susy_sig_ldp_1b"      ) ;
       mu_susy_sb_ldp_1b        = ws->function("mu_susy_sb_ldp_1b"       ) ;
       mu_susy_sig_2b           = ws->function("mu_susy_sig_2b"          ) ;
       mu_susy_sb_2b            = ws->function("mu_susy_sb_2b"           ) ;
       mu_susy_sig_sl_2b        = ws->function("mu_susy_sig_sl_2b"       ) ;
       mu_susy_sb_sl_2b         = ws->function("mu_susy_sb_sl_2b"        ) ;
       mu_susy_sig_ldp_2b       = ws->function("mu_susy_sig_ldp_2b"      ) ;
       mu_susy_sb_ldp_2b        = ws->function("mu_susy_sb_ldp_2b"       ) ;
       mu_susy_sig_3b           = ws->function("mu_susy_sig_3b"          ) ;
       mu_susy_sb_3b            = ws->function("mu_susy_sb_3b"           ) ;
       mu_susy_sig_sl_3b        = ws->function("mu_susy_sig_sl_3b"       ) ;
       mu_susy_sb_sl_3b         = ws->function("mu_susy_sb_sl_3b"        ) ;
       mu_susy_sig_ldp_3b       = ws->function("mu_susy_sig_ldp_3b"      ) ;
       mu_susy_sb_ldp_3b        = ws->function("mu_susy_sb_ldp_3b"       ) ;
       sf_mc                    = ws->function("sf_mc"                   ) ;
       sf_qcd_sig_1b            = ws->function("sf_qcd_sig_1b"           ) ;
       sf_qcd_sig_2b            = ws->function("sf_qcd_sig_2b"           ) ;
       sf_qcd_sig_3b            = ws->function("sf_qcd_sig_3b"           ) ;
       sf_qcd_sb_1b             = ws->function("sf_qcd_sb_1b"            ) ;
       sf_qcd_sb_2b             = ws->function("sf_qcd_sb_2b"            ) ;
       sf_qcd_sb_3b             = ws->function("sf_qcd_sb_3b"            ) ;
       sf_ttwj_sig_1b           = ws->function("sf_ttwj_sig_1b"          ) ;
       sf_ttwj_sig_2b           = ws->function("sf_ttwj_sig_2b"          ) ;
       sf_ttwj_sig_3b           = ws->function("sf_ttwj_sig_3b"          ) ;
       sf_ttwj_sb_2b            = ws->function("sf_ttwj_sb_2b"           ) ;
       sf_ttwj_sb_3b            = ws->function("sf_ttwj_sb_3b"           ) ;
       mu_qcd_sb_1b             = ws->function("mu_qcd_sb_1b"            ) ;


       if ( knn_sig_1b              == 0x0 ) { printf("\n\n\n *** can't find knn_sig_1b                 in workspace.  Quitting.\n\n\n") ; return ; }
       if ( knn_sig_2b              == 0x0 ) { printf("\n\n\n *** can't find knn_sig_2b                 in workspace.  Quitting.\n\n\n") ; return ; }
       if ( knn_sig_3b              == 0x0 ) { printf("\n\n\n *** can't find knn_sig_3b                 in workspace.  Quitting.\n\n\n") ; return ; }
       if ( knn_sb_1b               == 0x0 ) { printf("\n\n\n *** can't find knn_sb_1b                  in workspace.  Quitting.\n\n\n") ; return ; }
       if ( knn_sb_2b               == 0x0 ) { printf("\n\n\n *** can't find knn_sb_2b                  in workspace.  Quitting.\n\n\n") ; return ; }
       if ( knn_sb_3b               == 0x0 ) { printf("\n\n\n *** can't find knn_sb_3b                  in workspace.  Quitting.\n\n\n") ; return ; }
       if ( acc_ee_sig              == 0x0 ) { printf("\n\n\n *** can't find acc_ee_sig                 in workspace.  Quitting.\n\n\n") ; return ; }
       if ( acc_ee_sb               == 0x0 ) { printf("\n\n\n *** can't find acc_ee_sb                  in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fsig_ee                 == 0x0 ) { printf("\n\n\n *** can't find fsig_ee                    in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_ee                  == 0x0 ) { printf("\n\n\n *** can't find eff_ee                     in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_ee                   == 0x0 ) { printf("\n\n\n *** can't find sf_ee                      in workspace.  Quitting.\n\n\n") ; return ; }
       if ( acc_mm_sig              == 0x0 ) { printf("\n\n\n *** can't find acc_mm_sig                 in workspace.  Quitting.\n\n\n") ; return ; }
       if ( acc_mm_sb               == 0x0 ) { printf("\n\n\n *** can't find acc_mm_sb                  in workspace.  Quitting.\n\n\n") ; return ; }
       if ( fsig_mm                 == 0x0 ) { printf("\n\n\n *** can't find fsig_mm                    in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_mm                  == 0x0 ) { printf("\n\n\n *** can't find eff_mm                     in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_mm                   == 0x0 ) { printf("\n\n\n *** can't find sf_mm                      in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_ttwj_sig_ldp_1b      == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sig_ldp_1b         in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_ttwj_sb_ldp_1b       == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sb_ldp_1b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_znn_sig_ldp_1b       == 0x0 ) { printf("\n\n\n *** can't find mu_znn_sig_ldp_1b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_znn_sb_ldp_1b        == 0x0 ) { printf("\n\n\n *** can't find mu_znn_sb_ldp_1b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_ttwj_sig_ldp_2b      == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sig_ldp_2b         in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_ttwj_sb_ldp_2b       == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sb_ldp_2b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_znn_sig_ldp_2b       == 0x0 ) { printf("\n\n\n *** can't find mu_znn_sig_ldp_2b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_znn_sb_ldp_2b        == 0x0 ) { printf("\n\n\n *** can't find mu_znn_sb_ldp_2b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_ttwj_sig_ldp_3b      == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sig_ldp_3b         in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_ttwj_sb_ldp_3b       == 0x0 ) { printf("\n\n\n *** can't find mu_ttwj_sb_ldp_3b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_znn_sig_ldp_3b       == 0x0 ) { printf("\n\n\n *** can't find mu_znn_sig_ldp_3b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_znn_sb_ldp_3b        == 0x0 ) { printf("\n\n\n *** can't find mu_znn_sb_ldp_3b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sig_1b           == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sig_1b              in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sb_1b            == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sb_1b               in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sig_sl_1b        == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sig_sl_1b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sb_sl_1b         == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sb_sl_1b            in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sig_ldp_1b       == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sig_ldp_1b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sb_ldp_1b        == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sb_ldp_1b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sig_2b           == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sig_2b              in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sb_2b            == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sb_2b               in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sig_sl_2b        == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sig_sl_2b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sb_sl_2b         == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sb_sl_2b            in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sig_ldp_2b       == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sig_ldp_2b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sb_ldp_2b        == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sb_ldp_2b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sig_3b           == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sig_3b              in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sb_3b            == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sb_3b               in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sig_sl_3b        == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sig_sl_3b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sb_sl_3b         == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sb_sl_3b            in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sig_ldp_3b       == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sig_ldp_3b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( eff_sf_sb_ldp_3b        == 0x0 ) { printf("\n\n\n *** can't find eff_sf_sb_ldp_3b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sig_1b       == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sig_1b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sb_1b        == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sb_1b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sig_sl_1b    == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sig_sl_1b       in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sb_sl_1b     == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sb_sl_1b        in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sig_ldp_1b   == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sig_ldp_1b      in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sb_ldp_1b    == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sb_ldp_1b       in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sig_2b       == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sig_2b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sb_2b        == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sb_2b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sig_sl_2b    == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sig_sl_2b       in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sb_sl_2b     == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sb_sl_2b        in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sig_ldp_2b   == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sig_ldp_2b      in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sb_ldp_2b    == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sb_ldp_2b       in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sig_3b       == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sig_3b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sb_3b        == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sb_3b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sig_sl_3b    == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sig_sl_3b       in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sb_sl_3b     == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sb_sl_3b        in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sig_ldp_3b   == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sig_ldp_3b      in workspace.  Quitting.\n\n\n") ; return ; }
       if ( btageff_sf_sb_ldp_3b    == 0x0 ) { printf("\n\n\n *** can't find btageff_sf_sb_ldp_3b       in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sb_1b           == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sb_1b              in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sig_sl_1b       == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sig_sl_1b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sb_sl_1b        == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sb_sl_1b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sig_ldp_1b      == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sig_ldp_1b         in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sb_ldp_1b       == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sb_ldp_1b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sig_2b          == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sig_2b             in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sb_2b           == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sb_2b              in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sig_sl_2b       == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sig_sl_2b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sb_sl_2b        == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sb_sl_2b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sig_ldp_2b      == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sig_ldp_2b         in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sb_ldp_2b       == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sb_ldp_2b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sig_3b          == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sig_3b             in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sb_3b           == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sb_3b              in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sig_sl_3b       == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sig_sl_3b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sb_sl_3b        == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sb_sl_3b           in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sig_ldp_3b      == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sig_ldp_3b         in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_susy_sb_ldp_3b       == 0x0 ) { printf("\n\n\n *** can't find mu_susy_sb_ldp_3b          in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_mc                   == 0x0 ) { printf("\n\n\n *** can't find sf_mc                      in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_qcd_sig_1b           == 0x0 ) { printf("\n\n\n *** can't find sf_qcd_sig_1b              in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_qcd_sig_2b           == 0x0 ) { printf("\n\n\n *** can't find sf_qcd_sig_2b              in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_qcd_sig_3b           == 0x0 ) { printf("\n\n\n *** can't find sf_qcd_sig_3b              in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_qcd_sb_1b            == 0x0 ) { printf("\n\n\n *** can't find sf_qcd_sb_1b               in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_qcd_sb_2b            == 0x0 ) { printf("\n\n\n *** can't find sf_qcd_sb_2b               in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_qcd_sb_3b            == 0x0 ) { printf("\n\n\n *** can't find sf_qcd_sb_3b               in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_ttwj_sig_1b          == 0x0 ) { printf("\n\n\n *** can't find sf_ttwj_sig_1b             in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_ttwj_sig_2b          == 0x0 ) { printf("\n\n\n *** can't find sf_ttwj_sig_2b             in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_ttwj_sig_3b          == 0x0 ) { printf("\n\n\n *** can't find sf_ttwj_sig_3b             in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_ttwj_sb_2b           == 0x0 ) { printf("\n\n\n *** can't find sf_ttwj_sb_2b              in workspace.  Quitting.\n\n\n") ; return ; }
       if ( sf_ttwj_sb_3b           == 0x0 ) { printf("\n\n\n *** can't find sf_ttwj_sb_3b              in workspace.  Quitting.\n\n\n") ; return ; }
       if ( mu_qcd_sb_1b            == 0x0 ) { printf("\n\n\n *** can't find mu_qcd_sb_1b               in workspace.  Quitting.\n\n\n") ; return ; }







       //--- Get pointers to the observables.

       const RooArgSet* dsras = rds->get() ;
       TIterator* obsIter = dsras->createIterator() ;
       while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {

          if ( strcmp( obs->GetName(), "Nsig_1b"     ) == 0 ) { rrv_Nsig_1b      = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_1b"      ) == 0 ) { rrv_Nsb_1b       = obs ; }
          if ( strcmp( obs->GetName(), "Nsig_sl_1b"  ) == 0 ) { rrv_Nsig_sl_1b   = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_sl_1b"   ) == 0 ) { rrv_Nsb_sl_1b    = obs ; }
          if ( strcmp( obs->GetName(), "Nsig_ldp_1b" ) == 0 ) { rrv_Nsig_ldp_1b  = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_ldp_1b"  ) == 0 ) { rrv_Nsb_ldp_1b   = obs ; }

          if ( strcmp( obs->GetName(), "Nsig_2b"     ) == 0 ) { rrv_Nsig_2b      = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_2b"      ) == 0 ) { rrv_Nsb_2b       = obs ; }
          if ( strcmp( obs->GetName(), "Nsig_sl_2b"  ) == 0 ) { rrv_Nsig_sl_2b   = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_sl_2b"   ) == 0 ) { rrv_Nsb_sl_2b    = obs ; }
          if ( strcmp( obs->GetName(), "Nsig_ldp_2b" ) == 0 ) { rrv_Nsig_ldp_2b  = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_ldp_2b"  ) == 0 ) { rrv_Nsb_ldp_2b   = obs ; }

          if ( strcmp( obs->GetName(), "Nsig_3b"     ) == 0 ) { rrv_Nsig_3b      = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_3b"      ) == 0 ) { rrv_Nsb_3b       = obs ; }
          if ( strcmp( obs->GetName(), "Nsig_sl_3b"  ) == 0 ) { rrv_Nsig_sl_3b   = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_sl_3b"   ) == 0 ) { rrv_Nsb_sl_3b    = obs ; }
          if ( strcmp( obs->GetName(), "Nsig_ldp_3b" ) == 0 ) { rrv_Nsig_ldp_3b  = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_ldp_3b"  ) == 0 ) { rrv_Nsb_ldp_3b   = obs ; }

          if ( strcmp( obs->GetName(), "Nsig_ee"  ) == 0 ) { rrv_Nsig_ee   = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_ee"   ) == 0 ) { rrv_Nsb_ee    = obs ; }
          if ( strcmp( obs->GetName(), "Nsig_mm"  ) == 0 ) { rrv_Nsig_mm   = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_mm"   ) == 0 ) { rrv_Nsb_mm    = obs ; }

       }

       if ( rrv_Nsig_1b        == 0x0 ) { printf("\n\n\n *** can't find Nsig_1b       in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_1b         == 0x0 ) { printf("\n\n\n *** can't find Nsb_1b        in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsig_sl_1b     == 0x0 ) { printf("\n\n\n *** can't find Nsig_sl_1b    in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_sl_1b      == 0x0 ) { printf("\n\n\n *** can't find Nsb_sl_1b     in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsig_ldp_1b    == 0x0 ) { printf("\n\n\n *** can't find Nsig_ldp_1b   in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_ldp_1b     == 0x0 ) { printf("\n\n\n *** can't find Nsb_ldp_1b    in dataset.  Quitting.\n\n\n") ; return ; }

       if ( rrv_Nsig_2b        == 0x0 ) { printf("\n\n\n *** can't find Nsig_2b       in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_2b         == 0x0 ) { printf("\n\n\n *** can't find Nsb_2b        in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsig_sl_2b     == 0x0 ) { printf("\n\n\n *** can't find Nsig_sl_2b    in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_sl_2b      == 0x0 ) { printf("\n\n\n *** can't find Nsb_sl_2b     in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsig_ldp_2b    == 0x0 ) { printf("\n\n\n *** can't find Nsig_ldp_2b   in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_ldp_2b     == 0x0 ) { printf("\n\n\n *** can't find Nsb_ldp_2b    in dataset.  Quitting.\n\n\n") ; return ; }

       if ( rrv_Nsig_3b        == 0x0 ) { printf("\n\n\n *** can't find Nsig_3b       in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_3b         == 0x0 ) { printf("\n\n\n *** can't find Nsb_3b        in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsig_sl_3b     == 0x0 ) { printf("\n\n\n *** can't find Nsig_sl_3b    in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_sl_3b      == 0x0 ) { printf("\n\n\n *** can't find Nsb_sl_3b     in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsig_ldp_3b    == 0x0 ) { printf("\n\n\n *** can't find Nsig_ldp_3b   in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_ldp_3b     == 0x0 ) { printf("\n\n\n *** can't find Nsb_ldp_3b    in dataset.  Quitting.\n\n\n") ; return ; }

       if ( rrv_Nsig_ee    == 0x0 ) { printf("\n\n\n *** can't find Nsig_ee    in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_ee     == 0x0 ) { printf("\n\n\n *** can't find Nsb_ee     in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsig_mm    == 0x0 ) { printf("\n\n\n *** can't find Nsig_mm    in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_mm     == 0x0 ) { printf("\n\n\n *** can't find Nsb_mm     in dataset.  Quitting.\n\n\n") ; return ; }








       printf("\n\n\n === Model values for observables\n\n") ;

       printObservables() ;



      //--- save the actual values of the observables.

       saveObservables() ;











       //--- evaluate the test stat on the data: fit with susy floating.

       rrv_mu_susy_sig_1b->setVal( poiVal ) ;
       rrv_mu_susy_sig_1b->setConstant( kTRUE ) ;

       printf("\n\n\n ====== Fitting the data with susy fixed.\n\n") ;

       RooFitResult* dataFitResultSusyFixed = likelihood->fitTo(*rds, Save(true),Hesse(false),Minos(false),Strategy(1));
       int dataSusyFixedFitCovQual = dataFitResultSusyFixed->covQual() ;
       if ( dataSusyFixedFitCovQual < 2 ) { printf("\n\n\n *** Failed fit!  Cov qual %d.  Quitting.\n\n", dataSusyFixedFitCovQual ) ; return ; }
       double dataFitSusyFixedNll = dataFitResultSusyFixed->minNll() ;


       rrv_mu_susy_sig_1b->setVal( 0.0 ) ;
       rrv_mu_susy_sig_1b->setConstant( kFALSE ) ;

       printf("\n\n\n ====== Fitting the data with susy floating.\n\n") ;

       RooFitResult* dataFitResultSusyFloat = likelihood->fitTo(*rds, Save(true),Hesse(false),Minos(false),Strategy(1));
       int dataSusyFloatFitCovQual = dataFitResultSusyFloat->covQual() ;
       if ( dataSusyFloatFitCovQual < 2 ) { printf("\n\n\n *** Failed fit!  Cov qual %d.  Quitting.\n\n", dataSusyFloatFitCovQual ) ; return ; }
       double dataFitSusyFloatNll = dataFitResultSusyFloat->minNll() ;
       double data_maxL_mu_susy_sig_1b = rrv_mu_susy_sig_1b->getVal() ;

       double dataTestStat = 0. ;
       if ( data_maxL_mu_susy_sig_1b < poiVal ) {
          dataTestStat = 2.*( dataFitSusyFixedNll - dataFitSusyFloatNll) ;
       } else {
          printf("\n\n Setting data value of test stat to zero because mu_susy_sig_1b floated above hypo: %5.1f > %5.1f\n\n",
              data_maxL_mu_susy_sig_1b, poiVal ) ;
       }

       printf("\n\n\n Data value of test stat for mu_susy_sig_1b = %5.1f: %8.2f\n", poiVal, dataTestStat ) ;












       printf("\n\n\n === Nuisance parameters\n\n") ;

       {
          int npi(0) ;
          TIterator* npIter = nuisanceParameters->createIterator() ;
          while ( RooRealVar* np_rrv = (RooRealVar*) npIter->Next() ) {

             np_initial_val[npi] = np_rrv->getVal() ; //--- I am assuming that the order of the NPs in the iterator does not change.

             TString npname( np_rrv->GetName() ) ;
             npname.ReplaceAll("_prim","") ;
             RooAbsReal* np_rfv = ws->function( npname ) ;

             TString pdfname( np_rrv->GetName() ) ;
             pdfname.ReplaceAll("_prim","") ;
             pdfname.Prepend("pdf_") ;
             RooAbsPdf* np_pdf = ws->pdf( pdfname ) ;
             if ( np_pdf == 0x0 ) { printf("\n\n *** Can't find nuisance parameter pdf with name %s.\n\n", pdfname.Data() ) ; }

             if ( np_rfv != 0x0 ) {
                printf(" %20s : %8.2f , %20s, %8.2f\n", np_rrv->GetName(), np_rrv->getVal(), np_rfv->GetName(), np_rfv->getVal() ) ;
             } else {
                printf(" %20s : %8.2f\n", np_rrv->GetName(), np_rrv->getVal() ) ;
             }

             npi++ ;
          } // np_rrv iterator.

          np_count = npi ;

       }







       tt_dataTestStat = dataTestStat ;
       tt_hypo_mu_susy_sig_1b = poiVal ;














       printf("\n\n\n === Doing the toys\n\n") ;

       int nToyOK(0) ;
       int nToyWorseThanData(0) ;

       for ( int ti=0; ti<nToys; ti++ ) {

          printf("\n\n\n ======= Toy %4d\n\n\n", ti ) ;





          //--- 1) pick values for the nuisance parameters from the PDFs and fix them.

          {
             TIterator* npIter = nuisanceParameters->createIterator() ;
             while ( RooRealVar* np_rrv = (RooRealVar*) npIter->Next() ) {

                TString pdfname( np_rrv->GetName() ) ;
                pdfname.ReplaceAll("_prim","") ;
                pdfname.Prepend("pdf_") ;
                RooAbsPdf* np_pdf = ws->pdf( pdfname ) ;
                if ( np_pdf == 0x0 ) { printf("\n\n *** Can't find nuisance parameter pdf with name %s.\n\n", pdfname.Data() ) ; return ; }

                RooDataSet* nprds = np_pdf->generate( RooArgSet(*np_rrv) ,1) ;
                const RooArgSet* npdsras = nprds->get() ;
                TIterator* valIter = npdsras->createIterator() ;
                RooRealVar* val = (RooRealVar*) valIter->Next() ;

                //--- reset the value of the nuisance parameter and fix it for the toy model definition fit.
                np_rrv->setVal( val->getVal() ) ;
                np_rrv->setConstant( kTRUE ) ;


                TString npname( np_rrv->GetName() ) ;
                npname.ReplaceAll("_prim","") ;
                RooAbsReal* np_rfv = ws->function( npname ) ;

                if ( verbLevel > 0 ) {
                   if ( np_rfv != 0x0 ) {
                      printf(" %20s : %8.2f , %15s, %8.3f\n", val->GetName(), val->getVal(), np_rfv->GetName(), np_rfv->getVal() ) ;
               //  } else if ( strstr( npname.Data(), "eff_sf" ) != 0 ) {
               //     np_rfv = ws->function( "eff_sf_sig" ) ;
               //     RooAbsReal* np_rfv2 = ws->function( "eff_sf_sb" ) ;
               //     printf(" %20s : %8.2f , %15s, %8.3f , %15s, %8.3f\n", val->GetName(), val->getVal(), np_rfv->GetName(), np_rfv->getVal(), np_rfv2->GetName(), np_rfv2->getVal() ) ;
                   } else if ( strstr( npname.Data(), "sf_ll" ) != 0 ) {
                      np_rfv = ws->function( "sf_ee" ) ;
                      RooAbsReal* np_rfv2 = ws->function( "sf_mm" ) ;
                      printf(" %20s : %8.2f , %15s, %8.3f , %15s, %8.3f\n", val->GetName(), val->getVal(), np_rfv->GetName(), np_rfv->getVal(), np_rfv2->GetName(), np_rfv2->getVal() ) ;
                   } else {
                      printf(" %20s : %8.2f\n", val->GetName(), val->getVal() ) ;
                   }
                   cout << flush ;
                }

                delete nprds ;

             } // np_rrv iterator

             delete npIter ; // do I need to do this?  Is it safe?
          }






          //--- 2) Fit the dataset with these values for the nuisance parameters.

          if ( isBgonlyStudy ) {
            //-- fit with susy yield fixed to zero.
             rrv_mu_susy_sig_1b -> setVal( 0. ) ;
             if ( verbLevel > 0 ) { printf("\n Setting mu_susy_sig_1b to zero.\n\n") ; }
          } else {
            //-- fit with susy yield fixed to predicted value.
             rrv_mu_susy_sig_1b -> setVal( poiVal ) ;
             if ( verbLevel > 0 ) { printf("\n Setting mu_susy_sig_1b to %8.1f.\n\n", poiVal) ; }
          }
          rrv_mu_susy_sig_1b->setConstant( kTRUE ) ;

          if ( verbLevel > 0 ) {
             printf("\n\n") ;
             printf("  Fitting with these values for the observables to define the model for toy generation.\n") ;
             rds->printMultiline(cout, 1, kTRUE, "") ;
             printf("\n\n") ;
             printf("Before fit, before initialization: Current values for model predictions.\n\n") ;
             printObservables() ;
          }

          initializeFitpars() ;

          if ( verbLevel > 0 ) {
             printf("Before fit, after initialization: Current values for model predictions.\n\n") ;
             printObservables() ;
          }

          RooFitResult* toyModelDefinitionFitResult(0x0) ;
          if ( verbLevel < 2 ) {
             toyModelDefinitionFitResult = likelihood->fitTo(*rds, Save(true), PrintLevel(-1),Hesse(false),Minos(false),Strategy(1));
          } else {
             toyModelDefinitionFitResult = likelihood->fitTo(*rds, Save(true),Hesse(false),Minos(false),Strategy(1));
          }
          if ( verbLevel > 0 ) {
             printf("\n\n") ;
             printf("After fit: Current values for model predictions.\n\n") ;
             printObservables() ;
          }

          int toyModelDefFitCovQual = toyModelDefinitionFitResult->covQual() ;
          if ( verbLevel > 0 ) { printf("\n fit covariance matrix quality: %d\n\n", toyModelDefFitCovQual ) ; }
          if ( toyModelDefFitCovQual < 2 ) {
             printf("\n\n\n *** Bad toy model definition fit.  Cov qual %d.  Aborting this toy.\n\n\n", toyModelDefFitCovQual ) ;
             continue ;
          }

          delete toyModelDefinitionFitResult ;

          if ( verbLevel > 0 ) {
             printf("\n\n\n === Model values for observables.  These will be used to generate the toy dataset.\n\n") ;
             printObservables() ;
          }









          //--- 3) Generate a new set of observables based on this model.

          generateObservables() ;

          printf("\n\n\n   Generated dataset\n") ;
          rds->Print() ;
          rds->printMultiline(cout, 1, kTRUE, "") ;

          //--- Apparently, I need to make a new RooDataSet...  Resetting the
          //    values in the old one doesn't stick.  If you do likelihood->fitTo(*rds), it
          //    uses the original values, not the reset ones, in the fit.

          RooArgSet toyFitobservedParametersList ;

          toyFitobservedParametersList.add( *rrv_Nsig_1b        ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_1b         ) ;
          toyFitobservedParametersList.add( *rrv_Nsig_sl_1b     ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_sl_1b      ) ;
          toyFitobservedParametersList.add( *rrv_Nsig_ldp_1b    ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_ldp_1b     ) ;

          toyFitobservedParametersList.add( *rrv_Nsig_2b        ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_2b         ) ;
          toyFitobservedParametersList.add( *rrv_Nsig_sl_2b     ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_sl_2b      ) ;
          toyFitobservedParametersList.add( *rrv_Nsig_ldp_2b    ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_ldp_2b     ) ;

          toyFitobservedParametersList.add( *rrv_Nsig_3b        ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_3b         ) ;
          toyFitobservedParametersList.add( *rrv_Nsig_sl_3b     ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_sl_3b      ) ;
          toyFitobservedParametersList.add( *rrv_Nsig_ldp_3b    ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_ldp_3b     ) ;

          toyFitobservedParametersList.add( *rrv_Nsig_ee     ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_ee      ) ;
          toyFitobservedParametersList.add( *rrv_Nsig_mm     ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_mm      ) ;


          RooDataSet* toyFitdsObserved = new RooDataSet("toyfit_ra2b_observed_rds", "RA2b toy observed data values",
                                         toyFitobservedParametersList ) ;
          toyFitdsObserved->add( toyFitobservedParametersList ) ;





          //--- 4) Reset and free the nuisance parameters.

          {
             if ( verbLevel > 0 ) { printf("\n\n") ; }
             int npi(0) ;
             TIterator* npIter = nuisanceParameters->createIterator() ;
             while ( RooRealVar* np_rrv = (RooRealVar*) npIter->Next() ) {
                np_rrv -> setVal( np_initial_val[npi] ) ; // assuming that the order in the iterator does not change.
                np_rrv -> setConstant( kFALSE ) ;
                npi++ ;
                if ( verbLevel > 0 ) { printf("    reset %20s to %8.2f and freed it.\n", np_rrv->GetName() , np_rrv->getVal() ) ; }
             } // np_rrv iterator.
             if ( verbLevel > 0 ) { printf("\n\n") ; }
             delete npIter ; // do I need to do this?  Is it safe?
          }





          //--- 5a) Evaluate the test statistic: Fit with susy yield floating to get the absolute maximum log likelihood.

          if ( verbLevel > 0 ) { printf("\n\n  Evaluating the test statistic for this toy.  Fitting with susy floating.\n\n") ; }

          rrv_mu_susy_sig_1b->setVal( 0.0 ) ;
          rrv_mu_susy_sig_1b->setConstant( kFALSE ) ;

          if ( verbLevel > 0 ) {
             printf("\n toy dataset\n\n") ;
             toyFitdsObserved->printMultiline(cout, 1, kTRUE, "") ;
          }

     /////---- nfg.  Need to create a new dataset  ----------
     /////RooFitResult* maxLikelihoodFitResult = likelihood->fitTo(*rds, Save(true), PrintLevel(-1));
     /////RooFitResult* maxLikelihoodFitResult = likelihood->fitTo(*rds, Save(true));
     /////--------------

          if ( verbLevel > 0 ) {
             printf("\n\n") ;
             printf("Before fit, before initialization: Current values for model predictions.\n\n") ;
             printObservables() ;
          }

          initializeFitpars() ;

          if ( verbLevel > 0 ) {
             printf("Before fit, after initialization: Current values for model predictions.\n\n") ;
             printObservables() ;
          }

          RooFitResult* maxLikelihoodFitResult(0x0) ;
          if ( verbLevel < 2 ) {
             maxLikelihoodFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true), PrintLevel(-1),Hesse(false),Minos(false),Strategy(1));
          } else {
             maxLikelihoodFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true),Hesse(false),Minos(false),Strategy(1));
          }

          if ( verbLevel > 0 ) { printObservables() ; }

          int mlFitCovQual = maxLikelihoodFitResult->covQual() ;
          if ( verbLevel > 0 ) { printf("\n fit covariance matrix quality: %d , -log likelihood %f\n\n", mlFitCovQual, maxLikelihoodFitResult->minNll() ) ; }
          if ( mlFitCovQual < 2 ) {
             printf("\n\n\n *** Bad maximum likelihood fit (susy floating).  Cov qual %d.  Aborting this toy.\n\n\n", mlFitCovQual ) ;
             continue ;
          }
          double maxL_susyFloat = maxLikelihoodFitResult->minNll() ;
          double maxL_mu_susy_sig_1b = rrv_mu_susy_sig_1b->getVal() ;

          delete maxLikelihoodFitResult ;






          //--- 5b) Evaluate the test statistic: Fit with susy yield fixed to hypothesis value.
          //        This is only necessary if the maximum likelihood fit value of the susy yield
          //        is less than the hypothesis value (to get a one-sided limit).


          double testStat(0.0) ;
          double maxL_susyFixed(0.0) ;

          if ( maxL_mu_susy_sig_1b < poiVal ) {

             if ( verbLevel > 0 ) { printf("\n\n  Evaluating the test statistic for this toy.  Fitting with susy fixed to %8.2f.\n\n", poiVal ) ; }

             rrv_mu_susy_sig_1b->setVal( poiVal ) ;
             rrv_mu_susy_sig_1b->setConstant( kTRUE ) ;

             if ( verbLevel > 0 ) {
                printf("\n toy dataset\n\n") ;
                rds->printMultiline(cout, 1, kTRUE, "") ;
             }

         ////--------- nfg.  need to make a new dataset  ---------------
         ////RooFitResult* susyFixedFitResult = likelihood->fitTo(*rds, Save(true), PrintLevel(-1));
         ////RooFitResult* susyFixedFitResult = likelihood->fitTo(*rds, Save(true));
         ////-----------------------------

             if ( verbLevel > 0 ) {
                printf("\n\n") ;
                printf("Before fit, before initialization: Current values for model predictions.\n\n") ;
                printObservables() ;
             }

             initializeFitpars() ;

             if ( verbLevel > 0 ) {
                printf("Before fit, after initialization: Current values for model predictions.\n\n") ;
                printObservables() ;
             }

             RooFitResult* susyFixedFitResult(0x0) ;
             if ( verbLevel < 2 ) {
                susyFixedFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true), PrintLevel(-1),Hesse(false),Minos(false),Strategy(1));
             } else {
                susyFixedFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true),Hesse(false),Minos(false),Strategy(1));
             }

             if ( verbLevel > 0 ) { printObservables() ; }

             int susyFixedFitCovQual = susyFixedFitResult->covQual() ;
             if ( verbLevel > 0 ) { printf("\n fit covariance matrix quality: %d , -log likelihood %f\n\n", susyFixedFitCovQual, susyFixedFitResult->minNll()  ) ; }
             if ( susyFixedFitCovQual < 2 ) {
                printf("\n\n\n *** Bad maximum likelihood fit (susy fixed).  Cov qual %d.  Aborting this toy.\n\n\n", susyFixedFitCovQual ) ;
                continue ;
             }
             maxL_susyFixed = susyFixedFitResult->minNll() ;
             testStat = 2. * (maxL_susyFixed - maxL_susyFloat) ;


             delete susyFixedFitResult ;


          } else {

             if ( verbLevel > 0 ) { printf("\n\n  Floating value of susy yield greater than hypo value (%8.2f > %8.2f).  Setting test stat to zero.\n\n", maxL_mu_susy_sig_1b, poiVal ) ; }

             testStat = 0.0 ;

          }

          printf("   --- test stat for toy %4d : %8.2f\n", ti, testStat ) ;





          nToyOK++ ;

          if ( testStat >= dataTestStat ) { nToyWorseThanData++ ; }


          if ( makeTtree ) {

             tt_testStat = testStat ;

             tt_gen_Nsig_1b     = rrv_Nsig_1b->getVal() ;
             tt_gen_Nsb_1b      = rrv_Nsb_1b->getVal() ;
             tt_gen_Nsig_sl_1b  = rrv_Nsig_sl_1b->getVal() ;
             tt_gen_Nsb_sl_1b   = rrv_Nsb_sl_1b->getVal() ;
             tt_gen_Nsig_ldp_1b = rrv_Nsig_ldp_1b->getVal() ;
             tt_gen_Nsb_ldp_1b  = rrv_Nsb_ldp_1b->getVal() ;

             tt_gen_Nsig_2b     = rrv_Nsig_2b->getVal() ;
             tt_gen_Nsb_2b      = rrv_Nsb_2b->getVal() ;
             tt_gen_Nsig_sl_2b  = rrv_Nsig_sl_2b->getVal() ;
             tt_gen_Nsb_sl_2b   = rrv_Nsb_sl_2b->getVal() ;
             tt_gen_Nsig_ldp_2b = rrv_Nsig_ldp_2b->getVal() ;
             tt_gen_Nsb_ldp_2b  = rrv_Nsb_ldp_2b->getVal() ;

             tt_gen_Nsig_3b     = rrv_Nsig_3b->getVal() ;
             tt_gen_Nsb_3b      = rrv_Nsb_3b->getVal() ;
             tt_gen_Nsig_sl_3b  = rrv_Nsig_sl_3b->getVal() ;
             tt_gen_Nsb_sl_3b   = rrv_Nsb_sl_3b->getVal() ;
             tt_gen_Nsig_ldp_3b = rrv_Nsig_ldp_3b->getVal() ;
             tt_gen_Nsb_ldp_3b  = rrv_Nsb_ldp_3b->getVal() ;

             tt_gen_Nsig_ee  = rrv_Nsig_ee->getVal() ;
             tt_gen_Nsb_ee   = rrv_Nsb_ee->getVal() ;
             tt_gen_Nsig_mm  = rrv_Nsig_mm->getVal() ;
             tt_gen_Nsb_mm   = rrv_Nsb_mm->getVal() ;

             toytt->Fill() ;

          }





          //--- *) reset things for the next toy.

          resetObservables() ;

          delete toyFitdsObserved ;




       } // ti.

       wstf->Close() ;

       printf("\n\n\n") ;

       if ( nToyOK == 0 ) { printf("\n\n\n *** All toys bad !?!?!\n\n\n") ; return ; }

       double pValue = (1.0*nToyWorseThanData) / (1.0*nToyOK) ;

       if ( isBgonlyStudy ) {
          printf("\n\n\n p-value result, BG-only , poi=%3.0f : %4d / %4d = %6.3f\n\n\n\n", poiVal, nToyWorseThanData, nToyOK, pValue ) ;
       } else {
          printf("\n\n\n p-value result, S-plus-B, poi=%3.0f : %4d / %4d = %6.3f\n\n\n\n", poiVal, nToyWorseThanData, nToyOK, pValue ) ;
       }


       if ( makeTtree ) {
          printf("\n\n Writing TTree : %s : %s\n\n", ttname, tttitle ) ;
          ttfile->cd() ;
          toytt->Write() ;
       }


   } // ws_cls_hybrid1


  //==============================================================================================


   void printObservables() {

       printf("\n\n") ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_1b     ->GetName(), rfv_n_sig_1b     ->getVal(),  rrv_Nsig_1b     ->GetName(), rrv_Nsig_1b     ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_1b      ->GetName(), rfv_n_sb_1b      ->getVal(),  rrv_Nsb_1b      ->GetName(), rrv_Nsb_1b      ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_sl_1b  ->GetName(), rfv_n_sig_sl_1b  ->getVal(),  rrv_Nsig_sl_1b  ->GetName(), rrv_Nsig_sl_1b  ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_sl_1b   ->GetName(), rfv_n_sb_sl_1b   ->getVal(),  rrv_Nsb_sl_1b   ->GetName(), rrv_Nsb_sl_1b   ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_ldp_1b ->GetName(), rfv_n_sig_ldp_1b ->getVal(),  rrv_Nsig_ldp_1b ->GetName(), rrv_Nsig_ldp_1b ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_ldp_1b  ->GetName(), rfv_n_sb_ldp_1b  ->getVal(),  rrv_Nsb_ldp_1b  ->GetName(), rrv_Nsb_ldp_1b  ->getVal() ) ;

       printf("\n") ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_2b     ->GetName(), rfv_n_sig_2b     ->getVal(),  rrv_Nsig_2b     ->GetName(), rrv_Nsig_2b     ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_2b      ->GetName(), rfv_n_sb_2b      ->getVal(),  rrv_Nsb_2b      ->GetName(), rrv_Nsb_2b      ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_sl_2b  ->GetName(), rfv_n_sig_sl_2b  ->getVal(),  rrv_Nsig_sl_2b  ->GetName(), rrv_Nsig_sl_2b  ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_sl_2b   ->GetName(), rfv_n_sb_sl_2b   ->getVal(),  rrv_Nsb_sl_2b   ->GetName(), rrv_Nsb_sl_2b   ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_ldp_2b ->GetName(), rfv_n_sig_ldp_2b ->getVal(),  rrv_Nsig_ldp_2b ->GetName(), rrv_Nsig_ldp_2b ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_ldp_2b  ->GetName(), rfv_n_sb_ldp_2b  ->getVal(),  rrv_Nsb_ldp_2b  ->GetName(), rrv_Nsb_ldp_2b  ->getVal() ) ;

       printf("\n") ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_3b     ->GetName(), rfv_n_sig_3b     ->getVal(),  rrv_Nsig_3b     ->GetName(), rrv_Nsig_3b     ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_3b      ->GetName(), rfv_n_sb_3b      ->getVal(),  rrv_Nsb_3b      ->GetName(), rrv_Nsb_3b      ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_sl_3b  ->GetName(), rfv_n_sig_sl_3b  ->getVal(),  rrv_Nsig_sl_3b  ->GetName(), rrv_Nsig_sl_3b  ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_sl_3b   ->GetName(), rfv_n_sb_sl_3b   ->getVal(),  rrv_Nsb_sl_3b   ->GetName(), rrv_Nsb_sl_3b   ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_ldp_3b ->GetName(), rfv_n_sig_ldp_3b ->getVal(),  rrv_Nsig_ldp_3b ->GetName(), rrv_Nsig_ldp_3b ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_ldp_3b  ->GetName(), rfv_n_sb_ldp_3b  ->getVal(),  rrv_Nsb_ldp_3b  ->GetName(), rrv_Nsb_ldp_3b  ->getVal() ) ;

       printf("\n") ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_ee  ->GetName(), rfv_n_sig_ee  ->getVal(),  rrv_Nsig_ee  ->GetName(), rrv_Nsig_ee  ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_ee   ->GetName(), rfv_n_sb_ee   ->getVal(),  rrv_Nsb_ee   ->GetName(), rrv_Nsb_ee   ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_mm  ->GetName(), rfv_n_sig_mm  ->getVal(),  rrv_Nsig_mm  ->GetName(), rrv_Nsig_mm  ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_mm   ->GetName(), rfv_n_sb_mm   ->getVal(),  rrv_Nsb_mm   ->GetName(), rrv_Nsb_mm   ->getVal() ) ;

       printf("\n\n") ;

   }


  //==============================================================================================



   void saveObservables() {

       actual_Nsig_1b     =  rrv_Nsig_1b     ->getVal()  ;
       actual_Nsb_1b      =  rrv_Nsb_1b      ->getVal()  ;
       actual_Nsig_sl_1b  =  rrv_Nsig_sl_1b  ->getVal()  ;
       actual_Nsb_sl_1b   =  rrv_Nsb_sl_1b   ->getVal()  ;
       actual_Nsig_ldp_1b =  rrv_Nsig_ldp_1b ->getVal()  ;
       actual_Nsb_ldp_1b  =  rrv_Nsb_ldp_1b  ->getVal()  ;

       actual_Nsig_2b     =  rrv_Nsig_2b     ->getVal()  ;
       actual_Nsb_2b      =  rrv_Nsb_2b      ->getVal()  ;
       actual_Nsig_sl_2b  =  rrv_Nsig_sl_2b  ->getVal()  ;
       actual_Nsb_sl_2b   =  rrv_Nsb_sl_2b   ->getVal()  ;
       actual_Nsig_ldp_2b =  rrv_Nsig_ldp_2b ->getVal()  ;
       actual_Nsb_ldp_2b  =  rrv_Nsb_ldp_2b  ->getVal()  ;

       actual_Nsig_3b     =  rrv_Nsig_3b     ->getVal()  ;
       actual_Nsb_3b      =  rrv_Nsb_3b      ->getVal()  ;
       actual_Nsig_sl_3b  =  rrv_Nsig_sl_3b  ->getVal()  ;
       actual_Nsb_sl_3b   =  rrv_Nsb_sl_3b   ->getVal()  ;
       actual_Nsig_ldp_3b =  rrv_Nsig_ldp_3b ->getVal()  ;
       actual_Nsb_ldp_3b  =  rrv_Nsb_ldp_3b  ->getVal()  ;

       actual_Nsig_ee  =  rrv_Nsig_ee  ->getVal()  ;
       actual_Nsb_ee   =  rrv_Nsb_ee   ->getVal()  ;
       actual_Nsig_mm  =  rrv_Nsig_mm  ->getVal()  ;
       actual_Nsb_mm   =  rrv_Nsb_mm   ->getVal()  ;

   }


  //==============================================================================================



   void resetObservables() {

       rrv_Nsig_1b     ->setVal(actual_Nsig_1b    )  ;
       rrv_Nsb_1b      ->setVal(actual_Nsb_1b     )  ;
       rrv_Nsig_sl_1b  ->setVal(actual_Nsig_sl_1b )  ;
       rrv_Nsb_sl_1b   ->setVal(actual_Nsb_sl_1b  )  ;
       rrv_Nsig_ldp_1b ->setVal(actual_Nsig_ldp_1b)  ;
       rrv_Nsb_ldp_1b  ->setVal(actual_Nsb_ldp_1b )  ;

       rrv_Nsig_2b     ->setVal(actual_Nsig_2b    )  ;
       rrv_Nsb_2b      ->setVal(actual_Nsb_2b     )  ;
       rrv_Nsig_sl_2b  ->setVal(actual_Nsig_sl_2b )  ;
       rrv_Nsb_sl_2b   ->setVal(actual_Nsb_sl_2b  )  ;
       rrv_Nsig_ldp_2b ->setVal(actual_Nsig_ldp_2b)  ;
       rrv_Nsb_ldp_2b  ->setVal(actual_Nsb_ldp_2b )  ;

       rrv_Nsig_3b     ->setVal(actual_Nsig_3b    )  ;
       rrv_Nsb_3b      ->setVal(actual_Nsb_3b     )  ;
       rrv_Nsig_sl_3b  ->setVal(actual_Nsig_sl_3b )  ;
       rrv_Nsb_sl_3b   ->setVal(actual_Nsb_sl_3b  )  ;
       rrv_Nsig_ldp_3b ->setVal(actual_Nsig_ldp_3b)  ;
       rrv_Nsb_ldp_3b  ->setVal(actual_Nsb_ldp_3b )  ;

       rrv_Nsig_ee  ->setVal(actual_Nsig_ee )  ;
       rrv_Nsb_ee   ->setVal(actual_Nsb_ee  )  ;
       rrv_Nsig_mm  ->setVal(actual_Nsig_mm )  ;
       rrv_Nsb_mm   ->setVal(actual_Nsb_mm  )  ;

   }


  //==============================================================================================


   void generateObservables() {

       rrv_Nsig_1b     ->setVal( random_ng->Poisson( rfv_n_sig_1b    ->getVal() ) )  ;
       rrv_Nsb_1b      ->setVal( random_ng->Poisson( rfv_n_sb_1b     ->getVal() ) )  ;
       rrv_Nsig_sl_1b  ->setVal( random_ng->Poisson( rfv_n_sig_sl_1b ->getVal() ) )  ;
       rrv_Nsb_sl_1b   ->setVal( random_ng->Poisson( rfv_n_sb_sl_1b  ->getVal() ) )  ;
       rrv_Nsig_ldp_1b ->setVal( random_ng->Poisson( rfv_n_sig_ldp_1b->getVal() ) )  ;
       rrv_Nsb_ldp_1b  ->setVal( random_ng->Poisson( rfv_n_sb_ldp_1b ->getVal() ) )  ;

       rrv_Nsig_2b     ->setVal( random_ng->Poisson( rfv_n_sig_2b    ->getVal() ) )  ;
       rrv_Nsb_2b      ->setVal( random_ng->Poisson( rfv_n_sb_2b     ->getVal() ) )  ;
       rrv_Nsig_sl_2b  ->setVal( random_ng->Poisson( rfv_n_sig_sl_2b ->getVal() ) )  ;
       rrv_Nsb_sl_2b   ->setVal( random_ng->Poisson( rfv_n_sb_sl_2b  ->getVal() ) )  ;
       rrv_Nsig_ldp_2b ->setVal( random_ng->Poisson( rfv_n_sig_ldp_2b->getVal() ) )  ;
       rrv_Nsb_ldp_2b  ->setVal( random_ng->Poisson( rfv_n_sb_ldp_2b ->getVal() ) )  ;

       rrv_Nsig_3b     ->setVal( random_ng->Poisson( rfv_n_sig_3b    ->getVal() ) )  ;
       rrv_Nsb_3b      ->setVal( random_ng->Poisson( rfv_n_sb_3b     ->getVal() ) )  ;
       rrv_Nsig_sl_3b  ->setVal( random_ng->Poisson( rfv_n_sig_sl_3b ->getVal() ) )  ;
       rrv_Nsb_sl_3b   ->setVal( random_ng->Poisson( rfv_n_sb_sl_3b  ->getVal() ) )  ;
       rrv_Nsig_ldp_3b ->setVal( random_ng->Poisson( rfv_n_sig_ldp_3b->getVal() ) )  ;
       rrv_Nsb_ldp_3b  ->setVal( random_ng->Poisson( rfv_n_sb_ldp_3b ->getVal() ) )  ;

       rrv_Nsig_ee  ->setVal( random_ng->Poisson( rfv_n_sig_ee ->getVal() ) )  ;
       rrv_Nsb_ee   ->setVal( random_ng->Poisson( rfv_n_sb_ee  ->getVal() ) )  ;
       rrv_Nsig_mm  ->setVal( random_ng->Poisson( rfv_n_sig_mm ->getVal() ) )  ;
       rrv_Nsb_mm   ->setVal( random_ng->Poisson( rfv_n_sb_mm  ->getVal() ) )  ;

   }

  //==============================================================================================

   void initializeFitpars() {

     //----

       double ee_sig_denom = sf_ee->getVal()  *  acc_ee_sig->getVal()  *  eff_ee->getVal() ;
       double mm_sig_denom = sf_mm->getVal()  *  acc_mm_sig->getVal()  *  eff_mm->getVal() ;
       double initval_mu_znn_sig_1b = rrv_Nsig_1b->getVal() ;
       if ( ee_sig_denom > 0 && mm_sig_denom > 0 ) {
          initval_mu_znn_sig_1b = 0.5 * (  ( rrv_Nsig_ee->getVal()  *  fsig_ee->getVal()  *  znnoverll_bfratio->getVal()  *  knn_sig_1b->getVal() ) / ee_sig_denom
                                        +  ( rrv_Nsig_mm->getVal()  *  fsig_mm->getVal()  *  znnoverll_bfratio->getVal()  *  knn_sig_1b->getVal() ) / mm_sig_denom ) ;
       } else {
          printf("\n\n\n *** initializeFitpars :: zero denominator!\n\n\n") ;
       }

       fv_mu_znn_sig_1b -> setVal( initval_mu_znn_sig_1b ) ;

     //----


       double ee_sb_denom = sf_ee->getVal()  *  acc_ee_sb->getVal()  *  eff_ee->getVal() ;
       double mm_sb_denom = sf_mm->getVal()  *  acc_mm_sb->getVal()  *  eff_mm->getVal() ;
       double initval_mu_znn_sb_1b = rrv_Nsb_1b->getVal() ;
       if ( ee_sb_denom > 0 && mm_sb_denom > 0 ) {
          initval_mu_znn_sb_1b = 0.5 * (  ( rrv_Nsb_ee->getVal()  *  fsig_ee->getVal()  *  znnoverll_bfratio->getVal()  *  knn_sb_1b->getVal() ) / ee_sb_denom
                                        + ( rrv_Nsb_mm->getVal()  *  fsig_mm->getVal()  *  znnoverll_bfratio->getVal()  *  knn_sb_1b->getVal() ) / mm_sb_denom ) ;
       } else {
          printf("\n\n\n *** initializeFitpars :: zero denominator!\n\n\n") ;
       }

       fv_mu_znn_sb_1b -> setVal( initval_mu_znn_sb_1b ) ;

     //----

       double initval_mu_ttwj_sig_sl_1b =  rrv_Nsig_sl_1b->getVal()  -  btageff_sf_sig_sl_1b->getVal() * eff_sf_sig_sl_1b->getVal() * mu_susy_sig_sl_1b->getVal()  ;
       if ( initval_mu_ttwj_sig_sl_1b < 0 ) { initval_mu_ttwj_sig_sl_1b = 0. ; }
       fv_mu_ttwj_sig_sl_1b -> setVal(  initval_mu_ttwj_sig_sl_1b ) ;


     //----

       double initval_mu_ttwj_sb_sl_1b =  rrv_Nsb_sl_1b->getVal()  -  btageff_sf_sb_sl_1b->getVal() * eff_sf_sb_sl_1b->getVal() * mu_susy_sb_sl_1b->getVal()  ;
       if ( initval_mu_ttwj_sb_sl_1b < 0 ) { initval_mu_ttwj_sb_sl_1b = 0. ; }
       fv_mu_ttwj_sb_sl_1b -> setVal(  initval_mu_ttwj_sb_sl_1b ) ;


     //----

       double initval_mu_ttwj_sig_sl_2b =  rrv_Nsig_sl_2b->getVal()  -  btageff_sf_sig_sl_2b->getVal() * eff_sf_sig_sl_2b->getVal() * mu_susy_sig_sl_2b->getVal()  ;
       if ( initval_mu_ttwj_sig_sl_2b < 0 ) { initval_mu_ttwj_sig_sl_2b = 0. ; }
       fv_mu_ttwj_sig_sl_2b -> setVal(  initval_mu_ttwj_sig_sl_2b ) ;


     //----

       double initval_mu_ttwj_sb_sl_2b =  rrv_Nsb_sl_2b->getVal()  -  btageff_sf_sb_sl_2b->getVal() * eff_sf_sb_sl_2b->getVal() * mu_susy_sb_sl_2b->getVal()  ;
       if ( initval_mu_ttwj_sb_sl_2b < 0 ) { initval_mu_ttwj_sb_sl_2b = 0. ; }
       fv_mu_ttwj_sb_sl_2b -> setVal(  initval_mu_ttwj_sb_sl_2b ) ;


     //----

       double initval_mu_ttwj_sig_sl_3b =  rrv_Nsig_sl_3b->getVal()  -  btageff_sf_sig_sl_3b->getVal() * eff_sf_sig_sl_3b->getVal() * mu_susy_sig_sl_3b->getVal()  ;
       if ( initval_mu_ttwj_sig_sl_3b < 0 ) { initval_mu_ttwj_sig_sl_3b = 0. ; }
       fv_mu_ttwj_sig_sl_3b -> setVal(  initval_mu_ttwj_sig_sl_3b ) ;


     //----

       double initval_mu_ttwj_sb_sl_3b =  rrv_Nsb_sl_3b->getVal()  -  btageff_sf_sb_sl_3b->getVal() * eff_sf_sb_sl_3b->getVal() * mu_susy_sb_sl_3b->getVal()  ;
       if ( initval_mu_ttwj_sb_sl_3b < 0 ) { initval_mu_ttwj_sb_sl_3b = 0. ; }
       fv_mu_ttwj_sb_sl_3b -> setVal(  initval_mu_ttwj_sb_sl_3b ) ;












     //----

       double initval_mu_qcd_sb_ldp_1b =
                                     ( rrv_Nsb_ldp_1b->getVal()  -  btageff_sf_sb_ldp_1b->getVal() * eff_sf_sb_ldp_1b->getVal() * ( sf_mc->getVal() * ( mu_ttwj_sb_ldp_1b->getVal() + mu_znn_sb_ldp_1b->getVal() )
                                                                                        + mu_susy_sb_ldp_1b->getVal() ) ) ;

       if ( initval_mu_qcd_sb_ldp_1b < 0 ) { initval_mu_qcd_sb_ldp_1b = 0. ; }
       fv_mu_qcd_sb_ldp_1b -> setVal( initval_mu_qcd_sb_ldp_1b ) ;

     //----

       double initval_mu_qcd_sig_ldp_1b =
                                     ( rrv_Nsig_ldp_1b->getVal()  -  btageff_sf_sig_ldp_1b->getVal() * eff_sf_sig_ldp_1b->getVal() * ( sf_mc->getVal() * ( mu_ttwj_sig_ldp_1b->getVal() + mu_znn_sig_ldp_1b->getVal() )
                                                                                        + mu_susy_sig_ldp_1b->getVal() ) ) ;

       if ( initval_mu_qcd_sig_ldp_1b < 0 ) { initval_mu_qcd_sig_ldp_1b = 0. ; }
       fv_mu_qcd_sig_ldp_1b -> setVal( initval_mu_qcd_sig_ldp_1b ) ;


     //----

       double initval_mu_qcd_sb_ldp_2b =
                                     ( rrv_Nsb_ldp_2b->getVal()  -  btageff_sf_sb_ldp_2b->getVal() * eff_sf_sb_ldp_2b->getVal() * ( sf_mc->getVal() * ( mu_ttwj_sb_ldp_2b->getVal() + mu_znn_sb_ldp_2b->getVal() )
                                                                                        + mu_susy_sb_ldp_2b->getVal() ) ) ;

       if ( initval_mu_qcd_sb_ldp_2b < 0 ) { initval_mu_qcd_sb_ldp_2b = 0. ; }
       fv_mu_qcd_sb_ldp_2b -> setVal( initval_mu_qcd_sb_ldp_2b ) ;

     //----

       double initval_mu_qcd_sig_ldp_2b =
                                     ( rrv_Nsig_ldp_2b->getVal()  -  btageff_sf_sig_ldp_2b->getVal() * eff_sf_sig_ldp_2b->getVal() * ( sf_mc->getVal() * ( mu_ttwj_sig_ldp_2b->getVal() + mu_znn_sig_ldp_2b->getVal() )
                                                                                        + mu_susy_sig_ldp_2b->getVal() ) ) ;

       if ( initval_mu_qcd_sig_ldp_2b < 0 ) { initval_mu_qcd_sig_ldp_2b = 0. ; }
       fv_mu_qcd_sig_ldp_2b -> setVal( initval_mu_qcd_sig_ldp_2b ) ;


     //----

       double initval_mu_qcd_sb_ldp_3b =
                                     ( rrv_Nsb_ldp_3b->getVal()  -  btageff_sf_sb_ldp_3b->getVal() * eff_sf_sb_ldp_3b->getVal() * ( sf_mc->getVal() * ( mu_ttwj_sb_ldp_3b->getVal() + mu_znn_sb_ldp_3b->getVal() )
                                                                                        + mu_susy_sb_ldp_3b->getVal() ) ) ;

       if ( initval_mu_qcd_sb_ldp_3b < 0 ) { initval_mu_qcd_sb_ldp_3b = 0. ; }
       fv_mu_qcd_sb_ldp_3b -> setVal( initval_mu_qcd_sb_ldp_3b ) ;

     //----

       double initval_mu_qcd_sig_ldp_3b =
                                     ( rrv_Nsig_ldp_3b->getVal()  -  btageff_sf_sig_ldp_3b->getVal() * eff_sf_sig_ldp_3b->getVal() * ( sf_mc->getVal() * ( mu_ttwj_sig_ldp_3b->getVal() + mu_znn_sig_ldp_3b->getVal() )
                                                                                        + mu_susy_sig_ldp_3b->getVal() ) ) ;

       if ( initval_mu_qcd_sig_ldp_3b < 0 ) { initval_mu_qcd_sig_ldp_3b = 0. ; }
       fv_mu_qcd_sig_ldp_3b -> setVal( initval_mu_qcd_sig_ldp_3b ) ;


     //----











     //----

       if ( fv_mu_ttwj_sb_1b != 0x0 ) {
          double initval_mu_ttwj_sb_1b = rrv_Nsb_1b->getVal() - mu_qcd_sb_1b->getVal() - fv_mu_znn_sb_1b->getVal() - btageff_sf_sb_1b->getVal() * eff_sf_sb_1b->getVal() * mu_susy_sb_1b->getVal() ;
          if ( initval_mu_ttwj_sb_1b < 0. ) { initval_mu_ttwj_sb_1b = 0. ; }
          fv_mu_ttwj_sb_1b -> setVal( initval_mu_ttwj_sb_1b ) ;
       }


   }

  //==============================================================================================




