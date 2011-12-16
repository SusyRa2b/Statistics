//
//   Owen Long, UCR
//   Harrison Prosper, FSU
//   Sezen Sekmen, FSU
//
//


#include "ra2bRoostatsClass8.h"

#include <iostream>
#include <string.h>


#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TH2F.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TStringLong.h"

#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooTrace.h"
#include "RooUniform.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"

//////////#include "LikelihoodIntervalPlot.cxx"

  using namespace RooFit ;
  using namespace RooStats ;


  //=====================================================================================================


   ra2bRoostatsClass8::ra2bRoostatsClass8( bool ArgUseSigTtwjVar, bool ArgUseLdpVars ) {

      gStyle->SetOptStat(0) ;

      useSigTtwjVar = ArgUseSigTtwjVar ;
      useLdpVars = ArgUseLdpVars ;


     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  /// RooMsgService::instance().addStream(DEBUG,Topic(Tracing),OutputFile("debug.log")) ;
      printf("\n\n ==== RooFit output configuration ===============\n") ;
      RooMsgService::instance().Print("v") ;
      printf("\n\n ================================================\n") ;

      varsAtFitVals = false ;
      initialized = false ;




   }





  //===================================================================

   ra2bRoostatsClass8::~ra2bRoostatsClass8() { }



  //===================================================================

    bool ra2bRoostatsClass8::initialize( const char* infile ,
                                         const char* inputScanFile, double m0, double m12, bool isT1bbbb, double t1bbbbXsec ) {


       printf( "\n\n Opening input file : %s\n\n", infile ) ;

       FILE* infp ;
       if ( (infp=fopen( infile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
          return false ;
       }

       sprintf( initializeFile, "%s", infile ) ;


       int    Nsig                  ;
       int    Nsb                   ;
       int    Nsig_sl               ;
       int    Nsb_sl                ;
       int    Nsig_ldp              ;
       int    Nsb_ldp               ;
   /// int    Nhtonlytrig_lsb_0b      ;
   /// int    Nhtonlytrig_lsb_0b_ldp  ;
       int    Nsb_ee                  ;
       int    Nsig_ee                 ;
       int    Nsb_mm                  ;
       int    Nsig_mm                 ;
       float  Nttbarmc_sig_ldp      ;
       float  Nttbarmc_sb_ldp       ;
       float  NWJmc_sig_ldp         ;
       float  NWJmc_sb_ldp          ;
       float  NZnnmc_sig_ldp        ;
       float  NZnnmc_sb_ldp         ;
       float  NEwomc_sig_ldp        ;
       float  NEwomc_sb_ldp         ;
       float  NEwomc_sig  ;
       float  NEwomc_sb   ;


       float  sf_mc            ;
       float  sf_mc_err        ;
       float  sf_qcd_sb        ;
       float  sf_qcd_sb_err    ;
       float  sf_qcd_sig       ;
       float  sf_qcd_sig_err   ;
       float  sf_ttwj_sig      ;
       float  sf_ttwj_sig_err  ;
       float  sf_ee            ;
       float  sf_ee_err        ;
       float  sf_mm            ;
       float  sf_mm_err        ;

       float  Rlsb_passfail     ;
       float  Rlsb_passfail_err ;









       //--- read in description line.
       printf("\n\n") ;
       char c(0) ;
       while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
       printf("\n\n") ;


       char label[1000] ;

      //-----    The order here matters!

       fscanf( infp, "%s %d", label, &Nsig                  ) ;   printf( "%s %d\n", label, Nsig                  ) ;
       fscanf( infp, "%s %d", label, &Nsb                   ) ;   printf( "%s %d\n", label, Nsb                   ) ;
       fscanf( infp, "%s %d", label, &Nsig_sl               ) ;   printf( "%s %d\n", label, Nsig_sl               ) ;
       fscanf( infp, "%s %d", label, &Nsb_sl                ) ;   printf( "%s %d\n", label, Nsb_sl                ) ;
       fscanf( infp, "%s %d", label, &Nsig_ldp              ) ;   printf( "%s %d\n", label, Nsig_ldp              ) ;
       fscanf( infp, "%s %d", label, &Nsb_ldp               ) ;   printf( "%s %d\n", label, Nsb_ldp               ) ;
       fscanf( infp, "%s %d", label, &Nsb_ee                    ) ;   printf( "%s %d\n", label, Nsb_ee                    ) ;
       fscanf( infp, "%s %d", label, &Nsig_ee                   ) ;   printf( "%s %d\n", label, Nsig_ee                   ) ;
       fscanf( infp, "%s %d", label, &Nsb_mm                    ) ;   printf( "%s %d\n", label, Nsb_mm                    ) ;
       fscanf( infp, "%s %d", label, &Nsig_mm                   ) ;   printf( "%s %d\n", label, Nsig_mm                   ) ;
       fscanf( infp, "%s %g", label, &Rlsb_passfail         ) ;   printf( "%s %g\n", label, Rlsb_passfail         ) ;
       fscanf( infp, "%s %g", label, &Rlsb_passfail_err     ) ;   printf( "%s %g\n", label, Rlsb_passfail_err     ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sig_ldp      ) ;   printf( "%s %g\n", label, Nttbarmc_sig_ldp      ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sb_ldp       ) ;   printf( "%s %g\n", label, Nttbarmc_sb_ldp       ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig_ldp         ) ;   printf( "%s %g\n", label, NWJmc_sig_ldp         ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb_ldp          ) ;   printf( "%s %g\n", label, NWJmc_sb_ldp          ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig_ldp        ) ;   printf( "%s %g\n", label, NZnnmc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb_ldp         ) ;   printf( "%s %g\n", label, NZnnmc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig_ldp        ) ;   printf( "%s %g\n", label, NEwomc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb_ldp         ) ;   printf( "%s %g\n", label, NEwomc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig        ) ;   printf( "%s %g\n", label, NEwomc_sig        ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb         ) ;   printf( "%s %g\n", label, NEwomc_sb         ) ;
       fscanf( infp, "%s %g", label, &DataLumi                  ) ;   printf( "%s %g\n", label, DataLumi                  ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sig_mean              ) ;   printf( "%s %g\n", label, acc_ee_sig_mean               ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sig_err               ) ;   printf( "%s %g\n", label, acc_ee_sig_err                ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sb_mean               ) ;   printf( "%s %g\n", label, acc_ee_sb_mean                ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sb_err                ) ;   printf( "%s %g\n", label, acc_ee_sb_err                 ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sig_mean              ) ;   printf( "%s %g\n", label, acc_mm_sig_mean               ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sig_err               ) ;   printf( "%s %g\n", label, acc_mm_sig_err                ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sb_mean               ) ;   printf( "%s %g\n", label, acc_mm_sb_mean                ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sb_err                ) ;   printf( "%s %g\n", label, acc_mm_sb_err                 ) ;
       fscanf( infp, "%s %g", label, &eff_ee_mean               ) ;   printf( "%s %g\n", label, eff_ee_mean               ) ;
       fscanf( infp, "%s %g", label, &eff_ee_err                ) ;   printf( "%s %g\n", label, eff_ee_err                ) ;
       fscanf( infp, "%s %g", label, &eff_mm_mean               ) ;   printf( "%s %g\n", label, eff_mm_mean               ) ;
       fscanf( infp, "%s %g", label, &eff_mm_err                ) ;   printf( "%s %g\n", label, eff_mm_err                ) ;
       fscanf( infp, "%s %g", label, &Ztoll_lumi                ) ;   printf( "%s %g\n", label, Ztoll_lumi                ) ;
       fscanf( infp, "%s %g", label, &knn_ee_sig_mean              ) ;   printf( "%s %g\n", label, knn_ee_sig_mean              ) ;
       fscanf( infp, "%s %g", label, &knn_ee_sig_err               ) ;   printf( "%s %g\n", label, knn_ee_sig_err               ) ;
       fscanf( infp, "%s %g", label, &knn_ee_sb_mean               ) ;   printf( "%s %g\n", label, knn_ee_sb_mean               ) ;
       fscanf( infp, "%s %g", label, &knn_ee_sb_err                ) ;   printf( "%s %g\n", label, knn_ee_sb_err                ) ;
       fscanf( infp, "%s %g", label, &knn_mm_sig_mean              ) ;   printf( "%s %g\n", label, knn_mm_sig_mean              ) ;
       fscanf( infp, "%s %g", label, &knn_mm_sig_err               ) ;   printf( "%s %g\n", label, knn_mm_sig_err               ) ;
       fscanf( infp, "%s %g", label, &knn_mm_sb_mean               ) ;   printf( "%s %g\n", label, knn_mm_sb_mean               ) ;
       fscanf( infp, "%s %g", label, &knn_mm_sb_err                ) ;   printf( "%s %g\n", label, knn_mm_sb_err                ) ;
       fscanf( infp, "%s %g", label, &fsig_ee_mean              ) ;   printf( "fsig_ee_mean    %s %g\n", label, fsig_ee_mean              ) ;
       fscanf( infp, "%s %g", label, &fsig_ee_err               ) ;   printf( "fsig_ee_err     %s %g\n", label, fsig_ee_err               ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_mean              ) ;   printf( "fsig_mm_mean    %s %g\n", label, fsig_mm_mean              ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_err               ) ;   printf( "fsig_mm_err     %s %g\n", label, fsig_mm_err               ) ;
       fscanf( infp, "%s %g", label, &sf_mc                     ) ;   printf( "sf_mc           %s %g\n", label, sf_mc                     ) ;
       fscanf( infp, "%s %g", label, &sf_mc_err                 ) ;   printf( "sf_mc_err       %s %g\n", label, sf_mc_err                 ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb                 ) ;   printf( "sf_qcd_sb       %s %g\n", label, sf_qcd_sb                 ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb_err             ) ;   printf( "sf_qcd_sb_err   %s %g\n", label, sf_qcd_sb_err             ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig                ) ;   printf( "sf_qcd_sig      %s %g\n", label, sf_qcd_sig                ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig_err            ) ;   printf( "sf_qcd_sig_err  %s %g\n", label, sf_qcd_sig_err            ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig               ) ;   printf( "sf_ttwj_sig     %s %g\n", label, sf_ttwj_sig               ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig_err           ) ;   printf( "sf_ttwj_sig_err %s %g\n", label, sf_ttwj_sig_err           ) ;
       fscanf( infp, "%s %g", label, &sf_ee                     ) ;   printf( "sf_ee           %s %g\n", label, sf_ee                     ) ;
       fscanf( infp, "%s %g", label, &sf_ee_err                 ) ;   printf( "sf_ee_err       %s %g\n", label, sf_ee_err                 ) ;
       fscanf( infp, "%s %g", label, &sf_mm                     ) ;   printf( "sf_mm           %s %g\n", label, sf_mm                     ) ;
       fscanf( infp, "%s %g", label, &sf_mm_err                 ) ;   printf( "sf_mm_err       %s %g\n", label, sf_mm_err                 ) ;

       printf("\n Done reading in %s\n\n", infile ) ;
       fclose( infp ) ;

    //---- calculations for determining initial values for floating parameters.

     //-- Znunu stuff

       float initialval_znn_sig_ee(2.) ;
       float initialval_znn_sig_mm(2.) ;
       float initialval_znn_sig(2.) ;
       float initialval_znn_sb_ee(2.) ;
       float initialval_znn_sb_mm(2.) ;
       float initialval_znn_sb(2.) ;


       initialval_znn_sig_ee = (Nsig_ee) * ( 5.95 * DataLumi * fsig_ee_mean * knn_ee_sig_mean ) / ( acc_ee_sig_mean * eff_ee_mean * Ztoll_lumi ) ;
       initialval_znn_sb_ee  = (Nsb_ee ) * ( 5.95 * DataLumi * fsig_ee_mean * knn_ee_sb_mean  ) / ( acc_ee_sb_mean  * eff_ee_mean * Ztoll_lumi ) ;

       initialval_znn_sig_mm = (Nsig_mm) * ( 5.95 * DataLumi * fsig_mm_mean * knn_mm_sig_mean ) / ( acc_mm_sig_mean * eff_mm_mean * Ztoll_lumi ) ;
       initialval_znn_sb_mm  = (Nsb_mm ) * ( 5.95 * DataLumi * fsig_mm_mean * knn_mm_sb_mean  ) / ( acc_mm_sb_mean  * eff_mm_mean * Ztoll_lumi ) ;


       //-- really dumb ave.
       initialval_znn_sig = 0.5 * ( initialval_znn_sig_ee + initialval_znn_sig_mm ) ;
       initialval_znn_sb  = 0.5 * ( initialval_znn_sb_ee  + initialval_znn_sb_mm ) ;


  /////--- QCD lsb stuff

  ///  double initialval_qcd_lsb_0b     = Nhtonlytrig_lsb_0b     ;
  ///  double initialval_qcd_lsb_0b_ldp = Nhtonlytrig_lsb_0b_ldp ;


     //--- QCD ldp stuff

       double initialval_qcd_sig_ldp = Nsig_ldp - ( Nttbarmc_sig_ldp + NWJmc_sig_ldp + NZnnmc_sig_ldp + NEwomc_sig_ldp ) ;
       double initialval_qcd_sb_ldp  = Nsb_ldp  - ( Nttbarmc_sb_ldp  + NWJmc_sb_ldp  + NZnnmc_sb_ldp  + NEwomc_sb_ldp  ) ;


     //--- QCD sig and sb stuff

   /// if ( Nhtonlytrig_lsb_0b_ldp <= 0 ) {
   ///    printf("\n\n\n *** Nhtonlytrig_lsb_0b_ldp has a crazy value (%d).  I quit.\n\n\n", Nhtonlytrig_lsb_0b_ldp ) ;
   ///    return false ;
   /// }
   /// double Rldp = initialval_qcd_lsb_0b / initialval_qcd_lsb_0b_ldp ;
       double initialval_qcd_sig = Rlsb_passfail * initialval_qcd_sig_ldp ;
       double initialval_qcd_sb  = Rlsb_passfail * initialval_qcd_sb_ldp  ;


     //--- ttwj SL

       double initialval_ttwj_sig_sl = Nsig_sl ;
       double initialval_ttwj_sb_sl  = Nsb_sl ;


     //--- ttwj sig and sb

       if ( initialval_ttwj_sb_sl <= 0 ) {
          printf("\n\n\n *** initialval_ttwj_sb_sl has a crazy value (%.1f).  I quit.\n\n\n", initialval_ttwj_sb_sl ) ;
          return false ;
       }
       double initialval_ttwj_sb     = Nsb - ( initialval_qcd_sb + initialval_znn_sb ) ;
       double initialval_ttwj_sig    = initialval_ttwj_sb * ( initialval_ttwj_sig_sl / initialval_ttwj_sb_sl ) ;



       printf("\n\n\n --------- Observables and floating parameter initial values. ------------\n\n") ;

       printf("          |  Nobs   ||  ttwj  |  QCD  |  Znn  |\n") ;
       printf(" SIG      | %5d   || %5.1f | %5.1f | %5.1f |\n", Nsig, initialval_ttwj_sig, initialval_qcd_sig, initialval_znn_sig ) ;
       printf(" SB       | %5d   || %5.1f | %5.1f | %5.1f |\n", Nsb , initialval_ttwj_sb , initialval_qcd_sb , initialval_znn_sb  ) ;
       printf(" SIG,LDP  | %5d   ||*%5.1f | %5.1f |*%5.1f |\n", Nsig_ldp, (Nttbarmc_sig_ldp+NWJmc_sig_ldp+NEwomc_sig_ldp), initialval_qcd_sig_ldp, NZnnmc_sig_ldp ) ;
       printf(" SB ,LDP  | %5d   ||*%5.1f | %5.1f |*%5.1f |\n", Nsb_ldp , (Nttbarmc_sb_ldp +NWJmc_sb_ldp +NEwomc_sb_ldp) , initialval_qcd_sb_ldp , NZnnmc_sb_ldp  ) ;
       printf("\n * means fixed MC value.\n\n\n") ;




     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


       printf(" --- Defining observables.\n" ) ;


      rv_Nsig        = new RooRealVar( "Nsig"        , "Nsig"        , 0.0, 1000000. ) ;
      rv_Nsb         = new RooRealVar( "Nsb"         , "Nsb"         , 0.0, 1000000. ) ;

      rv_Nsig_sl     = new RooRealVar( "Nsig_sl"     , "Nsig_sl"     , 0.0, 1000000. ) ;
      rv_Nsb_sl      = new RooRealVar( "Nsb_sl"      , "Nsb_sl"      , 0.0, 1000000. ) ;

      rv_Nsig_ldp    = new RooRealVar( "Nsig_ldp"    , "Nsig_ldp"    , 0.0, 1000000. ) ;
      rv_Nsb_ldp     = new RooRealVar( "Nsb_ldp"     , "Nsb_ldp"     , 0.0, 1000000. ) ;

 ///  rv_Nlsb_0b     = new RooRealVar( "Nlsb_0b"     , "Nlsb_0b"     , 0.0, 1000000. ) ;
 ///  rv_Nlsb_0b_ldp = new RooRealVar( "Nlsb_0b_ldp" , "Nlsb_0b_ldp" , 0.0, 1000000. ) ;


      rv_Nsb_ee      = new RooRealVar( "Nsb_ee"      ,"Nsb_ee"       , 0., 10000. ) ;
      rv_Nsig_ee     = new RooRealVar( "Nsig_ee"     ,"Nsig_ee"      , 0., 10000. ) ;

      rv_Nsb_mm      = new RooRealVar( "Nsb_mm"      ,"Nsb_mm"       , 0., 10000. ) ;
      rv_Nsig_mm     = new RooRealVar( "Nsig_mm"     ,"Nsig_mm"      , 0., 10000. ) ;






      rv_Nsig        -> setVal( Nsig ) ;
      rv_Nsb         -> setVal( Nsb ) ;

      rv_Nsig_sl     -> setVal( Nsig_sl ) ;
      rv_Nsb_sl      -> setVal( Nsb_sl ) ;

      rv_Nsig_ldp    -> setVal( Nsig_ldp ) ;
      rv_Nsb_ldp     -> setVal( Nsb_ldp ) ;

  /// rv_Nlsb_0b     -> setVal( Nhtonlytrig_lsb_0b ) ;
  /// rv_Nlsb_0b_ldp -> setVal( Nhtonlytrig_lsb_0b_ldp ) ;


      rv_Nsb_ee      -> setVal( Nsb_ee ) ;
      rv_Nsig_ee     -> setVal( Nsig_ee ) ;

      rv_Nsb_mm      -> setVal( Nsb_mm ) ;
      rv_Nsig_mm     -> setVal( Nsig_mm ) ;








    //++++++++ Parameters of the likelihood +++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining parameters.\n" ) ;



    //____ Counts in SIG ______________________

      if ( useSigTtwjVar ) {
         rrv_mu_ttwj_sig = new RooRealVar( "mu_ttwj_sig"   , "mu_ttwj_sig"   , 0.0, 10000. ) ;
         rv_mu_ttwj_sig = rrv_mu_ttwj_sig ;
         rrv_mu_ttwj_sig   -> setVal( initialval_ttwj_sig ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sig  = new RooRealVar( "mu_qcd_sig"    , "mu_qcd_sig"    , 0.0, 200. ) ;
         rv_mu_qcd_sig = rrv_mu_qcd_sig ;
         rrv_mu_qcd_sig  -> setVal( initialval_qcd_sig ) ; //-- this is a starting value only.
      }


      //-- Note: Ewo is rfv

      rv_mu_znn_sig      = new RooRealVar( "mu_znn_sig"    , "mu_znn_sig"    , 0.0, 300. ) ;

      float maxSusySig = 4.0*Nsig ;
      rv_mu_susy_sig     = new RooRealVar( "mu_susy_sig"   , "mu_susy_sig"   , 0.0, maxSusySig ) ;


      rv_mu_znn_sig   -> setVal( initialval_znn_sig ) ;  //-- this is a starting value only.
      rv_mu_susy_sig    -> setVal( 0. ) ;  //-- this is a starting value only.





    //____ Counts in SB  ______________________


      if ( !useSigTtwjVar ) {
         rrv_mu_ttwj_sb  = new RooRealVar( "mu_ttwj_sb"    , "mu_ttwj_sb"    , 0.0, 10000. ) ;
         rv_mu_ttwj_sb = rrv_mu_ttwj_sb ;
         rrv_mu_ttwj_sb   -> setVal( initialval_ttwj_sb ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sb  = new RooRealVar( "mu_qcd_sb"    , "mu_qcd_sb"    , 0.0, 500. ) ;
         rv_mu_qcd_sb = rrv_mu_qcd_sb ;
         rrv_mu_qcd_sb  -> setVal( initialval_qcd_sb ) ; //-- this is a starting value only.
      }

      //-- Note: QCD is rfv
      //-- Note: Ewo is rfv
      //-- Note: SUSY is rfv


      rrv_mu_znn_sb       = new RooRealVar( "mu_znn_sb"     , "mu_znn_sb"     , 0.0, 350. ) ;

      rrv_mu_znn_sb   -> setVal( initialval_znn_sb ) ;  //-- this is a starting value only.

      rv_mu_znn_sb = rrv_mu_znn_sb ;
      //-- Note: Znn is rfv in Znn model 2.






    //____ Counts in SIG, SL  ______________________

      rv_mu_ttwj_sig_sl  = new RooRealVar( "mu_ttwj_sig_sl"    , "mu_ttwj_sig_sl"    , 0.0, 2500. ) ;

      //-- Note: QCD, Ewo, and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sig_sl  -> setVal( initialval_ttwj_sig_sl ) ;  //-- this is a starting value only.







    //____ Counts in SB, SL  ______________________

      rv_mu_ttwj_sb_sl  = new RooRealVar( "mu_ttwj_sb_sl"    , "mu_ttwj_sb_sl"    , 0.0, 3000. ) ;

      //-- Note: QCD, Ewo, and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sb_sl  -> setVal( initialval_ttwj_sb_sl ) ;  //-- this is a starting value only.







    //____ Counts in SIG, LDP  ______________________

      if ( useLdpVars ) {
         rrv_mu_qcd_sig_ldp  = new RooRealVar( "mu_qcd_sig_ldp"    , "mu_qcd_sig_ldp"    , 0.0, 3500. ) ;
         rv_mu_qcd_sig_ldp = rrv_mu_qcd_sig_ldp ;
         rrv_mu_qcd_sig_ldp  -> setVal( initialval_qcd_sig_ldp ) ; //-- this is a starting value only.
      }

      //-- Note: Ewo is assumed to be negligible and is not explicitly included.
      //-- Note: Znn is rfv (MC)
      //-- Note: ttwj is rfv (MC)
      //-- Note: SUSY is rfv








    //____ Counts in SB, LDP  ______________________

      if ( useLdpVars ) {
         rrv_mu_qcd_sb_ldp  = new RooRealVar( "mu_qcd_sb_ldp"    , "mu_qcd_sb_ldp"    , 0.0, 3000. ) ;
         rv_mu_qcd_sb_ldp = rrv_mu_qcd_sb_ldp ;
         rrv_mu_qcd_sb_ldp  -> setVal( initialval_qcd_sb_ldp ) ; //-- this is a starting value only.
      }


      //-- Note: Ewo is assumed to be negligible and is not explicitly included.
      //-- Note: Znn is rfv (MC)
      //-- Note: ttwj is rfv (MC)
      //-- Note: SUSY is rfv








 // //____ Counts in LSB, 0b  ______________________

 //   rv_mu_qcd_lsb_0b  = new RooRealVar( "mu_qcd_lsb_0b"    ,  "mu_qcd_lsb_0b" ,  0.    ,  10000. ) ;

 //   //-- Note: The 0btag LSB is assumed to be 100% QCD.

 //   rv_mu_qcd_lsb_0b  -> setVal( initialval_qcd_lsb_0b ) ;  //-- this is a starting value only.






 // //____ Counts in LSB, 0b, LDP  ______________________

 //   rv_mu_qcd_lsb_0b_ldp  = new RooRealVar( "mu_qcd_lsb_0b_ldp"    ,  "mu_qcd_lsb_0b_ldp",   0. ,  10000. ) ;

 //   //-- Note: The 0btag LSB, LDP is assumed to be 100% QCD.

 //   rv_mu_qcd_lsb_0b_ldp  -> setVal( initialval_qcd_lsb_0b_ldp ) ;  //-- this is a starting value only.






    //____ Counts in SB, ee  ______________________


      //-- Note: The Z to ee sample is assumed to be 100% Z to ee.
      //-- Note: zee is rfv






    //____ Counts in SIG, ee  ______________________


      //-- Note: The Z to ee sample is assumed to be 100% Z to ee.
      //-- Note: zee is rfv








    //____ Counts in SB, mm  ______________________


      //-- Note: The Z to mm sample is assumed to be 100% Z to mm.
      //-- Note: zmm is rfv






    //____ Counts in SIG, mm  ______________________


      //-- Note: The Z to mm sample is assumed to be 100% Z to mm.
      //-- Note: zmm is rfv






    //____ MC inputs _______________________________

      printf(" --- Defining MC parameters.\n" ) ;

     //-- SUSY

      rv_mu_susymc_sig      = new RooRealVar( "mu_susymc_sig"     , "mu_susymc_sig"     , 0.0, 100000. ) ;
      rv_mu_susymc_sb       = new RooRealVar( "mu_susymc_sb"      , "mu_susymc_sb"      , 0.0, 100000. ) ;
      rv_mu_susymc_sig_sl   = new RooRealVar( "mu_susymc_sig_sl"  , "mu_susymc_sig_sl"  , 0.0, 100000. ) ;
      rv_mu_susymc_sb_sl    = new RooRealVar( "mu_susymc_sb_sl"   , "mu_susymc_sb_sl"   , 0.0, 100000. ) ;
      rv_mu_susymc_sig_ldp  = new RooRealVar( "mu_susymc_sig_ldp" , "mu_susymc_sig_ldp" , 0.0, 100000. ) ;
      rv_mu_susymc_sb_ldp   = new RooRealVar( "mu_susymc_sb_ldp"  , "mu_susymc_sb_ldp"  , 0.0, 100000. ) ;

      rv_mu_susymc_sig     -> setVal( 0.1 ) ;
      rv_mu_susymc_sb      -> setVal( 0. ) ;
      rv_mu_susymc_sig_sl  -> setVal( 0. ) ;
      rv_mu_susymc_sb_sl   -> setVal( 0. ) ;
      rv_mu_susymc_sig_ldp -> setVal( 0. ) ;
      rv_mu_susymc_sb_ldp  -> setVal( 0. ) ;

      rv_mu_susymc_sig     -> setConstant(kTRUE) ;
      rv_mu_susymc_sb      -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_sl  -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_sl   -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_ldp -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_ldp  -> setConstant(kTRUE) ;



     //-- SIG, LDP

      rv_mu_ttbarmc_sig_ldp   = new RooRealVar( "mu_ttbarmc_sig_ldp" ,"mu_ttbarmc_sig_ldp" , 0., 1000. ) ;
      rv_mu_WJmc_sig_ldp      = new RooRealVar( "mu_WJmc_sig_ldp"    ,"mu_WJmc_sig_ldp"    , 0., 1000. ) ;
      rv_mu_Znnmc_sig_ldp     = new RooRealVar( "mu_Znnmc_sig_ldp"   ,"mu_Znnmc_sig_ldp"   , 0., 1000. ) ;
      rv_mu_Ewomc_sig_ldp     = new RooRealVar( "mu_Ewomc_sig_ldp"   ,"mu_Ewomc_sig_ldp"   , 0., 1000. ) ;

      rv_mu_ttbarmc_sig_ldp  -> setVal( Nttbarmc_sig_ldp ) ;
      rv_mu_WJmc_sig_ldp     -> setVal( NWJmc_sig_ldp ) ;
      rv_mu_Znnmc_sig_ldp    -> setVal( NZnnmc_sig_ldp ) ;
      rv_mu_Ewomc_sig_ldp    -> setVal( NEwomc_sig_ldp ) ;

      rv_mu_ttbarmc_sig_ldp  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sig_ldp     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sig_ldp    -> setConstant( kTRUE ) ;
      rv_mu_Ewomc_sig_ldp    -> setConstant( kTRUE ) ;


     //-- SB, LDP

      rv_mu_ttbarmc_sb_ldp   = new RooRealVar( "mu_ttbarmc_sb_ldp" ,"mu_ttbarmc_sb_ldp" , 0., 1000. ) ;
      rv_mu_WJmc_sb_ldp      = new RooRealVar( "mu_WJmc_sb_ldp"    ,"mu_WJmc_sb_ldp"    , 0., 1000. ) ;
      rv_mu_Znnmc_sb_ldp     = new RooRealVar( "mu_Znnmc_sb_ldp"   ,"mu_Znnmc_sb_ldp"   , 0., 1000. ) ;
      rv_mu_Ewomc_sb_ldp     = new RooRealVar( "mu_Ewomc_sb_ldp"   ,"mu_Ewomc_sb_ldp"   , 0., 1000. ) ;

      rv_mu_ttbarmc_sb_ldp  -> setVal( Nttbarmc_sb_ldp ) ;
      rv_mu_WJmc_sb_ldp     -> setVal( NWJmc_sb_ldp ) ;
      rv_mu_Znnmc_sb_ldp    -> setVal( NZnnmc_sb_ldp ) ;
      rv_mu_Ewomc_sb_ldp    -> setVal( NEwomc_sb_ldp ) ;

      rv_mu_ttbarmc_sb_ldp  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sb_ldp     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sb_ldp    -> setConstant( kTRUE ) ;
      rv_mu_Ewomc_sb_ldp    -> setConstant( kTRUE ) ;


     //-- SIG

      rv_mu_Ewomc_sig     = new RooRealVar( "mu_Ewomc_sig"   ,"mu_Ewomc_sig"   , 0., 1000. ) ;
      rv_mu_Ewomc_sig    -> setVal( NEwomc_sig ) ;
      rv_mu_Ewomc_sig    -> setConstant( kTRUE ) ;

     //-- SB

      rv_mu_Ewomc_sb     = new RooRealVar( "mu_Ewomc_sb"   ,"mu_Ewomc_sb"   , 0., 1000. ) ;
      rv_mu_Ewomc_sb    -> setVal( NEwomc_sb ) ;
      rv_mu_Ewomc_sb    -> setConstant( kTRUE ) ;





    //+++++++ Gaussian constraints ++++++++++++++++++++++++++++++++

      printf(" --- Defining Gaussian constraint and constant parameters.\n" ) ;

    //_______ Efficiency scale factor.  Applied to SUSY and all MC inputs _______________

    //   August 23, 2011:  switching to correlated log-normal PDFs
    //
    //
    //  double pmin, pmax ;


    //--- mean parameters.
      rv_mean_eff_sf_sig     = new RooRealVar( "mean_eff_sf_sig"    , "mean_eff_sf_sig", 0., 10. ) ;
      rv_mean_eff_sf_sb      = new RooRealVar( "mean_eff_sf_sb"     , "mean_eff_sf_sb", 0., 10. ) ;
      rv_mean_eff_sf_sig_sl  = new RooRealVar( "mean_eff_sf_sig_sl" , "mean_eff_sf_sig_sl", 0., 10. ) ;
      rv_mean_eff_sf_sb_sl   = new RooRealVar( "mean_eff_sf_sb_sl"  , "mean_eff_sf_sb_sl", 0., 10. ) ;
      rv_mean_eff_sf_sig_ldp = new RooRealVar( "mean_eff_sf_sig_ldp", "mean_eff_sf_sig_ldp", 0., 10. ) ;
      rv_mean_eff_sf_sb_ldp  = new RooRealVar( "mean_eff_sf_sb_ldp" , "mean_eff_sf_sb_ldp", 0., 10. ) ;

    //--- width parameters.
      rv_width_eff_sf_sig     = new RooRealVar( "width_eff_sf_sig"    , "width_eff_sf_sig", 0., 10. ) ;
      rv_width_eff_sf_sb      = new RooRealVar( "width_eff_sf_sb"     , "width_eff_sf_sb", 0., 10. ) ;
      rv_width_eff_sf_sig_sl  = new RooRealVar( "width_eff_sf_sig_sl" , "width_eff_sf_sig_sl", 0., 10. ) ;
      rv_width_eff_sf_sb_sl   = new RooRealVar( "width_eff_sf_sb_sl"  , "width_eff_sf_sb_sl", 0., 10. ) ;
      rv_width_eff_sf_sig_ldp = new RooRealVar( "width_eff_sf_sig_ldp", "width_eff_sf_sig_ldp", 0., 10. ) ;
      rv_width_eff_sf_sb_ldp  = new RooRealVar( "width_eff_sf_sb_ldp" , "width_eff_sf_sb_ldp", 0., 10. ) ;


      rv_mean_eff_sf_sig     -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb      -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_sl  -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_sl   -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_ldp -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_ldp  -> setVal( 1.00 ) ;

      rv_mean_eff_sf_sig     -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb      -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_sl  -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_sl   -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_ldp -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_ldp  -> setConstant( kTRUE ) ;


      //--- Initialize all width parameters to 15%.  They will be reset in the susy scan.
      rv_width_eff_sf_sig     -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb      -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_sl  -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_sl   -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_ldp -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_ldp  -> setVal( 0.15 ) ;

      rv_width_eff_sf_sig     -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb      -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_sl  -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_sl   -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_ldp -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_ldp  -> setConstant( kTRUE ) ;







   //
   // Owen : Sept 7, 2011: need to set these here so that it gets into the workspace.
   //
   //++++++++++
      setSusyScanPoint( inputScanFile,  m0,  m12,  isT1bbbb,  t1bbbbXsec ) ;
   //++++++++++










    //--- Systematics and other nuisance parameters
      char formula[1024];
      RooArgSet globalObservables ("globalObservables");
      RooArgSet allNuisances ("allNuisances");
      RooArgSet allNuisancePdfs ("allNuisancePdfs");

      RooRealVar Rlsb_passfail_prim ( "Rlsb_passfail_prim", "Rlsb_passfail_prim", 0, -5, 5);
      RooRealVar Rlsb_passfail_nom ( "Rlsb_passfail_nom", "Rlsb_passfail_nom", 0, -5, 5);
      RooGaussian pdf_Rlsb_passfail ("pdf_Rlsb_passfail" , "pdf_Rlsb_passfail", Rlsb_passfail_prim, Rlsb_passfail_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", Rlsb_passfail, exp(Rlsb_passfail_err/Rlsb_passfail));
      RooFormulaVar fv_Rlsb_passfail ("Rlsb_passfail", formula, RooArgList(Rlsb_passfail_prim));
      Rlsb_passfail_nom.setConstant();
      globalObservables.add (Rlsb_passfail_nom);
      allNuisances.add (Rlsb_passfail_prim);
      allNuisancePdfs.add (pdf_Rlsb_passfail);

      RooRealVar acc_ee_sig_prim ( "acc_ee_sig_prim", "acc_ee_sig_prim", 0, -5, 5);
      RooRealVar acc_ee_sig_nom ( "acc_ee_sig_nom", "acc_ee_sig_nom", 0, -5, 5);
      RooGaussian pdf_acc_ee_sig ("pdf_acc_ee_sig" , "pdf_acc_ee_sig", acc_ee_sig_prim, acc_ee_sig_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", acc_ee_sig_mean, exp(acc_ee_sig_err/acc_ee_sig_mean));
      RooFormulaVar fv_acc_ee_sig ("acc_ee_sig", formula, RooArgList(acc_ee_sig_prim));
      acc_ee_sig_nom.setConstant();
      globalObservables.add (acc_ee_sig_nom);
      allNuisances.add (acc_ee_sig_prim);
      allNuisancePdfs.add (pdf_acc_ee_sig);

      RooRealVar acc_ee_sb_prim ( "acc_ee_sb_prim", "acc_ee_sb_prim", 0, -5, 5);
      RooRealVar acc_ee_sb_nom ( "acc_ee_sb_nom", "acc_ee_sb_nom", 0, -5, 5);
      RooGaussian pdf_acc_ee_sb ("pdf_acc_ee_sb" , "pdf_acc_ee_sb", acc_ee_sb_prim, acc_ee_sb_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", acc_ee_sb_mean, exp(acc_ee_sb_err/acc_ee_sb_mean));
      RooFormulaVar fv_acc_ee_sb ("acc_ee_sb", formula, RooArgList(acc_ee_sb_prim));
      acc_ee_sb_nom.setConstant();
      globalObservables.add (acc_ee_sb_nom);
      allNuisances.add (acc_ee_sb_prim);
      allNuisancePdfs.add (pdf_acc_ee_sb);

      RooRealVar acc_mm_sig_prim ( "acc_mm_sig_prim", "acc_mm_sig_prim", 0, -5, 5);
      RooRealVar acc_mm_sig_nom ( "acc_mm_sig_nom", "acc_mm_sig_nom", 0, -5, 5);
      RooGaussian pdf_acc_mm_sig ("pdf_acc_mm_sig" , "pdf_acc_mm_sig", acc_mm_sig_prim, acc_mm_sig_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", acc_mm_sig_mean, exp(acc_mm_sig_err/acc_mm_sig_mean));
      RooFormulaVar fv_acc_mm_sig ("acc_mm_sig", formula, RooArgList(acc_mm_sig_prim));
      acc_mm_sig_nom.setConstant();
      globalObservables.add (acc_mm_sig_nom);
      allNuisances.add (acc_mm_sig_prim);
      allNuisancePdfs.add (pdf_acc_mm_sig);

      RooRealVar acc_mm_sb_prim ( "acc_mm_sb_prim", "acc_mm_sb_prim", 0, -5, 5);
      RooRealVar acc_mm_sb_nom ( "acc_mm_sb_nom", "acc_mm_sb_nom", 0, -5, 5);
      RooGaussian pdf_acc_mm_sb ("pdf_acc_mm_sb" , "pdf_acc_mm_sb", acc_mm_sb_prim, acc_mm_sb_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", acc_mm_sb_mean, exp(acc_mm_sb_err/acc_mm_sb_mean));
      RooFormulaVar fv_acc_mm_sb ("acc_mm_sb", formula, RooArgList(acc_mm_sb_prim));
      acc_mm_sb_nom.setConstant();
      globalObservables.add (acc_mm_sb_nom);
      allNuisances.add (acc_mm_sb_prim);
      allNuisancePdfs.add (pdf_acc_mm_sb);


      RooRealVar eff_ee_prim ( "eff_ee_prim", "eff_ee_prim", 0, -5, 5);
      RooRealVar eff_ee_nom ( "eff_ee_nom", "eff_ee_nom", 0, -5, 5);
      RooGaussian pdf_eff_ee ("pdf_eff_ee" , "pdf_eff_ee", eff_ee_prim, eff_ee_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", eff_ee_mean, exp(eff_ee_err/eff_ee_mean));
      RooFormulaVar fv_eff_ee ("eff_ee", formula, RooArgList(eff_ee_prim));
      eff_ee_nom.setConstant();
      globalObservables.add (eff_ee_nom);
      allNuisances.add (eff_ee_prim);
      allNuisancePdfs.add (pdf_eff_ee);

      RooRealVar eff_mm_prim ( "eff_mm_prim", "eff_mm_prim", 0, -5, 5);
      RooRealVar eff_mm_nom ( "eff_mm_nom", "eff_mm_nom", 0, -5, 5);
      RooGaussian pdf_eff_mm ("pdf_eff_mm" , "pdf_eff_mm", eff_mm_prim, eff_mm_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", eff_mm_mean, exp(eff_mm_err/eff_mm_mean));
      RooFormulaVar fv_eff_mm ("eff_mm", formula, RooArgList(eff_mm_prim));
      eff_mm_nom.setConstant();
      globalObservables.add (eff_mm_nom);
      allNuisances.add (eff_mm_prim);
      allNuisancePdfs.add (pdf_eff_mm);

      RooRealVar fsig_ee_prim ( "fsig_ee_prim", "fsig_ee_prim", 0, -5, 5);
      RooRealVar fsig_ee_nom ( "fsig_ee_nom", "fsig_ee_nom", 0, -5, 5);
      RooGaussian pdf_fsig_ee ("pdf_fsig_ee" , "pdf_fsig_ee", fsig_ee_prim, fsig_ee_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", fsig_ee_mean, exp(fsig_ee_err/fsig_ee_mean));
      RooFormulaVar fv_fsig_ee ("fsig_ee", formula, RooArgList(fsig_ee_prim));
      fsig_ee_nom.setConstant();
      globalObservables.add (fsig_ee_nom);
      allNuisances.add (fsig_ee_prim);
      allNuisancePdfs.add (pdf_fsig_ee);

      RooRealVar fsig_mm_prim ( "fsig_mm_prim", "fsig_mm_prim", 0, -5, 5);
      RooRealVar fsig_mm_nom ( "fsig_mm_nom", "fsig_mm_nom", 0, -5, 5);
      RooGaussian pdf_fsig_mm ("pdf_fsig_mm" , "pdf_fsig_mm", fsig_mm_prim, fsig_mm_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", fsig_mm_mean, exp(fsig_mm_err/fsig_mm_mean));
      RooFormulaVar fv_fsig_mm ("fsig_mm", formula, RooArgList(fsig_mm_prim));
      fsig_mm_nom.setConstant();
      globalObservables.add (fsig_mm_nom);
      allNuisances.add (fsig_mm_prim);
      allNuisancePdfs.add (pdf_fsig_mm);


    //--- Nov 24, 2011:  sf_ee and sf_mm are now derived from a common underlying Gaussian.

      RooRealVar sf_ll_prim ( "sf_ll_prim", "sf_ll_prim", 0, -5, 5);
      RooRealVar sf_ll_nom ( "sf_ll_nom", "sf_ll_nom", 0, -5, 5);
      RooGaussian pdf_sf_ll ("pdf_sf_ll" , "pdf_sf_ll", sf_ll_prim, sf_ll_nom, RooConst(1));
      sf_ll_nom.setConstant();
      globalObservables.add (sf_ll_nom);
      allNuisances.add (sf_ll_prim);
      allNuisancePdfs.add (pdf_sf_ll);

      sprintf (formula, "%f*pow(%f,@0)", sf_ee, exp(sf_ee_err/sf_ee));
      RooFormulaVar fv_sf_ee ("sf_ee", formula, RooArgList(sf_ll_prim));

      sprintf (formula, "%f*pow(%f,@0)", sf_mm, exp(sf_mm_err/sf_mm));
      RooFormulaVar fv_sf_mm ("sf_mm", formula, RooArgList(sf_ll_prim));





      RooRealVar sf_mc_prim ( "sf_mc_prim", "sf_mc_prim", 0, -5, 5);
      RooRealVar sf_mc_nom ( "sf_mc_nom", "sf_mc_nom", 0, -5, 5);
      RooGaussian pdf_sf_mc ("pdf_sf_mc" , "pdf_sf_mc", sf_mc_prim, sf_mc_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_mc, exp(sf_mc_err/sf_mc));
      RooFormulaVar fv_sf_mc ("sf_mc", formula, RooArgList(sf_mc_prim));
      sf_mc_nom.setConstant();
      globalObservables.add (sf_mc_nom);
      allNuisances.add (sf_mc_prim);
      allNuisancePdfs.add (pdf_sf_mc);

      RooRealVar sf_qcd_sb_prim ( "sf_qcd_sb_prim", "sf_qcd_sb_prim", 0, -5, 5);
      RooRealVar sf_qcd_sb_nom ( "sf_qcd_sb_nom", "sf_qcd_sb_nom", 0, -5, 5);
      RooGaussian pdf_sf_qcd_sb ("pdf_sf_qcd_sb" , "pdf_sf_qcd_sb", sf_qcd_sb_prim, sf_qcd_sb_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_qcd_sb, exp(sf_qcd_sb_err/sf_qcd_sb));
      RooFormulaVar fv_sf_qcd_sb ("sf_qcd_sb", formula, RooArgList(sf_qcd_sb_prim));
      sf_qcd_sb_nom.setConstant();
      globalObservables.add (sf_qcd_sb_nom);
      allNuisances.add (sf_qcd_sb_prim);
      allNuisancePdfs.add (pdf_sf_qcd_sb);

      RooRealVar sf_qcd_sig_prim ( "sf_qcd_sig_prim", "sf_qcd_sig_prim", 0, -5, 5);
      RooRealVar sf_qcd_sig_nom ( "sf_qcd_sig_nom", "sf_qcd_sig_nom", 0, -5, 5);
      RooGaussian pdf_sf_qcd_sig ("pdf_sf_qcd_sig" , "pdf_sf_qcd_sig", sf_qcd_sig_prim, sf_qcd_sig_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_qcd_sig, exp(sf_qcd_sig_err/sf_qcd_sig));
      RooFormulaVar fv_sf_qcd_sig ("sf_qcd_sig", formula, RooArgList(sf_qcd_sig_prim));
      sf_qcd_sig_nom.setConstant();
      globalObservables.add (sf_qcd_sig_nom);
      allNuisances.add (sf_qcd_sig_prim);
      allNuisancePdfs.add (pdf_sf_qcd_sig);

      RooRealVar sf_ttwj_sig_prim ( "sf_ttwj_sig_prim", "sf_ttwj_sig_prim", 0, -5, 5);
      RooRealVar sf_ttwj_sig_nom ( "sf_ttwj_sig_nom", "sf_ttwj_sig_nom", 0, -5, 5);
      RooGaussian pdf_sf_ttwj_sig ("pdf_sf_ttwj_sig" , "pdf_sf_ttwj_sig", sf_ttwj_sig_prim, sf_ttwj_sig_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_ttwj_sig, exp(sf_ttwj_sig_err/sf_ttwj_sig));
      RooFormulaVar fv_sf_ttwj_sig ("sf_ttwj_sig", formula, RooArgList(sf_ttwj_sig_prim));
      sf_ttwj_sig_nom.setConstant();
      globalObservables.add (sf_ttwj_sig_nom);
      allNuisances.add (sf_ttwj_sig_prim);
      allNuisancePdfs.add (pdf_sf_ttwj_sig);


   //--------------------------------------

      RooRealVar knn_ee_sig_prim ( "knn_ee_sig_prim", "knn_ee_sig_prim", 0, -5, 5);
      RooRealVar knn_ee_sig_nom ( "knn_ee_sig_nom", "knn_ee_sig_nom", 0, -5, 5);
      RooGaussian pdf_knn_ee_sig ("pdf_knn_ee_sig" , "pdf_knn_ee_sig", knn_ee_sig_prim, knn_ee_sig_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_ee_sig_mean, exp(knn_ee_sig_err/knn_ee_sig_mean));
      RooFormulaVar fv_knn_ee_sig ("knn_ee_sig", formula, RooArgList(knn_ee_sig_prim));
      knn_ee_sig_nom.setConstant();
      globalObservables.add (knn_ee_sig_nom);
      allNuisances.add (knn_ee_sig_prim);
      allNuisancePdfs.add (pdf_knn_ee_sig);

      RooRealVar knn_ee_sb_prim ( "knn_ee_sb_prim", "knn_ee_sb_prim", 0, -5, 5);
      RooRealVar knn_ee_sb_nom ( "knn_ee_sb_nom", "knn_ee_sb_nom", 0, -5, 5);
      RooGaussian pdf_knn_ee_sb ("pdf_knn_ee_sb" , "pdf_knn_ee_sb", knn_ee_sb_prim, knn_ee_sb_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_ee_sb_mean, exp(knn_ee_sb_err/knn_ee_sb_mean));
      RooFormulaVar fv_knn_ee_sb ("knn_ee_sb", formula, RooArgList(knn_ee_sb_prim));
      knn_ee_sb_nom.setConstant();
      globalObservables.add (knn_ee_sb_nom);
      allNuisances.add (knn_ee_sb_prim);
      allNuisancePdfs.add (pdf_knn_ee_sb);

      RooRealVar knn_mm_sig_prim ( "knn_mm_sig_prim", "knn_mm_sig_prim", 0, -5, 5);
      RooRealVar knn_mm_sig_nom ( "knn_mm_sig_nom", "knn_mm_sig_nom", 0, -5, 5);
      RooGaussian pdf_knn_mm_sig ("pdf_knn_mm_sig" , "pdf_knn_mm_sig", knn_mm_sig_prim, knn_mm_sig_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_mm_sig_mean, exp(knn_mm_sig_err/knn_mm_sig_mean));
      RooFormulaVar fv_knn_mm_sig ("knn_mm_sig", formula, RooArgList(knn_mm_sig_prim));
      knn_mm_sig_nom.setConstant();
      globalObservables.add (knn_mm_sig_nom);
      allNuisances.add (knn_mm_sig_prim);
      allNuisancePdfs.add (pdf_knn_mm_sig);

      RooRealVar knn_mm_sb_prim ( "knn_mm_sb_prim", "knn_mm_sb_prim", 0, -5, 5);
      RooRealVar knn_mm_sb_nom ( "knn_mm_sb_nom", "knn_mm_sb_nom", 0, -5, 5);
      RooGaussian pdf_knn_mm_sb ("pdf_knn_mm_sb" , "pdf_knn_mm_sb", knn_mm_sb_prim, knn_mm_sb_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_mm_sb_mean, exp(knn_mm_sb_err/knn_mm_sb_mean));
      RooFormulaVar fv_knn_mm_sb ("knn_mm_sb", formula, RooArgList(knn_mm_sb_prim));
      knn_mm_sb_nom.setConstant();
      globalObservables.add (knn_mm_sb_nom);
      allNuisances.add (knn_mm_sb_prim);
      allNuisancePdfs.add (pdf_knn_mm_sb);

    //--------------------------------------


      RooRealVar eff_sf_prim ( "eff_sf_prim", "eff_sf_prim", 0, -5, 5);
      RooRealVar eff_sf_nom ( "eff_sf_nom", "eff_sf_nom", 0, -5, 5);
      RooGaussian pdf_eff_sf_prim ("pdf_eff_sf_prim" , "master pdf_eff_sf_prim", eff_sf_prim, eff_sf_nom, RooConst(1));
      eff_sf_nom.setConstant();
      globalObservables.add (eff_sf_nom);
      allNuisances.add (eff_sf_prim);
      allNuisancePdfs.add (pdf_eff_sf_prim);

/*       RooRealVar xxx_prim ( "xxx_prim", "xxx_prim", 0, -5, 5); */
/*       RooRealVar xxx_nom ( "xxx_nom", "xxx_nom", 0, -5, 5); */
/*       RooGaussian pdf_xxx ("pdf_xxx" , "pdf_xxx", *xxx_prim, *xxx_nom, RooConst(1)); */
/*       sprintf (formula, "%f*pow(%f,@0)", xxx, exp(xxx_err/xxx)); */
/*       RooFormulaVar fv_xxx ("xxx", formula, RooArgList(xxx_prim)); */
/*       globalObservables.add (xxx_nom); */
/*       allNuisances.add (xxx_prim); */
/*       allNuisancePdfs.add (pdf_xxx); */



    //-- Z to nunu stuff



      rv_znnoverll_bfratio = new RooRealVar( "znnoverll_bfratio", "znnoverll_bfratio", 0.1, 10. ) ;
      rv_znnoverll_bfratio -> setVal( 5.95 ) ;
      rv_znnoverll_bfratio -> setConstant( kTRUE ) ;

      rv_dataoverll_lumiratio = new RooRealVar( "dataoverll_lumiratio", "dataoverll_lumiratio", 0.1, 10.0 ) ;
      rv_dataoverll_lumiratio  -> setVal( DataLumi / Ztoll_lumi ) ;
      rv_dataoverll_lumiratio  -> setConstant( kTRUE ) ;




     //+++++++++++++++++ Relationships between parameters ++++++++++++++++++++++++++++++++++++++++++++

       printf(" --- Defining relationships between parameters.\n" ) ;




    //-- ttwj

      if ( useSigTtwjVar ) {
         rfv_mu_ttwj_sb = new RooFormulaVar( "mu_ttwj_sb",
                                                       "mu_ttwj_sig * (1.0/sf_ttwj_sig) * (mu_ttwj_sb_sl/mu_ttwj_sig_sl)",
                                                       RooArgSet( *rv_mu_ttwj_sig, fv_sf_ttwj_sig, *rv_mu_ttwj_sb_sl, *rv_mu_ttwj_sig_sl ) ) ;
         rv_mu_ttwj_sb = rfv_mu_ttwj_sb ;
      } else {
         rfv_mu_ttwj_sig = new RooFormulaVar( "mu_ttwj_sig",
                                                       "mu_ttwj_sb * sf_ttwj_sig * (mu_ttwj_sig_sl/mu_ttwj_sb_sl)",
                                                       RooArgSet( *rv_mu_ttwj_sb, fv_sf_ttwj_sig, *rv_mu_ttwj_sig_sl, *rv_mu_ttwj_sb_sl ) ) ;
         rv_mu_ttwj_sig = rfv_mu_ttwj_sig ;
      }




      rv_mu_ttwj_sig_ldp = new RooFormulaVar( "mu_ttwj_sig_ldp",
                              "mu_ttbarmc_sig_ldp + mu_WJmc_sig_ldp",
                              RooArgSet( *rv_mu_ttbarmc_sig_ldp,
                                         *rv_mu_WJmc_sig_ldp ) ) ;


      rv_mu_ttwj_sb_ldp = new RooFormulaVar( "mu_ttwj_sb_ldp",
                              "mu_ttbarmc_sb_ldp + mu_WJmc_sb_ldp",
                              RooArgSet( *rv_mu_ttbarmc_sb_ldp,
                                         *rv_mu_WJmc_sb_ldp ) ) ;






    //-- QCD

   //-------------------------------
   ///if ( useLdpVars ) {

   ///   rfv_mu_qcd_sig = new RooFormulaVar( "mu_qcd_sig",
   ///                               "mu_qcd_sig_ldp * sf_qcd_sig * ( mu_qcd_lsb_0b / mu_qcd_lsb_0b_ldp )",
   ///                               RooArgSet( *rv_mu_qcd_sig_ldp, fv_sf_qcd_sig, *rv_mu_qcd_lsb_0b, *rv_mu_qcd_lsb_0b_ldp ) ) ;
   ///   rv_mu_qcd_sig = rfv_mu_qcd_sig ;

   ///   rfv_mu_qcd_sb = new RooFormulaVar( "mu_qcd_sb",
   ///                               "mu_qcd_sb_ldp * sf_qcd_sb * ( mu_qcd_lsb_0b / mu_qcd_lsb_0b_ldp )",
   ///                               RooArgSet( *rv_mu_qcd_sb_ldp, fv_sf_qcd_sb, *rv_mu_qcd_lsb_0b, *rv_mu_qcd_lsb_0b_ldp ) ) ;
   ///   rv_mu_qcd_sb = rfv_mu_qcd_sb ;

   ///} else {

   ///   rfv_mu_qcd_sig_ldp = new RooFormulaVar( "mu_qcd_sig_ldp",
   ///                               "mu_qcd_sig * (1.0/sf_qcd_sig) * ( mu_qcd_lsb_0b_ldp / mu_qcd_lsb_0b )",
   ///                               RooArgSet( *rv_mu_qcd_sig, fv_sf_qcd_sig, *rv_mu_qcd_lsb_0b_ldp, *rv_mu_qcd_lsb_0b ) ) ;
   ///   rv_mu_qcd_sig_ldp = rfv_mu_qcd_sig_ldp ;

   ///   rfv_mu_qcd_sb_ldp = new RooFormulaVar( "mu_qcd_sb_ldp",
   ///                               "mu_qcd_sb * (1.0/sf_qcd_sb) * ( mu_qcd_lsb_0b_ldp / mu_qcd_lsb_0b )",
   ///                               RooArgSet( *rv_mu_qcd_sb, fv_sf_qcd_sb, *rv_mu_qcd_lsb_0b_ldp, *rv_mu_qcd_lsb_0b ) ) ;
   ///   rv_mu_qcd_sb_ldp = rfv_mu_qcd_sb_ldp ;

   ///}
   //-------------------------------

      if ( useLdpVars ) {

         rfv_mu_qcd_sig = new RooFormulaVar( "mu_qcd_sig",
                                     "mu_qcd_sig_ldp * sf_qcd_sig * Rlsb_passfail",
                                     RooArgSet( *rv_mu_qcd_sig_ldp, fv_sf_qcd_sig, fv_Rlsb_passfail ) ) ;
         rv_mu_qcd_sig = rfv_mu_qcd_sig ;

         rfv_mu_qcd_sb = new RooFormulaVar( "mu_qcd_sb",
                                     "mu_qcd_sb_ldp * sf_qcd_sb * Rlsb_passfail",
                                     RooArgSet( *rv_mu_qcd_sb_ldp, fv_sf_qcd_sb, fv_Rlsb_passfail ) ) ;
         rv_mu_qcd_sb = rfv_mu_qcd_sb ;

      } else {

         rfv_mu_qcd_sig_ldp = new RooFormulaVar( "mu_qcd_sig_ldp",
                                     "mu_qcd_sig * (1.0/sf_qcd_sig) * ( 1.0 / Rlsb_passfail )",
                                     RooArgSet( *rv_mu_qcd_sig, fv_sf_qcd_sig, fv_Rlsb_passfail ) ) ;
         rv_mu_qcd_sig_ldp = rfv_mu_qcd_sig_ldp ;

         rfv_mu_qcd_sb_ldp = new RooFormulaVar( "mu_qcd_sb_ldp",
                                     "mu_qcd_sb * (1.0/sf_qcd_sb) * ( 1.0 / Rlsb_passfail )",
                                     RooArgSet( *rv_mu_qcd_sb, fv_sf_qcd_sb, fv_Rlsb_passfail ) ) ;
         rv_mu_qcd_sb_ldp = rfv_mu_qcd_sb_ldp ;

      }

   //-------------------------------



    //-- SUSY

      rv_mu_susy_sb = new RooFormulaVar( "mu_susy_sb",
                                        "mu_susymc_sb * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sb, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sig_sl = new RooFormulaVar( "mu_susy_sig_sl",
                                        "mu_susymc_sig_sl * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sig_sl, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sb_sl = new RooFormulaVar( "mu_susy_sb_sl",
                                        "mu_susymc_sb_sl * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sb_sl, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sig_ldp = new RooFormulaVar( "mu_susy_sig_ldp",
                                        "mu_susymc_sig_ldp * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sig_ldp, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sb_ldp = new RooFormulaVar( "mu_susy_sb_ldp",
                                        "mu_susymc_sb_ldp * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sb_ldp, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;





    //-- Z to nu nu

      rv_mu_znn_sig_ldp = new RooFormulaVar( "mu_znn_sig_ldp",
                              "mu_Znnmc_sig_ldp",
                              RooArgSet( *rv_mu_Znnmc_sig_ldp ) ) ;

      rv_mu_znn_sb_ldp = new RooFormulaVar( "mu_znn_sb_ldp",
                              "mu_Znnmc_sb_ldp",
                              RooArgSet( *rv_mu_Znnmc_sb_ldp ) ) ;


   //---------------------------------
   // if ( znnModel == 1 ) {

   //    rv_mu_zee_sb_ee = new RooFormulaVar( "mu_zee_sb_ee",
   //                                    "mu_znn_sb * sf_ee * ( acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
   //                                    RooArgSet( *rv_mu_znn_sb, fv_sf_ee, fv_acc_ee, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   //    rv_mu_zee_sig_ee = new RooFormulaVar( "mu_zee_sig_ee",
   //                                    "mu_znn_sig * sf_ee * ( acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
   //                                    RooArgSet( *rv_mu_znn_sig, fv_sf_ee, fv_acc_ee, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   //    rv_mu_zmm_sb_mm = new RooFormulaVar( "mu_zmm_sb_mm",
   //                                    "mu_znn_sb * sf_mm * ( acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
   //                                    RooArgSet( *rv_mu_znn_sb, fv_sf_mm, fv_acc_mm, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   //    rv_mu_zmm_sig_mm = new RooFormulaVar( "mu_zmm_sig_mm",
   //                                    "mu_znn_sig * sf_mm * ( acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
   //                                    RooArgSet( *rv_mu_znn_sig, fv_sf_mm, fv_acc_mm, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   // } else if ( znnModel == 2 ) {

   //    rfv_mu_znn_sb = new RooFormulaVar( "mu_znn_sb",
   //                                     "mu_znn_sig * ( knn_sb / knn_sig )",
   //                                     RooArgSet( *rv_mu_znn_sig, fv_knn_sb, fv_knn_sig ) ) ;

   //    rv_mu_znn_sb = rfv_mu_znn_sb ;

   //    rv_mu_zee_sigsb_ee = new RooFormulaVar( "mu_zee_sigsb_ee",
   //                                  "( mu_znn_sig / knn_sig ) * sf_ee * ( (acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
   //                                    RooArgSet( *rv_mu_znn_sig, fv_knn_sig, fv_sf_ee, fv_acc_ee, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   //    rv_mu_zmm_sigsb_mm = new RooFormulaVar( "mu_zmm_sigsb_mm",
   //                                  "( mu_znn_sig / knn_sig ) * sf_mm * ( (acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
   //                                    RooArgSet( *rv_mu_znn_sig, fv_knn_sig, fv_sf_mm, fv_acc_mm, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   // }
   //---------------------------------



      rv_mu_zee_sig_ee = new RooFormulaVar( "mu_zee_sig_ee",
                                    "( mu_znn_sig / knn_ee_sig ) * sf_ee * ( (acc_ee_sig * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                      RooArgSet( *rv_mu_znn_sig, fv_knn_ee_sig, fv_sf_ee, fv_acc_ee_sig, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

      rv_mu_zee_sb_ee = new RooFormulaVar( "mu_zee_sb_ee",
                                    "( mu_znn_sb / knn_ee_sb ) * sf_ee * ( (acc_ee_sb * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                      RooArgSet( *rv_mu_znn_sb, fv_knn_ee_sb, fv_sf_ee, fv_acc_ee_sb, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

      rv_mu_zmm_sig_mm = new RooFormulaVar( "mu_zmm_sig_mm",
                                    "( mu_znn_sig / knn_mm_sig ) * sf_mm * ( (acc_mm_sig * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                      RooArgSet( *rv_mu_znn_sig, fv_knn_mm_sig, fv_sf_mm, fv_acc_mm_sig, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

      rv_mu_zmm_sb_mm = new RooFormulaVar( "mu_zmm_sb_mm",
                                    "( mu_znn_sb / knn_mm_sb ) * sf_mm * ( (acc_mm_sb * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                      RooArgSet( *rv_mu_znn_sb, fv_knn_mm_sb, fv_sf_mm, fv_acc_mm_sb, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;










    //-- EWO

      rv_mu_ewo_sig     = new RooFormulaVar( "mu_ewo_sig"     , "mu_Ewomc_sig"     , RooArgSet( *rv_mu_Ewomc_sig     ) ) ;
      rv_mu_ewo_sb      = new RooFormulaVar( "mu_ewo_sb"      , "mu_Ewomc_sb"      , RooArgSet( *rv_mu_Ewomc_sb      ) ) ;
      rv_mu_ewo_sig_ldp = new RooFormulaVar( "mu_ewo_sig_ldp" , "mu_Ewomc_sig_ldp" , RooArgSet( *rv_mu_Ewomc_sig_ldp ) ) ;
      rv_mu_ewo_sb_ldp  = new RooFormulaVar( "mu_ewo_sb_ldp"  , "mu_Ewomc_sb_ldp"  , RooArgSet( *rv_mu_Ewomc_sb_ldp  ) ) ;




    //-- Parametric relations between correlated efficiency scale factors.



      rv_eff_sf_sig = new RooFormulaVar( "eff_sf_sig",
                                         "mean_eff_sf_sig * pow( exp( width_eff_sf_sig/mean_eff_sf_sig ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig, *rv_width_eff_sf_sig, *rv_mean_eff_sf_sig, eff_sf_prim ) ) ;

      rv_eff_sf_sb = new RooFormulaVar( "eff_sf_sb",
                                         "mean_eff_sf_sb * pow( exp( width_eff_sf_sb/mean_eff_sf_sb ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb, *rv_width_eff_sf_sb, *rv_mean_eff_sf_sb, eff_sf_prim ) ) ;

      rv_eff_sf_sig_sl = new RooFormulaVar( "eff_sf_sig_sl",
                                         "mean_eff_sf_sig_sl * pow( exp( width_eff_sf_sig_sl/mean_eff_sf_sig_sl ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_sl, *rv_width_eff_sf_sig_sl, *rv_mean_eff_sf_sig_sl, eff_sf_prim ) ) ;

      rv_eff_sf_sb_sl = new RooFormulaVar( "eff_sf_sb_sl",
                                         "mean_eff_sf_sb_sl * pow( exp( width_eff_sf_sb_sl/mean_eff_sf_sb_sl ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_sl, *rv_width_eff_sf_sb_sl, *rv_mean_eff_sf_sb_sl, eff_sf_prim ) ) ;

      rv_eff_sf_sig_ldp = new RooFormulaVar( "eff_sf_sig_ldp",
                                         "mean_eff_sf_sig_ldp * pow( exp( width_eff_sf_sig_ldp/mean_eff_sf_sig_ldp ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_ldp, *rv_width_eff_sf_sig_ldp, *rv_mean_eff_sf_sig_ldp, eff_sf_prim ) ) ;

      rv_eff_sf_sb_ldp = new RooFormulaVar( "eff_sf_sb_ldp",
                                         "mean_eff_sf_sb_ldp * pow( exp( width_eff_sf_sb_ldp/mean_eff_sf_sb_ldp ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_ldp, *rv_width_eff_sf_sb_ldp, *rv_mean_eff_sf_sb_ldp, eff_sf_prim ) ) ;


    //+++++++++++++ Expected counts for observables in terms of parameters ++++++++++++++++++

       printf(" --- Defining expected counts in terms of parameters.\n" ) ;


      rv_n_sig         = new RooFormulaVar( "n_sig",
                                     "mu_ttwj_sig + mu_qcd_sig + mu_znn_sig + eff_sf_sig*( mu_ewo_sig + mu_susy_sig)",
                                     RooArgSet( *rv_mu_ttwj_sig, *rv_mu_qcd_sig, *rv_mu_znn_sig, *rv_eff_sf_sig, *rv_mu_ewo_sig, *rv_mu_susy_sig ) ) ;

      rv_n_sb          = new RooFormulaVar( "n_sb",
                                     "mu_ttwj_sb  + mu_qcd_sb  + mu_znn_sb  + eff_sf_sb*( mu_ewo_sb  + mu_susy_sb )",
                                     RooArgSet( *rv_mu_ttwj_sb , *rv_mu_qcd_sb , *rv_mu_znn_sb , *rv_eff_sf_sb, *rv_mu_ewo_sb , *rv_mu_susy_sb  ) ) ;

      rv_n_sig_ldp     = new RooFormulaVar( "n_sig_ldp",
                                     "mu_qcd_sig_ldp + eff_sf_sig_ldp*( sf_mc * (mu_ttwj_sig_ldp + mu_znn_sig_ldp + mu_ewo_sig_ldp) + mu_susy_sig_ldp)",
                                     RooArgSet( *rv_mu_qcd_sig_ldp, *rv_eff_sf_sig_ldp, fv_sf_mc, *rv_mu_ttwj_sig_ldp, *rv_mu_znn_sig_ldp, *rv_mu_ewo_sig_ldp, *rv_mu_susy_sig_ldp ) ) ;

      rv_n_sb_ldp      = new RooFormulaVar( "n_sb_ldp",
                                     "mu_qcd_sb_ldp + eff_sf_sb_ldp*( sf_mc * (mu_ttwj_sb_ldp + mu_znn_sb_ldp + mu_ewo_sb_ldp) + mu_susy_sb_ldp)",
                                     RooArgSet( *rv_mu_qcd_sb_ldp, *rv_eff_sf_sb_ldp, fv_sf_mc, *rv_mu_ttwj_sb_ldp, *rv_mu_znn_sb_ldp, *rv_mu_ewo_sb_ldp, *rv_mu_susy_sb_ldp ) ) ;

      rv_n_sig_sl      = new RooFormulaVar( "n_sig_sl",
                                     "mu_ttwj_sig_sl + eff_sf_sig_sl*mu_susy_sig_sl",
                                     RooArgSet( *rv_mu_ttwj_sig_sl, *rv_eff_sf_sig_sl, *rv_mu_susy_sig_sl ) ) ;

      rv_n_sb_sl       = new RooFormulaVar( "n_sb_sl",
                                     "mu_ttwj_sb_sl + eff_sf_sb_sl*mu_susy_sb_sl",
                                     RooArgSet( *rv_mu_ttwj_sb_sl, *rv_eff_sf_sb_sl, *rv_mu_susy_sb_sl ) ) ;

 ///  rv_n_lsb_0b      = new RooFormulaVar( "n_lsb_0b",
 ///                                 "mu_qcd_lsb_0b",
 ///                                 RooArgSet( *rv_mu_qcd_lsb_0b ) ) ;

 ///  rv_n_lsb_0b_ldp  = new RooFormulaVar( "n_lsb_0b_ldp",
 ///                                 "mu_qcd_lsb_0b_ldp",
 ///                                 RooArgSet( *rv_mu_qcd_lsb_0b_ldp ) ) ;


      rv_n_sig_ee      = new RooFormulaVar( "n_sig_ee",
                                     "mu_zee_sig_ee / fsig_ee",
                                     RooArgSet( *rv_mu_zee_sig_ee, fv_fsig_ee ) ) ;

      rv_n_sb_ee       = new RooFormulaVar( "n_sb_ee",
                                     "mu_zee_sb_ee / fsig_ee",
                                     RooArgSet( *rv_mu_zee_sb_ee, fv_fsig_ee ) ) ;

      rv_n_sig_mm      = new RooFormulaVar( "n_sig_mm",
                                     "mu_zmm_sig_mm / fsig_mm",
                                     RooArgSet( *rv_mu_zmm_sig_mm, fv_fsig_mm ) ) ;

      rv_n_sb_mm       = new RooFormulaVar( "n_sb_mm",
                                     "mu_zmm_sb_mm / fsig_mm",
                                     RooArgSet( *rv_mu_zmm_sb_mm, fv_fsig_mm ) ) ;


   //++++++++++++ PDFs for the likelihood +++++++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining PDFs of the likelihood.\n" ) ;

      pdf_Nsig        = new RooPoisson( "pdf_Nsig"        , "Nsig Poisson PDF"        , *rv_Nsig        , *rv_n_sig ) ;
      pdf_Nsb         = new RooPoisson( "pdf_Nsb"         , "Nsb Poisson PDF"         , *rv_Nsb         , *rv_n_sb ) ;
      pdf_Nsig_ldp    = new RooPoisson( "pdf_Nsig_ldp"    , "Nsig_ldp Poisson PDF"    , *rv_Nsig_ldp    , *rv_n_sig_ldp ) ;
      pdf_Nsb_ldp     = new RooPoisson( "pdf_Nsb_ldp"     , "Nsb_ldp Poisson PDF"     , *rv_Nsb_ldp     , *rv_n_sb_ldp ) ;
      pdf_Nsig_sl     = new RooPoisson( "pdf_Nsig_sl"     , "Nsig_sl Poisson PDF"     , *rv_Nsig_sl     , *rv_n_sig_sl ) ;
      pdf_Nsb_sl      = new RooPoisson( "pdf_Nsb_sl"      , "Nsb_sl Poisson PDF"      , *rv_Nsb_sl      , *rv_n_sb_sl ) ;
 ///  pdf_Nlsb_0b     = new RooPoisson( "pdf_Nlsb_0b"     , "Nlsb_0b Poisson PDF"     , *rv_Nlsb_0b     , *rv_n_lsb_0b ) ;
 ///  pdf_Nlsb_0b_ldp = new RooPoisson( "pdf_Nlsb_0b_ldp" , "Nlsb_0b_ldp Poisson PDF" , *rv_Nlsb_0b_ldp , *rv_n_lsb_0b_ldp ) ;

      pdf_Nsig_ee     = new RooPoisson( "pdf_Nsig_ee"     , "Nsig_ee Poisson PDF"     , *rv_Nsig_ee     , *rv_n_sig_ee ) ;
      pdf_Nsb_ee      = new RooPoisson( "pdf_Nsb_ee"      , "Nsb_ee Poisson PDF"      , *rv_Nsb_ee      , *rv_n_sb_ee ) ;
      pdf_Nsig_mm     = new RooPoisson( "pdf_Nsig_mm"     , "Nsig_mm Poisson PDF"     , *rv_Nsig_mm     , *rv_n_sig_mm ) ;
      pdf_Nsb_mm      = new RooPoisson( "pdf_Nsb_mm"      , "Nsb_mm Poisson PDF"      , *rv_Nsb_mm      , *rv_n_sb_mm ) ;


      {
         RooArgSet pdflist ;
         pdflist.add( *pdf_Nsig        ) ;
         pdflist.add( *pdf_Nsb         ) ;
         pdflist.add( *pdf_Nsig_ldp    ) ;
         pdflist.add( *pdf_Nsb_ldp     ) ;
         pdflist.add( *pdf_Nsig_sl     ) ;
         pdflist.add( *pdf_Nsb_sl      ) ;
   ////  pdflist.add( *pdf_Nlsb_0b     ) ;
   ////  pdflist.add( *pdf_Nlsb_0b_ldp ) ;
         pdflist.add( *pdf_Nsig_ee     ) ;
         pdflist.add( *pdf_Nsb_ee      ) ;
         pdflist.add( *pdf_Nsig_mm     ) ;
         pdflist.add( *pdf_Nsb_mm      ) ;

	 pdflist.add(allNuisancePdfs);

         likelihood = new RooProdPdf( "likelihood", "ra2b likelihood", pdflist ) ;
      }


     //---- Define the list of observables.

       observedParametersList.add( *rv_Nsig        ) ;
       observedParametersList.add( *rv_Nsb         ) ;
       observedParametersList.add( *rv_Nsig_sl     ) ;
       observedParametersList.add( *rv_Nsb_sl      ) ;
       observedParametersList.add( *rv_Nsig_ldp    ) ;
       observedParametersList.add( *rv_Nsb_ldp     ) ;
  //// observedParametersList.add( *rv_Nlsb_0b     ) ;
  //// observedParametersList.add( *rv_Nlsb_0b_ldp ) ;
       observedParametersList.add( *rv_Nsb_ee      ) ;
       observedParametersList.add( *rv_Nsig_ee     ) ;
       observedParametersList.add( *rv_Nsb_mm      ) ;
       observedParametersList.add( *rv_Nsig_mm     ) ;



       dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
                                      observedParametersList ) ;
       dsObserved->add( observedParametersList ) ;


       RooWorkspace workspace ("ws") ;

       workspace.import(*dsObserved);

       // parameters of interest
       RooArgSet poi(*rv_mu_susy_sig, "poi");
       // flat prior for POI
       RooUniform signal_prior ("signal_prior","signal_prior",*rv_mu_susy_sig);

       //RooArgSet * newNuisances = likelihood->getParameters(*dsObserved);
       //RemoveConstantParameters(newNuisances);
       //newNuisances->remove(poi);

       // signal+background model
       ModelConfig sbModel ("SbModel");
       sbModel.SetWorkspace(workspace);
       sbModel.SetPdf(*likelihood);
       sbModel.SetParametersOfInterest(poi);
       //sbModel.SetPriorPdf(signal_prior);
       sbModel.SetNuisanceParameters(allNuisances);
       RooProdPdf nuisancePrior("nuisancePrior","nuisancePrior",allNuisancePdfs);
       sbModel.SetPriorPdf(nuisancePrior);
       //sbModel.SetNuisanceParameters(*newNuisances);
       sbModel.SetObservables(observedParametersList);
       sbModel.SetGlobalObservables(globalObservables);

       // find global maximum with the signal+background model
       // with conditional MLEs for nuisance parameters
       // and save the parameter point snapshot in the Workspace
       //  - safer to keep a default name because some RooStats calculators
       //    will anticipate it
       //RooAbsReal * pNll = sbModel.GetPdf()->createNLL(*dsObserved);
       //RooAbsReal * pProfile = pNll->createProfile(RooArgSet());
       //pProfile->getVal(); // this will do fit and set POI and nuisance parameters to fitted values
       //RooArgSet * pPoiAndNuisance = new RooArgSet();
       //pPoiAndNuisance->add(*sbModel.GetParametersOfInterest());
       //if(sbModel.GetNuisanceParameters()) pPoiAndNuisance->add(*sbModel.GetNuisanceParameters());
       //cout << "\nWill save these parameter points that correspond to the fit to data" << endl;
       //pPoiAndNuisance->Print("v");
       //sbModel.SetSnapshot(*pPoiAndNuisance);
       workspace.import (sbModel);

       //delete pProfile;
       //delete pNll;
       //delete pPoiAndNuisance;
       //delete newNuisances;


       // background-only model
       // use the same PDF as s+b, with xsec=0
       // POI value under the background hypothesis
       ModelConfig bModel (*(RooStats::ModelConfig *)workspace.obj("SbModel"));
       bModel.SetName("BModel");
       bModel.SetWorkspace(workspace);

       // Find a parameter point for generating pseudo-data
       // with the background-only data.
       // Save the parameter point snapshot in the Workspace
       //pNll = bModel.GetPdf()->createNLL(*dsObserved);
       //// bug discovered by Fedor on Sep 6th 2011:
       ////pProfile = pNll->createProfile(poi);
       ////((RooRealVar *)poi.first())->setVal(0.); // set signal = 0
       //pProfile = pNll->createProfile(*bModel.GetParametersOfInterest());
       //((RooRealVar *)(bModel.GetParametersOfInterest()->first()))->setVal(0.); // set signal = 0
       //pProfile->getVal(); // this will do fit and set nuisance parameters to profiled values
       //pPoiAndNuisance = new RooArgSet();
       //pPoiAndNuisance->add(*bModel.GetParametersOfInterest());
       //if(bModel.GetNuisanceParameters()) pPoiAndNuisance->add(*bModel.GetNuisanceParameters());
       //cout << "\nShould use these parameter points to generate pseudo data for bkg only" << endl;
       //pPoiAndNuisance->Print("v");
       ((RooRealVar *)(bModel.GetParametersOfInterest()->first()))->setVal(0.);
       bModel.SetSnapshot(*bModel.GetParametersOfInterest());
       workspace.import (bModel);

       //delete pProfile;
       //delete pNll;
       //delete pPoiAndNuisance;

       workspace.Print() ;
       workspace.writeToFile("ws.root");

       return true ;


    } // initialize.


   //===================================================================================================================================


    bool ra2bRoostatsClass8::setSusyScanPoint( const char* inputScanFile, double m0, double m12, bool isT1bbbb, double t1bbbbXsec ) {


       //--- Aug 15, 2011: updated to new format for AN, v3.


       printf("\n\n Opening SUSY scan input file : %s\n", inputScanFile ) ;

       FILE* infp ;
       if ( (infp=fopen( inputScanFile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", inputScanFile ) ;
          return false ;
       }

       double deltaM0(0.) ;
       double deltaM12(0.) ;

       if ( !isT1bbbb ) {
          deltaM0 = 20 ;
          deltaM12 = 20 ;
       } else {
          deltaM0 = 25 ;
          deltaM12 = 25 ;
       }

       bool found(false) ;

       //--- Loop over the scan points.
       while ( !feof( infp ) ) {

          float pointM0 ;
          float pointM12 ;



   //+++ ORL: Aug 14, 2011 ++++++++++++++++++++++++++
          float n_sig_raw ;
          float n_sb_raw ;
          float n_sig_sl_raw ;
          float n_sb_sl_raw ;
          float n_sig_ldp_raw ;
          float n_sb_ldp_raw ;

          float n_sig_correction ;
          float n_sb_correction ;
          float n_sig_sl_correction ;
          float n_sb_sl_correction ;
          float n_sig_ldp_correction ;
          float n_sb_ldp_correction ;

          float n_sig_error ;
          float n_sb_error ;
          float n_sig_sl_error ;
          float n_sb_sl_error ;
          float n_sig_ldp_error ;
          float n_sb_ldp_error ;

          int nGen ;

          fscanf( infp, "%f %f %d  %f %f %f %f %f %f  %f %f %f %f %f %f %f %f %f %f %f %f",
            &pointM0, &pointM12, &nGen,
            &n_sig_raw, &n_sb_raw, &n_sig_sl_raw, &n_sb_sl_raw, &n_sig_ldp_raw, &n_sb_ldp_raw,
            &n_sig_correction, &n_sb_correction, &n_sig_sl_correction, &n_sb_sl_correction, &n_sig_ldp_correction, &n_sb_ldp_correction,
            &n_sig_error, &n_sb_error, &n_sig_sl_error, &n_sb_sl_error, &n_sig_ldp_error, &n_sb_ldp_error ) ;

          if ( feof(infp) ) break ;
          if ( n_sig_raw < 0.00001 ) continue ;
   //--If you are asking for it, I'll assume it's good.  Josh is using 0 for ngen dummy in LM9.
   /////  if ( nGen != 10000 ) continue ; // get rid of bad scan points.

          if (    fabs( pointM0 - m0 ) <= deltaM0/2.
               && fabs( pointM12 - m12 ) <= deltaM12/2. ) {

              double nGenPerPoint = 10000 ; // for t1bbbb

             printf("\n\n") ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_error     = %6.1f %%\n", n_sig_error     ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_error      = %6.1f %%\n", n_sb_error      ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_sl_error  = %6.1f %%\n", n_sig_sl_error  ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_sl_error   = %6.1f %%\n", n_sb_sl_error   ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_ldp_error = %6.1f %%\n", n_sig_ldp_error ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_ldp_error  = %6.1f %%\n", n_sb_ldp_error  ) ;
             printf("\n\n") ;

      //--Include the stat error on the efficiency for t1bbbb.
             if ( isT1bbbb ) {

              //-- absolute raw eff
                 float n_sig_raw_eff     = n_sig_raw     / nGenPerPoint ;
                 float n_sb_raw_eff      = n_sb_raw      / nGenPerPoint ;
                 float n_sig_sl_raw_eff  = n_sig_sl_raw  / nGenPerPoint ;
                 float n_sb_sl_raw_eff   = n_sb_sl_raw   / nGenPerPoint ;
                 float n_sig_ldp_raw_eff = n_sig_ldp_raw / nGenPerPoint ;
                 float n_sb_ldp_raw_eff  = n_sb_ldp_raw  / nGenPerPoint ;


               //-- absolute stat err.
                 float n_sig_stat_error     =  sqrt(  n_sig_raw_eff    * ( 1.0 - n_sig_raw_eff     ) / nGenPerPoint ) ;
                 float n_sb_stat_error      =  sqrt(  n_sb_raw_eff     * ( 1.0 - n_sb_raw_eff      ) / nGenPerPoint ) ;
                 float n_sig_sl_stat_error  =  sqrt(  n_sig_sl_raw_eff * ( 1.0 - n_sig_sl_raw_eff  ) / nGenPerPoint ) ;
                 float n_sb_sl_stat_error   =  sqrt(  n_sb_sl_raw_eff  * ( 1.0 - n_sb_sl_raw_eff   ) / nGenPerPoint ) ;
                 float n_sig_ldp_stat_error =  sqrt(  n_sig_ldp_raw_eff* ( 1.0 - n_sig_ldp_raw_eff ) / nGenPerPoint ) ;
                 float n_sb_ldp_stat_error  =  sqrt(  n_sb_ldp_raw_eff * ( 1.0 - n_sb_ldp_raw_eff  ) / nGenPerPoint ) ;

               //-- relative stat err in percent.
                 if ( n_sig_raw_eff     > 0 ) { n_sig_stat_error     = 100.* n_sig_stat_error     / n_sig_raw_eff     ; } else { n_sig_stat_error     = 0. ; }
                 if ( n_sb_raw_eff      > 0 ) { n_sb_stat_error      = 100.* n_sb_stat_error      / n_sb_raw_eff      ; } else { n_sb_stat_error      = 0. ; }
                 if ( n_sig_sl_raw_eff  > 0 ) { n_sig_sl_stat_error  = 100.* n_sig_sl_stat_error  / n_sig_sl_raw_eff  ; } else { n_sig_sl_stat_error  = 0. ; }
                 if ( n_sb_sl_raw_eff   > 0 ) { n_sb_sl_stat_error   = 100.* n_sb_sl_stat_error   / n_sb_sl_raw_eff   ; } else { n_sb_sl_stat_error   = 0. ; }
                 if ( n_sig_ldp_raw_eff > 0 ) { n_sig_ldp_stat_error = 100.* n_sig_ldp_stat_error / n_sig_ldp_raw_eff ; } else { n_sig_ldp_stat_error = 0. ; }
                 if ( n_sb_ldp_raw_eff  > 0 ) { n_sb_ldp_stat_error  = 100.* n_sb_ldp_stat_error  / n_sb_ldp_raw_eff  ; } else { n_sb_ldp_stat_error  = 0. ; }

                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_stat_error     = %6.1f %%\n", n_sig_stat_error     ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_stat_error      = %6.1f %%\n", n_sb_stat_error      ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_sl_stat_error  = %6.1f %%\n", n_sig_sl_stat_error  ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_sl_stat_error   = %6.1f %%\n", n_sb_sl_stat_error   ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_ldp_stat_error = %6.1f %%\n", n_sig_ldp_stat_error ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_ldp_stat_error  = %6.1f %%\n", n_sb_ldp_stat_error  ) ;

               //-- total err in percent.
                 n_sig_error     = sqrt( pow( n_sig_error    , 2) + pow( n_sig_stat_error    , 2) ) ;
                 n_sb_error      = sqrt( pow( n_sb_error     , 2) + pow( n_sb_stat_error     , 2) ) ;
                 n_sig_sl_error  = sqrt( pow( n_sig_sl_error , 2) + pow( n_sig_sl_stat_error , 2) ) ;
                 n_sb_sl_error   = sqrt( pow( n_sb_sl_error  , 2) + pow( n_sb_sl_stat_error  , 2) ) ;
                 n_sig_ldp_error = sqrt( pow( n_sig_ldp_error, 2) + pow( n_sig_ldp_stat_error, 2) ) ;
                 n_sb_ldp_error  = sqrt( pow( n_sb_ldp_error , 2) + pow( n_sb_ldp_stat_error , 2) ) ;

                 printf("\n\n") ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_error     = %6.1f %%\n", n_sig_error     ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_error      = %6.1f %%\n", n_sb_error      ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_sl_error  = %6.1f %%\n", n_sig_sl_error  ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_sl_error   = %6.1f %%\n", n_sb_sl_error   ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_ldp_error = %6.1f %%\n", n_sig_ldp_error ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_ldp_error  = %6.1f %%\n", n_sb_ldp_error  ) ;
                 printf("\n\n") ;

             }

       //--- Not needed with log-normal
        ///  //-- enforce a maximum efficiency uncertainty (to avoid negative scale factors).
        ///  if ( n_sig_error     > 35. ) { n_sig_error     = 35. ; }
        ///  if ( n_sb_error      > 35. ) { n_sb_error      = 35. ; }
        ///  if ( n_sig_sl_error  > 35. ) { n_sig_sl_error  = 35. ; }
        ///  if ( n_sb_sl_error   > 35. ) { n_sb_sl_error   = 35. ; }
        ///  if ( n_sig_ldp_error > 35. ) { n_sig_ldp_error = 35. ; }
        ///  if ( n_sb_ldp_error  > 35. ) { n_sb_ldp_error  = 35. ; }

      //++++++++++++++++++++++++++++++++++++++++++++++++



             double setVal_n_sig(0.) ;
             double setVal_n_sb(0.) ;
             double setVal_n_sig_sl(0.) ;
             double setVal_n_sb_sl(0.) ;
             double setVal_n_sig_ldp(0.) ;
             double setVal_n_sb_ldp(0.) ;


             if ( !isT1bbbb ) {
                //-- tanb40
                setVal_n_sig     = n_sig_raw     * n_sig_correction     ;
                setVal_n_sb      = n_sb_raw      * n_sb_correction      ;
                setVal_n_sig_sl  = n_sig_sl_raw  * n_sig_sl_correction  ;
                setVal_n_sb_sl   = n_sb_sl_raw   * n_sb_sl_correction   ;
                setVal_n_sig_ldp = n_sig_ldp_raw * n_sig_ldp_correction ;
                setVal_n_sb_ldp  = n_sb_ldp_raw  * n_sb_ldp_correction  ;
             } else {
                //-- t1bbbb
                setVal_n_sig     = DataLumi * t1bbbbXsec * (( n_sig_raw     * n_sig_correction     )/ nGenPerPoint ) ;
                setVal_n_sb      = DataLumi * t1bbbbXsec * (( n_sb_raw      * n_sb_correction      )/ nGenPerPoint ) ;
                setVal_n_sig_sl  = DataLumi * t1bbbbXsec * (( n_sig_sl_raw  * n_sig_sl_correction  )/ nGenPerPoint ) ;
                setVal_n_sb_sl   = DataLumi * t1bbbbXsec * (( n_sb_sl_raw   * n_sb_sl_correction   )/ nGenPerPoint ) ;
                setVal_n_sig_ldp = DataLumi * t1bbbbXsec * (( n_sig_ldp_raw * n_sig_ldp_correction )/ nGenPerPoint ) ;
                setVal_n_sb_ldp  = DataLumi * t1bbbbXsec * (( n_sb_ldp_raw  * n_sb_ldp_correction  )/ nGenPerPoint ) ;
             }

             rv_mu_susymc_sig       -> setVal( setVal_n_sig      ) ;
             rv_mu_susymc_sb        -> setVal( setVal_n_sb       ) ;
             rv_mu_susymc_sig_sl    -> setVal( setVal_n_sig_sl   ) ;
             rv_mu_susymc_sb_sl     -> setVal( setVal_n_sb_sl    ) ;
             rv_mu_susymc_sig_ldp   -> setVal( setVal_n_sig_ldp  ) ;
             rv_mu_susymc_sb_ldp    -> setVal( setVal_n_sb_ldp   ) ;

             rv_width_eff_sf_sig     -> setVal( n_sig_error     / 100. ) ;
             rv_width_eff_sf_sb      -> setVal( n_sb_error      / 100. ) ;
             rv_width_eff_sf_sig_sl  -> setVal( n_sig_sl_error  / 100. ) ;
             rv_width_eff_sf_sb_sl   -> setVal( n_sb_sl_error   / 100. ) ;
             rv_width_eff_sf_sig_ldp -> setVal( n_sig_ldp_error / 100. ) ;
             rv_width_eff_sf_sb_ldp  -> setVal( n_sb_ldp_error  / 100. ) ;



             if ( !isT1bbbb ) {
                printf("\n\n Found point m0 = %4.0f,  m1/2 = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_sig ) ;
             } else {
                printf("\n\n Found point mGluino = %4.0f,  mLSP = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_sig ) ;
             }


             printf("\n\n") ;
             printf(" Setting susy N_sig     to  %7.1f\n", setVal_n_sig       ) ;
             printf(" Setting susy N_sb      to  %7.1f\n", setVal_n_sb        ) ;
             printf(" Setting susy N_sig_sl  to  %7.1f\n", setVal_n_sig_sl    ) ;
             printf(" Setting susy N_sb_sl   to  %7.1f\n", setVal_n_sb_sl     ) ;
             printf(" Setting susy N_sig_ldp to  %7.1f\n", setVal_n_sig_ldp   ) ;
             printf(" Setting susy N_sb_ldp  to  %7.1f\n", setVal_n_sb_ldp    ) ;
             printf("\n\n") ;

             found = true ;

             break ;

          } // point match?

       } // not eof ?

       fclose( infp ) ;

       if ( found ) {
          return true ;
       } else {
          printf("\n\n *** Point not found in scan.\n\n" ) ;
          return false ;
       }

    } // setSusyScanPoint.
