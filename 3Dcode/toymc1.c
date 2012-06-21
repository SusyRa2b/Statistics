#include "ra2bRoostatsClass3D_1.c"
#include "TRoot.h"
#include "TSystem.h"
#include "TRandom.h"

   static const int nBinsMET  = 4 ;
   static const int nBinsHT   = 4 ;

   static const int nBinsBjets = 3 ;    // this must always be 3


   //-- Global variables.

   float toy_mean_N_0lep[nBinsMET][nBinsHT][nBinsBjets] ;
   float toy_mean_N_1lep[nBinsMET][nBinsHT][nBinsBjets] ;
   float toy_mean_N_ldp [nBinsMET][nBinsHT][nBinsBjets] ;
   float toy_mean_N_Zee [nBinsMET][nBinsHT] ;
   float toy_mean_N_Zmm [nBinsMET][nBinsHT] ;

   float susy_0lep [nBinsMET][nBinsHT][nBinsBjets] ;
   float susy_1lep [nBinsMET][nBinsHT][nBinsBjets] ;
   float susy_ldp  [nBinsMET][nBinsHT][nBinsBjets] ;

   char datfile[10000] ;
   char susyfile[10000] ;
   char outputDir[10000] ;

   int   nFloatParInitVal ;
   char  floatParName[5000][100] ;
   float floatParInitVal[5000] ;

   float mgl, mlsp ;

   float nSusy0lep ;

   float susyPoi0lepRatio ;

   RooWorkspace* workspace ;

   RooRealVar* rrv_susy_poi ;
   RooAbsPdf* likelihood ;

   float true_susy_0lep ;
   float fit_susy_0lep ;
   float fit_susy_0lep_err ;
   float fit_susy_0lep_err_low ;
   float fit_susy_0lep_err_high ;
   float fit_susy_0lep_err_forpull ;
   float true_ttwj_0lep ;
   float fit_ttwj_0lep ;
   float true_qcd_0lep ;
   float fit_qcd_0lep ;
   float true_znn_0lep ;
   float fit_znn_0lep ;
   int   fit_covqual_susyfloat ;

   float fit_susy_signif ;

   float fit_susy_ul ;
   float fit_susy_ts_at_ul ;

   double minNllSusyFloat ;
   double susy_poi_atMinNll ;
   double susy_poi_plusErr ;


   bool doSignif ;
   bool doUL ;

   bool useExpected0lep ;




   //-- Prototypes

   bool readMCvals() ;
   bool addSusyToExpectedObs() ;
   RooDataSet* genToyData() ;
   bool processToyResult( int ti, RooFitResult* fitResult ) ;
   bool bookTree() ;
   float findUL( RooDataSet* toyds ) ;
   bool setReinitValues( RooFitResult* fitResult ) ;
   bool reinitFloatPars() ;


   TRandom* tran ;

   TTree* toytt ;
   TFile* ttfile ;

   //=====================================

   void toymc1( const char* input_datfile = "Input-met4-ht4.dat",
//                const char* input_susyfile = "Susy-mgl900-mlsp300-met4-ht4-v2.dat",
                const char* input_susyfile = "T1bbbb-met4-ht4-v2.dat",
                double input_mgl=900, double input_mlsp=300.,
                const char* input_deffdbtagfile = "dummy_DeffDbtag-met4-ht4-v2.dat",
                double input_nSusy0lep = 60.,
                const char* input_outputDir = "output-toymc1",
		int nToy = 10
                        ) {

       char command[10000] ;


       sprintf( datfile, "%s", input_datfile ) ;
       sprintf( susyfile, "%s", input_susyfile ) ;
       sprintf( outputDir, "%s", input_outputDir ) ;

       mgl  = input_mgl ;
       mlsp = input_mlsp ;

       nSusy0lep = input_nSusy0lep ;

       susyPoi0lepRatio = 0. ;
       true_susy_0lep = 0. ;

       doSignif = true ;
       doUL = false ;

       useExpected0lep = true ;


       tran = new TRandom(12345) ;



       //--- Create workspace.

       ra2bRoostatsClass3D_1 ra2b ;

       int qcdModelIndex = 1 ;
       ra2b.initialize( input_datfile, input_susyfile, mgl, mlsp, false, 0., input_deffdbtagfile, qcdModelIndex ) ;

       sprintf( command, "mkdir -p %s", outputDir ) ;
       gSystem->Exec( command ) ;
       sprintf( command, "mv ws.root %s/", outputDir ) ;
       gSystem->Exec( command ) ;

       char wsfilename[10000] ;
       sprintf( wsfilename, "%s/ws.root", outputDir ) ;
       TFile wsfile( wsfilename ) ;
       workspace = (RooWorkspace*) wsfile.Get("ws") ;
       if ( workspace == 0x0 ) {
          printf("\n\n *** Can't find the workspace.\n\n" ) ; return ;
       }






       //--- Access some workspace stuff needed later.
       rrv_susy_poi = workspace -> var( "mu_susy_M1_H1_1b" ) ;
       if ( rrv_susy_poi == 0x0 ) {
          printf("\n\n *** can't find susy poi mu_susy_M1_H1_1b.\n\n") ; return ;
       }

       ModelConfig* modelConfig = (ModelConfig*) workspace -> obj( "SbModel" ) ;
       if ( modelConfig == 0x0 ) { printf("\n\n *** can't find ModelConfig with name SbModel.\n\n") ; return ; }
       likelihood = modelConfig->GetPdf() ;








       //--- Do an initial fit to determine reasonable starting values for all floating parameters.

       RooDataSet* rdsMCvals = (RooDataSet*) workspace->obj( "ra2b_observed_rds" ) ;
       if ( rdsMCvals == 0x0 ) {
          printf("\n\n *** can't find dataset with name ra2b_observed_rds in workspace.\n\n") ;
          return ;
       }

       rrv_susy_poi -> setVal(0.) ;
       rrv_susy_poi -> setConstant( kTRUE ) ;

       RooFitResult* mcfitResult = likelihood->fitTo( *rdsMCvals, Save(true), PrintLevel(0) ) ;

       if ( ! setReinitValues( mcfitResult ) ) {
          printf("\n\n *** Problem in collecting inital vals for floating parameters from MC fit.\n\n") ; return ;
       }








       //--- Read in MC predictions from Input*.dat file and set the expected values for all observables.

       if ( ! readMCvals() ) {
          printf("\n\n *** Problem reading MC values.\n\n") ; return ;
       }







       //--- Add some SUSY, if requested.
       if ( ! addSusyToExpectedObs() ) {
          printf("\n\n *** Problem adding in some SUSY to expected observables.\n\n") ; return ;
       }













       //--- Create tree and define tree variables.
       bookTree() ;




       //--- Loop over toy experiments.

       for ( int ti=0; ti<nToy; ti++ ) {

          printf("\n\n\n\n ====== Begin toy experiment %d\n\n\n out of %d", ti, nToy ) ;




          //-- Do a fit with susy floating.  Get MINOS errors on susy yield.

          RooDataSet* toyds = genToyData() ;
          if ( toyds == 0x0 ) { printf("\n\n *** Fail!\n\n") ; return ; }
          toyds->printMultiline( cout, 1, kTRUE, "" ) ;

          RooArgSet ras ;
          ras.add( *rrv_susy_poi ) ;
          rrv_susy_poi->setConstant( kFALSE ) ;

          if ( ! reinitFloatPars() ) {
             printf("\n\n *** problem reinitializing floating pars.\n\n") ; return ;
          }

          RooFitResult* fitResult = likelihood -> fitTo( *toyds, Save(true), Hesse(true), Minos(ras), Strategy(1), PrintLevel(0) ) ;
          minNllSusyFloat = fitResult->minNll() ;

          processToyResult( ti, fitResult ) ;




          //-- Calculate significance from delta log likelihood, if requested.

          if ( doSignif ) {

             rrv_susy_poi->setVal(0.) ;
             rrv_susy_poi->setConstant( kTRUE ) ;
             RooFitResult* fitResult0susy = likelihood -> fitTo( *toyds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
             double minNll0Susy = fitResult0susy->minNll() ;

             double testStat = 2.*( minNll0Susy - minNllSusyFloat ) ;
             fit_susy_signif = 0. ;
             if ( testStat > 0 ) { fit_susy_signif = sqrt( testStat ) ; }
             printf("  test stat = %5.2f,  significance = %5.2f\n", testStat, fit_susy_signif ) ;

             delete fitResult0susy ;

             rrv_susy_poi->setConstant( kFALSE ) ;

          } // doSignif?





          //-- Calculate one-sided upper limit on susy signal, if requested.

          if ( doUL ) {

             findUL( toyds ) ;

          }




          //-- Clean up and save results for this toy in the TTree.

          delete fitResult ;
          delete toyds ;

          toytt -> Fill() ;

       } // ti.







       //--- Save output.

       ttfile->cd() ;
       toytt->Write() ;
       ttfile->Close() ;






   }


   //==================================================================================


   bool readMCvals() {

      //==== Data arrays.

      float N_0lep_input [nBinsMET][nBinsHT][nBinsBjets] ;
      float N_1lep       [nBinsMET][nBinsHT][nBinsBjets] ;
      float N_ldp        [nBinsMET][nBinsHT][nBinsBjets] ;
      float N_mc_ldp     [nBinsMET][nBinsHT][nBinsBjets] ;
      float N_Zee        [nBinsMET][nBinsHT] ;
      float N_Zmm        [nBinsMET][nBinsHT] ;

      float acc_Zee      [nBinsMET] ;
      float acc_Zmm      [nBinsMET] ;

      float Z_ee_eff(0.) ;
      float Z_mm_eff(0.) ;

      float knn[nBinsBjets] ;

      float Z_ee_pur(0.) ;
      float Z_mm_pur(0.) ;

      float R_passfail[nBinsHT] ;


      int nread ;
      char pname[1000] ;
      float pval ;
      char c ;

      FILE* infp ;



      //--- Go through the file again and read in data needed to compute the expected values
      //    for the 0lep observables.

      if ( (infp=fopen( datfile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", datfile ) ;
         return false ;
      }


      //--- read in description line.
      c = 0 ;
      while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }

      //--- Read in N_0lep lines
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               nread = fscanf( infp, "%s %g", pname, &pval ) ;
               ///// printf("  %d : %s %g\n", nread, pname, pval ) ;
               int metbin, htbin, nbbin ;
               int nmatch = sscanf( pname, "N_0lep_M%d_H%d_%db", &metbin, &htbin, &nbbin ) ;
               if ( nmatch == 3 ) {
                  if ( metbin == mbi+1 && htbin == hbi+1 && nbbin == bbi+1 ) {
                     N_0lep_input[mbi][hbi][bbi] = pval ;
                     printf(" Set N_0lep_input for mbi=%d, hbi=%d, bbi=%d to %g\n", mbi, hbi, bbi, pval ) ;
                  } else {
                     printf("\n\n *** input mismatch. expect %d, %d, %d.  found %d, %d, %d\n\n",
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return false ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return false ;
               }
            } // bbi.
         } // hbi
      } // mbi.

      //--- Read in N_1lep lines
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               nread = fscanf( infp, "%s %g", pname, &pval ) ;
               ///// printf("  %d : %s %g\n", nread, pname, pval ) ;
               int metbin, htbin, nbbin ;
               int nmatch = sscanf( pname, "N_1lep_M%d_H%d_%db", &metbin, &htbin, &nbbin ) ;
               if ( nmatch == 3 ) {
                  if ( metbin == mbi+1 && htbin == hbi+1 && nbbin == bbi+1 ) {
                     N_1lep[mbi][hbi][bbi] = pval ;
                     printf(" Set N_1lep for mbi=%d, hbi=%d, bbi=%d to %g\n", mbi, hbi, bbi, pval ) ;
                  } else {
                     printf("\n\n *** input mismatch. expect %d, %d, %d.  found %d, %d, %d\n\n",
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return false ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return false ;
               }
            } // bbi.
         } // hbi
      } // mbi.

      //--- Read in N_ldp lines
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               nread = fscanf( infp, "%s %g", pname, &pval ) ;
               ///// printf("  %d : %s %g\n", nread, pname, pval ) ;
               int metbin, htbin, nbbin ;
               int nmatch = sscanf( pname, "N_ldp_M%d_H%d_%db", &metbin, &htbin, &nbbin ) ;
               if ( nmatch == 3 ) {
                  if ( metbin == mbi+1 && htbin == hbi+1 && nbbin == bbi+1 ) {
                     N_ldp[mbi][hbi][bbi] = pval ;
                     printf(" Set N_ldp for mbi=%d, hbi=%d, bbi=%d to %g\n", mbi, hbi, bbi, pval ) ;
                  } else {
                     printf("\n\n *** input mismatch. expect %d, %d, %d.  found %d, %d, %d\n\n",
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return false ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return false ;
               }
            } // bbi.
         } // hbi
      } // mbi.

      //--- Read R_lsb lines
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
            nread = fscanf( infp, "%s %g", pname, &pval ) ;
            nread = fscanf( infp, "%s %g", pname, &pval ) ;
         } // bbi.
      } // hbi

      //--- Read in Zee lines
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            nread = fscanf( infp, "%s %g", pname, &pval ) ;
            ///// printf("  %d : %s %g\n", nread, pname, pval ) ;
            int metbin, htbin ;
            int nmatch = sscanf( pname, "N_Zee_M%d_H%d", &metbin, &htbin ) ;
            if ( nmatch == 2 ) {
               if ( metbin == mbi+1 && htbin == hbi+1 ) {
                  N_Zee[mbi][hbi] = pval ;
                  printf(" Set N_Zee for mbi=%d, hbi=%d to %g\n", mbi, hbi, pval ) ;
               }
            } else {
               printf("\n\n *** bad input line.\n\n") ; return false ;
            }
         } // hbi.
      } // mbi.

      //--- Read in Zmm lines
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            nread = fscanf( infp, "%s %g", pname, &pval ) ;
            ///// printf("  %d : %s %g\n", nread, pname, pval ) ;
            int metbin, htbin ;
            int nmatch = sscanf( pname, "N_Zmm_M%d_H%d", &metbin, &htbin ) ;
            if ( nmatch == 2 ) {
               if ( metbin == mbi+1 && htbin == hbi+1 ) {
                  N_Zmm[mbi][hbi] = pval ;
                  printf(" Set N_Zmm for mbi=%d, hbi=%d to %g\n", mbi, hbi, pval ) ;
               }
            } else {
               printf("\n\n *** bad input line.\n\n") ; return false ;
            }
         } // hbi.
      } // mbi.

      //--- Read in N_ttbarsingletopzjetsmc_ldp lines
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               nread = fscanf( infp, "%s %g", pname, &pval ) ;
               ///// printf("  %d : %s %g\n", nread, pname, pval ) ;
               int metbin, htbin, nbbin ;
               int nmatch = sscanf( pname, "N_ttbarsingletopzjetsmc_ldp_M%d_H%d_%db", &metbin, &htbin, &nbbin ) ;
               if ( nmatch == 3 ) {
                  if ( metbin == mbi+1 && htbin == hbi+1 && nbbin == bbi+1 ) {
                     N_mc_ldp[mbi][hbi][bbi] = pval ;
                     printf(" Set N_ttbarsingletopzjetsmc_ldp for mbi=%d, hbi=%d, bbi=%d to %g\n", mbi, hbi, bbi, pval ) ;
                  } else {
                     printf("\n\n *** input mismatch. expect %d, %d, %d.  found %d, %d, %d\n\n",
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return false ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return false ;
               }
            } // bbi.
         } // hbi
      } // mbi.

      //--- Read in N_WJmc_ldp lines
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               nread = fscanf( infp, "%s %g", pname, &pval ) ;
               ///// printf("  %d : %s %g\n", nread, pname, pval ) ;
               int metbin, htbin, nbbin ;
               int nmatch = sscanf( pname, "N_WJmc_ldp_M%d_H%d_%db", &metbin, &htbin, &nbbin ) ;
               if ( nmatch == 3 ) {
                  if ( metbin == mbi+1 && htbin == hbi+1 && nbbin == bbi+1 ) {
                     N_mc_ldp[mbi][hbi][bbi] += pval ;
                     printf(" Set N_WJmc_ldp for mbi=%d, hbi=%d, bbi=%d to %g\n", mbi, hbi, bbi, pval ) ;
                  } else {
                     printf("\n\n *** input mismatch. expect %d, %d, %d.  found %d, %d, %d\n\n",
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return false ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return false ;
               }
            } // bbi.
         } // hbi
      } // mbi.

      //--- Read in N_Znnmc_ldp lines
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               nread = fscanf( infp, "%s %g", pname, &pval ) ;
               ///// printf("  %d : %s %g\n", nread, pname, pval ) ;
               int metbin, htbin, nbbin ;
               int nmatch = sscanf( pname, "N_Znnmc_ldp_M%d_H%d_%db", &metbin, &htbin, &nbbin ) ;
               if ( nmatch == 3 ) {
                  if ( metbin == mbi+1 && htbin == hbi+1 && nbbin == bbi+1 ) {
                     N_mc_ldp[mbi][hbi][bbi] += pval ;
                     printf(" Set N_Znnmc_ldp for mbi=%d, hbi=%d, bbi=%d to %g\n", mbi, hbi, bbi, pval ) ;
                  } else {
                     printf("\n\n *** input mismatch. expect %d, %d, %d.  found %d, %d, %d\n\n",
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return false ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return false ;
               }
            } // bbi.
         } // hbi
      } // mbi.


      //--- Read acc_Zee lines
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
         int metbin ;
         int nmatch = sscanf( pname, "acc_Zee_M%d", &metbin ) ;
         if ( nmatch == 1 ) {
            acc_Zee[mbi] = pval ;
            printf(" Set acc_Zee for mbi=%d to %g\n", mbi, pval ) ;
         } else {
            printf("\n\n *** bad input line.\n\n") ; return false ;
         }
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
      } // mbi


      //--- Read acc_Zmm lines
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
         int metbin ;
         int nmatch = sscanf( pname, "acc_Zmm_M%d", &metbin ) ;
         if ( nmatch == 1 ) {
            acc_Zmm[mbi] = pval ;
            printf(" Set acc_Zmm for mbi=%d to %g\n", mbi, pval ) ;
         } else {
            printf("\n\n *** bad input line.\n\n") ; return false ;
         }
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
      } // mbi


      //--- Z_ee_eff
      nread = fscanf( infp, "%s %g", pname, &pval ) ;
      if ( strcmp( pname, "Z_ee_eff" ) == 0 ) {
         Z_ee_eff = pval ;
         printf(" Set Z_ee_eff to %g\n", pval ) ;
      } else {
         printf("\n\n *** bad input line.\n\n") ; return false ;
      }
      nread = fscanf( infp, "%s %g", pname, &pval ) ;


      //--- Z_mm_eff
      nread = fscanf( infp, "%s %g", pname, &pval ) ;
      if ( strcmp( pname, "Z_mm_eff" ) == 0 ) {
         Z_mm_eff = pval ;
         printf(" Set Z_mm_eff to %g\n", pval ) ;
      } else {
         printf("\n\n *** bad input line.\n\n") ; return false ;
      }
      nread = fscanf( infp, "%s %g", pname, &pval ) ;


      //--- knn
      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
         int nbbin ;
         int nmatch = sscanf( pname, "knn_%db", &nbbin ) ;
         if ( nmatch == 1 ) {
            if ( nbbin == bbi+1 ) {
               knn[bbi] = pval ;
               printf(" Set knn for bbi=%d to %g\n", bbi, pval ) ;
            } else {
               printf("\n\n *** bad input line.\n\n") ; return false ;
            }
         } else {
            printf("\n\n *** bad input line.\n\n") ; return false ;
         }
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
      } // bbi.


      //--- Z_ee_pur
      nread = fscanf( infp, "%s %g", pname, &pval ) ;
      if ( strcmp( pname, "Z_ee_pur" ) == 0 ) {
         Z_ee_pur = pval ;
         printf(" Set Z_ee_pur to %g\n", pval ) ;
      } else {
         printf("\n\n *** bad input line.\n\n") ; return false ;
      }
      nread = fscanf( infp, "%s %g", pname, &pval ) ;


      //--- Z_mm_pur
      nread = fscanf( infp, "%s %g", pname, &pval ) ;
      if ( strcmp( pname, "Z_mm_pur" ) == 0 ) {
         Z_mm_pur = pval ;
         printf(" Set Z_mm_pur to %g\n", pval ) ;
      } else {
         printf("\n\n *** bad input line.\n\n") ; return false ;
      }
      nread = fscanf( infp, "%s %g", pname, &pval ) ;

      fclose(infp) ;







      //--- Compute the expected values for the 0lep observables, not including susy.
      float Zee_factor[nBinsBjets] ;
      float Zmm_factor[nBinsBjets] ;

      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
         Zee_factor[bbi] = ( 5.94 * Z_ee_pur * knn[bbi] ) / ( acc_Zee[0] * Z_ee_eff ) ; // ignoring MET dependence of acceptance.
         Zmm_factor[bbi] = ( 5.94 * Z_mm_pur * knn[bbi] ) / ( acc_Zmm[0] * Z_mm_eff ) ; // ignoring MET dependence of acceptance.
         printf("  nB=%d Zll scale factors:  ee=%g, mm=%g\n", bbi, Zee_factor[bbi], Zmm_factor[bbi] ) ;
      } // bbi.

      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         R_passfail[hbi] = 0.07 ; // guess for now.
      }

      float R_ttwj_0lep_over_1lep = 1.15 ; // guess for now.

      true_ttwj_0lep = 0. ;
      true_qcd_0lep  = 0. ;
      true_znn_0lep  = 0. ;



      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               float exp_0lep_ttwj = N_1lep[mbi][hbi][bbi] * R_ttwj_0lep_over_1lep ;

               float exp_0lep_qcd =  ( N_ldp[mbi][hbi][bbi] - N_mc_ldp[mbi][hbi][bbi] ) * R_passfail[hbi] ;
               if ( exp_0lep_qcd < 0. ) { exp_0lep_qcd = 0. ; }

               float exp_0lep_znn = 0.5 * ( N_Zee[mbi][hbi] * Zee_factor[bbi] + N_Zmm[mbi][hbi] * Zmm_factor[bbi] ) ;

               float exp_0lep = exp_0lep_ttwj + exp_0lep_qcd + exp_0lep_znn ;

               if ( useExpected0lep ) {
                  toy_mean_N_0lep[mbi][hbi][bbi] = exp_0lep ;
               } else {
                  toy_mean_N_0lep[mbi][hbi][bbi] = N_0lep_input[mbi][hbi][bbi] ;
               }
               toy_mean_N_1lep[mbi][hbi][bbi] = N_1lep[mbi][hbi][bbi] ;
               toy_mean_N_ldp [mbi][hbi][bbi] = N_ldp [mbi][hbi][bbi] ;

               if ( useExpected0lep ) { // these are only well defined if using expected values for 0lep.
                  true_ttwj_0lep += exp_0lep_ttwj ;
                  true_qcd_0lep  += exp_0lep_qcd  ;
                  true_znn_0lep  += exp_0lep_znn  ;
               }

               printf(" 0lep: m,h,b (%d,%d,%d): ttwj=%5.1f, qcd=%5.1f, Znn=%5.1f,  expected total=%5.1f\n",
                    mbi, hbi, bbi, exp_0lep_ttwj, exp_0lep_qcd, exp_0lep_znn, toy_mean_N_0lep[mbi][hbi][bbi] ) ;

            } // bbi.
            toy_mean_N_Zee[mbi][hbi] = N_Zee[mbi][hbi] ;
            toy_mean_N_Zmm[mbi][hbi] = N_Zmm[mbi][hbi] ;
         } // hbi.
      } // mbi.


      return true ;


   } // readMCvals.

   //==================================================================================


   bool addSusyToExpectedObs() {

      //--- Get SUSY inputs, if requested.


         ifstream insusystr ;
         insusystr.open( susyfile ) ;
         if ( !insusystr.good() ) {
            printf("\n\n *** Problem opening susy input file: %s\n\n", susyfile ) ; return false ;
         }

         bool found(false) ;
         while ( insusystr.good() ) {

            float pointMgl ;
            float pointMlsp ;

            int nGen ;

            int ArraySize = 3 + 9*(nBinsMET*nBinsHT*nBinsBjets) ;
            double ArrayContent[ArraySize] ;

            for (int i = 0; i < ArraySize; ++ i) {
              insusystr >> ArrayContent[i];
            }

            pointMgl  = ArrayContent[0] ;
            pointMlsp = ArrayContent[1] ;

            if ( !( fabs(pointMgl - mgl) < 0.1 && fabs(pointMlsp - mlsp) < 0.1 ) ) { continue ; }

            nGen = (int)ArrayContent[2] ;

            int nBins = nBinsMET*nBinsHT*nBinsBjets*3 ;

            float input0lepSusyTotal(0.) ;

            for (int i = 0 ; i < nBinsMET ; i++) {
              for (int j = 0 ; j < nBinsHT ; j++) {
                for (int k = 0 ; k < nBinsBjets ; k++) {

                  float n_0l_raw  = ArrayContent[3 + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;
                  float n_1l_raw  = ArrayContent[4 + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;
                  float n_ldp_raw = ArrayContent[5 + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;

                  float n_0l_correction  = ArrayContent[3 + nBins + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;
                  float n_1l_correction  = ArrayContent[4 + nBins + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;
                  float n_ldp_correction = ArrayContent[5 + nBins + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;

                  susy_0lep[i][j][k] = n_0l_raw * n_0l_correction ;
                  susy_1lep[i][j][k] = n_1l_raw * n_1l_correction ;
                  susy_ldp[i][j][k]  = n_ldp_raw * n_ldp_correction ;

                  input0lepSusyTotal += n_0l_raw * n_0l_correction ;

                }
              }
            }

            susyPoi0lepRatio = input0lepSusyTotal / susy_0lep[0][0][0] ;

            printf("\n\n  pointMgl = %g , pointMlsp = %g,  total 0lep events = %g\n", pointMgl, pointMlsp, input0lepSusyTotal ) ;
            printf("   reset susy 0lep total to %g\n\n", nSusy0lep ) ;

            float newTotal(0.) ;
            for (int i = 0 ; i < nBinsMET ; i++) {
              for (int j = 0 ; j < nBinsHT ; j++) {
                for (int k = 0 ; k < nBinsBjets ; k++) {
                   susy_0lep[i][j][k] = (nSusy0lep/input0lepSusyTotal) * susy_0lep[i][j][k] ;
                   susy_1lep[i][j][k] = (nSusy0lep/input0lepSusyTotal) * susy_1lep[i][j][k] ;
                   susy_ldp[i][j][k]  = (nSusy0lep/input0lepSusyTotal) * susy_ldp[i][j][k] ;
                   printf(" %d,%d,%d  %8.2f\n", i,j,k, susy_0lep[i][j][k] ) ;
                   newTotal += susy_0lep[i][j][k] ;
                }
              }
            }
            printf("\n\n New total is %7.2f\n", newTotal ) ;

            printf("\n\n susy 0lepton 1,1,1 bin value for target is %9.4f\n\n", susy_0lep[0][0][0] ) ;


            found = true ;

            break ;


         } // reading file.

         if ( !found ) {
            printf("\n\n *** Didn't find susy point mgl=%g, mlsp=%g in input susy file %s\n\n", mgl, mlsp, susyfile ) ;
            return false ;
         }


        true_susy_0lep = nSusy0lep ;


        for (int i = 0 ; i < nBinsMET ; i++) {
          for (int j = 0 ; j < nBinsHT ; j++) {
            for (int k = 0 ; k < nBinsBjets ; k++) {
               printf("  0lep (%d,%d,%d) :  SM = %8.1f,  SUSY = %6.1f\n", i,j,k, toy_mean_N_0lep[i][j][k], susy_0lep[i][j][k] ) ;
               printf("  1lep (%d,%d,%d) :  SM = %8.1f,  SUSY = %6.1f\n", i,j,k, toy_mean_N_1lep[i][j][k], susy_1lep[i][j][k] ) ;
               printf("  ldp  (%d,%d,%d) :  SM = %8.1f,  SUSY = %6.1f\n", i,j,k, toy_mean_N_ldp [i][j][k], susy_ldp [i][j][k] ) ;
               toy_mean_N_0lep[i][j][k] += susy_0lep[i][j][k]  ;
               toy_mean_N_1lep[i][j][k] += susy_1lep[i][j][k]  ;
               toy_mean_N_ldp [i][j][k] += susy_ldp [i][j][k]   ;
               printf("\n") ;
            }
               printf(" ---- \n") ;
          }
               printf(" ======== \n") ;
        }




      return true ;

   } // addSusyToExpectedObs.


   //==================================================================================

   RooDataSet* genToyData() {

      RooArgSet obsList ;

      char oname[1000] ;

      int nsel(3) ;
      char selname[3][100] = { "0lep", "1lep", "ldp" } ;

      for ( int si=0; si<nsel; si++ ) {
         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
               for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

                  sprintf( oname, "N_%s_M%d_H%d_%db", selname[si], mbi+1, hbi+1, bbi+1 ) ;
                  RooRealVar* rrv = workspace -> var( oname ) ;
                  if ( rrv == 0x0 ) { printf("\n\n *** can't find variable with name %s\n", oname ) ; return 0x0 ; }

                  if ( si == 0 ) {
                     rrv -> setVal( tran->Poisson( toy_mean_N_0lep[mbi][hbi][bbi] ) ) ;
                  } else if ( si == 1 ) {
                     rrv -> setVal( tran->Poisson( toy_mean_N_1lep[mbi][hbi][bbi] ) ) ;
                  } else if ( si == 2 ) {
                     rrv -> setVal( tran->Poisson( toy_mean_N_ldp [mbi][hbi][bbi] ) ) ;
                  }

                  obsList.add( *rrv ) ;

               } // bbi.
            } // hbi.
         } // mbi.
      }

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            sprintf( oname, "N_Zee_M%d_H%d", mbi+1, hbi+1 ) ;
            RooRealVar* rrv = workspace -> var( oname ) ;
            if ( rrv == 0x0 ) { printf("\n\n *** can't find variable with name %s\n", oname ) ; return 0x0 ; }
            rrv -> setVal( tran->Poisson( toy_mean_N_Zee[mbi][hbi] ) ) ;
            obsList.add( *rrv ) ;
         } // hbi.
      } // mbi.

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            sprintf( oname, "N_Zmm_M%d_H%d", mbi+1, hbi+1 ) ;
            RooRealVar* rrv = workspace -> var( oname ) ;
            if ( rrv == 0x0 ) { printf("\n\n *** can't find variable with name %s\n", oname ) ; return 0x0 ; }
            rrv -> setVal( tran->Poisson( toy_mean_N_Zmm[mbi][hbi] ) ) ;
            obsList.add( *rrv ) ;
         } // hbi.
      } // mbi.

      RooDataSet* rds = new RooDataSet("toyfit_ra2b_observed_rds", "RA2b toy observed data values", obsList ) ;
      rds->add( obsList ) ;


      return rds ;

   } // genToyData.



   //==================================================================================

   bool processToyResult( int ti, RooFitResult* rfr ) {

      if ( rfr == 0x0 ) { return false ; }

      printf(" toy %4d : Fit nSusy 0lep : %4.1f +/- %4.1f (%4.1f, %4.1f)\n", ti,
          susyPoi0lepRatio * (rrv_susy_poi->getVal()),
          susyPoi0lepRatio * (rrv_susy_poi->getError()),
          susyPoi0lepRatio * (rrv_susy_poi->getErrorHi()),
          susyPoi0lepRatio * (rrv_susy_poi->getErrorLo())
          ) ;

      susy_poi_atMinNll = rrv_susy_poi->getVal() ;
      susy_poi_plusErr  = rrv_susy_poi->getErrorHi() ;

      fit_susy_0lep     = susyPoi0lepRatio * (rrv_susy_poi->getVal()) ;
      fit_susy_0lep_err = susyPoi0lepRatio * (rrv_susy_poi->getError()) ;
      fit_susy_0lep_err_low = susyPoi0lepRatio * (rrv_susy_poi->getErrorLo()) ;
      fit_susy_0lep_err_high = susyPoi0lepRatio * (rrv_susy_poi->getErrorHi()) ;
      if ( fit_susy_0lep > true_susy_0lep ) {
         fit_susy_0lep_err_forpull = fit_susy_0lep_err_low ;
      } else {
         fit_susy_0lep_err_forpull = fit_susy_0lep_err_high ;
      }

      fit_ttwj_0lep = 0. ;
      fit_qcd_0lep = 0. ;
      fit_znn_0lep = 0. ;

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               char vname[1000] ;
               RooAbsReal* rar ;

               sprintf( vname, "mu_ttwj_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               if ( mbi==0 && hbi==0 && bbi==0 ) {
                  rar = workspace -> var( vname ) ;
               } else {
                  rar = workspace -> function( vname ) ;
               }
               if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; return false ; }
               fit_ttwj_0lep += rar -> getVal() ;

               sprintf( vname, "mu_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               rar = workspace -> function( vname ) ;
               if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; return false ; }
               fit_qcd_0lep += rar -> getVal() ;

               sprintf( vname, "mu_znn_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               if ( bbi==0 ) {
                  rar = workspace -> var( vname ) ;
               } else {
                  rar = workspace -> function( vname ) ;
               }
               if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; return false ; }
               fit_znn_0lep += rar -> getVal() ;

            } // bbi.
         } // hbi.
      } // mbi.

      fit_covqual_susyfloat = rfr -> covQual() ;

      printf(" toy %4d : Fit total 0lep ttwj : %6.1f  (ave true %6.1f)\n", ti, fit_ttwj_0lep, true_ttwj_0lep ) ;
      printf(" toy %4d : Fit total 0lep qcd  : %6.1f  (ave true %6.1f)\n", ti, fit_qcd_0lep, true_qcd_0lep ) ;
      printf(" toy %4d : Fit total 0lep znn  : %6.1f  (ave true %6.1f)\n", ti, fit_znn_0lep, true_znn_0lep ) ;
      printf(" toy %4d : Fit covariance matrix quality: %d\n", ti, fit_covqual_susyfloat ) ;


      return true ;

   } // processToyResult.


   //==================================================================================


   bool bookTree() {

      char outfile[10000] ;
      sprintf( outfile, "%s/toy-results.root", outputDir ) ;
      ttfile = new TFile( outfile, "recreate" ) ;

      toytt = new TTree( "toytt", "Toy MC study" ) ;


      toytt -> Branch( "true_susy_0lep", &true_susy_0lep, "true_susy_0lep/F" ) ;
      toytt -> Branch( "fit_susy_0lep", &fit_susy_0lep, "fit_susy_0lep/F" ) ;
      toytt -> Branch( "fit_susy_0lep_err", &fit_susy_0lep_err, "fit_susy_0lep_err/F" ) ;
      toytt -> Branch( "fit_susy_0lep_err_low", &fit_susy_0lep_err_low, "fit_susy_0lep_err_low/F" ) ;
      toytt -> Branch( "fit_susy_0lep_err_high", &fit_susy_0lep_err_high, "fit_susy_0lep_err_high/F" ) ;
      toytt -> Branch( "fit_susy_0lep_err_forpull", &fit_susy_0lep_err_forpull, "fit_susy_0lep_err_forpull/F" ) ;
      toytt -> Branch( "true_ttwj_0lep", &true_ttwj_0lep, "true_ttwj_0lep/F" ) ;
      toytt -> Branch( "fit_ttwj_0lep", &fit_ttwj_0lep, "fit_ttwj_0lep/F" ) ;
      toytt -> Branch( "true_qcd_0lep", &true_qcd_0lep, "true_qcd_0lep/F" ) ;
      toytt -> Branch( "fit_qcd_0lep", &fit_qcd_0lep, "fit_qcd_0lep/F" ) ;
      toytt -> Branch( "true_znn_0lep", &true_znn_0lep, "true_znn_0lep/F" ) ;
      toytt -> Branch( "fit_znn_0lep", &fit_znn_0lep, "fit_znn_0lep/F" ) ;
      toytt -> Branch( "fit_covqual_susyfloat", &fit_covqual_susyfloat, "fit_covqual_susyfloat/I" ) ;

      if ( doSignif ) {
         toytt -> Branch( "fit_susy_signif", &fit_susy_signif, "fit_susy_signif/F" ) ;
      }

      if ( doUL ) {
         toytt -> Branch( "fit_susy_ul", &fit_susy_ul, "fit_susy_ul/F" ) ;
         toytt -> Branch( "fit_susy_ts_at_ul", &fit_susy_ts_at_ul, "fit_susy_ts_at_ul/F" ) ;
      }

      return true ;

   } // bookTree.


   //==================================================================================


   float findUL( RooDataSet* toyds ) {

      //-- Search for the 1-sided 95% CL UL in two steps of linear interpolation.
      //   Initially, use testStat at susy = +1.5 sigma and +2.0 sigma.

      float susy_poi_a = susy_poi_atMinNll + 1.5 * susy_poi_plusErr ;
      float susy_poi_b = susy_poi_atMinNll + 2.0 * susy_poi_plusErr ;

      rrv_susy_poi->setVal( susy_poi_a ) ;
      rrv_susy_poi->setConstant( kTRUE ) ;
      RooFitResult* fitResult_a = likelihood -> fitTo( *toyds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
      double minNll_a = fitResult_a->minNll() ;
      double testStat_a = 2.*( minNll_a - minNllSusyFloat ) ;
      delete fitResult_a ;

      rrv_susy_poi->setVal( susy_poi_b ) ;
      rrv_susy_poi->setConstant( kTRUE ) ;
      RooFitResult* fitResult_b = likelihood -> fitTo( *toyds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
      double minNll_b = fitResult_b->minNll() ;
      double testStat_b = 2.*( minNll_b - minNllSusyFloat ) ;
      delete fitResult_b ;




      double susy_poi_g1 = susy_poi_a + (2.70 - testStat_a) * (susy_poi_b-susy_poi_a)/(testStat_b-testStat_a) ;
      printf("  Test stat at +1.5 sigma = %5.2f, +2.0 sigma = %5.2f.  Will try %5.2f sigma.\n",
         testStat_a, testStat_b, (susy_poi_g1-susy_poi_atMinNll)/susy_poi_plusErr ) ;


      rrv_susy_poi->setVal( susy_poi_g1 ) ;
      rrv_susy_poi->setConstant( kTRUE ) ;
      RooFitResult* fitResult_g1 = likelihood -> fitTo( *toyds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
      double minNll_g1 = fitResult_g1->minNll() ;
      double testStat_g1 = 2.*( minNll_g1 - minNllSusyFloat ) ;
      delete fitResult_g1 ;
      printf("  Test stat at first guess is %5.2f\n", testStat_g1 ) ;





      double susy_poi_g2 = susy_poi_a + (2.70 - testStat_a) * (susy_poi_g1-susy_poi_a)/(testStat_g1-testStat_a) ;
      printf("  Test stat at +1.5 sigma = %5.2f, +%5.2f sigma = %5.2f.  Will try %5.2f sigma.\n",
         testStat_a,
         (susy_poi_g1-susy_poi_atMinNll)/susy_poi_plusErr,
         testStat_g1, (susy_poi_g2-susy_poi_atMinNll)/susy_poi_plusErr ) ;

      rrv_susy_poi->setVal( susy_poi_g2 ) ;
      rrv_susy_poi->setConstant( kTRUE ) ;
      RooFitResult* fitResult_g2 = likelihood -> fitTo( *toyds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
      double minNll_g2 = fitResult_g2->minNll() ;
      double testStat_g2 = 2.*( minNll_g2 - minNllSusyFloat ) ;
      delete fitResult_g2 ;
      printf("  Test stat at second guess is %5.2f\n", testStat_g2 ) ;

      fit_susy_ul = susyPoi0lepRatio * susy_poi_g2 ;
      fit_susy_ts_at_ul = testStat_g2 ;

      printf("  SUSY upper limit is : %5.2f 0lep events.\n", fit_susy_ul ) ;

      return fit_susy_ul ;


   } // findUL.

   //==================================================================================

   bool setReinitValues( RooFitResult* fitResult ) {

      const RooArgList floatPars = fitResult -> floatParsFinal() ;

      printf("\n\n ==== Will reinitialize to these values before every toy.\n\n") ;

      nFloatParInitVal = 0 ;

      TIterator* parIter = floatPars.createIterator() ;
      while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {

         sprintf( floatParName[nFloatParInitVal], "%s", par->GetName() ) ;
         floatParInitVal[nFloatParInitVal] = par->getVal() ;

         nFloatParInitVal++ ;

         printf(" %20s : %8.2f\n", par->GetName(), par->getVal() ) ;

      }
      printf("\n\n") ;

      return true ;

   } // setReinitValues.

   //==================================================================================

   bool reinitFloatPars() {

      for ( int pi=0; pi<nFloatParInitVal; pi++ ) {
         RooRealVar* par = workspace->var( floatParName[pi] ) ;
         if ( par == 0x0 ) {
            printf("\n\n *** can't find float par %s in workspace.\n\n", floatParName[pi] ) ;
            return false ;
         }
         par->setVal( floatParInitVal[pi] ) ;
      } // pi.
      rrv_susy_poi -> setVal( 0. ) ;

      return true ;

   } // reinitFloatPars.

   //==================================================================================




