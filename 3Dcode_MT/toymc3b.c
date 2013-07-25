#include "TROOT.h"
#include "TSystem.h"
#include "TRandom.h"
#include "ra2bRoostatsClass3D_3b.c"
#include "updateFileValue.c"
#include "getFileValue.c"


   // this has to be set by hand ...
   static const int nBinsVar1  = 3 ;
   static const int nBinsVar2   = 4 ;
   static const int nBinsBjets = 3 ;   


   //-- Global variables.

   char datfile[10000] ;
   char susyfile[10000] ;
   char deffdbtagfile[10000] ;
   char outputDir[10000] ;
   char mcvals_rootfile[10000] ;

   int   nFloatParInitVal ;
   char  floatParName[5000][100] ;
   float floatParInitVal[5000] ;

   float mgl, mlsp ;

   float nSusy0lep ;

   int qcdModelIndex ;

   RooWorkspace* workspace ;

   RooRealVar* rrv_susy_poi ;
   RooAbsPdf* likelihood ;

   RooRealVar* rrv_qcd_0lepLDP_ratio ;
   RooRealVar* rrv_qcd_0lepLDP_ratio_H1 ;
   RooRealVar* rrv_qcd_0lepLDP_ratio_H2 ;
   RooRealVar* rrv_qcd_0lepLDP_ratio_H3 ;
   RooRealVar* rrv_qcd_0lepLDP_ratio_H4 ;
   RooRealVar* rrv_SFqcd_met2 ;
   RooRealVar* rrv_SFqcd_met3 ;
   RooRealVar* rrv_SFqcd_met4 ;
   RooRealVar* rrv_SFqcd_nb2 ;
   RooRealVar* rrv_SFqcd_nb3 ;
   RooRealVar* rrv_SFqcd_nb4 ;

   RooRealVar* rrv_ttwj_0lep1lep_ratio ;
   RooRealVar* rrv_ttwj_1lepSig1lep_ratio ;

   float sf_ttwj[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float sf_qcd[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float ttwj_mc_ldpover0lep_ratio[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float  znn_mc_ldpover0lep_ratio[nBinsVar1][nBinsVar2][nBinsBjets] ;


   //--- Output ttree variables.

   int   Nobs_0lep[nBinsVar1][nBinsVar2][nBinsBjets] ;
   int   Nobs_1lepSig[nBinsVar1][nBinsVar2][nBinsBjets] ;
   int   Nobs_1lep[nBinsVar1][nBinsVar2][nBinsBjets] ;
   int   Nobs_ldp [nBinsVar1][nBinsVar2][nBinsBjets] ;
   int   Nobs_Zee [nBinsVar1][nBinsVar2] ;
   int   Nobs_Zmm [nBinsVar1][nBinsVar2] ;

   float toy_mean_N_0lep[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float toy_mean_N_1lepSig[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float toy_mean_N_1lep[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float toy_mean_N_ldp [nBinsVar1][nBinsVar2][nBinsBjets] ;
   float toy_mean_N_Zee [nBinsVar1][nBinsVar2] ;
   float toy_mean_N_Zmm [nBinsVar1][nBinsVar2] ;

   float fit_susy_0lep_err ;
   float fit_susy_0lep_err_low ;
   float fit_susy_0lep_err_high ;
   float fit_susy_0lep_err_forpull ;
   float fit_susy_0lep ;
   float fit_susy_0lep_wsfs ;
   float true_susy_0lep ;

   float fit_susy_0lep_3da[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_ttwj_0lep_3da[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_qcd_0lep_3da [nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_znn_0lep_3da [nBinsVar1][nBinsVar2][nBinsBjets] ;

   float fit_sf_ttwj_0lep_3da[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_sf_qcd_0lep_3da [nBinsVar1][nBinsVar2][nBinsBjets] ;

   float fit_susy_0lep_atUL_3da[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_ttwj_0lep_atUL_3da[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_qcd_0lep_atUL_3da [nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_znn_0lep_atUL_3da [nBinsVar1][nBinsVar2][nBinsBjets] ;

   float fit_sf_ttwj_0lep_atUL_3da[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_sf_qcd_0lep_atUL_3da [nBinsVar1][nBinsVar2][nBinsBjets] ;

   float fit_ttwj_0lep1lep_ratio_atUL;
   float fit_ttwj_1lepSig1lep_ratio_atUL;

   float fit_susy_0lep_at0susy_3da[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_ttwj_0lep_at0susy_3da[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_qcd_0lep_at0susy_3da [nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_znn_0lep_at0susy_3da [nBinsVar1][nBinsVar2][nBinsBjets] ;

   float fit_sf_ttwj_0lep_at0susy_3da[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float fit_sf_qcd_0lep_at0susy_3da [nBinsVar1][nBinsVar2][nBinsBjets] ;

   float fit_ttwj_0lep1lep_ratio_at0susy;
   float fit_ttwj_1lepSig1lep_ratio_at0susy;

   float fit_qcd_0lepLDP_ratio_at0susy ;
   float fit_qcd_0lepLDP_ratio_H1_at0susy ;
   float fit_qcd_0lepLDP_ratio_H2_at0susy ;
   float fit_qcd_0lepLDP_ratio_H3_at0susy ;
   float fit_qcd_0lepLDP_ratio_H4_at0susy ;
   float fit_SFqcd_met2_at0susy ;
   float fit_SFqcd_met3_at0susy ;
   float fit_SFqcd_met4_at0susy ;
   float fit_SFqcd_nb2_at0susy ;
   float fit_SFqcd_nb3_at0susy ;
   float fit_SFqcd_nb4_at0susy ;

   float fit_qcd_0lepLDP_ratio_atUL ;
   float fit_qcd_0lepLDP_ratio_H1_atUL ;
   float fit_qcd_0lepLDP_ratio_H2_atUL ;
   float fit_qcd_0lepLDP_ratio_H3_atUL ;
   float fit_qcd_0lepLDP_ratio_H4_atUL ;
   float fit_SFqcd_met2_atUL ;
   float fit_SFqcd_met3_atUL ;
   float fit_SFqcd_met4_atUL ;
   float fit_SFqcd_nb2_atUL ;
   float fit_SFqcd_nb3_atUL ;
   float fit_SFqcd_nb4_atUL ;

   float mcval_ttwj_0lep_3da[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float mcval_qcd_0lep_3da [nBinsVar1][nBinsVar2][nBinsBjets] ;
   float mcval_znn_0lep_3da [nBinsVar1][nBinsVar2][nBinsBjets] ;

   float susy_0lep   [nBinsVar1][nBinsVar2][nBinsBjets] ; // this is called mcval_susy_0lep_3da in the ttree.
   float susy_1lepSig[nBinsVar1][nBinsVar2][nBinsBjets] ;
   float susy_1lep   [nBinsVar1][nBinsVar2][nBinsBjets] ;
   float susy_ldp    [nBinsVar1][nBinsVar2][nBinsBjets] ;

   float fit_chi2_overall ;
   float fit_chi2_obs ;
   float fit_chi2_np ;
   float fit_chi2_prob ;

   float fit_qcd_0lepLDP_ratio ;
   float fit_qcd_0lepLDP_ratio_H1 ;
   float fit_qcd_0lepLDP_ratio_H2 ;
   float fit_qcd_0lepLDP_ratio_H3 ;
   float fit_qcd_0lepLDP_ratio_H4 ;
   float fit_SFqcd_met2 ;
   float fit_SFqcd_met3 ;
   float fit_SFqcd_met4 ;
   float fit_SFqcd_nb2 ;
   float fit_SFqcd_nb3 ;
   float fit_SFqcd_nb4 ;

   float fit_qcd_0lepLDP_ratio_err ;
   float fit_qcd_0lepLDP_ratio_H1_err ;
   float fit_qcd_0lepLDP_ratio_H2_err ;
   float fit_qcd_0lepLDP_ratio_H3_err ;
   float fit_qcd_0lepLDP_ratio_H4_err ;
   float fit_SFqcd_met2_err ;
   float fit_SFqcd_met3_err ;
   float fit_SFqcd_met4_err ;
   float fit_SFqcd_nb2_err ;
   float fit_SFqcd_nb3_err ;
   float fit_SFqcd_nb4_err ;

   float fit_ttwj_0lep1lep_ratio ;
   float fit_ttwj_0lep1lep_ratio_err ;

   float fit_ttwj_1lepSig1lep_ratio ;
   float fit_ttwj_1lepSig1lep_ratio_err ;


   int   fit_covqual_susyfloat ;

   float fit_susy_signif ;

   float fit_susy_ul ;
   float fit_susy_ts_at_ul ;


   //--- End output ttree variables.

   double initFit_R_ttwj_0lep_over_1lep ;
   double initFit_qcd_0lepLDP_ratio ;
   double initFit_qcd_0lepLDP_ratio_HT[20] ;
   double initFit_SFqcd_met[20] ;
   double initFit_SFqcd_nb[20] ;

   double minNllSusyFloat ;

   double susy_poi_atMinNll ;
   double susy_poi_plusErr ;


   bool doSignif ;
   bool doUL ;

   bool inputObservablesArePostTrigger ;


   bool useExpected0lep ;

   double trigeff_0lep[nBinsVar1][nBinsVar2] ;
   double trigeff_1lep[nBinsVar1][nBinsVar2] ;


   bool ignoreBin[nBinsVar1][nBinsVar2] ;
   bool blind0lepBin[nBinsVar1][nBinsVar2][nBinsBjets] ;
   int  num_ignoreBin(0) ;
   int  num_blind0lepBin(0) ;


   //-- Prototypes

   bool readMCvals() ;
   bool addSusyToExpectedObs() ;
   RooDataSet* genToyData( bool noFluctuations = false ) ;
   bool processToyResult( int ti, RooFitResult* fitResult ) ;
   bool bookTree() ;
   float findUL( RooDataSet* toyds ) ;
   bool setReinitValues( RooFitResult* fitResult ) ;
   bool reinitFloatPars() ;
   bool saveToyDatfile( int toyIndex, RooDataSet* toyds ) ;
   bool readAndSetMCVals() ;
   bool setSFVals() ;
   double getChi2Obs( const char* obsname, const char* modelname ) ;
   double getChi2GausNP( const char* npname ) ;
   double getChi2BetaNP( const char* npname ) ;
   bool getBetaModeRMS( const char* parName, double &mode, double &rms, double &alpha, double &beta ) ;


   TRandom* tran ;

   TTree* toytt ;
   TFile* ttfile ;

   char blindBinsList[10000] ;

   bool blindStudy ;

   //=====================================

   void toymc3b( const char* input_datfile = "datfiles/Input-met4-ht4-wsyst1.dat",
		 const char* input_susyfile = "datfiles/Susy-mgl900-mlsp300-met4-ht4.dat",
		 double input_mgl=900, double input_mlsp=300.,
		 const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met4-ht4.dat",
		 const char* input_lightmistagfile = "datfiles/dummy_DeffDbtag-met4-ht4.dat",
		 double input_nSusy0lep = 60.,
		 const char* input_outputDir = "output-toymc3b",
		 int nToy = 10,
		 const char* input_mcvals_rootfile = "rootfiles/gi-plots-met4-ht4-newMC.root",
		 bool input_useExpected0lep = false,
		 int input_qcdModelIndex = 3,
		 bool input_doSignif = false,
		 bool input_doUL = false,
		 const char* input_blindBinsList = "null",
		 bool input_inputObservablesArePostTrigger = true,
		 bool constrainBjetShape = false,
		 bool floatSLSigRatios = false,
		 const char* systfile1 = "systFile1.txt",
		 const char* systfile2 = "systFile2.txt",
		 const char* systfile3 = "systFile3.txt",
		 const char* wjets_syst = "wjetsfile.txt",
		 const char* singletop_syst = "singletopfile.txt"
		 ) {
     
       char command[10000] ;

       sprintf( datfile, "%s", input_datfile ) ;
       sprintf( susyfile, "%s", input_susyfile ) ;
       sprintf( outputDir, "%s", input_outputDir ) ;
       sprintf( deffdbtagfile, "%s", input_deffdbtagfile ) ;
       sprintf( mcvals_rootfile, "%s", input_mcvals_rootfile ) ;
       sprintf( blindBinsList, "%s", input_blindBinsList ) ;


      //-- Hardwire in to ignore the highest Var1 bin in the lowest Var2 bin.
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
            ignoreBin[mbi][hbi] = false ;
         } // hbi
      } // mbi
      ignoreBin[3][0] = true ; 

      printf("\n\n *** Ignoring these bins in the analysis.\n\n") ;
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
            if ( ignoreBin[mbi][hbi] ) {
               num_ignoreBin ++ ;
               printf("  Var1 %d, Var2 %d\n", mbi+1, hbi+1 ) ;
            }
         } // hbi
      } // mbi
      printf("\n\n\n") ;


     //--- read in blind bins, if file given.

      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               blind0lepBin[mbi][hbi][bbi] = false ;
            } // bbi.
         } // hbi.
      } // mbi.

      blindStudy = false ;

      if ( strcmp( blindBinsList, "null" ) != 0 ) {
         sprintf( command, "ls %s >& /dev/null", blindBinsList ) ;
         int returnstat = gSystem->Exec( command ) ;
         if ( returnstat != 0 ) {
            printf("\n\n *** Given non-null blindBinsList file, but it doesn't exist: %s\n\n", blindBinsList ) ;
            return ;
         }
         FILE* bbl_file ;
         if ( (bbl_file = fopen( blindBinsList, "r"))==NULL ) {
            printf("\n\n *** Problem opening blindBinsList file: %s\n\n", blindBinsList ) ;
            return ;
         }
         blindStudy = true ;
         printf("\n\n Reading blindBinsList file: %s\n\n", blindBinsList ) ;
         while ( ! feof(bbl_file) ) {
            char binlabel[1000] ;
            fscanf( bbl_file, "%s", binlabel ) ;
            if ( feof(bbl_file) ) break ;
            int bmb(-1), bhb(-1), bbb(-1) ;
            int rv = sscanf( binlabel, "N_0lep_M%d_H%d_%db", &bmb, &bhb, &bbb ) ;
            if ( rv == 0 ) {
               printf("\n\n *** bad bin label format: %s\n\n", binlabel ) ;
               return ;
            }
            if ( bmb < 1 || bhb < 1 || bbb < 1 || bmb > nBinsVar1 || bhb > nBinsVar2 || bbb > nBinsBjets ) {
               printf("\n\n *** bin index out of range: Var1=%d, Var2=%d, nb=%d\n\n", bmb, bhb, bbb ) ;
               return ;
            }
            num_blind0lepBin++ ;
            blind0lepBin[bmb-1][bhb-1][bbb-1] = true ;
            printf("   bin label %s : will blind Var1=%d, Var2=%d, nb=%d 0lep bin.\n", binlabel, bmb, bhb, bbb ) ;
         } // still reading?
      } // process blindBinsList file?





       if ( strcmp( blindBinsList, "null") == 0 ) {
          blindStudy = false ;
       } else {
          printf("\n\n *** blind bins list file set to %s.  Will fix signal to zero in fits.\n", blindBinsList ) ;
          blindStudy = true ;
       }

       sprintf( command, "mkdir -p %s", outputDir ) ;
       gSystem->Exec( command ) ;

       sprintf( command, "cp %s %s", datfile, outputDir ) ;
       gSystem->Exec( command ) ;

       sprintf( command, "cp %s %s", susyfile, outputDir ) ;
       gSystem->Exec( command ) ;

       sprintf( command, "cp %s %s", deffdbtagfile, outputDir ) ;
       gSystem->Exec( command ) ;

       sprintf( command, "cp %s %s", mcvals_rootfile, outputDir ) ;
       gSystem->Exec( command ) ;


       inputObservablesArePostTrigger = input_inputObservablesArePostTrigger ;


       mgl  = input_mgl ;
       mlsp = input_mlsp ;

       nSusy0lep = input_nSusy0lep ;

       true_susy_0lep = 0. ;

       doSignif = input_doSignif ;
       doUL = input_doUL ;

       useExpected0lep = input_useExpected0lep ;


       tran = new TRandom(12345) ;
       //tran = new TRandom() ;
       //tran->SetSeed(0);


       qcdModelIndex = input_qcdModelIndex ;



       //--- Create workspace.

       ra2bRoostatsClass3D_3b ra2b ;

       char wsfilename[10000] ;
       sprintf( wsfilename, "%s/ws.root", outputDir ) ;
       ra2b.initialize( input_datfile, input_susyfile, mgl, mlsp, false, 0., input_deffdbtagfile, input_lightmistagfile, qcdModelIndex, wsfilename, 
			blindBinsList, constrainBjetShape, floatSLSigRatios, systfile1, systfile2, systfile3, wjets_syst, singletop_syst ) ;
       TFile wsfile( wsfilename ) ;
       workspace = (RooWorkspace*) wsfile.Get("ws") ;
       if ( workspace == 0x0 ) {
          printf("\n\n *** Can't find the workspace.\n\n" ) ; return ;
       }








       //--- Read in ttwj scale factors (used if useExpected0lep is true).

       if ( !setSFVals() ) {
          printf("\n\n *** Problem reading in scale factors from %s.\n\n", datfile ) ;
       }







       //--- Access some workspace stuff needed later.

       rrv_susy_poi = workspace -> var( "mu_susy_all0lep" ) ;
       if ( rrv_susy_poi == 0x0 ) {
          printf("\n\n *** can't find susy poi mu_susy_all0lep.\n\n") ; return ;
       }

       if ( qcdModelIndex < 2 || qcdModelIndex > 5 ) {
          printf("\n\n Unsupported QCD model index %d.\n\n", qcdModelIndex ) ;
          return ;
       }
       if ( qcdModelIndex == 2 || qcdModelIndex == 4 || qcdModelIndex == 5 ) {
          if ( nBinsVar2 >=1 ) {
             rrv_qcd_0lepLDP_ratio_H1 = workspace -> var( "qcd_0lepLDP_ratio_H1" ) ;
             if ( rrv_qcd_0lepLDP_ratio_H1 == 0x0 ) {
                printf("\n\n *** can't find qcd_0lepLDP_ratio_H1.\n\n") ; return ;
             }
          }
          if ( nBinsVar2 >=2 ) {
             rrv_qcd_0lepLDP_ratio_H2 = workspace -> var( "qcd_0lepLDP_ratio_H2" ) ;
             if ( rrv_qcd_0lepLDP_ratio_H2 == 0x0 ) {
                printf("\n\n *** can't find qcd_0lepLDP_ratio_H2.\n\n") ; return ;
             }
          }
          if ( nBinsVar2 >=3 ) {
             rrv_qcd_0lepLDP_ratio_H3 = workspace -> var( "qcd_0lepLDP_ratio_H3" ) ;
             if ( rrv_qcd_0lepLDP_ratio_H3 == 0x0 ) {
                printf("\n\n *** can't find qcd_0lepLDP_ratio_H3.\n\n") ; return ;
             }
          }
          if ( nBinsVar2 >=4 ) {
             rrv_qcd_0lepLDP_ratio_H4 = workspace -> var( "qcd_0lepLDP_ratio_H4" ) ;
             if ( rrv_qcd_0lepLDP_ratio_H4 == 0x0 ) {
                printf("\n\n *** can't find qcd_0lepLDP_ratio_H4. nBinsVar2=%d\n\n", nBinsVar2) ; return ;
             }
          }
       } 
       if ( qcdModelIndex == 3 ) {
          rrv_qcd_0lepLDP_ratio = workspace -> var( "qcd_0lepLDP_ratio" ) ;
          if ( rrv_qcd_0lepLDP_ratio == 0x0 ) {
             printf("\n\n *** can't find qcd_0lepLDP_ratio.\n\n") ; return ;
          }
       }
       if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {
          if ( nBinsVar1 >=2 ) {
             rrv_SFqcd_met2 = workspace -> var( "SFqcd_met2" ) ;
             if ( rrv_SFqcd_met2 == 0x0 ) {
	       printf("\n\n *** can't find SFqcd_met2.\n\n") ; return ;
             }
          }
          if ( nBinsBjets >=2 ) {
             rrv_SFqcd_nb2 = workspace -> var( "SFqcd_nb2" ) ;
             if ( rrv_SFqcd_nb2 == 0x0 ) {
                printf("\n\n *** can't find SFqcd_nb2.\n\n") ; return ;
             }
          }
          if ( qcdModelIndex == 4 ) {
             if ( nBinsVar1 >=3 ) {
                rrv_SFqcd_met3 = workspace -> var( "SFqcd_met3" ) ;
                if ( rrv_SFqcd_met3 == 0x0 ) {
		  printf("\n\n *** can't find SFqcd_met3.\n\n") ; //return ;
                }
             }
             if ( nBinsVar1 >=4 ) {
                rrv_SFqcd_met4 = workspace -> var( "SFqcd_met4" ) ;
                if ( rrv_SFqcd_met4 == 0x0 ) {
		  printf("\n\n *** can't find SFqcd_met4.\n\n") ; //return ;
                }
             }
             if ( nBinsBjets >=3 ) {
                rrv_SFqcd_nb3 = workspace -> var( "SFqcd_nb3" ) ;
                if ( rrv_SFqcd_nb3 == 0x0 ) {
		  printf("\n\n *** can't find SFqcd_nb3.\n\n") ; //return ;
                }
             }
	     //             if ( nBinsBjets >=4 ) {
             //   rrv_SFqcd_nb4 = workspace -> var( "SFqcd_nb4" ) ;
             //   if ( rrv_SFqcd_nb4 == 0x0 ) {
             //      printf("\n\n *** can't find SFqcd_nb4.\n\n") ; return ;
             //   }
             //}
          }
       }


       rrv_ttwj_0lep1lep_ratio = workspace -> var( "ttwj_0lep1lep_ratio" ) ;
       if ( rrv_ttwj_0lep1lep_ratio == 0x0 ) {
          printf("\n\n *** can't find rrv_ttwj_0lep1lep_ratio.\n\n") ; return ;
       }

       rrv_ttwj_1lepSig1lep_ratio = workspace -> var( "ttwj_slSigsl_ratio" ) ;
       if ( rrv_ttwj_1lepSig1lep_ratio == 0x0 ) {
          printf("\n\n *** can't find rrv_ttwj_1lepSig1lep_ratio.\n\n") ; return ;
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









       //--- Do an another initial fit to determine reasonable starting values for all floating parameters.
       //    Do it from dataset at MC values with no fluctuations, including susy.
       
       RooDataSet* toydsNofluctuations = genToyData(true) ;

       rrv_susy_poi -> setVal( nSusy0lep ) ;
       rrv_susy_poi -> setConstant( kFALSE ) ;

       printf("\n\n\n ===== Starting fit to determine reinitialization values...\n\n\n") ;

       RooFitResult* mcfitResult2 = likelihood->fitTo( *toydsNofluctuations, Save(true), PrintLevel(0) ) ;

       if ( ! setReinitValues( mcfitResult2 ) ) {
          printf("\n\n *** Problem in collecting inital vals for floating parameters from MC fit.\n\n") ; return ;
       }








       //--- set MC vals to be saved in output ttree.

       if ( !useExpected0lep ) {
          if ( !readAndSetMCVals() ) {
             printf("\n\n *** Problem reading MC vals from file %s.\n\n", mcvals_rootfile ) ;
          }
       }










       //--- Create tree and define tree variables.
       bookTree() ;




       //--- Loop over toy experiments.

       for ( int ti=0; ti<nToy; ti++ ) {

          printf("\n\n\n\n ====== Begin toy experiment %d out of %d\n\n\n", ti, nToy ) ;




          //-- Do a fit with susy floating.  Get MINOS errors on susy yield.

          RooDataSet* toyds = genToyData() ;
          if ( toyds == 0x0 ) { printf("\n\n *** Fail!\n\n") ; return ; }
          toyds->printMultiline( cout, 1, kTRUE, "" ) ;



          saveToyDatfile( ti, toyds ) ;



          RooArgSet ras ;
          ras.add( *rrv_susy_poi ) ;
          if ( !blindStudy ) {
             rrv_susy_poi->setConstant( kFALSE ) ;
          } else {
             rrv_susy_poi->setVal( 0.0 ) ;
             rrv_susy_poi->setConstant( kTRUE ) ;
          }

          if ( ! reinitFloatPars() ) {
             printf("\n\n *** problem reinitializing floating pars.\n\n") ; return ;
          }

          RooFitResult* fitResult ;
          if ( !blindStudy ) {
             fitResult = likelihood -> fitTo( *toyds, Save(true), Hesse(true), Minos(ras), Strategy(1), PrintLevel(0) ) ;
          } else {
             fitResult = likelihood -> fitTo( *toyds, Save(true), Hesse(true), Strategy(1), PrintLevel(0) ) ;
          }
          minNllSusyFloat = fitResult->minNll() ;

          processToyResult( ti, fitResult ) ;




          //-- Calculate significance from delta log likelihood, if requested.

          if ( doSignif && !blindStudy ) {
             double susy_yield = rrv_susy_poi->getVal();
             rrv_susy_poi->setVal(0.) ;
             rrv_susy_poi->setConstant( kTRUE ) ;
             RooFitResult* fitResult0susy = likelihood -> fitTo( *toyds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
             double minNll0Susy = fitResult0susy->minNll() ;

             double testStat = 2.*( minNll0Susy - minNllSusyFloat ) ;
             fit_susy_signif = 0. ;
             if ( testStat > 0 ) { 
               if (susy_yield>0) fit_susy_signif = sqrt( testStat ) ;
               else  fit_susy_signif = -sqrt( testStat ) ;
             }
             printf("  test stat = %5.2f,  significance = %5.2f\n", testStat, fit_susy_signif ) ;

             fit_ttwj_0lep1lep_ratio_at0susy = rrv_ttwj_0lep1lep_ratio->getVal();

             if ( qcdModelIndex == 2 || qcdModelIndex == 4 || qcdModelIndex == 5 ) {
                if ( nBinsVar2 >= 1 ) { fit_qcd_0lepLDP_ratio_H1_at0susy     = rrv_qcd_0lepLDP_ratio_H1 -> getVal() ; }
                if ( nBinsVar2 >= 2 ) { fit_qcd_0lepLDP_ratio_H2_at0susy     = rrv_qcd_0lepLDP_ratio_H2 -> getVal() ; }
                if ( nBinsVar2 >= 3 ) { fit_qcd_0lepLDP_ratio_H3_at0susy     = rrv_qcd_0lepLDP_ratio_H3 -> getVal() ; }
                if ( nBinsVar2 >= 4 ) { fit_qcd_0lepLDP_ratio_H4_at0susy     = rrv_qcd_0lepLDP_ratio_H4 -> getVal() ; }
             }
             if ( qcdModelIndex == 3 ) { fit_qcd_0lepLDP_ratio_at0susy     = rrv_qcd_0lepLDP_ratio -> getVal() ; }
             if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {
                if ( nBinsVar1 >=2 ) { fit_SFqcd_met2_at0susy     = rrv_SFqcd_met2 -> getVal() ; }
                if ( nBinsBjets >=2 ) { fit_SFqcd_nb2_at0susy     = rrv_SFqcd_nb2 -> getVal() ; }
                if ( qcdModelIndex == 4 ) {
		  //if ( nBinsVar1 >=3 ) { fit_SFqcd_met3_at0susy     = rrv_SFqcd_met3 -> getVal() ; }
		  //if ( nBinsVar1 >=4 ) { fit_SFqcd_met4_at0susy     = rrv_SFqcd_met4 -> getVal() ; }
		  //if ( nBinsBjets >=3 ) { fit_SFqcd_nb3_at0susy     = rrv_SFqcd_nb3 -> getVal() ; }
                   if ( nBinsBjets >=4 ) { fit_SFqcd_nb4_at0susy     = rrv_SFqcd_nb4 -> getVal() ; }
                }
             }


             for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
                for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
                   for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
		      //skip highest Var1, lowest Var2 bin, which is omitted from the analysis
		      if (mbi==3&&hbi==0) {
	                 fit_susy_0lep_at0susy_3da[mbi][hbi][bbi] = 0.0;
	                 fit_ttwj_0lep_at0susy_3da[mbi][hbi][bbi] = 0.0;
	                 fit_qcd_0lep_at0susy_3da[mbi][hbi][bbi] = 0.0;
	                 fit_znn_0lep_at0susy_3da[mbi][hbi][bbi] = 0.0;
	                 fit_sf_ttwj_0lep_at0susy_3da[mbi][hbi][bbi] = 1.0;
	                 fit_sf_qcd_0lep_at0susy_3da[mbi][hbi][bbi] = 1.0;
		      } else {

                         char vname[1000] ;
                         RooAbsReal* rar ;
                         RooAbsReal* effsf  ;
                         RooAbsReal* btagsf  ;


                         sprintf( vname, "mu_susy_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                         rar = workspace -> function( vname ) ;
                         if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; }
                         sprintf( vname, "btageff_sf_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                         btagsf = (RooAbsReal*) workspace -> obj( vname ) ;
                         if ( btagsf == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ;  }
                         sprintf( vname, "eff_sf_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                         effsf = (RooAbsReal*) workspace -> obj( vname ) ;
                         if ( effsf == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; }
                         double nsusy = ( btagsf -> getVal() ) * ( effsf -> getVal() ) * ( rar -> getVal() ) ;

                         fit_susy_0lep_wsfs += nsusy ;

                         fit_susy_0lep_at0susy_3da[mbi][hbi][bbi] = nsusy ;

                         sprintf( vname, "mu_ttwj_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                         rar = (RooAbsReal*) workspace -> obj( vname ) ;
                         if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ;  }
                         double nttwj = rar -> getVal() ;
                         fit_ttwj_0lep_at0susy_3da[mbi][hbi][bbi] = nttwj ;


                         sprintf( vname, "mu_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                         rar = workspace -> function( vname ) ;
                         if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ;}
                         double nqcd = rar -> getVal() ;
                         fit_qcd_0lep_at0susy_3da[mbi][hbi][bbi] = nqcd ;

                         sprintf( vname, "mu_znn_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                         rar = (RooAbsReal*) workspace -> obj( vname ) ;
                         if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; }
                         double nznn = rar -> getVal() ;
                         fit_znn_0lep_at0susy_3da[mbi][hbi][bbi] = nznn ;

                         sprintf( vname, "sf_ttwj_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                         rar = (RooAbsReal*) workspace -> obj( vname ) ;
                         if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ;  }
                         fit_sf_ttwj_0lep_at0susy_3da[mbi][hbi][bbi] = rar -> getVal() ;

                         sprintf( vname, "sf_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                         rar = (RooAbsReal*) workspace -> obj( vname ) ;
                         if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; }
                         fit_sf_qcd_0lep_at0susy_3da[mbi][hbi][bbi] = rar -> getVal() ; 	      

                      }
		   }
                }
             }


             delete fitResult0susy ;

             rrv_susy_poi->setConstant( kFALSE ) ;

          } // doSignif?





          //-- Calculate one-sided upper limit on susy signal, if requested.

          if ( doUL && !blindStudy ) {

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

      float N_0lep_input [nBinsVar1][nBinsVar2][nBinsBjets] ;
      float N_1lepSig    [nBinsVar1][nBinsVar2][nBinsBjets] ;
      float N_1lep       [nBinsVar1][nBinsVar2][nBinsBjets] ;
      float N_ldp        [nBinsVar1][nBinsVar2][nBinsBjets] ;
      float N_mc_ldp     [nBinsVar1][nBinsVar2][nBinsBjets] ;
      float N_Zee        [nBinsVar1][nBinsVar2] ;
      float N_Zmm        [nBinsVar1][nBinsVar2] ;

      float acc_Zee      [nBinsVar1] ;
      float acc_Zmm      [nBinsVar1] ;

      float Z_ee_eff ;
      float Z_mm_eff ;

      float knn[nBinsVar1][nBinsBjets] ;

      float Z_ee_pur(0.) ;
      float Z_mm_pur(0.) ;

      // float R_passfail[nBinsVar2] ;


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
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
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

      //--- Read in N_1lepSig lines
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               nread = fscanf( infp, "%s %g", pname, &pval ) ;
               ///// printf("  %d : %s %g\n", nread, pname, pval ) ;
               int metbin, htbin, nbbin ;
               int nmatch = sscanf( pname, "N_1lepSig_M%d_H%d_%db", &metbin, &htbin, &nbbin ) ;
               if ( nmatch == 3 ) {
                  if ( metbin == mbi+1 && htbin == hbi+1 && nbbin == bbi+1 ) {
                     N_1lepSig[mbi][hbi][bbi] = pval ;
                     printf(" Set N_1lepSig for mbi=%d, hbi=%d, bbi=%d to %g\n", mbi, hbi, bbi, pval ) ;
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
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
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
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
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
      for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
         for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
            nread = fscanf( infp, "%s %g", pname, &pval ) ;
            nread = fscanf( infp, "%s %g", pname, &pval ) ;
         } // bbi.
      } // hbi

      //--- Read in Zee lines
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
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
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
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
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
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
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
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
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
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
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
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
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
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
         if (bbi==0) {
	    for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
	       nread = fscanf( infp, "%s %g", pname, &pval ) ;
               int nbbin ;
               int metbin ;
               int nmatch = sscanf( pname, "knn_%db_M%d", &nbbin, &metbin ) ;
               if ( nmatch == 2 ) {
	           if ( nbbin == bbi+1 ) {
                      knn[mbi][bbi] = pval ;
                      printf(" Set knn for bbi=%d and mbi=%d to %g\n", bbi, mbi, pval ) ;
                   } else {
                      printf(" Expecting nB=%d but found %d.  Label is %s, val is %g, nmatch=%d\n",
                           bbi+1, nbbin, pname, pval, nmatch ) ;
                      return false ;
                   }
               } else {
                  printf("\n\n *** bad input line.  Label = %s, val = %g\n\n", pname, pval) ;
                  return false ;
               }
               nread = fscanf( infp, "%s %g", pname, &pval ) ;  
	    }
	 } else {
	    nread = fscanf( infp, "%s %g", pname, &pval ) ;
            int nbbin ;
            int nmatch = sscanf( pname, "knn_%db", &nbbin ) ;
            if ( nmatch == 1 ) {
               if ( nbbin == bbi+1 ) {
           	  knn[0][bbi] = pval ;
           	  printf(" Set knn for bbi=%d to %g\n", bbi, pval ) ;
               } else {
           	  printf("\n\n *** bad input line.\n\n") ; return false ;
               }
            } else {
               printf("\n\n *** bad input line.\n\n") ; return false ;
            }
            nread = fscanf( infp, "%s %g", pname, &pval ) ;
	 }
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
         Zee_factor[bbi] = ( 5.94 * Z_ee_pur * knn[0][bbi] ) / ( acc_Zee[0] * Z_ee_eff ) ; // ignoring Var1 dependence of acceptance, knn and Var2 dependence of eff
         Zmm_factor[bbi] = ( 5.94 * Z_mm_pur * knn[0][bbi] ) / ( acc_Zmm[0] * Z_mm_eff ) ; // ignoring Var1 dependence of acceptance, knn and Var2 dependence of eff
         printf("  nB=%d Zll scale factors:  ee=%g, mm=%g\n", bbi, Zee_factor[bbi], Zmm_factor[bbi] ) ;
      } // bbi.




      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               float exp_0lep_ttwj(0.) ;
               if ( inputObservablesArePostTrigger ) {
                  exp_0lep_ttwj = sf_ttwj[mbi][hbi][bbi] * (N_1lep[mbi][hbi][bbi]/trigeff_1lep[mbi][hbi]) * initFit_R_ttwj_0lep_over_1lep ;
               } else {
                  exp_0lep_ttwj = sf_ttwj[mbi][hbi][bbi] * N_1lep[mbi][hbi][bbi] * initFit_R_ttwj_0lep_over_1lep ;
               }

               float exp_0lep_znn = 0.5 * ( N_Zee[mbi][hbi] * Zee_factor[bbi] + N_Zmm[mbi][hbi] * Zmm_factor[bbi] ) ;

               double qcdratio(0.) ;
               if ( qcdModelIndex == 2 ) {
                  qcdratio = initFit_qcd_0lepLDP_ratio_HT[hbi] ;
               } else if ( qcdModelIndex == 3 ) {
                  qcdratio = initFit_qcd_0lepLDP_ratio ;
               } else if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {
                  qcdratio = initFit_qcd_0lepLDP_ratio_HT[hbi] * initFit_SFqcd_met[mbi] * initFit_SFqcd_nb[bbi] ;
               }

               float exp_0lep_qcd(0.) ;
               if ( inputObservablesArePostTrigger ) {
               // exp_0lep_qcd =  sf_qcd[mbi][hbi][bbi] * ( N_ldp[mbi][hbi][bbi] / trigeff_0lep[mbi][hbi]
               //     - exp_0lep_ttwj* ttwj_mc_ldpover0lep_ratio[mbi][hbi][bbi]
               //     - exp_0lep_znn *  znn_mc_ldpover0lep_ratio[mbi][hbi][bbi] )
               //     * qcdratio ;
                  exp_0lep_qcd =  sf_qcd[mbi][hbi][bbi] * (
                  (1./trigeff_0lep[mbi][hbi])*(N_ldp[mbi][hbi][bbi]
                      - trigeff_1lep[mbi][hbi]*( exp_0lep_ttwj* ttwj_mc_ldpover0lep_ratio[mbi][hbi][bbi]
                                               + exp_0lep_znn *  znn_mc_ldpover0lep_ratio[mbi][hbi][bbi] ) ) )
                      * qcdratio ;
               } else {
                  exp_0lep_qcd =  sf_qcd[mbi][hbi][bbi] * ( N_ldp[mbi][hbi][bbi]
                      - exp_0lep_ttwj* ttwj_mc_ldpover0lep_ratio[mbi][hbi][bbi]
                      - exp_0lep_znn *  znn_mc_ldpover0lep_ratio[mbi][hbi][bbi] )
                      * qcdratio ;
               }

               if ( exp_0lep_qcd < 0. ) { exp_0lep_qcd = 0. ; }

               /// float exp_0lep = trigeff_0lep[mbi][hbi] * ( exp_0lep_ttwj + exp_0lep_qcd + exp_0lep_znn ) ;
               float exp_0lep = trigeff_0lep[mbi][hbi] * ( exp_0lep_qcd ) + trigeff_1lep[mbi][hbi] * ( exp_0lep_ttwj + exp_0lep_znn ) ;

               if ( useExpected0lep ) {
                  toy_mean_N_0lep[mbi][hbi][bbi] = exp_0lep ;
                  mcval_ttwj_0lep_3da[mbi][hbi][bbi] = exp_0lep_ttwj ;
                  mcval_qcd_0lep_3da [mbi][hbi][bbi] = exp_0lep_qcd ;
                  mcval_znn_0lep_3da [mbi][hbi][bbi] = exp_0lep_znn ;
               } else {
                  toy_mean_N_0lep[mbi][hbi][bbi] = N_0lep_input[mbi][hbi][bbi] ;
               }
               toy_mean_N_1lepSig[mbi][hbi][bbi] = N_1lepSig[mbi][hbi][bbi] ;
               toy_mean_N_1lep   [mbi][hbi][bbi] = N_1lep[mbi][hbi][bbi] ;
               toy_mean_N_ldp    [mbi][hbi][bbi] = N_ldp [mbi][hbi][bbi] ;


               printf("\n 0lep: m,h,b (%d,%d,%d): ttwj=%5.1f, qcd=%5.1f, Znn=%5.1f,  expected total=%5.1f\n",
                    mbi, hbi, bbi, exp_0lep_ttwj, exp_0lep_qcd, exp_0lep_znn, toy_mean_N_0lep[mbi][hbi][bbi] ) ;
               printf("   ttwj : trig_eff * sf * N1lep * Rttwj = %5.2f * %5.2f * %6.1f * %5.2f = %6.1f\n",
		      trigeff_1lep[mbi][hbi], sf_ttwj[mbi][hbi][bbi],  N_1lep[mbi][hbi][bbi], initFit_R_ttwj_0lep_over_1lep,
		      trigeff_1lep[mbi][hbi]*exp_0lep_ttwj ) ;
               printf("   qcd  : trig_eff * sf * (Nldp - ttwj * RMCttwj - Znn * RMCznn) * Rqcd = %5.2f * %5.2f * (%6.1f - %6.1f * %4.2f - %6.1f * %4.2f) * %5.3f = %6.1f\n",
		      trigeff_0lep[mbi][hbi], sf_qcd[mbi][hbi][bbi], N_ldp[mbi][hbi][bbi], exp_0lep_ttwj, ttwj_mc_ldpover0lep_ratio[mbi][hbi][bbi],
		      exp_0lep_znn, znn_mc_ldpover0lep_ratio[mbi][hbi][bbi], qcdratio,
		      trigeff_0lep[mbi][hbi]*exp_0lep_qcd ) ;
               printf("   Znn  : 0.5 * (Nee * Fee + Nmm * Fmm)   =   0.5 * (%5.1f * %5.2f + %5.1f * %5.2f)   =   %6.1f\n",
                    N_Zee[mbi][hbi], Zee_factor[bbi], N_Zmm[mbi][hbi], Zmm_factor[bbi], exp_0lep_znn ) ;

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

            int ArraySize = 3 + 2*4*(nBinsVar1*nBinsVar2*nBinsBjets) ;
            double ArrayContent[ArraySize] ;

            for (int i = 0; i < ArraySize; ++ i) {
              insusystr >> ArrayContent[i];
            }

            pointMgl  = ArrayContent[0] ;
            pointMlsp = ArrayContent[1] ;

            if ( !( fabs(pointMgl - mgl) < 0.1 && fabs(pointMlsp - mlsp) < 0.1 ) ) { continue ; }

            nGen = (int)ArrayContent[2] ;

            /// int nBins = nBinsVar1*nBinsVar2*nBinsBjets*3 ;

            float input0lepSusyTotal(0.) ;

            for (int i = 0 ; i < nBinsVar1 ; i++) {
              for (int j = 0 ; j < nBinsVar2 ; j++) {
                for (int k = 0 ; k < nBinsBjets ; k++) {


                  float n_0l_raw    =  ArrayContent[3 + i*(nBinsVar2*nBinsBjets) + j*(nBinsBjets) + k] ;
                  float n_1lSig_raw =  ArrayContent[3 + (nBinsVar1*nBinsVar2*nBinsBjets) + i*(nBinsVar2*nBinsBjets) + j*(nBinsBjets) + k] ;
                  float n_1l_raw    =  ArrayContent[3 + 2*(nBinsVar1*nBinsVar2*nBinsBjets) + i*(nBinsVar2*nBinsBjets) + j*(nBinsBjets) + k] ;
                  float n_ldp_raw   =  ArrayContent[3 + 3*(nBinsVar1*nBinsVar2*nBinsBjets) + i*(nBinsVar2*nBinsBjets) + j*(nBinsBjets) + k] ;

                  float n_0l_correction    =  1. ;
                  float n_1lSig_correction =  1. ;
                  float n_1l_correction    =  1. ;
                  float n_ldp_correction   =  1. ;

                  susy_0lep[i][j][k]    = n_0l_raw * n_0l_correction ;
                  susy_1lepSig[i][j][k] = n_1lSig_raw * n_1l_correction ;
                  susy_1lep[i][j][k]    = n_1l_raw * n_1l_correction ;
                  susy_ldp[i][j][k]     = n_ldp_raw * n_ldp_correction ;

                  input0lepSusyTotal += n_0l_raw * n_0l_correction ;

                  /* here's how to read in the errors, but not sure what to do with them yet
		  n_0l_error[i][j][k]    = ArrayContent[3 + 4*(nBinsVar1*nBinsVar2*nBinsBtag) + i*(nBinsVar2*nBinsBtag) + j*(nBinsBtag) + k] ;
		  n_1lSig_error[i][j][k] = ArrayContent[3 + 5*(nBinsVar1*nBinsVar2*nBinsBtag) + i*(nBinsVar2*nBinsBtag) + j*(nBinsBtag) + k] ;
		  n_1l_error[i][j][k]    = ArrayContent[3 + 6*(nBinsVar1*nBinsVar2*nBinsBtag) + i*(nBinsVar2*nBinsBtag) + j*(nBinsBtag) + k] ;
		  n_ldp_error[i][j][k]   = ArrayContent[3 + 7*(nBinsVar1*nBinsVar2*nBinsBtag) + i*(nBinsVar2*nBinsBtag) + j*(nBinsBtag) + k] ;
	    
		  if ( n_0l_raw[i][j][k] > 0.00001 ) {
		    n_0l_error[i][j][k] = 100 * n_0l_error[i][j][k] / n_0l_raw[i][j][k] ;
		  }
		  else n_0l_error[i][j][k] = 0.1 ;

		  if ( n_1lSig_raw[i][j][k] > 0.00001 ) {
		    n_1lSig_error[i][j][k] = 100 * n_1lSig_error[i][j][k] / n_1lSig_raw[i][j][k] ;
		  }
		  else n_1lSig_error[i][j][k] = 0.1 ;

		  if ( n_1l_raw[i][j][k] > 0.00001 ) {
		    n_1l_error[i][j][k] = 100 * n_1l_error[i][j][k] / n_1l_raw[i][j][k] ;
		  }
		  else n_1l_error[i][j][k] = 0.1 ;
	    
		  if ( n_ldp_raw[i][j][k] > 0.00001 ) {
		    n_ldp_error[i][j][k] = 100 * n_ldp_error[i][j][k] / n_ldp_raw[i][j][k] ;
		  }
		  else n_ldp_error[i][j][k] = 0.1 ;
                  */

                }
              }
            }


            printf("\n\n  pointMgl = %g , pointMlsp = %g,  total 0lep events = %g\n", pointMgl, pointMlsp, input0lepSusyTotal ) ;
            printf("   reset susy 0lep total to %g\n\n", nSusy0lep ) ;

            float newTotal(0.) ;
            for (int i = 0 ; i < nBinsVar1 ; i++) {
              for (int j = 0 ; j < nBinsVar2 ; j++) {
                for (int k = 0 ; k < nBinsBjets ; k++) {
                   susy_0lep[i][j][k]    = (nSusy0lep/input0lepSusyTotal) * susy_0lep[i][j][k] ;
                   susy_1lepSig[i][j][k] = (nSusy0lep/input0lepSusyTotal) * susy_1lepSig[i][j][k] ;
                   susy_1lep[i][j][k]    = (nSusy0lep/input0lepSusyTotal) * susy_1lep[i][j][k] ;
                   susy_ldp[i][j][k]     = (nSusy0lep/input0lepSusyTotal) * susy_ldp[i][j][k] ;
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


        for (int i = 0 ; i < nBinsVar1 ; i++) {
          for (int j = 0 ; j < nBinsVar2 ; j++) {
            for (int k = 0 ; k < nBinsBjets ; k++) {
               printf("  0lep (%d,%d,%d)    :  SM = %8.1f,  SUSY = %6.1f\n", i,j,k, toy_mean_N_0lep[i][j][k],    susy_0lep[i][j][k] ) ;
               printf("  1lepSig (%d,%d,%d) :  SM = %8.1f,  SUSY = %6.1f\n", i,j,k, toy_mean_N_1lepSig[i][j][k], susy_1lepSig[i][j][k] ) ;
               printf("  1lep (%d,%d,%d)    :  SM = %8.1f,  SUSY = %6.1f\n", i,j,k, toy_mean_N_1lep[i][j][k],    susy_1lep[i][j][k] ) ;
               printf("  ldp  (%d,%d,%d)    :  SM = %8.1f,  SUSY = %6.1f\n", i,j,k, toy_mean_N_ldp [i][j][k],    susy_ldp [i][j][k] ) ;
            // toy_mean_N_0lep[i][j][k]    += susy_0lep[i][j][k]  ;
            // toy_mean_N_1lepSig[i][j][k] += susy_1lepSig[i][j][k]  ;
            // toy_mean_N_1lep[i][j][k]    += susy_1lep[i][j][k]  ;
            // toy_mean_N_ldp [i][j][k]    += susy_ldp [i][j][k]   ;
               toy_mean_N_0lep[i][j][k]    += trigeff_1lep[i][j] * susy_0lep[i][j][k]  ;
               toy_mean_N_1lepSig[i][j][k] += trigeff_1lep[i][j] * susy_1lepSig[i][j][k]  ;
               toy_mean_N_1lep[i][j][k]    += trigeff_1lep[i][j] * susy_1lep[i][j][k]  ;
               toy_mean_N_ldp [i][j][k]    += trigeff_1lep[i][j] * susy_ldp [i][j][k]   ;
               printf("\n") ;
            }
               printf(" ---- \n") ;
          }
               printf(" ======== \n") ;
        }





      return true ;

   } // addSusyToExpectedObs.


   //==================================================================================

   RooDataSet* genToyData( bool noFluctuations ) {

      RooArgSet obsList ;

      char oname[1000] ;

      int nsel(4) ;
      char selname[4][100] = { "0lep", "1lepSig", "1lep", "ldp" } ;

      for ( int si=0; si<nsel; si++ ) {
         for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
            for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
               for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

                  sprintf( oname, "N_%s_M%d_H%d_%db", selname[si], mbi+1, hbi+1, bbi+1 ) ;
                  RooRealVar* rrv = workspace -> var( oname ) ;
                  if ( rrv == 0x0 ) {
                     printf("\n\n *** can't find variable with name %s\n", oname ) ; cout << flush ;
                     if ( si == 0 ) Nobs_0lep[mbi][hbi][bbi] = 0 ;
                     if ( si == 1 ) Nobs_1lepSig[mbi][hbi][bbi] = 0 ;
                     if ( si == 2 ) Nobs_1lep[mbi][hbi][bbi] = 0 ;
                     if ( si == 3 ) Nobs_ldp [mbi][hbi][bbi] = 0 ;
                     continue ;
                  }

                  if ( si == 0 ) {
                     int nobs ;
                     if ( noFluctuations ) {
                        nobs = TMath::Nint( toy_mean_N_0lep[mbi][hbi][bbi] ) ;
                     } else {
                        nobs = tran->Poisson( toy_mean_N_0lep[mbi][hbi][bbi] ) ;
                     }

		     printf("\nDebug: toy_mean_N_0lep[%d][%d][%d] = %f", mbi+1, hbi+1, bbi+1, toy_mean_N_0lep[mbi][hbi][bbi] ) ; 

                     rrv -> setVal( nobs  ) ;
                     Nobs_0lep[mbi][hbi][bbi] = nobs ;
                  } else if ( si == 1 ) {
                     int nobs ;
                     if ( noFluctuations ) {
                        nobs = TMath::Nint( toy_mean_N_1lepSig[mbi][hbi][bbi] ) ;
                     } else {
                        nobs = tran->Poisson( toy_mean_N_1lepSig[mbi][hbi][bbi] ) ;
                     }

		     printf("\nDebug: toy_mean_N_1lepSig[%d][%d][%d] = %f", mbi+1, hbi+1, bbi+1, toy_mean_N_1lepSig[mbi][hbi][bbi] ) ; 

                     rrv -> setVal( nobs ) ;
                     Nobs_1lepSig[mbi][hbi][bbi] = nobs ;                  
		  } else if ( si == 2 ) {
                     int nobs ;
                     if ( noFluctuations ) {
                        nobs = TMath::Nint( toy_mean_N_1lep[mbi][hbi][bbi] ) ;
                     } else {
                        nobs = tran->Poisson( toy_mean_N_1lep[mbi][hbi][bbi] ) ;
                     }

		     printf("\nDebug: toy_mean_N_1lep[%d][%d][%d] = %f", mbi+1, hbi+1, bbi+1, toy_mean_N_1lep[mbi][hbi][bbi] ) ; 

                     rrv -> setVal( nobs ) ;
                     Nobs_1lep[mbi][hbi][bbi] = nobs ;
                  } else if ( si == 3 ) {
                     int nobs ;
                     if ( noFluctuations ) {
                        nobs = TMath::Nint( toy_mean_N_ldp [mbi][hbi][bbi] ) ;
                     } else {
                        nobs = tran->Poisson( toy_mean_N_ldp [mbi][hbi][bbi] ) ;
                     }

		     printf("\nDebug: toy_mean_N_ldp[%d][%d][%d] = %f", mbi+1, hbi+1, bbi+1, toy_mean_N_ldp[mbi][hbi][bbi] ) ; 

                     rrv -> setVal( nobs ) ;
                     Nobs_ldp[mbi][hbi][bbi] = nobs ;
                  }

                  obsList.add( *rrv ) ;

               } // bbi.
            } // hbi.
         } // mbi.
      }

      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
            sprintf( oname, "N_Zee_M%d_H%d", mbi+1, hbi+1 ) ;
            RooRealVar* rrv = workspace -> var( oname ) ;
            if ( rrv == 0x0 ) {
               printf("\n\n *** can't find variable with name %s\n", oname ) ;
               Nobs_Zee[mbi][hbi] = 0 ;
               continue ;
            }
            int nobs ;
            if ( noFluctuations ) {
               nobs = TMath::Nint( toy_mean_N_Zee[mbi][hbi] ) ;
            } else {
               nobs = tran->Poisson( toy_mean_N_Zee[mbi][hbi] ) ;
            }
            rrv -> setVal( nobs ) ;
            Nobs_Zee[mbi][hbi] = nobs ;
            obsList.add( *rrv ) ;
         } // hbi.
      } // mbi.

      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
            sprintf( oname, "N_Zmm_M%d_H%d", mbi+1, hbi+1 ) ;
            RooRealVar* rrv = workspace -> var( oname ) ;
            if ( rrv == 0x0 ) {
               printf("\n\n *** can't find variable with name %s\n", oname ) ;
               Nobs_Zmm[mbi][hbi] = 0 ;
               continue ;
            }
            int nobs ;
            if ( noFluctuations ) {
               nobs = TMath::Nint( toy_mean_N_Zmm[mbi][hbi] ) ;
            } else {
               nobs = tran->Poisson( toy_mean_N_Zmm[mbi][hbi] ) ;
            }
            rrv -> setVal( nobs ) ;
            Nobs_Zmm[mbi][hbi] = nobs ;
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

      fit_chi2_overall = 0. ;
      fit_chi2_obs = 0. ;
      fit_chi2_np = 0. ;
      fit_chi2_prob = 0. ;

      printf("\n\n") ;
      printf(" toy %4d : Fit nSusy 0lep : %4.1f +/- %4.1f (%4.1f, %4.1f)\n", ti,
           (rrv_susy_poi->getVal()),
           (rrv_susy_poi->getError()),
           (rrv_susy_poi->getErrorHi()),
           (rrv_susy_poi->getErrorLo())
          ) ;

      susy_poi_atMinNll = rrv_susy_poi->getVal() ;
      susy_poi_plusErr  = rrv_susy_poi->getErrorHi() ;

      fit_susy_0lep     =  (rrv_susy_poi->getVal()) ;
      fit_susy_0lep_err =  (rrv_susy_poi->getError()) ;
      fit_susy_0lep_err_low =  (rrv_susy_poi->getErrorLo()) ;
      fit_susy_0lep_err_high =  (rrv_susy_poi->getErrorHi()) ;


      if ( fit_susy_0lep > true_susy_0lep ) {
         fit_susy_0lep_err_forpull = fit_susy_0lep_err_low ;
      } else {
         fit_susy_0lep_err_forpull = fit_susy_0lep_err_high ;
      }

      if ( qcdModelIndex == 2 || qcdModelIndex == 4 || qcdModelIndex == 5 ) {
         if ( nBinsVar2 >= 1 ) {
            fit_qcd_0lepLDP_ratio_H1     = rrv_qcd_0lepLDP_ratio_H1 -> getVal() ;
            fit_qcd_0lepLDP_ratio_H1_err = rrv_qcd_0lepLDP_ratio_H1 -> getError() ;
         }
         if ( nBinsVar2 >= 2 ) {
            fit_qcd_0lepLDP_ratio_H2     = rrv_qcd_0lepLDP_ratio_H2 -> getVal() ;
            fit_qcd_0lepLDP_ratio_H2_err = rrv_qcd_0lepLDP_ratio_H2 -> getError() ;
         }
         if ( nBinsVar2 >= 3 ) {
            fit_qcd_0lepLDP_ratio_H3     = rrv_qcd_0lepLDP_ratio_H3 -> getVal() ;
            fit_qcd_0lepLDP_ratio_H3_err = rrv_qcd_0lepLDP_ratio_H3 -> getError() ;
         }
         if ( nBinsVar2 >= 4 ) {
            fit_qcd_0lepLDP_ratio_H4     = rrv_qcd_0lepLDP_ratio_H4 -> getVal() ;
            fit_qcd_0lepLDP_ratio_H4_err = rrv_qcd_0lepLDP_ratio_H4 -> getError() ;
         }
      } 
      if ( qcdModelIndex == 3 ) {
         fit_qcd_0lepLDP_ratio     = rrv_qcd_0lepLDP_ratio -> getVal() ;
         fit_qcd_0lepLDP_ratio_err = rrv_qcd_0lepLDP_ratio -> getError() ;
      }
      if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {
         if ( nBinsVar1 >=2 ) {
            fit_SFqcd_met2     = rrv_SFqcd_met2 -> getVal() ;
            fit_SFqcd_met2_err = rrv_SFqcd_met2 -> getError() ;
         }
         if ( nBinsBjets >=2 ) {
            fit_SFqcd_nb2     = rrv_SFqcd_nb2 -> getVal() ;
            fit_SFqcd_nb2_err = rrv_SFqcd_nb2 -> getError() ;
         }
      }

      /*
      if ( qcdModelIndex == 4 ) {
         if ( nBinsVar1 >=3 ) {
            fit_SFqcd_met3     = rrv_SFqcd_met3 -> getVal() ;
            fit_SFqcd_met3_err = rrv_SFqcd_met3 -> getError() ;
         }
         if ( nBinsVar1 >=4 ) {
            fit_SFqcd_met4     = rrv_SFqcd_met4 -> getVal() ;
            fit_SFqcd_met4_err = rrv_SFqcd_met4 -> getError() ;
         }
         if ( nBinsBjets >=3 ) {
            fit_SFqcd_nb3     = rrv_SFqcd_nb3 -> getVal() ;
            fit_SFqcd_nb3_err = rrv_SFqcd_nb3 -> getError() ;
         }
         if ( nBinsBjets >=4 ) {
            fit_SFqcd_nb4     = rrv_SFqcd_nb4 -> getVal() ;
            fit_SFqcd_nb4_err = rrv_SFqcd_nb4 -> getError() ;
         }
      }
      */

      fit_ttwj_0lep1lep_ratio     = rrv_ttwj_0lep1lep_ratio -> getVal() ;
      fit_ttwj_0lep1lep_ratio_err = rrv_ttwj_0lep1lep_ratio -> getError() ;

      fit_ttwj_1lepSig1lep_ratio     = rrv_ttwj_1lepSig1lep_ratio -> getVal() ;
      fit_ttwj_1lepSig1lep_ratio_err = rrv_ttwj_1lepSig1lep_ratio -> getError() ;

      fit_susy_0lep_wsfs = 0. ;
      double fit_ttwj_0lep = 0. ;
      double fit_qcd__0lep = 0. ;
      double fit_znn__0lep = 0. ;

      double fit_susy_0lep_1b = 0. ;
      double fit_ttwj_0lep_1b = 0. ;
      double fit_qcd__0lep_1b = 0. ;
      double fit_znn__0lep_1b = 0. ;

      double fit_susy_0lep_2b = 0. ;
      double fit_ttwj_0lep_2b = 0. ;
      double fit_qcd__0lep_2b = 0. ;
      double fit_znn__0lep_2b = 0. ;

      double fit_susy_0lep_3b = 0. ;
      double fit_ttwj_0lep_3b = 0. ;
      double fit_qcd__0lep_3b = 0. ;
      double fit_znn__0lep_3b = 0. ;

      double fit_susy_0lep_nb_hm1 = 0. ;
      double fit_ttwj_0lep_nb_hm1 = 0. ;
      double fit_qcd__0lep_nb_hm1 = 0. ;
      double fit_znn__0lep_nb_hm1 = 0. ;

      double fit_susy_0lep_1b_hm1 = 0. ;
      double fit_ttwj_0lep_1b_hm1 = 0. ;
      double fit_qcd__0lep_1b_hm1 = 0. ;
      double fit_znn__0lep_1b_hm1 = 0. ;

      double fit_susy_0lep_2b_hm1 = 0. ;
      double fit_ttwj_0lep_2b_hm1 = 0. ;
      double fit_qcd__0lep_2b_hm1 = 0. ;
      double fit_znn__0lep_2b_hm1 = 0. ;

      double fit_susy_0lep_3b_hm1 = 0. ;
      double fit_ttwj_0lep_3b_hm1 = 0. ;
      double fit_qcd__0lep_3b_hm1 = 0. ;
      double fit_znn__0lep_3b_hm1 = 0. ;

      double fit_susy_0lep_nb_hh1 = 0. ;
      double fit_ttwj_0lep_nb_hh1 = 0. ;
      double fit_qcd__0lep_nb_hh1 = 0. ;
      double fit_znn__0lep_nb_hh1 = 0. ;

      double fit_susy_0lep_1b_hh1 = 0. ;
      double fit_ttwj_0lep_1b_hh1 = 0. ;
      double fit_qcd__0lep_1b_hh1 = 0. ;
      double fit_znn__0lep_1b_hh1 = 0. ;

      double fit_susy_0lep_2b_hh1 = 0. ;
      double fit_ttwj_0lep_2b_hh1 = 0. ;
      double fit_qcd__0lep_2b_hh1 = 0. ;
      double fit_znn__0lep_2b_hh1 = 0. ;

      double fit_susy_0lep_3b_hh1 = 0. ;
      double fit_ttwj_0lep_3b_hh1 = 0. ;
      double fit_qcd__0lep_3b_hh1 = 0. ;
      double fit_znn__0lep_3b_hh1 = 0. ;

      double fit_susy_0lep_nb_hm2_hh2 = 0. ;
      double fit_ttwj_0lep_nb_hm2_hh2 = 0. ;
      double fit_qcd__0lep_nb_hm2_hh2 = 0. ;
      double fit_znn__0lep_nb_hm2_hh2 = 0. ;

      double fit_susy_0lep_1b_hm2_hh2 = 0. ;
      double fit_ttwj_0lep_1b_hm2_hh2 = 0. ;
      double fit_qcd__0lep_1b_hm2_hh2 = 0. ;
      double fit_znn__0lep_1b_hm2_hh2 = 0. ;

      double fit_susy_0lep_2b_hm2_hh2 = 0. ;
      double fit_ttwj_0lep_2b_hm2_hh2 = 0. ;
      double fit_qcd__0lep_2b_hm2_hh2 = 0. ;
      double fit_znn__0lep_2b_hm2_hh2 = 0. ;

      double fit_susy_0lep_3b_hm2_hh2 = 0. ;
      double fit_ttwj_0lep_3b_hm2_hh2 = 0. ;
      double fit_qcd__0lep_3b_hm2_hh2 = 0. ;
      double fit_znn__0lep_3b_hm2_hh2 = 0. ;

      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
            char vname[1000] ;
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               if ( ignoreBin[mbi][hbi] ) continue ;

               RooAbsReal* rar ;
               RooAbsReal* effsf  ;
               RooAbsReal* btagsf  ;


               double nsusy(0.) ;
               sprintf( vname, "mu_susy_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               rar = workspace -> function( vname ) ;
               if ( rar == 0x0 ) { printf("\n\n *** processToyResult: missing var %s\n\n", vname ) ; }
               sprintf( vname, "btageff_sf_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               btagsf = (RooAbsReal*) workspace -> obj( vname ) ;
               if ( btagsf == 0x0 ) { printf("\n\n *** processToyResult: missing var %s\n\n", vname ) ; }
               sprintf( vname, "eff_sf_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               effsf = (RooAbsReal*) workspace -> obj( vname ) ;
               if ( effsf == 0x0 ) { printf("\n\n *** processToyResult: missing var %s\n\n", vname ) ; }
               if ( rar != 0x0 && btagsf != 0x0 && effsf != 0x0 ) {
                  nsusy = ( btagsf -> getVal() ) * ( effsf -> getVal() ) * ( rar -> getVal() ) ;
               }

               fit_susy_0lep_wsfs += nsusy ;

               fit_susy_0lep_3da[mbi][hbi][bbi] = nsusy ;

               if ( bbi==0 ) { fit_susy_0lep_1b += nsusy ; }
               if ( bbi==1 ) { fit_susy_0lep_2b += nsusy ; }
               if ( bbi==2 ) { fit_susy_0lep_3b += nsusy ; }
               if ( mbi==(nBinsVar1-1) ) {
                  fit_susy_0lep_nb_hm1 += nsusy ;
                  if ( bbi==0 ) { fit_susy_0lep_1b_hm1 += nsusy ; }
                  if ( bbi==1 ) { fit_susy_0lep_2b_hm1 += nsusy ; }
                  if ( bbi==2 ) { fit_susy_0lep_3b_hm1 += nsusy ; }
               }
               if ( hbi==(nBinsVar2-1) ) {
                  fit_susy_0lep_nb_hh1 += nsusy ;
                  if ( bbi==0 ) { fit_susy_0lep_1b_hh1 += nsusy ; }
                  if ( bbi==1 ) { fit_susy_0lep_2b_hh1 += nsusy ; }
                  if ( bbi==2 ) { fit_susy_0lep_3b_hh1 += nsusy ; }
               }
               if ( mbi>=(nBinsVar1-2) && hbi>=(nBinsVar2-2) ) {
                  fit_susy_0lep_nb_hm2_hh2 += nsusy ;
                  if ( bbi==0 ) { fit_susy_0lep_1b_hm2_hh2 += nsusy ; }
                  if ( bbi==1 ) { fit_susy_0lep_2b_hm2_hh2 += nsusy ; }
                  if ( bbi==2 ) { fit_susy_0lep_3b_hm2_hh2 += nsusy ; }
               }



               double nttwj(0.) ;
               sprintf( vname, "mu_ttwj_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               rar = (RooAbsReal*) workspace -> obj( vname ) ;
               if ( rar == 0x0 ) {
                  printf("\n\n *** missing var %s\n\n", vname ) ;
               } else {
                  nttwj = rar -> getVal() ;
               }
               fit_ttwj_0lep_3da[mbi][hbi][bbi] = nttwj ;
               fit_ttwj_0lep += nttwj  ;
               if ( bbi==0 ) { fit_ttwj_0lep_1b +=  nttwj  ; }
               if ( bbi==1 ) { fit_ttwj_0lep_2b +=  nttwj  ; }
               if ( bbi==2 ) { fit_ttwj_0lep_3b +=  nttwj  ; }
               if ( mbi==(nBinsVar1-1) ) {
                  fit_ttwj_0lep_nb_hm1 += nttwj ;
                  if ( bbi==0 ) { fit_ttwj_0lep_1b_hm1 += nttwj ; }
                  if ( bbi==1 ) { fit_ttwj_0lep_2b_hm1 += nttwj ; }
                  if ( bbi==2 ) { fit_ttwj_0lep_3b_hm1 += nttwj ; }
               }
               if ( hbi==(nBinsVar2-1) ) {
                  fit_ttwj_0lep_nb_hh1 += nttwj ;
                  if ( bbi==0 ) { fit_ttwj_0lep_1b_hh1 += nttwj ; }
                  if ( bbi==1 ) { fit_ttwj_0lep_2b_hh1 += nttwj ; }
                  if ( bbi==2 ) { fit_ttwj_0lep_3b_hh1 += nttwj ; }
               }
               if ( mbi>=(nBinsVar1-2) && hbi>=(nBinsVar2-2) ) {
                  fit_ttwj_0lep_nb_hm2_hh2 += nttwj ;
                  if ( bbi==0 ) { fit_ttwj_0lep_1b_hm2_hh2 += nttwj ; }
                  if ( bbi==1 ) { fit_ttwj_0lep_2b_hm2_hh2 += nttwj ; }
                  if ( bbi==2 ) { fit_ttwj_0lep_3b_hm2_hh2 += nttwj ; }
               }

               double nqcd(0.) ;
               sprintf( vname, "mu_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               rar = workspace -> function( vname ) ;
               if ( rar == 0x0 ) {
                  printf("\n\n *** missing var %s\n\n", vname ) ;
               } else {
                  nqcd = rar -> getVal() ;
               }
               fit_qcd_0lep_3da[mbi][hbi][bbi] = nqcd ;
               fit_qcd__0lep += nqcd  ;
               if ( bbi==0 ) { fit_qcd__0lep_1b +=  nqcd   ; }
               if ( bbi==1 ) { fit_qcd__0lep_2b +=  nqcd   ; }
               if ( bbi==2 ) { fit_qcd__0lep_3b +=  nqcd   ; }
               if ( mbi==(nBinsVar1-1) ) {
                  fit_qcd__0lep_nb_hm1 += nqcd ;
                  if ( bbi==0 ) { fit_qcd__0lep_1b_hm1 += nqcd ; }
                  if ( bbi==1 ) { fit_qcd__0lep_2b_hm1 += nqcd ; }
                  if ( bbi==2 ) { fit_qcd__0lep_3b_hm1 += nqcd ; }
               }
               if ( hbi==(nBinsVar2-1) ) {
                  fit_qcd__0lep_nb_hh1 += nqcd ;
                  if ( bbi==0 ) { fit_qcd__0lep_1b_hh1 += nqcd ; }
                  if ( bbi==1 ) { fit_qcd__0lep_2b_hh1 += nqcd ; }
                  if ( bbi==2 ) { fit_qcd__0lep_3b_hh1 += nqcd ; }
               }
               if ( mbi>=(nBinsVar1-2) && hbi>=(nBinsVar2-2) ) {
                  fit_qcd__0lep_nb_hm2_hh2 += nqcd ;
                  if ( bbi==0 ) { fit_qcd__0lep_1b_hm2_hh2 += nqcd ; }
                  if ( bbi==1 ) { fit_qcd__0lep_2b_hm2_hh2 += nqcd ; }
                  if ( bbi==2 ) { fit_qcd__0lep_3b_hm2_hh2 += nqcd ; }
               }



               double nznn(0.) ;
               sprintf( vname, "mu_znn_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               rar = (RooAbsReal*) workspace -> obj( vname ) ;
               if ( rar == 0x0 ) {
                  printf("\n\n *** missing var %s\n\n", vname ) ;
               } else {
                  nznn = rar -> getVal() ;
               }
               fit_znn_0lep_3da[mbi][hbi][bbi] = nznn ;
               fit_znn__0lep += nznn  ;
               if ( bbi==0 ) { fit_znn__0lep_1b +=  nznn   ; }
               if ( bbi==1 ) { fit_znn__0lep_2b +=  nznn   ; }
               if ( bbi==2 ) { fit_znn__0lep_3b +=  nznn   ; }
               if ( mbi==(nBinsVar1-1) ) {
                  fit_znn__0lep_nb_hm1 += nznn ;
                  if ( bbi==0 ) { fit_znn__0lep_1b_hm1 += nznn ; }
                  if ( bbi==1 ) { fit_znn__0lep_2b_hm1 += nznn ; }
                  if ( bbi==2 ) { fit_znn__0lep_3b_hm1 += nznn ; }
               }
               if ( hbi==(nBinsVar2-1) ) {
                  fit_znn__0lep_nb_hh1 += nznn ;
                  if ( bbi==0 ) { fit_znn__0lep_1b_hh1 += nznn ; }
                  if ( bbi==1 ) { fit_znn__0lep_2b_hh1 += nznn ; }
                  if ( bbi==2 ) { fit_znn__0lep_3b_hh1 += nznn ; }
               }
               if ( mbi>=(nBinsVar1-2) && hbi>=(nBinsVar2-2) ) {
                  fit_znn__0lep_nb_hm2_hh2 += nznn ;
                  if ( bbi==0 ) { fit_znn__0lep_1b_hm2_hh2 += nznn ; }
                  if ( bbi==1 ) { fit_znn__0lep_2b_hm2_hh2 += nznn ; }
                  if ( bbi==2 ) { fit_znn__0lep_3b_hm2_hh2 += nznn ; }
               }

               sprintf( vname, "sf_ttwj_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               rar = (RooAbsReal*) workspace -> obj( vname ) ;
               if ( rar == 0x0 ) {
                  printf("\n\n *** missing var %s\n\n", vname ) ;
                  fit_sf_ttwj_0lep_3da[mbi][hbi][bbi] = 1. ;
               } else {
                  fit_sf_ttwj_0lep_3da[mbi][hbi][bbi] = rar -> getVal() ;
               }

               sprintf( vname, "sf_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               rar = (RooAbsReal*) workspace -> obj( vname ) ;
               if ( rar == 0x0 ) {
                  printf("\n\n *** missing var %s\n\n", vname ) ;
                  fit_sf_qcd_0lep_3da[mbi][hbi][bbi] = 1. ;
               } else {
                  fit_sf_qcd_0lep_3da[mbi][hbi][bbi] = rar -> getVal() ;
               }





               //++++ compute chi2 (below here).

               char obsname[1000] ;
               char modelname[1000] ;

               //-- 0lep observable
               if ( !blind0lepBin[mbi][hbi][bbi] ) {
                  sprintf( obsname,   "N_0lep_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  sprintf( modelname, "n_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  fit_chi2_obs += getChi2Obs( obsname, modelname ) ;
               }

               //-- 1lepSig observable
               sprintf( obsname,   "N_1lepSig_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               sprintf( modelname, "n_slSig_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               fit_chi2_obs += getChi2Obs( obsname, modelname ) ;

               //-- 1lep observable
               sprintf( obsname,   "N_1lep_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               sprintf( modelname, "n_sl_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               fit_chi2_obs += getChi2Obs( obsname, modelname ) ;

               //-- LDP observable
               sprintf( obsname,   "N_ldp_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               sprintf( modelname, "n_ldp_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               fit_chi2_obs += getChi2Obs( obsname, modelname ) ;




               //-- Zll observables
               if ( bbi == 0 ) {

                  sprintf( obsname,   "N_Zee_M%d_H%d", mbi+1, hbi+1 ) ;
                  sprintf( modelname, "n_ee_M%d_H%d", mbi+1, hbi+1 ) ;
                  fit_chi2_obs += getChi2Obs( obsname, modelname ) ;

                  sprintf( obsname,   "N_Zmm_M%d_H%d", mbi+1, hbi+1 ) ;
                  sprintf( modelname, "n_mm_M%d_H%d", mbi+1, hbi+1 ) ;
                  fit_chi2_obs += getChi2Obs( obsname, modelname ) ;

               }



               if ( !blind0lepBin[mbi][hbi][bbi] ) {


                 //-- ttwj and QCD SFs

                  sprintf( vname, "sf_ttwj_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  fit_chi2_np += getChi2GausNP( vname ) ;

                  sprintf( vname, "sf_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  fit_chi2_np += getChi2GausNP( vname ) ;



                 //-- Signal efficiency uncertainties (bin-by-bin MC stat err).

                  sprintf( vname, "eff_sf_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  fit_chi2_np += getChi2GausNP( vname ) ;

                  sprintf( vname, "eff_sf_ldp_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  fit_chi2_np += getChi2GausNP( vname ) ;

                  sprintf( vname, "eff_sf_slSig_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  fit_chi2_np += getChi2GausNP( vname ) ;

                  sprintf( vname, "eff_sf_sl_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  fit_chi2_np += getChi2GausNP( vname ) ;



               }




            } // bbi.


            //-- trigger efficiencies.

            sprintf( vname, "trigeff_M%d_H%d", mbi+1, hbi+1 ) ;
            fit_chi2_np += getChi2BetaNP( vname ) ;

            sprintf( vname, "trigeff_sl_M%d_H%d", mbi+1, hbi+1 ) ;
            fit_chi2_np += getChi2BetaNP( vname ) ;


         } // hbi.
      } // mbi.

      //-- more nuisance parameter chi2 contributions.
      char npname[1000] ;

      ////// sprintf( npname, "eff_sf" ) ;             /// obsolete
      ////// fit_chi2_np += getChi2GausNP( npname ) ;  /// obsolete

      sprintf( npname, "rar_vv_sf" ) ;
      fit_chi2_np += getChi2GausNP( npname ) ;


      sprintf( npname, "SFqcd_met3" ) ;
      fit_chi2_np += getChi2GausNP( npname ) ;

      sprintf( npname, "SFqcd_met4" ) ;
      fit_chi2_np += getChi2GausNP( npname ) ;

      sprintf( npname, "SFqcd_nb3" ) ;
      fit_chi2_np += getChi2GausNP( npname ) ;


      sprintf( npname, "btageff_sf" ) ;
      fit_chi2_np += getChi2GausNP( npname ) ;

      sprintf( npname, "sf_mc" ) ;
      fit_chi2_np += getChi2GausNP( npname ) ;

      sprintf( npname, "sf_ll" ) ;
      fit_chi2_np += getChi2GausNP( npname ) ;

      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         sprintf( npname, "knn_1b_M%d", mbi+1 ) ;
         fit_chi2_np += getChi2GausNP( npname ) ;
      } // mbi

      sprintf( npname, "knn_2b" ) ;
      fit_chi2_np += getChi2GausNP( npname ) ;

      sprintf( npname, "knn_3b" ) ;
      fit_chi2_np += getChi2GausNP( npname ) ;

      sprintf( npname, "eff_Zee" ) ;
      fit_chi2_np += getChi2BetaNP( npname ) ;
      
      sprintf( npname, "eff_Zmm" ) ;
      fit_chi2_np += getChi2BetaNP( npname ) ;
      
      sprintf( npname, "pur_Zee" ) ;
      fit_chi2_np += getChi2BetaNP( npname ) ;

      sprintf( npname, "pur_Zmm" ) ;
      fit_chi2_np += getChi2BetaNP( npname ) ;

      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {

         sprintf( npname, "acc_Zee_M%d", mbi+1 ) ;
         fit_chi2_np += getChi2BetaNP( npname ) ;

         sprintf( npname, "acc_Zmm_M%d", mbi+1 ) ;
         fit_chi2_np += getChi2BetaNP( npname ) ;

      } // mbi.

      fit_chi2_overall = fit_chi2_obs + fit_chi2_np ;

      //-- 3 is ZL,SL,LDP, 2 is ee,mm.
      int nobs = 3 * ( ( nBinsVar1*nBinsVar2 - num_ignoreBin ) * nBinsBjets ) + 2 * ( nBinsVar1*nBinsVar2 - num_ignoreBin ) - num_blind0lepBin ;
      //-- 2 is ttwj,QCD, 8 at end is ttwj ZL/SL + 6 QCD model pars + signal yield.
      int nfloat = 2 * ( ( nBinsVar1*nBinsVar2 - num_ignoreBin ) * nBinsBjets) + (nBinsVar1*nBinsVar2 - num_ignoreBin ) + 8 ;

      int ndof = nobs - nfloat ;
      printf(" Nobs = %d, Nfloat = %d\n", nobs, nfloat ) ;

      fit_chi2_prob = TMath::Prob( fit_chi2_overall, ndof ) ;
      printf(" toy %4d : chi2 / ndof = %8.2f / %d,  prob = %7.4f\n", ti, fit_chi2_overall, ndof, fit_chi2_prob ) ;

      //-- end nuisance parameter chi2 contributions.


      fit_covqual_susyfloat = rfr -> covQual() ;

      printf(" toy %4d : Fit total 0lep ttwj : %6.1f\n", ti, fit_ttwj_0lep ) ;
      printf(" toy %4d : Fit total 0lep qcd  : %6.1f\n", ti, fit_qcd__0lep ) ;
      printf(" toy %4d : Fit total 0lep znn  : %6.1f\n", ti, fit_znn__0lep ) ;
      printf(" toy %4d : Fit covariance matrix quality: %d\n", ti, fit_covqual_susyfloat ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 1b 0lep susy : %6.1f  \n", ti, fit_susy_0lep_1b ) ;
      printf(" toy %4d : Fit 1b 0lep ttwj : %6.1f  \n", ti, fit_ttwj_0lep_1b ) ;
      printf(" toy %4d : Fit 1b 0lep qcd  : %6.1f  \n", ti, fit_qcd__0lep_1b ) ;
      printf(" toy %4d : Fit 1b 0lep znn  : %6.1f  \n", ti, fit_znn__0lep_1b ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 2b 0lep susy : %6.1f  \n", ti, fit_susy_0lep_2b ) ;
      printf(" toy %4d : Fit 2b 0lep ttwj : %6.1f  \n", ti, fit_ttwj_0lep_2b ) ;
      printf(" toy %4d : Fit 2b 0lep qcd  : %6.1f  \n", ti, fit_qcd__0lep_2b ) ;
      printf(" toy %4d : Fit 2b 0lep znn  : %6.1f  \n", ti, fit_znn__0lep_2b ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 3b 0lep susy : %6.1f  \n", ti, fit_susy_0lep_3b ) ;
      printf(" toy %4d : Fit 3b 0lep ttwj : %6.1f  \n", ti, fit_ttwj_0lep_3b ) ;
      printf(" toy %4d : Fit 3b 0lep qcd  : %6.1f  \n", ti, fit_qcd__0lep_3b ) ;
      printf(" toy %4d : Fit 3b 0lep znn  : %6.1f  \n", ti, fit_znn__0lep_3b ) ;



      printf("\n\n") ;
      printf(" toy %4d : Fit nb 0lep susy, highest Var2 bin  : %6.1f  \n", ti, fit_susy_0lep_nb_hh1 ) ;
      printf(" toy %4d : Fit nb 0lep ttwj, highest Var2 bin  : %6.1f  \n", ti, fit_ttwj_0lep_nb_hh1 ) ;
      printf(" toy %4d : Fit nb 0lep qcd,  highest Var2 bin  : %6.1f  \n", ti, fit_qcd__0lep_nb_hh1 ) ;
      printf(" toy %4d : Fit nb 0lep znn,  highest Var2 bin  : %6.1f  \n", ti, fit_znn__0lep_nb_hh1 ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 1b 0lep susy, highest Var2 bin  : %6.1f  \n", ti, fit_susy_0lep_1b_hh1 ) ;
      printf(" toy %4d : Fit 1b 0lep ttwj, highest Var2 bin  : %6.1f  \n", ti, fit_ttwj_0lep_1b_hh1 ) ;
      printf(" toy %4d : Fit 1b 0lep qcd,  highest Var2 bin  : %6.1f  \n", ti, fit_qcd__0lep_1b_hh1 ) ;
      printf(" toy %4d : Fit 1b 0lep znn,  highest Var2 bin  : %6.1f  \n", ti, fit_znn__0lep_1b_hh1 ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 2b 0lep susy, highest Var2 bin  : %6.1f  \n", ti, fit_susy_0lep_2b_hh1 ) ;
      printf(" toy %4d : Fit 2b 0lep ttwj, highest Var2 bin  : %6.1f  \n", ti, fit_ttwj_0lep_2b_hh1 ) ;
      printf(" toy %4d : Fit 2b 0lep qcd,  highest Var2 bin  : %6.1f  \n", ti, fit_qcd__0lep_2b_hh1 ) ;
      printf(" toy %4d : Fit 2b 0lep znn,  highest Var2 bin  : %6.1f  \n", ti, fit_znn__0lep_2b_hh1 ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 3b 0lep susy, highest Var2 bin  : %6.1f  \n", ti, fit_susy_0lep_3b_hh1 ) ;
      printf(" toy %4d : Fit 3b 0lep ttwj, highest Var2 bin  : %6.1f  \n", ti, fit_ttwj_0lep_3b_hh1 ) ;
      printf(" toy %4d : Fit 3b 0lep qcd,  highest Var2 bin  : %6.1f  \n", ti, fit_qcd__0lep_3b_hh1 ) ;
      printf(" toy %4d : Fit 3b 0lep znn,  highest Var2 bin  : %6.1f  \n", ti, fit_znn__0lep_3b_hh1 ) ;


      printf("\n\n") ;
      printf(" toy %4d : Fit nb 0lep susy, highest Var1 bin : %6.1f  \n", ti, fit_susy_0lep_nb_hm1 ) ;
      printf(" toy %4d : Fit nb 0lep ttwj, highest Var1 bin : %6.1f  \n", ti, fit_ttwj_0lep_nb_hm1 ) ;
      printf(" toy %4d : Fit nb 0lep qcd,  highest Var1 bin : %6.1f  \n", ti, fit_qcd__0lep_nb_hm1 ) ;
      printf(" toy %4d : Fit nb 0lep znn,  highest Var1 bin : %6.1f  \n", ti, fit_znn__0lep_nb_hm1 ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 1b 0lep susy, highest Var1 bin : %6.1f  \n", ti, fit_susy_0lep_1b_hm1 ) ;
      printf(" toy %4d : Fit 1b 0lep ttwj, highest Var1 bin : %6.1f  \n", ti, fit_ttwj_0lep_1b_hm1 ) ;
      printf(" toy %4d : Fit 1b 0lep qcd,  highest Var1 bin : %6.1f  \n", ti, fit_qcd__0lep_1b_hm1 ) ;
      printf(" toy %4d : Fit 1b 0lep znn,  highest Var1 bin : %6.1f  \n", ti, fit_znn__0lep_1b_hm1 ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 2b 0lep susy, highest Var1 bin : %6.1f  \n", ti, fit_susy_0lep_2b_hm1 ) ;
      printf(" toy %4d : Fit 2b 0lep ttwj, highest Var1 bin : %6.1f  \n", ti, fit_ttwj_0lep_2b_hm1 ) ;
      printf(" toy %4d : Fit 2b 0lep qcd,  highest Var1 bin : %6.1f  \n", ti, fit_qcd__0lep_2b_hm1 ) ;
      printf(" toy %4d : Fit 2b 0lep znn,  highest Var1 bin : %6.1f  \n", ti, fit_znn__0lep_2b_hm1 ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 3b 0lep susy, highest Var1 bin : %6.1f  \n", ti, fit_susy_0lep_3b_hm1 ) ;
      printf(" toy %4d : Fit 3b 0lep ttwj, highest Var1 bin : %6.1f  \n", ti, fit_ttwj_0lep_3b_hm1 ) ;
      printf(" toy %4d : Fit 3b 0lep qcd,  highest Var1 bin : %6.1f  \n", ti, fit_qcd__0lep_3b_hm1 ) ;
      printf(" toy %4d : Fit 3b 0lep znn,  highest Var1 bin : %6.1f  \n", ti, fit_znn__0lep_3b_hm1 ) ;


      printf("\n\n") ;
      printf(" toy %4d : Fit nb 0lep susy, highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_susy_0lep_nb_hm2_hh2 ) ;
      printf(" toy %4d : Fit nb 0lep ttwj, highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_ttwj_0lep_nb_hm2_hh2 ) ;
      printf(" toy %4d : Fit nb 0lep qcd,  highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_qcd__0lep_nb_hm2_hh2 ) ;
      printf(" toy %4d : Fit nb 0lep znn,  highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_znn__0lep_nb_hm2_hh2 ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 1b 0lep susy, highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_susy_0lep_1b_hm2_hh2 ) ;
      printf(" toy %4d : Fit 1b 0lep ttwj, highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_ttwj_0lep_1b_hm2_hh2 ) ;
      printf(" toy %4d : Fit 1b 0lep qcd,  highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_qcd__0lep_1b_hm2_hh2 ) ;
      printf(" toy %4d : Fit 1b 0lep znn,  highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_znn__0lep_1b_hm2_hh2 ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 2b 0lep susy, highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_susy_0lep_2b_hm2_hh2 ) ;
      printf(" toy %4d : Fit 2b 0lep ttwj, highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_ttwj_0lep_2b_hm2_hh2 ) ;
      printf(" toy %4d : Fit 2b 0lep qcd,  highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_qcd__0lep_2b_hm2_hh2 ) ;
      printf(" toy %4d : Fit 2b 0lep znn,  highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_znn__0lep_2b_hm2_hh2 ) ;

      printf("\n") ;
      printf(" toy %4d : Fit 3b 0lep susy, highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_susy_0lep_3b_hm2_hh2 ) ;
      printf(" toy %4d : Fit 3b 0lep ttwj, highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_ttwj_0lep_3b_hm2_hh2 ) ;
      printf(" toy %4d : Fit 3b 0lep qcd,  highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_qcd__0lep_3b_hm2_hh2 ) ;
      printf(" toy %4d : Fit 3b 0lep znn,  highest two Var1 and Var2 bins : %6.1f  \n", ti, fit_znn__0lep_3b_hm2_hh2 ) ;


      printf("\n") ;
      if ( qcdModelIndex == 2 || qcdModelIndex == 4 || qcdModelIndex == 5 ) {
         if ( nBinsVar2>=1 ) { printf(" toy %4d : Fit QCD 0lep/LDP ratio, H1 : (%5.1f +/- %4.1f)%%\n", ti, 100*fit_qcd_0lepLDP_ratio_H1, 100*fit_qcd_0lepLDP_ratio_H1_err) ; }
         if ( nBinsVar2>=2 ) { printf(" toy %4d : Fit QCD 0lep/LDP ratio, H2 : (%5.1f +/- %4.1f)%%\n", ti, 100*fit_qcd_0lepLDP_ratio_H2, 100*fit_qcd_0lepLDP_ratio_H2_err) ; }
         if ( nBinsVar2>=3 ) { printf(" toy %4d : Fit QCD 0lep/LDP ratio, H3 : (%5.1f +/- %4.1f)%%\n", ti, 100*fit_qcd_0lepLDP_ratio_H3, 100*fit_qcd_0lepLDP_ratio_H3_err) ; }
         if ( nBinsVar2>=4 ) { printf(" toy %4d : Fit QCD 0lep/LDP ratio, H4 : (%5.1f +/- %4.1f)%%\n", ti, 100*fit_qcd_0lepLDP_ratio_H4, 100*fit_qcd_0lepLDP_ratio_H4_err) ; }
      }
      if ( qcdModelIndex == 3 ) {
         printf(" toy %4d : Fit QCD 0lep/LDP ratio : (%5.1f +/- %4.1f)%%\n", ti, 100*fit_qcd_0lepLDP_ratio, 100*fit_qcd_0lepLDP_ratio_err) ;
      }
      if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {
         if ( nBinsVar1>=2 ) { printf(" toy %4d : Fit QCD SFqcd_met2 : %5.2f +/- %5.2f\n", ti, fit_SFqcd_met2, fit_SFqcd_met2_err) ; }
         if ( nBinsBjets>=2 ) { printf(" toy %4d : Fit QCD SFqcd_nb2 : %5.2f +/- %5.2f\n", ti, fit_SFqcd_nb2, fit_SFqcd_nb2_err) ; }
      }
      if ( qcdModelIndex == 4 ) {
         if ( nBinsVar1>=3 ) { printf(" toy %4d : Fit QCD SFqcd_met3 : %5.2f +/- %5.2f\n", ti, fit_SFqcd_met3, fit_SFqcd_met3_err) ; }
         if ( nBinsVar1>=4 ) { printf(" toy %4d : Fit QCD SFqcd_met4 : %5.2f +/- %5.2f\n", ti, fit_SFqcd_met4, fit_SFqcd_met4_err) ; }
         if ( nBinsBjets>=3 ) { printf(" toy %4d : Fit QCD SFqcd_nb3 : %5.2f +/- %5.2f\n", ti, fit_SFqcd_nb3, fit_SFqcd_nb3_err) ; }
         if ( nBinsBjets>=4 ) { printf(" toy %4d : Fit QCD SFqcd_nb4 : %5.2f +/- %5.2f\n", ti, fit_SFqcd_nb4, fit_SFqcd_nb4_err) ; }
      }

      printf(" toy %4d : Fit ttwj ZL/SL ratio : ( %6.3f +/- %5.3f)\n", ti, fit_ttwj_0lep1lep_ratio, fit_ttwj_0lep1lep_ratio_err ) ;

      printf("\n\n") ;


      return true ;

   } // processToyResult.


   //==================================================================================


   bool bookTree() {

      char branchstring[1000] ;

      char outfile[10000] ;
      sprintf( outfile, "%s/toy-results.root", outputDir ) ;
      ttfile = new TFile( outfile, "recreate" ) ;

      toytt = new TTree( "toytt", "Toy MC study" ) ;


      toytt -> Branch( "fit_susy_0lep_err", &fit_susy_0lep_err, "fit_susy_0lep_err/F" ) ;
      toytt -> Branch( "fit_susy_0lep_err_low", &fit_susy_0lep_err_low, "fit_susy_0lep_err_low/F" ) ;
      toytt -> Branch( "fit_susy_0lep_err_high", &fit_susy_0lep_err_high, "fit_susy_0lep_err_high/F" ) ;
      toytt -> Branch( "fit_susy_0lep_err_forpull", &fit_susy_0lep_err_forpull, "fit_susy_0lep_err_forpull/F" ) ;
      toytt -> Branch( "fit_susy_0lep", &fit_susy_0lep, "fit_susy_0lep/F" ) ;
      toytt -> Branch( "fit_susy_0lep_wsfs", &fit_susy_0lep_wsfs, "fit_susy_0lep_wsfs/F" ) ;



      sprintf( branchstring, "fit_susy_0lep_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "fit_susy_0lep_3da", &fit_susy_0lep_3da, branchstring ) ;

      sprintf( branchstring, "fit_ttwj_0lep_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "fit_ttwj_0lep_3da", &fit_ttwj_0lep_3da, branchstring ) ;

      sprintf( branchstring, "fit_qcd_0lep_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "fit_qcd_0lep_3da" , &fit_qcd_0lep_3da , branchstring ) ;

      sprintf( branchstring, "fit_znn_0lep_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "fit_znn_0lep_3da" , &fit_znn_0lep_3da , branchstring ) ;


      sprintf( branchstring, "fit_sf_ttwj_0lep_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "fit_sf_ttwj_0lep_3da", &fit_sf_ttwj_0lep_3da, branchstring ) ;

      sprintf( branchstring, "fit_sf_qcd_0lep_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "fit_sf_qcd_0lep_3da" , &fit_sf_qcd_0lep_3da , branchstring ) ;



      sprintf( branchstring, "mcval_susy_0lep_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "mcval_susy_0lep_3da", &susy_0lep          , branchstring ) ;

      sprintf( branchstring, "mcval_ttwj_0lep_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "mcval_ttwj_0lep_3da", &mcval_ttwj_0lep_3da, branchstring ) ;

      sprintf( branchstring, "mcval_qcd_0lep_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "mcval_qcd_0lep_3da" , &mcval_qcd_0lep_3da , branchstring ) ;

      sprintf( branchstring, "mcval_znn_0lep_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "mcval_znn_0lep_3da" , &mcval_znn_0lep_3da , branchstring ) ;



      if ( qcdModelIndex == 2 || qcdModelIndex == 4 || qcdModelIndex == 5 ) {
         if ( nBinsVar2 >=1 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H1", &fit_qcd_0lepLDP_ratio_H1, "fit_qcd_0lepLDP_ratio_H1/F" ) ; }
         if ( nBinsVar2 >=2 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H2", &fit_qcd_0lepLDP_ratio_H2, "fit_qcd_0lepLDP_ratio_H2/F" ) ; }
         if ( nBinsVar2 >=3 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H3", &fit_qcd_0lepLDP_ratio_H3, "fit_qcd_0lepLDP_ratio_H3/F" ) ; }
         if ( nBinsVar2 >=4 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H4", &fit_qcd_0lepLDP_ratio_H4, "fit_qcd_0lepLDP_ratio_H4/F" ) ; }

         if ( nBinsVar2 >=1 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H1_err", &fit_qcd_0lepLDP_ratio_H1_err, "fit_qcd_0lepLDP_ratio_H1_err/F" ) ; }
         if ( nBinsVar2 >=2 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H2_err", &fit_qcd_0lepLDP_ratio_H2_err, "fit_qcd_0lepLDP_ratio_H2_err/F" ) ; }
         if ( nBinsVar2 >=3 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H3_err", &fit_qcd_0lepLDP_ratio_H3_err, "fit_qcd_0lepLDP_ratio_H3_err/F" ) ; }
         if ( nBinsVar2 >=4 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H4_err", &fit_qcd_0lepLDP_ratio_H4_err, "fit_qcd_0lepLDP_ratio_H4_err/F" ) ; }
      }
      if ( qcdModelIndex == 3 ) {
         toytt -> Branch( "fit_qcd_0lepLDP_ratio", &fit_qcd_0lepLDP_ratio, "fit_qcd_0lepLDP_ratio/F" ) ;
         toytt -> Branch( "fit_qcd_0lepLDP_ratio_err", &fit_qcd_0lepLDP_ratio_err, "fit_qcd_0lepLDP_ratio_err/F" ) ;
      }
      if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {
         if ( nBinsVar1 >=2 ) { toytt -> Branch( "fit_SFqcd_met2", &fit_SFqcd_met2, "fit_SFqcd_met2/F" ) ; }
         if ( nBinsVar1 >=2 ) { toytt -> Branch( "fit_SFqcd_met2_err", &fit_SFqcd_met2_err, "fit_SFqcd_met2_err/F" ) ; }
         if ( nBinsBjets >=2 ) { toytt -> Branch( "fit_SFqcd_nb2", &fit_SFqcd_nb2, "fit_SFqcd_nb2/F" ) ; }
         if ( nBinsBjets >=2 ) { toytt -> Branch( "fit_SFqcd_nb2_err", &fit_SFqcd_nb2_err, "fit_SFqcd_nb2_err/F" ) ; }
      }
      if ( qcdModelIndex == 4 ) {
         if ( nBinsVar1 >=3 ) { toytt -> Branch( "fit_SFqcd_met3", &fit_SFqcd_met3, "fit_SFqcd_met3/F" ) ; }
         if ( nBinsVar1 >=4 ) { toytt -> Branch( "fit_SFqcd_met4", &fit_SFqcd_met4, "fit_SFqcd_met4/F" ) ; }
         if ( nBinsVar1 >=3 ) { toytt -> Branch( "fit_SFqcd_met3_err", &fit_SFqcd_met3_err, "fit_SFqcd_met3_err/F" ) ; }
         if ( nBinsVar1 >=4 ) { toytt -> Branch( "fit_SFqcd_met4_err", &fit_SFqcd_met4_err, "fit_SFqcd_met4_err/F" ) ; }
         if ( nBinsBjets >=3 ) { toytt -> Branch( "fit_SFqcd_nb3", &fit_SFqcd_nb3, "fit_SFqcd_nb3/F" ) ; }
         if ( nBinsBjets >=4 ) { toytt -> Branch( "fit_SFqcd_nb4", &fit_SFqcd_nb4, "fit_SFqcd_nb4/F" ) ; }
         if ( nBinsBjets >=3 ) { toytt -> Branch( "fit_SFqcd_nb3_err", &fit_SFqcd_nb3_err, "fit_SFqcd_nb3_err/F" ) ; }
         if ( nBinsBjets >=4 ) { toytt -> Branch( "fit_SFqcd_nb4_err", &fit_SFqcd_nb4_err, "fit_SFqcd_nb4_err/F" ) ; }
      }


      toytt -> Branch( "fit_ttwj_0lep1lep_ratio"     , &fit_ttwj_0lep1lep_ratio    , "fit_ttwj_0lep1lep_ratio/F"     ) ;
      toytt -> Branch( "fit_ttwj_0lep1lep_ratio_err" , &fit_ttwj_0lep1lep_ratio_err, "fit_ttwj_0lep1lep_ratio_err/F" ) ;

      toytt -> Branch( "fit_ttwj_1lepSig1lep_ratio"     , &fit_ttwj_1lepSig1lep_ratio    , "fit_ttwj_1lepSig1lep_ratio/F"     ) ;
      toytt -> Branch( "fit_ttwj_1lepSig1lep_ratio_err" , &fit_ttwj_1lepSig1lep_ratio_err, "fit_ttwj_1lepSig1lep_ratio_err/F" ) ;

      toytt -> Branch( "true_susy_0lep", &true_susy_0lep, "true_susy_0lep/F" ) ;

      toytt -> Branch( "fit_covqual_susyfloat", &fit_covqual_susyfloat, "fit_covqual_susyfloat/I" ) ;

      if ( doSignif ) {
         toytt -> Branch( "fit_susy_signif", &fit_susy_signif, "fit_susy_signif/F" ) ;
      
      	 sprintf( branchstring, "fit_susy_0lep_at0susy_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_susy_0lep_at0susy_3da", &fit_susy_0lep_at0susy_3da, branchstring ) ;

      	 sprintf( branchstring, "fit_ttwj_0lep_at0susy_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_ttwj_0lep_at0susy_3da", &fit_ttwj_0lep_at0susy_3da, branchstring ) ;

      	 sprintf( branchstring, "fit_qcd_0lep_at0susy_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_qcd_0lep_at0susy_3da" , &fit_qcd_0lep_at0susy_3da , branchstring ) ;

      	 sprintf( branchstring, "fit_znn_0lep_at0susy_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_znn_0lep_at0susy_3da" , &fit_znn_0lep_at0susy_3da , branchstring ) ;


      	 sprintf( branchstring, "fit_sf_ttwj_0lep_at0susy_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_sf_ttwj_0lep_at0susy_3da", &fit_sf_ttwj_0lep_at0susy_3da, branchstring ) ;

      	 sprintf( branchstring, "fit_sf_qcd_0lep_at0susy_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_sf_qcd_0lep_at0susy_3da" , &fit_sf_qcd_0lep_at0susy_3da , branchstring ) ;

         toytt -> Branch( "fit_ttwj_0lep1lep_ratio_at0susy", &fit_ttwj_0lep1lep_ratio_at0susy, "fit_ttwj_0lep1lep_ratio_at0susy/F" ) ;

         toytt -> Branch( "fit_ttwj_1lepSig1lep_ratio_at0susy", &fit_ttwj_1lepSig1lep_ratio_at0susy, "fit_ttwj_1lepSig1lep_ratio_at0susy/F" ) ;

         if ( qcdModelIndex == 2 || qcdModelIndex == 4 || qcdModelIndex == 5 ) {
            if ( nBinsVar2 >=1 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H1_at0susy", &fit_qcd_0lepLDP_ratio_H1_at0susy, "fit_qcd_0lepLDP_ratio_H1_at0susy/F" ) ; }
            if ( nBinsVar2 >=2 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H2_at0susy", &fit_qcd_0lepLDP_ratio_H2_at0susy, "fit_qcd_0lepLDP_ratio_H2_at0susy/F" ) ; }
            if ( nBinsVar2 >=3 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H3_at0susy", &fit_qcd_0lepLDP_ratio_H3_at0susy, "fit_qcd_0lepLDP_ratio_H3_at0susy/F" ) ; }
            if ( nBinsVar2 >=4 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H4_at0susy", &fit_qcd_0lepLDP_ratio_H4_at0susy, "fit_qcd_0lepLDP_ratio_H4_at0susy/F" ) ; }
         }
         if ( qcdModelIndex == 3 ) {
            toytt -> Branch( "fit_qcd_0lepLDP_ratio_at0susy", &fit_qcd_0lepLDP_ratio_at0susy, "fit_qcd_0lepLDP_ratio_at0susy/F" ) ;
         }
         if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {
            if ( nBinsVar1 >=2 ) { toytt -> Branch( "fit_SFqcd_met2_at0susy", &fit_SFqcd_met2_at0susy, "fit_SFqcd_met2_at0susy/F" ) ; }
            if ( nBinsBjets >=2 ) { toytt -> Branch( "fit_SFqcd_nb2_at0susy", &fit_SFqcd_nb2_at0susy, "fit_SFqcd_nb2_at0susy/F" ) ; }
         }
         if ( qcdModelIndex == 4 ) {
            if ( nBinsVar1 >=3 ) { toytt -> Branch( "fit_SFqcd_met3_at0susy", &fit_SFqcd_met3_at0susy, "fit_SFqcd_met3_at0susy/F" ) ; }
            if ( nBinsVar1 >=4 ) { toytt -> Branch( "fit_SFqcd_met4_at0susy", &fit_SFqcd_met4_at0susy, "fit_SFqcd_met4_at0susy/F" ) ; }
            if ( nBinsBjets >=3 ) { toytt -> Branch( "fit_SFqcd_nb3_at0susy", &fit_SFqcd_nb3_at0susy, "fit_SFqcd_nb3_at0susy/F" ) ; }
            if ( nBinsBjets >=4 ) { toytt -> Branch( "fit_SFqcd_nb4_at0susy", &fit_SFqcd_nb4_at0susy, "fit_SFqcd_nb4_at0susy/F" ) ; }
         }

      } // doSignif?

      if ( doUL ) {
         toytt -> Branch( "fit_susy_ul", &fit_susy_ul, "fit_susy_ul/F" ) ;
         toytt -> Branch( "fit_susy_ts_at_ul", &fit_susy_ts_at_ul, "fit_susy_ts_at_ul/F" ) ;

      	 sprintf( branchstring, "fit_susy_0lep_atUL_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_susy_0lep_atUL_3da", &fit_susy_0lep_atUL_3da, branchstring ) ;

      	 sprintf( branchstring, "fit_ttwj_0lep_atUL_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_ttwj_0lep_atUL_3da", &fit_ttwj_0lep_atUL_3da, branchstring ) ;

      	 sprintf( branchstring, "fit_qcd_0lep_atUL_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_qcd_0lep_atUL_3da" , &fit_qcd_0lep_atUL_3da , branchstring ) ;

      	 sprintf( branchstring, "fit_znn_0lep_atUL_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_znn_0lep_atUL_3da" , &fit_znn_0lep_atUL_3da , branchstring ) ;


      	 sprintf( branchstring, "fit_sf_ttwj_0lep_atUL_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_sf_ttwj_0lep_atUL_3da", &fit_sf_ttwj_0lep_atUL_3da, branchstring ) ;

      	 sprintf( branchstring, "fit_sf_qcd_0lep_atUL_3da[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      	 toytt -> Branch( "fit_sf_qcd_0lep_atUL_3da" , &fit_sf_qcd_0lep_atUL_3da , branchstring ) ;

         toytt -> Branch( "fit_ttwj_0lep1lep_ratio_atUL", &fit_ttwj_0lep1lep_ratio_atUL, "fit_ttwj_0lep1lep_ratio_atUL/F" ) ;
         toytt -> Branch( "fit_ttwj_1lepSig1lep_ratio_atUL", &fit_ttwj_1lepSig1lep_ratio_atUL, "fit_ttwj_1lepSig1lep_ratio_atUL/F" ) ;

         if ( qcdModelIndex == 2 || qcdModelIndex == 4 || qcdModelIndex == 5 ) {
            if ( nBinsVar2 >=1 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H1_atUL", &fit_qcd_0lepLDP_ratio_H1_atUL, "fit_qcd_0lepLDP_ratio_H1_atUL/F" ) ; }
            if ( nBinsVar2 >=2 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H2_atUL", &fit_qcd_0lepLDP_ratio_H2_atUL, "fit_qcd_0lepLDP_ratio_H2_atUL/F" ) ; }
            if ( nBinsVar2 >=3 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H3_atUL", &fit_qcd_0lepLDP_ratio_H3_atUL, "fit_qcd_0lepLDP_ratio_H3_atUL/F" ) ; }
            if ( nBinsVar2 >=4 ) { toytt -> Branch( "fit_qcd_0lepLDP_ratio_H4_atUL", &fit_qcd_0lepLDP_ratio_H4_atUL, "fit_qcd_0lepLDP_ratio_H4_atUL/F" ) ; }
         }
         if ( qcdModelIndex == 3 ) {
            toytt -> Branch( "fit_qcd_0lepLDP_ratio_atUL", &fit_qcd_0lepLDP_ratio_atUL, "fit_qcd_0lepLDP_ratio_atUL/F" ) ;
         }
         if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {
            if ( nBinsVar1 >=2 ) { toytt -> Branch( "fit_SFqcd_met2_atUL", &fit_SFqcd_met2_atUL, "fit_SFqcd_met2_atUL/F" ) ; }
            if ( nBinsBjets >=2 ) { toytt -> Branch( "fit_SFqcd_nb2_atUL", &fit_SFqcd_nb2_atUL, "fit_SFqcd_nb2_atUL/F" ) ; }
         }
         if ( qcdModelIndex == 4 ) {
            if ( nBinsVar1 >=3 ) { toytt -> Branch( "fit_SFqcd_met3_atUL", &fit_SFqcd_met3_atUL, "fit_SFqcd_met3_atUL/F" ) ; }
            if ( nBinsVar1 >=4 ) { toytt -> Branch( "fit_SFqcd_met4_atUL", &fit_SFqcd_met4_atUL, "fit_SFqcd_met4_atUL/F" ) ; }
            if ( nBinsBjets >=3 ) { toytt -> Branch( "fit_SFqcd_nb3_atUL", &fit_SFqcd_nb3_atUL, "fit_SFqcd_nb3_atUL/F" ) ; }
            if ( nBinsBjets >=4 ) { toytt -> Branch( "fit_SFqcd_nb4_atUL", &fit_SFqcd_nb4_atUL, "fit_SFqcd_nb4_atUL/F" ) ; }
         }
      } // doUL?


      sprintf( branchstring, "Nobs_0lep[%d][%d][%d]/I", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "Nobs_0lep", Nobs_0lep, branchstring ) ;

      sprintf( branchstring, "Nobs_1lepSig[%d][%d][%d]/I", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "Nobs_1lepSig", Nobs_1lepSig, branchstring ) ;

      sprintf( branchstring, "Nobs_1lep[%d][%d][%d]/I", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "Nobs_1lep", Nobs_1lep, branchstring ) ;

      sprintf( branchstring, "Nobs_ldp[%d][%d][%d]/I", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "Nobs_ldp" , Nobs_ldp , branchstring ) ;

      sprintf( branchstring, "Nobs_Zee[%d][%d]/I", nBinsVar1, nBinsVar2 ) ;
      toytt -> Branch( "Nobs_Zee" , Nobs_Zee , branchstring ) ;

      sprintf( branchstring, "Nobs_Zmm[%d][%d]/I", nBinsVar1, nBinsVar2 ) ;
      toytt -> Branch( "Nobs_Zmm" , Nobs_Zmm , branchstring ) ;



      sprintf( branchstring, "toy_mean_N_0lep[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "toy_mean_N_0lep", toy_mean_N_0lep, branchstring ) ;

      sprintf( branchstring, "toy_mean_N_1lepSig[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "toy_mean_N_1lepSig", toy_mean_N_1lepSig, branchstring ) ;

      sprintf( branchstring, "toy_mean_N_1lep[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "toy_mean_N_1lep", toy_mean_N_1lep, branchstring ) ;

      sprintf( branchstring, "toy_mean_N_ldp[%d][%d][%d]/F", nBinsVar1, nBinsVar2, nBinsBjets ) ;
      toytt -> Branch( "toy_mean_N_ldp", toy_mean_N_ldp, branchstring ) ;

      sprintf( branchstring, "toy_mean_N_Zee[%d][%d]/F", nBinsVar1, nBinsVar2 ) ;
      toytt -> Branch( "toy_mean_N_Zee", toy_mean_N_Zee, branchstring ) ;

      sprintf( branchstring, "toy_mean_N_Zmm[%d][%d]/F", nBinsVar1, nBinsVar2 ) ;
      toytt -> Branch( "toy_mean_N_Zmm", toy_mean_N_Zmm, branchstring ) ;





      toytt -> Branch( "fit_chi2_overall", &fit_chi2_overall, "fit_chi2_overall/F" ) ;
      toytt -> Branch( "fit_chi2_obs"    , &fit_chi2_obs    , "fit_chi2_obs/F"     ) ;
      toytt -> Branch( "fit_chi2_np"     , &fit_chi2_np     , "fit_chi2_np/F"      ) ;
      toytt -> Branch( "fit_chi2_prob"   , &fit_chi2_prob   , "fit_chi2_prob/F"    ) ;






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

      fit_susy_ul = susy_poi_g2 ;
      fit_susy_ts_at_ul = testStat_g2 ;

      printf("  SUSY upper limit is : %5.2f 0lep events.\n", fit_susy_ul ) ;


      fit_ttwj_0lep1lep_ratio_atUL = rrv_ttwj_0lep1lep_ratio->getVal();

      if ( qcdModelIndex == 2 || qcdModelIndex == 4 || qcdModelIndex == 5 ) {
         if ( nBinsVar2 >= 1 ) { fit_qcd_0lepLDP_ratio_H1_atUL     = rrv_qcd_0lepLDP_ratio_H1 -> getVal() ; }
         if ( nBinsVar2 >= 2 ) { fit_qcd_0lepLDP_ratio_H2_atUL     = rrv_qcd_0lepLDP_ratio_H2 -> getVal() ; }
         if ( nBinsVar2 >= 3 ) { fit_qcd_0lepLDP_ratio_H3_atUL     = rrv_qcd_0lepLDP_ratio_H3 -> getVal() ; }
         if ( nBinsVar2 >= 4 ) { fit_qcd_0lepLDP_ratio_H4_atUL     = rrv_qcd_0lepLDP_ratio_H4 -> getVal() ; }
      }
      if ( qcdModelIndex == 3 ) { fit_qcd_0lepLDP_ratio_atUL     = rrv_qcd_0lepLDP_ratio -> getVal() ; }
      if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {
         if ( nBinsVar1 >=2 ) { fit_SFqcd_met2_atUL     = rrv_SFqcd_met2 -> getVal() ; }
         if ( nBinsBjets >=2 ) { fit_SFqcd_nb2_atUL     = rrv_SFqcd_nb2 -> getVal() ; }
      }
      if ( qcdModelIndex == 4 ) {
         if ( nBinsVar1 >=3 ) { fit_SFqcd_met3_atUL     = rrv_SFqcd_met3 -> getVal() ; }
         if ( nBinsVar1 >=4 ) { fit_SFqcd_met4_atUL     = rrv_SFqcd_met4 -> getVal() ; }
         if ( nBinsBjets >=3 ) { fit_SFqcd_nb3_atUL     = rrv_SFqcd_nb3 -> getVal() ; }
         if ( nBinsBjets >=4 ) { fit_SFqcd_nb4_atUL     = rrv_SFqcd_nb4 -> getVal() ; }
      }

      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
	    for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               if (mbi==3 && hbi==0) {
	          fit_susy_0lep_atUL_3da[mbi][hbi][bbi] = 0.0;
	          fit_ttwj_0lep_atUL_3da[mbi][hbi][bbi] = 0.0;
	          fit_qcd_0lep_atUL_3da[mbi][hbi][bbi] = 0.0;
	          fit_znn_0lep_atUL_3da[mbi][hbi][bbi] = 0.0;
	          fit_sf_ttwj_0lep_atUL_3da[mbi][hbi][bbi] = 1.0;
	          fit_sf_qcd_0lep_atUL_3da[mbi][hbi][bbi] = 1.0;
	       } else {

                  char vname[1000] ;
                  RooAbsReal* rar ;
                  RooAbsReal* effsf  ;
                  RooAbsReal* btagsf  ;


                  sprintf( vname, "mu_susy_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  rar = workspace -> function( vname ) ;
                  if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; }
                  sprintf( vname, "btageff_sf_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  btagsf = (RooAbsReal*) workspace -> obj( vname ) ;
                  if ( btagsf == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; }
                  sprintf( vname, "eff_sf_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  effsf = (RooAbsReal*) workspace -> obj( vname ) ;
                  if ( effsf == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; }
                  double nsusy = ( btagsf -> getVal() ) * ( effsf -> getVal() ) * ( rar -> getVal() ) ;

                  fit_susy_0lep_wsfs += nsusy ;

                  fit_susy_0lep_atUL_3da[mbi][hbi][bbi] = nsusy ;

                  sprintf( vname, "mu_ttwj_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  rar = (RooAbsReal*) workspace -> obj( vname ) ;
                  if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; }
                  double nttwj = rar -> getVal() ;
                  fit_ttwj_0lep_atUL_3da[mbi][hbi][bbi] = nttwj ;


                  sprintf( vname, "mu_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  rar = workspace -> function( vname ) ;
                  if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; }
                  double nqcd = rar -> getVal() ;
                  fit_qcd_0lep_atUL_3da[mbi][hbi][bbi] = nqcd ;

                  sprintf( vname, "mu_znn_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  rar = (RooAbsReal*) workspace -> obj( vname ) ;
                  if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; }
                  double nznn = rar -> getVal() ;
                  fit_znn_0lep_atUL_3da[mbi][hbi][bbi] = nznn ;

                  sprintf( vname, "sf_ttwj_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  rar = (RooAbsReal*) workspace -> obj( vname ) ;
                  if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ;  }
                  fit_sf_ttwj_0lep_atUL_3da[mbi][hbi][bbi] = rar -> getVal() ;

                  sprintf( vname, "sf_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  rar = (RooAbsReal*) workspace -> obj( vname ) ;
                  if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; }
                  fit_sf_qcd_0lep_atUL_3da[mbi][hbi][bbi] = rar -> getVal() ;

               }
            }
         }
      }

      return fit_susy_ul ;


   } // findUL.

   //==================================================================================

   bool setReinitValues( RooFitResult* fitResult ) {

      const RooArgList floatPars = fitResult -> floatParsFinal() ;

      printf("\n\n ==== Will reinitialize to these values before every toy.\n\n") ;

      nFloatParInitVal = 0 ;

      initFit_R_ttwj_0lep_over_1lep = -1. ;
      initFit_qcd_0lepLDP_ratio = -1. ;
      for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
         initFit_qcd_0lepLDP_ratio_HT[hbi] = -1. ;
      } // hbi.
      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         initFit_SFqcd_met[mbi] = 1. ;
      } // mbi.
      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
         initFit_SFqcd_nb[bbi] = 1. ;
      } // bbi.

      TIterator* parIter = floatPars.createIterator() ;
      while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {

         sprintf( floatParName[nFloatParInitVal], "%s", par->GetName() ) ;
         floatParInitVal[nFloatParInitVal] = par->getVal() ;

         nFloatParInitVal++ ;

         printf(" %20s : %8.2f\n", par->GetName(), par->getVal() ) ;

         if ( strcmp( par->GetName(), "ttwj_0lep1lep_ratio" ) == 0 ) {
            printf("    Setting initFit_R_ttwj_0lep_over_1lep to %5.3f\n", par->getVal() ) ;
            initFit_R_ttwj_0lep_over_1lep = par->getVal() ;
         }
         if ( qcdModelIndex == 3 ) {
            if ( strcmp( par->GetName(), "qcd_0lepLDP_ratio" ) == 0 ) {
               printf("    Setting initFit_qcd_0lepLDP_ratio to %5.3f\n", par->getVal() ) ;
               initFit_qcd_0lepLDP_ratio = par->getVal() ;
            }
         } else if ( qcdModelIndex == 2 ) {
            for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
               char pname[100] ;
               sprintf( pname, "qcd_0lepLDP_ratio_H%d", hbi+1 ) ;
               if ( strcmp( par->GetName(), pname ) == 0 ) {
                  printf("    Setting initFit_%s to %5.3f\n", pname, par->getVal() ) ;
                  initFit_qcd_0lepLDP_ratio_HT[hbi] = par->getVal() ;
               }
            } // hbi.
         } else if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {
            for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
               char pname[100] ;
               sprintf( pname, "qcd_0lepLDP_ratio_H%d", hbi+1 ) ;
               if ( strcmp( par->GetName(), pname ) == 0 ) {
                  printf("    Setting initFit_%s to %5.3f\n", pname, par->getVal() ) ;
                  initFit_qcd_0lepLDP_ratio_HT[hbi] = par->getVal() ;
               }
            } // hbi.
            for ( int mbi=1; mbi<nBinsVar1; mbi++ ) { //-- this is ok for model 5
               char pname[100] ;
               sprintf( pname, "SFqcd_met%d", mbi+1 ) ;
               if ( strcmp( par->GetName(), pname ) == 0 ) {
                  printf("    Setting initFit_%s to %5.3f\n", pname, par->getVal() ) ;
                  initFit_SFqcd_met[mbi] = par->getVal() ;
               }
            } // mbi.
            for ( int bbi=1; bbi<nBinsBjets; bbi++ ) { //-- this is ok for model 5
               char pname[100] ;
               sprintf( pname, "SFqcd_nb%d", bbi+1 ) ;
               if ( strcmp( par->GetName(), pname ) == 0 ) {
                  printf("    Setting initFit_%s to %5.3f\n", pname, par->getVal() ) ;
                  initFit_SFqcd_nb[bbi] = par->getVal() ;
               }
            } // bbi.
         } // QCD model 4 or 5?

      } // par iterator.
      printf("\n\n") ;

      if ( initFit_R_ttwj_0lep_over_1lep < 0 ) {

      // printf("\n\n *** Did not find floating parameter ttwj_0lep1lep_ratio.  Can't continue.\n\n") ;
      // return false ;

        //zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
         printf("\n\n *** Did not find floating parameter ttwj_0lep1lep_ratio.  Assuming its fixed.\n\n") ;
         RooAbsReal* rar = (RooAbsReal*) workspace->obj( "ttwj_0lep1lep_ratio" ) ;
         if ( rar == 0x0 ) {
            printf("\n\n *** Did not find fixed parameter ttwj_0lep1lep_ratio.  Can't continue.\n\n") ;
         } else {
            initFit_R_ttwj_0lep_over_1lep = rar->getVal() ;
         }
        //zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

      }
      if ( qcdModelIndex == 3 ) {
         if ( initFit_qcd_0lepLDP_ratio < 0 ) {
            printf("\n\n *** Did not find floating parameter qcd_0lepLDP_ratio.  Can't continue.\n\n") ;
            return false ;
         }
      } else if ( qcdModelIndex == 2 ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
            if ( initFit_qcd_0lepLDP_ratio_HT[hbi] < 0 ) {
               printf("\n\n *** Did not find floating parameter initFit_qcd_0lepLDP_ratio_H%d.  Can't continue.\n\n", hbi+1 ) ;
               return false ;
            }
         } // hbi.
      } else if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
            if ( initFit_qcd_0lepLDP_ratio_HT[hbi] < 0 ) {
               printf("\n\n *** Did not find floating parameter initFit_qcd_0lepLDP_ratio_H%d.  Can't continue.\n\n", hbi+1 ) ;
               return false ;
            }
         } // hbi.
         for ( int mbi=1; mbi<nBinsVar1; mbi++ ) {
            if ( qcdModelIndex == 5 && mbi>1 ) continue ;
            if ( initFit_SFqcd_met[mbi] == 1 ) {
               printf("\n\n *** Did not find floating parameter SFqcd_met%d.  Can't continue.\n\n", mbi+1 ) ;
               //return false ;
            }
         } // mbi.
         for ( int bbi=1; bbi<nBinsBjets; bbi++ ) {
            if ( qcdModelIndex == 5 && bbi>1 ) continue ;
            if ( initFit_SFqcd_nb[bbi] == 1 ) {
               printf("\n\n *** Did not find floating parameter SFqcd_nb%d.  Can't continue.\n\n", bbi+1 ) ;
               //return false ;
            }
         } // bbi.
      }

      return true ;

   } // setReinitValues.

   //==================================================================================

   bool reinitFloatPars() {

      for ( int pi=0; pi<nFloatParInitVal; pi++ ) {
         RooRealVar* par = workspace->var( floatParName[pi] ) ;
         if ( par == 0x0 ) {
            printf("\n\n *** can't find float par %s in workspace.\n\n", floatParName[pi] ) ;
            //return false ;
         }
         else { par->setVal( floatParInitVal[pi] ) ; }
      } // pi.
      rrv_susy_poi -> setVal( 0. ) ;

      return true ;

   } // reinitFloatPars.

   //==================================================================================

   bool saveToyDatfile( int toyIndex, RooDataSet* toyds ) {

      TString toyoutfilename( datfile ) ;
      char replacestring[1000] ;
      sprintf( replacestring, "-toy%04d.dat", toyIndex ) ;
      toyoutfilename.ReplaceAll( ".dat", replacestring ) ;
      toyoutfilename.ReplaceAll( "datfiles", outputDir ) ;

      printf("\n\n Saving toy dataset in %s\n\n", toyoutfilename.Data() ) ;
      char command[1000] ;
      sprintf( command, "cp %s %s", datfile, toyoutfilename.Data() );
      gSystem->Exec( command ) ;

      const RooArgSet* dsras = toyds->get() ;
      TIterator* obsIter = dsras->createIterator() ;

      TString ufvname = outputDir ;
      ufvname += ".txt";

      while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {
	updateFileValue( toyoutfilename, obs->GetName(), obs->getVal(), ufvname ) ;
      }

      return true ;

   } // saveToyDatfile

   //==================================================================================

   bool readAndSetMCVals() {


      //--- Note: the truth susy values do NOT come from the mcvals root file.  They come from the SUSY .dat file
      //          and the ttree variables are set in addSusyToExpectedObs.



      printf("\n\n Reading in MC values from %s\n\n", mcvals_rootfile ) ;

      TFile mcfile( mcvals_rootfile, "READ" ) ;
      if ( ! mcfile.IsOpen() ) {
         return false ;
      }

      TH1F* httwj[nBinsBjets] ;
      TH1F* hqcd[nBinsBjets] ;
      TH1F* hznn[nBinsBjets] ;

      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

         char hname[1000] ;
         TH1F* hist ;

         sprintf( hname, "hmctruth_ttwj_0lep_%db", bbi+1 ) ;
         hist = (TH1F*) mcfile.Get( hname ) ;
         if ( hist == 0x0 ) { printf("\n\n *** Missing MC histogram %s\n\n", hname ) ; return false ; }
         httwj[bbi] = hist ;

         sprintf( hname, "hmctruth_qcd_0lep_%db", bbi+1 ) ;
         hist = (TH1F*) mcfile.Get( hname ) ;
         if ( hist == 0x0 ) { printf("\n\n *** Missing MC histogram %s\n\n", hname ) ; return false ; }
         hqcd[bbi] = hist ;

         sprintf( hname, "hmctruth_znn_0lep_%db", bbi+1 ) ;
         hist = (TH1F*) mcfile.Get( hname ) ;
         if ( hist == 0x0 ) { printf("\n\n *** Missing MC histogram %s\n\n", hname ) ; return false ; }
         hznn[bbi] = hist ;

      } // bbi.

      int mcnBinsVar1(0), mcnBinsVar2(0) ;

      char binlabel[1000] ;
      sprintf( binlabel, "%s", httwj[0] -> GetXaxis() -> GetBinLabel( httwj[0]->GetNbinsX() - 1 ) ) ;
      sscanf( binlabel, "0lep_M%d_H%d_1b", &mcnBinsVar1, &mcnBinsVar2 ) ;
      if ( mcnBinsVar1 != nBinsVar1 || mcnBinsVar2 != nBinsVar2 ) {
         printf("\n\n *** Incompatible Var1,Var2 binning.\n") ;
         printf("  toymc setting   : Var1=%d, Var2=%d\n", nBinsVar1, nBinsVar2) ;
         printf("  mcvals rootfile : Var1=%d, Var2=%d\n", mcnBinsVar1, mcnBinsVar2) ;
         printf("\n\n") ;
         return false ;
      }


      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
         for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               int hbin = 1 + hbi + mbi*(nBinsVar2+1) + 1 ;

               mcval_ttwj_0lep_3da[mbi][hbi][bbi] = httwj[bbi] -> GetBinContent( hbin ) ;
               mcval_qcd_0lep_3da [mbi][hbi][bbi] = hqcd[bbi]  -> GetBinContent( hbin ) ;
               mcval_znn_0lep_3da [mbi][hbi][bbi] = hznn[bbi]  -> GetBinContent( hbin ) ;


            } // bbi.
         } // hbi.
      } // mbi.




      return true ;

   } // readAndSetMCVals


   //==================================================================================


   bool setSFVals() {


      for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
	for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {

            char ename[1000] ;
            float trigeff ;

            sprintf( ename, "trigeff_val_0L_M%d_H%d", mbi+1, hbi+1 ) ;
            if ( !getFileValue( datfile, ename, trigeff ) ) {
               return false ;
            }
            trigeff_0lep[mbi][hbi] = trigeff ;

            sprintf( ename, "trigeff_val_1L_M%d_H%d", mbi+1, hbi+1 ) ;
            if ( !getFileValue( datfile, ename, trigeff ) ) {
               return false ;
            }
            trigeff_1lep[mbi][hbi] = trigeff ;

            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               float val ;
               char sfname[10000] ;

               sprintf( sfname, "sf_ttwj_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               if ( !getFileValue( datfile, sfname, val ) ) {
                  return false ;
               }
               sf_ttwj[mbi][hbi][bbi] = val ;

               sprintf( sfname, "sf_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               if ( !getFileValue( datfile, sfname, val ) ) {
                  return false ;
               }
               sf_qcd[mbi][hbi][bbi] = val ;



               sprintf( sfname, "ttwj_mc_ldpover0lep_ratio_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               if ( !getFileValue( datfile, sfname, val ) ) {
                  return false ;
               }
               ttwj_mc_ldpover0lep_ratio[mbi][hbi][bbi] = val ;

               sprintf( sfname, "znn_mc_ldpover0lep_ratio_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
               if ( !getFileValue( datfile, sfname, val ) ) {
                  return false ;
               }
               znn_mc_ldpover0lep_ratio[mbi][hbi][bbi] = val ;

            } // bbi.
         } // hbi.
      } // mbi.


      return true ;

   } // setSFVals



   //==================================================================================

   double getChi2Obs( const char* obsname, const char* modelname ) {

      double obsval, modelval, chi(0.) ;
      RooAbsReal* rar ;

      rar = (RooAbsReal*) workspace -> obj( obsname ) ;
      if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", obsname ) ; return 0. ; }
      obsval = rar -> getVal() ;

      rar = (RooAbsReal*) workspace -> obj( modelname ) ;
      if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", modelname ) ; return 0. ; }
      modelval = rar -> getVal() ;

      if ( obsval > 0 ) {
         chi = (obsval - modelval) / sqrt(obsval) ;
         fit_chi2_obs += chi*chi ;
      }

      return chi*chi ;

   } // getChi2Obs.


   //==================================================================================

   double getChi2GausNP( const char* npname ) {

      double npval, npmean, npsigma, chi(0.) ;
      char vname[1000] ;
      RooAbsReal* rar ;

      rar = (RooAbsReal*) workspace -> obj( npname ) ;
      if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", npname ) ; return 0. ; }
      npval = rar -> getVal() ;

      sprintf( vname, "mean_%s", npname ) ;
      rar = (RooAbsReal*) workspace -> obj( vname ) ;
      if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; return 0. ; }
      npmean = rar -> getVal() ;

      sprintf( vname, "sigma_%s", npname ) ;
      rar = (RooAbsReal*) workspace -> obj( vname ) ;
      if ( rar == 0x0 ) { printf("\n\n *** missing var %s\n\n", vname ) ; return 0. ; }
      npsigma = rar -> getVal() ;

      if ( npsigma > 0. ) {
         chi = (npval-npmean)/npsigma ;
         fit_chi2_np += chi*chi ;
      }

      return chi*chi ;

   } // getChi2GausNP.


   //==================================================================================

   double getChi2BetaNP( const char* npname ) {

      double  chi(0.) ;

      double mode, rms, alpha, beta ;
      getBetaModeRMS( npname, mode, rms, alpha, beta ) ;
      RooAbsReal* np = (RooAbsReal*) workspace->obj( npname ) ;
      if ( np == 0x0 ) {
         printf("\n\n *** missing nuisance parameter? %s\n\n", npname ) ;
         return 0. ;
      }
      double npVal = np->getVal() ;
      chi = (npVal-mode)/rms ;

      return chi*chi ;

   } // getChi2GausNP.


   //==================================================================================


  bool getBetaModeRMS( const char* parName, double &mode, double &rms, double &alpha, double &beta ) {

     mode = 1.0 ;
     rms = 0.0 ;

     char varname[1000] ;

     sprintf( varname, "alpha_%s", parName ) ;
     RooAbsReal* alphaPar = (RooAbsReal*) workspace->obj( varname ) ;
     if ( alphaPar == 0x0 ) {
        printf("\n\n *** getNPModeRMS : can't find alpha for %s\n\n", parName ) ;
        return false ;
     }
     alpha = alphaPar->getVal() ;

     sprintf( varname, "beta_%s", parName ) ;
     RooAbsReal* betaPar = (RooAbsReal*) workspace->obj( varname ) ;
     if ( betaPar == 0x0 ) {
        printf("\n\n *** getNPModeRMS : can't find beta for %s\n\n", parName ) ;
        return false ;
     }
     beta = betaPar->getVal() ;

     mode = (alpha - 1.)/(alpha+beta -2.) ;

     rms = sqrt( alpha * beta / ( pow(alpha + beta,2) * (alpha + beta + 1) ) ) ;

     return true ;

  } // getNPModeRMS.


//==========================================================================================







