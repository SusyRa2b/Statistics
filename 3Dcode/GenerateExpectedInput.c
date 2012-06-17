
#include <fstream>
#include <stdio.h>



   void GenerateExpectedInput( const char* infilename = "Input-met4-ht4.dat",
                               const char* outfilename = "expected-inputs.dat",
                               double mgl = 900., double mlsp = 300., double nsusy0leptarget = 50.,
                               const char* insusyfilename = "Susy-mgl900-mlsp300-met4-ht4-v2.dat" ) {



      printf("\n\n Will generate %s from %s\n", outfilename, infilename ) ;
      if ( mgl>0. && mlsp>0. ) {
         printf("    Including %5.1f events from mgl=%.0f, mlsp=%.0f from %s\n",
           nsusy0leptarget, mgl, mlsp, insusyfilename ) ;
      }

     //============================================

      char pname[1000] ;
      float pval ;

      //--- Read in the file and determine the binning configuration.

      FILE* infp ;
      if ( (infp=fopen( infilename,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infilename ) ;
         return ;
      }

      //--- read in description line.
      printf("\n\n --- Header line:\n") ;
      char c(0) ;
      while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
      printf("\n\n") ;

      int nBinsMET(0), nBinsHT(0), nBinsBjets(0) ;
      int nread(1) ;
      while ( nread>0 ) {
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
         /// printf("  %d : %s %g\n", nread, pname, pval ) ;
         int metbin, htbin, nbbin ;
         int nmatch = sscanf( pname, "N_0lep_M%d_H%d_%db", &metbin, &htbin, &nbbin ) ;
         if ( nmatch > 0 ) {
            /// printf("  nmatch=%d : Observable: metbin=%d, htbin=%d, nbbin=%d\n", nmatch, metbin, htbin, nbbin ) ;
            if ( metbin > nBinsMET   ) { nBinsMET   = metbin ; }
            if ( htbin  > nBinsHT    ) { nBinsHT    = htbin ; }
            if ( nbbin  > nBinsBjets ) { nBinsBjets = nbbin ; }
         }
      }
      fclose(infp) ;

      printf("\n  File has this configuration:  nBinsMET = %d ,  nBinsHT = %d ,  nBinsBjets = %d\n\n",
           nBinsMET, nBinsHT, nBinsBjets ) ;


      //==== Set up data arrays.

      int*** N_0lep = new int**[nBinsMET] ;
      int*** N_0lep_input = new int**[nBinsMET] ;
      int*** N_1lep = new int**[nBinsMET] ;
      int*** N_ldp  = new int**[nBinsMET] ;
      float*** N_mc_ldp = new float**[nBinsMET] ;
      float*** susy_0lep = new float**[nBinsMET] ;
      float*** susy_1lep = new float**[nBinsMET] ;
      float*** susy_ldp  = new float**[nBinsMET] ;
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         N_0lep[mbi] = new int*[nBinsHT] ;
         N_0lep_input[mbi] = new int*[nBinsHT] ;
         N_1lep[mbi] = new int*[nBinsHT] ;
         N_ldp[mbi]  = new int*[nBinsHT] ;
         N_mc_ldp[mbi] = new float*[nBinsHT] ;
         susy_0lep[mbi] = new float*[nBinsHT] ;
         susy_1lep[mbi] = new float*[nBinsHT] ;
         susy_ldp[mbi]  = new float*[nBinsHT] ;
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            N_0lep[mbi][hbi] = new int[nBinsBjets] ;
            N_0lep_input[mbi][hbi] = new int[nBinsBjets] ;
            N_1lep[mbi][hbi] = new int[nBinsBjets] ;
            N_ldp[mbi][hbi]  = new int[nBinsBjets] ;
            N_mc_ldp[mbi][hbi] = new float[nBinsBjets] ;
            susy_0lep[mbi][hbi] = new float[nBinsBjets] ;
            susy_1lep[mbi][hbi] = new float[nBinsBjets] ;
            susy_ldp[mbi][hbi]  = new float[nBinsBjets] ;
         }
      }

      int** N_Zee = new int*[nBinsMET] ;
      int** N_Zmm = new int*[nBinsMET] ;
      for ( int mbi=0; mbi<nBinsHT; mbi++ ) {
         N_Zee[mbi] = new int[nBinsHT] ;
         N_Zmm[mbi] = new int[nBinsHT] ;
      }

      float* acc_Zee = new float[nBinsMET] ;
      float* acc_Zmm = new float[nBinsMET] ;

      float Z_ee_eff(0.) ;
      float Z_mm_eff(0.) ;

      float* knn = new float[nBinsBjets] ;

      float Z_ee_pur(0.) ;
      float Z_mm_pur(0.) ;

      float* R_passfail = new float[nBinsHT] ;


     //============================================


      //--- Go through the file again and read in data needed to compute the expected values
      //    for the 0lep observables.

      if ( (infp=fopen( infilename,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infilename ) ;
         return ;
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
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return ;
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
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return ;
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
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return ;
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
               printf("\n\n *** bad input line.\n\n") ; return ;
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
               printf("\n\n *** bad input line.\n\n") ; return ;
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
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return ;
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
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return ;
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
                        mbi+1, hbi+1, bbi+1,  metbin, htbin, nbbin ) ; return ;
                  }
               } else {
                  printf("\n\n *** bad input line.\n\n") ; return ;
               }
            } // bbi.
         } // hbi
      } // mbi.


      //--- Read acc_Zee lines
      for ( int mbi=0; mbi<nBinsHT; mbi++ ) {
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
         int metbin ;
         int nmatch = sscanf( pname, "acc_Zee_M%d", &metbin ) ;
         if ( nmatch == 1 ) {
            acc_Zee[mbi] = pval ;
            printf(" Set acc_Zee for mbi=%d to %g\n", mbi, pval ) ;
         } else {
            printf("\n\n *** bad input line.\n\n") ; return ;
         }
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
      } // mbi


      //--- Read acc_Zmm lines
      for ( int mbi=0; mbi<nBinsHT; mbi++ ) {
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
         int metbin ;
         int nmatch = sscanf( pname, "acc_Zmm_M%d", &metbin ) ;
         if ( nmatch == 1 ) {
            acc_Zmm[mbi] = pval ;
            printf(" Set acc_Zmm for mbi=%d to %g\n", mbi, pval ) ;
         } else {
            printf("\n\n *** bad input line.\n\n") ; return ;
         }
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
      } // mbi


      //--- Z_ee_eff
      nread = fscanf( infp, "%s %g", pname, &pval ) ;
      if ( strcmp( pname, "Z_ee_eff" ) == 0 ) {
         Z_ee_eff = pval ;
         printf(" Set Z_ee_eff to %g\n", pval ) ;
      } else {
         printf("\n\n *** bad input line.\n\n") ; return ;
      }
      nread = fscanf( infp, "%s %g", pname, &pval ) ;


      //--- Z_mm_eff
      nread = fscanf( infp, "%s %g", pname, &pval ) ;
      if ( strcmp( pname, "Z_mm_eff" ) == 0 ) {
         Z_mm_eff = pval ;
         printf(" Set Z_mm_eff to %g\n", pval ) ;
      } else {
         printf("\n\n *** bad input line.\n\n") ; return ;
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
               printf("\n\n *** bad input line.\n\n") ; return ;
            }
         } else {
            printf("\n\n *** bad input line.\n\n") ; return ;
         }
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
      } // bbi.


      //--- Z_ee_pur
      nread = fscanf( infp, "%s %g", pname, &pval ) ;
      if ( strcmp( pname, "Z_ee_pur" ) == 0 ) {
         Z_ee_pur = pval ;
         printf(" Set Z_ee_pur to %g\n", pval ) ;
      } else {
         printf("\n\n *** bad input line.\n\n") ; return ;
      }
      nread = fscanf( infp, "%s %g", pname, &pval ) ;


      //--- Z_mm_pur
      nread = fscanf( infp, "%s %g", pname, &pval ) ;
      if ( strcmp( pname, "Z_mm_pur" ) == 0 ) {
         Z_mm_pur = pval ;
         printf(" Set Z_mm_pur to %g\n", pval ) ;
      } else {
         printf("\n\n *** bad input line.\n\n") ; return ;
      }
      nread = fscanf( infp, "%s %g", pname, &pval ) ;

      fclose(infp) ;


     //======================================================

      //--- Get SUSY inputs, if requested.

      if ( mgl>0. && mlsp>0. ) {

         ifstream insusystr ;
         insusystr.open( insusyfilename ) ;
         if ( !insusystr.good() ) {
            printf("\n\n *** Problem opening susy input file: %s\n\n", insusyfilename ) ; return ;
         }

         bool found(false) ;
         while ( insusystr.good() ) {

            float pointMgl ;
            float pointMlsp ;

            int nGen ;

            int ArraySize = 3 + 9*(nBinsMET*nBinsHT*nBinsBjets) ;
            double ArrayContent[ArraySize] ;

            for (int i = 0; infp && i < ArraySize; ++ i) {
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


            printf("\n\n  pointMgl = %g , pointMlsp = %g,  total 0lep events = %g\n", pointMgl, pointMlsp, input0lepSusyTotal ) ;
            printf("   reset susy 0lep total to %g\n\n", nsusy0leptarget ) ;

            float newTotal(0.) ;
            for (int i = 0 ; i < nBinsMET ; i++) {
              for (int j = 0 ; j < nBinsHT ; j++) {
                for (int k = 0 ; k < nBinsBjets ; k++) {
                   susy_0lep[i][j][k] = (nsusy0leptarget/input0lepSusyTotal) * susy_0lep[i][j][k] ;
                   susy_1lep[i][j][k] = (nsusy0leptarget/input0lepSusyTotal) * susy_1lep[i][j][k] ;
                   susy_ldp[i][j][k]  = (nsusy0leptarget/input0lepSusyTotal) * susy_ldp[i][j][k] ;
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
            printf("\n\n *** Didn't find susy point mgl=%g, mlsp=%g in input susy file %s\n\n", mgl, mlsp, insusyfilename ) ;
            return ;
         }

      } else {

            for (int i = 0 ; i < nBinsMET ; i++) {
              for (int j = 0 ; j < nBinsHT ; j++) {
                for (int k = 0 ; k < nBinsBjets ; k++) {
                   susy_0lep[i][j][k] = 0. ;
                   susy_1lep[i][j][k] = 0. ;
                   susy_ldp[i][j][k]  = 0. ;
                }
              }
            }

      }




     //======================================================

      //--- Compute the expected values for the 0lep observables, not including susy.
      float* Zee_factor = new float[nBinsBjets] ;
      float* Zmm_factor = new float[nBinsBjets] ;

      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
         Zee_factor[bbi] = ( 5.94 * Z_ee_pur * knn[bbi] ) / ( acc_Zee[0] * Z_ee_eff ) ; //*** ignoring MET dependence of acceptance.
         Zmm_factor[bbi] = ( 5.94 * Z_mm_pur * knn[bbi] ) / ( acc_Zmm[0] * Z_mm_eff ) ; //*** ignoring MET dependence of acceptance.
         printf("  nB=%d Zll scale factors:  ee=%g, mm=%g\n", bbi, Zee_factor[bbi], Zmm_factor[bbi] ) ;
      } // bbi.

      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         R_passfail[hbi] = 0.07 ; // guess for now.
      }

      float R_ttwj_0lep_over_1lep = 1.15 ; // guess for now.

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               float exp_0lep_ttwj = N_1lep[mbi][hbi][bbi] * R_ttwj_0lep_over_1lep ;

               float exp_0lep_qcd =  ( N_ldp[mbi][hbi][bbi] - N_mc_ldp[mbi][hbi][bbi] ) * R_passfail[hbi] ;
               if ( exp_0lep_qcd < 0. ) { exp_0lep_qcd = 0. ; }

               float exp_0lep_znn = 0.5 * ( N_Zee[mbi][hbi] * Zee_factor[bbi] + N_Zmm[mbi][hbi] * Zmm_factor[bbi] ) ;

               float exp_0lep = exp_0lep_ttwj + exp_0lep_qcd + exp_0lep_znn ;

               N_0lep[mbi][hbi][bbi] = (int) exp_0lep ;

               int Ndiff = N_0lep[mbi][hbi][bbi] - N_0lep_input[mbi][hbi][bbi] ;
               float diffPercent(0.) ;
               if ( N_0lep_input[mbi][hbi][bbi] > 0 ) { diffPercent = 100.* (Ndiff) / N_0lep_input[mbi][hbi][bbi] ; }

               printf(" m,h,b (%d,%d,%d): ttwj=%5.1f, qcd=%5.1f, Znn=%5.1f,  expected total=%5.1f,  N_0lep set to %4d, inupt val %4d  (diff %4d, %5.1f %%)\n",
                    mbi, hbi, bbi, exp_0lep_ttwj, exp_0lep_qcd, exp_0lep_znn, exp_0lep, N_0lep[mbi][hbi][bbi], N_0lep_input[mbi][hbi][bbi],
                    Ndiff, diffPercent ) ;

            } // bbi.
         } // hbi.
      } // mbi.

    //====================================

     //--- One more pass through the input file.  Write output file with 0lep observables set to expected
     //     values and susy added to 0lep, ldp, and 1lep observables, if requested.

      if ( (infp=fopen( infilename,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infilename ) ;
         return ;
      }

      FILE* outfp ;
      if ( (outfp=fopen( outfilename, "w"))==NULL ) {
         printf("\n\n *** Problem opening output file: %s\n\n", outfilename ) ;
         return ;
      }

      c = 0 ;
      while ( c!=10  ) {
         c = fgetc( infp ) ;
         fprintf(outfp, "%c", c ) ;
      }

      //--- N_0lep
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               nread = fscanf( infp, "%s %g", pname, &pval ) ;
               fprintf( outfp, "%s  %d\n", pname, (int)(N_0lep[mbi][hbi][bbi]+susy_0lep[mbi][hbi][bbi]) ) ;
            } // bbi.
         } // hbi.
      } // mbi.

      //--- N_1lep
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               nread = fscanf( infp, "%s %g", pname, &pval ) ;
               fprintf( outfp, "%s  %d\n", pname, (int)(N_1lep[mbi][hbi][bbi]+susy_1lep[mbi][hbi][bbi]) ) ;
            } // bbi.
         } // hbi.
      } // mbi.

      //--- N_ldp
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               nread = fscanf( infp, "%s %g", pname, &pval ) ;
               fprintf( outfp, "%s  %d\n", pname, (int)(N_ldp[mbi][hbi][bbi]+susy_ldp[mbi][hbi][bbi]) ) ;
            } // bbi.
         } // hbi.
      } // mbi.

      //--- the rest of the file stays the same.

      nread = 1 ;
      while ( nread>0 ) {
         nread = fscanf( infp, "%s %g", pname, &pval ) ;
         fprintf( outfp, "%s %g\n", pname , pval ) ;
      }

      fclose(infp) ;
      fclose(outfp) ;

      printf("\n\n Created output file %s\n\n", outfilename ) ;

   }










