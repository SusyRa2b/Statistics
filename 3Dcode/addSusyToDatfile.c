#include <iostream>
#include <fstream>
#include "getFileValue.c"
#include "updateFileValue.c"
#include "TSystem.h"


    void addSusyToDatfile( const char* indatfile,
                           const char* insusyfile,
                           float mgl = 1200,
                           float mlsp = 300,
                           float nsusy0lep = -1
                         ) {

      const int nBinsBjets = 3 ;
      const int nBinsMET = 4 ;
      const int nBinsHT = 4 ;

      printf("\n\n Opening SUSY scan input file : %s\n", insusyfile ) ;

      ifstream infp ;
      infp.open(insusyfile) ;
      if ( !infp.good() ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", insusyfile ) ;
         return ;
      }

      while ( infp.good() ) {

         int ArraySize = 3 + 9*nBinsMET*nBinsHT*nBinsBjets ;
         double ArrayContent[ArraySize] ;
         for (int i = 0; infp && i < ArraySize; ++ i) {
           infp >> ArrayContent[i];
         }

         float this_mgl  = ArrayContent[0] ;
         float this_mlsp = ArrayContent[1] ;
         if ( !( fabs(this_mgl-mgl)<1. && fabs(this_mlsp-mlsp)<1.) ) continue ;

         printf("\n\n Found mgl=%.0f, mlsp=%.0f\n\n", this_mgl, this_mlsp ) ;

         float n_0lep[nBinsMET][nBinsHT][nBinsBjets] ;
         float n_1lep[nBinsMET][nBinsHT][nBinsBjets] ;
         float n_ldp [nBinsMET][nBinsHT][nBinsBjets] ;
         float sum0lep(0.) ;
         for (int i = 0 ; i < nBinsMET ; i++) {
           for (int j = 0 ; j < nBinsHT ; j++) {
             for (int k = 0 ; k < nBinsBjets ; k++) {

               float n_0lep_raw  = ArrayContent[3 + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;
               float n_1lep_raw  = ArrayContent[4 + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;
               float n_ldp_raw = ArrayContent[5 + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;

               int nBins = nBinsMET*nBinsHT*nBinsBjets*3 ;

               float n_0lep_correction  = ArrayContent[3 + nBins + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;
               float n_1lep_correction  = ArrayContent[4 + nBins + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;
               float n_ldp_correction = ArrayContent[5 + nBins + i*(nBinsHT*nBinsBjets*3) + j*(nBinsBjets*3) + k*3] ;

               n_0lep[i][j][k] = n_0lep_raw * n_0lep_correction ;
               n_1lep[i][j][k] = n_1lep_raw * n_1lep_correction ;
               n_ldp [i][j][k] = n_ldp_raw  * n_ldp_correction  ;

               sum0lep += n_0lep[i][j][k] ;

             }
           }
         }

         printf("\n\n Total 0lep susy events, without rescaling: %.1f\n\n", sum0lep ) ;

         TString newfile( indatfile ) ;
         char newend[1000] ;
         sprintf( newend, "-wsusy-mgl%.0f-mlsp%.0f-%.0fevts.dat", mgl, mlsp, nsusy0lep ) ;
         newfile.ReplaceAll( ".dat", newend ) ;

         printf("\n\n Creating new dat file : %s\n\n", newfile.Data() ) ;

         char command[10000] ;
         sprintf( command, "cp %s %s", indatfile, newfile.Data() ) ;
         gSystem -> Exec( command ) ;

         float rescalefactor(1.) ;

         if ( nsusy0lep > 0 ) { rescalefactor = nsusy0lep / sum0lep ; }

         for (int i = 0 ; i < nBinsMET ; i++) {
           for (int j = 0 ; j < nBinsHT ; j++) {
             for (int k = 0 ; k < nBinsBjets ; k++) {

                float rs_0lep = rescalefactor * n_0lep[i][j][k] ;
                float rs_1lep = rescalefactor * n_1lep[i][j][k] ;
                float rs_ldp  = rescalefactor * n_ldp [i][j][k] ;

                printf(" m,h,b : %d,%d,%d :  0lep = %5.1f,  1lep = %5.1f,  LDP = %5.1f\n", i+1, j+1, k+1,
                     rs_0lep, rs_1lep, rs_ldp ) ;

                char pname[100] ;
                float oldval ;
                float newval ;

                sprintf( pname, "N_0lep_M%d_H%d_%db", i+1, j+1, k+1 ) ;
                getFileValue( newfile.Data(), pname, oldval ) ;
                newval = oldval + rs_0lep ;
                updateFileValue( newfile.Data(), pname, newval ) ;

                sprintf( pname, "N_1lep_M%d_H%d_%db", i+1, j+1, k+1 ) ;
                getFileValue( newfile.Data(), pname, oldval ) ;
                newval = oldval + rs_1lep ;
                updateFileValue( newfile.Data(), pname, newval ) ;

                sprintf( pname, "N_ldp_M%d_H%d_%db", i+1, j+1, k+1 ) ;
                getFileValue( newfile.Data(), pname, oldval ) ;
                newval = oldval + rs_ldp ;
                updateFileValue( newfile.Data(), pname, newval ) ;

                printf("\n") ;

             }
           }
         }

         return ;

      } // read while


    }

