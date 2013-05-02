
#include "frac_arrays.h"

#include "../getFileValue.h"



   void fill_frac_arrays( const char* frac_file = "../outputfiles/fitresults-ws-data-unblind-susyFixed0.0.txt"
                         ) {

      char tagstring[1000] ;
      int trow, tcol ;

      printf("\n\n") ;


   //--------


     //-- =2b

      printf("\n\n =2b, fractions\n\n") ;

      for ( int ci=0; ci<ncomps; ci++ ) {

         sprintf( tagstring, "frac_zl_%s_M4_H2_2b", comp_fracfile_name[ci] ) ; tcol = 0 ;
         getFileValueWithError( frac_file, tagstring, frac_2b_val[ci][tcol], frac_2b_err[ci][tcol] ) ;

         sprintf( tagstring, "frac_zl_%s_M4_H3_2b", comp_fracfile_name[ci] ) ; tcol = 1 ;
         getFileValueWithError( frac_file, tagstring, frac_2b_val[ci][tcol], frac_2b_err[ci][tcol] ) ;

         sprintf( tagstring, "frac_zl_%s_M4_H4_2b", comp_fracfile_name[ci] ) ; tcol = 2 ;
         getFileValueWithError( frac_file, tagstring, frac_2b_val[ci][tcol], frac_2b_err[ci][tcol] ) ;

         sprintf( tagstring, "htsum_zl_frac_%s_M4_2b", comp_fracfile_name[ci] ) ; tcol = 3 ;
         getFileValueWithError( frac_file, tagstring, frac_2b_val[ci][tcol], frac_2b_err[ci][tcol] ) ;

      } // ci.


      printf("\n\n") ;


     //-- >=3b

      printf("\n\n >=3b, fractions\n\n") ;

      for ( int ci=0; ci<ncomps; ci++ ) {


         sprintf( tagstring, "frac_zl_%s_M2_H1_3b", comp_fracfile_name[ci] ) ; trow = 0; tcol = 0 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "frac_zl_%s_M2_H2_3b", comp_fracfile_name[ci] ) ; trow = 0; tcol = 1 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "frac_zl_%s_M2_H3_3b", comp_fracfile_name[ci] ) ; trow = 0; tcol = 2 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "frac_zl_%s_M2_H4_3b", comp_fracfile_name[ci] ) ; trow = 0; tcol = 3 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "htsum_zl_frac_%s_M2_3b", comp_fracfile_name[ci] ) ; trow = 0; tcol = 4 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;


         sprintf( tagstring, "frac_zl_%s_M3_H1_3b", comp_fracfile_name[ci] ) ; trow = 1; tcol = 0 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "frac_zl_%s_M3_H2_3b", comp_fracfile_name[ci] ) ; trow = 1; tcol = 1 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "frac_zl_%s_M3_H3_3b", comp_fracfile_name[ci] ) ; trow = 1; tcol = 2 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "frac_zl_%s_M3_H4_3b", comp_fracfile_name[ci] ) ; trow = 1; tcol = 3 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "htsum_zl_frac_%s_M3_3b", comp_fracfile_name[ci] ) ; trow = 1; tcol = 4 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;


         frac_3b_val[ci][2][0] = 0. ; frac_3b_err[ci][2][0] = 0. ;

         sprintf( tagstring, "frac_zl_%s_M4_H2_3b", comp_fracfile_name[ci] ) ; trow = 2; tcol = 1 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "frac_zl_%s_M4_H3_3b", comp_fracfile_name[ci] ) ; trow = 2; tcol = 2 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "frac_zl_%s_M4_H4_3b", comp_fracfile_name[ci] ) ; trow = 2; tcol = 3 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "htsum_zl_frac_%s_M4_3b", comp_fracfile_name[ci] ) ; trow = 2; tcol = 4 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;


         sprintf( tagstring, "metsum_zl_frac_%s_H1_3b", comp_fracfile_name[ci] ) ; trow = 3; tcol = 0 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "metsum_zl_frac_%s_H2_3b", comp_fracfile_name[ci] ) ; trow = 3; tcol = 1 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "metsum_zl_frac_%s_H3_3b", comp_fracfile_name[ci] ) ; trow = 3; tcol = 2 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "metsum_zl_frac_%s_H4_3b", comp_fracfile_name[ci] ) ; trow = 3; tcol = 3 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;

         sprintf( tagstring, "methtsum_zl_frac_%s_3b", comp_fracfile_name[ci] ) ; trow = 3; tcol = 4 ;
         getFileValueWithError( frac_file, tagstring, frac_3b_val[ci][trow][tcol], frac_3b_err[ci][trow][tcol] ) ;



      } // ci.




      printf("\n\n\n\n\n +++ =2b ++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n") ;

      for ( int ci=0; ci<ncomps; ci++ ) {

         printf("\n\n ==== %s ===============================================\n\n", comp_out_name[ci] ) ;

         for ( int coli=0; coli<4; coli++ ) {
            printf(" %5.2f +/- %5.2f  |  ", frac_2b_val[ci][coli], frac_2b_err[ci][coli] ) ;
         } // coli.
         printf("\n") ;

      } // ci.




      printf("\n\n\n\n\n +++ >=3b ++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n") ;

      for ( int ci=0; ci<ncomps; ci++ ) {

         printf("\n\n ==== %s ===============================================\n\n", comp_out_name[ci] ) ;

         for ( int rowi=0; rowi<4; rowi++ ) {
            for ( int coli=0; coli<5; coli++ ) {
               printf(" %5.2f +/- %5.2f  |  ", frac_3b_val[ci][rowi][coli], frac_3b_err[ci][rowi][coli] ) ;
            } // coli.
            printf("\n") ;
         } // rowi.

      } // ci.



      printf("\n\n") ;


   }

