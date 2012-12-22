
#include "arrays.h"

#include "getScanDataLine.c"
#include "../getFileValue.c"


   void fill_arrays( const char* scan_results_file = "scan-results.txt", const char* datfile = "../datfiles_18fb/data-vals-unblind.dat" ) {

      char tagstring[1000] ;
      int trow, tcol ;

      printf("\n\n") ;





     //--- =2b, unbiased fit.

      sprintf( tagstring, "scan-hb-n_M4_H2_2b.log" ) ; tcol = 0 ;
      getScanDataLine( scan_results_file, tagstring, ub_2b_val[tcol], ub_2b_p1s[tcol], ub_2b_m1s[tcol], ub_2b_p2s[tcol], ub_2b_m2s[tcol] ) ;

      sprintf( tagstring, "scan-hb-n_M4_H3_2b.log" ) ; tcol = 1 ;
      getScanDataLine( scan_results_file, tagstring, ub_2b_val[tcol], ub_2b_p1s[tcol], ub_2b_m1s[tcol], ub_2b_p2s[tcol], ub_2b_m2s[tcol] ) ;

      sprintf( tagstring, "scan-hb-n_M4_H4_2b.log" ) ; tcol = 2 ;
      getScanDataLine( scan_results_file, tagstring, ub_2b_val[tcol], ub_2b_p1s[tcol], ub_2b_m1s[tcol], ub_2b_p2s[tcol], ub_2b_m2s[tcol] ) ;

      sprintf( tagstring, "scan-hb-bs-n_M4_H234_2b.log" ) ; tcol = 3 ;
      getScanDataLine( scan_results_file, tagstring, ub_2b_val[tcol], ub_2b_p1s[tcol], ub_2b_m1s[tcol], ub_2b_p2s[tcol], ub_2b_m2s[tcol] ) ;



      printf("\n\n") ;

      for ( int ci=0; ci<4; ci++ ) {
         if ( ub_2b_val[ci] > 50 ) {
            printf("  %5.0f ", ub_2b_val[ci] ) ;
         } else {
            printf("  %5.1f ", ub_2b_val[ci] ) ;
         }
      } // ci.
      printf("\n") ;

      printf("\n\n") ;





     //--- =2b, full fit.

      sprintf( tagstring, "scan-ff-susyfixed-n_M4_H2_2b.log" ) ; tcol = 0 ;
      getScanDataLine( scan_results_file, tagstring, ff_2b_val[tcol], ff_2b_p1s[tcol], ff_2b_m1s[tcol], ff_2b_p2s[tcol], ff_2b_m2s[tcol] ) ;

      sprintf( tagstring, "scan-ff-susyfixed-n_M4_H3_2b.log" ) ; tcol = 1 ;
      getScanDataLine( scan_results_file, tagstring, ff_2b_val[tcol], ff_2b_p1s[tcol], ff_2b_m1s[tcol], ff_2b_p2s[tcol], ff_2b_m2s[tcol] ) ;

      sprintf( tagstring, "scan-ff-susyfixed-n_M4_H4_2b.log" ) ; tcol = 2 ;
      getScanDataLine( scan_results_file, tagstring, ff_2b_val[tcol], ff_2b_p1s[tcol], ff_2b_m1s[tcol], ff_2b_p2s[tcol], ff_2b_m2s[tcol] ) ;

      sprintf( tagstring, "scan-ff-bs-n_M4_H234_2b.log" ) ; tcol = 3 ;
      getScanDataLine( scan_results_file, tagstring, ff_2b_val[tcol], ff_2b_p1s[tcol], ff_2b_m1s[tcol], ff_2b_p2s[tcol], ff_2b_m2s[tcol] ) ;



      printf("\n\n") ;

      for ( int ci=0; ci<4; ci++ ) {
         if ( ff_2b_val[ci] > 50 ) {
            printf("  %5.0f ", ff_2b_val[ci] ) ;
         } else {
            printf("  %5.1f ", ff_2b_val[ci] ) ;
         }
      } // ci.
      printf("\n") ;

      printf("\n\n") ;




     //-- =2b, data obs

      sprintf( tagstring, "N_0lep_M4_H2_2b" ) ; tcol = 0 ;
      getFileValue( datfile, tagstring, data_2b_obs[tcol] ) ;

      sprintf( tagstring, "N_0lep_M4_H3_2b" ) ; tcol = 1 ;
      getFileValue( datfile, tagstring, data_2b_obs[tcol] ) ;

      sprintf( tagstring, "N_0lep_M4_H4_2b" ) ; tcol = 2 ;
      getFileValue( datfile, tagstring, data_2b_obs[tcol] ) ;

      data_2b_obs[3] = data_2b_obs[0] + data_2b_obs[1] + data_2b_obs[2] ;


      printf("\n\n") ;

      for ( int ci=0; ci<4; ci++ ) {
         printf("  %5.0f ", data_2b_obs[ci] ) ;
      } // ci.
      printf("\n") ;

      printf("\n\n") ;





     //--- >=3b, unbiased fit.

      sprintf( tagstring, "scan-hb-n_M2_H1_3b.log" ) ; trow = 0 ; tcol = 0 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-n_M2_H2_3b.log" ) ; trow = 0 ; tcol = 1 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-n_M2_H3_3b.log" ) ; trow = 0 ; tcol = 2 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-n_M2_H4_3b.log" ) ; trow = 0 ; tcol = 3 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-bs-n_M2_H1234_3b.log" ) ; trow = 0 ; tcol = 4 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;



      sprintf( tagstring, "scan-hb-n_M3_H1_3b.log" ) ; trow = 1 ; tcol = 0 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-n_M3_H2_3b.log" ) ; trow = 1 ; tcol = 1 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-n_M3_H3_3b.log" ) ; trow = 1 ; tcol = 2 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-n_M3_H4_3b.log" ) ; trow = 1 ; tcol = 3 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-bs-n_M3_H1234_3b.log" ) ; trow = 1 ; tcol = 4 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;




      sprintf( tagstring, "scan-hb-n_M4_H2_3b.log" ) ; trow = 2 ; tcol = 1 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-n_M4_H3_3b.log" ) ; trow = 2 ; tcol = 2 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-n_M4_H4_3b.log" ) ; trow = 2 ; tcol = 3 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-bs-n_M4_H234_3b.log" ) ; trow = 2 ; tcol = 4 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;



      sprintf( tagstring, "scan-hb-bs-n_M23_H1_3b.log" ) ; trow = 3 ; tcol = 0 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-bs-n_M234_H2_3b.log" ) ; trow = 3 ; tcol = 1 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-bs-n_M234_H3_3b.log" ) ; trow = 3 ; tcol = 2 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-hb-bs-n_M234_H4_3b.log" ) ; trow = 3 ; tcol = 3 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;



      sprintf( tagstring, "scan-hb-bs-n_M234_H1234_3b.log" ) ; trow = 3 ; tcol = 4 ;
      getScanDataLine( scan_results_file, tagstring, ub_3b_val[trow][tcol], ub_3b_p1s[trow][tcol], ub_3b_m1s[trow][tcol], ub_3b_p2s[trow][tcol], ub_3b_m2s[trow][tcol] ) ;


      printf("\n\n") ;

      for ( int ri=0; ri<4; ri++ ) {
         for ( int ci=0; ci<5; ci++ ) {
            if ( ub_3b_val[ri][ci] > 50 ) {
               printf("  %5.0f ", ub_3b_val[ri][ci] ) ;
            } else {
               printf("  %5.1f ", ub_3b_val[ri][ci] ) ;
            }
         } // ci.
         printf("\n") ;
      } // ri

      printf("\n\n") ;





     //--- >=3b, full fit.

      sprintf( tagstring, "scan-ff-susyfixed-n_M2_H1_3b.log" ) ; trow = 0 ; tcol = 0 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-susyfixed-n_M2_H2_3b.log" ) ; trow = 0 ; tcol = 1 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-susyfixed-n_M2_H3_3b.log" ) ; trow = 0 ; tcol = 2 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-susyfixed-n_M2_H4_3b.log" ) ; trow = 0 ; tcol = 3 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-bs-n_M2_H1234_3b.log" ) ; trow = 0 ; tcol = 4 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;



      sprintf( tagstring, "scan-ff-susyfixed-n_M3_H1_3b.log" ) ; trow = 1 ; tcol = 0 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-susyfixed-n_M3_H2_3b.log" ) ; trow = 1 ; tcol = 1 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-susyfixed-n_M3_H3_3b.log" ) ; trow = 1 ; tcol = 2 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-susyfixed-n_M3_H4_3b.log" ) ; trow = 1 ; tcol = 3 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-bs-n_M3_H1234_3b.log" ) ; trow = 1 ; tcol = 4 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;




      sprintf( tagstring, "scan-ff-susyfixed-n_M4_H2_3b.log" ) ; trow = 2 ; tcol = 1 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-susyfixed-n_M4_H3_3b.log" ) ; trow = 2 ; tcol = 2 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-susyfixed-n_M4_H4_3b.log" ) ; trow = 2 ; tcol = 3 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-bs-n_M4_H234_3b.log" ) ; trow = 2 ; tcol = 4 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;



      sprintf( tagstring, "scan-ff-bs-n_M23_H1_3b.log" ) ; trow = 3 ; tcol = 0 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-bs-n_M234_H2_3b.log" ) ; trow = 3 ; tcol = 1 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-bs-n_M234_H3_3b.log" ) ; trow = 3 ; tcol = 2 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;

      sprintf( tagstring, "scan-ff-bs-n_M234_H4_3b.log" ) ; trow = 3 ; tcol = 3 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;



      sprintf( tagstring, "scan-ff-bs-n_M234_H1234_3b.log" ) ; trow = 3 ; tcol = 4 ;
      getScanDataLine( scan_results_file, tagstring, ff_3b_val[trow][tcol], ff_3b_p1s[trow][tcol], ff_3b_m1s[trow][tcol], ff_3b_p2s[trow][tcol], ff_3b_m2s[trow][tcol] ) ;


      printf("\n\n") ;

      for ( int ri=0; ri<4; ri++ ) {
         for ( int ci=0; ci<5; ci++ ) {
            if ( ff_3b_val[ri][ci] > 50 ) {
               printf("  %5.0f ", ff_3b_val[ri][ci] ) ;
            } else {
               printf("  %5.1f ", ff_3b_val[ri][ci] ) ;
            }
         } // ci.
         printf("\n") ;
      } // ri

      printf("\n\n") ;






     //-- >=3b, data obs

      sprintf( tagstring, "N_0lep_M2_H1_3b" ) ; trow = 0 ; tcol = 0 ;
      getFileValue( datfile, tagstring, data_3b_obs[trow][tcol] ) ;

      sprintf( tagstring, "N_0lep_M2_H2_3b" ) ; trow = 0 ; tcol = 1 ;
      getFileValue( datfile, tagstring, data_3b_obs[trow][tcol] ) ;

      sprintf( tagstring, "N_0lep_M2_H3_3b" ) ; trow = 0 ; tcol = 2 ;
      getFileValue( datfile, tagstring, data_3b_obs[trow][tcol] ) ;

      sprintf( tagstring, "N_0lep_M2_H4_3b" ) ; trow = 0 ; tcol = 3 ;
      getFileValue( datfile, tagstring, data_3b_obs[trow][tcol] ) ;

      data_3b_obs[0][4] = data_3b_obs[0][0] + data_3b_obs[0][1] + data_3b_obs[0][2] + data_3b_obs[0][3] ;


      sprintf( tagstring, "N_0lep_M3_H1_3b" ) ; trow = 1 ; tcol = 0 ;
      getFileValue( datfile, tagstring, data_3b_obs[trow][tcol] ) ;

      sprintf( tagstring, "N_0lep_M3_H2_3b" ) ; trow = 1 ; tcol = 1 ;
      getFileValue( datfile, tagstring, data_3b_obs[trow][tcol] ) ;

      sprintf( tagstring, "N_0lep_M3_H3_3b" ) ; trow = 1 ; tcol = 2 ;
      getFileValue( datfile, tagstring, data_3b_obs[trow][tcol] ) ;

      sprintf( tagstring, "N_0lep_M3_H4_3b" ) ; trow = 1 ; tcol = 3 ;
      getFileValue( datfile, tagstring, data_3b_obs[trow][tcol] ) ;

      data_3b_obs[1][4] = data_3b_obs[1][0] + data_3b_obs[1][1] + data_3b_obs[1][2] + data_3b_obs[1][3] ;



      sprintf( tagstring, "N_0lep_M4_H2_3b" ) ; trow = 2 ; tcol = 1 ;
      getFileValue( datfile, tagstring, data_3b_obs[trow][tcol] ) ;

      sprintf( tagstring, "N_0lep_M4_H3_3b" ) ; trow = 2 ; tcol = 2 ;
      getFileValue( datfile, tagstring, data_3b_obs[trow][tcol] ) ;

      sprintf( tagstring, "N_0lep_M4_H4_3b" ) ; trow = 2 ; tcol = 3 ;
      getFileValue( datfile, tagstring, data_3b_obs[trow][tcol] ) ;

      data_3b_obs[2][4] = data_3b_obs[2][1] + data_3b_obs[2][2] + data_3b_obs[2][3] ;

      data_3b_obs[3][0] = data_3b_obs[0][0] + data_3b_obs[1][0] ;
      data_3b_obs[3][1] = data_3b_obs[0][1] + data_3b_obs[1][1] + data_3b_obs[2][1] ;
      data_3b_obs[3][2] = data_3b_obs[0][2] + data_3b_obs[1][2] + data_3b_obs[2][2] ;
      data_3b_obs[3][3] = data_3b_obs[0][3] + data_3b_obs[1][3] + data_3b_obs[2][3] ;
      data_3b_obs[3][4] = data_3b_obs[0][4] + data_3b_obs[1][4] + data_3b_obs[2][4] ;



      printf("\n\n") ;

      for ( int ri=0; ri<4; ri++ ) {
         for ( int ci=0; ci<5; ci++ ) {
            printf("  %5.0f ", data_3b_obs[ri][ci] ) ;
         } // ci.
         printf("\n") ;
      } // ri

      printf("\n\n") ;




   }

