
//    float ub_3b_val[4][5] ;
//    float ub_3b_p1s[4][5] ;
//    float ub_3b_m1s[4][5] ;
//    float ub_3b_p2s[4][5] ;
//    float ub_3b_m2s[4][5] ;

//    float ff_3b_val[4][5] ;
//    float ff_3b_p1s[4][5] ;
//    float ff_3b_m1s[4][5] ;
//    float ff_3b_p2s[4][5] ;
//    float ff_3b_m2s[4][5] ;

#include "fill_arrays.c"


   void generate_table1( const char* scan_results_file = "scan-results.txt" ) {


      fill_arrays( scan_results_file ) ;


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




   }

