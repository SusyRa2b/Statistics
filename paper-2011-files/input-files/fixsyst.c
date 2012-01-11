
#include <stdio.h>

  void fixsyst( const char* allsystfname = "likelihood-newfit-oldsusysyst-LM9-HT400-SIGMET250.txt" ,
                const char* btageffsystfname = "likelihood-newfit-syst-btageff-LM9-HT400-SIGMET250.txt",
                const char* outfname = "likelihood-newfit-syst-allbutbtageff-LM9-HT400-SIGMET250.txt"
              ) {

      FILE* allsystfile(0x0) ;
      FILE* btageffsystfile(0x0) ;

      FILE* outfile(0x0) ;

      if ( (allsystfile = fopen( allsystfname, "r" ))==0x0 ) {
         printf("\n\n\n *** problem opening allsystfile : %s\n\n", allsystfname ) ;
         return ;
      }

      if ( (btageffsystfile = fopen( btageffsystfname, "r" ))==0x0 ) {
         printf("\n\n\n *** problem opening btageffsystfile : %s\n\n", btageffsystfname ) ;
         return ;
      }

      if ( (outfile = fopen( outfname, "w" ))==0x0 ) {
         printf("\n\n\n *** problem opening outfile : %s\n\n", outfname ) ;
         return ;
      }

      while ( !feof( allsystfile ) ) {

         float allsyst_m1, allsyst_m2 ;
         float allsyst_sig_1b, allsyst_sb_1b, allsyst_sig_sl_1b, allsyst_sb_sl_1b, allsyst_sig_ldp_1b, allsyst_sb_ldp_1b ;
         float allsyst_sig_2b, allsyst_sb_2b, allsyst_sig_sl_2b, allsyst_sb_sl_2b, allsyst_sig_ldp_2b, allsyst_sb_ldp_2b ;
         float allsyst_sig_3b, allsyst_sb_3b, allsyst_sig_sl_3b, allsyst_sb_sl_3b, allsyst_sig_ldp_3b, allsyst_sb_ldp_3b ;

         fscanf( allsystfile, "%f %f   %f %f %f %f %f %f    %f %f %f %f %f %f    %f %f %f %f %f %f     ",
                 &allsyst_m1, &allsyst_m2,
                &allsyst_sig_1b, &allsyst_sb_1b, &allsyst_sig_sl_1b, &allsyst_sb_sl_1b, &allsyst_sig_ldp_1b, &allsyst_sb_ldp_1b,
                &allsyst_sig_2b, &allsyst_sb_2b, &allsyst_sig_sl_2b, &allsyst_sb_sl_2b, &allsyst_sig_ldp_2b, &allsyst_sb_ldp_2b,
                &allsyst_sig_3b, &allsyst_sb_3b, &allsyst_sig_sl_3b, &allsyst_sb_sl_3b, &allsyst_sig_ldp_3b, &allsyst_sb_ldp_3b
               ) ;


         float btageffsyst_m1, btageffsyst_m2 ;
         float btageffsyst_sig_1b, btageffsyst_sb_1b, btageffsyst_sig_sl_1b, btageffsyst_sb_sl_1b, btageffsyst_sig_ldp_1b, btageffsyst_sb_ldp_1b ;
         float btageffsyst_sig_2b, btageffsyst_sb_2b, btageffsyst_sig_sl_2b, btageffsyst_sb_sl_2b, btageffsyst_sig_ldp_2b, btageffsyst_sb_ldp_2b ;
         float btageffsyst_sig_3b, btageffsyst_sb_3b, btageffsyst_sig_sl_3b, btageffsyst_sb_sl_3b, btageffsyst_sig_ldp_3b, btageffsyst_sb_ldp_3b ;

         fscanf( btageffsystfile, "%f %f   %f %f %f %f %f %f    %f %f %f %f %f %f    %f %f %f %f %f %f     ",
                 &btageffsyst_m1, &btageffsyst_m2,
                &btageffsyst_sig_1b, &btageffsyst_sb_1b, &btageffsyst_sig_sl_1b, &btageffsyst_sb_sl_1b, &btageffsyst_sig_ldp_1b, &btageffsyst_sb_ldp_1b,
                &btageffsyst_sig_2b, &btageffsyst_sb_2b, &btageffsyst_sig_sl_2b, &btageffsyst_sb_sl_2b, &btageffsyst_sig_ldp_2b, &btageffsyst_sb_ldp_2b,
                &btageffsyst_sig_3b, &btageffsyst_sb_3b, &btageffsyst_sig_sl_3b, &btageffsyst_sb_sl_3b, &btageffsyst_sig_ldp_3b, &btageffsyst_sb_ldp_3b
               ) ;




         printf("\n allsyst m1=%.1f, m2=%.1f\n", allsyst_m1, allsyst_m2 ) ;

         if ( btageffsyst_m1!= allsyst_m1 || btageffsyst_m2 != allsyst_m2 ) {
            printf(" inconsistent m1 and m2 in btageffsyst file: %.1f, %.1f\n\n", btageffsyst_m1, btageffsyst_m2 ) ;
         }

         float allbutsyst_m1, allbutsyst_m2 ;
         float allbutsyst_sig_1b, allbutsyst_sb_1b, allbutsyst_sig_sl_1b, allbutsyst_sb_sl_1b, allbutsyst_sig_ldp_1b, allbutsyst_sb_ldp_1b ;
         float allbutsyst_sig_2b, allbutsyst_sb_2b, allbutsyst_sig_sl_2b, allbutsyst_sb_sl_2b, allbutsyst_sig_ldp_2b, allbutsyst_sb_ldp_2b ;
         float allbutsyst_sig_3b, allbutsyst_sb_3b, allbutsyst_sig_sl_3b, allbutsyst_sb_sl_3b, allbutsyst_sig_ldp_3b, allbutsyst_sb_ldp_3b ;

         allbutsyst_m1 = allsyst_m1 ;
         allbutsyst_m2 = allsyst_m2 ;

         allbutsyst_sig_1b     = sqrt( allsyst_sig_1b*allsyst_sig_1b - 100*100*btageffsyst_sig_1b*btageffsyst_sig_1b ) ;
         allbutsyst_sb_1b      = sqrt( allsyst_sb_1b*allsyst_sb_1b - 100*100*btageffsyst_sb_1b*btageffsyst_sb_1b ) ;
         allbutsyst_sig_sl_1b  = sqrt( allsyst_sig_sl_1b*allsyst_sig_sl_1b - 100*100*btageffsyst_sig_sl_1b*btageffsyst_sig_sl_1b ) ;
         allbutsyst_sb_sl_1b   = sqrt( allsyst_sb_sl_1b*allsyst_sb_sl_1b - 100*100*btageffsyst_sb_sl_1b*btageffsyst_sb_sl_1b ) ;
         allbutsyst_sig_ldp_1b = sqrt( allsyst_sig_ldp_1b*allsyst_sig_ldp_1b - 100*100*btageffsyst_sig_ldp_1b*btageffsyst_sig_ldp_1b ) ;
         allbutsyst_sb_ldp_1b  = sqrt( allsyst_sb_ldp_1b*allsyst_sb_ldp_1b - 100*100*btageffsyst_sb_ldp_1b*btageffsyst_sb_ldp_1b ) ;

         allbutsyst_sig_2b     = sqrt( allsyst_sig_2b*allsyst_sig_2b - 100*100*btageffsyst_sig_2b*btageffsyst_sig_2b ) ;
         allbutsyst_sb_2b      = sqrt( allsyst_sb_2b*allsyst_sb_2b - 100*100*btageffsyst_sb_2b*btageffsyst_sb_2b ) ;
         allbutsyst_sig_sl_2b  = sqrt( allsyst_sig_sl_2b*allsyst_sig_sl_2b - 100*100*btageffsyst_sig_sl_2b*btageffsyst_sig_sl_2b ) ;
         allbutsyst_sb_sl_2b   = sqrt( allsyst_sb_sl_2b*allsyst_sb_sl_2b - 100*100*btageffsyst_sb_sl_2b*btageffsyst_sb_sl_2b ) ;
         allbutsyst_sig_ldp_2b = sqrt( allsyst_sig_ldp_2b*allsyst_sig_ldp_2b - 100*100*btageffsyst_sig_ldp_2b*btageffsyst_sig_ldp_2b ) ;
         allbutsyst_sb_ldp_2b  = sqrt( allsyst_sb_ldp_2b*allsyst_sb_ldp_2b - 100*100*btageffsyst_sb_ldp_2b*btageffsyst_sb_ldp_2b ) ;

         allbutsyst_sig_3b     = sqrt( allsyst_sig_3b*allsyst_sig_3b - 100*100*btageffsyst_sig_3b*btageffsyst_sig_3b ) ;
         allbutsyst_sb_3b      = sqrt( allsyst_sb_3b*allsyst_sb_3b - 100*100*btageffsyst_sb_3b*btageffsyst_sb_3b ) ;
         allbutsyst_sig_sl_3b  = sqrt( allsyst_sig_sl_3b*allsyst_sig_sl_3b - 100*100*btageffsyst_sig_sl_3b*btageffsyst_sig_sl_3b ) ;
         allbutsyst_sb_sl_3b   = sqrt( allsyst_sb_sl_3b*allsyst_sb_sl_3b - 100*100*btageffsyst_sb_sl_3b*btageffsyst_sb_sl_3b ) ;
         allbutsyst_sig_ldp_3b = sqrt( allsyst_sig_ldp_3b*allsyst_sig_ldp_3b - 100*100*btageffsyst_sig_ldp_3b*btageffsyst_sig_ldp_3b ) ;
         allbutsyst_sb_ldp_3b  = sqrt( allsyst_sb_ldp_3b*allsyst_sb_ldp_3b - 100*100*btageffsyst_sb_ldp_3b*btageffsyst_sb_ldp_3b ) ;

         printf("\n") ;
         printf(" sig_1b     : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sig_1b, 100*btageffsyst_sig_1b, allbutsyst_sig_1b ) ;
         printf(" sb_1b      : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sb_1b, 100*btageffsyst_sb_1b, allbutsyst_sb_1b ) ;
         printf(" sig_sl_1b  : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sig_sl_1b, 100*btageffsyst_sig_sl_1b, allbutsyst_sig_sl_1b ) ;
         printf(" sb_sl_1b   : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sb_sl_1b, 100*btageffsyst_sb_sl_1b, allbutsyst_sb_sl_1b ) ;
         printf(" sig_ldp_1b : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sig_ldp_1b, 100*btageffsyst_sig_ldp_1b, allbutsyst_sig_ldp_1b ) ;
         printf(" sb_ldp_1b  : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sb_ldp_1b, 100*btageffsyst_sb_ldp_1b, allbutsyst_sb_ldp_1b ) ;

         printf("\n") ;
         printf(" sig_2b     : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sig_2b, 100*btageffsyst_sig_2b, allbutsyst_sig_2b ) ;
         printf(" sb_2b      : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sb_2b, 100*btageffsyst_sb_2b, allbutsyst_sb_2b ) ;
         printf(" sig_sl_2b  : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sig_sl_2b, 100*btageffsyst_sig_sl_2b, allbutsyst_sig_sl_2b ) ;
         printf(" sb_sl_2b   : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sb_sl_2b, 100*btageffsyst_sb_sl_2b, allbutsyst_sb_sl_2b ) ;
         printf(" sig_ldp_2b : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sig_ldp_2b, 100*btageffsyst_sig_ldp_2b, allbutsyst_sig_ldp_2b ) ;
         printf(" sb_ldp_2b  : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sb_ldp_2b, 100*btageffsyst_sb_ldp_2b, allbutsyst_sb_ldp_2b ) ;

         printf("\n") ;
         printf(" sig_3b     : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sig_3b, 100*btageffsyst_sig_3b, allbutsyst_sig_3b ) ;
         printf(" sb_3b      : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sb_3b, 100*btageffsyst_sb_3b, allbutsyst_sb_3b ) ;
         printf(" sig_sl_3b  : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sig_sl_3b, 100*btageffsyst_sig_sl_3b, allbutsyst_sig_sl_3b ) ;
         printf(" sb_sl_3b   : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sb_sl_3b, 100*btageffsyst_sb_sl_3b, allbutsyst_sb_sl_3b ) ;
         printf(" sig_ldp_3b : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sig_ldp_3b, 100*btageffsyst_sig_ldp_3b, allbutsyst_sig_ldp_3b ) ;
         printf(" sb_ldp_3b  : all = %6.3f,  btageff = %6.3f,  all-but-btageff = %6.3f\n", allsyst_sb_ldp_3b, 100*btageffsyst_sb_ldp_3b, allbutsyst_sb_ldp_3b ) ;

         fprintf( outfile, "%.0f %.0f   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f    %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f    %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f",
                 allbutsyst_m1, allbutsyst_m2,
                allbutsyst_sig_1b, allbutsyst_sb_1b, allbutsyst_sig_sl_1b, allbutsyst_sb_sl_1b, allbutsyst_sig_ldp_1b, allbutsyst_sb_ldp_1b,
                allbutsyst_sig_2b, allbutsyst_sb_2b, allbutsyst_sig_sl_2b, allbutsyst_sb_sl_2b, allbutsyst_sig_ldp_2b, allbutsyst_sb_ldp_2b,
                allbutsyst_sig_3b, allbutsyst_sb_3b, allbutsyst_sig_sl_3b, allbutsyst_sb_sl_3b, allbutsyst_sig_ldp_3b, allbutsyst_sb_ldp_3b
               ) ;

         if ( feof( allsystfile ) ) break ;

      }

      fclose( outfile ) ;

  }

