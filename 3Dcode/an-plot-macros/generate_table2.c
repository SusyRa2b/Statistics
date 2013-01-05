

#include "fill_frac_arrays.c"

#include <stdio.h>

   void fitfrac_3b_table( FILE* outfp, const char* table_title ) ;
   void fitfrac_2b_table( FILE* outfp, const char* table_title ) ;

   //---------------

   void generate_table2( const char* frac_file_unbiased = "../outputfiles/fitresults-halfblind-ws-data-unblind-susyFixed0.0.txt",
                         const char* frac_file_fullfit  = "../outputfiles/fitresults-ws-data-unblind-susyFixed0.0.txt"
                         ) {




      FILE* outfp ;
      if ( (outfp=fopen( "tables2.tex", "w" ))==NULL ) {
         printf("\n\n *** can't open output file.\n\n") ; return ;
      }


      fprintf( outfp, "\\documentclass[11pt]{article}\n" ) ;
      fprintf( outfp, "\\def\\znunu {$Z \\to \\nu \\bar{\\nu}$\\ }\n") ;
      fprintf( outfp, "\\def\\vsmtvs {\\rule[-0.3cm]{0cm}{0.75cm}}\n" ) ;
      fprintf( outfp, "\\begin{document}\n" ) ;


      fprintf( outfp, "\n\n\n" ) ;

      fill_frac_arrays( frac_file_unbiased ) ;

      fitfrac_2b_table(  outfp, "Unbiased fit, SM-only" ) ;
      fitfrac_3b_table(  outfp, "Unbiased fit, SM-only" ) ;

      fprintf( outfp, "\n\n\n" ) ;

      fprintf( outfp, "    \\pagebreak\n") ;

      fprintf( outfp, "\n\n\n" ) ;


      fill_frac_arrays( frac_file_fullfit ) ;

      fitfrac_2b_table(  outfp, "Full fit, SM-only" ) ;
      fitfrac_3b_table(  outfp, "Full fit, SM-only" ) ;


      fprintf( outfp, "\\end{document}\n" ) ;

      fclose( outfp ) ;

   }


  //================================================================================================================================

   void fitfrac_3b_table( FILE* outfp, const char* table_title ) {

      fprintf( outfp, "\n\n\n" ) ;
      fprintf( outfp, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n") ;
      fprintf( outfp, "    \\begin{tabular}{|c||c|c|c|c||c|}\n" ) ;

      for ( int ci=0; ci<3; ci++ ) {

         fprintf( outfp, "    \\hline\n" ) ;
         fprintf( outfp, "    \\hline\n" ) ;
         fprintf( outfp, "    \\multicolumn{6}{|c|}{ %s fraction, %s \\vsmtvs } \\\\ \n", comp_out_name[ci], table_title ) ;
         fprintf( outfp, "    \\hline\n" ) ;
         fprintf( outfp, "      & HT1 & HT2 & HT3 & HT4 & HT1-4 \\\\ \n" ) ;
         fprintf( outfp, "    \\hline\n" ) ;
         fprintf( outfp, "    \\hline\n" ) ;
         fprintf( outfp, "       $N_b\\ge3$, MET2   &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f     \\vsmtvs \\\\ \n",
           frac_3b_val[ci][0][0], frac_3b_err[ci][0][0],
           frac_3b_val[ci][0][1], frac_3b_err[ci][0][1],
           frac_3b_val[ci][0][2], frac_3b_err[ci][0][2],
           frac_3b_val[ci][0][3], frac_3b_err[ci][0][3],
           frac_3b_val[ci][0][4], frac_3b_err[ci][0][4]
           ) ;
         fprintf( outfp, "       $N_b\\ge3$, MET3   &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f     \\vsmtvs \\\\ \n",
           frac_3b_val[ci][1][0], frac_3b_err[ci][1][0],
           frac_3b_val[ci][1][1], frac_3b_err[ci][1][1],
           frac_3b_val[ci][1][2], frac_3b_err[ci][1][2],
           frac_3b_val[ci][1][3], frac_3b_err[ci][1][3],
           frac_3b_val[ci][1][4], frac_3b_err[ci][1][4]
           ) ;
         fprintf( outfp, "       $N_b\\ge3$, MET4   &     ---        &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f     \\vsmtvs \\\\ \n",
           frac_3b_val[ci][2][1], frac_3b_err[ci][2][1],
           frac_3b_val[ci][2][2], frac_3b_err[ci][2][2],
           frac_3b_val[ci][2][3], frac_3b_err[ci][2][3],
           frac_3b_val[ci][2][4], frac_3b_err[ci][2][4]
           ) ;
         fprintf( outfp, "    \\hline\n" ) ;
         fprintf( outfp, "    \\hline\n" ) ;
         fprintf( outfp, "       $N_b\\ge3$, MET2-4 &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f     \\vsmtvs \\\\ \n",
           frac_3b_val[ci][3][0], frac_3b_err[ci][3][0],
           frac_3b_val[ci][3][1], frac_3b_err[ci][3][1],
           frac_3b_val[ci][3][2], frac_3b_err[ci][3][2],
           frac_3b_val[ci][3][3], frac_3b_err[ci][3][3],
           frac_3b_val[ci][3][4], frac_3b_err[ci][3][4]
           ) ;

      } // ci.



      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "    \\end{tabular}\n" ) ;
      fprintf( outfp, "    \\label{tab:fracs-nb3}\n") ;
      fprintf( outfp, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n") ;
      fprintf( outfp, "\n\n\n" ) ;



   } // fitfrac_3b_table.

  //================================================================================================================================


   void fitfrac_2b_table( FILE* outfp, const char* table_title ) {

      fprintf( outfp, "\n\n\n" ) ;
      fprintf( outfp, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n") ;
      fprintf( outfp, "    \\begin{tabular}{|l||c|c|c||c|}\n") ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "    \\multicolumn{5}{|c|}{ Component fractions, %s  \\vsmtvs } \\\\ \n",  table_title ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "       &  HT2 & HT3 & HT4 & HT2-4 \\vsmtvs \\\\ \n" ) ;
      fprintf( outfp, "    \\hline\n") ;

      for ( int ci=0; ci<3; ci++ ) {

         fprintf( outfp, "          %s fraction       &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f    &  %5.2f $\\pm$ %5.2f      \\vsmtvs \\\\ \n",
           comp_out_name[ci],
           frac_2b_val[ci][0], frac_2b_err[ci][0],
           frac_2b_val[ci][1], frac_2b_err[ci][1],
           frac_2b_val[ci][2], frac_2b_err[ci][2],
           frac_2b_val[ci][3], frac_2b_err[ci][3]
           ) ;
         fprintf( outfp, "    \\hline\n") ;

      } // ci.
      fprintf( outfp, "    \\end{tabular}\n" ) ;
      fprintf( outfp, "    \\label{tab:fracs-nb2}\n") ;
      fprintf( outfp, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n") ;

      fprintf( outfp, "\n\n\n" ) ;

   } // fitfrac_2b_table.

  //================================================================================================================================





