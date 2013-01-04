

#include "fill_arrays.c"

#include <stdio.h>

   void generate_table1( const char* scan_results_file = "scan-results-19fb.txt" ) {


      fill_arrays( scan_results_file ) ;


      FILE* outfp ;
      if ( (outfp=fopen( "tables.tex", "w" ))==NULL ) {
         printf("\n\n *** can't open output file.\n\n") ; return ;
      }








      fprintf( outfp, "\\documentclass[11pt]{article}\n" ) ;
      fprintf( outfp, "\\def\\vsmtvs {\\rule[-0.3cm]{0cm}{0.75cm}}\n" ) ;
      fprintf( outfp, "\\begin{document}\n" ) ;

      fprintf( outfp, "\n\n\n" ) ;

      fprintf( outfp, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n") ;
      fprintf( outfp, "    \\begin{tabular}{|l||c|c|c||c|}\n") ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "       &  HT2 & HT3 & HT4 & HT2-4 \\vsmtvs \\\\ \n" ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "          Observed events                     &  %.0f  & %.0f  & %.0f  &  %.0f  \\vsmtvs \\\\ \n",
        data_2b_obs[0], data_2b_obs[1], data_2b_obs[2], data_2b_obs[3] ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "          Unbiased SM background predictions  & $%.1f\\ ^{+%.1f}_{-%.1f}$    & $%.1f\\ ^{+%.1f}_{-%.1f}$    & $%.1f\\ ^{+%.1f}_{-%.1f}$    & $%.0f\\ ^{+%.0f}_{-%.0f}$  \\vsmtvs \\\\ \n",
               ub_2b_val[0], ub_2b_p1s[0], ub_2b_m1s[0],
               ub_2b_val[1], ub_2b_p1s[1], ub_2b_m1s[1],
               ub_2b_val[2], ub_2b_p1s[2], ub_2b_m1s[2],
               ub_2b_val[3], ub_2b_p1s[3], ub_2b_m1s[3]
               ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "          SM background from full fit  & $%.1f\\ ^{+%.1f}_{-%.1f}$    & $%.1f\\ ^{+%.1f}_{-%.1f}$    & $%.1f\\ ^{+%.1f}_{-%.1f}$    & $%.0f\\ ^{+%.0f}_{-%.0f}$  \\vsmtvs \\\\ \n",
               ff_2b_val[0], ff_2b_p1s[0], ff_2b_m1s[0],
               ff_2b_val[1], ff_2b_p1s[1], ff_2b_m1s[1],
               ff_2b_val[2], ff_2b_p1s[2], ff_2b_m1s[2],
               ff_2b_val[3], ff_2b_p1s[3], ff_2b_m1s[3]
               ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "          MC value  &    $%.1f \\pm %.1f$  &       $%.1f \\pm %.1f$  &    $%.1f \\pm %.1f$  &    $%.0f \\pm %.0f$  \\vsmtvs \\\\ \n",
               mc_2b_val[0], mc_2b_err[0],
               mc_2b_val[1], mc_2b_err[1],
               mc_2b_val[2], mc_2b_err[2],
               mc_2b_val[3], mc_2b_err[3]
               ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "    \\end{tabular}\n" ) ;
      fprintf( outfp, "    \\label{tab:pred-and-observed-nb2}\n") ;
      fprintf( outfp, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n") ;

      fprintf( outfp, "\n\n\n" ) ;



      fprintf( outfp, "    \\begin{tabular}{|c||c|c|c|c||c|}\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\multicolumn{6}{|c|}{ Observed number of events} \\\\ \n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "      & HT1 & HT2 & HT3 & HT4 & HT1-4 \\\\ \n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET2  &  %.0f  &  %.0f  &  %.0f  &  %.0f  &  %.0f  \\vsmtvs \\\\ \n",
           data_3b_obs[0][0],
           data_3b_obs[0][1],
           data_3b_obs[0][2],
           data_3b_obs[0][3],
           data_3b_obs[0][4] ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET3  &  %.0f  &  %.0f  &  %.0f  &  %.0f  &  %.0f  \\vsmtvs \\\\ \n",
           data_3b_obs[1][0],
           data_3b_obs[1][1],
           data_3b_obs[1][2],
           data_3b_obs[1][3],
           data_3b_obs[1][4] ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET4  &  ---  &  %.0f  &  %.0f  &  %.0f  &  %.0f  \\vsmtvs \\\\ \n",
           data_3b_obs[2][1],
           data_3b_obs[2][2],
           data_3b_obs[2][3],
           data_3b_obs[2][4] ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET2-4  &  %.0f  &  %.0f  &  %.0f  &  %.0f  &  %.0f  \\vsmtvs \\\\ \n",
           data_3b_obs[3][0],
           data_3b_obs[3][1],
           data_3b_obs[3][2],
           data_3b_obs[3][3],
           data_3b_obs[3][4] ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "    \\multicolumn{6}{|c|}{ Unbiased background predictions, SM-only} \\\\ \n" ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "      & HT1 & HT2 & HT3 & HT4 & HT1-4 \\\\ \n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET2  &  $%.0f\\ ^{+%.0f}_{-%.0f}$    &  $%.0f\\ ^{+%.0f}_{-%.0f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.0f\\ ^{+%.0f}_{-%.0f}$   \\vsmtvs \\\\ \n"
          ,ub_3b_val[0][0], ub_3b_p1s[0][0], ub_3b_m1s[0][0]
          ,ub_3b_val[0][1], ub_3b_p1s[0][1], ub_3b_m1s[0][1]
          ,ub_3b_val[0][2], ub_3b_p1s[0][2], ub_3b_m1s[0][2]
          ,ub_3b_val[0][3], ub_3b_p1s[0][3], ub_3b_m1s[0][3]
          ,ub_3b_val[0][4], ub_3b_p1s[0][4], ub_3b_m1s[0][4]
           ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET3  &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$   \\vsmtvs \\\\ \n"
          ,ub_3b_val[1][0], ub_3b_p1s[1][0], ub_3b_m1s[1][0]
          ,ub_3b_val[1][1], ub_3b_p1s[1][1], ub_3b_m1s[1][1]
          ,ub_3b_val[1][2], ub_3b_p1s[1][2], ub_3b_m1s[1][2]
          ,ub_3b_val[1][3], ub_3b_p1s[1][3], ub_3b_m1s[1][3]
          ,ub_3b_val[1][4], ub_3b_p1s[1][4], ub_3b_m1s[1][4]
           ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET4  &  ---    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$   \\vsmtvs \\\\ \n"
          ,ub_3b_val[2][1], ub_3b_p1s[2][1], ub_3b_m1s[2][1]
          ,ub_3b_val[2][2], ub_3b_p1s[2][2], ub_3b_m1s[2][2]
          ,ub_3b_val[2][3], ub_3b_p1s[2][3], ub_3b_m1s[2][3]
          ,ub_3b_val[2][4], ub_3b_p1s[2][4], ub_3b_m1s[2][4]
           ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET2-4  &  $%.0f\\ ^{+%.0f}_{-%.0f}$    &  $%.0f\\ ^{+%.0f}_{-%.0f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.0f\\ ^{+%.0f}_{-%.0f}$   \\vsmtvs \\\\ \n"
          ,ub_3b_val[3][0], ub_3b_p1s[3][0], ub_3b_m1s[3][0]
          ,ub_3b_val[3][1], ub_3b_p1s[3][1], ub_3b_m1s[3][1]
          ,ub_3b_val[3][2], ub_3b_p1s[3][2], ub_3b_m1s[3][2]
          ,ub_3b_val[3][3], ub_3b_p1s[3][3], ub_3b_m1s[3][3]
          ,ub_3b_val[3][4], ub_3b_p1s[3][4], ub_3b_m1s[3][4]
           ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "    \\multicolumn{6}{|c|}{ Full fit results, SM-only} \\\\ \n" ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "      & HT1 & HT2 & HT3 & HT4 & HT1-4 \\\\ \n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET2  &  $%.0f\\ ^{+%.0f}_{-%.0f}$    &  $%.0f\\ ^{+%.0f}_{-%.0f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.0f\\ ^{+%.0f}_{-%.0f}$   \\vsmtvs \\\\ \n"
          ,ff_3b_val[0][0], ff_3b_p1s[0][0], ff_3b_m1s[0][0]
          ,ff_3b_val[0][1], ff_3b_p1s[0][1], ff_3b_m1s[0][1]
          ,ff_3b_val[0][2], ff_3b_p1s[0][2], ff_3b_m1s[0][2]
          ,ff_3b_val[0][3], ff_3b_p1s[0][3], ff_3b_m1s[0][3]
          ,ff_3b_val[0][4], ff_3b_p1s[0][4], ff_3b_m1s[0][4]
           ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET3  &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$   \\vsmtvs \\\\ \n"
          ,ff_3b_val[1][0], ff_3b_p1s[1][0], ff_3b_m1s[1][0]
          ,ff_3b_val[1][1], ff_3b_p1s[1][1], ff_3b_m1s[1][1]
          ,ff_3b_val[1][2], ff_3b_p1s[1][2], ff_3b_m1s[1][2]
          ,ff_3b_val[1][3], ff_3b_p1s[1][3], ff_3b_m1s[1][3]
          ,ff_3b_val[1][4], ff_3b_p1s[1][4], ff_3b_m1s[1][4]
           ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET4  &  ---    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$   \\vsmtvs \\\\ \n"
          ,ff_3b_val[2][1], ff_3b_p1s[2][1], ff_3b_m1s[2][1]
          ,ff_3b_val[2][2], ff_3b_p1s[2][2], ff_3b_m1s[2][2]
          ,ff_3b_val[2][3], ff_3b_p1s[2][3], ff_3b_m1s[2][3]
          ,ff_3b_val[2][4], ff_3b_p1s[2][4], ff_3b_m1s[2][4]
           ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET2-4  &  $%.0f\\ ^{+%.0f}_{-%.0f}$    &  $%.0f\\ ^{+%.0f}_{-%.0f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.1f\\ ^{+%.1f}_{-%.1f}$    &  $%.0f\\ ^{+%.0f}_{-%.0f}$   \\vsmtvs \\\\ \n"
          ,ff_3b_val[3][0], ff_3b_p1s[3][0], ff_3b_m1s[3][0]
          ,ff_3b_val[3][1], ff_3b_p1s[3][1], ff_3b_m1s[3][1]
          ,ff_3b_val[3][2], ff_3b_p1s[3][2], ff_3b_m1s[3][2]
          ,ff_3b_val[3][3], ff_3b_p1s[3][3], ff_3b_m1s[3][3]
          ,ff_3b_val[3][4], ff_3b_p1s[3][4], ff_3b_m1s[3][4]
           ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "    \\multicolumn{6}{|c|}{ Monte Carlo, SM-only} \\\\ \n" ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "      & HT1 & HT2 & HT3 & HT4 & HT1-4 \\\\ \n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET2  &  $%.0f \\pm %.0f$    &  $%.0f \\pm %.0f$  &  $%.0f \\pm %.0f$  &  $%.0f \\pm %.0f$  &  $%.0f \\pm %.0f$ \\vsmtvs \\\\ \n"
          ,mc_3b_val[0][0], mc_3b_err[0][0]
          ,mc_3b_val[0][1], mc_3b_err[0][1]
          ,mc_3b_val[0][2], mc_3b_err[0][2]
          ,mc_3b_val[0][3], mc_3b_err[0][3]
          ,mc_3b_val[0][4], mc_3b_err[0][4]
           ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET3  &  $%.1f \\pm %.1f$    &  $%.1f \\pm %.1f$  &  $%.1f \\pm %.1f$  &  $%.1f \\pm %.1f$  &  $%.1f \\pm %.1f$ \\vsmtvs \\\\ \n"
          ,mc_3b_val[1][0], mc_3b_err[1][0]
          ,mc_3b_val[1][1], mc_3b_err[1][1]
          ,mc_3b_val[1][2], mc_3b_err[1][2]
          ,mc_3b_val[1][3], mc_3b_err[1][3]
          ,mc_3b_val[1][4], mc_3b_err[1][4]
           ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET4  &  ---    &  $%.1f \\pm %.1f$  &  $%.1f \\pm %.1f$  &  $%.1f \\pm %.1f$  &  $%.1f \\pm %.1f$ \\vsmtvs \\\\ \n"
          ,mc_3b_val[2][1], mc_3b_err[2][1]
          ,mc_3b_val[2][2], mc_3b_err[2][2]
          ,mc_3b_val[2][3], mc_3b_err[2][3]
          ,mc_3b_val[2][4], mc_3b_err[2][4]
           ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "    \\hline\n" ) ;
      fprintf( outfp, "       $N_b\\ge3$, MET2-4  &  $%.0f \\pm %.0f$    &  $%.0f \\pm %.0f$  &  $%.0f \\pm %.0f$  &  $%.0f \\pm %.0f$  &  $%.0f \\pm %.0f$ \\vsmtvs \\\\ \n"
          ,mc_3b_val[3][0], mc_3b_err[3][0]
          ,mc_3b_val[3][1], mc_3b_err[3][1]
          ,mc_3b_val[3][2], mc_3b_err[3][2]
          ,mc_3b_val[3][3], mc_3b_err[3][3]
          ,mc_3b_val[3][4], mc_3b_err[3][4]
           ) ;
      fprintf( outfp, "    \\hline\n") ;
      fprintf( outfp, "    \\end{tabular}\n" ) ;
      fprintf( outfp, "    \\label{tab:pred-and-observed-nb3}\n") ;
      fprintf( outfp, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n") ;

      fprintf( outfp, "\n\n\n" ) ;










      fprintf( outfp, "\\end{document}\n" ) ;

      fclose( outfp ) ;

   }

