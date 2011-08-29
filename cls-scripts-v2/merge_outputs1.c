
#include <stdio.h>

   void merge_outputs1( const char* selection = "ge1btight" ) {



      char infile[10000] ;
      char outfile[10000] ;
      char c(0) ;
      FILE* infp ;
      FILE* outfp ;


      int nPoints(0) ;

      int   mg[1000] ;
      int   mlsp[1000] ;
      float eff[1000] ;
      float plnsigul[1000] ;
      float plxsecul[1000] ;
      float clsxsec1[1000] ;
      float clsxsec2[1000] ;
      float clsxsec3[1000] ;
      float clsxsec4[1000] ;
      float clsxsec5[1000] ;
      float clsxsec6[1000] ;
      float clsxsec7[1000] ;
      float clsxsec8[1000] ;
      float clsval1[1000] ;
      float clsval2[1000] ;
      float clsval3[1000] ;
      float clsval4[1000] ;
      float clsval5[1000] ;
      float clsval6[1000] ;
      float clsval7[1000] ;
      float clsval8[1000] ;

      for ( int i=0; i<1000; i++ ) {
         mg[i] = 0 ;
         mlsp[i] = 0 ;
         eff[i] = 0. ;
         plnsigul[i] = 0. ;
         plxsecul[i] = 0. ;
         clsxsec1[i] = 999. ;
         clsxsec2[i] = 999. ;
         clsxsec3[i] = 999. ;
         clsxsec4[i] = 999. ;
         clsxsec5[i] = 999. ;
         clsxsec6[i] = 999. ;
         clsxsec7[i] = 999. ;
         clsxsec8[i] = 999. ;
         clsval1[i] = 99. ;
         clsval2[i] = 99. ;
         clsval3[i] = 99. ;
         clsval4[i] = 99. ;
         clsval5[i] = 99. ;
         clsval6[i] = 99. ;
         clsval7[i] = 99. ;
         clsval8[i] = 99. ;
      }




      //--- read in the output of the PL analysis: t1bbbb-scan-points-50GevBins-ge1btight.txt

      sprintf( infile, "t1bbbb-scan-points-50GevBins-%s.txt", selection ) ;

      if ( (infp=fopen( infile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
         return ;
      }

      //-- read in header (1st line).
      while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }

      //-- read in file.
      while ( !feof( infp ) ) {

         int tmg, tmlsp ;
         float teff, tplnsigul, tplxsecul ;
         fscanf( infp, "%d %d %f %f %f", &tmg, &tmlsp, &teff, &tplnsigul, &tplxsecul ) ;
         if ( feof(infp) ) break ;

         int pi = nPoints ;

         mg[pi] = tmg ;
         mlsp[pi] = tmlsp ;
         eff[pi] = teff ;
         plnsigul[pi] = tplnsigul ;
         plxsecul[pi] = tplxsecul ;
         nPoints++ ;

      } // eof?

      fclose(infp) ;







      //--- Iteration 1 ======================================================================

      sprintf( infile, "output-%s-t1bbbb-iter1.txt", selection ) ;

      if ( (infp=fopen( infile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
         return ;
      }

      //-- read in file.
      while ( !feof( infp ) ) {

         int tmg, tmlsp ;
         float tcls, tpspb, tpbgo, txsec ;
         fscanf( infp, "%d %d %f %f %f %f", &tmg, &tmlsp, &tcls, &tpspb, &tpbgo, &txsec ) ;
         if ( feof(infp) ) break ;

         int thisInd(-1) ;
         for ( int i=0; i<nPoints; i++ ) {
            if ( tmg == mg[i] && tmlsp == mlsp[i] ) { thisInd = i  ; break ; }
         }
         if ( thisInd<0 ) {
            printf("\n\n *** point not found!  mg=%d, mlsp=%d\n", tmg, tmlsp ) ;
         } else {
            clsxsec1[thisInd] = txsec ;
            clsval1[thisInd] = tcls ;
         }

      } // eof?

      fclose(infp) ;











      //--- Iteration 2 ======================================================================

      sprintf( infile, "output-%s-t1bbbb-iter2.txt", selection ) ;

      if ( (infp=fopen( infile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
         return ;
      }

      //-- read in file.
      while ( !feof( infp ) ) {

         int tmg, tmlsp ;
         float tcls, tpspb, tpbgo, txsec ;
         fscanf( infp, "%d %d %f %f %f %f", &tmg, &tmlsp, &tcls, &tpspb, &tpbgo, &txsec ) ;
         if ( feof(infp) ) break ;

         int thisInd(-1) ;
         for ( int i=0; i<nPoints; i++ ) {
            if ( tmg == mg[i] && tmlsp == mlsp[i] ) { thisInd = i  ; break ; }
         }
         if ( thisInd<0 ) {
            printf("\n\n *** point not found!  mg=%d, mlsp=%d\n", tmg, tmlsp ) ;
         } else {
            clsxsec2[thisInd] = txsec ;
            clsval2[thisInd] = tcls ;
         }

      } // eof?

      fclose(infp) ;








      //--- Iteration 3 ======================================================================

      sprintf( infile, "output-%s-t1bbbb-iter3.txt", selection ) ;

      if ( (infp=fopen( infile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
         return ;
      }

      //-- read in file.
      while ( !feof( infp ) ) {

         int tmg, tmlsp ;
         float tcls, tpspb, tpbgo, txsec ;
         fscanf( infp, "%d %d %f %f %f %f", &tmg, &tmlsp, &tcls, &tpspb, &tpbgo, &txsec ) ;
         if ( feof(infp) ) break ;

         int thisInd(-1) ;
         for ( int i=0; i<nPoints; i++ ) {
            if ( tmg == mg[i] && tmlsp == mlsp[i] ) { thisInd = i  ; break ; }
         }
         if ( thisInd<0 ) {
            printf("\n\n *** point not found!  mg=%d, mlsp=%d\n", tmg, tmlsp ) ;
         } else {
            clsxsec3[thisInd] = txsec ;
            clsval3[thisInd] = tcls ;
         }

      } // eof?

      fclose(infp) ;










      //--- Iteration 4 ======================================================================

      sprintf( infile, "output-%s-t1bbbb-iter4.txt", selection ) ;

      if ( (infp=fopen( infile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
         return ;
      }

      //-- read in file.
      while ( !feof( infp ) ) {

         int tmg, tmlsp ;
         float tcls, tpspb, tpbgo, txsec ;
         fscanf( infp, "%d %d %f %f %f %f", &tmg, &tmlsp, &tcls, &tpspb, &tpbgo, &txsec ) ;
         if ( feof(infp) ) break ;

         int thisInd(-1) ;
         for ( int i=0; i<nPoints; i++ ) {
            if ( tmg == mg[i] && tmlsp == mlsp[i] ) { thisInd = i  ; break ; }
         }
         if ( thisInd<0 ) {
            printf("\n\n *** point not found!  mg=%d, mlsp=%d\n", tmg, tmlsp ) ;
         } else {
            clsxsec4[thisInd] = txsec ;
            clsval4[thisInd] = tcls ;
         }

      } // eof?

      fclose(infp) ;










      //--- Iteration 5 ======================================================================

      sprintf( infile, "output-%s-t1bbbb-iter5.txt", selection ) ;

      if ( (infp=fopen( infile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
         return ;
      }

      //-- read in file.
      while ( !feof( infp ) ) {

         int tmg, tmlsp ;
         float tcls, tpspb, tpbgo, txsec ;
         fscanf( infp, "%d %d %f %f %f %f", &tmg, &tmlsp, &tcls, &tpspb, &tpbgo, &txsec ) ;
         if ( feof(infp) ) break ;

         int thisInd(-1) ;
         for ( int i=0; i<nPoints; i++ ) {
            if ( tmg == mg[i] && tmlsp == mlsp[i] ) { thisInd = i  ; break ; }
         }
         if ( thisInd<0 ) {
            printf("\n\n *** point not found!  mg=%d, mlsp=%d\n", tmg, tmlsp ) ;
         } else {
            clsxsec5[thisInd] = txsec ;
            clsval5[thisInd] = tcls ;
         }

      } // eof?

      fclose(infp) ;










      //--- Iteration 6 ======================================================================

      sprintf( infile, "output-%s-t1bbbb-iter6.txt", selection ) ;

      if ( (infp=fopen( infile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
         return ;
      }

      //-- read in file.
      while ( !feof( infp ) ) {

         int tmg, tmlsp ;
         float tcls, tpspb, tpbgo, txsec ;
         fscanf( infp, "%d %d %f %f %f %f", &tmg, &tmlsp, &tcls, &tpspb, &tpbgo, &txsec ) ;
         if ( feof(infp) ) break ;

         int thisInd(-1) ;
         for ( int i=0; i<nPoints; i++ ) {
            if ( tmg == mg[i] && tmlsp == mlsp[i] ) { thisInd = i  ; break ; }
         }
         if ( thisInd<0 ) {
            printf("\n\n *** point not found!  mg=%d, mlsp=%d\n", tmg, tmlsp ) ;
         } else {
            clsxsec6[thisInd] = txsec ;
            clsval6[thisInd] = tcls ;
         }

      } // eof?

      fclose(infp) ;






      //--- Iteration 7 ======================================================================

      sprintf( infile, "output-%s-t1bbbb-iter7.txt", selection ) ;

      if ( (infp=fopen( infile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
         return ;
      }

      //-- read in file.
      while ( !feof( infp ) ) {

         int tmg, tmlsp ;
         float tcls, tpspb, tpbgo, txsec ;
         fscanf( infp, "%d %d %f %f %f %f", &tmg, &tmlsp, &tcls, &tpspb, &tpbgo, &txsec ) ;
         if ( feof(infp) ) break ;

         int thisInd(-1) ;
         for ( int i=0; i<nPoints; i++ ) {
            if ( tmg == mg[i] && tmlsp == mlsp[i] ) { thisInd = i  ; break ; }
         }
         if ( thisInd<0 ) {
            printf("\n\n *** point not found!  mg=%d, mlsp=%d\n", tmg, tmlsp ) ;
         } else {
            clsxsec7[thisInd] = txsec ;
            clsval7[thisInd] = tcls ;
         }

      } // eof?

      fclose(infp) ;






      //--- Iteration 8 ======================================================================

      sprintf( infile, "output-%s-t1bbbb-iter8.txt", selection ) ;

      if ( (infp=fopen( infile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
         return ;
      }

      //-- read in file.
      while ( !feof( infp ) ) {

         int tmg, tmlsp ;
         float tcls, tpspb, tpbgo, txsec ;
         fscanf( infp, "%d %d %f %f %f %f", &tmg, &tmlsp, &tcls, &tpspb, &tpbgo, &txsec ) ;
         if ( feof(infp) ) break ;

         int thisInd(-1) ;
         for ( int i=0; i<nPoints; i++ ) {
            if ( tmg == mg[i] && tmlsp == mlsp[i] ) { thisInd = i  ; break ; }
         }
         if ( thisInd<0 ) {
            printf("\n\n *** point not found!  mg=%d, mlsp=%d\n", tmg, tmlsp ) ;
         } else {
            clsxsec8[thisInd] = txsec ;
            clsval8[thisInd] = tcls ;
         }

      } // eof?

      fclose(infp) ;











      //--- output ============================================================================


      sprintf( outfile, "t1bbbb-all-output-%s.txt", selection ) ;

      if ( (outfp=fopen( outfile,"w"))==NULL ) {
         printf("\n\n *** Problem opening output file: %s.\n\n", outfile ) ;
         return ;
      }

      printf("\n\n") ;
      printf(        "mg:mlsp:Eff:PLNsigUL:PLXsecUL:clsXsec1:clsXsec2:clsXsec3:clsXsec4:clsXsec5:clsXsec6:clsXsec7:clsXsec8:clsVal1:clsVal2:clsVal3:clsVal4:clsVal5:clsVal6:clsVal7:clsVal8\n") ;
      fprintf(outfp, "mg:mlsp:Eff:PLNsigUL:PLXsecUL:clsXsec1:clsXsec2:clsXsec3:clsXsec4:clsXsec5:clsXsec6:clsXsec7:clsXsec8:clsVal1:clsVal2:clsVal3:clsVal4:clsVal5:clsVal6:clsVal7:clsVal8\n") ;
      for ( int i=0; i<nPoints; i++ ) {
       //if ( clsval1[i] > 9. ) continue ;
         printf(        "  %7d %7d  %8.4f  %8.3f %8.3f   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
             mg[i], mlsp[i], eff[i], plnsigul[i], plxsecul[i], clsxsec1[i], clsxsec2[i], clsxsec3[i], clsxsec4[i], clsxsec5[i], clsxsec6[i], clsxsec7[i], clsxsec8[i], clsval1[i], clsval2[i], clsval3[i], clsval4[i], clsval5[i], clsval6[i], clsval7[i], clsval8[i]  ) ;
         fprintf(outfp, "  %7d %7d  %8.4f  %8.3f %8.3f   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
             mg[i], mlsp[i], eff[i], plnsigul[i], plxsecul[i], clsxsec1[i], clsxsec2[i], clsxsec3[i], clsxsec4[i], clsxsec5[i], clsxsec6[i], clsxsec7[i], clsxsec8[i], clsval1[i], clsval2[i], clsval3[i], clsval4[i], clsval5[i], clsval6[i], clsval7[i], clsval8[i]  ) ;
      } // i.
      printf("\n\n") ;


      fclose( outfp ) ;

      printf("\n\n Created output file : %s\n\n\n\n", outfile ) ;


   }

