
#include "updateFileValue.c"


   void makeblind( const char* datfile, const char* blindBinsList="macros1/blind-bins-list.txt" ) {

      FILE* bbl_file ;
      if ( (bbl_file = fopen( blindBinsList, "r"))==NULL ) {
         printf("\n\n *** Problem opening blindBinsList file: %s\n\n", blindBinsList ) ;
         return ;
      }
      printf("\n\n Reading blindBinsList file: %s\n\n", blindBinsList ) ;
      while ( ! feof(bbl_file) ) {
         char binlabel[1000] ;
         fscanf( bbl_file, "%s", binlabel ) ;
         if ( feof(bbl_file) ) break ;
         int bmb(-1), bhb(-1), bbb(-1) ;
         int rv = sscanf( binlabel, "N_0lep_M%d_H%d_%db", &bmb, &bhb, &bbb ) ;
         if ( rv == 0 ) {
            printf("\n\n *** bad bin label format: %s\n\n", binlabel ) ;
            return  ;
         }
         if ( bmb < 1 || bhb < 1 || bbb < 1 || bmb > 4 || bhb > 4 || bbb > 3 ) {
            printf("\n\n *** bin index out of range: met=%d, ht=%d, nb=%d\n\n", bmb, bhb, bbb ) ;
            return  ;
         }
         printf("   bin label %s : will blind met=%d, ht=%d, nb=%d 0lep bin.\n", binlabel, bmb, bhb, bbb ) ;
         updateFileValue( datfile, binlabel, 0 ) ;
      } // still reading?


   } // makeblind


