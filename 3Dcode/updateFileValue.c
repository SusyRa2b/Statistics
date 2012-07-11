
#include "TSystem.h"

   bool updateFileValue( const char* inFile,
                         const char* parameterName,
                         double newValue, 
			 const char* ufvname) {


      FILE* infp ;
      if ( (infp=fopen( inFile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", inFile ) ;
         return false ;
      }

      FILE* outfp ;
      if ( (outfp=fopen( ufvname,"w"))==NULL ) {
         printf("\n\n *** Problem opening output file.\n\n" ) ;
         return false ;
      }

      //--- read in description line.
      char c(0) ;
      while ( c!=10  ) {
         c = fgetc( infp ) ;
         fprintf( outfp, "%c", c ) ;
      }

      //--- parameters.
      while ( !feof(infp) ) {
         char label[1000] ;
         float value ;
         fscanf( infp, "%s %g", label, &value ) ;
         if ( strcmp( label, parameterName ) == 0 ) {
            printf(" Found %s.  Changing value from %g to %g\n", parameterName, value, newValue ) ;
            value = newValue ;
         }
         if ( !feof(infp) ) {
            fprintf( outfp, "%s   %g\n", label, value ) ;
         }
      }

      fclose( infp ) ;
      fclose( outfp ) ;

      char command[10000] ;
      sprintf( command, "mv %s  %s\n", ufvname, inFile ) ;
      gSystem->Exec( command ) ;

      return true ;

   }


