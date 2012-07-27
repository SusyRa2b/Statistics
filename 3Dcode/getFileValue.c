
#include "TSystem.h"

   bool getFileValue( const char* inFile,
                      const char* parameterName,
                      float& returnValue ) {


      returnValue = 1.0 ;

      FILE* infp ;
      if ( (infp=fopen( inFile,"r"))==NULL ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", inFile ) ;
         return false ;
      }

      //--- read in description line.
      char c(0) ;
      while ( c!=10  ) {
         c = fgetc( infp ) ;
      }

      //--- parameters.
      while ( !feof(infp) ) {
         char label[1000] ;
         float value ;
         fscanf( infp, "%s %g", label, &value ) ;
         if ( strcmp( label, parameterName ) == 0 ) {
            printf(" Found %s.  Value is %g\n", parameterName, value ) ;
            returnValue = value ;
            return true ;
         }
      }

      fclose( infp ) ;

      printf("\n\n *** Could not find parameter %s in file %s.\n\n", parameterName, inFile ) ;

      return false ;

   }


