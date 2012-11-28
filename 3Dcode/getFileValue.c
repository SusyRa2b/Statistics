
#include "TSystem.h"
#include <iostream>

   bool getFileValue( const char* inFile,
                      const char* parameterName,
                      float& returnValue ) {


      returnValue = 1.0 ;

    //-- New way that doesn't use fopen.
    //-- Nov 27, 2012: Include a blank space at the end to avoid multiple
    //                 matches in grep when target is a substring found in other lines.

      char command[10000] ;
      sprintf( command, "grep \"%s \" %s\n", parameterName, inFile ) ;
      TString commandOutput = gSystem->GetFromPipe( command ) ;

      /// printf( " Output of command is : %s\n", commandOutput.Data() ) ;

      char label[1000] ;
      float value ;
      sscanf( commandOutput.Data(), "%s %g", label, &value ) ;
      if ( strcmp( label, parameterName ) == 0 ) {
         printf(" Found %s.  Value is %g\n", parameterName, value ) ;
         returnValue = value ;
         return true ;
      }

      printf("\n\n *** Could not find parameter %s in file %s.\n\n", parameterName, inFile ) ;

      return false ;

   }


