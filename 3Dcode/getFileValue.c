
#include "TSystem.h"
#include <iostream>

   bool getFileValue( const char* inFile,
                      const char* parameterName,
                      float& returnValue ) {


      returnValue = 1.0 ;

    //-- New way that doesn't use fopen.

      char command[10000] ;
      sprintf( command, "grep %s %s\n", parameterName, inFile ) ;
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

   //--- old way below here.  Runs into too many open files, even though
   //    I close them with fclose.  Don't know why...

 //// FILE* infp ;
 //// if ( (infp=fopen( inFile,"r"))==NULL ) {
 ////    printf("\n\n *** getFileValue: Problem opening input file: %s.\n\n", inFile ) ;
 ////    char command[10000] ;
 ////    sprintf( command, "ls -l %s\n", inFile ) ;
 ////    gSystem->Exec( command ) ;
 ////    cout << flush ;
 ////    return false ;
 //// }

 //// //--- read in description line.
 //// char c(0) ;
 //// while ( c!=10  ) {
 ////    c = fgetc( infp ) ;
 //// }

 //// //--- parameters.
 //// while ( !feof(infp) ) {
 ////    char label[1000] ;
 ////    float value ;
 ////    fscanf( infp, "%s %g", label, &value ) ;
 ////    if ( strcmp( label, parameterName ) == 0 ) {
 ////       printf(" Found %s.  Value is %g\n", parameterName, value ) ;
 ////       returnValue = value ;
 ////       return true ;
 ////    }
 //// }

 //// fclose( infp ) ;

      printf("\n\n *** Could not find parameter %s in file %s.\n\n", parameterName, inFile ) ;

      return false ;

   }


