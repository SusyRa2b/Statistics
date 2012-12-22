
#include "TSystem.h"
#include <iostream>

   bool getScanDataLine( const char* inFile,
                      const char* tagString,
                      float& val, float& p1sig, float& m1sig, float& p2sig, float& m2sig ) {

      val = -1. ;
      p1sig = 0. ;
      m1sig = 0. ;
      p2sig = 0. ;
      m2sig = 0. ;

      char command[10000] ;
      sprintf( command, "grep \"%s\" %s\n", tagString, inFile ) ;
      TString commandOutput = gSystem->GetFromPipe( command ) ;

      printf( " Output of command is : %s\n", commandOutput.Data() ) ;

      char fname[1000] ;
      char varname[1000] ;
      char formatstring[1000] ;
      sscanf( commandOutput.Data(), "%s %s %s %g %g %g %g %g",
             fname, varname, formatstring, &val, &p1sig, &m1sig, &p2sig, &m2sig ) ;

      return true ;


   }


