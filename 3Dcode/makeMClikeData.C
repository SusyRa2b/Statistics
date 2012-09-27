#include "TLine.h"
#include "TString.h"

#include "updateFileValue.c"
#include "getFileValue.c"

#include <iostream>
#include <fstream>

using std::cout ;							    
using std::endl ;							    

void makeMClikeData( const char* datfile = "Input-met4-ht4-wsyst1-testing.dat" ) {  

   const int nBinsMET = 4 ;
   const int nBinsHT = 4 ;
   const int nBinsBjets = 3;
   const int nSelTrig = 2;
   char seltrigname[2][100] = { "0L", "1L" } ;
   int nSel = 3;
   char selname[3][100] = { "0lep", "1lep", "ldp" } ;

   float trigeff[nSelTrig][nBinsMET][nBinsHT];

   for ( int si=0; si<nSelTrig; si++ ) {
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         char pname[1000] ;
         sprintf( pname, "trigeff_val_%s_M%d_H%d", seltrigname[si], mbi+1, hbi+1 ) ;
         getFileValue( datfile, pname, trigeff[si][mbi][hbi] ) ;         
	 }
      }
   }

   for ( int si=0; si<nSel; si++ ) {
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               int selbin = 0;
	       if (si==1) selbin = 1;
               char pname[1000] ;
               sprintf( pname, "N_%s_M%d_H%d_%db", selname[si], mbi+1, hbi+1, bbi+1 ) ;
               float oldvalue;
	       getFileValue( datfile, pname, oldvalue ) ;   
               updateFileValue( datfile, pname, oldvalue*trigeff[selbin][mbi][hbi] ) ;
	    }
	 }
      }
   }
  
}
