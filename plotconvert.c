//
//  $Id: plotconvert.c 30 2011-02-09 02:02:19Z owenl $
//
//

#include "TH1.h"
#include "RooPlot.h"
#include "RooHist.h"
#include <stdio.h>


    TH1F* plotconvert( RooPlot* rp, const char* hname ) {
       if ( rp==0 ) return 0x0 ;
       double xl = rp->GetXaxis()->GetXmin() ;
       double xh = rp->GetXaxis()->GetXmax() ;
       int nbins = rp->GetNbinsX() ;
       printf(" plotconvert: GetNbinsX returns %d\n", nbins ) ;
       int nbins2 = rp->getHist()->GetN() ;
       printf(" plotconvert: GetN returns %d\n", nbins2 ) ;
       TH1F* hret = new TH1F( hname, "", nbins, xl, xh ) ;
       for ( int i=1; i<=nbins; i++ ) {
          double x, y ;
          rp->getHist()->GetPoint(i,x,y) ;
      //  printf(" bin %3d : %7.0f\n", i, y ) ;
          hret->SetBinContent( i, y ) ;
       } // i
       return hret ;
    }


