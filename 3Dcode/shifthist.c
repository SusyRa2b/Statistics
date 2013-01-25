
#include "TH1F.h"


   TH1F* shifthist( TH1F* hin, double delta = 0.1 ) {

      if ( hin == 0x0 ) {
         printf("\n\n *** bad input hist: %p\n\n", hin ) ;
         return 0x0 ;
      }

      int nbins = hin -> GetNbinsX() ;
      double xlow  = hin -> GetXaxis() -> GetBinLowEdge( 1 ) ;
      double xhigh = hin -> GetXaxis() -> GetBinUpEdge( nbins ) ;

      printf(" %s : nbins=%d, xlow=%6.2f, xhigh=%6.2f\n", hin->GetName(), nbins, xlow, xhigh ) ;

      char hname[1000] ;
      sprintf( hname, "%s_shift", hin->GetName() ) ;

      TH1F* retval = new TH1F( hname, hin->GetTitle(), nbins, xlow-delta, xhigh-delta ) ;

      for ( int bi=1; bi<=nbins; bi++ ) {
         retval -> GetXaxis() -> SetBinLabel( bi, hin->GetXaxis()->GetBinLabel( bi ) ) ;
         retval -> SetBinContent( bi, hin->GetBinContent(bi) ) ;
         retval -> SetBinError( bi, hin->GetBinError(bi) ) ;
      }

      retval -> SetMarkerStyle( hin->GetMarkerStyle() ) ;
      retval -> SetLineColor( hin->GetLineColor() ) ;
      retval -> SetLineWidth( hin->GetLineWidth() ) ;

      retval -> SetLabelSize( hin->GetLabelSize(),"x") ;
      retval -> GetXaxis() -> LabelsOption("v") ;

      return retval ;

   }


