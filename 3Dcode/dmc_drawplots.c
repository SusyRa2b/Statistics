
#include "TH1F.h"
#include "THStack.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TText.h"

#include <iostream>



    int nComps(3) ;
    char compname[3][100] ;

    bool savePdf ;

  //----------
  // prototypes

   void drawSet( const char* hname_base, const char* xtitle ) ;
   void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;

  //----------

    void dmc_drawplots( const char* infile = "rootfiles/dmc_plots1.root", bool arg_savePdf = false ) {

       savePdf = arg_savePdf ;

       gDirectory->Delete("h*") ;

       gStyle -> SetOptStat(0) ;

       sprintf( compname[0], "data" ) ;
       sprintf( compname[1], "ttbar" ) ;
       sprintf( compname[2], "wjets" ) ;


       loadHist( infile ) ;

       gStyle->SetOptLogy(1) ;

       drawSet( "h_ht_sl_all", "HT (GeV)" ) ; // draw first one twice.
       drawSet( "h_ht_sl_all", "HT (GeV)" ) ; // draw first one twice.
       drawSet( "h_ht_sl_nb1", "HT (GeV)" ) ;
       drawSet( "h_ht_sl_nb2", "HT (GeV)" ) ;
       drawSet( "h_ht_sl_nb3", "HT (GeV)" ) ;

       drawSet( "h_met_sl_all", "MET (GeV)" ) ;
       drawSet( "h_met_sl_nb1", "MET (GeV)" ) ;
       drawSet( "h_met_sl_nb2", "MET (GeV)" ) ;
       drawSet( "h_met_sl_nb3", "MET (GeV)" ) ;

       drawSet( "h_nb_sl_all", "N btags" ) ;




       gStyle->SetOptLogy(0) ;

       drawSet( "h_ht_sl_all", "HT (GeV)" ) ;
       drawSet( "h_ht_sl_nb1", "HT (GeV)" ) ;
       drawSet( "h_ht_sl_nb2", "HT (GeV)" ) ;
       drawSet( "h_ht_sl_nb3", "HT (GeV)" ) ;

       drawSet( "h_met_sl_all", "MET (GeV)" ) ;
       drawSet( "h_met_sl_nb1", "MET (GeV)" ) ;
       drawSet( "h_met_sl_nb2", "MET (GeV)" ) ;
       drawSet( "h_met_sl_nb3", "MET (GeV)" ) ;

       drawSet( "h_nb_sl_all", "N btags" ) ;

    } // dmc_drawplots

  //--------------------------------------------------------


   void drawSet( const char* hname_base, const char* xtitle ) {

      bool islogy = gStyle->GetOptLogy() ;

      char cname[1000] ;
      if ( islogy ) {
         sprintf( cname, "can_logy_%s", hname_base ) ;
      } else {
         sprintf( cname, "can_%s", hname_base ) ;
      }
      TCanvas* dmccan = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( dmccan == 0x0 ) {
         dmccan = new TCanvas( cname, hname_base, 600, 750 ) ;
      }
      dmccan->Clear() ;

      char hname[1000] ;
      sprintf( hname, "%s_%s", hname_base, "data" ) ;
      TH1F* hdata = (TH1F*) gDirectory->FindObject( hname ) ;
      if ( hdata == 0x0 ) {
         printf("\n\n *** drawSet: can't find data hist with name %s\n\n", hname ) ; return ;
      }

      sprintf( hname, "%s_mcstack", "" ) ;
      THStack* hmcstack = new THStack() ;

      sprintf( hname, "%s_mcsum", "" ) ;
      TH1F* hmcsum = (TH1F*) hdata->Clone( hname ) ;
      hmcsum -> Reset() ;

      hdata -> SetLineWidth(2) ;
      hdata -> SetMarkerStyle(20) ;

      TLegend* legend = new TLegend( 0.80, 0.80, 0.95, 0.95 ) ;
      legend->SetFillColor(kWhite) ;

      for ( int ci=1; ci<nComps; ci++ ) {

         sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
         TH1F* hmc = (TH1F*) gDirectory->FindObject( hname ) ;
         if ( hmc == 0x0 ) { printf("\n\n *** drawSet: missing MC hist %s\n", hname ) ; return ; }
         hmcsum -> Add( hmc ) ;
         hmcstack -> Add( hmc ) ;
         legend -> AddEntry( hmc, compname[ci] ) ;

      }

      sprintf( hname, "%s_diff", hname_base ) ;
      TH1F* hdiff = (TH1F*) hdata->Clone( hname ) ;
      hdiff->Reset() ;

      for ( int bi=1; bi<hdata->GetNbinsX(); bi++ ) {
         double data = hdata->GetBinContent(bi) ;
         double data_err = hdata->GetBinError(bi) ;
         double mc = hmcsum->GetBinContent(bi) ;
         double mc_err = hmcsum->GetBinError(bi) ;
         double val = 0. ;
         double err = 0. ;
         if ( mc > 0 ) {
            val = data / mc ;
            double errsq(0.) ;
            if ( mc > 0 ) { errsq += pow(mc_err/mc,2) ; }
            if ( data > 0 ) { errsq += pow(data_err/data,2) ; }
            err = sqrt(errsq) ;
         }
         hdiff->SetBinContent(bi,val) ;
         hdiff->SetBinError(bi,err) ;
      } // bi.

      hdiff->SetMinimum(0.3) ;
      hdiff->SetMaximum(1.7) ;




      double hmax = hdata->GetBinContent( hdata->GetMaximumBin() ) ;
      if ( islogy ) {
         hmax = 3*hmax ;
      } else {
         hmax = 1.2*hmax ;
      }
      hdata->SetMaximum(hmax) ;
      hdata->SetXTitle( xtitle ) ;

      char padname[1000] ;
      sprintf( padname, "tp_%s", hname_base ) ;
      TPad* toppad = new TPad( padname, padname, 0.02, 0.3, 0.98, 0.98 ) ;
      toppad->Draw() ;
      sprintf( padname, "bp_%s", hname_base ) ;
      TPad* bottompad = new TPad( padname, padname, 0.02, 0.02, 0.98, 0.28 ) ;
      bottompad->Draw() ;

      hmcsum->SetMarkerStyle(0) ;

      toppad->cd() ;
      gStyle->SetOptTitle(0) ;
      hdata->DrawCopy() ;
      hmcstack->Draw("samehist") ;
      hmcsum->Draw("esame") ;
      hdata->DrawCopy("same") ;
      legend->Draw() ;
      toppad->Update() ;
      toppad->Draw() ;
      TText* title = new TText() ;
      title->DrawTextNDC( 0.05, 0.95, hdata->GetTitle() ) ;

      dmccan->Update() ;
      dmccan->Draw() ;


      bottompad->cd() ;
      hdiff->UseCurrentStyle() ;
      gPad->SetLogy(0) ;
      gStyle->SetNdivisions(404,"y") ;
      gStyle->SetOptTitle(0) ;
      gStyle->SetLabelSize(0.10,"x") ;
      gStyle->SetLabelSize(0.10,"y") ;
      gStyle->SetTitleSize(0.11,"y") ;
      gStyle->SetTitleOffset(0.4,"y") ;
      hdiff->SetLineWidth(2) ;
      hdiff->SetMarkerStyle(20) ;
      hdiff->Draw() ;
      hdiff->SetYTitle("Data/MC") ;
      TLine* line = new TLine() ;
      line->SetLineStyle(2) ;
      double xl = hdata->GetBinLowEdge(1) ;
      double xh = hdata->GetBinLowEdge( hdata->GetNbinsX() ) + hdata->GetBinWidth(1) ;
      line->DrawLine( xl, 1., xh, 1. ) ;

      dmccan->Update() ;
      dmccan->Draw() ;

      if ( savePdf ) {
         char filename[10000] ;
         if ( islogy ) {
            sprintf( filename, "outputfiles/%s_logy.pdf", hname_base ) ;
         } else {
            sprintf( filename, "outputfiles/%s.pdf", hname_base ) ;
         }
         dmccan->SaveAs( filename ) ;
      }

   } // drawSet




   //=======================================================================================
void loadHist(const char* filename, const char* pfx, const char* pat, Bool_t doAdd, Double_t scaleFactor)
{
  cout << " Reading histograms from file: " << filename << endl << flush ;
  TFile inf(filename) ;
  //inf.ReadAll() ;
  TList* list = inf.GetListOfKeys() ;
  TIterator* iter = list->MakeIterator();

  TRegexp re(pat,kTRUE) ;
  std::cout << "pat = " << pat << std::endl ;

  gDirectory->cd("Rint:") ;

  TObject* obj ;
  TKey* key ;
  std::cout << "doAdd = " << (doAdd?"T":"F") << std::endl ;
  std::cout << "loadHist: reading." ;
  while((key=(TKey*)iter->Next())) {
   
    Int_t ridx = TString(key->GetName()).Index(re) ;    
    if (ridx==-1) {
      continue ;
    }

    obj = inf.Get(key->GetName()) ;
    TObject* clone ;
    if (pfx) {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	oldObj = gDirectory->Get(Form("%s_%s",pfx,obj->GetName())) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }
      if (oldObj) {
	clone = oldObj ;
        if ( scaleFactor > 0 ) {
           ((TH1*)clone)->Sumw2() ;
           ((TH1*)clone)->Add((TH1*)obj, scaleFactor) ;
        } else {
           ((TH1*)clone)->Add((TH1*)obj) ;
        }
      } else {
	clone = obj->Clone(Form("%s_%s",pfx,obj->GetName())) ;
      }


    } else {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	oldObj = gDirectory->Get(key->GetName()) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }

      if (oldObj) {
	clone = oldObj ;
        if ( scaleFactor > 0 ) {
           ((TH1*)clone)->Sumw2() ;
           ((TH1*)clone)->Add((TH1*)obj, scaleFactor) ;
        } else {
           ((TH1*)clone)->Add((TH1*)obj) ;
        }
      } else {
	clone = obj->Clone() ;
      }
    }
    if ( scaleFactor > 0 && !doAdd ) {
       ((TH1*) clone)->Sumw2() ;
       ((TH1*) clone)->Scale(scaleFactor) ;
    }
    if (!gDirectory->GetList()->FindObject(clone)) {
      gDirectory->Append(clone) ;
    }
    std::cout << "." ;
    std::cout.flush() ;
  }
  std::cout << std::endl;
  inf.Close() ;
  delete iter ;
}

//==========================================================================================
