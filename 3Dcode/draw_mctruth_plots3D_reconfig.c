
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TAxis.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TText.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TDirectory.h"

#include <string.h>
#include <sstream>

#include <iostream>

 using std::cout ;


   void loadHist1(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;

//----------------

   void draw_mctruth_plots3D_reconfig( const char* histfile = "rootfiles/gi-plots-met4-ht4.root",
                              const char* primaryGroupingVar = "met",
                              const char* secondaryGroupingVar = "ht",
                              bool logy=false,
                              bool doNorm=false,
                              double normmax=2.0,
                              int metgroupzoom=1,
                              bool splitttwj=false ) {





     if ( doNorm ) { logy = false ; }

     gStyle->SetOptStat(0) ;
     gStyle->SetPadTopMargin(0.03) ;
     gStyle->SetPadBottomMargin(0.30) ;
     gStyle->SetPadRightMargin(0.10) ;
     gStyle->SetTitleX(0.95) ;
     gStyle->SetTitleAlign(33) ;

     gDirectory->Delete("h*") ;

     loadHist1( histfile ) ;
     gDirectory->ls() ;

     if ( logy ) {
        gStyle->SetOptLogy(1) ;
     } else {
        gStyle->SetOptLogy(0) ;
     }




     int ncomp(7) ;
     char compname[7][100] = { "allsm", "susy", "ttwj", "qcd", "znn", "ttbar", "wjets" } ;

     int nsel(3) ;
     char selname[3][100] = { "0lep", "1lep", "ldp" } ;

     //int nbtagmult(3) ;
     char btagmultname[3][100] = { "1b", "2b", "3b" } ;

     int compcolor[7] ;
     compcolor[0] = -1 ;
     compcolor[1] = 6 ;
     compcolor[2] = kBlue-9 ;
     compcolor[3] = 2 ;
     compcolor[4] = kGreen-3 ;
     compcolor[5] = kBlue-9 ;
     compcolor[6] = kBlue-7 ;



     //-- figure out the binning from a couple of histograms in the input file.

     int nBinsBjets(3) ;
     TH1F* hmctruth_qcd_lsb_pass = (TH1F*) gDirectory->FindObject("hmctruth_qcd_lsb_pass") ;
     if ( hmctruth_qcd_lsb_pass == 0x0 ) { printf("\n\n *** Can't find hmctruth_qcd_lsb_pass.\n\n") ; return ; }
     int nBinsHT = (hmctruth_qcd_lsb_pass->GetNbinsX() - 1) / nBinsBjets - 1 ;
     printf("\n\n Number of HT bins : %d\n\n", nBinsHT ) ;
     TH1F* hmctruth_qcd_0lep_1b  = (TH1F*) gDirectory->FindObject("hmctruth_qcd_0lep_1b") ;
     if ( hmctruth_qcd_0lep_1b == 0 ) { printf("\n\n*** hmctruth_qcd_0lep_1b missing.\n\n") ; return ; }
     int nBinsMET = (hmctruth_qcd_0lep_1b->GetNbinsX()-1)/(nBinsHT+1) ;
     printf("\n\n Number of MET bins : %d\n\n", nBinsMET ) ;





     //-- read in the input histograms and save pointers.

     TH1F* hinput_[3][3][7] ;
     for ( int si=0; si<nsel; si++ ) {
        for ( int bbi=0 ; bbi<nBinsBjets; bbi++ ) {
           for ( int ci=0; ci<ncomp; ci++ ) {

              char hname[1000] ;

              sprintf( hname, "hmctruth_%s_%s_%s", compname[ci], selname[si], btagmultname[bbi] ) ;
              hinput_[si][bbi][ci] = (TH1F*) gDirectory->FindObject( hname ) ;
              if ( hinput_[si][bbi][ci] == 0x0 ) { printf("\n\n *** can't find histogram %s\n\n", hname ) ; return ; }

           } // ci.
        } // bbi.
     } // si.



     TH1F* hmctruth_fit_zee_1b  =  (TH1F*) gDirectory->FindObject("hmctruth_fit_zee_1b") ;
     TH1F* hmctruth_fit_zmm_1b  =  (TH1F*) gDirectory->FindObject("hmctruth_fit_zmm_1b") ;







     //-- create new histograms with alternate grouping.

      int nbv1(0), nbv2(0), nbv3(0) ;

      char thirdGroupingVar[10] ;

      if ( strcmp( primaryGroupingVar, "ht" ) == 0 && strcmp( secondaryGroupingVar, "met" ) == 0 ) {
         nbv1 = nBinsHT ;
         nbv2 = nBinsMET ;
         nbv3 = nBinsBjets ;
         sprintf( thirdGroupingVar, "nbjet" ) ;
      }
      if ( strcmp( primaryGroupingVar, "ht" ) == 0 && strcmp( secondaryGroupingVar, "nbjet" ) == 0 ) {
         nbv1 = nBinsHT ;
         nbv2 = nBinsBjets ;
         nbv3 = nBinsMET ;
         sprintf( thirdGroupingVar, "met" ) ;
      }
      if ( strcmp( primaryGroupingVar, "met" ) == 0 && strcmp( secondaryGroupingVar, "ht" ) == 0 ) {
         nbv1 = nBinsMET ;
         nbv2 = nBinsHT ;
         nbv3 = nBinsBjets ;
         sprintf( thirdGroupingVar, "nbjet" ) ;
      }
      if ( strcmp( primaryGroupingVar, "met" ) == 0 && strcmp( secondaryGroupingVar, "nbjet" ) == 0 ) {
         nbv1 = nBinsMET ;
         nbv2 = nBinsBjets ;
         nbv3 = nBinsHT ;
         sprintf( thirdGroupingVar, "ht" ) ;
      }
      if ( strcmp( primaryGroupingVar, "nbjet" ) == 0 && strcmp( secondaryGroupingVar, "ht" ) == 0 ) {
         nbv1 = nBinsBjets ;
         nbv2 = nBinsHT ;
         nbv3 = nBinsMET ;
         sprintf( thirdGroupingVar, "met" ) ;
      }
      if ( strcmp( primaryGroupingVar, "nbjet" ) == 0 && strcmp( secondaryGroupingVar, "met" ) == 0 ) {
         nbv1 = nBinsBjets ;
         nbv2 = nBinsMET ;
         nbv3 = nBinsHT ;
         sprintf( thirdGroupingVar, "ht" ) ;
      }

      cout << "\n Checking for canvas.\n\n" << flush ;

      TCanvas* cmctruth = (TCanvas*) gDirectory->FindObject("cmctruth") ;
      if ( cmctruth==0x0 ) {
         printf("\n\n Creating cmctruth canvas.\n\n") ;
         cmctruth = new TCanvas("cmctruth","RA2b fit quality", 850, 1000 ) ;
      } else {
         printf("\n\n Found existing cmctruth canvas.\n\n") ;
         cmctruth->Clear() ;
      }
      cmctruth->Divide( nbv3, 4 ) ;



      int padIndex(1) ;




      if ( nBinsMET > 10 || nBinsHT > 10 ) { printf("\n\n *** too many bins.\n\n") ; return ; }
      TH1F*    hnew_[3][10][7] ; // first index selection, second is variable bin, third is component index.
      THStack* hstack_[3][10] ;

      int mbi(0), hbi(0), bbi(0) ;

      for ( int si=0; si<nsel; si++ ) {

         for ( int v3bi=0; v3bi<nbv3; v3bi++ ) {


            if ( strcmp( thirdGroupingVar, "met"   ) == 0x0 ) { mbi = v3bi ; }
            if ( strcmp( thirdGroupingVar, "ht"    ) == 0x0 ) { hbi = v3bi ; }
            if ( strcmp( thirdGroupingVar, "nbjet" ) == 0x0 ) { bbi = v3bi ; }

            int nbins_new = 1 + (nbv1+1)*nbv2 ;

            char hname[1000] ;
            char htitle[1000] ;


            sprintf( htitle, "%s, %s %d", selname[si], thirdGroupingVar, v3bi+1 ) ;

            printf( "%s\n", htitle ) ;

            for ( int ci=0; ci<ncomp; ci++ ) {

               sprintf( hname, "hnew_%s_%s_%s%d", compname[ci], selname[si], thirdGroupingVar, v3bi+1 ) ;
               hnew_[si][v3bi][ci] = new TH1F( hname, htitle, nbins_new, 0.5, nbins_new+0.5 ) ;
               if ( compcolor[ci] > 0 ) { hnew_[si][v3bi][ci] -> SetFillColor( compcolor[ci] ) ; }

               //-- fill the new histogram.

               for ( int v2bi=0; v2bi<nbv2; v2bi++ ) {

                  if ( strcmp( secondaryGroupingVar, "met"   ) == 0x0 ) { mbi = v2bi ; }
                  if ( strcmp( secondaryGroupingVar, "ht"    ) == 0x0 ) { hbi = v2bi ; }
                  if ( strcmp( secondaryGroupingVar, "nbjet" ) == 0x0 ) { bbi = v2bi ; }

                  for ( int v1bi=0; v1bi<nbv1; v1bi++ ) {

                     if ( strcmp( primaryGroupingVar, "met"   ) == 0x0 ) { mbi = v1bi ; }
                     if ( strcmp( primaryGroupingVar, "ht"    ) == 0x0 ) { hbi = v1bi ; }
                     if ( strcmp( primaryGroupingVar, "nbjet" ) == 0x0 ) { bbi = v1bi ; }

                     int newhbin = 1 + (nbv1+1)*v2bi + v1bi + 1 ;
                     int oldhbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

                     hnew_[si][v3bi][ci] -> SetBinContent( newhbin, hinput_[si][bbi][ci]->GetBinContent( oldhbin ) ) ;
                     hnew_[si][v3bi][ci] -> SetBinError(   newhbin, hinput_[si][bbi][ci]->GetBinError( oldhbin ) ) ;
                     hnew_[si][v3bi][ci] -> GetXaxis() -> SetBinLabel( newhbin, hinput_[si][bbi][ci] -> GetXaxis() -> GetBinLabel( oldhbin ) ) ;

                  } // v1bi

               } // v2bi

               hnew_[si][v3bi][ci] -> SetLabelSize(0.055,"x") ;
               hnew_[si][v3bi][ci] -> GetXaxis() -> LabelsOption("v") ;


            } // ci.

            if ( doNorm ) {
               for ( int v2bi=0; v2bi<nbv2; v2bi++ ) {
                  for ( int v1bi=0; v1bi<nbv1; v1bi++ ) {
                     int newhbin = 1 + (nbv1+1)*v2bi + v1bi + 1 ;
                     double modelTotal(0.) ;
                     for ( int ci=1; ci<=4; ci++ ) {
                        modelTotal += hnew_[si][v3bi][ci] -> GetBinContent( newhbin ) ;
                     } // ci.
                     for ( int ci=0; ci<ncomp; ci++ ) {
                        if ( modelTotal > 0. ) {
                           hnew_[si][v3bi][ci] -> SetBinContent( newhbin, (hnew_[si][v3bi][ci] -> GetBinContent( newhbin ))/modelTotal ) ;
                           hnew_[si][v3bi][ci] -> SetBinError(   newhbin, (hnew_[si][v3bi][ci] -> GetBinError( newhbin ))/modelTotal ) ;
                        } else {
                           hnew_[si][v3bi][ci] -> SetBinContent( newhbin, 0. ) ;
                           hnew_[si][v3bi][ci] -> SetBinError( newhbin, 0. ) ;
                        }
                        hnew_[si][v3bi][ci] -> SetMaximum( normmax ) ;
                     } // ci.
                  } // v1bi
               } // v2bi
            }

            sprintf( hname, "hmctruth_fit_%s_%s%d", selname[si], thirdGroupingVar, v3bi+1 ) ;
            sprintf( htitle, " " ) ;
            hstack_[si][v3bi] = new THStack( hname, htitle ) ;

            if ( splitttwj ) {
               hstack_[si][v3bi] -> Add( hnew_[si][v3bi][4] ) ;
               hstack_[si][v3bi] -> Add( hnew_[si][v3bi][3] ) ;
               hstack_[si][v3bi] -> Add( hnew_[si][v3bi][6] ) ;
               hstack_[si][v3bi] -> Add( hnew_[si][v3bi][5] ) ;
               hstack_[si][v3bi] -> Add( hnew_[si][v3bi][1] ) ;
            } else {
               hstack_[si][v3bi] -> Add( hnew_[si][v3bi][4] ) ;
               hstack_[si][v3bi] -> Add( hnew_[si][v3bi][3] ) ;
               hstack_[si][v3bi] -> Add( hnew_[si][v3bi][2] ) ;
               hstack_[si][v3bi] -> Add( hnew_[si][v3bi][1] ) ;
            }

            cmctruth -> cd( padIndex++ ) ;
            gPad->SetTicks(1,0) ;

            hnew_[si][v3bi][0] -> SetLineWidth(2) ;

            hnew_[si][v3bi][0] -> Draw("histpe") ;
            hstack_[si][v3bi] -> Draw("histsame") ;
            hnew_[si][v3bi][0] -> Draw("samee") ;
            gPad->SetGridy(1) ;

         } // v3bi

      } // si.


     cmctruth->cd(padIndex++) ;
     hmctruth_fit_zee_1b->SetFillColor(kGreen-3) ;
     hmctruth_fit_zee_1b->SetLineWidth(1) ;
     hmctruth_fit_zee_1b->Draw("hist") ;
     hmctruth_fit_zee_1b->Draw("esame") ;
     gPad->SetGridy(1) ;

     cmctruth->cd(padIndex++) ;
     hmctruth_fit_zmm_1b->SetFillColor(kGreen-3) ;
     hmctruth_fit_zmm_1b->SetLineWidth(1) ;
     hmctruth_fit_zmm_1b->Draw("histe") ;
     hmctruth_fit_zmm_1b->Draw("esame") ;
     gPad->SetGridy(1) ;



     printf("\n Making legend...\n") ; cout << flush ;
     TLegend* legend = new TLegend(0.4,0.35,0.7,0.85) ;
     if ( splitttwj ) {
        legend->AddEntry( hnew_[0][0][1], compname[1] ) ;
        legend->AddEntry( hnew_[0][0][5], compname[5] ) ;
        legend->AddEntry( hnew_[0][0][6], compname[6] ) ;
        legend->AddEntry( hnew_[0][0][3], compname[3] ) ;
        legend->AddEntry( hnew_[0][0][4], compname[4] ) ;
     } else {
        legend->AddEntry( hnew_[0][0][1], compname[1] ) ;
        legend->AddEntry( hnew_[0][0][2], compname[2] ) ;
        legend->AddEntry( hnew_[0][0][3], compname[3] ) ;
        legend->AddEntry( hnew_[0][0][4], compname[4] ) ;
     }

     printf("\n\n Done making legend.\n\n\n") ; cout << flush ;

     cmctruth->cd(padIndex++);
     legend->Draw() ;

     cmctruth->Update() ;




     gSystem->Exec("mkdir -p outputfiles") ;

     TString saveFile( histfile ) ;
     saveFile.ReplaceAll( "rootfiles", "outputfiles" ) ;
     char newstring[1000] ;
     sprintf( newstring, "-%s-%s-%s.pdf", primaryGroupingVar, secondaryGroupingVar, thirdGroupingVar ) ;
     saveFile.ReplaceAll( ".root", newstring ) ;

     cmctruth->SaveAs( saveFile ) ;














//   if ( metgroupzoom>1 && !logy ) {

//      int nhistbins = hmctruth_ttwj_0lep_1b->GetNbinsX() ;

//      int zoomRefBin = 1 + (nBinsHT+1)*(metgroupzoom-1) + 1 ;

//      printf("\n\n Number of HT bins=%d, met group ref bin=%d\n\n", nBinsHT, zoomRefBin ) ;

//      double maxSF(1.3) ;

//      hmctruth_allsm_0lep_1b->SetMaximum( maxSF*(hmctruth_allsm_0lep_1b->GetBinContent( zoomRefBin )) ) ;
//      hmctruth_allsm_0lep_2b->SetMaximum( maxSF*(hmctruth_allsm_0lep_2b->GetBinContent( zoomRefBin )) ) ;
//      hmctruth_allsm_0lep_3b->SetMaximum( maxSF*(hmctruth_allsm_0lep_3b->GetBinContent( zoomRefBin )) ) ;

//      hmctruth_allsm_1lep_1b->SetMaximum( maxSF*(hmctruth_allsm_1lep_1b->GetBinContent( zoomRefBin )) ) ;
//      hmctruth_allsm_1lep_2b->SetMaximum( maxSF*(hmctruth_allsm_1lep_2b->GetBinContent( zoomRefBin )) ) ;
//      hmctruth_allsm_1lep_3b->SetMaximum( maxSF*(hmctruth_allsm_1lep_3b->GetBinContent( zoomRefBin )) ) ;

//      hmctruth_allsm_ldp_1b->SetMaximum( maxSF*(hmctruth_allsm_ldp_1b->GetBinContent( zoomRefBin+2 )) ) ;
//      hmctruth_allsm_ldp_2b->SetMaximum( maxSF*(hmctruth_allsm_ldp_2b->GetBinContent( zoomRefBin+2 )) ) ;
//      hmctruth_allsm_ldp_3b->SetMaximum( maxSF*(hmctruth_allsm_ldp_3b->GetBinContent( zoomRefBin+2 )) ) ;

//      hmctruth_fit_zee_1b->SetMaximum( maxSF*(hmctruth_fit_zee_1b->GetBinContent( zoomRefBin )) ) ;
//      hmctruth_fit_zmm_1b->SetMaximum( maxSF*(hmctruth_fit_zmm_1b->GetBinContent( zoomRefBin )) ) ;

//   }











   }



//==========================================================================================


void loadHist1(const char* filename, const char* pfx, const char* pat, Bool_t doAdd, Double_t scaleFactor)
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

