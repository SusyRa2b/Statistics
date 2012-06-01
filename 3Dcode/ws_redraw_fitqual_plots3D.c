
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


void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;

//----------------

   void ws_redraw_fitqual_plots3D( const char* histfile = "fitqual-hists-ws-met3-ht3-v1.root", bool logy=false ) {

     gStyle->SetOptStat(0) ;
     gStyle->SetPadTopMargin(0.03) ;
     gStyle->SetPadBottomMargin(0.30) ;
     gStyle->SetPadRightMargin(0.10) ;
     gStyle->SetTitleX(0.95) ;
     gStyle->SetTitleAlign(33) ;

     gDirectory->Delete("*") ;

     loadHist( histfile ) ;
     gDirectory->ls() ;

     if ( logy ) {
        gStyle->SetOptLogy(1) ;
     } else {
        gStyle->SetOptLogy(0) ;
     }


     TH1F* hfitqual_data_0lep_1b = (TH1F*) gDirectory->FindObject("hfitqual_data_0lep_1b") ;
     TH1F* hfitqual_susy_0lep_1b = (TH1F*) gDirectory->FindObject("hfitqual_susy_0lep_1b") ;
     TH1F* hfitqual_ttwj_0lep_1b = (TH1F*) gDirectory->FindObject("hfitqual_ttwj_0lep_1b") ;
     TH1F* hfitqual_qcd_0lep_1b  = (TH1F*) gDirectory->FindObject("hfitqual_qcd_0lep_1b") ;
     TH1F* hfitqual_znn_0lep_1b  = (TH1F*) gDirectory->FindObject("hfitqual_znn_0lep_1b") ;
     if ( hfitqual_data_0lep_1b == 0 ) { printf("\n\n*** hfitqual_data_0lep_1b missing.\n\n") ; return ; }
     if ( hfitqual_susy_0lep_1b == 0 ) { printf("\n\n*** hfitqual_susy_0lep_1b missing.\n\n") ; return ; }
     if ( hfitqual_ttwj_0lep_1b == 0 ) { printf("\n\n*** hfitqual_ttwj_0lep_1b missing.\n\n") ; return ; }
     if ( hfitqual_qcd_0lep_1b == 0 ) { printf("\n\n*** hfitqual_qcd_0lep_1b missing.\n\n") ; return ; }
     if ( hfitqual_znn_0lep_1b == 0 ) { printf("\n\n*** hfitqual_znn_0lep_1b missing.\n\n") ; return ; }

     TH1F* hfitqual_data_0lep_2b = (TH1F*) gDirectory->FindObject("hfitqual_data_0lep_2b") ;
     TH1F* hfitqual_susy_0lep_2b = (TH1F*) gDirectory->FindObject("hfitqual_susy_0lep_2b") ;
     TH1F* hfitqual_ttwj_0lep_2b = (TH1F*) gDirectory->FindObject("hfitqual_ttwj_0lep_2b") ;
     TH1F* hfitqual_qcd_0lep_2b  = (TH1F*) gDirectory->FindObject("hfitqual_qcd_0lep_2b") ;
     TH1F* hfitqual_znn_0lep_2b  = (TH1F*) gDirectory->FindObject("hfitqual_znn_0lep_2b") ;
     if ( hfitqual_data_0lep_2b == 0 ) { printf("\n\n*** hfitqual_data_0lep_2b missing.\n\n") ; return ; }
     if ( hfitqual_susy_0lep_2b == 0 ) { printf("\n\n*** hfitqual_susy_0lep_2b missing.\n\n") ; return ; }
     if ( hfitqual_ttwj_0lep_2b == 0 ) { printf("\n\n*** hfitqual_ttwj_0lep_2b missing.\n\n") ; return ; }
     if ( hfitqual_qcd_0lep_2b == 0 ) { printf("\n\n*** hfitqual_qcd_0lep_2b missing.\n\n") ; return ; }
     if ( hfitqual_znn_0lep_2b == 0 ) { printf("\n\n*** hfitqual_znn_0lep_2b missing.\n\n") ; return ; }

     TH1F* hfitqual_data_0lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_data_0lep_3b") ;
     TH1F* hfitqual_susy_0lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_0lep_3b") ;
     TH1F* hfitqual_ttwj_0lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_0lep_3b") ;
     TH1F* hfitqual_qcd_0lep_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_0lep_3b") ;
     TH1F* hfitqual_znn_0lep_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_0lep_3b") ;
     if ( hfitqual_data_0lep_3b == 0 ) { printf("\n\n*** hfitqual_data_0lep_3b missing.\n\n") ; return ; }
     if ( hfitqual_susy_0lep_3b == 0 ) { printf("\n\n*** hfitqual_susy_0lep_3b missing.\n\n") ; return ; }
     if ( hfitqual_ttwj_0lep_3b == 0 ) { printf("\n\n*** hfitqual_ttwj_0lep_3b missing.\n\n") ; return ; }
     if ( hfitqual_qcd_0lep_3b == 0 ) { printf("\n\n*** hfitqual_qcd_0lep_3b missing.\n\n") ; return ; }
     if ( hfitqual_znn_0lep_3b == 0 ) { printf("\n\n*** hfitqual_znn_0lep_3b missing.\n\n") ; return ; }



     TH1F* hfitqual_data_1lep_1b =  (TH1F*) gDirectory->FindObject("hfitqual_data_1lep_1b") ;
     TH1F* hfitqual_susy_1lep_1b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_1lep_1b") ;
     TH1F* hfitqual_ttwj_1lep_1b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_1lep_1b") ;
     TH1F* hfitqual_qcd_1lep_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_1lep_1b") ;
     TH1F* hfitqual_znn_1lep_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_1lep_1b") ;

     TH1F* hfitqual_data_1lep_2b =  (TH1F*) gDirectory->FindObject("hfitqual_data_1lep_2b") ;
     TH1F* hfitqual_susy_1lep_2b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_1lep_2b") ;
     TH1F* hfitqual_ttwj_1lep_2b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_1lep_2b") ;
     TH1F* hfitqual_qcd_1lep_2b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_1lep_2b") ;
     TH1F* hfitqual_znn_1lep_2b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_1lep_2b") ;

     TH1F* hfitqual_data_1lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_data_1lep_3b") ;
     TH1F* hfitqual_susy_1lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_1lep_3b") ;
     TH1F* hfitqual_ttwj_1lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_1lep_3b") ;
     TH1F* hfitqual_qcd_1lep_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_1lep_3b") ;
     TH1F* hfitqual_znn_1lep_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_1lep_3b") ;



     TH1F* hfitqual_data_ldp_1b =  (TH1F*) gDirectory->FindObject("hfitqual_data_ldp_1b") ;
     TH1F* hfitqual_susy_ldp_1b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_ldp_1b") ;
     TH1F* hfitqual_ttwj_ldp_1b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_ldp_1b") ;
     TH1F* hfitqual_qcd_ldp_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_ldp_1b") ;
     TH1F* hfitqual_znn_ldp_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_ldp_1b") ;

     TH1F* hfitqual_data_ldp_2b =  (TH1F*) gDirectory->FindObject("hfitqual_data_ldp_2b") ;
     TH1F* hfitqual_susy_ldp_2b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_ldp_2b") ;
     TH1F* hfitqual_ttwj_ldp_2b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_ldp_2b") ;
     TH1F* hfitqual_qcd_ldp_2b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_ldp_2b") ;
     TH1F* hfitqual_znn_ldp_2b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_ldp_2b") ;

     TH1F* hfitqual_data_ldp_3b =  (TH1F*) gDirectory->FindObject("hfitqual_data_ldp_3b") ;
     TH1F* hfitqual_susy_ldp_3b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_ldp_3b") ;
     TH1F* hfitqual_ttwj_ldp_3b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_ldp_3b") ;
     TH1F* hfitqual_qcd_ldp_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_ldp_3b") ;
     TH1F* hfitqual_znn_ldp_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_ldp_3b") ;


     TH1F* hfitqual_data_zee_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_data_zee_1b") ;
     TH1F* hfitqual_data_zmm_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_data_zmm_1b") ;
     TH1F* hfitqual_fit_zee_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_fit_zee_1b") ;
     TH1F* hfitqual_fit_zmm_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_fit_zmm_1b") ;

     TH1F* hfitqual_np   =  (TH1F*) gDirectory->FindObject("hfitqual_np") ;


     printf("\n\n Making stacks...\n") ; cout << flush ;

     THStack* hfitqual_fit_0lep_1b = new THStack( "hfitqual_fit_0lep_1b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_0lep_2b = new THStack( "hfitqual_fit_0lep_2b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_0lep_3b = new THStack( "hfitqual_fit_0lep_3b", "RA2b likelihood fit results, fit" ) ;

     THStack* hfitqual_fit_1lep_1b = new THStack( "hfitqual_fit_1lep_1b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_1lep_2b = new THStack( "hfitqual_fit_1lep_2b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_1lep_3b = new THStack( "hfitqual_fit_1lep_3b", "RA2b likelihood fit results, fit" ) ;

     THStack* hfitqual_fit_ldp_1b  = new THStack( "hfitqual_fit_ldp_1b",  "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_ldp_2b  = new THStack( "hfitqual_fit_ldp_2b",  "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_ldp_3b  = new THStack( "hfitqual_fit_ldp_3b",  "RA2b likelihood fit results, fit" ) ;

     hfitqual_fit_0lep_1b->Add( hfitqual_znn_0lep_1b ) ;
     hfitqual_fit_0lep_1b->Add( hfitqual_qcd_0lep_1b ) ;
     hfitqual_fit_0lep_1b->Add( hfitqual_ttwj_0lep_1b ) ;
     hfitqual_fit_0lep_1b->Add( hfitqual_susy_0lep_1b ) ;

     hfitqual_fit_0lep_2b->Add( hfitqual_znn_0lep_2b ) ;
     hfitqual_fit_0lep_2b->Add( hfitqual_qcd_0lep_2b ) ;
     hfitqual_fit_0lep_2b->Add( hfitqual_ttwj_0lep_2b ) ;
     hfitqual_fit_0lep_2b->Add( hfitqual_susy_0lep_2b ) ;

     hfitqual_fit_0lep_3b->Add( hfitqual_znn_0lep_3b ) ;
     hfitqual_fit_0lep_3b->Add( hfitqual_qcd_0lep_3b ) ;
     hfitqual_fit_0lep_3b->Add( hfitqual_ttwj_0lep_3b ) ;
     hfitqual_fit_0lep_3b->Add( hfitqual_susy_0lep_3b ) ;



     hfitqual_fit_1lep_1b->Add( hfitqual_znn_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_qcd_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_ttwj_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_susy_1lep_1b ) ;

     hfitqual_fit_1lep_2b->Add( hfitqual_znn_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_qcd_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_ttwj_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_susy_1lep_2b ) ;

     hfitqual_fit_1lep_3b->Add( hfitqual_znn_1lep_3b ) ;
     hfitqual_fit_1lep_3b->Add( hfitqual_qcd_1lep_3b ) ;
     hfitqual_fit_1lep_3b->Add( hfitqual_ttwj_1lep_3b ) ;
     hfitqual_fit_1lep_3b->Add( hfitqual_susy_1lep_3b ) ;




     hfitqual_fit_ldp_1b->Add( hfitqual_znn_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_qcd_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_ttwj_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_susy_ldp_1b ) ;

     hfitqual_fit_ldp_2b->Add( hfitqual_znn_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_qcd_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_ttwj_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_susy_ldp_2b ) ;

     hfitqual_fit_ldp_3b->Add( hfitqual_znn_ldp_3b ) ;
     hfitqual_fit_ldp_3b->Add( hfitqual_qcd_ldp_3b ) ;
     hfitqual_fit_ldp_3b->Add( hfitqual_ttwj_ldp_3b ) ;
     hfitqual_fit_ldp_3b->Add( hfitqual_susy_ldp_3b ) ;



     printf("\n\n Done making stacks.\n\n") ; cout << flush ;

     printf("\n Making legend...\n") ; cout << flush ;
     TLegend* legend = new TLegend(0.4,0.35,0.7,0.85) ;

     legend->AddEntry( hfitqual_data_0lep_1b, "data" ) ;
     legend->AddEntry( hfitqual_susy_0lep_1b, "SUSY" ) ;
     legend->AddEntry( hfitqual_ttwj_0lep_1b, "ttwj" ) ;
     legend->AddEntry( hfitqual_qcd_0lep_1b,  "QCD" ) ;
     legend->AddEntry( hfitqual_znn_0lep_1b,  "Znunu" ) ;
     legend->AddEntry( hfitqual_np,           "Eff PG" ) ;

     printf("\n\n Done making legend.\n\n\n") ; cout << flush ;

     TCanvas* cfitqual = new TCanvas("cfitqual","RA2b fit quality", 850, 1000 ) ;

     cfitqual->Divide(3,4);

     gPad->SetTicks(1,0) ;


     printf(" pad 1\n") ; cout << flush ;
     cfitqual->cd(1);
     hfitqual_data_0lep_1b->Draw("histpe") ;
     hfitqual_fit_0lep_1b->Draw("same") ;
     hfitqual_data_0lep_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 2\n") ; cout << flush ;
     cfitqual->cd(2);
     hfitqual_data_0lep_2b->Draw("histpe") ;
     hfitqual_fit_0lep_2b->Draw("same") ;
     hfitqual_data_0lep_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 3\n") ; cout << flush ;
     cfitqual->cd(3);
     hfitqual_data_0lep_3b->Draw("histpe") ;
     hfitqual_fit_0lep_3b->Draw("same") ;
     hfitqual_data_0lep_3b->Draw("same") ;
     gPad->SetGridy(1) ;
     


     printf(" pad 4\n") ; cout << flush ;
     cfitqual->cd(4);
     hfitqual_data_1lep_1b->Draw("histpe") ;
     hfitqual_fit_1lep_1b->Draw("same") ;
     hfitqual_data_1lep_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 5\n") ; cout << flush ;
     cfitqual->cd(5);
     hfitqual_data_1lep_2b->Draw("histpe") ;
     hfitqual_fit_1lep_2b->Draw("same") ;
     hfitqual_data_1lep_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 6\n") ; cout << flush ;
     cfitqual->cd(6);
     hfitqual_data_1lep_3b->Draw("histpe") ;
     hfitqual_fit_1lep_3b->Draw("same") ;
     hfitqual_data_1lep_3b->Draw("same") ;
     gPad->SetGridy(1) ;




     
     printf(" pad 7\n") ; cout << flush ;
     cfitqual->cd(7);
     hfitqual_data_ldp_1b->Draw("histpe") ;
     hfitqual_fit_ldp_1b->Draw("same") ;
     hfitqual_data_ldp_1b->Draw("same") ;
     gPad->SetGridy(1) ;

     printf(" pad 8\n") ; cout << flush ;
     cfitqual->cd(8);
     hfitqual_data_ldp_2b->Draw("histpe") ;
     hfitqual_fit_ldp_2b->Draw("same") ;
     hfitqual_data_ldp_2b->Draw("same") ;
     gPad->SetGridy(1) ;

     printf(" pad 9\n") ; cout << flush ;
     cfitqual->cd(9);
     hfitqual_data_ldp_3b->Draw("histpe") ;
     hfitqual_fit_ldp_3b->Draw("same") ;
     hfitqual_data_ldp_3b->Draw("same") ;
     gPad->SetGridy(1) ;




     cfitqual->cd(10) ;
     hfitqual_data_zee_1b->Draw("histpe") ;
     hfitqual_fit_zee_1b->Draw("same") ;
     hfitqual_data_zee_1b->Draw("same") ;
     gPad->SetGridy(1) ;

     cfitqual->cd(11) ;
     hfitqual_data_zmm_1b->Draw("histpe") ;
     hfitqual_fit_zmm_1b->Draw("same") ;
     hfitqual_data_zmm_1b->Draw("same") ;
     gPad->SetGridy(1) ;




     cfitqual->cd(12);
     legend->Draw() ;

     cfitqual->Update() ;



     //--- Efficiency scale factor primary Gaussian value.

     hfitqual_np->SetMinimum(-5.) ;
     hfitqual_np->SetMaximum( 5.) ;
     
     hfitqual_np->SetNdivisions(101,"x") ;
     hfitqual_np->SetNdivisions(101,"y") ;
     hfitqual_np->SetLabelOffset(99,"y") ;
     

     TPad* tp = new TPad("tp","tp",0.09,0.,0.18,1.0) ;

     tp->SetRightMargin(0.4) ;

     tp->Draw() ;
     tp->cd() ;
     hfitqual_np->SetLabelSize(0.5,"x") ;
     TAxis *xaxis ;
     xaxis = hfitqual_np->GetXaxis() ;
     xaxis->SetBinLabel(1,"Eff PG") ;
     hfitqual_np->GetXaxis()->LabelsOption("v") ;
     hfitqual_np->Draw() ;

     cfitqual->Update() ;


     TGaxis* axis = new TGaxis() ;
     axis->SetLabelOffset(0.1) ;
     axis->SetLabelSize(0.30) ;
     axis->SetTickSize(0.2) ;
     axis->DrawAxis( 1.0, -5., 1.0, 5., -5., 5., 510, "+LS") ;
     
     cfitqual->Update() ;

     cfitqual->SaveAs("fitqual.gif") ;




   }



//==========================================================================================


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

