

#include "TH1F.h"
#include "TStyle.h"
#include "TChain.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TLegend.h"


   void toyplotsge1btight() {

      gStyle->SetPadBottomMargin(0.15) ;
      gStyle->SetPadLeftMargin(0.15) ;
      gStyle->SetOptTitle(0) ;


      TCanvas* c1 = new TCanvas("c1","c1") ;

      TArrow* arrow = new TArrow() ;
      arrow->SetLineWidth(3) ;


      TChain toyns("toytt") ;
      toyns.Add("output-files/toy-data-ge1btight.root") ;


     //--- extract true values from TTree.

      double fit_tru_ttwj ;
      TBranch* b_fit_tru_ttwj ;
      toyns.SetBranchAddress("fit_tru_ttwj", &fit_tru_ttwj, &b_fit_tru_ttwj ) ;
      toyns.GetEntry() ;
      printf("\n\n true value for ttwj : %8.2f\n", fit_tru_ttwj ) ;

      double fit_tru_qcd ;
      TBranch* b_fit_tru_qcd ;
      toyns.SetBranchAddress("fit_tru_qcd", &fit_tru_qcd, &b_fit_tru_qcd ) ;
      toyns.GetEntry() ;
      printf("\n\n true value for qcd : %8.2f\n", fit_tru_qcd ) ;

      double fit_tru_znn ;
      TBranch* b_fit_tru_znn ;
      toyns.SetBranchAddress("fit_tru_znn", &fit_tru_znn, &b_fit_tru_znn ) ;
      toyns.GetEntry() ;
      printf("\n\n true value for znn : %8.2f\n", fit_tru_znn ) ;


     //============ fit values

      gStyle->SetOptStat("emr") ;

      TH1F* httwjfit = new TH1F("httwjfit","ttwj",40, 0., 40. ) ;
      TH1F* hqcdfit  = new TH1F("hqcdfit" ,"qcd" ,40, 0., 15. ) ;
      TH1F* hznnfit  = new TH1F("hznnfit" ,"znn" ,40, 0., 25. ) ;

      toyns.Draw("fit_val_ttwj>>httwjfit","fit_cov_qual==3") ;
      toyns.Draw("fit_val_qcd>>hqcdfit","fit_cov_qual==3") ;
      toyns.Draw("fit_val_znn>>hznnfit","fit_cov_qual==3") ;

      httwjfit->SetLineWidth(2) ;
      hqcdfit->SetLineWidth(2) ;
      hznnfit->SetLineWidth(2) ;

      httwjfit->SetFillColor(11) ;
      hqcdfit->SetFillColor(11) ;
      hznnfit->SetFillColor(11) ;

      httwjfit->SetXTitle("Fit ttwj SIG events") ;
      hqcdfit->SetXTitle("Fit QCD SIG events") ;
      hznnfit->SetXTitle("Fit Znn SIG events") ;

      httwjfit->SetYTitle("Toy experiments") ;
      hqcdfit->SetYTitle("Toy experiments") ;
      hznnfit->SetYTitle("Toy experiments") ;

      httwjfit->Draw() ;
      arrow->DrawArrow(fit_tru_ttwj, 0.4*(httwjfit->GetMaximum()), fit_tru_ttwj, 0) ;
      c1->SaveAs("output-files/toymc-ttwj-sig-fit-ge1btight.png") ;
      hqcdfit->Draw() ;
      arrow->DrawArrow(fit_tru_qcd, 0.4*(hqcdfit->GetMaximum()), fit_tru_qcd, 0) ;
      c1->SaveAs("output-files/toymc-qcd-sig-fit-ge1btight.png") ;
      hznnfit->Draw() ;
      arrow->DrawArrow(fit_tru_znn, 0.4*(hznnfit->GetMaximum()), fit_tru_znn, 0) ;
      c1->SaveAs("output-files/toymc-znn-sig-fit-ge1btight.png") ;


     //============ fit uncertainty

      TH1F* httwjerr = new TH1F("httwjerr","ttwj",40, 0., 15. ) ;
      TH1F* hqcderr  = new TH1F("hqcderr" ,"qcd" ,40, 0., 10. ) ;
      TH1F* hznnerr  = new TH1F("hznnerr" ,"znn" ,40, 0., 15. ) ;

      toyns.Draw("fit_err_ttwj>>httwjerr","fit_cov_qual==3") ;
      toyns.Draw("fit_err_qcd>>hqcderr","fit_cov_qual==3") ;
      toyns.Draw("fit_err_znn>>hznnerr","fit_cov_qual==3") ;

      httwjerr->SetLineWidth(2) ;
      hqcderr->SetLineWidth(2) ;
      hznnerr->SetLineWidth(2) ;

      httwjerr->SetFillColor(11) ;
      hqcderr->SetFillColor(11) ;
      hznnerr->SetFillColor(11) ;

      httwjerr->SetXTitle("err ttwj SIG events") ;
      hqcderr->SetXTitle("err QCD SIG events") ;
      hznnerr->SetXTitle("err Znn SIG events") ;

      httwjerr->SetYTitle("Toy experiments") ;
      hqcderr->SetYTitle("Toy experiments") ;
      hznnerr->SetYTitle("Toy experiments") ;

      httwjerr->Draw() ;
      c1->SaveAs("output-files/toymc-ttwj-sig-err-ge1btight.png") ;
      hqcderr->Draw() ;
      c1->SaveAs("output-files/toymc-qcd-sig-err-ge1btight.png") ;
      hznnerr->Draw() ;
      c1->SaveAs("output-files/toymc-znn-sig-err-ge1btight.png") ;






   }



