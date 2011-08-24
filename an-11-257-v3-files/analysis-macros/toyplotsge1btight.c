

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


      TChain toyns("tt_toy_nosusyfit") ;
      toyns.Add("an-11-257-v3-files/output-files/toyge1btight-smonly-mctest.root") ;

      TChain bgo("tt_cls_bgonly") ;
      TChain spb("tt_cls_splusb") ;

      bgo.Add("an-11-257-v3-files/output-files/toyge1btight-smonly-mctest.root") ;
      spb.Add("an-11-257-v3-files/output-files/toyge1btight-smonly-mctest.root") ;



     //============ fit values

      gStyle->SetOptStat("emr") ;

      TH1F* httwjfit = new TH1F("httwjfit","ttwj",40, 0., 45. ) ;
      TH1F* hqcdfit  = new TH1F("hqcdfit" ,"qcd" ,40, 0., 10. ) ;
      TH1F* hznnfit  = new TH1F("hznnfit" ,"znn" ,40, 0., 20. ) ;

      toyns.Draw("ttwj_sig_fit>>httwjfit","") ;
      toyns.Draw("qcd_sig_fit>>hqcdfit","") ;
      toyns.Draw("znn_sig_fit>>hznnfit","") ;

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
      arrow->DrawArrow(19.6,40,19.6,0) ;
      c1->SaveAs("an-11-257-v3-files/output-files/toymc-ttwj-sig-fit-ge1btight.png") ;
      hqcdfit->Draw() ;
      arrow->DrawArrow(1.3,60,1.3,0) ;
      c1->SaveAs("an-11-257-v3-files/output-files/toymc-qcd-sig-fit-ge1btight.png") ;
      hznnfit->Draw() ;
      arrow->DrawArrow(4.25,40,4.25,0) ;
      c1->SaveAs("an-11-257-v3-files/output-files/toymc-znn-sig-fit-ge1btight.png") ;


     //============ fit uncertainty

      TH1F* httwjerr = new TH1F("httwjerr","ttwj",40, 0., 10. ) ;
      TH1F* hqcderr  = new TH1F("hqcderr" ,"qcd" ,40, 0., 8. ) ;
      TH1F* hznnerr  = new TH1F("hznnerr" ,"znn" ,40, 0., 10. ) ;

      toyns.Draw("ttwj_sig_err>>httwjerr","") ;
      toyns.Draw("qcd_sig_err>>hqcderr","") ;
      toyns.Draw("znn_sig_err>>hznnerr","") ;

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
      c1->SaveAs("an-11-257-v3-files/output-files/toymc-ttwj-sig-err-ge1btight.png") ;
      hqcderr->Draw() ;
      c1->SaveAs("an-11-257-v3-files/output-files/toymc-qcd-sig-err-ge1btight.png") ;
      hznnerr->Draw() ;
      c1->SaveAs("an-11-257-v3-files/output-files/toymc-znn-sig-err-ge1btight.png") ;


     //============ q value distributions.

      gStyle->SetOptStat(0) ;

      TH1F* hbgo = new TH1F("hbgo","BG only", 50, 0., 15.) ;
      TH1F* hspb = new TH1F("hspb","SIG + BG", 50, 0., 15.) ;

      hbgo->SetLineWidth(2) ;
      hspb->SetLineWidth(2) ;
      hbgo->SetLineColor(2) ;
      hspb->SetLineColor(4) ;

      bgo.Draw("testStat>>hbgo","") ;
      spb.Draw("testStat>>hspb","") ;

      hspb->SetXTitle("q value") ;
      hspb->SetYTitle("Toy experiments") ;

      hspb->Draw() ;
      hbgo->Draw("same") ;


    //+++ X value below is hardwired to 2.44.
    //    If inputs are changed, grep the string "Data value of test statistic"
    //    in the runtoyge1btight log file to find the new value.
    //
      arrow->DrawArrow(1.13,150,1.13,0.,0,">") ;

      TLegend* legend = new TLegend(0.4,0.5,0.8,0.8) ;
      legend->AddEntry( hbgo, "BG-only hypothesis") ;
      legend->AddEntry( hspb, "SUSY+BG hypothesis") ;
      legend->SetFillColor(kWhite) ;

      legend->Draw() ;


      c1->SaveAs("an-11-257-v3-files/output-files/toymc-qvalue-distributions-ge1btight.png") ;




   }



