
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"


#include <iostream>

   void saveHist(const char* filename, const char* pat) ;
   void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;
   TH1F* bookHist(const char* hname, const char* htitle, int nBinsMET, int nBinsHT, int sampleIndex ) ;

   //-------

   void qcd_study1( bool fillHists=false, bool savePlots=false ) {

      TLine* line = new TLine() ;

      gStyle->SetPalette(1) ;
      gStyle->SetPadBottomMargin(0.20) ;

      gDirectory->Delete("h*") ;

      int nBinsBjets = 3 ;

      const int nBinsMET = 3 ;
      const int nBinsHT  = 3 ;
      float Mbins[nBinsMET+1] = { 125, 200,  350, 99999. } ;
      float Hbins[nBinsHT+1]  = { 400, 600, 1000, 99999. } ;

      TH2F* hdummy = new TH2F("hdummy","",2, Mbins[0], 1500., 2, Hbins[0], 1500. ) ;

      const int nQcdSamples(9) ;

      char inputfile[9][1000] = {
        "files15fb_8TeV/nn-qcd-0120-to-0170.root"
       ,"files15fb_8TeV/nn-qcd-0170-to-0300.root"
       ,"files15fb_8TeV/nn-qcd-0300-to-0470.root"
       ,"files15fb_8TeV/nn-qcd-0470-to-0600.root"
       ,"files15fb_8TeV/nn-qcd-0600-to-0800.root"
       ,"files15fb_8TeV/nn-qcd-0800-to-1000.root"
       ,"files15fb_8TeV/nn-qcd-1000-to-1400.root"
       ,"files15fb_8TeV/nn-qcd-1400-to-1800.root"
       ,"files15fb_8TeV/nn-qcd-1800.root"
      } ;

      char samplename[9][100] = {
        "qcd_0120_to_0170"
       ,"qcd_0170_to_0300"
       ,"qcd_0300_to_0470"
       ,"qcd_0470_to_0600"
       ,"qcd_0600_to_0800"
       ,"qcd_0800_to_1000"
       ,"qcd_1000_to_1400"
       ,"qcd_1400_to_1800"
       ,"qcd_1800_to_9999"
      } ;

 //   int samplecolor[9] = {
 //      kOrange+8,
 //      kRed,
 //      KMagenta+2,
 //      KBlue+1,
 //      KAzure+2,
 //      KGreen+2,
 //      kPink-2,
 //      kViolet+2,
 //      kOrange+2
 //   } ;

      int samplecolor[9] = {
         800+8,
         632,
         616+2,
         600+1,
         860+2,
         416+2,
         900-2,
         880+2,
         800+2
      } ;

      TH2F*   h0lep[nQcdSamples][nBinsBjets] ;
      TH2F*   hldp [nQcdSamples][nBinsBjets] ;

      TCanvas* cqcd = (TCanvas*) gDirectory->FindObject("cqcd") ;
      if ( cqcd == 0x0 ) {
         cqcd = new TCanvas("cqcd", "qcd study", 700, 900 ) ;
      }


      if ( fillHists ) {

         TChain* qcdch[nQcdSamples] ;

         printf("\n\n") ;
         for ( int si=0; si<nQcdSamples; si++ ) {

            qcdch[si] = new TChain("tree") ;
            printf(" %2d : connecting to %s\n", si, inputfile[si] ) ;
            qcdch[si] -> Add( inputfile[si] ) ;

            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               char hname[1000] ;
               char htitle[1000] ;
               sprintf( hname, "h_0lep_%db_%s", bbi+1, samplename[si] ) ;
               sprintf( htitle, "QCD 0lep yield, nb=%d, %s", bbi+1, samplename[si] ) ;
               printf("         booking hist %s : %s\n", hname, htitle ) ;
               h0lep[si][bbi] = new TH2F( hname, htitle, nBinsMET, Mbins, nBinsHT, Hbins ) ;
               h0lep[si][bbi] -> Sumw2() ;
               sprintf( hname, "h_ldp_%db_%s", bbi+1, samplename[si] ) ;
               sprintf( htitle, "QCD  LDP yield, nb=%d, %s", bbi+1, samplename[si] ) ;
               printf("         booking hist %s  : %s\n", hname, htitle ) ;
               hldp [si][bbi] = new TH2F( hname, htitle, nBinsMET, Mbins, nBinsHT, Hbins ) ;
               hldp [si][bbi] -> Sumw2() ;

            } // bbi.

         } // si.
         printf("\n\n") ;





         char basecuts[10000] ;
         sprintf( basecuts, "pt_1st_leadJet>50&&pt_2nd_leadJet>50&&pt_3rd_leadJet>50&&HT>400&&MET>125&&nMu==0&&nEl==0&&nB>=1&&nJets>=3" ) ;

         char bcut[3][100] = { "nB==1", "nB==2", "nB>=3" } ;


         for ( int si=0; si<nQcdSamples; si++ ) {

            printf(" %2d : %s : 0lep\n", si, samplename[si] ) ; cout << flush ;
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               char arg1[1000] ;

               char cuts0lep[10000] ;
               sprintf( cuts0lep, "(%s)&&(minDelPhiN>4&&%s)", basecuts, bcut[bbi] ) ;
               printf("     %db, 0lep cuts : %s\n", bbi+1, cuts0lep ) ;
               sprintf( arg1, "HT:MET>>h_0lep_%db_%s", bbi+1, samplename[si] ) ;
               qcdch[si] -> Draw( arg1, cuts0lep ) ;
               hdummy->Draw() ;
               h0lep[si][bbi]->Draw("samecolz") ;
               cqcd->Update() ; cqcd->Draw() ;


               char cutsldp[10000] ;
               sprintf( cutsldp, "(%s)&&(minDelPhiN<=4&&%s)", basecuts, bcut[bbi] ) ;
               printf("     %db, ldp  cuts : %s\n", bbi+1, cutsldp  ) ;
               sprintf( arg1, "HT:MET>>h_ldp_%db_%s", bbi+1, samplename[si] ) ;
               qcdch[si] -> Draw( arg1, cutsldp, "colz" ) ;
               hdummy->Draw() ;
               hldp[si][bbi]->Draw("samecolz") ;
               cqcd->Update() ; cqcd->Draw() ;

            // hldp[si][bbi]->Print("all") ;
            // getchar() ;

            } // bbi.

         } // si.

         saveHist( "rootfiles/qcd-study1.root", "h*" ) ;

      }

      //--------

      gStyle->SetOptStat(0) ;

      gDirectory->Delete("h*") ;
      loadHist( "rootfiles/qcd-study1.root" ) ;

      for ( int si=0; si<nQcdSamples; si++ ) {
         for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

            char hname[1000] ;
            sprintf( hname, "h_0lep_%db_%s", bbi+1, samplename[si] ) ;
            printf("  loading %s\n", hname ) ;
            h0lep[si][bbi] = (TH2F*) gDirectory->FindObject( hname ) ;
            if ( h0lep[si][bbi] == 0x0 ) { printf("\n\n *** %s missing.\n\n", hname ) ; return ; }
            sprintf( hname, "h_ldp_%db_%s", bbi+1, samplename[si] ) ;
            printf("  loading %s\n", hname ) ;
            hldp [si][bbi] = (TH2F*) gDirectory->FindObject( hname ) ;
            if ( hldp [si][bbi] == 0x0 ) { printf("\n\n *** %s missing.\n\n", hname ) ; return ; }

         } // bbi.
      } // si.

      int nbins = nBinsMET*(nBinsHT+1) + 1 ;

      TH1F* hflat_0lep[nQcdSamples][nBinsBjets] ;
      TH1F* hflat_ldp [nQcdSamples][nBinsBjets] ;
      TH1F* hflat_0lepldp_ratio[nQcdSamples][nBinsBjets] ;

      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
         for ( int si=0; si<nQcdSamples; si++ ) {

            char hname[1000] ;
            char htitle[1000] ;

            sprintf( hname, "hflat_0lep_%db_%s", bbi+1, samplename[si] ) ;
            sprintf( htitle, "QCD 0lep events, nb=%d, %s", bbi+1, samplename[si] ) ;
            hflat_0lep[si][bbi] = bookHist( hname, htitle, nBinsMET, nBinsHT, si ) ;
            hflat_0lep[si][bbi]->SetFillColor(11) ;
            hldp[si][bbi]->Print("all") ;

            sprintf( hname, "hflat_ldp_%db_%s", bbi+1, samplename[si] ) ;
            sprintf( htitle, "QCD LDP events, nb=%d, %s", bbi+1, samplename[si] ) ;
            hflat_ldp[si][bbi] = bookHist( hname, htitle, nBinsMET, nBinsHT, si ) ;
            hflat_ldp[si][bbi]->SetFillColor(11) ;

            sprintf( hname, "hflat_0lepldp_ratio_%db_%s", bbi+1, samplename[si] ) ;
            sprintf( htitle, "QCD 0lep/LDP ratio, nb=%d, %s", bbi+1, samplename[si] ) ;
            hflat_0lepldp_ratio[si][bbi] = bookHist( hname, htitle, nBinsMET, nBinsHT, si ) ;

            for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
               for ( int hbi=0; hbi<nBinsHT; hbi++ ) {

                  int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

                  float n0lep     = h0lep[si][bbi] -> GetBinContent( mbi+1, hbi+1 ) ;
                  float n0lep_err = h0lep[si][bbi] -> GetBinError( mbi+1, hbi+1 ) ;

                  float nldp      = hldp [si][bbi] -> GetBinContent( mbi+1, hbi+1 ) ;
                  float nldp_err  = hldp [si][bbi] -> GetBinError( mbi+1, hbi+1 ) ;

                  hflat_0lep[si][bbi] -> SetBinContent( histbin, n0lep ) ;
                  hflat_0lep[si][bbi] -> SetBinError( histbin, n0lep_err ) ;

                  hflat_ldp[si][bbi]  -> SetBinContent( histbin, nldp ) ;
                  hflat_ldp[si][bbi]  -> SetBinError( histbin, nldp_err ) ;

                  double ratio(0.) ;
                  if ( nldp > 0 ) { ratio = n0lep / nldp ; }
                  double ratio_err(0.) ;
                  if ( n0lep_err > 0. && nldp_err > 0. ) {
                     ratio_err = ratio * sqrt( pow(n0lep_err/n0lep,2) + pow(nldp_err/nldp,2) ) ;
                  }

                  hflat_0lepldp_ratio[si][bbi] -> SetBinContent( histbin, ratio ) ;
                  hflat_0lepldp_ratio[si][bbi] -> SetBinError( histbin, ratio_err ) ;


               } // hbi.
            } // mbi.

            cqcd->Clear() ;
            cqcd->Divide(1,3) ;

            cqcd->cd(1) ;
            hflat_0lep[si][bbi] -> Draw("histe") ;
            hflat_0lep[si][bbi] -> Draw("samee") ;

            cqcd->cd(2) ;
            hflat_ldp[si][bbi] -> Draw("histe") ;
            hflat_ldp[si][bbi] -> Draw("samee") ;

            cqcd->cd(3) ;
            hflat_0lepldp_ratio[si][bbi]->SetMinimum(-0.1) ;
            hflat_0lepldp_ratio[si][bbi]->SetMaximum(0.6) ;
            hflat_0lepldp_ratio[si][bbi]->SetMarkerStyle(20) ;
            hflat_0lepldp_ratio[si][bbi]->SetLineWidth(2) ;
            hflat_0lepldp_ratio[si][bbi]->Draw() ;
            gPad->SetGridx(1) ;
            gPad->SetGridy(1) ;
            line->DrawLine(0.5,0,nbins+0.5,0) ;


            cqcd->Update() ; cqcd->Draw() ;

            // char a = getchar() ;
            // if ( a == 'q') { return ; }

         } // si.
      } // bbi.




      TH1F* hflat_0lepldp_ratio_ave[nBinsBjets] ;

      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

         char hname[1000] ;
         char htitle[1000] ;

         sprintf( hname, "hflat_0lepldp_ratio_ave_%db", bbi+1 ) ;
         sprintf( htitle, "QCD 0lep/LDP average ratio, nb=%d", bbi+1 ) ;
         hflat_0lepldp_ratio_ave[bbi] = bookHist( hname, htitle, nBinsMET, nBinsHT, 5 ) ;

         hflat_0lepldp_ratio_ave[bbi] -> SetMarkerStyle(24) ;
         hflat_0lepldp_ratio_ave[bbi] -> SetMarkerSize(2.0) ;
         hflat_0lepldp_ratio_ave[bbi] -> SetLineWidth(3) ;

         for ( int binind=1; binind<=nbins; binind++ ) {

            double wvsum(0.) ;
            double wsum(0.) ;

            for ( int si=0; si<nQcdSamples; si++ ) {

               double val = hflat_0lepldp_ratio[si][bbi] -> GetBinContent( binind ) ;
               double err = hflat_0lepldp_ratio[si][bbi] -> GetBinError( binind ) ;
               if ( err > 0 ) {
                  double weight = 1./(err*err) ;
                  wvsum += val*weight ;
                  wsum  += weight ;
               }
            } // si.

            if ( wsum > 0. ) {
               hflat_0lepldp_ratio_ave[bbi] -> SetBinContent( binind, wvsum/wsum ) ;
               hflat_0lepldp_ratio_ave[bbi] -> SetBinError( binind, 1./sqrt(wsum) ) ;
            }

         } // binind

      } // bbi.







     //---

      TCanvas* cqcd2 = (TCanvas*) gDirectory->FindObject("cqcd2") ;
      if ( cqcd2 == 0x0 ) {
         cqcd2 = new TCanvas("cqcd2", "qcd study", 1700, 600 ) ;
      }

      TLegend* l2 = new TLegend( 0.79, 0.70,  0.99, 0.99 ) ;
      l2->AddEntry( hflat_0lepldp_ratio[0][0], "120 to 170" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[1][0], "170 to 300" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[2][0], "300 to 470" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[3][0], "470 to 600" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[4][0], "600 to 800" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[5][0], "800 to 1000" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[6][0], "1000 to 1400" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[7][0], "1400 to 1800" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[8][0], "> 1800" ) ;
      l2->AddEntry( hflat_0lepldp_ratio_ave[0], "Average ratio" ) ;


      cqcd2 -> Clear() ;
      cqcd2 -> Divide(3,1) ;

      cqcd2 -> cd(1) ;
      hflat_0lepldp_ratio[0][0] -> SetTitle("QCD 0lep/LDP ratio, nb=1") ;
      for ( int si=0; si<nQcdSamples; si++ ) {
         hflat_0lepldp_ratio[si][0] -> SetMarkerStyle(20+si) ;
         hflat_0lepldp_ratio[si][0] -> SetLineColor(samplecolor[si]) ;
         hflat_0lepldp_ratio[si][0] -> SetMarkerColor(samplecolor[si]) ;
         if ( si == 0 ) { hflat_0lepldp_ratio[si][0] -> DrawCopy() ; } else { hflat_0lepldp_ratio[si][0] -> Draw("same") ; }
      }
      hflat_0lepldp_ratio_ave[0]->Draw("same") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      l2->Draw() ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd2 -> cd(2) ;
      hflat_0lepldp_ratio[0][1] -> SetTitle("QCD 0lep/LDP ratio, nb=2") ;
      for ( int si=0; si<nQcdSamples; si++ ) {
         hflat_0lepldp_ratio[si][1] -> SetMarkerStyle(20+si) ;
         hflat_0lepldp_ratio[si][1] -> SetLineColor(samplecolor[si]) ;
         hflat_0lepldp_ratio[si][1] -> SetMarkerColor(samplecolor[si]) ;
         if ( si == 0 ) { hflat_0lepldp_ratio[si][1] -> DrawCopy() ; } else { hflat_0lepldp_ratio[si][1] -> Draw("same") ; }
      }
      hflat_0lepldp_ratio_ave[1]->Draw("same") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      l2->Draw() ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd2 -> cd(3) ;
      hflat_0lepldp_ratio[0][2] -> SetTitle("QCD 0lep/LDP ratio, nb>=3") ;
      for ( int si=0; si<nQcdSamples; si++ ) {
         hflat_0lepldp_ratio[si][2] -> SetMarkerStyle(20+si) ;
         hflat_0lepldp_ratio[si][2] -> SetLineColor(samplecolor[si]) ;
         hflat_0lepldp_ratio[si][2] -> SetMarkerColor(samplecolor[si]) ;
         if ( si == 0 ) { hflat_0lepldp_ratio[si][2] -> DrawCopy() ; } else { hflat_0lepldp_ratio[si][2] -> Draw("same") ; }
      }
      hflat_0lepldp_ratio_ave[2]->Draw("same") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      l2->Draw() ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;


      cqcd2->Update() ; cqcd2->Draw() ;

      if ( savePlots ) {
         cqcd2 -> SaveAs("outputfiles/qcd-study-allonone.png") ;
         cqcd2 -> SaveAs("outputfiles/qcd-study-allonone.pdf") ;
      }

     //---

      TCanvas* cqcd3 = (TCanvas*) gDirectory->FindObject("cqcd3") ;
      if ( cqcd3 == 0x0 ) {
         cqcd3 = new TCanvas("cqcd3", "qcd study", 1200, 900 ) ;
      }

      char oldtitle[1000] ;
      sprintf( oldtitle, "QCD 0lep/LDP ratio, nb=1, %s", samplename[0] ) ;
      hflat_0lepldp_ratio[0][0] -> SetTitle(oldtitle) ;
      sprintf( oldtitle, "QCD 0lep/LDP ratio, nb=2, %s", samplename[0] ) ;
      hflat_0lepldp_ratio[0][1] -> SetTitle(oldtitle) ;
      sprintf( oldtitle, "QCD 0lep/LDP ratio, nb>=3, %s", samplename[0] ) ;
      hflat_0lepldp_ratio[0][2] -> SetTitle(oldtitle) ;

      cqcd3->Clear() ;
      cqcd3->Divide(3,3) ;

      for ( int si=0; si<nQcdSamples; si++ ) {
         TLegend* l3 = new TLegend( 0.65, 0.78, 0.99, 0.92) ;
         l3->AddEntry( hflat_0lepldp_ratio[si][0], samplename[si] ) ;
         l3->AddEntry( hflat_0lepldp_ratio_ave[0], "Average") ;
         for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

            cqcd3->cd(1+bbi) ;
            hflat_0lep[si][bbi] -> Draw("histe") ;
            hflat_0lep[si][bbi] -> Draw("samee") ;

            cqcd3->cd(4+bbi) ;
            hflat_ldp[si][bbi] -> Draw("histe") ;
            hflat_ldp[si][bbi] -> Draw("samee") ;


            cqcd3->cd(7+bbi) ;
            gPad->SetGridx(1) ;
            gPad->SetGridy(1) ;
            hflat_0lepldp_ratio[si][bbi]->Draw() ;
            hflat_0lepldp_ratio_ave[bbi]->Draw("same") ;
            line->DrawLine(0.5,0,nbins+0.5,0) ;
            l3->Draw() ;

         } // bbi.

         cqcd3->Update() ; cqcd3->Draw() ;

         if ( savePlots ) {
            char filename[1000] ;
            sprintf( filename, "outputfiles/qcd-study-%s.pdf", samplename[si] ) ;
            cqcd3 -> SaveAs( filename ) ;
            sprintf( filename, "outputfiles/qcd-study-%s.png", samplename[si] ) ;
            cqcd3 -> SaveAs( filename ) ;
         }
         char a = getchar() ;
         if ( a == 'q') { return ; }

      } // si.


   } // qcd_study.

  //==========================================================================================
  
  
    TH1F* bookHist(const char* hname, const char* htitle, int nBinsMET, int nBinsHT, int sampleIndex ) {
  
       int nbins = nBinsMET*(nBinsHT+1) + 1 ;
  
       TH1F* retVal = new TH1F( hname, htitle, nbins, 0.5+0.05*(sampleIndex-5), nbins+0.5+0.05*(sampleIndex-5) ) ;
       TAxis* xaxis = retVal->GetXaxis() ;
  
       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
             int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
             char binlabel[1000] ;
             sprintf( binlabel, "M%d_H%d", mbi+1, hbi+1 ) ;
             xaxis->SetBinLabel( histbin, binlabel ) ;
          } // hbi.
       } // mbi.
  
       retVal->SetLabelSize(0.055,"x") ;
       xaxis->LabelsOption("v") ;
  
       return retVal ;
  
    }
  
  //==========================================================================================


void saveHist(const char* filename, const char* pat)
{

  cout << "\n\n Saving histograms matching " << pat << " in file " << filename << "\n\n" << flush ;

  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();

  TRegexp re(pat,kTRUE) ;

  TFile outf(filename,"RECREATE") ;
  TObject* obj ;
  while((obj=iter->Next())) {
    if (TString(obj->GetName()).Index(re)>=0) {
      obj->Write() ;
      std::cout << "." ;
    }
  }
  std::cout << std::endl ;
  outf.Close() ;

  delete iter ;
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


