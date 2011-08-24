
#include "TH2.h"
#include "TStyle.h"
#include "TLine.h"
#include "TPad.h"
#include "TText.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TPaletteAxis.h"
#include <iostream>
#include "histio.c"

  using namespace std ;

    void t1bbbb_pl_plots() {

        gStyle->SetPadGridX(1) ;
        gStyle->SetPadGridY(1) ;
        gStyle->SetOptStat(0.) ;
        gStyle->SetOptTitle(0) ;
        gStyle->SetPadRightMargin(0.20) ;
        gStyle->SetPadLeftMargin(0.15) ;
        gStyle->SetPadBottomMargin(0.15) ;
        gStyle->SetTitleOffset(1.4,"x") ;
        gStyle->SetTitleOffset(1.4,"y") ;
        gStyle->SetBarOffset(0.1) ;

        TText* text = new TText()  ;

        TCanvas* cxsec = new TCanvas("cxsec","cxsec") ;

        gPad->SetLogz(1) ;

        TPaletteAxis* paxis(0x0) ;


        loadHist("an-11-257-v3-files/output-files/an-scanplot-unblind-t1bbbb-withcontam-ge1b-loose.root","ge1bloose" ) ;
        loadHist("an-11-257-v3-files/output-files/an-scanplot-unblind-t1bbbb-withcontam-ge1b-tight.root","ge1btight" ) ;
        loadHist("an-11-257-v3-files/output-files/an-scanplot-unblind-t1bbbb-withcontam-ge2b-loose.root","ge2bloose" ) ;
        loadHist("an-11-257-v3-files/output-files/an-scanplot-unblind-t1bbbb-withcontam-ge2b-tight.root","ge2btight" ) ;



    //--------
        TH2F* hxsecul_ge1bloose = (TH2F*) gDirectory->FindObject("ge1bloose_hsusyscanXsecul") ;

        hxsecul_ge1bloose->SetXTitle("Gluino mass (GeV)") ;
        hxsecul_ge1bloose->SetYTitle("LSP mass (GeV)") ;
        hxsecul_ge1bloose->UseCurrentStyle() ;
        hxsecul_ge1bloose->SetMinimum(0.010) ;
        hxsecul_ge1bloose->SetMaximum(150) ;

        paxis = (TPaletteAxis*)hxsecul_ge1bloose->GetListOfFunctions()->FindObject("palette");
        paxis->SetX1NDC(0.82) ;
        paxis->SetX2NDC(0.86) ;
        paxis->SetY1NDC(0.15) ;
        paxis->SetY2NDC(0.90) ;
        hxsecul_ge1bloose->DrawCopy("colz") ;

        text->SetTextAngle(90.) ;
        text->SetTextSize(0.04) ;
        text->DrawTextNDC( 0.95, 0.15, "Maximum cross section [pb](PL 95% UL)" ) ;

        text->SetTextAngle(0.) ;
        text->SetTextSize(0.045) ;
        text->DrawTextNDC( 0.15, 0.92, "T1bbbb: ge1bloose, PL cross section upper limit") ;


        cxsec->SaveAs("an-11-257-v3-files/output-files/an-t1bbbb-pl-xsecul-ge1bloose.png") ;


    //--------

        TH2F* hxsecul_ge1btight = (TH2F*) gDirectory->FindObject("ge1btight_hsusyscanXsecul") ;

        hxsecul_ge1btight->SetXTitle("Gluino mass (GeV)") ;
        hxsecul_ge1btight->SetYTitle("LSP mass (GeV)") ;
        hxsecul_ge1btight->UseCurrentStyle() ;
        hxsecul_ge1btight->SetMinimum(0.010) ;
        hxsecul_ge1btight->SetMaximum(150) ;

        paxis = (TPaletteAxis*)hxsecul_ge1btight->GetListOfFunctions()->FindObject("palette");
        paxis->SetX1NDC(0.82) ;
        paxis->SetX2NDC(0.86) ;
        paxis->SetY1NDC(0.15) ;
        paxis->SetY2NDC(0.90) ;
        hxsecul_ge1btight->DrawCopy("colz") ;

        text->SetTextAngle(90.) ;
        text->SetTextSize(0.04) ;
        text->DrawTextNDC( 0.95, 0.15, "Maximum cross section [pb](PL 95% UL)" ) ;

        text->SetTextAngle(0.) ;
        text->SetTextSize(0.045) ;
        text->DrawTextNDC( 0.15, 0.92, "T1bbbb: ge1btight, PL cross section upper limit") ;


        cxsec->SaveAs("an-11-257-v3-files/output-files/an-t1bbbb-pl-xsecul-ge1btight.png") ;


    //--------

        TH2F* hxsecul_ge2bloose = (TH2F*) gDirectory->FindObject("ge2bloose_hsusyscanXsecul") ;

        hxsecul_ge2bloose->SetXTitle("Gluino mass (GeV)") ;
        hxsecul_ge2bloose->SetYTitle("LSP mass (GeV)") ;
        hxsecul_ge2bloose->UseCurrentStyle() ;
        hxsecul_ge2bloose->SetMinimum(0.010) ;
        hxsecul_ge2bloose->SetMaximum(150) ;

        paxis = (TPaletteAxis*)hxsecul_ge2bloose->GetListOfFunctions()->FindObject("palette");
        paxis->SetX1NDC(0.82) ;
        paxis->SetX2NDC(0.86) ;
        paxis->SetY1NDC(0.15) ;
        paxis->SetY2NDC(0.90) ;
        hxsecul_ge2bloose->DrawCopy("colz") ;

        text->SetTextAngle(90.) ;
        text->SetTextSize(0.04) ;
        text->DrawTextNDC( 0.95, 0.15, "Maximum cross section [pb](PL 95% UL)" ) ;

        text->SetTextAngle(0.) ;
        text->SetTextSize(0.045) ;
        text->DrawTextNDC( 0.15, 0.92, "T1bbbb: ge2bloose, PL cross section upper limit") ;


        cxsec->SaveAs("an-11-257-v3-files/output-files/an-t1bbbb-pl-xsecul-ge2bloose.png") ;


    //--------

        TH2F* hxsecul_ge2btight = (TH2F*) gDirectory->FindObject("ge2btight_hsusyscanXsecul") ;


        hxsecul_ge2btight->SetXTitle("Gluino mass (GeV)") ;
        hxsecul_ge2btight->SetYTitle("LSP mass (GeV)") ;
        hxsecul_ge2btight->UseCurrentStyle() ;
        hxsecul_ge2btight->SetMinimum(0.010) ;
        hxsecul_ge2btight->SetMaximum(150) ;

        paxis = (TPaletteAxis*)hxsecul_ge2btight->GetListOfFunctions()->FindObject("palette");
        paxis->SetX1NDC(0.82) ;
        paxis->SetX2NDC(0.86) ;
        paxis->SetY1NDC(0.15) ;
        paxis->SetY2NDC(0.90) ;
        hxsecul_ge2btight->DrawCopy("colz") ;

        text->SetTextAngle(90.) ;
        text->SetTextSize(0.04) ;
        text->DrawTextNDC( 0.95, 0.15, "Maximum cross section [pb](PL 95% UL)" ) ;

        text->SetTextAngle(0.) ;
        text->SetTextSize(0.045) ;
        text->DrawTextNDC( 0.15, 0.92, "T1bbbb: ge2btight, PL cross section upper limit") ;


        cxsec->SaveAs("an-11-257-v3-files/output-files/an-t1bbbb-pl-xsecul-ge2btight.png") ;













        TCanvas* ceff = new TCanvas("ceff","ceff") ;

        gPad->SetLogz(0) ;

    //--------
        TH2F* heff_ge1bloose = (TH2F*) gDirectory->FindObject("ge1bloose_hsusyscanEfficiency") ;

        cout << "heff_ge1bloose pointer: " <<  heff_ge1bloose << endl ;

        heff_ge1bloose->SetXTitle("Gluino mass (GeV)") ;
        heff_ge1bloose->SetYTitle("LSP mass (GeV)") ;
        heff_ge1bloose->UseCurrentStyle() ;
        heff_ge1bloose->SetMinimum(0.0) ;
        heff_ge1bloose->SetMaximum(0.6) ;

        heff_ge1bloose->DrawCopy("colz") ;
        heff_ge1bloose->DrawCopy("colz") ;

        text->SetTextAngle(90.) ;
        text->SetTextSize(0.04) ;
        text->DrawTextNDC( 0.95, 0.15, "Efficiency" ) ;

        text->SetTextAngle(0.) ;
        text->SetTextSize(0.045) ;
        text->DrawTextNDC( 0.15, 0.92, "T1bbbb: ge1bloose, Efficiency") ;

        ceff->SaveAs("an-11-257-v3-files/output-files/an-t1bbbb-eff-ge1bloose.png") ;

    //--------
        TH2F* heff_ge1btight = (TH2F*) gDirectory->FindObject("ge1btight_hsusyscanEfficiency") ;

        cout << "heff_ge1btight pointer: " <<  heff_ge1btight << endl ;

        heff_ge1btight->SetXTitle("Gluino mass (GeV)") ;
        heff_ge1btight->SetYTitle("LSP mass (GeV)") ;
        heff_ge1btight->UseCurrentStyle() ;
        heff_ge1btight->SetMinimum(0.0) ;
        heff_ge1btight->SetMaximum(0.6) ;

        heff_ge1btight->DrawCopy("colz") ;
        heff_ge1btight->DrawCopy("colz") ;

        text->SetTextAngle(90.) ;
        text->SetTextSize(0.04) ;
        text->DrawTextNDC( 0.95, 0.15, "Efficiency" ) ;

        text->SetTextAngle(0.) ;
        text->SetTextSize(0.045) ;
        text->DrawTextNDC( 0.15, 0.92, "T1bbbb: ge1btight, Efficiency") ;

        ceff->SaveAs("an-11-257-v3-files/output-files/an-t1bbbb-eff-ge1btight.png") ;

    //--------
        TH2F* heff_ge2bloose = (TH2F*) gDirectory->FindObject("ge2bloose_hsusyscanEfficiency") ;

        cout << "heff_ge2bloose pointer: " <<  heff_ge2bloose << endl ;

        heff_ge2bloose->SetXTitle("Gluino mass (GeV)") ;
        heff_ge2bloose->SetYTitle("LSP mass (GeV)") ;
        heff_ge2bloose->UseCurrentStyle() ;
        heff_ge2bloose->SetMinimum(0.0) ;
        heff_ge2bloose->SetMaximum(0.6) ;

        heff_ge2bloose->DrawCopy("colz") ;
        heff_ge2bloose->DrawCopy("colz") ;

        text->SetTextAngle(90.) ;
        text->SetTextSize(0.04) ;
        text->DrawTextNDC( 0.95, 0.15, "Efficiency" ) ;

        text->SetTextAngle(0.) ;
        text->SetTextSize(0.045) ;
        text->DrawTextNDC( 0.15, 0.92, "T1bbbb: ge2bloose, Efficiency") ;

        ceff->SaveAs("an-11-257-v3-files/output-files/an-t1bbbb-eff-ge2bloose.png") ;

    //--------
        TH2F* heff_ge2btight = (TH2F*) gDirectory->FindObject("ge2btight_hsusyscanEfficiency") ;

        cout << "heff_ge2btight pointer: " <<  heff_ge2btight << endl ;

        heff_ge2btight->SetXTitle("Gluino mass (GeV)") ;
        heff_ge2btight->SetYTitle("LSP mass (GeV)") ;
        heff_ge2btight->UseCurrentStyle() ;
        heff_ge2btight->SetMinimum(0.0) ;
        heff_ge2btight->SetMaximum(0.6) ;

        heff_ge2btight->DrawCopy("colz") ;
        heff_ge2btight->DrawCopy("colz") ;

        text->SetTextAngle(90.) ;
        text->SetTextSize(0.04) ;
        text->DrawTextNDC( 0.95, 0.15, "Efficiency" ) ;

        text->SetTextAngle(0.) ;
        text->SetTextSize(0.045) ;
        text->DrawTextNDC( 0.15, 0.92, "T1bbbb: ge2btight, Efficiency") ;

        ceff->SaveAs("an-11-257-v3-files/output-files/an-t1bbbb-eff-ge2btight.png") ;

    }




