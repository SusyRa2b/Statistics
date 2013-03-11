#include "TROOT.h"
#include "TStyle.h"

//CMS style
TStyle *theCmsStyle_ =0;
void initCmsStyle() {

  //check if the style is already defined
  if (theCmsStyle_==0 && gROOT->GetStyle("CMS")==0) {
    theCmsStyle_ = new TStyle("CMS","Style for P-TDR");


    // For the canvas:
    theCmsStyle_->SetCanvasBorderMode(0);
    theCmsStyle_->SetCanvasColor(kWhite);
    theCmsStyle_->SetCanvasDefH(600); //Height of canvas
    theCmsStyle_->SetCanvasDefW(600); //Width of canvas
    theCmsStyle_->SetCanvasDefX(0);   //POsition on screen
    theCmsStyle_->SetCanvasDefY(0);
    
    // For the Pad:
    theCmsStyle_->SetPadBorderMode(0);
    // theCmsStyle_->SetPadBorderSize(Width_t size = 1);
    theCmsStyle_->SetPadColor(kWhite);
    theCmsStyle_->SetPadGridX(false);
    theCmsStyle_->SetPadGridY(false);
    theCmsStyle_->SetGridColor(0);
    theCmsStyle_->SetGridStyle(3);
    theCmsStyle_->SetGridWidth(1);
    
    // For the frame:
    theCmsStyle_->SetFrameBorderMode(0);
    theCmsStyle_->SetFrameBorderSize(1);
    theCmsStyle_->SetFrameFillColor(0);
    theCmsStyle_->SetFrameFillStyle(0);
    theCmsStyle_->SetFrameLineColor(1);
    theCmsStyle_->SetFrameLineStyle(1);
    theCmsStyle_->SetFrameLineWidth(1);
    
    // For the histo:
    // theCmsStyle_->SetHistFillColor(1);
    // theCmsStyle_->SetHistFillStyle(0);
    theCmsStyle_->SetHistLineColor(1);
    theCmsStyle_->SetHistLineStyle(0);
    theCmsStyle_->SetHistLineWidth(1);
    // theCmsStyle_->SetLegoInnerR(Float_t rad = 0.5);
  // theCmsStyle_->SetNumberContours(Int_t number = 20);
    
    theCmsStyle_->SetEndErrorSize(2);
    //  theCmsStyle_->SetErrorMarker(20);
    theCmsStyle_->SetErrorX(0.);
    
    theCmsStyle_->SetMarkerStyle(20);
    
    //For the fit/function:
    theCmsStyle_->SetOptFit(1);
    theCmsStyle_->SetFitFormat("5.4g");
    theCmsStyle_->SetFuncColor(2);
    theCmsStyle_->SetFuncStyle(1);
    theCmsStyle_->SetFuncWidth(1);
    
    //For the date:
    theCmsStyle_->SetOptDate(0);
    // theCmsStyle_->SetDateX(Float_t x = 0.01);
    // theCmsStyle_->SetDateY(Float_t y = 0.01);
    
    // For the statistics box:
    theCmsStyle_->SetOptFile(0);
    theCmsStyle_->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
    theCmsStyle_->SetStatColor(kWhite);
    theCmsStyle_->SetStatFont(42);
    theCmsStyle_->SetStatFontSize(0.025);
    theCmsStyle_->SetStatTextColor(1);
    theCmsStyle_->SetStatFormat("6.4g");
    theCmsStyle_->SetStatBorderSize(1);
    theCmsStyle_->SetStatH(0.1);
    theCmsStyle_->SetStatW(0.15);
    // theCmsStyle_->SetStatStyle(Style_t style = 1001);
    // theCmsStyle_->SetStatX(Float_t x = 0);
    // theCmsStyle_->SetStatY(Float_t y = 0);
    
    // Margins:
    theCmsStyle_->SetPadTopMargin(0.05);
    theCmsStyle_->SetPadBottomMargin(0.13);
    theCmsStyle_->SetPadLeftMargin(0.16);
    theCmsStyle_->SetPadRightMargin(0.02);
    
    // For the Global title:
    
    theCmsStyle_->SetOptTitle(0);
    theCmsStyle_->SetTitleFont(42);
    theCmsStyle_->SetTitleColor(1);
    theCmsStyle_->SetTitleTextColor(1);
    theCmsStyle_->SetTitleFillColor(10);
    theCmsStyle_->SetTitleFontSize(0.05);
    // theCmsStyle_->SetTitleH(0); // Set the height of the title box
    // theCmsStyle_->SetTitleW(0); // Set the width of the title box
    // theCmsStyle_->SetTitleX(0); // Set the position of the title box
    // theCmsStyle_->SetTitleY(0.985); // Set the position of the title box
    // theCmsStyle_->SetTitleStyle(Style_t style = 1001);
    // theCmsStyle_->SetTitleBorderSize(2);
    
    // For the axis titles:
    
    theCmsStyle_->SetTitleColor(1, "XYZ");
    theCmsStyle_->SetTitleFont(42, "XYZ");
    theCmsStyle_->SetTitleSize(0.06, "XYZ");
    // theCmsStyle_->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
    // theCmsStyle_->SetTitleYSize(Float_t size = 0.02);
    theCmsStyle_->SetTitleXOffset(0.9);
    theCmsStyle_->SetTitleYOffset(1.25);
    // theCmsStyle_->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
    
    // For the axis labels:
    
    theCmsStyle_->SetLabelColor(1, "XYZ");
    theCmsStyle_->SetLabelFont(42, "XYZ");
    theCmsStyle_->SetLabelOffset(0.007, "XYZ");
    theCmsStyle_->SetLabelSize(0.05, "XYZ");
    
    // For the axis:
    
    theCmsStyle_->SetAxisColor(1, "XYZ");
    theCmsStyle_->SetStripDecimals(kTRUE);
    theCmsStyle_->SetTickLength(0.03, "XYZ");
    theCmsStyle_->SetNdivisions(510, "XYZ");
    theCmsStyle_->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    theCmsStyle_->SetPadTickY(1);
    
    // Change for log plots:
    theCmsStyle_->SetOptLogx(0);
    theCmsStyle_->SetOptLogy(0);
    theCmsStyle_->SetOptLogz(0);
    
    // Postscript options:
    theCmsStyle_->SetPaperSize(20.,20.);
    // theCmsStyle_->SetLineScalePS(Float_t scale = 3);
    // theCmsStyle_->SetLineStyleString(Int_t i, const char* text);
    // theCmsStyle_->SetHeaderPS(const char* header);
    // theCmsStyle_->SetTitlePS(const char* pstitle);
    
    // theCmsStyle_->SetBarOffset(Float_t baroff = 0.5);
    // theCmsStyle_->SetBarWidth(Float_t barwidth = 0.5);
    // theCmsStyle_->SetPaintTextFormat(const char* format = "g");
    // theCmsStyle_->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    // theCmsStyle_->SetTimeOffset(Double_t toffset);
    // theCmsStyle_->SetHistMinimumZero(kTRUE);
    
    theCmsStyle_->cd(); 
    //end CMS style
    
  }

}
