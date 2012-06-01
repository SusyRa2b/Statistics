
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"

#include <iostream>
#include <sstream>

  using namespace RooFit;
  using namespace RooStats;


// simple minded PL scan (integrating SUSY signal over all the bins)

void ws_simple_PLscan( TString wsfile = "ws.root" ) {


  // hardcode here the number of bins of the analysis
  
  const int nBinsMET  = 2 ;
  const int nBinsHT   = 1 ;
  const int nBinsBtag = 3 ;

  TString sMbins[nBinsMET];
  TString sHbins[nBinsHT];
  TString sBbins[3] = {"_1b","_2b","_3b"};
  
  for (int i = 0 ; i < nBinsMET ; i++) {
    TString base = "_M";
    stringstream sbin;
    sbin << i+1;
    base += sbin.str();
    sMbins[i] = base;
  }
  
  for (int j = 0 ; j < nBinsHT ; j++) {
    TString base = "_H";
    stringstream sbin;
    sbin << j+1;
    base += sbin.str();
    sHbins[j] = base;
  }
  

  TFile* wstf = new TFile( wsfile ) ;
  RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
  ws->Print() ;

  RooDataSet* rds = (RooDataSet*) ws->obj( "ra2b_observed_rds" ) ;
  cout << "\n\n\n  ===== RooDataSet ====================\n\n" << endl ;
  rds->Print() ;
  rds->printMultiline(cout, 1, kTRUE, "") ;


  // use the SUSY yield in the (0,0,0) bin as poi:

  TString parName = "mu_susy_M1_H1_1b" ;

  cout << "\n\n\n  ===== Grabbing " << parName << " rrv ====================\n\n" ;
  RooRealVar* rrv_par = ws->var( parName ) ;
  if ( rrv_par == 0x0 ) {
    cout << "\n\n\n *** can't find " << parName << " in workspace.  Quitting.\n\n\n" ; return ;
  } else {
    cout << " current value is : " << rrv_par->getVal() << endl ; cout << flush ;
  }


  cout << "\n\n\n  ===== Grabbing likelihood pdf ====================\n\n" ;
  RooAbsPdf* likelihood = ws->pdf("likelihood") ;
  if ( likelihood == 0x0 ) {
    cout << "\n\n\n *** can't find likelihood pdf in workspace.  Quitting.\n\n\n" ; return ;
  } else {
    cout << "\n\n likelihood pdf: \n\n" ;
    likelihood->Print() ;
  }


  cout << "\n\n\n  ===== Doing a fit ====================\n\n" ;

  rrv_par->setConstant(kFALSE);
  likelihood->fitTo( *rds ) ;
  double mlValue = rrv_par->getVal() ;

  RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true) ) ;

  double logLikelihoodSusyFloat = fitResult->minNll() ;
  double logLikelihoodSusyFixed ;

  cout << "\n\nlogLikelihoodSusyFloat = " << logLikelihoodSusyFloat << endl ;

  // just print out the value that the fit prefers for the poi:
  cout << "\n\n  Maximum likelihood value of " << parName << ": " << mlValue << " +/- " <<  rrv_par->getError() << "\n\n" ;

  double testStatVal(0.) ;

  cout << "\n\n\n  ===== Mini scan to find 0.95 CL limit ====================\n\n" ;

  double poival = 0.05 ;

  while ( (2.70 - testStatVal) > 0.01 ) {

    poival = poival*( 1 + 0.15*(2.70 - testStatVal )) ;

    rrv_par->setVal(poival);
    rrv_par->setConstant(kTRUE);

    fitResult = likelihood->fitTo( *rds, Save(true) ) ;
    logLikelihoodSusyFixed = fitResult->minNll() ;
    testStatVal = 2.*(logLikelihoodSusyFixed - logLikelihoodSusyFloat) ;

    cout << "\n\nlogLikelihoodSusyFixed = " << logLikelihoodSusyFixed << endl ;
    cout << " poival = " << poival << "     testStat = " << testStatVal << endl ;

  }

  // ok, found the poi for which the testStat is 2.70
  double poiLimit = poival ;

  cout << "\n\n\n  95% CL limit for poi = " << poiLimit << "\n\n" ;


  // now run a scan between 0 and 120% of the 95% CL limit

  cout << "\n\n\n  ===== Drawing profile likelihood plot ====================\n\n" ;


  // find the fraction of SUSY signal in the (0,0,0) bin

  double mu_susy_tot = 0. ;
  
  for ( int i = 0 ; i < nBinsMET ; i++ ) {
    for ( int j = 0 ; j < nBinsHT ; j++ ) {
      for ( int k = 0 ; k < nBinsBtag ; k++ ) {
	
	TString MuSusyString = "mu_susy" ;
	MuSusyString += sMbins[i]+sHbins[j]+sBbins[k] ;
	
	mu_susy_tot += ((RooRealVar*) ws->obj(MuSusyString))->getVal() ;
	
      }
    }
  }

  double poiFrac = poival / mu_susy_tot ;


  const int nPoints = 12 ; 
  
  double poiScan[nPoints+1];
  double muSusyScan[nPoints+1];
  double testStatScan[nPoints+1] ;
 
  for ( int p = 0 ; p <= nPoints ; p++ ) {

    poiScan[p] = p*poiLimit/10 ;
    muSusyScan[p] = poiScan[p]/poiFrac ;

    rrv_par->setVal(poiScan[p]);
    rrv_par->setConstant(kTRUE);

    fitResult = likelihood->fitTo( *rds, Save(true) ) ;
    logLikelihoodSusyFixed = fitResult->minNll() ;

    testStatScan[p] = 2.*(logLikelihoodSusyFixed - logLikelihoodSusyFloat) ;

    cout << "mu_susy_tot = " << muSusyScan[p] << "     testStat = " << testStatScan[p] << endl ;

  }


  // draw graph

  double totLimit = poiLimit / poiFrac ;

  TCanvas *c0 = new TCanvas("c0","simple PL scan",700,600);

  TLine *lv = new TLine(0,2.70,totLimit,2.70);
  lv->SetLineColor(kGreen+3);
  lv->SetLineWidth(2);

  TLine *lh = new TLine(totLimit,0,totLimit,2.70);
  lh->SetLineColor(kGreen+3);
  lh->SetLineWidth(2);

  TGraph *ScanPlot = new TGraph(nPoints+1,muSusyScan,testStatScan) ;

  stringstream MetBins ; MetBins << nBinsMET ;
  stringstream HtBins  ; HtBins  << nBinsHT ;

  TString grTitle = "PL scan - nBins MET = " ;
  grTitle += MetBins.str() ;
  grTitle += ", nBins HT = " ;
  grTitle += HtBins.str() ;

  ScanPlot->SetTitle(grTitle);
  ScanPlot->GetXaxis()->SetTitle("mu_susy_sig");
  ScanPlot->GetYaxis()->SetTitle("test statistic");

  ScanPlot->SetMarkerStyle(20);
  ScanPlot->SetMarkerColor(2);
  ScanPlot->SetLineColor(2);
  ScanPlot->SetLineWidth(2);

  c0->cd();
  ScanPlot->Draw("APL");
  lv->Draw("same");
  lh->Draw("same");
  c0->cd(0);

  c0->SaveAs("simple_PLscan.gif");

  return ;

}
