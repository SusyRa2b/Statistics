{

  //  gROOT->ProcessLine(".x setup.C");

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
  

  
  //////////////////////////////////////////////////////////////////////////
  // THIS GREATLY SIMPLIFIES THE TAU->HAD PORTION

  gROOT->ProcessLine(".L metReweightingBuilderSIMPLETAU.C++");

  TString wspacefile("testoutputMORE_xsec1.root");



  // SETUP WORKSPACE IN FRESH FILE
  RooWorkspace*  wspace = new RooWorkspace("wspace");


  // SETUP SIGNAL CROSS-SECTION
  RooRealVar signalCrossSection("signalCrossSection","signalCrossSection",0.01,0.,10.);
  wspace->import(signalCrossSection);
  wspace->defineSet("namesfordata","");
  wspace->defineSet("nuisances","");
  //wspace->writeToFile( wspacefile.Data(), true ); // RECREATE ROOT FILE

  /*
  // SET UP TRIGGER EFFICIENCIES
  RooAbsArg* zeroLeptonTriggerEfficiency =
    getGaussianConstraint(ws, "zeroLeptonTriggerEfficiency_", abcd.zeroLeptonTriggerEfficiencyName, //FIXME change to beta after synch
			  abcd.zeroLeptonTriggerEfficiency, abcd.zeroLeptonTriggerEfficiencyError,
			  names.observables,names.nuisances);
  ws.import(*zeroLeptonTriggerEfficiency, RecycleConflictNodes());

  RooAbsArg* oneLeptonTriggerEfficiency =
    getGaussianConstraint(ws, "oneLeptonTriggerEfficiency_", abcd.oneLeptonTriggerEfficiencyName, //FIXME change to beta after synch
			  abcd.oneLeptonTriggerEfficiency, abcd.oneLeptonTriggerEfficiencyError,
			  names.observables,names.nuisances);
  ws.import(*oneLeptonTriggerEfficiency, RecycleConflictNodes());
  */


    
  //buildMRLikelihood( wspacefile.Data(), "testSimpleSetupFileHT3.txt" );
  buildMRLikelihood( *wspace, wspacefile.Data(), "testDataMRSetupFileSIMPLE.txt", true );

  wspace->writeToFile( outputFile.Data(), true ); 


   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetHistFillColor(kWhite);
  gStyle->SetCanvasColor(kWhite);     // background is no longer mouse-dropping white
  gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
  gStyle->SetPadBorderMode(0);


  gROOT->SetStyle("Plain");      // Switches off the ROOT default style
  //  gPad->UseCurrentStyle();       // this makes everything black and white,
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);


  
  TFile *test = new TFile( wspacefile.Data(), "READ" );
  RooWorkspace *wspace2 = (RooWorkspace*) test->Get("wspace");
  wspace2->autoImportClassCode(true);


  /*

  //  (wspace2->obj("zeroLepton_bin333_SignalYield")).Print();

  RooArgSet parofinterest(*wspace->set("poi")) ;
  //  RooArgSet parofinterest(*wspace2->set("oneTightMu_bin333_Theta4_TopWJetsYield"));
 
  ProfileLikelihoodCalculator plc(*wspace2->data("dataset"), *wspace2->pdf("model"), 
				  parofinterest);
				  //*wspace2->var("oneTightMu_bin333_Theta4_TopWJetsYield"));

  plc.SetConfidenceLevel(0.95);
  LikelihoodInterval* interval = plc.GetInterval() ;
  LikelihoodIntervalPlot*  lplot = new LikelihoodIntervalPlot(interval);

  double profileLikelihoodUpperLimit = interval->UpperLimit(*wspace2->var("signalCrossSection");
  double profileLikelihoodLowerLimit = interval->LowerLimit(*wspace2->var("signalCrossSection");
  firstPOI->setRange(profileLikelihoodLowerLimit,profileLikelihoodUpperLimit);

  TCanvas* canvas1 = new TCanvas("Result", "Result", 10, 10, 500, 500);
  lplot->Draw();
  canvas1->Update();
  canvas1->SaveAs("PLC_xsec2_MORE.eps");

  */



  RooAbsReal* nll = (wspace2->pdf("model")).createNLL(*wspace2->data("dataset"));
  RooAbsReal* pll1 = nll->createProfile(*wspace2->var("signalCrossSection"));
  //RooAbsReal* pll1 = nll->createProfile(*wspace2->var("oneTightMu_bin333_Theta4_TopWJetsYield"));
  //RooAbsReal* pll2 = nll->createProfile(*wspace2->var("oneLooseLep_bin333_Theta4_TopWJetsYield"));
  //RooAbsReal* pll3 = nll->createProfile(*wspace2->var("twoTightMu_bin333_TopWJetsYield"));
  //RooAbsReal* pll4 = nll->createProfile(*wspace2->var("twoLooseLep_bin333_TopWJetsYield"));

  //  (*wspace2->var("twoTightMu_bin333_TopWJetsYield")).Print();
  //  RooPlot* frame1 = signalCrossSection.frame(0,5);//yield.frame();
  //  RooPlot* frame1 = (*wspace2->var("oneLooseLep_bin301_Theta4_TopWJetsYield")).frame(1,20);//yield.frame();
  
  RooPlot* frame1 = (*wspace2->var("signalCrossSection")).frame(0,5);//yield.frame();
  //RooPlot* frame1 = (*wspace2->var("oneTightMu_bin333_Theta4_TopWJetsYield")).frame(0,2);//yield.frame();
  //RooPlot* frame2 = (*wspace2->var("oneLooseLep_bin333_Theta4_TopWJetsYield")).frame(0,2);//yield.frame();
  //RooPlot* frame3 = (*wspace2->var("twoTightMu_bin333_TopWJetsYield")).frame(0,2);//yield.frame();
  //RooPlot* frame4 = (*wspace2->var("twoLooseLep_bin333_TopWJetsYield")).frame(0,2);//yield.frame();
  //  nll->plotOn(frame1) ;


  TCanvas *c1 =  new TCanvas();
  pll1->plotOn(frame1) ;
  frame1->Draw();
  //  c1->SaveAs("PLL_xsec2_1tmu4.eps");
  c1->SaveAs("PLL_xsec1_XSEC.eps");

  /*
  TCanvas *c2 =  new TCanvas();
  pll2->plotOn(frame2) ;
  frame2->Draw();
  c2->SaveAs("PLL_xsec2_1lep4.eps");

  TCanvas *c3 =  new TCanvas();
  pll3->plotOn(frame3) ;
  frame3->Draw();
  c3->SaveAs("PLL_xsec2_2tmu.eps");

  TCanvas *c4 =  new TCanvas();
  pll4->plotOn(frame4) ;
  frame4->Draw();
  c4->SaveAs("PLL_xsec2_2lep.eps");

  */


  // crude plotting of all NBs/MET for a slice of HT

  TH1F *bkg_pol_0L = new TH1F("bkg_pol_0L","",20,0,20);
  TH1F *bkg_tau_0L = new TH1F("bkg_tau_0L","",20,0,20);
  TH1F *bkg_dil_0L = new TH1F("bkg_dil_0L","",20,0,20);
  TH1F *background_1L = new TH1F("background_1L","",60,0,60);

  TH1F *signal_0L = new TH1F("signal_0L","",20,0,20);
  TH1F *signal_1L = new TH1F("signal_1L","",60,0,60);

  TH1F *observed_0L = new TH1F("observed_0L","",20,0,20);
  TH1F *observed_1L = new TH1F("observed_1L","",60,0,60);

  TH1F *truesig_0L = new TH1F("truesig_0L","",20,0,20);
  TH1F *truesig_1L = new TH1F("truesig_1L","",60,0,60);
  TH1F *truebkg_0L = new TH1F("truebkg_0L","",20,0,20);
  TH1F *truebkg_1L = new TH1F("truebkg_1L","",60,0,60);


  cout << " Before met and nb loops " << endl; 

  for( int metbin=0; metbin<4; metbin++ ){

    for( int nb=1; nb<4; nb++  ){

      //      if( (metbin==0 && nb==1) || (metbin==3 && nb==3)  ){
	//if( (metbin==0 && nb==1) ){

    TString zeroLeptonName("zeroLepton_bin3");
    zeroLeptonName+=metbin;
    zeroLeptonName+=nb;
    TString oneTightMuName("oneTightMu_bin3");
    oneTightMuName+=metbin;
    oneTightMuName+=nb;
    TString oneLooseLepName("oneLooseLep_bin3");
    oneLooseLepName+=metbin;
    oneLooseLepName+=nb;
    TString twoTightMuName("twoTightMu_bin3");
    twoTightMuName+=metbin;
    twoTightMuName+=nb;
    TString twoLooseLepName("twoLooseLep_bin3");
    twoLooseLepName+=metbin;
    twoLooseLepName+=nb;

    cout << " After names: " << zeroLeptonName << "  " << oneTightMuName << "  " << oneLooseLepName << "  " 
	 << twoTightMuName << "  " << twoLooseLepName << endl;

    double tofill = 0.;
    int fillbin = 0;

    ///////////////////////////////////////////////////////////////////////////////////
    // SL observation

    TString name1(oneTightMuName+"_Theta1_Count");
    TString name2(oneTightMuName+"_Theta2_Count");
    TString name3(oneTightMuName+"_Theta3_Count");
    TString name4(oneTightMuName+"_Theta4_Count");
    TString name5(oneTightMuName+"_Theta5_Count");
    
    tofill = (wspace2->var(name1.Data()))->getVal() + (wspace2->var(name2.Data()))->getVal() 
      + (wspace2->var(name3.Data()))->getVal() + (wspace2->var(name4.Data()))->getVal() + (wspace2->var(name5.Data()))->getVal();

    cout << " observed: " << tofill << "  " << (wspace2->var(name1.Data()))->getVal() << "  " 
	 << (wspace2->var(name2.Data()))->getVal() << "  " << (wspace2->var(name3.Data()))->getVal() << "  " 
	 << (wspace2->var(name4.Data()))->getVal() <<"  " << (wspace2->var(name5.Data()))->getVal() <<endl;
    fillbin = (metbin*15 + 1) + (nb-1)*5;
    //    if( metbin==3 ) fillbin = fillbin-5;
    observed_1L->Fill( fillbin, tofill );

    //////
    TString name1a(oneLooseLepName+"_Theta1_Count");
    TString name2a(oneLooseLepName+"_Theta2_Count");
    TString name3a(oneLooseLepName+"_Theta3_Count");
    TString name4a(oneLooseLepName+"_Theta4_Count");
    TString name5a(oneLooseLepName+"_Theta5_Count");

    tofill = (wspace2->var(name1a.Data()))->getVal() + (wspace2->var(name2a.Data()))->getVal()
      + (wspace2->var(name3a.Data()))->getVal() + (wspace2->var(name4a.Data()))->getVal() + (wspace2->var(name5a.Data()))->getVal();
    fillbin = (metbin*15 + 1) + (nb-1)*5 + 1;
    //if( metbin==3 ) fillbin = fillbin-5;
    observed_1L->Fill( fillbin, tofill );

    //////
    TString name6(twoTightMuName+"_Count");
    TString name7(twoLooseLepName+"_Count");

    tofill = (wspace2->var(name6.Data()))->getVal();
    fillbin = (metbin*15 + 1) + (nb-1)*5 + 2;
    //if( metbin==3 ) fillbin = fillbin-5;
    observed_1L->Fill( fillbin, tofill );

    tofill = (wspace2->var(name7.Data()))->getVal();
    fillbin = (metbin*15 + 1) + (nb-1)*5 + 3;
    //if( metbin==3 ) fillbin = fillbin-5;
    observed_1L->Fill( fillbin, tofill );

    // end of SL observation

    cout << " after SL observation " << endl; 

    ///////////////////////////////////////////////////////////////////////////////////
    // SL signal

    TString names1(oneTightMuName+"_Theta1_SignalYield");
    TString names2(oneTightMuName+"_Theta2_SignalYield");
    TString names3(oneTightMuName+"_Theta3_SignalYield");
    TString names4(oneTightMuName+"_Theta4_SignalYield");
    TString names5(oneTightMuName+"_Theta5_SignalYield");

    cout << names1.Data() << "  " << names2.Data() << " " <<names3.Data() << " " <<names4.Data() << " " <<names5.Data() << endl;
    
    tofill = (wspace2->function(names1.Data()))->getVal() + (wspace2->function(names2.Data()))->getVal() 
      + (wspace2->function(names3.Data()))->getVal() + (wspace2->function(names4.Data()))->getVal() + (wspace2->function(names5.Data()))->getVal();
    fillbin = (metbin*15 + 1) + (nb-1)*5;
    //if( metbin==3 ) fillbin = fillbin-5;
    signal_1L->Fill( fillbin, tofill );
    //if( metbin==0 ) truesig_1L->Fill( fillbin+5, 0.3*5. );
    //if( metbin==3 ) truesig_1L->Fill( fillbin+5, 2.6*5. );


    cout << "  tight mu sig " << endl; 
    //////
    TString names1a(oneLooseLepName+"_Theta1_SignalYield");
    TString names2a(oneLooseLepName+"_Theta2_SignalYield");
    TString names3a(oneLooseLepName+"_Theta3_SignalYield");
    TString names4a(oneLooseLepName+"_Theta4_SignalYield");
    TString names5a(oneLooseLepName+"_Theta5_SignalYield");

    tofill = (wspace2->function(names1a.Data()))->getVal() + (wspace2->function(names2a.Data()))->getVal()
      + (wspace2->function(names3a.Data()))->getVal() + (wspace2->function(names4a.Data()))->getVal() + (wspace2->function(names5a.Data()))->getVal();
    fillbin = (metbin*15 + 1) + (nb-1)*5 + 1;
    //if( metbin==3 ) fillbin = fillbin-5;
    signal_1L->Fill( fillbin, tofill );
    //if( metbin==0 ) truesig_1L->Fill( fillbin+5, 0.2*5. );
    //if( metbin==3 ) truesig_1L->Fill( fillbin+5, 2.6*5. );

    cout << "  loose lep sig " << endl; 
    //////
    TString names6(twoTightMuName+"_SignalYield");
    TString names7(twoLooseLepName+"_SignalYield");

    tofill = (wspace2->function(names6.Data()))->getVal();
    fillbin = (metbin*15 + 1) + (nb-1)*5 + 2;
    //if( metbin==3 ) fillbin = fillbin-5;
    signal_1L->Fill( fillbin, tofill );
    //if( metbin==0 ) truesig_1L->Fill( fillbin+5, 0. );
    //if( metbin==3 ) truesig_1L->Fill( fillbin+5, 0.4*5. );

    tofill = (wspace2->function(names7.Data()))->getVal();
    fillbin = (metbin*15 + 1) + (nb-1)*5 + 3;
    //if( metbin==3 ) fillbin = fillbin-5;
    signal_1L->Fill( fillbin, tofill );
    //if( metbin==0 ) truesig_1L->Fill( fillbin+5, 0. );
    //if( metbin==3 ) truesig_1L->Fill( fillbin+5, 0.6*5. );

    cout << "after SL signal " << endl;

    ///////////////////////////////////////////////////////////////////////////////////
    // SL background 

    TString nameb1(oneTightMuName+"_Theta1_TopWJetsYield");
    TString nameb2(oneTightMuName+"_Theta2_TopWJetsYield");
    TString nameb3(oneTightMuName+"_Theta3_TopWJetsYield");
    TString nameb4(oneTightMuName+"_Theta4_TopWJetsYield");
    TString nameb5(oneTightMuName+"_Theta5_TopWJetsYield");
    
    tofill = (wspace2->var(nameb1.Data()))->getVal() + (wspace2->var(nameb2.Data()))->getVal() 
      + (wspace2->var(nameb3.Data()))->getVal() + (wspace2->var(nameb4.Data()))->getVal() + (wspace2->var(nameb5.Data()))->getVal();
    fillbin = (metbin*15 + 1) + (nb-1)*5;
    //if( metbin==3 ) fillbin = fillbin-5;
    background_1L->Fill( fillbin, tofill );
    //if( metbin==0 ) truebkg_1L->Fill( fillbin+5, 11 );
    //if( metbin==3 ) truebkg_1L->Fill( fillbin+5, 0 );

    //////
    TString nameb1a(oneLooseLepName+"_Theta1_TopWJetsYield");
    TString nameb2a(oneLooseLepName+"_Theta2_TopWJetsYield");
    TString nameb3a(oneLooseLepName+"_Theta3_TopWJetsYield");
    TString nameb4a(oneLooseLepName+"_Theta4_TopWJetsYield");
    TString nameb5a(oneLooseLepName+"_Theta5_TopWJetsYield");

    tofill = (wspace2->var(nameb1a.Data()))->getVal() + (wspace2->var(nameb2a.Data()))->getVal()
      + (wspace2->var(nameb3a.Data()))->getVal() + (wspace2->var(nameb4a.Data()))->getVal() + (wspace2->var(nameb5a.Data()))->getVal();
    fillbin = (metbin*15 + 1) + (nb-1)*5 + 1;
    //if( metbin==3 ) fillbin = fillbin-5;
    background_1L->Fill( fillbin, tofill );
    //if( metbin==0 ) truebkg_1L->Fill( fillbin+5, 18 );
    //if( metbin==3 ) truebkg_1L->Fill( fillbin+5, 0 );

    //////
    TString nameb6(twoTightMuName+"_TopWJetsYield");
    TString nameb7(twoLooseLepName+"_TopWJetsYield");

    tofill = (wspace2->var(nameb6.Data()))->getVal();
    fillbin = (metbin*15 + 1) + (nb-1)*5 + 2;
    //if( metbin==3 ) fillbin = fillbin-5;
    background_1L->Fill( fillbin, tofill );
    //if( metbin==0 ) truebkg_1L->Fill( fillbin+5, 1 );
    //if( metbin==3 ) truebkg_1L->Fill( fillbin+5, 0 );

    tofill = (wspace2->var(nameb7.Data()))->getVal();
    fillbin = (metbin*15 + 1) + (nb-1)*5 + 3;
    //if( metbin==3 ) fillbin = fillbin-5;
    background_1L->Fill( fillbin, tofill );
    //if( metbin==0 ) truebkg_1L->Fill( fillbin+5, 1 );
    //if( metbin==3 ) truebkg_1L->Fill( fillbin+5, 0 );

    // end of SL observation

    cout << " after SL background " << endl;

    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    // 0L observation

    TString nameo01(zeroLeptonName+"_Count");
    tofill = (wspace2->var(nameo01.Data()))->getVal();
    fillbin = (metbin*5) + nb;
    //if( metbin==3 ) fillbin = fillbin-5;
    observed_0L->Fill( fillbin, tofill );

    // 0L signal

    TString names01(zeroLeptonName+"_SignalYield");
    tofill = (wspace2->function(names01.Data()))->getVal();
    fillbin = (metbin*5) + nb;
    //if( metbin==3 ) fillbin = fillbin-5;
    signal_0L->Fill( fillbin, tofill );
    //if( metbin==0 ) truesig_0L->Fill( fillbin+5, 3.8*5. );
    //if( metbin==3 ) truesig_0L->Fill( fillbin+5, 19.6*5. );


    cout << " after 0L signal and observation " << endl;
    // 0L background

    TString nameb01(zeroLeptonName+"_Theta1_TopWJetsYield");
    TString nameb02(zeroLeptonName+"_Theta2_TopWJetsYield");
    TString nameb03(zeroLeptonName+"_Theta3_TopWJetsYield");
    TString nameb04(zeroLeptonName+"_Theta4_TopWJetsYield");
    TString nameb05(zeroLeptonName+"_Theta5_TopWJetsYield");
    tofill = (wspace2->function(nameb01.Data()))->getVal() + (wspace2->function(nameb02.Data()))->getVal() 
      + (wspace2->function(nameb03.Data()))->getVal() + (wspace2->function(nameb04.Data()))->getVal() 
      + (wspace2->function(nameb05.Data()))->getVal() ;
    fillbin = (metbin*5) + nb;
    //if( metbin==3 ) fillbin = fillbin-5;
    bkg_pol_0L->Fill( fillbin, tofill );

    //if( metbin==0 ) truebkg_0L->Fill( fillbin+5, 67 );
    //if( metbin==3 ) truebkg_0L->Fill( fillbin+5, 3 );

    TString nameb11(zeroLeptonName+"_1Tau_TopWJetsYield");
    tofill = (wspace2->function(nameb11.Data()))->getVal();
    fillbin = (metbin*5) + nb;
    //if( metbin==3 ) fillbin = fillbin-5;
    bkg_tau_0L->Fill( fillbin, tofill );

    TString nameb21(zeroLeptonName+"_2Tau_TopWJetsYield");
    TString nameb22(zeroLeptonName+"_Dilep_TopWJetsYield");
    tofill = (wspace2->function(nameb21.Data()))->getVal() + (wspace2->function(nameb22.Data()))->getVal();
    fillbin = (metbin*5) + nb;
    //if( metbin==3 ) fillbin = fillbin-5;
    bkg_dil_0L->Fill( fillbin, tofill );


    cout << " after 0L background " << metbin << "  " << nb << endl;

    ///////////////////////////////////////////////////////////////////////////////////

    //     }// select a few bins

    }
  }

  TCanvas *c1 = new TCanvas();

  bkg_pol_0L->SetFillColor(kRed-5);
  bkg_tau_0L->SetFillColor(kBlue-5);
  bkg_dil_0L->SetFillColor(kGreen-5);
  signal_0L->SetFillColor(kMagenta);
  observed_0L->SetMarkerStyle(20);

  THStack *st1 = new THStack();
  st1->Add(bkg_pol_0L);
  st1->Add(bkg_tau_0L);
  st1->Add(bkg_dil_0L);
  st1->Add(signal_0L);
  observed_0L->Draw("elp");
  st1->Draw("histsame");
  observed_0L->Draw("elpsame");
/*
  truebkg_0L->SetFillColor(kBlue+1);
  truesig_0L->SetFillColor(kMagenta+1);
  THStack *st1a = new THStack();
  st1a->Add(truebkg_0L);
  st1a->Add(truesig_0L);
  st1a->Draw("histsame");
*/
  c1->SaveAs("BINS_HT3_xsec1_0L.eps");

  //////

  TCanvas *c2 = new TCanvas();

  background_1L->SetFillColor(kBlue-5);
  signal_1L->SetFillColor(kMagenta);
  observed_1L->SetMarkerStyle(20);

  THStack *st2 = new THStack();
  st2->Add(background_1L);
  st2->Add(signal_1L);
  observed_1L->Draw("elp");
  st2->Draw("histsame");
  observed_1L->Draw("elpsame");
/*
  truebkg_1L->SetFillColor(kBlue+1);
  truesig_1L->SetFillColor(kMagenta+1);
  THStack *st2a = new THStack();
  st2a->Add(truebkg_1L);
  st2a->Add(truesig_1L);
  st2a->Draw("histsame");
*/
  c2->SaveAs("BINS_HT3_xsec1_1L.eps");





  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	



}
