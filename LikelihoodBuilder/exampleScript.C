exampleScript()
{
  gSystem->CompileMacro("betaHelperFunctions.h"      ,"kO") ;
  gSystem->CompileMacro("RooNormalFromFlatPdf.cxx"      ,"kO") ;
  gSystem->CompileMacro("RooBetaInverseCDF.cxx"      ,"kO") ;
  gSystem->CompileMacro("RooBetaPrimeInverseCDF.cxx" ,"kO") ;
  gSystem->CompileMacro("RooCorrelatedBetaGeneratorHelper.cxx"  ,"kO") ;
  gSystem->CompileMacro("RooCorrelatedBetaPrimeGeneratorHelper.cxx"  ,"kO") ;
  gSystem->CompileMacro("rooFitBetaHelperFunctions.h","kO") ;

  TFile betaTest("betaTest.root","RECREATE");
  betaTest.cd();
  
  RooWorkspace workspace("workspace");
  TString correlatedName("testVariable");
  TString observables("observables");
  TString nuisances("nuisances");

  RooAbsArg* betaOne = getCorrelatedBetaConstraint(workspace,"betaOne","",
						   0.5 , 0.1 ,
						   observables, nuisances,
						   correlatedName );

  printf("\n\n *** constraint name is %s from betaOne and %s\n\n", betaOne->GetName(), correlatedName.Data() ) ;

  RooAbsArg* betaTwo = getCorrelatedBetaConstraint(workspace,"betaTwo","",
						   0 , 0 ,
						   observables, nuisances,
						   correlatedName );

  RooAbsArg* betaThree = getCorrelatedBetaConstraint(workspace,"betaThree","",
						     0.2 , 0.01 ,
						     observables, nuisances,
						     correlatedName );

  RooAbsArg* betaFour = getCorrelatedBetaConstraint(workspace,"betaFour","",
						    0.7 , 0.1 ,
						    observables, nuisances,
						    correlatedName );

  RooAbsArg* betaFourC = getCorrelatedBetaConstraint(workspace,"betaFourC","",
						    0.7 , 0.1 ,
						    observables, nuisances,
						    correlatedName, kTRUE );

  RooAbsArg* betaPrimeOne = getCorrelatedBetaPrimeConstraint(workspace,"betaPrimeOne","",
							     1.0 , 0.5 ,
							     observables, nuisances,
							     correlatedName );

  RooAbsArg* betaPrimeOneC = getCorrelatedBetaPrimeConstraint(workspace,"betaPrimeOneC","",
							     1.0 , 0.5 ,
							     observables, nuisances,
							     correlatedName, kTRUE );

  RooAbsArg* betaPrimeTwo = getCorrelatedBetaPrimeConstraint(workspace,"betaPrimeTwo","",
							     0.7 , 0.5 ,
							     observables, nuisances,
							     correlatedName );

  RooAbsArg* betaPrimeThree = getCorrelatedBetaPrimeConstraint(workspace,"betaPrimeThree","",
							       0.1 , 0.05 ,
							       observables, nuisances,
							       correlatedName );

  RooAbsArg* betaPrimeFour = getCorrelatedBetaPrimeConstraint(workspace,"betaPrimeFour","",
							      7 , 1 ,
							      observables, nuisances,
							      correlatedName );

  RooRealVar* correlatedParameter = workspace.var(correlatedName);

  RooAbsPdf* normalFromFlat = workspace.pdf(correlatedName+"_Constraint");

  RooDataSet* data = normalFromFlat->generate(RooArgSet(*correlatedParameter),1e5);

  data->addColumn(*normalFromFlat);

  data->addColumn(*betaOne);
  data->addColumn(*betaTwo);
  data->addColumn(*betaThree);
  data->addColumn(*betaFour);
  data->addColumn(*betaFourC);
  
  data->addColumn(*betaPrimeOne);
  data->addColumn(*betaPrimeTwo);
  data->addColumn(*betaPrimeThree);
  data->addColumn(*betaPrimeFour);
  data->addColumn(*betaPrimeOneC);

  data->Print("v");

  workspace.Print() ;

  //Setup Plotting Kluges:

  RooRealVar normalPlotter  (correlatedName+"_Constraint" , correlatedName+"_Constraint"  ,0,1);
  RooPlot* normalPlot = normalPlotter.frame();
  data->plotOn(normalPlot);

  RooRealVar betaOnePlotter  ("betaOne_BetaInverseCDF"  ,"betaOne_BetaInverseCDF"  ,0,1);
  RooRealVar betaTwoPlotter  ("betaTwo_BetaInverseCDF"  ,"betaTwo_BetaInverseCDF"  ,0,1);
  RooRealVar betaThreePlotter("betaThree_BetaInverseCDF","betaThree_BetaInverseCDF",0,1);
  RooRealVar betaFourPlotter ("betaFour_BetaInverseCDF" ,"betaFour_BetaInverseCDF" ,0,1);
  RooRealVar betaFourCPlotter ("betaFourC_BetaInverseCDF" ,"betaFourC_BetaInverseCDF" ,0,1);

  RooRealVar betaPrimeOnePlotter  ("betaPrimeOne_BetaPrimeInverseCDF"  ,"betaPrimeOne_BetaPrimeInverseCDF"  ,0,4);
  RooRealVar betaPrimeOneCPlotter  ("betaPrimeOneC_BetaPrimeInverseCDF"  ,"betaPrimeOneC_BetaPrimeInverseCDF"  ,0,4);
  RooRealVar betaPrimeTwoPlotter  ("betaPrimeTwo_BetaPrimeInverseCDF"  ,"betaPrimeTwo_BetaPrimeInverseCDF"  ,0,4);
  RooRealVar betaPrimeThreePlotter("betaPrimeThree_BetaPrimeInverseCDF","betaPrimeThree_BetaPrimeInverseCDF",0,0.3);
  RooRealVar betaPrimeFourPlotter ("betaPrimeFour_BetaPrimeInverseCDF" ,"betaPrimeFour_BetaPrimeInverseCDF" ,4,12);

  RooPlot* betaOnePlot   = betaOnePlotter  .frame();
  RooPlot* betaTwoPlot   = betaTwoPlotter  .frame();
  RooPlot* betaThreePlot = betaThreePlotter.frame();
  RooPlot* betaFourPlot  = betaFourPlotter .frame();
  RooPlot* betaFourCPlot  = betaFourCPlotter .frame();

  data->plotOn(betaOnePlot  );
  data->plotOn(betaTwoPlot  );
  data->plotOn(betaThreePlot);
  data->plotOn(betaFourPlot );
  data->plotOn(betaFourCPlot );

  RooPlot* betaPrimeOnePlot   = betaPrimeOnePlotter  .frame();
  RooPlot* betaPrimeOneCPlot   = betaPrimeOneCPlotter  .frame();
  RooPlot* betaPrimeTwoPlot   = betaPrimeTwoPlotter  .frame();
  RooPlot* betaPrimeThreePlot = betaPrimeThreePlotter.frame();
  RooPlot* betaPrimeFourPlot  = betaPrimeFourPlotter .frame();

  data->plotOn(betaPrimeOnePlot  );
  data->plotOn(betaPrimeOneCPlot  );
  data->plotOn(betaPrimeTwoPlot  );
  data->plotOn(betaPrimeThreePlot);
  data->plotOn(betaPrimeFourPlot );

  TCanvas* underlyingVariable = new TCanvas("underlyingVariable","underlyingVariable",800,800);
  underlyingVariable->Divide(2,2);
  underlyingVariable->cd(1);
  RooPlot* underlyingPlot   = correlatedParameter->frame();
  data->plotOn(underlyingPlot);
  underlyingPlot->Draw();
  underlyingVariable->cd(2);
  normalPlot->Draw();
  underlyingVariable->cd(3);
  TH2F* underlying = data->createHistogram(*correlatedParameter,normalPlotter,50,50);
  underlying->Draw("col");
  TH2F* legoUnderlying = (TH2F*)underlying->Clone();
  underlyingVariable->cd(4);
  legoUnderlying->Draw("lego");

  underlyingVariable->SaveAs("underlyingVariable.pdf");
  
  TCanvas* betaCanvas = new TCanvas("betaCanvas","betaCanvas",800,800);
  
  betaCanvas->Divide(3,2);
  
  betaCanvas->cd(1);
  betaOnePlot->Draw();
  betaCanvas->cd(2);
  betaTwoPlot->Draw();
  betaCanvas->cd(3);
  betaThreePlot->Draw();
  betaCanvas->cd(4);
  betaFourPlot->Draw();
  betaCanvas->cd(5);
  betaFourCPlot->Draw();

  betaCanvas->SaveAs("betaVariables.pdf");

  TCanvas* betaPrimeCanvas = new TCanvas("betaPrimeCanvas","betaPrimeCanvas",1200,800);
  
  betaPrimeCanvas->Divide(3,2);
  
  betaPrimeCanvas->cd(1);
  betaPrimeOnePlot->Draw();
  betaPrimeCanvas->cd(2);
  betaPrimeTwoPlot->Draw();
  betaPrimeCanvas->cd(3);
  betaPrimeThreePlot->Draw();
  betaPrimeCanvas->cd(4);
  betaPrimeFourPlot->Draw();
  betaPrimeCanvas->cd(5);
  betaPrimeOneCPlot->Draw();

  betaPrimeCanvas->SaveAs("betaPrimeVariables.pdf");
  
  TCanvas* betaCorrelationsCanvas = new TCanvas("betaCorrelationsCanvas","betaCorrelationsCanvas",1600,800);
  
  betaCorrelationsCanvas->Divide(4,2);

  TH2F* oneTwo = data->createHistogram(betaOnePlotter,betaTwoPlotter,30,30);
  TH2F* oneThree = data->createHistogram(betaOnePlotter,betaThreePlotter,30,30);
  TH2F* oneFour = data->createHistogram(betaOnePlotter,betaFourPlotter,30,30);
  TH2F* twoThree = data->createHistogram(betaTwoPlotter,betaThreePlotter,30,30);
  TH2F* twoFour = data->createHistogram(betaTwoPlotter,betaFourPlotter,30,30);
  TH2F* threeFour = data->createHistogram(betaThreePlotter,betaFourPlotter,30,30);
  TH2F* twoFourC = data->createHistogram(betaTwoPlotter,betaFourCPlotter,30,30);
  TH2F* fourFourC = data->createHistogram(betaFourPlotter,betaFourCPlotter,30,30);

  betaCorrelationsCanvas->cd(1);
  oneTwo->DrawCopy("lego");
  betaCorrelationsCanvas->cd(2);
  oneThree->DrawCopy("lego");
  betaCorrelationsCanvas->cd(3);
  oneFour->DrawCopy("lego");
  betaCorrelationsCanvas->cd(4);
  twoThree->DrawCopy("lego");
  betaCorrelationsCanvas->cd(5);
  twoFour->DrawCopy("lego");
  betaCorrelationsCanvas->cd(6);
  threeFour->DrawCopy("lego");
  betaCorrelationsCanvas->cd(7);
  twoFourC->DrawCopy("lego");
  betaCorrelationsCanvas->cd(8);
  fourFourC->DrawCopy("lego");

  betaCorrelationsCanvas->SaveAs("betaCorrelations.pdf");

  TCanvas* betaPrimeCorrelationsCanvas = new TCanvas("betaPrimeCorrelationsCanvas","betaPrimeCorrelationsCanvas",1600,800);
  
  betaPrimeCorrelationsCanvas->Divide(4,2);

  TH2F* oneTwo = data->createHistogram(betaPrimeOnePlotter,betaPrimeTwoPlotter,30,30);
  TH2F* oneThree = data->createHistogram(betaPrimeOnePlotter,betaPrimeThreePlotter,30,30);
  TH2F* oneFour = data->createHistogram(betaPrimeOnePlotter,betaPrimeFourPlotter,30,30);
  TH2F* twoThree = data->createHistogram(betaPrimeTwoPlotter,betaPrimeThreePlotter,30,30);
  TH2F* twoFour = data->createHistogram(betaPrimeTwoPlotter,betaPrimeFourPlotter,30,30);
  TH2F* threeFour = data->createHistogram(betaPrimeThreePlotter,betaPrimeFourPlotter,30,30);
  TH2F* oneOneC = data->createHistogram(betaPrimeOnePlotter,betaPrimeOneCPlotter,30,30);

  betaPrimeCorrelationsCanvas->cd(1);
  oneTwo->DrawCopy("lego");
  betaPrimeCorrelationsCanvas->cd(2);
  oneThree->DrawCopy("lego");
  betaPrimeCorrelationsCanvas->cd(3);
  oneFour->DrawCopy("lego");
  betaPrimeCorrelationsCanvas->cd(4);
  twoThree->DrawCopy("lego");
  betaPrimeCorrelationsCanvas->cd(5);
  twoFour->DrawCopy("lego");
  betaPrimeCorrelationsCanvas->cd(6);
  threeFour->DrawCopy("lego");
  betaPrimeCorrelationsCanvas->cd(7);
  oneOneC->DrawCopy("lego");

  betaPrimeCorrelationsCanvas->SaveAs("betaPrimeCorrelations.pdf");

  RooProdPdf totalPdf("totalPdf","totalPdf",workspace.allPdfs());
  totalPdf.Print("v");

  RooArgSet* observableSet = workspace.set("observables");

  observableSet->Print();

  RooDataSet* allDataOne = totalPdf.generate(*observableSet,1);
  allDataOne->Print("v");

  correlatedParameter->setVal(0.25);

  RooDataSet* allDataTwo = totalPdf.generate(*observableSet,1);
  allDataTwo->Print("v");

  correlatedParameter->setVal(0.75);

  RooDataSet* allDataThree = totalPdf.generate(*observableSet,1);
  allDataThree->Print("v");

  //Testing for extreme values!

  for(int i = 0; i< 101; i++)
    {
      correlatedParameter->setVal((double)i/100.);
      cout << "Correlation parameter has value of " << correlatedParameter->getVal();
      cout << " and the pdf has an unnormalized value of " << normalFromFlat->getVal() << endl;
    }


}
