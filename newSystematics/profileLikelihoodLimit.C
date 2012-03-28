//Limit Calculator.  
//First calculates limit using profile likelihood calculator.  
//Uses this result to set the scan region for the Frequentist Calculator
//Outputs a result in a text file with the format:
//m0 m12 asymptoticLimit toyLimit toyLimit_err expectedToyLimit expectedToyLimit_errHigh expectedToyLimit_errLow 

#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooProfileLL.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooCurve.h"
#include "RooHist.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooStats/ProfileLikelihoodTestStat.h"
#include "ProfileLikelihoodTestStat_New.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "DebuggingToyMCSampler.h"
#include "HybridToyMCSampler.h"

#include <iostream>
#include <string>
#include <fstream>

using namespace RooFit;
using namespace RooStats;

vector<TObject*> makeStatDistPlot( HypoTestResult* theResult )
{
  vector<TObject*> theObjs;
  if(!theResult) return theObjs;
  TH1F* nullHypBelow(0),*nullHypAbove(0),*altHypBelow(0),*altHypAbove(0);

  const SamplingDistribution *null = theResult->GetNullDistribution();
  const SamplingDistribution *alt = theResult->GetAltDistribution();

  vector<Double_t> nullDist = null->GetSamplingDistribution();
  vector<Double_t> altDist = alt->GetSamplingDistribution();

  vector<Double_t> nullWeightDist = null->GetSampleWeights();
  vector<Double_t> altWeightDist = alt->GetSampleWeights();

  cout << "Array sizes:" << endl;
  cout << nullDist.size() << endl;
  cout << altDist.size() << endl;

  //if( nullDist.size() == 0 || altDist.size() == 0 ) return theObjs ;

  Double_t nullmin(0) ;
  if( nullDist.size() != 0 ) nullmin = *(min_element(nullDist.begin(), nullDist.end())) ;
  Double_t nullmax(0) ;
  if( nullDist.size() != 0 ) nullmax = *(max_element(nullDist.begin(), nullDist.end())) ;
  
  cout << "The nullmax is " << nullmax << endl;
  cout << "The nullmin is " << nullmin << endl;

  if(nullmin != nullmin || nullmax!= nullmax || nullmax == numeric_limits<Double_t>::infinity())
    {
      nullmin = 1000.;
      nullmax = 0.;
      for(vector<Double_t>::iterator thisPoint = nullDist.begin() ; thisPoint != nullDist.end() ; thisPoint++ )
	{
	  if(*thisPoint != *thisPoint) continue ;
	  if(*thisPoint == numeric_limits<Double_t>::infinity() ) continue;
	  if(*thisPoint == -numeric_limits<Double_t>::infinity() ) continue;
	  if(*thisPoint>nullmax) nullmax = *thisPoint ;
	  
	  if(*thisPoint<nullmin) nullmin = *thisPoint ;
	}
    }

  Double_t altmin(0) ;
  if( altDist.size() != 0 ) altmin = *(min_element(altDist.begin(), altDist.end())) ;
  Double_t altmax(0) ;
  if( altDist.size() != 0 ) altmax = *(max_element(altDist.begin(), altDist.end())) ;

  cout << "The altmax is " << altmax << endl;
  cout << "The altmin is " << altmin << endl;

  if(altmin != altmin || altmax != altmax || altmax == numeric_limits<Double_t>::infinity() )
    {
      altmin = 1000.;
      altmax = 0.;
      for(vector<Double_t>::iterator thisPoint = altDist.begin() ; thisPoint != altDist.end() ; thisPoint++ )
	{
	  if(*thisPoint != *thisPoint) continue ;
	  if(*thisPoint == numeric_limits<Double_t>::infinity() ) continue;
	  if(*thisPoint == -numeric_limits<Double_t>::infinity() ) continue;
	  if(*thisPoint>altmax) altmax = *thisPoint ;
	  if(*thisPoint<altmin) altmin = *thisPoint ;
	}
    }

  double xmin = min( nullmin , altmin );
  double xmax = max( nullmax , altmax );
  double xmid = theResult->HasTestStatisticData() ? theResult->GetTestStatisticData() : -1. ;
  if( xmid != xmid || xmid == numeric_limits<Double_t>::infinity() ) xmid = -1 ;
  

  if(xmin == xmax)
    {
      xmin = xmax - 1.0;
      xmax = xmax + 1.0;
    }

  nullHypBelow = new TH1F("nullHypBelow","nullHypBelow",100,xmin,xmax);
  nullHypAbove = new TH1F("nullHypAbove","nullHypAbove",100,xmin,xmax);
  altHypBelow =  new TH1F("altHypBelow", "altHypBelow" ,100,xmin,xmax);
  altHypAbove =  new TH1F("altHypAbove", "altHypAbove" ,100,xmin,xmax);

  theObjs.push_back((TObject*)nullHypBelow);
  theObjs.push_back((TObject*)nullHypAbove);
  theObjs.push_back((TObject*)altHypBelow);
  theObjs.push_back((TObject*)altHypAbove);

  nullHypAbove->SetFillStyle(3005);
  altHypAbove->SetFillStyle(3004);
  nullHypAbove->SetFillColor(kRed);
  altHypAbove->SetFillColor(kBlue);
  
  nullHypBelow->SetLineColor(kRed);
  nullHypAbove->SetLineColor(kRed);
  altHypBelow->SetLineColor(kBlue);
  altHypAbove->SetLineColor(kBlue);

  nullHypBelow->SetLineWidth(1);
  nullHypAbove->SetLineWidth(1);
  altHypBelow->SetLineWidth(1);
  altHypAbove->SetLineWidth(1);

  for(vector<Double_t>::iterator thisPoint = nullDist.begin() ; thisPoint != nullDist.end() ; thisPoint++ )
    {
      if(*thisPoint != *thisPoint) continue ;
      if(*thisPoint == numeric_limits<Double_t>::infinity() ) continue;
      if(*thisPoint == -numeric_limits<Double_t>::infinity() ) continue;
      if(*thisPoint < xmid) nullHypBelow->Fill(*thisPoint);
      else nullHypAbove->Fill(*thisPoint);
    }

  for(vector<Double_t>::iterator thisPoint = altDist.begin() ; thisPoint != altDist.end() ; thisPoint++ )
    {
      if(*thisPoint != *thisPoint) continue ;
      if(*thisPoint == numeric_limits<Double_t>::infinity() ) continue;
      if(*thisPoint == -numeric_limits<Double_t>::infinity() ) continue;
      if(*thisPoint < xmid) altHypBelow->Fill(*thisPoint);
      else altHypAbove->Fill(*thisPoint);
    }

  if( nullHypBelow->Integral() + nullHypAbove->Integral() !=  altHypBelow->Integral() + altHypAbove->Integral() )
    {
      altHypBelow->Scale( ( nullHypBelow->Integral() + nullHypAbove->Integral() ) / ( altHypBelow->Integral() + altHypAbove->Integral() ) );
      altHypAbove->Scale( ( nullHypBelow->Integral() + nullHypAbove->Integral() ) / ( altHypBelow->Integral() + altHypAbove->Integral() ) );
    }

  double ymax = pow( max( max(nullHypBelow->GetMaximum(),nullHypAbove->GetMaximum()) , max(altHypBelow->GetMaximum(),altHypAbove->GetMaximum()) ) , 1.1 );

  TGraph* theTestStatistic = new TGraph(2);
  theTestStatistic->SetPoint(0,xmid,0);
  theTestStatistic->SetPoint(1,xmid,ymax);
  theTestStatistic->SetLineColor(kBlack);
  theTestStatistic->SetLineWidth(1);  

  theObjs.push_back((TObject*)theTestStatistic);

  altHypBelow->Draw("HIST");
  altHypAbove->Draw("HIST SAME");
  nullHypBelow->Draw("HIST SAME");
  nullHypAbove->Draw("HIST SAME");
  theTestStatistic->Draw("SAME");

  TLegend* theLegend = new TLegend(0.70,0.95-0.2*0.66,0.95,0.95);
  theLegend->AddEntry(altHypBelow,"B_model","l");
  theLegend->AddEntry(nullHypBelow,"S+B_model","l");
  theLegend->AddEntry(theTestStatistic,"test statistic data","l");
  theLegend->Draw();

  theObjs.push_back((TObject*)theLegend);

  cout << "The xmax is " << xmax << endl;
  cout << "The xmin is " << xmin << endl;
  cout << "The ymax is " << ymax << endl;
  cout << "The xmid is " << xmid << endl;


  return theObjs;
}

void profileLikelihoodLimit(const char * fileName =0,
			    const char * wsName = "combined",
			    const char * modelSBName = "ModelConfig",
			    const char * modelBName = "",
			    const char * dataName = "obsData",    
			    const char * modelName = "",
			    bool isMeasured = true,
			    int m0 = 0,
			    int m12 = 0,
			    int nToys = 0,
			    int nToysMax = 0,
			    int nScanPoints = 10,
			    bool isFrequentist = true,
			    bool isHybrid = false,
			    bool doPlots = false,
			    bool makePlots = false,
			    bool makeMorePlots = false)
{

  TFile *file = TFile::Open(fileName);
  // if input file was specified byt not found, quit
  if(!file){
    cout <<"file not found" << endl;
    return;
  } 

  // get the workspace out of the file
  RooWorkspace* w = (RooWorkspace*) file->Get(wsName);
  if(!w){
    cout <<"workspace not found" << endl;
    return;
  }

  // get the modelConfig out of the file
  ModelConfig* mc = (ModelConfig*) w->obj(modelSBName);

  // get the modelConfig out of the file
  RooAbsData* data = w->data(dataName);
  RooAbsPdf* pdf = mc->GetPdf();

  // make sure ingredients are found
  if(!data || !mc){
    w->Print();
    cout << "data or ModelConfig was not found" <<endl;
    return;
  }

  double tempMax = ((RooRealVar*)data->get()->find("Nsig"))->getVal() > 0 ? ((RooRealVar*)data->get()->find("Nsig"))->getVal() + 3*sqrt(((RooRealVar*)data->get()->find("Nsig"))->getVal()) : 5 ;

  RooRealVar* firstPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
  firstPOI->setRange(0, tempMax);

  // create and use the ProfileLikelihoodCalculator
  // to find and plot the 95% confidence interval
  // on the parameter of interest as specified
  // in the model config

  //ProfileLikelihoodCalculator pl(*data,*mc);
  ProfileLikelihoodCalculator pl(*data,*pdf, RooArgSet(*firstPOI));
  pl.SetTestSize( 0.025 ) ;
  //pl.SetConfidenceLevel(0.975); // 95% one sided limit
  LikelihoodInterval* interval = pl.GetInterval();
  //LikelihoodInterval* interval = new LikelihoodInterval("theInterval",profile,&bestFit,bestFitSnapshot);
  ///interval->SetConfidenceLevel(0.975);

  //// print out the iterval on the first Parameter of Interest

  cout << "\n95% interval on " <<firstPOI->GetName()<<" is : ["<<
    interval->LowerLimit(*firstPOI) << ", "<<
    interval->UpperLimit(*firstPOI) <<"] "<<endl;

  double asymptoticLimit = interval->UpperLimit(*firstPOI);
  //double asymptoticLimit = profile->getVal();

  TCanvas* c ;

  double mediumTempMax = ((RooRealVar*)data->get()->find("Nsig"))->getVal() > 0 ? ((RooRealVar*)data->get()->find("Nsig"))->getVal() + 5*sqrt(((RooRealVar*)data->get()->find("Nsig"))->getVal()) : 7 ; 

  if(doPlots){
    firstPOI->setRange(0, mediumTempMax );
    c = new TCanvas("c");
    cout << "making a plot of the profile likelihood function ....(if it is taking a lot of time use less points or the TF1 drawing option)\n";
    LikelihoodIntervalPlot plot(interval);
    plot.SetNPoints(50);  // do not use too many points, it could become very slow for some models
    plot.SetMaximum(4); 
    plot.SetRange(0.,1.3*asymptoticLimit);
    plot.Draw("");  // use option TF1 if too slow (plot.Draw("tf1")
    c->SaveAs(TString(modelName)+"_profileLikelihoodScan.pdf");
    delete c;
  }

  double bigTempMax = ((RooRealVar*)data->get()->find("Nsig"))->getVal() > 0 ? ((RooRealVar*)data->get()->find("Nsig"))->getVal() + 10*sqrt(((RooRealVar*)data->get()->find("Nsig"))->getVal()) : 10 ; 

  firstPOI->setRange(0, bigTempMax );

  //Now set up the FrequentistCalculator so that it can look around the asymptotic limit

  TFile *fDebug;

  if(doPlots)
    {
      fDebug = new TFile("debuggingFile.root","RECREATE");
      fDebug->cd();
    }

  ModelConfig* bModel = (ModelConfig*) w->obj(modelBName);
  ModelConfig* sbModel = (ModelConfig*) w->obj(modelSBName);

  RooArgSet dataSetVars(*sbModel->GetPdf()->getVariables());
  RemoveConstantParameters(&dataSetVars);
  RooRealVar poiHyp(*firstPOI,(TString(firstPOI->GetName())+"_hypothesis"));
  dataSetVars.add(poiHyp);
  RooRealVar poiFit(*firstPOI,(TString(firstPOI->GetName())+"_fitted"));
  dataSetVars.add(poiFit);
  dataSetVars.Print("v");

  RooDataSet genDebugging("genDebugging","genDebugging",dataSetVars);
  RooDataSet trueDebugging("trueDebugging","trueDebugging",dataSetVars);
  RooDataSet uncondDebugging("uncondDebugging","uncondDebugging",dataSetVars);
  RooDataSet condDebugging("condDebugging","condDebugging",dataSetVars);

  ProfileLikelihoodTestStat_New profll(*sbModel->GetPdf());
  //profll.SetDoProfllCheck(false);
  profll.SetOneSided(1);
  profll.SetReuseNLL(true);

  DebuggingToyMCSampler poiSampler(profll,nToys);  

  cout << "limit scan on data" << endl;

  FrequentistCalculator frequentist_calculator(*data, *bModel, *sbModel,&poiSampler);
  ((ToyMCSampler*)frequentist_calculator.GetTestStatSampler())->SetNEventsPerToy(1);
  ((ToyMCSampler*)frequentist_calculator.GetTestStatSampler())->SetTestStatistic(&profll);
  ((ToyMCSampler*)frequentist_calculator.GetTestStatSampler())->SetUseMultiGen(true);
  ((ToyMCSampler*)frequentist_calculator.GetTestStatSampler())->SetMaxToys(nToysMax);

  ((DebuggingToyMCSampler*)frequentist_calculator.GetTestStatSampler())->SetFitToData(data);
  if(doPlots) ((DebuggingToyMCSampler*)frequentist_calculator.GetTestStatSampler())->SetDebuggingData(&genDebugging);
  if(doPlots) profll.SetDebuggingData(&uncondDebugging , &condDebugging, &genDebugging , &trueDebugging );
  //profll.SetFluctuatedNuisanceParameters(*sbModel->GetNuisanceParameters());

  //((ToyMCSampler_New*)frequentist_calculator.GetTestStatSampler())->SetFitToData(theData);
  frequentist_calculator.SetToys(nToys,nToys);

  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);

  HypoTestInverter calc_frequentist_calculator(frequentist_calculator);

  calc_frequentist_calculator.SetConfidenceLevel(0.95);
  calc_frequentist_calculator.SetVerbose(true);
  calc_frequentist_calculator.UseCLs(true);
  double poimin = 0.5*asymptoticLimit;
  double poimax = 2.0*asymptoticLimit;
  if(!(asymptoticLimit > 0.&& asymptoticLimit < tempMax) )
    {
      poimin = 0;
      poimax = mediumTempMax;
    }

  cout << "The min is " << poimin << endl;
  cout << "The max is " << poimax << endl;
  cout << "There are " << nScanPoints << " scan points." << endl;

  calc_frequentist_calculator.SetFixedScan( nScanPoints , poimin , poimax) ;

  cout << "setup all done" << endl;

  HypoTestInverterResult* res_frequentist_calculator;
  if(isFrequentist) res_frequentist_calculator = calc_frequentist_calculator.GetInterval();

  cout << "done with the interval calculation" << endl; 
  
  if( isFrequentist )for(int j = 0 ; j < res_frequentist_calculator->ArraySize() ; j++)
    {
        cout << "CLb(" << j << ")          : " <<res_frequentist_calculator->CLb(j)           << endl;
  	cout << "CLbError(" << j << ")     : " <<res_frequentist_calculator->CLbError(j)      << endl;
  	cout << "CLs(" << j << ")          : " <<res_frequentist_calculator->CLs(j) 	      << endl;
  	cout << "CLsError(" << j << ")     : " <<res_frequentist_calculator->CLsError(j)      << endl;
  	cout << "CLsplusb(" << j << ")     : " <<res_frequentist_calculator->CLsplusb(j)      << endl;
  	cout << "CLsplusbError(" << j << "): " <<res_frequentist_calculator->CLsplusbError(j) << endl;
    }
  
  cout << "getting relevant data information:" << endl;
  
  if(isFrequentist) res_frequentist_calculator->SetInterpolationOption(HypoTestInverterResult::kLinear);
  
  double upperLimit_linear             = !isFrequentist ? 0 : res_frequentist_calculator->UpperLimit();
  double upperLimitError_linear        = !isFrequentist ? 0 : res_frequentist_calculator->UpperLimitEstimatedError();
  double expectedUpperLimit_linear     = !isFrequentist ? 0 : res_frequentist_calculator->GetExpectedUpperLimit(0);
  double expectedUpperLimitHigh_linear = !isFrequentist ? 0 : res_frequentist_calculator->GetExpectedUpperLimit(1);
  double expectedUpperLimitLow_linear  = !isFrequentist ? 0 : res_frequentist_calculator->GetExpectedUpperLimit(-1);

  HypoTestInverterResult* thisResult;
  vector<vector<TObject*> > objPointerVectors;

  if(makePlots && isFrequentist)
    {
      HypoTestInverterPlot plot_frequentist_calculator( "profll" , "profll" , res_frequentist_calculator ) ;
      
      c = new TCanvas("plot") ;
      c->cd();
      c->SetLogy(0);
      plot_frequentist_calculator.Draw("CLB 2CL");
      TString limitPlotName("frequentist_limit_");
      limitPlotName+=modelName;
      limitPlotName+=".pdf";
      c->SaveAs(limitPlotName) ;
      delete c;
      if(makeMorePlots)
  	{
  	  thisResult = res_frequentist_calculator;
  	  c = new TCanvas("plot" , "plot" , 2400, 600*TMath::Ceil(thisResult->ArraySize()/2.));
  	  c->SetLogy(1);
  	  c->Divide( 2, TMath::Ceil(thisResult->ArraySize()/2.));
  	  for (int i=0; i<thisResult->ArraySize(); i++) {
  	    cout << "Drawing plot number " << i << endl;
  	    c->cd(i+1)->SetLogy(1);
  	    objPointerVectors.push_back(makeStatDistPlot(thisResult->GetResult(i)));
  	  }
  	  TString distPlotName("frequentist_dist_");
  	  distPlotName+=modelName;
  	  distPlotName+=".pdf";
  	  c->SaveAs(distPlotName);
  	  delete c;
  	}
    }

  //for(vector<vector<TObject*> >::iterator thisVector = objPointerVectors.begin(); thisVector != objPointerVectors.end() ; thisVector++ )
  //  {
  //    for(vector<TObject*>::iterator thisObj = thisVector->begin(); thisObj != thisVector->end(); thisObj++)
  //	{
  //	  cout << *thisObj << endl;
  //	  if( *thisObj ) delete *thisObj;
  //	}
  //  }

  cout << "getting CLs+b" << endl;

  if(isFrequentist) res_frequentist_calculator->UseCLs(false);

  double clsplusb_upperLimit_linear             =  !isFrequentist ? 0 : res_frequentist_calculator->UpperLimit();
  double clsplusb_upperLimitError_linear        =  !isFrequentist ? 0 : res_frequentist_calculator->UpperLimitEstimatedError();
  double clsplusb_expectedUpperLimit_linear     =  !isFrequentist ? 0 : res_frequentist_calculator->GetExpectedUpperLimit(0);
  double clsplusb_expectedUpperLimitHigh_linear =  !isFrequentist ? 0 : res_frequentist_calculator->GetExpectedUpperLimit(1);
  double clsplusb_expectedUpperLimitLow_linear  =  !isFrequentist ? 0 : res_frequentist_calculator->GetExpectedUpperLimit(-1);

  ///// --- BEGIN HYBRID STUFF

  RooDataSet genHybridDebugging("genHybridDebugging","genHybridDebugging",dataSetVars);
  RooDataSet trueHybridDebugging("trueHybridDebugging","trueHybridDebugging",dataSetVars);
  RooDataSet uncondHybridDebugging("uncondHybridDebugging","uncondHybridDebugging",dataSetVars);
  RooDataSet condHybridDebugging("condHybridDebugging","condHybridDebugging",dataSetVars);

  ProfileLikelihoodTestStat_New hybridProfll(*sbModel->GetPdf());
  //hybridProfll.SetDoHybridProfllCheck(false);
  hybridProfll.SetOneSided(1);
  hybridProfll.SetReuseNLL(true);

  HybridToyMCSampler hybridPoiSampler(hybridProfll,nToys);  

  cout << "limit scan on data" << endl;

  RooAbsPdf* nuisancePrior = w->pdf("nuisancePrior");
  RooArgSet nonPoissonNuisances(  *w->set("allNonPoissonNuisances") );

  HybridCalculator hybrid_calculator(*data, *bModel, *sbModel,&hybridPoiSampler);
  ((ToyMCSampler*)hybrid_calculator.GetTestStatSampler())->SetNEventsPerToy(1);
  ((ToyMCSampler*)hybrid_calculator.GetTestStatSampler())->SetTestStatistic(&hybridProfll);
  ((ToyMCSampler*)hybrid_calculator.GetTestStatSampler())->SetUseMultiGen(true);
  ((ToyMCSampler*)hybrid_calculator.GetTestStatSampler())->SetMaxToys(nToysMax);
  ((HybridToyMCSampler*)hybrid_calculator.GetTestStatSampler())->SetSampledNuisances(&nonPoissonNuisances);

  hybrid_calculator.ForcePriorNuisanceAlt(*nuisancePrior);
  hybrid_calculator.ForcePriorNuisanceNull(*nuisancePrior);

  ((HybridToyMCSampler*)hybrid_calculator.GetTestStatSampler())->SetFitToData(data);
  if(doPlots) ((HybridToyMCSampler*)hybrid_calculator.GetTestStatSampler())->SetDebuggingData(&genHybridDebugging);
  if(doPlots) hybridProfll.SetDebuggingData(&uncondHybridDebugging , &condHybridDebugging, &genHybridDebugging , &trueHybridDebugging );

  hybrid_calculator.SetToys(nToys,nToys);

  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);

  HypoTestInverter calc_hybrid_calculator(hybrid_calculator);

  calc_hybrid_calculator.SetConfidenceLevel(0.95);
  calc_hybrid_calculator.SetVerbose(true);
  calc_hybrid_calculator.UseCLs(true);
  calc_hybrid_calculator.SetFixedScan( nScanPoints , poimin , poimax) ;

  cout << "hybrid setup all done" << endl;

  HypoTestInverterResult* res_hybrid_calculator ;
  if (isHybrid) res_hybrid_calculator = calc_hybrid_calculator.GetInterval();

  cout << "done with the hybrid interval calculation" << endl; 

  if( isHybrid ) for(int j = 0 ; j < res_hybrid_calculator->ArraySize() ; j++)
    {
        cout << "CLb(" << j << ")          : " <<res_hybrid_calculator->CLb(j)           << endl;
	cout << "CLbError(" << j << ")     : " <<res_hybrid_calculator->CLbError(j)      << endl;
	cout << "CLs(" << j << ")          : " <<res_hybrid_calculator->CLs(j) 	      << endl;
	cout << "CLsError(" << j << ")     : " <<res_hybrid_calculator->CLsError(j)      << endl;
	cout << "CLsplusb(" << j << ")     : " <<res_hybrid_calculator->CLsplusb(j)      << endl;
	cout << "CLsplusbError(" << j << "): " <<res_hybrid_calculator->CLsplusbError(j) << endl;
    }
  
  cout << "getting relevant data information:" << endl;

  if(isHybrid) res_hybrid_calculator->SetInterpolationOption(HypoTestInverterResult::kLinear);
  
  double hybrid_upperLimit_linear             = !isHybrid ? 0 : res_hybrid_calculator->UpperLimit();
  double hybrid_upperLimitError_linear        = !isHybrid ? 0 : res_hybrid_calculator->UpperLimitEstimatedError();
  double hybrid_expectedUpperLimit_linear     = !isHybrid ? 0 : res_hybrid_calculator->GetExpectedUpperLimit(0);
  double hybrid_expectedUpperLimitHigh_linear = !isHybrid ? 0 : res_hybrid_calculator->GetExpectedUpperLimit(1);
  double hybrid_expectedUpperLimitLow_linear  = !isHybrid ? 0 : res_hybrid_calculator->GetExpectedUpperLimit(-1);

  HypoTestInverterResult* thisHybridResult;

  if(makePlots && isHybrid)
    {
      HypoTestInverterPlot plot_hybrid_calculator( "hybridProfll" , "hybridProfll" , res_hybrid_calculator ) ;
      
      c = new TCanvas("plot") ;
      c->cd();
      c->SetLogy(0);
      plot_hybrid_calculator.Draw("CLB 2CL");
      TString limitPlotName("hybrid_limit_");
      limitPlotName+=modelName;
      limitPlotName+=".pdf";
      c->SaveAs(limitPlotName) ;
      delete c;
      if(makeMorePlots)
	{
	  thisHybridResult = res_hybrid_calculator;
	  c = new TCanvas("plot" , "plot" , 2400, 600*TMath::Ceil(thisHybridResult->ArraySize()/2.));
	  c->SetLogy(1);
	  c->Divide( 2, TMath::Ceil(thisHybridResult->ArraySize()/2.));
	  for (int i=0; i<thisHybridResult->ArraySize(); i++) {
	    cout << "Drawing plot number " << i << endl;
	    c->cd(i+1)->SetLogy(1);
	    objPointerVectors.push_back(makeStatDistPlot(thisHybridResult->GetResult(i)));
	  }
	  TString distPlotName("hybrid_dist_");
	  distPlotName+=modelName;
	  distPlotName+=".pdf";
	  c->SaveAs(distPlotName);
	  delete c;
	}
    }
  for(vector<vector<TObject*> >::iterator thisVector = objPointerVectors.begin(); thisVector != objPointerVectors.end() ; thisVector++ )
    {
      for(vector<TObject*>::iterator thisObj = thisVector->begin(); thisObj != thisVector->end(); thisObj++)
	{
	  if(*thisObj != NULL) delete *thisObj;
	  *thisObj = NULL;
	}
    }

  cout << "getting CLs+b" << endl;

  if(isHybrid) res_hybrid_calculator->UseCLs(false);

  double clsplusb_hybrid_upperLimit_linear             = !isHybrid ? 0 : res_hybrid_calculator->UpperLimit();
  double clsplusb_hybrid_upperLimitError_linear        = !isHybrid ? 0 : res_hybrid_calculator->UpperLimitEstimatedError();
  double clsplusb_hybrid_expectedUpperLimit_linear     = !isHybrid ? 0 : res_hybrid_calculator->GetExpectedUpperLimit(0);
  double clsplusb_hybrid_expectedUpperLimitHigh_linear = !isHybrid ? 0 : res_hybrid_calculator->GetExpectedUpperLimit(1);
  double clsplusb_hybrid_expectedUpperLimitLow_linear  = !isHybrid ? 0 : res_hybrid_calculator->GetExpectedUpperLimit(-1);

  ///// --- END HYBRID STUFF

  
  //Saving things to a text file:
  TString outfileName("results/");
  outfileName+=modelName;
  outfileName+="_";
  outfileName+=m0;
  outfileName+="_";
  outfileName+=m12;
  if(isMeasured) outfileName+="_measured";
  else outfileName+="_expected";
  if(isHybrid) outfileName+="_hybrid";
  if(isFrequentist) outfileName+="_frequentist";
  outfileName+=".dat";

  ofstream outfile;  outfile.open(outfileName.Data());
  outfile << m0 ;  outfile << " ";
  outfile << m12;  outfile << " ";
  outfile << asymptoticLimit;
  if(isFrequentist){
    outfile << " ";
    outfile << upperLimit_linear;  outfile << " ";
    outfile << upperLimitError_linear;  outfile << " ";
    outfile << expectedUpperLimit_linear;  outfile << " ";
    outfile << expectedUpperLimitHigh_linear;  outfile << " ";
    outfile << expectedUpperLimitLow_linear;  outfile << " ";
    outfile << clsplusb_upperLimit_linear;  outfile << " ";
    outfile << clsplusb_upperLimitError_linear;  outfile << " ";
    outfile << clsplusb_expectedUpperLimit_linear;  outfile << " ";
    outfile << clsplusb_expectedUpperLimitHigh_linear;  outfile << " ";
    outfile << clsplusb_expectedUpperLimitLow_linear;
  }
  if(isHybrid){
    outfile << " ";
    outfile << hybrid_upperLimit_linear;  outfile << " ";
    outfile << hybrid_upperLimitError_linear;  outfile << " ";
    outfile << hybrid_expectedUpperLimit_linear;  outfile << " ";
    outfile << hybrid_expectedUpperLimitHigh_linear;  outfile << " ";
    outfile << hybrid_expectedUpperLimitLow_linear;  outfile << " ";
    outfile << clsplusb_hybrid_upperLimit_linear;  outfile << " ";
    outfile << clsplusb_hybrid_upperLimitError_linear;  outfile << " ";
    outfile << clsplusb_hybrid_expectedUpperLimit_linear;  outfile << " ";
    outfile << clsplusb_hybrid_expectedUpperLimitHigh_linear;  outfile << " ";
    outfile << clsplusb_hybrid_expectedUpperLimitLow_linear;
  }
  outfile << "\n";
  outfile.close();

  //outfileName=TString("results/TGraph.");
  //outfileName+=modelName;
  //outfileName+="_";
  //outfileName+=m0;
  //outfileName+="_";
  //outfileName+=m12;
  //if(isMeasured) outfileName+="_measured";
  //else outfileName+="_expected";
  //outfileName+=".dat";

  //ofstream outfileTGraph;  outfileTGraph.open(outfileName.Data());
  //outfileTGraph << m0 ;  outfileTGraph << " ";
  //outfileTGraph << m12;  outfileTGraph << " ";
  //outfileTGraph << upperLimit_linear;  outfileTGraph << "\n";

  //outfileTGraph.close();
  
  if(doPlots){

    //vector<TString> rangeNames;
    vector<pair< double, double > > ranges;

    double width = (poimax - poimin)/(nScanPoints-1);

    cout << "The min is " << poimin << endl;
    cout << "The max is " << poimax << endl;
    cout << "There are " << nScanPoints << " scan points." << endl;


    for(int iPoint = 0 ; iPoint < nScanPoints ; iPoint++)
      {
	double mid = poimin + iPoint*width;
	double low = mid - width/4;
	double high = mid + width/4;    
	TString rangeName("");
	rangeName+=iPoint;
	firstPOI->setRange(rangeName,low,high);
	poiHyp.setRange(rangeName+"poiHyp",low,high);
	ranges.push_back(make_pair(low,high));
	//rangeNames.push_back(rangeName);
      }
    
    RooPlot* signalFrame;
    RooPlot* backgroundFrame;
    TCanvas d("c","c",600,800);
    d.Divide(1,2);

    if(isFrequentist){
   
    signalFrame=firstPOI->frame(AutoRange(trueDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    trueDebugging.plotOn(signalFrame,Name("true"));
    backgroundFrame=poiHyp.frame(AutoRange(trueDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    trueDebugging.plotOn(backgroundFrame,Name("true"));
    d.cd(1); 
    signalFrame->Draw();
    d.cd(2);
    backgroundFrame->Draw();
    d.SaveAs("genTest.pdf");
    
    delete signalFrame;
    delete backgroundFrame;
    
    signalFrame=firstPOI->frame(AutoRange(uncondDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    uncondDebugging.plotOn(signalFrame,Name("true"));
    backgroundFrame=poiHyp.frame(AutoRange(uncondDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    uncondDebugging.plotOn(backgroundFrame,Name("true"));
    d.cd(1);
    signalFrame->Draw();
    d.cd(2);
    backgroundFrame->Draw();
    d.SaveAs("uncondTest.pdf");
    
    delete signalFrame;
    delete backgroundFrame;
    
    signalFrame=firstPOI->frame(AutoRange(condDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    condDebugging.plotOn(signalFrame,Name("true"));
    backgroundFrame=poiHyp.frame(AutoRange(condDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    condDebugging.plotOn(backgroundFrame,Name("true"));
    d.cd(1);
    signalFrame->Draw();
    d.cd(2);
    backgroundFrame->Draw();
    d.SaveAs("condTest.pdf");
    
    delete signalFrame;
    delete backgroundFrame;
    
    //vector<TString>::iterator rangeName = rangeNames.begin() ;
    for(vector<pair<double,double> >::iterator thisRange = ranges.begin() ; thisRange != ranges.end() /*|| rangeName!=rangeNames.end()*/ ; thisRange++/*, rangeName++*/)
      {
	
	TString trueCut = (TString(firstPOI->GetName())+">="+=thisRange->first)+(TString("&&")+firstPOI->GetName()+"<"+=thisRange->second);
	//TString zeroCut = (TString(firstPOI->GetName())+">="+=ranges.begin()->first)+(TString("&&")+firstPOI->GetName()+"<"+=ranges.begin()->second);
	TString zeroCut = (TString(firstPOI->GetName())+">="+=0)+(TString("&&")+firstPOI->GetName()+"<"+=ranges.begin()->first);
	TString hypoCut = (TString(poiHyp.GetName())+">="+=thisRange->first)+(TString("&&")+poiHyp.GetName()+"<"+=thisRange->second);
	TString sigCut = trueCut + "&&" + hypoCut;
	TString bgCut = zeroCut+"&&"+hypoCut;
	TString fileNamer = TString("mu_susy_")+=(thisRange->first+0.5*(thisRange->second-thisRange->first));
	
	
	RooAbsData* trueDebuggingSignal ;
	RooAbsData* uncondDebuggingSignal ;
	RooAbsData* condDebuggingSignal ;
	
	RooAbsData* trueDebuggingBackground ;
	RooAbsData* uncondDebuggingBackground ;
	RooAbsData* condDebuggingBackground ;
	
	trueDebuggingSignal = trueDebugging.reduce(SelectVars(dataSetVars),Cut(sigCut));
	uncondDebuggingSignal = uncondDebugging.reduce(SelectVars(dataSetVars),Cut(sigCut));
	condDebuggingSignal = condDebugging.reduce(SelectVars(dataSetVars),Cut(sigCut));
	
	trueDebuggingBackground = trueDebugging.reduce(SelectVars(dataSetVars),Cut(bgCut));
	uncondDebuggingBackground = uncondDebugging.reduce(SelectVars(dataSetVars),Cut(bgCut));
	condDebuggingBackground = condDebugging.reduce(SelectVars(dataSetVars),Cut(bgCut));
	
	
	RooLinkedListIter it = dataSetVars.iterator();
	RooRealVar* tmpPar  = NULL;
	while((tmpPar = (RooRealVar*)it.Next())){
	  double lo,hi,lowest,highest;
	  trueDebuggingSignal->getRange(*tmpPar,lowest,highest);
	  uncondDebuggingSignal->getRange(*tmpPar,lo,hi);
	  lowest = lowest<lo?lowest:lo;
	  highest = highest>hi?highest:hi;
	  condDebuggingSignal->getRange(*tmpPar,lo,hi);
	  lowest = lowest<lo?lowest:lo;
	  highest = highest>hi?highest:hi;
	  trueDebuggingBackground->getRange(*tmpPar,lo,hi);
	  lowest = lowest<lo?lowest:lo;
	  highest = highest>hi?highest:hi;
	  uncondDebuggingBackground->getRange(*tmpPar,lo,hi);
	  lowest = lowest<lo?lowest:lo;
	  highest = highest>hi?highest:hi;
	  condDebuggingBackground->getRange(*tmpPar,lo,hi);
	  lowest = lowest<lo?lowest:lo;
	  highest = highest>hi?highest:hi;
	  
	  lowest = lowest - 0.05*(highest - lowest);
	  highest = highest + 0.05*(highest - lowest);
	  
	  
	  backgroundFrame=tmpPar->frame(Range(lowest,highest),Bins(nToys/25),Title(TString("s+b Sample Distribution of ")+ tmpPar->GetName()));
	  trueDebuggingBackground->plotOn(backgroundFrame,Name("true"));
	  uncondDebuggingBackground->plotOn(backgroundFrame,LineColor(kRed),MarkerColor(kRed),Name("uncondFit"));
	  condDebuggingBackground->plotOn(backgroundFrame,LineColor(kBlue),MarkerColor(kBlue),Name("condFit"));
	  TLegend* backgroundLegend = new TLegend(.6,.75,1.,1.,TString("Background Only Distribution of ")+tmpPar->GetName());
	  backgroundFrame->getHist("true")->SetFillColor(kWhite);
	  backgroundFrame->getHist("uncondFit")->SetFillColor(kWhite);
	  backgroundFrame->getHist("condFit")->SetFillColor(kWhite);
	  backgroundFrame->getHist("true")->SetFillStyle(0);
	  backgroundFrame->getHist("uncondFit")->SetFillStyle(0);
	  backgroundFrame->getHist("condFit")->SetFillStyle(0);
	  backgroundLegend->AddEntry(backgroundFrame->findObject("true"),"True Distribution");
	  backgroundLegend->AddEntry(backgroundFrame->findObject("uncondFit"),"Fitted Distribution with Floating SUSY");
	  backgroundLegend->AddEntry(backgroundFrame->findObject("condFit"),"Fitted Distribution with Fixed SUSY");
	  backgroundLegend->SetFillColor(kWhite);
	  backgroundLegend->SetFillStyle(0);
	  
	  backgroundFrame->addObject(backgroundLegend);
	  
	  signalFrame=tmpPar->frame(Range(lowest,highest),Bins(nToys/25),Title(TString("s+b Sample Distribution of ")+ tmpPar->GetName()));
	  trueDebuggingSignal->plotOn(signalFrame,Name("true"));
	  uncondDebuggingSignal->plotOn(signalFrame,LineColor(kRed),MarkerColor(kRed),Name("uncondFit"));
	  condDebuggingSignal->plotOn(signalFrame,LineColor(kBlue),MarkerColor(kBlue),Name("condFit"));
	  TLegend* signalLegend = new TLegend(.6,.75,1.,1.,TString("Signal + Background Distribution of ")+tmpPar->GetName());
	  signalFrame->getHist("true")->SetFillColor(kWhite);
	  signalFrame->getHist("uncondFit")->SetFillColor(kWhite);
	  signalFrame->getHist("condFit")->SetFillColor(kWhite);
	  signalFrame->getHist("true")->SetFillStyle(0);
	  signalFrame->getHist("uncondFit")->SetFillStyle(0);
	  signalFrame->getHist("condFit")->SetFillStyle(0);
	  signalLegend->AddEntry(signalFrame->findObject("true"),"True Distribution");
	  signalLegend->AddEntry(signalFrame->findObject("uncondFit"),"Fitted Distribution with Floating SUSY");
	  signalLegend->AddEntry(signalFrame->findObject("condFit"),"Fitted Distribution with Fixed SUSY");
	  signalFrame->addObject(signalLegend);
	  signalLegend->SetFillColor(kWhite);
	  signalLegend->SetFillStyle(0);
	  d.cd(1);
	  backgroundFrame->Draw();
	  d.cd(2);
	  signalFrame->Draw();
	  d.SaveAs(TString("debugging_plots/")+fileNamer+"_"+tmpPar->GetName()+"_dist.pdf");
	}
	//cout << "Doing something with " << tmpPar->GetName() << " in range " << "Signal Yield Hypothesis of " << *rangeName << endl;
	delete backgroundFrame;
	delete signalFrame;
	delete trueDebuggingSignal;
	delete uncondDebuggingSignal;
	delete condDebuggingSignal;	 
	delete trueDebuggingBackground;
	delete uncondDebuggingBackground;
	delete condDebuggingBackground;
      } 
    }
    
    if(isHybrid){
    signalFrame=firstPOI->frame(AutoRange(trueHybridDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    trueHybridDebugging.plotOn(signalFrame,Name("true"));
    backgroundFrame=poiHyp.frame(AutoRange(trueHybridDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    trueHybridDebugging.plotOn(backgroundFrame,Name("true"));
    d.cd(1); 
    signalFrame->Draw();
    d.cd(2);
    backgroundFrame->Draw();
    d.SaveAs("genTest_hybrid.pdf");
    
    delete signalFrame;
    delete backgroundFrame;
    
    signalFrame=firstPOI->frame(AutoRange(uncondHybridDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    uncondHybridDebugging.plotOn(signalFrame,Name("true"));
    backgroundFrame=poiHyp.frame(AutoRange(uncondHybridDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    uncondHybridDebugging.plotOn(backgroundFrame,Name("true"));
    d.cd(1);
    signalFrame->Draw();
    d.cd(2);
    backgroundFrame->Draw();
    d.SaveAs("uncondTest_hybrid.pdf");
    
    delete signalFrame;
    delete backgroundFrame;
    
    signalFrame=firstPOI->frame(AutoRange(condHybridDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    condHybridDebugging.plotOn(signalFrame,Name("true"));
    backgroundFrame=poiHyp.frame(AutoRange(condHybridDebugging),Bins(100),Title(TString("s+b Sample Distribution of  for ")));
    condHybridDebugging.plotOn(backgroundFrame,Name("true"));
    d.cd(1);
    signalFrame->Draw();
    d.cd(2);
    backgroundFrame->Draw();
    d.SaveAs("condTest_hybrid.pdf");
    
    delete signalFrame;
    delete backgroundFrame;
    
    //vector<TString>::iterator rangeName = rangeNames.begin() ;
    for(vector<pair<double,double> >::iterator thisRange = ranges.begin() ; thisRange != ranges.end() /*|| rangeName!=rangeNames.end()*/ ; thisRange++/*, rangeName++*/)
      {
	
	TString trueCut = (TString(firstPOI->GetName())+">="+=thisRange->first)+(TString("&&")+firstPOI->GetName()+"<"+=thisRange->second);
	//TString zeroCut = (TString(firstPOI->GetName())+">="+=ranges.begin()->first)+(TString("&&")+firstPOI->GetName()+"<"+=ranges.begin()->second);
	TString zeroCut = (TString(firstPOI->GetName())+">="+=0)+(TString("&&")+firstPOI->GetName()+"<"+=ranges.begin()->first);
	TString hypoCut = (TString(poiHyp.GetName())+">="+=thisRange->first)+(TString("&&")+poiHyp.GetName()+"<"+=thisRange->second);
	TString sigCut = trueCut + "&&" + hypoCut;
	TString bgCut = zeroCut+"&&"+hypoCut;
	TString fileNamer = TString("mu_susy_hybrid")+=(thisRange->first+0.5*(thisRange->second-thisRange->first));
	
	
	RooAbsData* trueHybridDebuggingSignal ;
	RooAbsData* uncondHybridDebuggingSignal ;
	RooAbsData* condHybridDebuggingSignal ;
	
	RooAbsData* trueHybridDebuggingBackground ;
	RooAbsData* uncondHybridDebuggingBackground ;
	RooAbsData* condHybridDebuggingBackground ;
	
	trueHybridDebuggingSignal = trueHybridDebugging.reduce(SelectVars(dataSetVars),Cut(sigCut));
	uncondHybridDebuggingSignal = uncondHybridDebugging.reduce(SelectVars(dataSetVars),Cut(sigCut));
	condHybridDebuggingSignal = condHybridDebugging.reduce(SelectVars(dataSetVars),Cut(sigCut));
	
	trueHybridDebuggingBackground = trueHybridDebugging.reduce(SelectVars(dataSetVars),Cut(bgCut));
	uncondHybridDebuggingBackground = uncondHybridDebugging.reduce(SelectVars(dataSetVars),Cut(bgCut));
	condHybridDebuggingBackground = condHybridDebugging.reduce(SelectVars(dataSetVars),Cut(bgCut));
	
	
	RooLinkedListIter it = dataSetVars.iterator();
	RooRealVar* tmpPar  = NULL;
	while((tmpPar = (RooRealVar*)it.Next())){
	  double lo,hi,lowest,highest;
	  trueHybridDebuggingSignal->getRange(*tmpPar,lowest,highest);
	  uncondHybridDebuggingSignal->getRange(*tmpPar,lo,hi);
	  lowest = lowest<lo?lowest:lo;
	  highest = highest>hi?highest:hi;
	  condHybridDebuggingSignal->getRange(*tmpPar,lo,hi);
	  lowest = lowest<lo?lowest:lo;
	  highest = highest>hi?highest:hi;
	  trueHybridDebuggingBackground->getRange(*tmpPar,lo,hi);
	  lowest = lowest<lo?lowest:lo;
	  highest = highest>hi?highest:hi;
	  uncondHybridDebuggingBackground->getRange(*tmpPar,lo,hi);
	  lowest = lowest<lo?lowest:lo;
	  highest = highest>hi?highest:hi;
	  condHybridDebuggingBackground->getRange(*tmpPar,lo,hi);
	  lowest = lowest<lo?lowest:lo;
	  highest = highest>hi?highest:hi;
	  
	  lowest = lowest - 0.05*(highest - lowest);
	  highest = highest + 0.05*(highest - lowest);
	  
	  
	  backgroundFrame=tmpPar->frame(Range(lowest,highest),Bins(nToys/25),Title(TString("s+b Sample Distribution of ")+ tmpPar->GetName()));
	  trueHybridDebuggingBackground->plotOn(backgroundFrame,Name("true"));
	  uncondHybridDebuggingBackground->plotOn(backgroundFrame,LineColor(kRed),MarkerColor(kRed),Name("uncondFit"));
	  condHybridDebuggingBackground->plotOn(backgroundFrame,LineColor(kBlue),MarkerColor(kBlue),Name("condFit"));
	  TLegend* backgroundLegend = new TLegend(.6,.75,1.,1.,TString("Background Only Distribution of ")+tmpPar->GetName());
	  backgroundFrame->getHist("true")->SetFillColor(kWhite);
	  backgroundFrame->getHist("uncondFit")->SetFillColor(kWhite);
	  backgroundFrame->getHist("condFit")->SetFillColor(kWhite);
	  backgroundFrame->getHist("true")->SetFillStyle(0);
	  backgroundFrame->getHist("uncondFit")->SetFillStyle(0);
	  backgroundFrame->getHist("condFit")->SetFillStyle(0);
	  backgroundLegend->AddEntry(backgroundFrame->findObject("true"),"True Distribution");
	  backgroundLegend->AddEntry(backgroundFrame->findObject("uncondFit"),"Fitted Distribution with Floating SUSY");
	  backgroundLegend->AddEntry(backgroundFrame->findObject("condFit"),"Fitted Distribution with Fixed SUSY");
	  backgroundLegend->SetFillColor(kWhite);
	  backgroundLegend->SetFillStyle(0);
	  
	  backgroundFrame->addObject(backgroundLegend);
	  
	  signalFrame=tmpPar->frame(Range(lowest,highest),Bins(nToys/25),Title(TString("s+b Sample Distribution of ")+ tmpPar->GetName()));
	  trueHybridDebuggingSignal->plotOn(signalFrame,Name("true"));
	  uncondHybridDebuggingSignal->plotOn(signalFrame,LineColor(kRed),MarkerColor(kRed),Name("uncondFit"));
	  condHybridDebuggingSignal->plotOn(signalFrame,LineColor(kBlue),MarkerColor(kBlue),Name("condFit"));
	  TLegend* signalLegend = new TLegend(.6,.75,1.,1.,TString("Signal + Background Distribution of ")+tmpPar->GetName());
	  signalFrame->getHist("true")->SetFillColor(kWhite);
	  signalFrame->getHist("uncondFit")->SetFillColor(kWhite);
	  signalFrame->getHist("condFit")->SetFillColor(kWhite);
	  signalFrame->getHist("true")->SetFillStyle(0);
	  signalFrame->getHist("uncondFit")->SetFillStyle(0);
	  signalFrame->getHist("condFit")->SetFillStyle(0);
	  signalLegend->AddEntry(signalFrame->findObject("true"),"True Distribution");
	  signalLegend->AddEntry(signalFrame->findObject("uncondFit"),"Fitted Distribution with Floating SUSY");
	  signalLegend->AddEntry(signalFrame->findObject("condFit"),"Fitted Distribution with Fixed SUSY");
	  signalFrame->addObject(signalLegend);
	  signalLegend->SetFillColor(kWhite);
	  signalLegend->SetFillStyle(0);
	  d.cd(1);
	  backgroundFrame->Draw();
	  d.cd(2);
	  signalFrame->Draw();
	  d.SaveAs(TString("debugging_plots/")+fileNamer+"_"+tmpPar->GetName()+"_dist.pdf");
	}
	//cout << "Doing something with " << tmpPar->GetName() << " in range " << "Signal Yield Hypothesis of " << *rangeName << endl;
	delete backgroundFrame;
	delete signalFrame;
	delete trueHybridDebuggingSignal;
	delete uncondHybridDebuggingSignal;
	delete condHybridDebuggingSignal;	 
	delete trueHybridDebuggingBackground;
	delete uncondHybridDebuggingBackground;
	delete condHybridDebuggingBackground;
      } 
    }

    fDebug->Write();
    fDebug->Close();       
  }
}
