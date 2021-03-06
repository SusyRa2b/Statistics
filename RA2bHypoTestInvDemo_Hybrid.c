// Standard tutorial macro for performing an inverted  hypothesis test 
//
// This macro will perform a scan of tehe p-values for computing the limit
// 

#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH2F.h"

#include "RooStats/HybridCalculator.h"
#include "RooStats/HypoTestPlot.h"


//-------------
//#include "HypoTestCalculatorGeneric.h"
//#include "RooStats/ToyMCSampler.h"
#include "HypoTestCalculatorGeneric.cxx"
//#include "ToyMCSampler.cxx"
#include "HybridToyMCSampler.h"
//-------------

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

#include "RooStats/FrequentistCalculator.h"

using namespace RooFit;
using namespace RooStats;


bool plotHypoTestResult = true; 
bool useProof = false;
bool optimize = false;
bool writeResult = false;
int nworkers = 1;


// internal routine to run the inverter
HypoTestInverterResult * RunInverter(RooWorkspace * w, const char * modelSBName, const char * modelBName, const char * dataName,
                                     int type,  int testStatType, int npoints, double poimin, double poimax, int ntoys, bool useCls );

vector<TObject*> makeStatDistPlot( HypoTestResult* theResult );


void RA2bHypoTestInvDemo(const char * fileName =0,
			 const char * wsName = "combined",
			 const char * modelSBName = "ModelConfig",
			 const char * modelBName = "",
			 const char * dataName = "obsData",                 
			 int calculatorType = 0,
			 int testStatType = 3, 
			 bool useCls = true ,  
			 int npoints = 5,   
			 double poimin = 0,  
			 double poimax = 5, 
			 int ntoys=1000,
			 int mgl = -1,
			 int mlsp = -1,
			 const char * outFileName = "test")    
{
/*

   Other Parameter to pass in tutorial
   apart from standard for filename, ws, modelconfig and data

    type = 0 Freq calculator 
    type = 1 Hybrid 

    testStatType = 0 LEP
                 = 1 Tevatron 
                 = 2 Profile Likelihood
                 = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)

    useCLs          scan for CLs (otherwise for CLs+b)    

    npoints:        number of points to scan , for autoscan set npoints = -1 

    poimin,poimax:  min/max value to scan in case of fixed scans 
                    (if min >= max, try to find automatically)                           

    ntoys:         number of toys to use 

    extra options are available as global paramters of the macro. They are: 

    plotHypoTestResult   plot result of tests at each point (TS distributions) 
    useProof = true;
    writeResult = true;
    nworkers = 4;


   */

   if (fileName==0) { 
      fileName = "results/example_combined_GaussExample_model.root";
      std::cout << "Use standard file generated with HistFactory :" << fileName << std::endl;
   }
   TFile * file = new TFile(fileName); 

   RooWorkspace * w = dynamic_cast<RooWorkspace*>( file->Get(wsName) );

   HypoTestInverterResult * r = 0; 
   std::cout << w << "\t" << fileName << std::endl;
   if (w != NULL) {
      r = RunInverter(w, modelSBName, modelBName, dataName, calculatorType, testStatType, npoints, poimin, poimax,  ntoys, useCls );    
      if (!r) { 
         std::cerr << "Error running the HypoTestInverter - Exit " << std::endl;
         return;          
      }
   }
   else 
   { 
      // case workspace is not present look for the inverter result
      std::cout << "Reading an HypoTestInverterResult with name " << wsName << " from file " << fileName << std::endl;
      r = dynamic_cast<HypoTestInverterResult*>( file->Get(wsName) ); //
      if (!r) { 
         std::cerr << "File " << fileName << " does not contain a workspace or an HypoTestInverterResult - Exit " 
                   << std::endl;
         file->ls();
         return; 
      }
   }		
      		


   printf("\n\n") ;
   HypoTestResult* htr = r->GetResult(0) ;
   printf("  Data value for test stat : %7.3f\n", htr->GetTestStatisticData() ) ;
   printf("  CLsplusb : %9.4f\n", r->CLsplusb(0) ) ;
   printf("  CLb      : %9.4f\n", r->CLb(0) ) ;
   printf("  CLs      : %9.4f\n", r->CLs(0) ) ;
   printf("\n\n") ;
   cout << flush ;

   double upperLimit = r->UpperLimit();
   double ulError = r->UpperLimitEstimatedError();


   std::cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << std::endl;
 
   const int nEntries = r->ArraySize();


   const char *  typeName = (calculatorType == 0) ? "Frequentist" : "Hybrid";
   const char * resultName = (w) ? w->GetName() : r->GetName();
   TString plotTitle = TString::Format("%s CL Scan for workspace %s",typeName,resultName);
   HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot",plotTitle,r);
   TCanvas* c1 = new TCanvas() ;
   plot->Draw("CLb 2CL");  // plot all and Clb
   c1->Update() ;
   c1->SaveAs("cls-canv1.png") ;
   c1->SaveAs("cls-canv1.pdf") ;

   if (plotHypoTestResult) { 
     vector<vector<TObject*> > objPointerVectors;
      TCanvas * c2 = new TCanvas();
      c2->Divide( 2, TMath::Ceil(nEntries/2));
      for (int i=0; i<nEntries; i++) {
         //c2->cd(i+1);
         //SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
         //pl->SetLogYaxis(true);
         //pl->Draw();
	c2->cd(i+1)->SetLogy(1);
	objPointerVectors.push_back(makeStatDistPlot(r->GetResult(i)));
      }
      c2->Update() ;
      c2->SaveAs("cls-canv2.png") ;
      c2->SaveAs("cls-canv2.pdf") ;

      for(vector<vector<TObject*> >::iterator thisVector = objPointerVectors.begin(); thisVector != objPointerVectors.end() ; thisVector++ )
	{
	  for(vector<TObject*>::iterator thisObj = thisVector->begin(); thisObj != thisVector->end(); thisObj++)
	    {
	      if(! *thisObj ) delete *thisObj;
	    }
	}
      
   }


   std::cout << " expected limit (median) " <<  r->GetExpectedUpperLimit(0) << std::endl;
   std::cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << std::endl;
   std::cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << std::endl;


   // save 2d histograms bin to file

   TH2F *result = new TH2F("result","result",22,100,1200,23,50,1200); 
   TH2F *exp_res = new TH2F("exp_res","exp_res",22,100,1200,23,50,1200); 
   TH2F *exp_res_minus = new TH2F("exp_res_minus","exp_res_minus",22,100,1200,23,50,1200); 
   TH2F *exp_res_plus = new TH2F("exp_res_plus","exp_res_plus",22,100,1200,23,50,1200); 

   result->Fill(mgl,mlsp,upperLimit);
   exp_res->Fill(mgl,mlsp,r->GetExpectedUpperLimit(0));
   exp_res_minus->Fill(mgl,mlsp,r->GetExpectedUpperLimit(-1));
   exp_res_plus->Fill(mgl,mlsp,r->GetExpectedUpperLimit(1));


   TFile *f = new TFile(outFileName,"RECREATE");
   f->cd();

   result->Write();
   exp_res->Write();
   exp_res_minus->Write();
   exp_res_plus->Write();

   f->Close();


   if (w != NULL && writeResult) {

      // write to a file the results
      const char *  calcType = (calculatorType == 0) ? "Freq" : "Hybr";
      const char *  limitType = (useCls) ? "CLs" : "Cls+b";
      const char * scanType = (npoints < 0) ? "auto" : "grid";
      TString resultFileName = TString::Format("%s_%s_%s_ts%d_",calcType,limitType,scanType,testStatType);      
      resultFileName += fileName;
      
      TFile * fileOut = new TFile(resultFileName,"RECREATE");
      r->Write();
      fileOut->Close();                                                                     
   }   

}

//==================================================================================================================

// internal routine to run the inverter
HypoTestInverterResult *  RunInverter(RooWorkspace * w, const char * modelSBName, const char * modelBName, 
                                      const char * dataName, int type,  int testStatType, 
                                      int npoints, double poimin, double poimax, 
                                      int ntoys, bool useCls ) 
{

   std::cout << "Running HypoTestInverter on the workspace " << w->GetName() << std::endl;

   w->Print();


   RooAbsData * data = w->data(dataName); 
   if (!data) { 
      Error("RA2bHypoTestDemo","Not existing data %s",dataName);
      return 0;
   }
   else 
      std::cout << "Using data set " << dataName << std::endl;

   
   // get models from WS
   // get the modelConfig out of the file
   ModelConfig* bModel = (ModelConfig*) w->obj(modelBName);
   ModelConfig* sbModel = (ModelConfig*) w->obj(modelSBName);

   if (!sbModel) {
      Error("RA2bHypoTestDemo","Not existing ModelConfig %s",modelSBName);
      return 0;
   }
   // check the model 
   if (!sbModel->GetPdf()) { 
      Error("RA2bHypoTestDemo","Model %s has no pdf ",modelSBName);
      return 0;
   }
   if (!sbModel->GetParametersOfInterest()) {
      Error("RA2bHypoTestDemo","Model %s has no poi ",modelSBName);
      return 0;
   }
   if (!sbModel->GetParametersOfInterest()) {
      Error("RA2bHypoTestInvDemo","Model %s has no poi ",modelSBName);
      return 0;
   }
   if (!sbModel->GetSnapshot() ) { 
      Info("RA2bHypoTestInvDemo","Model %s has no snapshot  - make one using model poi",modelSBName);
      sbModel->SetSnapshot( *sbModel->GetParametersOfInterest() );
   }


   if (!bModel || bModel == sbModel) {
      Info("RA2bHypoTestInvDemo","The background model %s does not exist",modelBName);
      Info("RA2bHypoTestInvDemo","Copy it from ModelConfig %s and set POI to zero",modelSBName);
      bModel = (ModelConfig*) sbModel->Clone();
      bModel->SetName(TString(modelSBName)+TString("_with_poi_0"));      
      RooRealVar * var = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
      if (!var) return 0;
      double oldval = var->getVal();
      var->setVal(0);
      bModel->SetSnapshot( RooArgSet(*var)  );
      var->setVal(oldval);
   }
   else { 
      if (!bModel->GetSnapshot() ) { 
         Info("RA2bHypoTestInvDemo","Model %s has no snapshot  - make one using model poi and 0 values ",modelBName);
         RooRealVar * var = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
         if (var) { 
            double oldval = var->getVal();
            var->setVal(0);
            bModel->SetSnapshot( RooArgSet(*var)  );
            var->setVal(oldval);
         }
         else { 
            Error("RA2bHypoTestInvDemo","Model %s has no valid poi",modelBName);
            return 0;
         }         
      }
   }


   SimpleLikelihoodRatioTestStat slrts(*sbModel->GetPdf(),*bModel->GetPdf());
   if (sbModel->GetSnapshot()) slrts.SetNullParameters(*sbModel->GetSnapshot());
   if (bModel->GetSnapshot()) slrts.SetAltParameters(*bModel->GetSnapshot());

   // ratio of profile likelihood - need to pass snapshot for the alt
   RatioOfProfiledLikelihoodsTestStat 
      ropl(*sbModel->GetPdf(), *bModel->GetPdf(), bModel->GetSnapshot());
   ropl.SetSubtractMLE(false);
   
   //MyProfileLikelihoodTestStat profll(*sbModel->GetPdf());
   ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
   if (testStatType == 3) profll.SetOneSided(1);
   if (optimize) profll.SetReuseNLL(true);

   TestStatistic * testStat = &slrts;
   if (testStatType == 1) testStat = &ropl;
   if (testStatType == 2 || testStatType == 3) testStat = &profll;
  
   HybridToyMCSampler poiSampler(*testStat,ntoys);
   HypoTestCalculatorGeneric *  hc = 0;
   if (type == 0) hc = new FrequentistCalculator(*data, *bModel, *sbModel,&poiSampler);
   else hc = new HybridCalculator(*data, *bModel, *sbModel, &poiSampler);

   ToyMCSampler *toymcs = (ToyMCSampler*)hc->GetTestStatSampler();
   //=== DEBUG
   ///// toymcs->SetWS( w ) ;
   //=== DEBUG
   toymcs->SetNEventsPerToy(1);
   toymcs->SetTestStatistic(testStat);
   if (optimize) toymcs->SetUseMultiGen(true);


   if (type == 1) { 
      HybridCalculator *hhc = (HybridCalculator*) hc;
      hhc->SetToys(ntoys,ntoys); 

      // check for nuisance prior pdf 
      if (bModel->GetPriorPdf() && sbModel->GetPriorPdf() ) {
         hhc->ForcePriorNuisanceAlt(*bModel->GetPriorPdf());
         hhc->ForcePriorNuisanceNull(*sbModel->GetPriorPdf());
	 ((HybridToyMCSampler*)hhc->GetTestStatSampler())->SetFitToData(data);

      }
      else {
         if (bModel->GetNuisanceParameters() || sbModel->GetNuisanceParameters() ) {
            Error("RA2bHypoTestInvDemo","Cannnot run Hybrid calculator because no prior on the nuisance parameter is specified");
            return 0;
         }
      }
   } 
   else 
      ((FrequentistCalculator*) hc)->SetToys(ntoys,ntoys); 

   // Get the result
   RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);


   TStopwatch tw; tw.Start(); 
   const RooArgSet * poiSet = sbModel->GetParametersOfInterest();
   RooRealVar *poi = (RooRealVar*)poiSet->first();

   // fit the data first
   sbModel->GetPdf()->fitTo(*data);
   double poihat  = poi->getVal();


   HypoTestInverter calc(*hc);
   calc.SetConfidenceLevel(0.95);

   calc.UseCLs(useCls);
   calc.SetVerbose(true);

   // can speed up using proof-lite
   if (useProof && nworkers > 1) { 
      ProofConfig pc(*w, nworkers, "", kFALSE);
      toymcs->SetProofConfig(&pc);    // enable proof
   }


   printf(" npoints = %d, poimin = %7.2f, poimax = %7.2f\n\n", npoints, poimin, poimax ) ;
   cout << flush ;

   if ( npoints==1 ) {

      std::cout << "Evaluating one point : " << poimax << std::endl;
      calc.RunOnePoint(poimax);

   } else if (npoints > 0) {
      if (poimin >= poimax) { 
         // if no min/max given scan between MLE and +4 sigma 
         poimin = int(poihat);
         poimax = int(poihat +  4 * poi->getError());
      }
      std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
      calc.SetFixedScan(npoints,poimin,poimax);
   }
   else { 
      //poi->setMax(10*int( (poihat+ 10 *poi->getError() )/10 ) );
      std::cout << "Doing an  automatic scan  in interval : " << poi->getMin() << " , " << poi->getMax() << std::endl;
   }

   //cout << "Setting New Nuisance Parameters for Test" << endl;
   //RooArgSet * newNuisances = sbModel->GetPdf()->getParameters(*data);
   //RemoveConstantParameters(newNuisances);
   //newNuisances->remove(*poiSet);
   //sbModel->SetNuisanceParameters(*newNuisances);
   //bModel->SetNuisanceParameters(*newNuisances);

   cout << "\n\n right before calc.GetInterval(), ntoys = " << ntoys << " \n\n" << flush ;
   HypoTestInverterResult * r = calc.GetInterval();


   return r; 
}

void ReadResult(const char * fileName, const char * resultName="", bool useCLs=true) { 
   // read a previous stored result from a file given the result name

   RA2bHypoTestInvDemo(fileName, resultName,"","","",0,0,useCLs);
}

int main() {
   RA2bHypoTestInvDemo();
}

vector<TObject*> makeStatDistPlot( HypoTestResult* theResult )
{
  vector<TObject*> theObjs;
  if(!theResult) return theObjs;
  TH1F* nullHypBelow(0),*nullHypAbove(0),*altHypBelow(0),*altHypAbove(0);

  const SamplingDistribution *null = theResult->GetNullDistribution();
  const SamplingDistribution *alt = theResult->GetAltDistribution();

  vector<Double_t> nullDist = null->GetSamplingDistribution();
  vector<Double_t> altDist = alt->GetSamplingDistribution();

  //cout << "Array sizes:" << endl;
  //cout << nullDist.size() << endl;
  //cout << altDist.size() << endl;

  //if( nullDist.size() == 0 || altDist.size() == 0 ) return theObjs ;

  Double_t nullmin(0) ;
  if( nullDist.size() != 0 ) nullmin = *(min_element(nullDist.begin(), nullDist.end())) ;
  Double_t nullmax(0) ;
  if( nullDist.size() != 0 ) nullmax = *(max_element(nullDist.begin(), nullDist.end())) ;
  
  //cout << "The nullmax is " << nullmax << endl;
  //cout << "The nullmin is " << nullmin << endl;

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

  vector<Double_t> nullWeightDist = null->GetSampleWeights();
  bool nullSizeComp = nullWeightDist.size() == nullDist.size();
  vector<Double_t> altWeightDist = alt->GetSampleWeights();
  bool altSizeComp = altWeightDist.size() == altDist.size();

  vector<Double_t>::iterator thisWeight;

  if(nullSizeComp) thisWeight = nullWeightDist.begin();

  for(vector<Double_t>::iterator thisPoint = nullDist.begin() ; thisPoint != nullDist.end() ; thisPoint++ )
    {
      double weight=1.;
      if(nullSizeComp){ 
	weight = *thisWeight;
	thisWeight++;
      }
      if(*thisPoint != *thisPoint) continue ;
      if(*thisPoint == numeric_limits<Double_t>::infinity() ) continue;
      if(*thisPoint == -numeric_limits<Double_t>::infinity() ) continue;
      if(*thisPoint < xmid) nullHypBelow->Fill(*thisPoint);
      else nullHypAbove->Fill(*thisPoint,weight);
    }

  if(altSizeComp) thisWeight = altWeightDist.begin();

  for(vector<Double_t>::iterator thisPoint = altDist.begin() ; thisPoint != altDist.end() ; thisPoint++ )
    {
      double weight=1.;
      if(altSizeComp){ 
	weight = *thisWeight;
	thisWeight++;
      }
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

  //cout << "The xmax is " << xmax << endl;
  //cout << "The xmin is " << xmin << endl;
  //cout << "The ymax is " << ymax << endl;
  //cout << "The xmid is " << xmid << endl;


  return theObjs;
}
