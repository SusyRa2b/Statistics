// @(#)root/roostats:$Id: HybridToyMCSampler.cxx 42273 2011-11-28 14:28:37Z moneta $
// Author: Lucas Winstrom and Sven Kreiss    June 2010
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "HybridToyMCSampler.h"

#ifndef ROO_MSG_SERVICE
#include "RooMsgService.h"
#endif

#ifndef ROO_DATA_HIST
#include "RooDataHist.h"
#endif

#ifndef ROO_REAL_VAR
#include "RooRealVar.h"
#endif

#include "RooStats/RooStatsUtils.h"

#include "TCanvas.h"
#include "RooPlot.h"
#include "RooRandom.h"

#include "RooStudyManager.h"
#include "RooStats/ToyMCStudy.h"
#include "RooSimultaneous.h"

#include "RooProfileLL.h"
#include "RooMinuit.h"

#include "TMath.h"
#include "TTree.h"

#include "RooMinimizer.h"
#include "Math/MinimizerOptions.h"

using namespace RooFit;


ClassImp(RooStats::HybridToyMCSampler)

namespace RooStats {

class NuisanceParametersSamplerAndFitter {
   // Helper for HybridToyMCSampler. Handles all of the nuisance parameter related
   // functions. Once instantiated, it gives a new nuisance parameter point
   // at each call to nextPoint(...).

   public:
       NuisanceParametersSamplerAndFitter(RooAbsPdf *prior=NULL, RooAbsPdf *fullPdf=NULL, RooAbsData* data=NULL , const RooArgSet *poi=NULL , const RooArgSet *parameters=NULL , Int_t nToys=1000, Bool_t asimov=kFALSE) :
         fPrior(prior),
	 fPdf((RooAbsPdf*)fullPdf->cloneTree()),
	 fData(data),
	 fPOI((RooArgSet*)poi->snapshot()),
         fParams(parameters),
         fNToys(nToys),
         fExpected(asimov),
         fPoints(NULL),
         fIndex(0)
      {
	if(prior) Refresh();
      }
      virtual ~NuisanceParametersSamplerAndFitter() {
         if(fPoints) delete fPoints;
	 delete fPdf;
	 delete fPOI;
      }

      void NextPoint(RooArgSet& nuisPoint, Double_t& weight) {
         // Assigns new nuisance parameter point to members of nuisPoint.
         // nuisPoint can be more objects than just the nuisance
         // parameters.

         // check whether to get new set of nuisanceParPoints
         if (fIndex >= fNToys) {
            Refresh();
            fIndex = 0;
         }


	 //New Code for Hybrid Hybrid Calculator
	 
	 if(fData)
	   {
	     RooArgSet * allVariables = fPdf->getParameters(*fData);
	     RemoveConstantParameters(allVariables);
	     *allVariables=*fPoints->get(fIndex++);
	     RooArgSet nuisancePlusPoi;

	     RooLinkedListIter it = fParams->iterator();
	     RooAbsArg*  tmpPar  = NULL;
	     RooRealVar* tmpParA = NULL;
	     while((tmpPar = (RooRealVar*)it.Next())){
	       tmpParA =  dynamic_cast<RooRealVar*>( allVariables->find(tmpPar->GetName()) );
	       if (tmpParA) nuisancePlusPoi.add(*tmpParA);
	       //cout << "Adding variable " << tmpParA->GetName() << " to nuisancePlusPoi" << endl;
	     }

	     it=fPOI->iterator();
	     while((tmpPar = (RooRealVar*)it.Next())){
	       tmpParA =  dynamic_cast<RooRealVar*>( allVariables->find(tmpPar->GetName()) );
	       if (tmpParA) nuisancePlusPoi.add(*tmpParA);
	       //cout << "Adding variable " << tmpParA->GetName() << " to nuisancePlusPoi" << endl;
	     }

	     nuisancePlusPoi=nuisPoint;
	     nuisancePlusPoi=*fPoints->get(fIndex++);

	     it = nuisancePlusPoi.iterator();
	     while((tmpParA = (RooRealVar*)it.Next())){
	       tmpParA->setConstant();
	       //cout << "Setting variable " << tmpParA->GetName() << " constant and to value " << tmpParA->getVal() << endl;
	     }

	     RooAbsCollection* allVarSnap = allVariables->snapshot();
	     //RooArgSet * allNuisance = fPrior->getVariables();
	     //RemoveConstantParameters(allNuisance);
	     //RooArgSet nuisancePlusPoi(*allNuisance);


	     //it = nuisancePlusPoi.iterator();
	     //while((tmpParA = (RooRealVar*)it.Next())){
	     //  tmpParA->setConstant(kFALSE);
	     //}
	     RooAbsReal* nll = fPdf->createNLL(*fData);
	     int status =0;
	     MinNLL(nll,status);

	     it.Reset();
	     while((tmpParA = (RooRealVar*)it.Next())){
	       tmpParA->setConstant();
	       //cout << "Setting variable " << tmpParA->GetName() << " not constant and to value " << tmpParA->getVal() << endl;
	     }
	     
	     //nuisancePlusPoi.Print("v");
	     //RooProfileLL* profile = (RooProfileLL*) nll->createProfile(nuisancePlusPoi);
	     //profile->setEvalErrorLoggingMode(RooAbsReal::CountErrors);
	     //profile->minuit()->setWarnLevel(-1);
	     //profile->getVal(); // this will do fit and set nuisance parameters to profiled values
	     // add nuisance parameters to parameter point
	     //cout << "========= First Values ==============" << endl;
	     //allVarSnap->Print("v");
	     //cout << "=====================================" << endl;
	     //cout << "========= Second Values =============" << endl;
	     //allVariables->Print("v");
	     //cout << "=====================================" << endl;
	     
	     nuisPoint=*allVariables;
	     weight = fPoints->weight();

	     delete nll;
	     delete allVariables;
	     delete allVarSnap;
	     
	     if(status != 0) NextPoint(nuisPoint, weight);

	   }
	 //End New Code
	 else{
         // get value
         nuisPoint =  *fPoints->get(fIndex++);
         weight = fPoints->weight();
	 }
         // check whether result will have any influence
         if(fPoints->weight() == 0.0) {
            oocoutI((TObject*)NULL,Generation) << "Weight 0 encountered. Skipping." << endl;
            NextPoint(nuisPoint, weight);
         }
      }


   protected:

      void Refresh() {
         // Creates the initial set of nuisance parameter points. It also refills the
         // set with new parameter points if called repeatedly. This helps with
         // adaptive sampling as the required number of nuisance parameter points
         // might increase during the run.

         if (!fPrior || !fParams) return;

         if (fPoints) delete fPoints;

         if (fExpected) {
            // UNDER CONSTRUCTION
            oocoutI((TObject*)NULL,InputArguments) << "Using expected nuisance parameters." << endl;

            int nBins = fNToys;

            // From FeldmanCousins.cxx:
            // set nbins for the POI
            TIter it2 = fParams->createIterator();
            RooRealVar *myarg2;
            while ((myarg2 = dynamic_cast<RooRealVar*>(it2.Next()))) {
              myarg2->setBins(nBins);
            }


            fPoints = fPrior->generate(
               *fParams,
               AllBinned(),
               ExpectedData(),
               NumEvents(1) // for Asimov set, this is only a scale factor
            );


            if(fPoints->numEntries() != fNToys) {
               fNToys = fPoints->numEntries();
               oocoutI((TObject*)NULL,InputArguments) <<
                  "Adjusted number of toys to number of bins of nuisance parameters: " << fNToys << endl;
            }

/*
            // check
            TCanvas *c1 = new TCanvas;
            RooPlot *p = dynamic_cast<RooRealVar*>(fParams->first())->frame();
            fPoints->plotOn(p);
            p->Draw();
            for(int x=0; x < fPoints->numEntries(); x++) {
               fPoints->get(x)->Print("v");
               cout << fPoints->weight() << endl;
            }
*/
         }else{
            oocoutI((TObject*)NULL,InputArguments) << "Using randomized nuisance parameters." << endl;
            fPoints = fPrior->generate(*fParams, fNToys);
         }
      }

  void MinNLL(RooAbsReal* NLL,int& status) {
    //find minimum of NLL using RooMinimizer
    
    TString fMinimizer;
    Int_t fStrategy;
    Double_t fTolerance; 
    Int_t fPrintLevel;
    fMinimizer=::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    fStrategy=::ROOT::Math::MinimizerOptions::DefaultStrategy();
    // avoid default tolerance to be too small (1. is default in RooMinimizer)
    fTolerance=TMath::Max(1.,::ROOT::Math::MinimizerOptions::DefaultTolerance());
    fPrintLevel=::ROOT::Math::MinimizerOptions::DefaultPrintLevel();

   RooMinimizer minim(*NLL);
   minim.setStrategy(fStrategy);
   minim.setEps(fTolerance);
   //LM: RooMinimizer.setPrintLevel has +1 offset - so subtruct  here -1 + an extra -1 
   int level = (fPrintLevel == 0) ? -1 : fPrintLevel -2;
   minim.setPrintLevel(level);
   minim.setEps(fTolerance);
   // this cayses a memory leak
   minim.optimizeConst(true); 
   for (int tries = 0, maxtries = 4; tries <= maxtries; ++tries) {
      //	 status = minim.minimize(fMinimizer, ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
      status = minim.minimize(fMinimizer, "Minimize");
      if (status == 0) {  
         break;
      } else {
         if (tries > 1) {
            printf("    ----> Doing a re-scan first\n");
            minim.minimize(fMinimizer,"Scan");
         }
         if (tries > 2) {
            printf("    ----> trying with strategy = 1\n");
            minim.setStrategy(1);
         }
      }
   }
}

   private:
      RooAbsPdf *fPrior;           // prior for nuisance parameters
      RooAbsPdf *fPdf;            // full Pdf parameters
      RooAbsData *fData;
      const RooArgSet *fPOI;
      const RooArgSet *fParams;    // nuisance parameters
      Int_t fNToys;
      Bool_t fExpected;

      RooAbsData *fPoints;         // generated nuisance parameter points
      Int_t fIndex;                // current index in fPoints array
};


HybridToyMCSampler::~HybridToyMCSampler() {

   if(fNuisanceParametersSamplerAndFitter) delete fNuisanceParametersSamplerAndFitter;
}

//SamplingDistribution* HybridToyMCSampler::GetSamplingDistribution(RooArgSet& paramPointIn) {
//   // Use for serial and parallel runs.
//  //cout << "The paramPointIn is"<<endl;
//  //paramPointIn.Print("v");
//  
//  //Modification to get a new variable definition for each parameter point for sampling.
//
//  RooArgSet* allVars = NULL;
//  RooArgSet* allVarsSnap = NULL;
//  RooArgSet* paramPointInSnap = NULL;
//  bool doTest(false);
//  if(doTest && fFitToData)
//    {
//      //paramPointIn.Print("v");  
//      allVars = fPdf->getVariables();
//      allVarsSnap = (RooArgSet*) allVars->snapshot();
//      *allVars=paramPointIn;
//      RooLinkedListIter it = fNullPOI->iterator();
//      RooAbsArg*  tmpPar  = NULL;
//      RooRealVar* tmpParA = NULL;
//      while((tmpPar = (RooRealVar*)it.Next())){
//	tmpParA =  dynamic_cast<RooRealVar*>( allVars->find(tmpPar->GetName()) );
//	if (tmpParA) tmpParA->setConstant();
//	//cout << "Setting variable " << tmpParA->GetName() << " constant and to value " << tmpParA->getVal() << endl;
//      }
//      fPdf->fitTo(*fFitToData);
//      it.Reset();
//      while((tmpPar = (RooRealVar*)it.Next())){
//	tmpParA =  dynamic_cast<RooRealVar*>( allVars->find(tmpPar->GetName()) );
//	if (tmpParA) tmpParA->setConstant(kFALSE);
//    //cout << "Setting variable " << tmpParA->GetName() << " not constant with value " << tmpParA->getVal() << endl;
//      }
//      paramPointInSnap =  (RooArgSet*) paramPointIn.snapshot();
//      paramPointIn = *allVars ;
//      //paramPointIn.Print("v");
//    }
//  SamplingDistribution *result = ToyMCSampler::GetSamplingDistribution(paramPointIn) ;
//
//  if(allVarsSnap)
//    {
//      *allVars = *allVarsSnap;
//      delete allVarsSnap;
//    }
//  if(paramPointInSnap)
//    {
//      paramPointIn = *paramPointInSnap;
//    }
//
//   return result;
//}

RooAbsData* HybridToyMCSampler::GenerateToyData(RooArgSet& paramPoint, double& weight) const {
   // This method generates a toy data set for the given parameter point taking
   // global observables into account.
   // The values of the generated global observables remain in the pdf's variables.
   // They have to have those values for the subsequent evaluation of the
   // test statistics.

   if(!fObservables) {
      ooccoutE((TObject*)NULL,InputArguments) << "Observables not set." << endl;
      return NULL;
   }

   if(fImportanceDensity) {
      oocoutW((TObject*)NULL,InputArguments) << "HybridToyMCSampler: importance density given but ignored for generating toys." << endl;
   }

   // assign input paramPoint
   RooArgSet* allVars = fPdf->getVariables();
   *allVars = paramPoint;


   //create nuisance parameter points
   //if(fPriorNuisance)
   //  {
   //    RooArgSet* tempNuisanceVars = fPriorNuisance->getVariables();
   //    if(fNuisancePars) delete fNuisancePars;
   //    fNuisancePars = tempNuisanceVars;
   //  }
   //fNuisancePars->Print();
   if(!fNuisanceParametersSamplerAndFitter && fPriorNuisance && fNuisancePars)
     fNuisanceParametersSamplerAndFitter = new NuisanceParametersSamplerAndFitter(fPriorNuisance, fPdf, fFitToData, fNullPOI, fNuisancePars, fNToys, fExpectedNuisancePar);


   // generate global observables
   RooArgSet observables(*fObservables);
   if(fGlobalObservables  &&  fGlobalObservables->getSize()) {
      observables.remove(*fGlobalObservables);
      GenerateGlobalObservables();
   }

   // save values to restore later.
   // but this must remain after(!) generating global observables
   const RooArgSet* saveVars = (const RooArgSet*)allVars->snapshot();

   if(fNuisanceParametersSamplerAndFitter) { // use nuisance parameters?
      // Construct a set of nuisance parameters that has the parameters
      // in the input paramPoint removed. Therefore, no parameter in
      // paramPoint is randomized.
      // Therefore when a parameter is given (should be held fixed),
      // but is also in the list of nuisance parameters, the parameter
      // will be held fixed. This is useful for debugging to hold single
      // parameters fixed although under "normal" circumstances it is
      // randomized.
      //RooArgSet nuisanceVars(*allVars);
      //RemoveConstantParameters(&nuisanceVars);
      //RooArgSet nonNuisanceVars(nuisanceVars);
      //nonNuisanceVars.remove(*fNuisancePars,kFALSE,kTRUE);
      //nuisanceVars.remove(nonNuisanceVars,kFALSE,kTRUE);
      //RooArgSet* nuisanceSnap = (RooArgSet*)nuisanceVars.snapshot();
      //RooArgSet* nonNuisanceSnap = (RooArgSet*)nonNuisanceVars.snapshot();
      RooArgSet allVarsMinusParamPoint(*allVars);
      allVarsMinusParamPoint.remove(paramPoint, kFALSE, kTRUE); // match by name
      //allVarsMinusParamPoint.Print("v");
      fNuisanceParametersSamplerAndFitter->NextPoint(allVarsMinusParamPoint, weight);
      //allVarsMinusParamPoint.Print("v");
      //nonNuisanceSnap->Print("v");
      //nonNuisanceVars.Print("v");
      //nuisanceSnap->Print("v");
      //nuisanceVars.Print("v");

   }else{
      weight = -1.0;
   }

   RooAbsData *data = Generate(*fPdf, observables);

   // We generated the data with the randomized nuisance parameter (if hybrid)
   // but now we are setting the nuisance parameters back to where they were.
   *allVars = *saveVars;
   delete allVars;
   delete saveVars;

   return data;
}

} // end namespace RooStats
