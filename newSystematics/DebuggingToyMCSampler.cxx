// @(#)root/roostats:$Id: ToyMCSampler.cxx,v 1.7 2011/11/30 14:40:49 owen Exp $
// Author: Sven Kreiss    June 2010
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//===== DEBUG
//#include "RooStats/DebuggingToyMCSampler.h"
#include "DebuggingToyMCSampler.h"
#include "TFile.h"
#include "TSystem.h"
//===== DEBUG

#ifndef ROO_MSG_SERVICE
#include "RooMsgService.h"
#endif

#ifndef ROO_DATA_HIST
#include "RooDataHist.h"
#endif

#ifndef ROO_REAL_VAR
#include "RooRealVar.h"
#endif

#include "TCanvas.h"
#include "RooPlot.h"
#include "RooRandom.h"

#include "RooStudyManager.h"
#include "RooStats/ToyMCStudy.h"
#include "RooSimultaneous.h"

#include "TMath.h"
#include "TTree.h"


ClassImp(RooStats::DebuggingToyMCSampler)

namespace RooStats {

  //  class NuisanceParametersSampler {
  //     // Helper for DebuggingToyMCSampler. Handles all of the nuisance parameter related
  //     // functions. Once instantiated, it gives a new nuisance parameter point
  //     // at each call to nextPoint(...).
  //  
  //     public:
  //        NuisanceParametersSampler(RooAbsPdf *prior=NULL, const RooArgSet *parameters=NULL, Int_t nToys=1000, Bool_t asimov=kFALSE) :
  //           fPrior(prior),
  //           fParams(parameters),
  //           fNToys(nToys),
  //           fExpected(asimov),
  //           fPoints(NULL),
  //           fIndex(0)
  //        {
  //           if(prior) Refresh();
  //        }
  //        virtual ~NuisanceParametersSampler() {
  //           if(fPoints) delete fPoints;
  //        }
  //  
  //        void NextPoint(RooArgSet& nuisPoint, Double_t& weight) {
  //           // Assigns new nuisance parameter point to members of nuisPoint.
  //           // nuisPoint can be more objects than just the nuisance
  //           // parameters.
  //  
  //           // check whether to get new set of nuisanceParPoints
  //           if (fIndex >= fNToys) {
  //              Refresh();
  //              fIndex = 0;
  //           }
  //  
  //           // get value
  //           nuisPoint =  *fPoints->get(fIndex++);
  //           weight = fPoints->weight();
  //  
  //           // check whether result will have any influence
  //           if(fPoints->weight() == 0.0) {
  //              oocoutI((TObject*)NULL,Generation) << "Weight 0 encountered. Skipping." << endl;
  //              NextPoint(nuisPoint, weight);
  //           }
  //        }
  //  
  //  
  //     protected:
  //  
  //        void Refresh() {
  //           // Creates the initial set of nuisance parameter points. It also refills the
  //           // set with new parameter points if called repeatedly. This helps with
  //           // adaptive sampling as the required number of nuisance parameter points
  //           // might increase during the run.
  //  
  //           if (!fPrior || !fParams) return;
  //  
  //           if (fPoints) delete fPoints;
  //  
  //           if (fExpected) {
  //              // UNDER CONSTRUCTION
  //              oocoutI((TObject*)NULL,InputArguments) << "Using expected nuisance parameters." << endl;
  //  
  //              int nBins = fNToys;
  //  
  //              // From FeldmanCousins.cxx:
  //              // set nbins for the POI
  //              TIter it2 = fParams->createIterator();
  //              RooRealVar *myarg2;
  //              while ((myarg2 = dynamic_cast<RooRealVar*>(it2.Next()))) {
  //                myarg2->setBins(nBins);
  //              }
  //  
  //  
  //              fPoints = fPrior->generateBinned(
  //                 *fParams,
  //                 RooFit::ExpectedData(),
  //                 RooFit::NumEvents(1) // for Asimov set, this is only a scale factor
  //              );
  //              if(fPoints->numEntries() != fNToys) {
  //                 fNToys = fPoints->numEntries();
  //                 oocoutI((TObject*)NULL,InputArguments) <<
  //                    "Adjusted number of toys to number of bins of nuisance parameters: " << fNToys << endl;
  //              }
  //  
  //  /*
  //              // check
  //              TCanvas *c1 = new TCanvas;
  //              RooPlot *p = dynamic_cast<RooRealVar*>(fParams->first())->frame();
  //              fPoints->plotOn(p);
  //              p->Draw();
  //              for(int x=0; x < fPoints->numEntries(); x++) {
  //                 fPoints->get(x)->Print("v");
  //                 cout << fPoints->weight() << endl;
  //              }
  //  */
  //  
  //           }else{
  //              oocoutI((TObject*)NULL,InputArguments) << "Using randomized nuisance parameters." << endl;
  //  
  //              fPoints = fPrior->generate(*fParams, fNToys);
  //           }
  //        }
  //  
  //  
  //     private:
  //        RooAbsPdf *fPrior;           // prior for nuisance parameters
  //        const RooArgSet *fParams;    // nuisance parameters
  //        Int_t fNToys;
  //        Bool_t fExpected;
  //  
  //        RooAbsData *fPoints;         // generated nuisance parameter points
  //        Int_t fIndex;                // current index in fPoints array
  //  };




DebuggingToyMCSampler::~DebuggingToyMCSampler() {
   if(fNuisanceParametersSampler) delete fNuisanceParametersSampler;
   ClearCache();
}



RooAbsData* DebuggingToyMCSampler::GenerateToyData(RooArgSet& paramPoint, double& weight, RooAbsPdf& pdf) const {
   // This method generates a toy data set for the given parameter point taking
   // global observables into account.
   // The values of the generated global observables remain in the pdf's variables.
   // They have to have those values for the subsequent evaluation of the
   // test statistics.

   if(!fObservables) {
      ooccoutE((TObject*)NULL,InputArguments) << "Observables not set." << endl;
      return NULL;
   }

   // assign input paramPoint
   RooArgSet* allVars = fPdf->getVariables();
   *allVars = paramPoint;


   // create nuisance parameter points
   if(!fNuisanceParametersSampler && fPriorNuisance && fNuisancePars)
      fNuisanceParametersSampler = new NuisanceParametersSampler(fPriorNuisance, fNuisancePars, fNToys, fExpectedNuisancePar);


   // generate global observables
   RooArgSet observables(*fObservables);
   if(fGlobalObservables  &&  fGlobalObservables->getSize()) {
      observables.remove(*fGlobalObservables);
      GenerateGlobalObservables(pdf);
   }

   // save values to restore later.
   // but this must remain after(!) generating global observables
   const RooArgSet* saveVars = (const RooArgSet*)allVars->snapshot();

   if(fNuisanceParametersSampler) { // use nuisance parameters?
      // Construct a set of nuisance parameters that has the parameters
      // in the input paramPoint removed. Therefore, no parameter in
      // paramPoint is randomized.
      // Therefore when a parameter is given (should be held fixed),
      // but is also in the list of nuisance parameters, the parameter
      // will be held fixed. This is useful for debugging to hold single
      // parameters fixed although under "normal" circumstances it is
      // randomized.
      RooArgSet allVarsMinusParamPoint(*allVars);
      allVarsMinusParamPoint.remove(paramPoint, kFALSE, kTRUE); // match by name

      // get nuisance parameter point and weight
      fNuisanceParametersSampler->NextPoint(allVarsMinusParamPoint, weight);
   }else{
      weight = -1.0;
   }

   RooAbsData *data = Generate(pdf, observables);

   if(fDebuggingData){
     RooArgSet saveSet(*allVars);
     saveSet = *data->get();
     fDebuggingData->add(saveSet);
   }

   // We generated the data with the randomized nuisance parameter (if hybrid)
   // but now we are setting the nuisance parameters back to where they were.
   *allVars = *saveVars;
   delete allVars;
   delete saveVars;

   return data;
}

}
