// @(#)root/roostats:$Id: ToyMCSampler.cxx 43199 2012-03-01 20:17:42Z moneta $
// Author: Sven Kreiss    June 2010
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

#include "TCanvas.h"
#include "RooPlot.h"
#include "RooRandom.h"

#include "RooStudyManager.h"
#include "RooStats/ToyMCStudy.h"
#include "RooSimultaneous.h"

#include "TMath.h"
#include "RooProfileLL.h"
#include "RooMinuit.h"

#include "RooStats/RooStatsUtils.h"
#include "TTree.h"

#include "RooMinimizer.h"
#include "Math/MinimizerOptions.h"


using namespace RooFit;
using namespace RooStats;

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
	 if (fIndex >= fNToys || fIndex >= fPoints->numEntries() ) {
	    cout << "fIndex refreshed" << endl;
            Refresh();
            fIndex = 0;
         }


	 //New Code for Hybrid Hybrid Calculator
	 
	 if(fData)
	   {

	     //Setting up custom original values to help with psychotic limit setting behavior

	     RooArgSet everyPossibleThing;
	     fPdf->treeNodeServerList(&everyPossibleThing);
	     fPdf->branchNodeServerList(&everyPossibleThing);
	     fPdf->leafNodeServerList(&everyPossibleThing);

	     ////everyPossibleThing.Print("v");
	     //
	     //double mu_Znn_sig = 0.5*(everyPossibleThing.getRealValue("Nsb_ee") * everyPossibleThing.getRealValue("fsig_ee") * everyPossibleThing.getRealValue("knn_ee_sig") * everyPossibleThing.getRealValue("znnoverll_bfratio") /( everyPossibleThing.getRealValue("acc_ee_sig") * everyPossibleThing.getRealValue("eff_ee") ) + everyPossibleThing.getRealValue("Nsb_mm") * everyPossibleThing.getRealValue("fsig_mm") * everyPossibleThing.getRealValue("knn_mm_sig") * everyPossibleThing.getRealValue("znnoverll_bfratio") /( everyPossibleThing.getRealValue("acc_mm_sig") * everyPossibleThing.getRealValue("eff_mm") ) );
	     //
	     //double mu_Znn_sb = 0.5*(everyPossibleThing.getRealValue("Nsb_ee") * everyPossibleThing.getRealValue("fsb_ee") * everyPossibleThing.getRealValue("knn_ee_sb") * everyPossibleThing.getRealValue("znnoverll_bfratio") /( everyPossibleThing.getRealValue("acc_ee_sb") * everyPossibleThing.getRealValue("eff_ee") ) + everyPossibleThing.getRealValue("Nsb_mm") * everyPossibleThing.getRealValue("fsb_mm") * everyPossibleThing.getRealValue("knn_mm_sb") * everyPossibleThing.getRealValue("znnoverll_bfratio") /( everyPossibleThing.getRealValue("acc_mm_sb") * everyPossibleThing.getRealValue("eff_mm") ) );
	     //
	     //double mu_qcd_sig = max(everyPossibleThing.getRealValue("Nsig_ldp") - (everyPossibleThing.getRealValue("eff_sf_sig_ldp") * ( everyPossibleThing.getRealValue("sf_mc") * (everyPossibleThing.getRealValue("mu_ttwj_sig_ldp") + everyPossibleThing.getRealValue("mu_znn_sig_ldp") + everyPossibleThing.getRealValue("mu_ewo_sig_ldp") ) + everyPossibleThing.getRealValue("mu_susy_sig_ldp") ) ) , 0. ) * everyPossibleThing.getRealValue("sf_qcd_sig") * everyPossibleThing.getRealValue("Rlsig_passfail");
	     //
	     //double mu_qcd_sb = max(everyPossibleThing.getRealValue("Nsb_ldp") - (everyPossibleThing.getRealValue("eff_sf_sb_ldp") * ( everyPossibleThing.getRealValue("sf_mc") * (everyPossibleThing.getRealValue("mu_ttwj_sb_ldp") + everyPossibleThing.getRealValue("mu_znn_sb_ldp") + everyPossibleThing.getRealValue("mu_ewo_sb_ldp") ) + everyPossibleThing.getRealValue("mu_susy_sb_ldp") ) ) , 0. ) * everyPossibleThing.getRealValue("sf_qcd_sb") * everyPossibleThing.getRealValue("Rlsb_passfail");
	     //
	     //
	     //double mu_ttwj_sig_sl = max(everyPossibleThing.getRealValue("Nsig_sl") - everyPossibleThing.getRealValue("eff_sf_sig_sl") * everyPossibleThing.getRealValue("mu_susy_sig_sl") , 0. );
	     //
	     //double mu_ttwj_sb_sl = max(everyPossibleThing.getRealValue("Nsb_sl") - everyPossibleThing.getRealValue("eff_sf_sb_sl") * everyPossibleThing.getRealValue("mu_susy_sb_sl") , 0. );
	     //
	     //double mu_ttwj_sb = max(everyPossibleThing.getRealValue("Nsb") - ( mu_qcd_sb + mu_Znn_sb + everyPossibleThing.getRealValue("eff_sf_sb")*( everyPossibleThing.getRealValue("mu_ewo_sb") + everyPossibleThing.getRealValue("mu_susy_sb")) ) , 0. );
	     //
	     //((RooRealVar*)everyPossibleThing.find("mu_znn_sig"    ))->setVal(mu_Znn_sig    );
	     //((RooRealVar*)everyPossibleThing.find("mu_znn_sb"     ))->setVal(mu_Znn_sb     );
	     //((RooRealVar*)everyPossibleThing.find("mu_qcd_sig"    ))->setVal(mu_qcd_sig    );
	     //((RooRealVar*)everyPossibleThing.find("mu_qcd_sb"     ))->setVal(mu_qcd_sb     );
	     //((RooRealVar*)everyPossibleThing.find("mu_ttwj_sig_sl"))->setVal(mu_ttwj_sig_sl);
	     //((RooRealVar*)everyPossibleThing.find("mu_ttwj_sb_sl" ))->setVal(mu_ttwj_sb_sl );
	     //((RooRealVar*)everyPossibleThing.find("mu_ttwj_sb"    ))->setVal(mu_ttwj_sb    );

	     //end of set up
	     //This still sometimes fails to find a good fit.  This means that the generated distributions will potentially look strange!

	     RooArgSet * allVariables = fPdf->getParameters(*fData);
	     RemoveConstantParameters(allVariables);
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
	     //allVariables->find("Rlsb_passfail_prim")->Print();
	     nuisancePlusPoi=*fPoints->get(fIndex++);
	     //allVariables->find("Rlsb_passfail_prim")->Print();

	     it = nuisancePlusPoi.iterator();
	     while((tmpParA = (RooRealVar*)it.Next())){
	       tmpParA->setConstant();
	       //cout << "Setting variable " << tmpParA->GetName() << " constant and to value " << tmpParA->getVal() << endl;
	     }

	     //RooAbsCollection* allVarSnap = allVariables->snapshot();
	     //RooArgSet * allNuisance = fPrior->getVariables();
	     //RemoveConstantParameters(allNuisance);
	     //RooArgSet nuisancePlusPoi(*allNuisance);


	     //it = nuisancePlusPoi.iterator();
	     //while((tmpParA = (RooRealVar*)it.Next())){
	     //  tmpParA->setConstant(kFALSE);
	     //}

	     RooAbsReal* nll = fPdf->createNLL(*fData);
	     //fPdf->findServer("pdf_Nsb_ldp")->Print("v");
	     //fPdf->getVariables()->Print("v");
	     //nll->getVariables()->Print("v");
	     //nll->getVariables()->find("mu_ttbarmc_sb_ldp")->Print("v");
	     //everyPossibleThing.Print("v");
	     int status =0;
	     MinNLL(nll,status,everyPossibleThing);

	     it.Reset();
	     while((tmpParA = (RooRealVar*)it.Next())){
	       tmpParA->setConstant(kFALSE);
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

	     if(nll) delete nll;
	     if(allVariables) delete allVariables;
	     //delete allVarSnap;
	     
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

  void MinNLL(RooAbsReal* NLL,int& status ,RooArgSet everyPossibleThing) {
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
   //minim.minimize(fMinimizer,"Scan");
   for (int tries = 0, maxtries = 3; tries <= maxtries; ++tries) {
      //	 status = minim.minimize(fMinimizer, ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
      status = minim.minimize(fMinimizer, "Minimize");
      double val = NLL->getVal();
      if(val != val || val == numeric_limits<double>::infinity() || val == -numeric_limits<double>::infinity() ) {
	status = -1;
	if(tries<=1) tries = 2;
	//everyPossibleThing.find("pdf_Nsb_ldp")->Print();
	//everyPossibleThing.find("Nsb_ldp")->Print();
	//everyPossibleThing.find("n_sb_ldp")->Print();
	//everyPossibleThing.find("mu_qcd_sb_ldp")->Print();
	//everyPossibleThing.find("mu_qcd_sb")->Print();
	//cout << "calculated value: " << max( everyPossibleThing.getRealValue("Nsb_ldp") - (everyPossibleThing.getRealValue("eff_sf_sb_ldp") * ( everyPossibleThing.getRealValue("sf_mc") * (everyPossibleThing.getRealValue("mu_ttwj_sb_ldp") + everyPossibleThing.getRealValue("mu_znn_sb_ldp") + everyPossibleThing.getRealValue("mu_ewo_sb_ldp") ) + everyPossibleThing.getRealValue("mu_susy_sb_ldp") ) ) , 0. ) * everyPossibleThing.getRealValue("sf_qcd_sb") * everyPossibleThing.getRealValue("Rlsb_passfail") << endl;
	//everyPossibleThing.find("mu_znn_sb")->Print();
	//everyPossibleThing.find("eff_sf_sb")->Print();
	//everyPossibleThing.find("mu_ewo_sb")->Print();
	//everyPossibleThing.find("mu_susy_sb")->Print();
	
      }
      if (status == 0) {  
         break;
      } else {
         if (tries > 1) {
            printf("    HybridToyMCSampler ----> Doing a re-scan first\n");
	    RooLinkedListIter jt = everyPossibleThing.iterator();
	    int whichBit = 0;
	    RooRealVar* tmpPar(NULL);
	    while((tmpPar = (RooRealVar*)jt.Next())){
	      whichBit++;
	      if (!tmpPar) {
		cout << "I can't find number " << whichBit << endl;
		return;
	      }
	       //cout << "Adding variable " << tmpParA->GetName() << " to nuisancePlusPoi" << endl;
	     }
            minim.minimize(fMinimizer,"Scan");
         }
         if (tries > 2) {
            printf("    HybridToyMCSampler ----> trying with strategy = 1\n");
	    RooLinkedListIter jt = everyPossibleThing.iterator();
	    int whichBit = 0;
	    RooRealVar* tmpPar(NULL);
	    while((tmpPar = (RooRealVar*)jt.Next())){
	      whichBit++;
	      if (!tmpPar) {
		cout << "I can't find number " << whichBit << endl;
		return;
	      }
	       //cout << "Adding variable " << tmpParA->GetName() << " to nuisancePlusPoi" << endl;
	     }
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

RooAbsData* HybridToyMCSampler::GenerateToyData(RooArgSet& paramPoint, double& weight, RooAbsPdf& pdf) const {
   // This method generates a toy data set for the given parameter point taking
   // global observables into account.
   // The values of the generated global observables remain in the pdf's variables.
   // They have to have those values for the subsequent evaluation of the
   // test statistics.

  //cout << "doing a generation step" << endl;

   if(!fObservables) {
      ooccoutE((TObject*)NULL,InputArguments) << "Observables not set." << endl;
      return NULL;
   }

   // assign input paramPoint
   RooArgSet* allVars = fPdf->getVariables();
   *allVars = paramPoint;


   // create nuisance parameter points
   if(!fNuisanceParametersSamplerAndFitter && fPriorNuisance && fNuisancePars)
      fNuisanceParametersSamplerAndFitter = new NuisanceParametersSamplerAndFitter(fPriorNuisance, fPdf, fFitToData, fParametersForTestStat, fSampledNuisances, fNToys, fExpectedNuisancePar);


   // generate global observables
   RooArgSet observables(*fObservables);
   if(fGlobalObservables  &&  fGlobalObservables->getSize()) {
      observables.remove(*fGlobalObservables);
      GenerateGlobalObservables(pdf);
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
      RooArgSet allVarsMinusParamPoint(*allVars);
      allVarsMinusParamPoint.remove(paramPoint, kFALSE, kTRUE); // match by name

      // get nuisance parameter point and weight
      fNuisanceParametersSamplerAndFitter->NextPoint(allVarsMinusParamPoint, weight);
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


} // end namespace RooStats
