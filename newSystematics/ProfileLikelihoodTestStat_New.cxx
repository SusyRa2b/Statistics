// @(#)root/roostats:$Id: ProfileLikelihoodTestStat_New.cxx 42339 2011-11-30 23:54:18Z moneta $
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
// Additional Contributions: Giovanni Petrucciani 
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ProfileLikelihoodTestStat_New.h"

Bool_t RooStats::ProfileLikelihoodTestStat_New::fgAlwaysReuseNll = kTRUE ;

Double_t RooStats::ProfileLikelihoodTestStat_New::EvaluateProfileLikelihood(int type, RooAbsData& data, RooArgSet& paramsOfInterest) {
        // interna function to evaluate test statistics
        // can do depending on type: 
        // type  = 0 standard evaluation, type = 1 find only unconditional NLL minimum, type = 2 conditional MLL

       if (!&data) {
	 cout << "problem with data" << endl;
	 return -1 ;
       }


       //data.Print("V");
       
       TStopwatch tsw; 
       tsw.Start();

       double initial_mu_value  = 0;
       RooRealVar* firstPOI = dynamic_cast<RooRealVar*>( paramsOfInterest.first());       
       if (firstPOI) initial_mu_value = firstPOI->getVal();
       //paramsOfInterest.getRealValue(firstPOI->GetName());

       RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();

       if (fPrintLevel < 3) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

       // simple
       Bool_t reuse=false;//(fReuseNll || fgAlwaysReuseNll) ;
       
       Bool_t created(kFALSE) ;
       if (!reuse || fNll==0) {
          RooArgSet* allParams = fPdf->getParameters(data);
          RooStats::RemoveConstantParameters(allParams);

	  RooArgSet everyPossibleThing;
	  fPdf->treeNodeServerList(&everyPossibleThing);
	  fPdf->branchNodeServerList(&everyPossibleThing);
	  fPdf->leafNodeServerList(&everyPossibleThing);
	  
	  RooLinkedListIter jt = everyPossibleThing.iterator();
	  int whichBit = 0;
	  RooRealVar* tmpPar(NULL);
	  while((tmpPar = (RooRealVar*)jt.Next())){
	    whichBit++;
	    if (!tmpPar) {
	      cout << "ProfileLikelihoodTestStat: I can't find variable number: " << whichBit << endl;
	      return -1;
	    }
	  }

          // need to call constrain for RooSimultaneous until stripDisconnected problem fixed
          fNll = (RooNLLVar*) fPdf->createNLL(data, RooFit::CloneData(kFALSE),RooFit::Constrain(*allParams));

          //	 fNll = (RooNLLVar*) fPdf->createNLL(data, RooFit::CloneData(kFALSE));
          //	 fProfile = (RooProfileLL*) fNll->createProfile(paramsOfInterest);
          created = kTRUE ;
          delete allParams;
          //cout << "creating profile LL " << fNll << " " << fProfile << " data = " << &data << endl ;
       }
       if (reuse && !created) {
          //cout << "reusing profile LL " << fNll << " new data = " << &data << endl ;
          fNll->setData(data,kFALSE) ;
       }

       // make sure we set the variables attached to this nll
       RooArgSet* attachedSet = fNll->getVariables();

       if(fFluctuatedNuisanceParameters && fDebuggingGenData){
	 RooLinkedListIter it = fFluctuatedNuisanceParameters->iterator();
	 RooRealVar* tmpPar = NULL, *tmpParA=NULL, *tmpParB=NULL;
	 while((tmpPar = (RooRealVar*)it.Next())){
	   tmpParA = dynamic_cast<RooRealVar*>( attachedSet->find(tmpPar->GetName()));
	   tmpParB = dynamic_cast<RooRealVar*>(((RooArgSet*)fDebuggingGenData->get())->find(tmpPar->GetName()));
	   if (tmpParA && tmpParB) tmpParA->setVal(0.0/*tmpParB->getVal()*/) , tmpParA->setConstant() ;
	 }
       }


       if (fPrintLevel > 0)
	 { 
	   RooLinkedListIter it = paramsOfInterest.iterator();
	   RooRealVar* tmpPar = NULL, *tmpParA=NULL;
	   while((tmpPar = (RooRealVar*)it.Next())){
             tmpParA =  dynamic_cast<RooRealVar*>( attachedSet->find(tmpPar->GetName()));
             if (tmpParA){
	       if( tmpPar->getVal() != tmpParA->getVal()){
	       cout << "POI value of " << tmpPar->GetName() << " is: " << tmpPar->getVal() << endl;
	       cout << "NOM value of " << tmpPar->GetName() << " is: " << tmpParA->getVal() << endl;
	       }
	     }
	   }
	 }
       
       *attachedSet = paramsOfInterest;

       RooArgSet* origAttachedSet = (RooArgSet*) attachedSet->snapshot();

       ///////////////////////////////////////////////////////////////////////
       // New profiling based on RooMinimizer (allows for Minuit2)
       // based on major speed increases seen by CMS for complex problems

	 //Fill some debugging information
       RooArgSet* trueSavedSet(NULL);
       RooArgSet* uncondSavedSet(NULL);
       RooArgSet* condSavedSet(NULL); 

       // other order
       // get the numerator
       RooArgSet* snap =  (RooArgSet*)paramsOfInterest.snapshot();

       tsw.Stop(); 
       double createTime = tsw.CpuTime();
       tsw.Start();

       // get the denominator
       double uncondML = 0;
       double fit_favored_mu = 0;
       int statusD = 0;
       if (type != 2) {
	 uncondML = GetMinNLL(statusD);
	 if (uncondML != uncondML || uncondML == numeric_limits<Double_t>::infinity() || uncondML == -numeric_limits<Double_t>::infinity() ) uncondML = GetMinNLL(statusD);
	 
	 // get best fit value for one-sided interval 
	 if (firstPOI) fit_favored_mu = attachedSet->getRealValue(firstPOI->GetName()) ;

	 //Fill some debugging information
	   
	 if(fDebuggingUncondData){
	   RooArgSet* saveSet= (RooArgSet*)fPdf->getVariables()->snapshot();
	   *saveSet=*data.get();
	   RooRealVar* poiPtr=NULL;
	   RooRealVar* poiFittedPtr=NULL;
	   if (firstPOI) poiPtr=new RooRealVar(*firstPOI,TString(firstPOI->GetName())+"_hypothesis") ;
	   if (firstPOI) poiPtr->setVal(initial_mu_value);
	   if (firstPOI) poiFittedPtr=new RooRealVar(*firstPOI,TString(firstPOI->GetName())+"_fitted") ;
	   if (firstPOI) poiFittedPtr->setVal(fit_favored_mu);
	   saveSet->addOwned(*poiPtr);
	   saveSet->addOwned(*poiFittedPtr);
	   //RemoveConstantParameters(saveSet);
	   //saveSet->Print("v");
	   //if(fDebuggingGenData) cout << "have the debug thing" << endl;
	   //fDebuggingGenData->Print("v");
	   RooArgSet* trueSet= (RooArgSet*)fDebuggingGenData->get()->snapshot();
	   ((RooRealVar*)trueSet->find(TString(firstPOI->GetName())+"_hypothesis"))->setVal(initial_mu_value);
	   ((RooRealVar*)trueSet->find(TString(firstPOI->GetName())+"_fitted"))->setVal(initial_mu_value);
	   double truePoiValue = ((RooRealVar*)fDebuggingGenData->get()->find(TString(firstPOI->GetName())))->getVal();
	   ((RooRealVar*)saveSet->find(firstPOI->GetName()))->setVal(truePoiValue);
	   //saveSet->find(TString(firstPOI->GetName())+"_fitted")->Print("v");
	   //fDebuggingGenData->get()->find(TString(firstPOI->GetName())+"_hypothesis")->Print();
	   //fDebuggingTrueData->add(*trueSet);
	   //fDebuggingUncondData->add(*saveSet);
	   trueSavedSet = trueSet;
	   uncondSavedSet = saveSet;
	   //cout << "------------------------------" << endl;
	   //cout << "input from function" << endl;
	   //saveSet->find(TString("Nsig"))->Print();
	   //cout << "input from debugging tool" << endl;
	   //trueSet->find(TString("Nsig"))->Print();
	   //cout << "------------------------------" << endl;
	   //delete trueSet;
	   //delete saveSet;
	 }

       }
       tsw.Stop();
       double fitTime1  = tsw.CpuTime();
          
       double ret = 0; 
       int statusN = 0;
       tsw.Start();

       double condML = 0; 

       bool doConditionalFit = (type != 1); 

       // skip the conditional ML (the numerator) only when fit value is smaller than test value
       if (fOneSided &&  fit_favored_mu > initial_mu_value) { 
          doConditionalFit = false; 
          condML = uncondML;
       }

       if (doConditionalFit) {  


          //       cout <<" reestablish snapshot"<<endl;
          *attachedSet = *origAttachedSet;

 
          // set the POI to constant
          RooLinkedListIter it = paramsOfInterest.iterator();
          RooRealVar* tmpPar = NULL, *tmpParA=NULL;
          while((tmpPar = (RooRealVar*)it.Next())){
             tmpParA =  dynamic_cast<RooRealVar*>( attachedSet->find(tmpPar->GetName()));
             if (tmpParA) tmpParA->setConstant();
	  }


          // check if there are non-const parameters so it is worth to do the minimization
          RooArgSet allParams(*attachedSet); 
          RooStats::RemoveConstantParameters(&allParams);
          
          // in case no nuisance parameters are present
          // no need to minimize just evaluate the nll
          if (allParams.getSize() == 0 ) {
             condML = fNll->getVal(); 
          }
          else {              
	    condML = GetMinNLL(statusN);
	    if(condML != condML|| condML == numeric_limits<Double_t>::infinity() || condML == -numeric_limits<Double_t>::infinity() ) condML = GetMinNLL(statusN);
          }
       }

       //Fill some debugging information
       if(fDebuggingCondData){
	 RooArgSet* saveSet= (RooArgSet*)fPdf->getVariables()->snapshot();
	 *saveSet=*data.get();
	 RooRealVar* poiPtr=NULL;
	 RooRealVar* poiFittedPtr=NULL;
	 if (firstPOI) poiPtr=new RooRealVar(*firstPOI,TString(firstPOI->GetName())+"_hypothesis") ;
	 if (firstPOI) poiPtr->setVal(initial_mu_value);
	 if (firstPOI) poiFittedPtr=new RooRealVar(*firstPOI,TString(firstPOI->GetName())+"_fitted") ;
	 if (firstPOI) poiFittedPtr->setVal(0);
	 saveSet->addOwned(*poiPtr);
	 saveSet->addOwned(*poiFittedPtr);
	 //((RooRealVar*)fDebuggingGenData->get()->find(TString(firstPOI->GetName())+"_hypothesis"))->setVal(initial_mu_value);
	 double truePoiValue = ((RooRealVar*)fDebuggingGenData->get()->find(TString(firstPOI->GetName())))->getVal();
	 ((RooRealVar*)saveSet->find(firstPOI->GetName()))->setVal(truePoiValue);
	 condSavedSet = saveSet;
	 //fDebuggingCondData->add(*saveSet);
	 //delete saveSet;
       }

       //case where profll < 0 must be a numeric effect rather than a real effect 
       //re-run in this case to correct any errors:

       if(type != 2 && condML<uncondML)
	 {
	   if (fPrintLevel > 0) { 
	     cout << "fixing negative profile likelihood value" << endl;
	     cout << "current value : " << condML-uncondML << endl;
	   }
	   RooLinkedListIter it = paramsOfInterest.iterator();
	   RooRealVar* tmpPar = NULL, *tmpParA=NULL;
	   while((tmpPar = (RooRealVar*)it.Next())){
             tmpParA =  dynamic_cast<RooRealVar*>( attachedSet->find(tmpPar->GetName()));
             if (tmpParA) tmpParA->setConstant(kFALSE);
	   }
	   uncondML = GetMinNLL(statusD);
	   if (uncondML != uncondML || uncondML == numeric_limits<Double_t>::infinity() || uncondML == -numeric_limits<Double_t>::infinity() ) uncondML = GetMinNLL(statusD);
	   if (firstPOI) fit_favored_mu = attachedSet->getRealValue(firstPOI->GetName()) ;
	   RooArgSet allParams(*attachedSet); 
	   RooStats::RemoveConstantParameters(&allParams);
	   if ((fOneSided && fit_favored_mu > initial_mu_value) || allParams.getSize() == 0 ) { 
	     condML = uncondML;
	   }
	   else
	     {
	       it.Reset();
	       while((tmpPar = (RooRealVar*)it.Next())){
		 tmpParA =  dynamic_cast<RooRealVar*>( snap->find(tmpPar->GetName()));
		 double val = 0;
		 if (tmpParA) val = tmpParA->getVal();
		 tmpParA =  dynamic_cast<RooRealVar*>( attachedSet->find(tmpPar->GetName()));
		 if (tmpParA) tmpParA->setVal(val);
		 if (tmpParA) tmpParA->setConstant();
	       }       
	       condML = GetMinNLL(statusN);
	       if(condML != condML || condML == numeric_limits<Double_t>::infinity() || condML == -numeric_limits<Double_t>::infinity() ) condML = GetMinNLL(statusN);
	     }
	   if (fPrintLevel > 0) { 
	     cout << "new value     : " << condML-uncondML << endl ;
	   }
	   if(condML<uncondML) condML = uncondML;
	 }

       tsw.Stop();
       double fitTime2 = tsw.CpuTime();

       if (fPrintLevel > 0) { 
          std::cout << "EvaluateProfileLikelihood - ";
          if (type <= 1)  
             std::cout << "mu hat = " << fit_favored_mu  <<  " uncond ML = " << uncondML; 
          if (type != 1) 
             std::cout << " cond ML = " << condML;
          if (type == 0)
             std::cout << " pll =  " << condML-uncondML; 
          std::cout << " time (create/fit1/2) " << createTime << " , " << fitTime1 << " , " << fitTime2  
                    << std::endl;
       }

       //(fDoProfllCheck)
       //{
       //  double profileLLCheck = 0;
       //  if(type == 0)
       //    {
       //      if(fOneSided && fit_favored_mu > initial_mu_value) profileLLCheck = 0;
       //      else {
       //	 *attachedSet = *origAttachedSet;
       //	 RooMinuit minim(*fNll);
       //	 minim.setPrintLevel(-999) ;
       //	 minim.zeroEvalCount() ;
       //	 minim.migrad();
       //	 double firstProfileLLCheck = fNll->getVal() ;
       //	 *attachedSet = *origAttachedSet;
       //	 RooLinkedListIter it = paramsOfInterest.iterator();
       //	 RooRealVar* tmpPar = NULL, *tmpParA=NULL;
       //	 while((tmpPar = (RooRealVar*)it.Next())){
       //	   tmpParA =  dynamic_cast<RooRealVar*>( attachedSet->find(tmpPar->GetName()));
       //	   if (tmpParA) tmpParA->setConstant();
       //	 }
       //	 minim.migrad();
       //	 double secondProfileLLCheck = fNll->getVal();
       //	 if (fPrintLevel > 0) { 
       //	   cout << "Looking at variable: " << firstPOI->GetName() << endl;
       //	   cout << "The value of the variable is : " << initial_mu_value << endl ;
       //	   cout << "The best fit from Minuit is: " << firstProfileLLCheck << endl ;
       //	   cout << "The best constrained fit from Minuit is: " << fNll->getVal() << endl ;
       //	   cout << "The profile likelihood is :" << secondProfileLLCheck << endl;
       //	 }
       //	 *attachedSet = *origAttachedSet;
       //	 fProfile = (RooProfileLL*) fNll->createProfile(paramsOfInterest) ;
       //	 profileLLCheck = fProfile->getVal();
       //	 if(fPrintLevel > 0)
       //	   {
       //	     cout << "uncond by hand is : " << firstProfileLLCheck << endl;
       //	     cout << "cond by hand is   : " << secondProfileLLCheck << endl;
       //	     cout << "The by hand is    : " << secondProfileLLCheck - firstProfileLLCheck << endl;
       //	     cout << "The RooFit Tool is: "<< profileLLCheck << endl ;
       //	   }
       //	 if(profileLLCheck < 0)
       //	   {
       //	     it.Reset();
       //	     while((tmpPar = (RooRealVar*)it.Next())){
       //	       tmpParA =  dynamic_cast<RooRealVar*>( attachedSet->find(tmpPar->GetName()));
       //	       if (tmpParA) tmpParA->setConstant();
       //	     }
       //	     minim.migrad();
       //	     it.Reset();
       //	     while((tmpPar = (RooRealVar*)it.Next())){
       //	       tmpParA =  dynamic_cast<RooRealVar*>( attachedSet->find(tmpPar->GetName()));
       //	       if (tmpParA) tmpParA->setConstant(kFALSE);
       //	     }		     
       //	     fProfile->validateAbsMin() ;
       //	     profileLLCheck = fProfile->getVal();
       //	     if(profileLLCheck < 0) profileLLCheck = 0;
       //	   }
       //	 delete fProfile;
       //	 if (fPrintLevel > 0) { 
       //	   cout << "The best fit from RooFit Tool is: " << profileLLCheck << endl ;
       //	   if( profileLLCheck - (condML-uncondML)  > 1E-3 || profileLLCheck - (condML-uncondML) < -1E-3 )
       //	     {
       //	       cout << "difference in profile Likelihood Methods!" << endl;
       //	       cout << "Standard Value is: " << condML-uncondML << endl;
       //	       cout << "Check Value is:    " << profileLLCheck << endl;
       //	       cout << "Standard/Check is: " << (condML-uncondML)/profileLLCheck << endl;
       //	       cout << "condML is   : " << condML << endl ; 
       //	       cout << "uncondML is : " << uncondML << endl ; 
       //	     }
       //	   else 
       //	     {
       //	       cout << "no difference in profile Likelihood Methods!" << endl;
       //	       cout << "Standard Value is: " << condML-uncondML << endl;
       //	       cout << "Check Value is:    " << profileLLCheck << endl;
       //	       cout << "Standard/Check is: " << (condML-uncondML)/profileLLCheck << endl;
       //	       cout << "condML is   : " << condML << endl ; 
       //	       cout << "uncondML is : " << uncondML << endl ; 
       //	     }
       //	 }
       //	 condML = profileLLCheck ;
       //	 uncondML = 0 ;
       //      }
       //    }
       //}

       // need to restore the values ?
       *attachedSet = *origAttachedSet;

       if(condML-uncondML != condML-uncondML || statusN!=0 || statusD!=0)
	 {
	   cout << "uncondML: " << uncondML << endl;
	   cout << "condML: " << condML << endl;
	   cout << "fitted mu " << fit_favored_mu << endl;
	   data.get()->Print("v");
	   if (firstPOI) attachedSet->find(firstPOI->GetName())->Print() ;
	   
	   if( trueSavedSet   ) delete trueSavedSet   ;
	   if( uncondSavedSet ) delete uncondSavedSet ;
	   if( condSavedSet   ) delete condSavedSet   ;
	 }
       else 
	 {
	   //trueSavedSet->find("Nsig")->Print();
	   //uncondSavedSet->find("Nsig")->Print();
	   //condSavedSet->find("Nsig")->Print();
	   //trueSavedSet->find("mu_susy_sig")->Print();
	   //uncondSavedSet->find("mu_susy_sig")->Print();
	   //condSavedSet->find("mu_susy_sig")->Print();
	   //trueSavedSet->find("mu_susy_sig_fitted")->Print();
	   //uncondSavedSet->find("mu_susy_sig_fitted")->Print();
	   //condSavedSet->find("mu_susy_sig_fitted")->Print();
	   //trueSavedSet->find("mu_ttwj_sb")->Print();
	   //uncondSavedSet->find("mu_ttwj_sb")->Print();
	   //condSavedSet->find("mu_ttwj_sb")->Print();

	   if( trueSavedSet   ) fDebuggingTrueData->add(  *trueSavedSet  )  ;
	   if( uncondSavedSet ) fDebuggingUncondData->add(*uncondSavedSet) ;
	   if( condSavedSet   ) fDebuggingCondData->add(  *condSavedSet  ) ;

	   if( trueSavedSet   ) delete trueSavedSet   ;
	   if( uncondSavedSet ) delete uncondSavedSet ;
	   if( condSavedSet   ) delete condSavedSet   ;
	 }

       if(fFluctuatedNuisanceParameters && fDebuggingGenData){
	 RooLinkedListIter it = fFluctuatedNuisanceParameters->iterator();
	 RooRealVar* tmpPar = NULL, *tmpParA=NULL;
	 while((tmpPar = (RooRealVar*)it.Next())){
	   tmpParA = dynamic_cast<RooRealVar*>( attachedSet->find(tmpPar->GetName()));
	   if (tmpParA) tmpParA->setConstant(kFALSE) ;
	 }
       }

       if(attachedSet) delete attachedSet;
       if(origAttachedSet) delete origAttachedSet;
       if(snap) delete snap;

       if (!reuse) {
	 delete fNll;
	 fNll = 0; 
	 // delete fProfile;
	 fProfile = 0 ;
       }

       RooMsgService::instance().setGlobalKillBelow(msglevel);

       if(statusN!=0 || statusD!=0 || condML-uncondML != condML-uncondML)
	 ret= -1; // indicate failed fit

       if (type == 1) return uncondML;
       if (type == 2) return condML;
       return condML-uncondML;
             
     }     

double RooStats::ProfileLikelihoodTestStat_New::GetMinNLL(int& status) {
   //find minimum of NLL using RooMinimizer

   RooMinimizer minim(*fNll);
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
      double val = fNll->getVal();
      if(val != val || val == numeric_limits<double>::infinity() || val == -numeric_limits<double>::infinity() ) status = -1;
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

  //RooMinuit minim(*fNll);
  //minim.setPrintLevel(-999) ;
  //minim.zeroEvalCount() ;
  //minim.migrad();
  double val =  fNll->getVal();
  //minim.optimizeConst(false); 
  
  return val;
}
