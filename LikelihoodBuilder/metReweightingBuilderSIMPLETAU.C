#include <iostream>
#include <fstream>
#include <string.h>
#include <complex>
#include <map>

#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TH2F.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TStringLong.h"
#include "TString.h"
#include "TPRegexp.h"

#include "RooArgSet.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooTrace.h"
#include "RooUniform.h"
#include "RooAddition.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooRatio.h"
#include "RooAddition.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooPlot.h"

#include "RooCorrelatedBetaGeneratorHelper.h"
#include "RooCorrelatedBetaPrimeGeneratorHelper.h"
#include "betaHelperFunctions.h"
#include "RooBetaInverseCDF.h"
#include "RooBetaPrimeInverseCDF.h"
#include "rooFitBetaHelperFunctions.h"
#include "RooNormalFromFlatPdf.h"


#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "TMath.h"

using namespace RooFit ;
using namespace RooStats ;



//////////////////////////////////////////////////////////////////////////////////////////////////////

void makePrediction( RooWorkspace& wspace, TString thisBin ){


      ////////////////////////////////////
      ////////////////////////////////////
      // PUTTING TOGETHER THE 0L PREDICTION FOR THIS BIN

      TString zeroLeptonName("zeroLepton_");
      zeroLeptonName+=thisBin;
      //////
      //cout << " BEFORE GETTING PREDICTION COMPONENTS " << endl;

      TString PolarizationName1(zeroLeptonName);
      PolarizationName1.Append("_Theta1_TopWJetsYield");
    
      TString PolarizationName2(zeroLeptonName);
      PolarizationName2.Append("_Theta2_TopWJetsYield");
        
      TString PolarizationName3(zeroLeptonName);
      PolarizationName3.Append("_Theta3_TopWJetsYield");
        
      TString PolarizationName4(zeroLeptonName);
      PolarizationName4.Append("_Theta4_TopWJetsYield");
        
      TString PolarizationName5(zeroLeptonName);
      PolarizationName5.Append("_Theta5_TopWJetsYield");
    
      //    cout << " GOT POLARIZATION PREDICTIONS" << endl;
      //////
      TString TauHadName1(zeroLeptonName);
      TauHadName1.Append("_1Tau_TopWJetsYield");
      
      TString TauHadName2(zeroLeptonName);
      TauHadName2.Append("_2Tau_TopWJetsYield");
      
      TString TauHadName3(zeroLeptonName);
      TauHadName3.Append("_MuTau_TopWJetsYield");
      
      TString TauHadName4(zeroLeptonName);
      TauHadName4.Append("_ETau_TopWJetsYield");
      
      //    cout << " GOT TAU->HAD PREDICTIONS" << endl;
      //////
      TString DilepName1(zeroLeptonName);
      DilepName1.Append("_MuMu_TopWJetsYield");
      
      TString DilepName2(zeroLeptonName);
      DilepName2.Append("_EMu_TopWJetsYield");
      
      TString DilepName3(zeroLeptonName);
      DilepName3.Append("_EE_TopWJetsYield");
      
      //    cout << " GOT DILEPTON PREDICTIONS" << endl;
      //////
      
      //cout << " GOT ALL COMPONENTS OF 0L PREDICTION, ADDING THEM... " << endl;
      
      RooAddition zeroLeptonTopWJetsPolarizationYield(zeroLeptonName+"_TopWJetsPolarizationYield",zeroLeptonName+"_TopWJetsPolarizationYield",
						      RooArgSet( *wspace.arg(PolarizationName1.Data()),
								 *wspace.arg(PolarizationName2.Data()),
								 *wspace.arg(PolarizationName3.Data()),
								 *wspace.arg(PolarizationName4.Data()),
								 *wspace.arg(PolarizationName5.Data()) ));
      
      //cout << " after adding the pol " << endl;

      RooAddition zeroLeptonTopWJetsTauHadYield(zeroLeptonName+"_TopWJetsTauHadYield",zeroLeptonName+"_TopWJetsTauHadYield",
						RooArgSet( *wspace.arg(TauHadName1.Data()),
							   *wspace.arg(TauHadName2.Data()),
							   *wspace.arg(TauHadName3.Data()),
							   *wspace.arg(TauHadName4.Data()) ));
    
      //cout << " after adding the tau " << endl;

      RooAddition zeroLeptonTopWJetsDileptonYield(zeroLeptonName+"_TopWJetsDileptonYield",zeroLeptonName+"_TopWJetsDileptonYield",
						  RooArgSet( *wspace.arg(DilepName1.Data()),
							     *wspace.arg(DilepName2.Data()),
							     *wspace.arg(DilepName3.Data()) )); 

      /*
      cout << " pol pieces: ";
      (*wspace.function(PolarizationName2.Data())).Print();
      (*wspace.function(PolarizationName3.Data())).Print();
      (*wspace.function(PolarizationName4.Data())).Print();
      (*wspace.function(PolarizationName5.Data())).Print();
      cout << endl;
      
      cout << " tau pieces: ";
      (*wspace.arg(TauHadName1.Data())).Print();
      (*wspace.arg(TauHadName2.Data())).Print();
      (*wspace.arg(TauHadName3.Data())).Print();
      (*wspace.arg(TauHadName4.Data())).Print();
      cout << endl;

      cout << " dilep pieces: ";
      (*wspace.arg(DilepName1.Data())).Print(); 
      (*wspace.arg(DilepName2.Data())).Print(); 
      (*wspace.arg(DilepName3.Data())).Print(); 
      cout << endl;

      cout << " BACKGROUND PREDICTION PIECES:  " << zeroLeptonTopWJetsPolarizationYield.getValV() << "  " 
	   << zeroLeptonTopWJetsTauHadYield.getValV() << "  " << zeroLeptonTopWJetsDileptonYield.getValV() << endl;
      */

      /////////////////////////////////////////////////////////////////////////////////////
      // GET AND APPLY SIGNAL FRACTION FOR THIS 0L BIN, ADD TO ALL PREDICTIONS
      
      //cout << " GET 0L COUNTS AND SIG FRAC " << endl;


      TString zeroLeptonCountName(zeroLeptonName);
      zeroLeptonCountName.Append("_Count");

      TString zeroLeptonSignalYieldName(zeroLeptonName);
      zeroLeptonSignalYieldName.Append("_SignalYield");
      


      RooAddition zeroLeptonYieldSum(zeroLeptonName+"_YieldSum",zeroLeptonName+"_YieldSum",
				     RooArgSet( *wspace.arg(zeroLeptonSignalYieldName.Data()), 
						zeroLeptonTopWJetsPolarizationYield,
						zeroLeptonTopWJetsTauHadYield,
						zeroLeptonTopWJetsDileptonYield 
						));
      
      /*
      cout << " SIGNAL YIELD:  " << zeroLeptonSignalYieldName.Data() << "  " ;
      (*wspace.arg(zeroLeptonSignalYieldName.Data())).Print();

      cout << " COUNTS:  " << zeroLeptonCountName.Data() << "  "
	   << (*wspace.var(zeroLeptonCountName.Data())).getValV() << "  " << endl;

      cout << " TOTAL PREDICTED IN BIN:  ";
      zeroLeptonYieldSum.Print();
      */

      /////////////////////////////////////////////////////////////////////////////////////
      // SET CONSTRAINTS FOR THIS 0L BIN
      /****** SHOULD LEAVE THIS OUT IN COMPLETE LIKELIHOOD!!! ******/
      //cout << " SETTING 0L CONSTRAINTS " << endl;
      
      RooPoisson zeroLeptonConstraint(zeroLeptonName+"_Constraint",zeroLeptonName+"_Constraint",
				    *wspace.var(zeroLeptonCountName.Data()),zeroLeptonYieldSum);
      
      wspace.import( zeroLeptonConstraint,RecycleConflictNodes() );

      /*
      cout << endl << endl;
      zeroLeptonConstraint.Print();
      cout << endl << endl;
      */
      
}// END OF PREDICTIONS


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////



void makePolarizationConstraintsPredictions( RooWorkspace& wspace, TString binname ){

  //////////////////////////////////////////////////////////////////
  // LOOP OVER STUFF IN DTHETA BINS!!!
  //////////////////////////////////////////////////////////////////

  for( int dtheta=1; dtheta<6; dtheta++ ){

    //cout << " in dtheta bin " << endl;
    TString ending1("_Theta");
    ending1+=dtheta;

    TString zeroLeptonName("zeroLepton_");
    zeroLeptonName+=binname;
    zeroLeptonName.Append(ending1);
    TString oneLeptonName("oneLepton_");
    oneLeptonName+=binname;
    oneLeptonName.Append(ending1);
 
    TString oneTightMuName("oneTightMu_");
    oneTightMuName+=binname;
    oneTightMuName.Append(ending1);
    TString oneLooseMuName("oneLooseMu_");
    oneLooseMuName+=binname;
    oneLooseMuName.Append(ending1);
    TString oneLooseEName("oneLooseE_");
    oneLooseEName+=binname;
    oneLooseEName.Append(ending1);

    ////////////////////////////////////////////////////
    // GET CONTROL SAMPLE COUNTS FOR THIS DTHETA BIN 
    
    //TString ending("_Theta");
    //ending+=dtheta;
    //ending.Append("_Count");
    TString ending("_Count");
    
    //   cout << " getting counts for pol method " << endl;

    TString oneTightMuCountName(oneTightMuName);
    oneTightMuCountName.Append(ending);
    TString oneLooseMuCountName(oneLooseMuName);
    oneLooseMuCountName.Append(ending);
    TString oneLooseECountName(oneLooseEName);
    oneLooseECountName.Append(ending);

    //    cout << oneTightMuCountName << "  " << oneLooseMuCountName << "  " << oneLooseECountName << endl;

    RooRealVar oneTightMuCount = *wspace.var(oneTightMuCountName.Data());
    RooRealVar oneLooseMuCount = *wspace.var(oneLooseMuCountName.Data());
    RooRealVar oneLooseECount = *wspace.var(oneLooseECountName.Data());

    ////////////////////////////////////////////////////
    // GET SIGNAL YIELDS FOR THIS DTHETA BIN 
 
    //    cout << " getting signal yields for pol method " << endl;

    TString oneTightMuSignalYieldName(oneTightMuName);
    oneTightMuSignalYieldName.Append("_SignalYield");
    TString oneLooseMuSignalYieldName(oneLooseMuName);
    oneLooseMuSignalYieldName.Append("_SignalYield");
    TString oneLooseESignalYieldName(oneLooseEName);
    oneLooseESignalYieldName.Append("_SignalYield");

    //    cout << oneTightMuSignalYieldName << "  " << oneLooseMuSignalYieldName << "  " << oneLooseESignalYieldName << endl;
    /*
    RooAbsArg oneTightMuSignalYield = *wspace.arg(oneTightMuSignalYieldName.Data());
    RooAbsArg oneLooseMuSignalYield = *wspace.arg(oneLooseMuSignalYieldName.Data());
    RooAbsArg oneLooseESignalYield = *wspace.arg(oneLooseESignalYieldName.Data());
    */
    //    cout << " getting scale factor for pol method " << endl;

    ////////////////////////////////////////////////////
    // GET THE 1L->0L SCALE FACTOR FOR THIS DTHETA BIN 
    
    TString oneLeptonScaleFactorName(oneLeptonName);
    oneLeptonScaleFactorName.Append("_ScaleFactor");

    ////////////////////////////////////////////////////
    // INITIALIZE TTBAR COMPONENT OF THIS DTHETA BIN TO START WITH OBSERVED CONTENT
    // EACH YIELD IS A SEPARATE NUISANCE PARAMETER

    RooRealVar oneTightMuTopWJetsYield(oneTightMuName+"_TopWJetsYield",oneTightMuName+"_TopWJetsYield",
				       oneTightMuCount.getVal(),0.000001,20000.);
    wspace.import(oneTightMuTopWJetsYield,RecycleConflictNodes());
    wspace.extendSet("nuisances",oneTightMuTopWJetsYield.GetName());

    RooRealVar oneLooseMuTopWJetsYield(oneLooseMuName+"_TopWJetsYield",oneLooseMuName+"_TopWJetsYield",
				       oneLooseMuCount.getVal(),0.000001,20000.);
    wspace.import(oneLooseMuTopWJetsYield,RecycleConflictNodes());
    wspace.extendSet("nuisances",oneLooseMuTopWJetsYield.GetName());

    RooRealVar oneLooseETopWJetsYield(oneLooseEName+"_TopWJetsYield",oneLooseEName+"_TopWJetsYield",
				      oneLooseECount.getVal(),0.000001,20000.);
    wspace.import(oneLooseETopWJetsYield,RecycleConflictNodes());
    wspace.extendSet("nuisances",oneLooseETopWJetsYield.GetName());


    ////////////////////////////////////////////////////
    // COMBINE SIGNAL AND TTWT COMPONENTS OF THIS DTHETA BIN AND SET CONTROL SAMPLE CONSTRAINTS
    
    /****** NOTE: INSERT SINGLE LEPTON EFFICIENCIES HERE??? ******/

    RooAddition oneTightMuYieldSum(oneTightMuName+"_YieldSum",oneTightMuName+"_YieldSum",
				   RooArgSet(*wspace.arg(oneTightMuSignalYieldName.Data()),oneTightMuTopWJetsYield));
    RooAddition oneLooseMuYieldSum(oneLooseMuName+"_YieldSum",oneLooseMuName+"_YieldSum",
				   RooArgSet(*wspace.arg(oneLooseMuSignalYieldName.Data()),oneLooseMuTopWJetsYield));
    RooAddition oneLooseEYieldSum(oneLooseEName+"_YieldSum",oneLooseEName+"_YieldSum",
				  RooArgSet(*wspace.arg(oneLooseESignalYieldName.Data()),oneLooseETopWJetsYield));
   

    RooPoisson oneTightMuConstraint(oneTightMuName+"_Constraint",oneTightMuName+"_Constraint",oneTightMuCount,oneTightMuYieldSum);
    wspace.import( oneTightMuConstraint,RecycleConflictNodes() );
    RooPoisson oneLooseMuConstraint(oneLooseMuName+"_Constraint",oneLooseMuName+"_Constraint",oneLooseMuCount,oneLooseMuYieldSum);
    wspace.import( oneLooseMuConstraint,RecycleConflictNodes() );
    RooPoisson oneLooseEConstraint(oneLooseEName+"_Constraint",oneLooseEName+"_Constraint",oneLooseECount,oneLooseEYieldSum);
    wspace.import( oneLooseEConstraint,RecycleConflictNodes() );


    ////////////////////////////////////////////////////
    // CONSTRUCT THE 0L BACKGROUND PREDICTION FROM THIS DTHETA BIN
    // ****** DONT FORGET ABOUT TRIGGER INEFFICIENCIES! ****** //
 
    RooAddition oneLeptonTopWJetsYield(oneLeptonName+"_TopWJetsYield",oneLeptonName+"_TopWJetsYield",
				       RooArgSet(oneTightMuTopWJetsYield,oneLooseMuTopWJetsYield,oneLooseETopWJetsYield));
    RooProduct  zeroLeptonTopWJetsYield(zeroLeptonName+"_TopWJetsYield",zeroLeptonName+"_TopWJetsYield",
					RooArgSet(*wspace.arg(oneLeptonScaleFactorName.Data()),oneLeptonTopWJetsYield));
    wspace.import(zeroLeptonTopWJetsYield,RecycleConflictNodes());

  }// END OF LOOP OVER DTHETA BINS
  
}// END OF BIG BIN IN NBS/HT/MET


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

void makeDileptonConstraintsPredictions( RooWorkspace& wspace, TString binname,  TString binname_outside ){


  TString twoTightMuName("twoTightMu_");
  twoTightMuName+=binname;
  TString twoTightMuLooseMuName("twoTightMuLooseMu_");
  twoTightMuLooseMuName+=binname;
  TString twoLooseMuName("twoLooseMu_");
  twoLooseMuName+=binname;
  TString twoTightMuLooseEName("twoTightMuLooseE_");
  twoTightMuLooseEName+=binname;
  TString twoLooseMuLooseEName("twoLooseMuLooseE_");
  twoLooseMuLooseEName+=binname;
  TString twoLooseEName("twoLooseE_");
  twoLooseEName+=binname;
   
  /////////////////////////////////////////////////////////////////////////////////////
  // GET DILEPTON COUNTS
  
  TString twoTightMuCountName(twoTightMuName);
  twoTightMuCountName.Append("_Count");
  TString twoTightMuLooseMuCountName(twoTightMuLooseMuName);
  twoTightMuLooseMuCountName.Append("_Count");
  TString twoLooseMuCountName(twoLooseMuName);
  twoLooseMuCountName.Append("_Count");
  TString twoTightMuLooseECountName(twoTightMuLooseEName);
  twoTightMuLooseECountName.Append("_Count");
  TString twoLooseMuLooseECountName(twoLooseMuLooseEName);
  twoLooseMuLooseECountName.Append("_Count");
  TString twoLooseECountName(twoLooseEName);
  twoLooseECountName.Append("_Count");

  RooRealVar twoTightMuCount = *wspace.var(twoTightMuCountName.Data());
  RooRealVar twoTightMuLooseMuCount = *wspace.var(twoTightMuLooseMuCountName.Data());
  RooRealVar twoLooseMuCount = *wspace.var(twoLooseMuCountName.Data());
  RooRealVar twoTightMuLooseECount = *wspace.var(twoTightMuLooseECountName.Data());
  RooRealVar twoLooseMuLooseECount = *wspace.var(twoLooseMuLooseECountName.Data());
  RooRealVar twoLooseECount = *wspace.var(twoLooseECountName.Data());

  //  cout << " END OF GETTING DILEPTON COUNTS " << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // GET DILEPTON SCALE FACTORS
  // SHARED ACROSS HT/MET BINS 

  TString MuMuScaleFactorName("MuMu_");
  MuMuScaleFactorName+=binname_outside;
  MuMuScaleFactorName.Append("_ScaleFactor");
  TString EMuScaleFactorName("EMu_");
  EMuScaleFactorName+=binname_outside;
  EMuScaleFactorName.Append("_ScaleFactor");
  TString EEScaleFactorName("EE_");
  EEScaleFactorName+=binname_outside;
  EMuScaleFactorName.Append("_ScaleFactor");


  /////////////////////////////////////////////////////////////////////////////////////
  // GET DILEPTON SIGNAL YIELDS

  TString twoTightMuSignalYieldName(twoTightMuName);
  twoTightMuSignalYieldName.Append("_SignalYield");
  TString twoTightMuLooseMuSignalYieldName(twoTightMuLooseMuName);
  twoTightMuLooseMuSignalYieldName.Append("_SignalYield");
  TString twoLooseMuSignalYieldName(twoLooseMuName);
  twoLooseMuSignalYieldName.Append("_SignalYield");
  TString twoTightMuLooseESignalYieldName(twoTightMuLooseEName);
  twoTightMuLooseESignalYieldName.Append("_SignalYield");
  TString twoLooseMuLooseESignalYieldName(twoLooseMuLooseEName);
  twoLooseMuLooseESignalYieldName.Append("_SignalYield");
  TString twoLooseESignalYieldName(twoLooseEName);
  twoLooseESignalYieldName.Append("_SignalYield");


  //  cout << " END OF GETTING DILEPTON SIGNAL YIELDS " << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // SET BACKGROUND CONTRIBUTIONS TO DILEPTON SAMPLES 
  // IMPORT BACKGROUND CONTRIBUTIONS TO DILEPTON SAMPLES TO WORKSPACE AS NUISANCE PARAMS

  // ****** TRIGGER INEFFICIENCIES MUST BE PUT IN SOMEWHERE FOR BACKGROUND AND SIGNAL ******//
  
  RooRealVar twoTightMuTopWJetsYield(twoTightMuName+"_TopWJetsYield",twoTightMuName+"_TopWJetsYield",
				     twoTightMuCount.getVal(),0.,1000.);
  wspace.import(twoTightMuTopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoTightMuTopWJetsYield.GetName());

  RooRealVar twoTightMuLooseMuTopWJetsYield(twoTightMuLooseMuName+"_TopWJetsYield",twoTightMuLooseMuName+"_TopWJetsYield",
					    twoTightMuLooseMuCount.getVal(),0.,1000.);
  wspace.import(twoTightMuLooseMuTopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoTightMuLooseMuTopWJetsYield.GetName());

  RooRealVar twoLooseMuTopWJetsYield(twoLooseMuName+"_TopWJetsYield",twoLooseMuName+"_TopWJetsYield",
				     twoLooseMuCount.getVal(),0.,1000.);
  wspace.import(twoLooseMuTopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoLooseMuTopWJetsYield.GetName());

  RooRealVar twoTightMuLooseETopWJetsYield(twoTightMuLooseEName+"_TopWJetsYield",twoTightMuLooseEName+"_TopWJetsYield",
					   twoTightMuLooseECount.getVal(),0.,1000.);
  wspace.import(twoTightMuLooseETopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoTightMuLooseETopWJetsYield.GetName());

  RooRealVar twoLooseMuLooseETopWJetsYield(twoLooseMuLooseEName+"_TopWJetsYield",twoLooseMuLooseEName+"_TopWJetsYield",
					   twoLooseMuLooseECount.getVal(),0.,1000.);
  wspace.import(twoLooseMuLooseETopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoLooseMuLooseETopWJetsYield.GetName());

  RooRealVar twoLooseETopWJetsYield(twoLooseEName+"_TopWJetsYield",twoLooseEName+"_TopWJetsYield",
				    twoLooseECount.getVal(),0.,1000.);
  wspace.import(twoLooseETopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoLooseETopWJetsYield.GetName());


  /////////////////////////////////////////////////////////////////////////////////////
  // PRODUCE SIGNAL+BACKGROUND DILEPTON CONSTRAINTS (NOT EVER BINNED IN DTHETA) 
  // CONSTRAINTS ARE DONE SEPARATELY FOR EACH DILEPTON TYPE
 
  RooAddition twoTightMuYieldSum(twoTightMuName+"_YieldSum",twoTightMuName+"_YieldSum",
				 RooArgSet(*wspace.arg(twoTightMuSignalYieldName.Data()),twoTightMuTopWJetsYield));
  RooAddition twoTightMuLooseMuYieldSum(twoTightMuLooseMuName+"_YieldSum",twoTightMuLooseMuName+"_YieldSum",
					RooArgSet(*wspace.arg(twoTightMuLooseMuSignalYieldName.Data()),twoTightMuLooseMuTopWJetsYield));
  RooAddition twoLooseMuYieldSum(twoLooseMuName+"_YieldSum",twoLooseMuName+"_YieldSum",
				 RooArgSet(*wspace.arg(twoLooseMuSignalYieldName.Data()),twoLooseMuTopWJetsYield));
  RooAddition twoTightMuLooseEYieldSum(twoTightMuLooseEName+"_YieldSum",twoTightMuLooseEName+"_YieldSum",
				       RooArgSet(*wspace.arg(twoTightMuLooseESignalYieldName.Data()),twoTightMuLooseETopWJetsYield));
  RooAddition twoLooseMuLooseEYieldSum(twoLooseMuLooseEName+"_YieldSum",twoLooseMuLooseEName+"_YieldSum",
				       RooArgSet(*wspace.arg(twoLooseMuLooseESignalYieldName.Data()),twoLooseMuLooseETopWJetsYield));
  RooAddition twoLooseEYieldSum(twoLooseEName+"_YieldSum",twoLooseEName+"_YieldSum",
				RooArgSet(*wspace.arg(twoLooseESignalYieldName.Data()),twoLooseETopWJetsYield));
  

  RooPoisson twoTightMuConstraint(twoTightMuName+"_Constraint",twoTightMuName+"_Constraint",twoTightMuCount,twoTightMuYieldSum);
  wspace.import( twoTightMuConstraint,RecycleConflictNodes() );

  RooPoisson twoTightMuLooseMuConstraint(twoTightMuLooseMuName+"_Constraint",twoTightMuLooseMuName+"_Constraint",twoTightMuLooseMuCount,twoTightMuLooseMuYieldSum);
  wspace.import( twoTightMuLooseMuConstraint,RecycleConflictNodes() );

  RooPoisson twoLooseMuConstraint(twoLooseMuName+"_Constraint",twoLooseMuName+"_Constraint",twoLooseMuCount,twoLooseMuYieldSum);
  wspace.import( twoLooseMuConstraint,RecycleConflictNodes() );

  RooPoisson twoTightMuLooseEConstraint(twoTightMuLooseEName+"_Constraint",twoTightMuLooseEName+"_Constraint",twoTightMuLooseECount,twoTightMuLooseEYieldSum);
  wspace.import( twoTightMuLooseEConstraint,RecycleConflictNodes() );

  RooPoisson twoLooseMuLooseEConstraint(twoLooseMuLooseEName+"_Constraint",twoLooseMuLooseEName+"_Constraint",twoLooseMuLooseECount,twoLooseMuLooseEYieldSum);
  wspace.import( twoLooseMuLooseEConstraint,RecycleConflictNodes() );

  RooPoisson twoLooseEConstraint(twoLooseEName+"_Constraint",twoLooseEName+"_Constraint",twoLooseECount,twoLooseEYieldSum);
  wspace.import( twoLooseEConstraint,RecycleConflictNodes() );



  /////////////////////////////////////////////////////////////////////////////////////
  // PRODUCE 0L PREDICTIONS FOR MUMU, EMU, EE DILEPTONS (NOT EVER BINNED IN DTHETA) 
  // NOTE THAT DILEPTON TYPE BINS ARE COMBINED BEFORE SCALE FACTOR IS APPLIED

  TString MuMuYieldName("MuMu_");
  MuMuYieldName+=binname;
  TString EMuYieldName("EMu_");
  EMuYieldName+=binname;

  RooAddition MuMuTopWJetsYield(MuMuYieldName+"_TopWJetsYield",MuMuYieldName+"_TopWJetsYield",
				RooArgSet(twoTightMuTopWJetsYield,twoTightMuLooseMuTopWJetsYield,twoLooseMuTopWJetsYield));
  RooAddition EMuTopWJetsYield(EMuYieldName+"_TopWJetsYield",EMuYieldName+"_TopWJetsYield",
				RooArgSet(twoTightMuLooseETopWJetsYield,twoLooseMuLooseETopWJetsYield));


  RooProduct  zeroLeptonTopWJetsYieldMuMu("zeroLepton_"+binname+"_MuMu_TopWJetsYield","zeroLepton_"+binname+"_MuMu_TopWJetsYield",
					  RooArgSet( *wspace.arg(MuMuScaleFactorName.Data()), MuMuTopWJetsYield ));
  wspace.import(zeroLeptonTopWJetsYieldMuMu,RecycleConflictNodes());
 
  RooProduct  zeroLeptonTopWJetsYieldEMu("zeroLepton_"+binname+"_EMu_TopWJetsYield","zeroLepton_"+binname+"_EMu_TopWJetsYield",
					 RooArgSet( *wspace.arg(EMuScaleFactorName.Data()), EMuTopWJetsYield ));
  wspace.import(zeroLeptonTopWJetsYieldEMu,RecycleConflictNodes());
 
  RooProduct  zeroLeptonTopWJetsYieldEE("zeroLepton_"+binname+"_EE_TopWJetsYield","zeroLepton_"+binname+"_EE_TopWJetsYield",
					RooArgSet( *wspace.arg(EEScaleFactorName.Data()), twoLooseETopWJetsYield ));
  wspace.import(zeroLeptonTopWJetsYieldEE,RecycleConflictNodes());


}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
// MAKE TAU->HAD AND DITAU PREDICTIONS FOR THIS NB/MET/HT BIN 

void makeTauHadBinPrediction( RooWorkspace& wspace, TString binname, TString binname_outside ){

  //  cout << " IN TAU PREDICTION" << endl;
  
  TString oneTightMuName("oneTightMu_");
  oneTightMuName+=binname;
  TString twoTightMuName("twoTightMu_");
  twoTightMuName+=binname;
  TString twoTightMuLooseMuName("twoTightMuLooseMu_");
  twoTightMuLooseMuName+=binname;
  TString twoTightMuLooseEName("twoTightMuLooseE_");
  twoTightMuLooseEName+=binname;
  
  ///////////////////////////////////////////////////////////////////////////////////
  // THE 1 TIGHT MU PIECES ARE ALREADY INITIALIZED AND CONSTRAINED IN BINS OF DTHETA
  
  //  cout << endl << endl << " GETTING DTHETA-BINNED TTBAR YIELDS " << endl << endl;
  
  TString oneTightMuTopWJetsYieldName1(oneTightMuName+"_Theta1_TopWJetsYield");	
  //cout << (*wspace.var(oneTightMuTopWJetsYieldName1.Data())).getValV() << endl;;
  
  TString oneTightMuTopWJetsYieldName2(oneTightMuName+"_Theta2_TopWJetsYield");	
  //cout << (*wspace.var(oneTightMuTopWJetsYieldName2.Data())).getValV() << endl;;	
  
  TString oneTightMuTopWJetsYieldName3(oneTightMuName+"_Theta3_TopWJetsYield");	
  //cout << (*wspace.var(oneTightMuTopWJetsYieldName3.Data())).getValV() << endl;;	
  
  TString oneTightMuTopWJetsYieldName4(oneTightMuName+"_Theta4_TopWJetsYield");	
  //cout << (*wspace.var(oneTightMuTopWJetsYieldName4.Data())).getValV() << endl;;	

  TString oneTightMuTopWJetsYieldName5(oneTightMuName+"_Theta5_TopWJetsYield");	
  //cout << (*wspace.var(oneTightMuTopWJetsYieldName5.Data())).getValV() << endl;
  
  /////////////////////////////////////////////////////////////////////////////////////
  // PRODUCE DITAU PREDICTIONS FROM "THIS BIN" TO "BINNAME" (NOT EVER BINNED IN DTHETA) 

  //  cout << " GETTING DILEPTON YIELDS " << endl;
  
  TString twoTightMuTopWJetsYieldName(twoTightMuName+"_TopWJetsYield");	
  //cout << (*wspace.var(twoTightMuTopWJetsYieldName.Data())).getValV() << endl;;
  
  TString twoTightMuLooseMuTopWJetsYieldName(twoTightMuLooseMuName+"_TopWJetsYield");	
  //cout << (*wspace.var(twoTightMuLooseMuTopWJetsYieldName.Data())).getValV() << endl;;
  
  TString twoTightMuLooseETopWJetsYieldName(twoTightMuLooseEName+"_TopWJetsYield");	
  //cout << (*wspace.var(twoTightMuLooseETopWJetsYieldName.Data())).getValV() << endl;;
  

  /////////////////////////////////////////////////////////////////////////////////////
  // SUM MU SAMPLES

  RooAddition oneTauInput( "oneTauInput_"+binname,"oneTauInput_"+binname,
			   RooArgSet(*wspace.var(oneTightMuTopWJetsYieldName1.Data()),*wspace.var(oneTightMuTopWJetsYieldName2.Data()),
				     *wspace.var(oneTightMuTopWJetsYieldName4.Data()),*wspace.var(oneTightMuTopWJetsYieldName3.Data()),
				     *wspace.var(oneTightMuTopWJetsYieldName5.Data())) );
  
  //RooRealVar twoTauInput( "twoTauInput_"+binname,"twoTauInput_"+binname, twoTightMuTopWJetsYield );
  
  RooAddition twoMuTauInput( "twoMuTauInput_"+binname,"twoMuTauInput_"+binname,
			     RooArgSet(*wspace.var(twoTightMuTopWJetsYieldName.Data()),
				       *wspace.var(twoTightMuLooseMuTopWJetsYieldName.Data())) );
  
  //RooRealVar twoETauInput( "twoETauInput_"+binname,"twoETauInput_"+binname, *twoTightMuLooseETopWJetsYield );
  
  
  /////////////////////////////////////////////////////////////////////////////////////
  // GET COMMON MU->TAU SCALE FACTORS
    
  TString oneTauSFName("oneTau_");
  oneTauSFName+=binname_outside;
  oneTauSFName.Append("_ScaleFactor");
  TString twoTauSFName("twoTau_");
  twoTauSFName+=binname_outside;
  twoTauSFName.Append("_ScaleFactor");
  TString twoMuTauSFName("twoMuTau_");
  twoMuTauSFName+=binname_outside;
  twoMuTauSFName.Append("_ScaleFactor");
  TString twoETauSFName("twoETau_");
  twoETauSFName+=binname_outside;
  twoETauSFName.Append("_ScaleFactor");
  
  /*
  cout << " tauhad MC Scale Factors: " << endl;
  (*wspace.arg(oneTauSFName.Data())).Print();
  (*wspace.arg(twoTauSFName.Data())).Print();
  (*wspace.arg(twoMuTauSFName.Data())).Print();
  (*wspace.arg(twoETauSFName.Data())).Print();
  */

  //cout << " APPLYING SCALE FACTORS TO SUMS " << endl;

  RooProduct  zeroLeptonTopWJetsYield1Tau("zeroLepton_"+binname+"_1Tau_TopWJetsYield","zeroLepton_"+binname+"_1Tau_TopWJetsYield",
					  RooArgSet( *wspace.arg(oneTauSFName.Data()), oneTauInput ));
  
  RooProduct  zeroLeptonTopWJetsYield2Tau("zeroLepton_"+binname+"_2Tau_TopWJetsYield","zeroLepton_"+binname+"_2Tau_TopWJetsYield",
					  RooArgSet( *wspace.arg(twoTauSFName.Data()), *wspace.var(twoTightMuTopWJetsYieldName.Data()) ));
  
  RooProduct  zeroLeptonTopWJetsYieldMuTau("zeroLepton_"+binname+"_MuTau_TopWJetsYield","zeroLepton_"+binname+"_MuTau_TopWJetsYield",
					   RooArgSet( *wspace.arg(twoMuTauSFName.Data()), twoMuTauInput ));

  RooProduct  zeroLeptonTopWJetsYieldETau("zeroLepton_"+binname+"_ETau_TopWJetsYield","zeroLepton_"+binname+"_ETau_TopWJetsYield",
					  RooArgSet( *wspace.arg(twoETauSFName.Data()), *wspace.var(twoTightMuLooseETopWJetsYieldName.Data()) ));  
  
  /////////////////////////////////////////////////////////////////////////////////////
  // ADD TAU->HAD AND DITAU CONTRIBUTIONS TO 0L PREDICTION TO WORKSPACE
  
  //cout << " IMPORTING PREDICTIONS " << endl;
  /*
  zeroLeptonTopWJetsYield1Tau.Print();
  zeroLeptonTopWJetsYield2Tau.Print();
  zeroLeptonTopWJetsYieldMuTau.Print();
  zeroLeptonTopWJetsYieldETau.Print();
  */
  wspace.import(zeroLeptonTopWJetsYield1Tau,RecycleConflictNodes());
  wspace.import(zeroLeptonTopWJetsYield2Tau,RecycleConflictNodes());
  wspace.import(zeroLeptonTopWJetsYieldMuTau,RecycleConflictNodes());
  wspace.import(zeroLeptonTopWJetsYieldETau,RecycleConflictNodes());
  
  
}// END OF TAU->HAD FUNCTION FOR ONE BIN



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


void buildMRLikelihood( TString outputFile, TString setupFileName ) 
{

  ///////////////////////////////////////////////////

  int error = 0;
  //  TFile *test = new TFile(outputFile.Data(),"RECREATE");

  ofstream test;
  test.open(outputFile.Data(),ios::app);
  

  ///////////////////////////////////////////////////

  RooWorkspace* wspace = new RooWorkspace("wspace");

  RooRealVar signalCrossSection("signalCrossSection","signalCrossSection",0.1,0.,10.);
  wspace->import(signalCrossSection);

  wspace->defineSet("namesfordata","");
  wspace->defineSet("nuisances","");
  vector<TString> allbinnames;
  allbinnames.clear();

  // DEFINE BINS AND LOCATION OF BIN INPUTS
  vector<TString> binnames;
  map<TString,TString> binnamesoutside;
  map<TString,TString> binContentFileNames;
  map<TString,TString> binPolarizationScaleFactorFileNames;
  map<TString,TString> binTauHadScaleFactorFileNames;
  map<TString,TString> binSignalFracFileNames;


  ////////////////////////////////////
  ////// READ IN FILE WITH LIST OF BIN NAMES AND THE CORRESPONDING INPUT FILES
  ////// FOR CONTENT AND SCALE FACTORS (repeated scale factors will be specified there)

  ifstream setupFile;

  string fileLine;
  TString index, outsidename, fileName1, fileName2, fileName3, fileName4;

  setupFile.open(setupFileName.Data(),fstream::in);


  while(!setupFile.eof())
    {
      getline(setupFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken listOfFiles(thisLine," ");
      listOfFiles.NextToken();
      index = listOfFiles;
      listOfFiles.NextToken();
      fileName1 = listOfFiles;
      listOfFiles.NextToken();
      outsidename = listOfFiles;
      listOfFiles.NextToken();
      fileName2 = listOfFiles;
      listOfFiles.NextToken();
      fileName3 = listOfFiles;
      listOfFiles.NextToken();
      fileName4 = listOfFiles;

      //cout << index << " : " << fileName << endl;
      if(index == "") continue;
    
      binnames.push_back(index);
      binContentFileNames[index] = fileName1;
      binnamesoutside[index] = outsidename;
      binPolarizationScaleFactorFileNames[index] = fileName2;
      binTauHadScaleFactorFileNames[index] = fileName3;
      binSignalFracFileNames[index] = fileName4;

    }
  
  setupFile.close();
  

  ////////////////////////////////////
  ////// READ IN SIGNAL YIELDS FROM OUTSIDE, ADD TO WORKSPACE

  for( int i=0; i<binnames.size(); i++ )
    {

      TString thisBin = binnames.at(i);

      TString signalfracFileName = binSignalFracFileNames[thisBin];
      ifstream signalfracFile;
      signalfracFile.open(signalfracFileName.Data(),fstream::in);
      double count = 0.;

      //////

      TString zeroLeptonName("zeroLepton_");
      zeroLeptonName+=thisBin;

      TString oneTightMuName("oneTightMu_");
      oneTightMuName+=thisBin;
      TString oneLooseMuName("oneLooseMu_");
      oneLooseMuName+=thisBin;
      TString oneLooseEName("oneLooseE_");
      oneLooseEName+=thisBin;

      TString twoTightMuName("twoTightMu_");
      twoTightMuName+=thisBin;
      TString twoTightMuLooseMuName("twoTightMuLooseMu_");
      twoTightMuLooseMuName+=thisBin;
      TString twoLooseMuName("twoLooseMu_");
      twoLooseMuName+=thisBin;
      TString twoTightMuLooseEName("twoTightMuLooseE_");
      twoTightMuLooseEName+=thisBin;
      TString twoLooseMuLooseEName("twoLooseMuLooseE_");
      twoLooseMuLooseEName+=thisBin;
      TString twoLooseEName("twoLooseE_");
      twoLooseEName+=thisBin;

      ////// 

      signalfracFile>>count;
      count = count/10.;

      RooRealVar acount0(zeroLeptonName+"_SignalFrac",zeroLeptonName+"_SignalFrac",count); 
      RooProduct zeroLepton(zeroLeptonName+"_SignalYield",zeroLeptonName+"_SignalYield",RooArgSet(acount0,signalCrossSection));
      wspace->import(zeroLepton,RecycleConflictNodes());
      
      TString outputthis(zeroLeptonName);
      outputthis.Append("_SignalYield");
      //cout << endl << endl << " ZERO LEPTON SIGNAL YIELD MADE HERE AS " << outputthis.Data() << "  " << zeroLepton.getValV() << endl << endl;

      //////
      
      for( int j=1; j<6; j++ ){
	
	TString oneTightMuThetaName(oneTightMuName+"_Theta");
	oneTightMuThetaName+=j;
	TString oneLooseMuThetaName(oneLooseMuName+"_Theta");
	oneLooseMuThetaName+=j;
	TString oneLooseEThetaName(oneLooseEName+"_Theta");
	oneLooseEThetaName+=j;

	signalfracFile>>count;
	count = count/10.;
	RooRealVar acount1(oneTightMuThetaName+"_SignalFrac",oneTightMuThetaName+"_SignalFrac",count); 
	RooProduct oneTightMuTheta(oneTightMuThetaName+"_SignalYield",oneTightMuThetaName+"_SignalYield",RooArgSet(acount1,signalCrossSection));
	wspace->import(oneTightMuTheta,RecycleConflictNodes());
		
	signalfracFile>>count;
	count = count/10.;
	RooRealVar acount2(oneLooseMuThetaName+"_SignalFrac",oneLooseMuThetaName+"_SignalFrac",count); 
	RooProduct oneLooseMuTheta(oneLooseMuThetaName+"_SignalYield",oneLooseMuThetaName+"_SignalYield",RooArgSet(acount2,signalCrossSection));
	wspace->import(oneLooseMuTheta,RecycleConflictNodes());
		
	signalfracFile>>count;
	count = count/10.;
	RooRealVar acount3(oneLooseEThetaName+"_SignalFrac",oneLooseEThetaName+"_SignalFrac",count); 
	RooProduct oneLooseETheta(oneLooseEThetaName+"_SignalYield",oneLooseEThetaName+"_SignalYield",RooArgSet(acount3,signalCrossSection));
	wspace->import(oneLooseETheta,RecycleConflictNodes());
		
      }
      
      signalfracFile>>count;
      count = count/10.;
      RooRealVar acount4(twoTightMuName+"_SignalFrac",twoTightMuName+"_SignalFrac",count); 
      RooProduct twoTightMu(twoTightMuName+"_SignalYield",twoTightMuName+"_SignalYield",RooArgSet(acount4,signalCrossSection));
      wspace->import(twoTightMu,RecycleConflictNodes());
      
      signalfracFile>>count;
      count = count/10.;
      RooRealVar acount5(twoTightMuLooseMuName+"_SignalFrac",twoTightMuLooseMuName+"_SignalFrac",count); 
      RooProduct twoTightMuLooseMu(twoTightMuLooseMuName+"_SignalYield",twoTightMuLooseMuName+"_SignalYield",RooArgSet(acount5,signalCrossSection));
      wspace->import(twoTightMuLooseMu,RecycleConflictNodes());
      
      signalfracFile>>count;
      count = count/10.;
      RooRealVar acount6(twoLooseMuName+"_SignalFrac",twoLooseMuName+"_SignalFrac",count); 
      RooProduct twoLooseMu(twoLooseMuName+"_SignalYield",twoLooseMuName+"_SignalYield",RooArgSet(acount6,signalCrossSection));
      wspace->import(twoLooseMu,RecycleConflictNodes());
      
      signalfracFile>>count;
      count = count/10.;
      RooRealVar acount7(twoTightMuLooseEName+"_SignalFrac",twoTightMuLooseEName+"_SignalFrac",count); 
      RooProduct twoTightMuLooseE(twoTightMuLooseEName+"_SignalYield",twoTightMuLooseEName+"_SignalYield",RooArgSet(acount7,signalCrossSection));
      wspace->import(twoTightMuLooseE,RecycleConflictNodes());
      
      signalfracFile>>count;
      count = count/10.;
      RooRealVar acount8(twoLooseMuLooseEName+"_SignalFrac",twoLooseMuLooseEName+"_SignalFrac",count); 
      RooProduct twoLooseMuLooseE(twoLooseMuLooseEName+"_SignalYield",twoLooseMuLooseEName+"_SignalYield",RooArgSet(acount8,signalCrossSection));
      wspace->import(twoLooseMuLooseE,RecycleConflictNodes());
      
      signalfracFile>>count;
      count = count/10.;
      RooRealVar acount9(twoLooseEName+"_SignalFrac",twoLooseEName+"_SignalFrac",count); 
      RooProduct twoLooseE(twoLooseEName+"_SignalYield",twoLooseEName+"_SignalYield",RooArgSet(acount9,signalCrossSection));
      wspace->import(twoLooseE,RecycleConflictNodes());

          
    }// end of loop over bins setting 



  ////////////////////////////////////
  ////// READ IN CONTENT OF ALL BINS, ADD TO WORKSPACE

  //  for(vector<TString>::iterator thisBin = binnames.begin(); thisBin != binnames.end() ; thisBin++)
  for( int i=0; i<binnames.size(); i++ )
    {

      TString thisBin = binnames.at(i);

      TString contentFileName = binContentFileNames[thisBin];
      ifstream contentFile;
      contentFile.open(contentFileName.Data(),fstream::in);
      double count = 0.;

      //////

      TString zeroLeptonName("zeroLepton_");
      zeroLeptonName+=thisBin;

      TString oneTightMuName("oneTightMu_");
      oneTightMuName+=thisBin;
      TString oneLooseMuName("oneLooseMu_");
      oneLooseMuName+=thisBin;
      TString oneLooseEName("oneLooseE_");
      oneLooseEName+=thisBin;

      TString twoTightMuName("twoTightMu_");
      twoTightMuName+=thisBin;
      TString twoTightMuLooseMuName("twoTightMuLooseMu_");
      twoTightMuLooseMuName+=thisBin;
      TString twoLooseMuName("twoLooseMu_");
      twoLooseMuName+=thisBin;
      TString twoTightMuLooseEName("twoTightMuLooseE_");
      twoTightMuLooseEName+=thisBin;
      TString twoLooseMuLooseEName("twoLooseMuLooseE_");
      twoLooseMuLooseEName+=thisBin;
      TString twoLooseEName("twoLooseE_");
      twoLooseEName+=thisBin;

      ////// 

      contentFile>>count;
      TString name1(zeroLeptonName+"_SignalFrac");
      //count = count + 10.0*(wspace->var(name1.Data())->getValV());
      RooRealVar zeroLeptonCount(zeroLeptonName+"_Count",zeroLeptonName+"_Count",count);
      zeroLeptonCount.setConstant();
      wspace->import(zeroLeptonCount,RecycleConflictNodes());
      wspace->extendSet("namesfordata",zeroLeptonCount.GetName());
      
      //////
      
      for( int j=1; j<6; j++ ){
	
	TString oneTightMuThetaName(oneTightMuName+"_Theta");
	oneTightMuThetaName+=j;
	TString oneLooseMuThetaName(oneLooseMuName+"_Theta");
	oneLooseMuThetaName+=j;
	TString oneLooseEThetaName(oneLooseEName+"_Theta");
	oneLooseEThetaName+=j;
	
	contentFile>>count;
	TString name2(oneTightMuThetaName+"_SignalFrac");
	//count = count + 10.0*(wspace->var(name2.Data())->getValV());
 	RooRealVar oneTightMuThetaCount(oneTightMuThetaName+"_Count",oneTightMuThetaName+"_Count",count);
	oneTightMuThetaCount.setConstant();
	wspace->import(oneTightMuThetaCount,RecycleConflictNodes());
	wspace->extendSet("namesfordata",oneTightMuThetaCount.GetName());
	
	contentFile>>count;
	TString name3(oneLooseMuThetaName+"_SignalFrac");
	//count = count + 10.0*(wspace->var(name3.Data())->getValV());
	RooRealVar oneLooseMuThetaCount(oneLooseMuThetaName+"_Count",oneLooseMuThetaName+"_Count",count);
	oneLooseMuThetaCount.setConstant();
	wspace->import(oneLooseMuThetaCount,RecycleConflictNodes());
	wspace->extendSet("namesfordata",oneLooseMuThetaCount.GetName());
	
	contentFile>>count;
	TString name4(oneLooseEThetaName+"_SignalFrac");
	//count = count + 10.0*(wspace->var(name4.Data())->getValV());
	RooRealVar oneLooseEThetaCount(oneLooseEThetaName+"_Count",oneLooseEThetaName+"_Count",count);
	oneLooseEThetaCount.setConstant();
	wspace->import(oneLooseEThetaCount,RecycleConflictNodes());
	wspace->extendSet("namesfordata",oneLooseEThetaCount.GetName());
	
      }
      
      contentFile>>count;
      TString name5(twoTightMuName+"_SignalFrac");
      //count = count + 10.0*(wspace->var(name5.Data())->getValV());
      RooRealVar twoTightMuCount(twoTightMuName+"_Count",twoTightMuName+"_Count",count);
      twoTightMuCount.setConstant();
      wspace->import(twoTightMuCount,RecycleConflictNodes());
      wspace->extendSet("namesfordata",twoTightMuCount.GetName());
      
      contentFile>>count;
      TString name6(twoTightMuLooseMuName+"_SignalFrac");
      //count = count + 10.0*(wspace->var(name6.Data())->getValV());
      RooRealVar twoTightMuLooseMuCount(twoTightMuLooseMuName+"_Count",twoTightMuLooseMuName+"_Count",count);
      twoTightMuLooseMuCount.setConstant();
      wspace->import(twoTightMuLooseMuCount,RecycleConflictNodes());
      wspace->extendSet("namesfordata",twoTightMuLooseMuCount.GetName());
      
      contentFile>>count;
      TString name7(twoLooseMuName+"_SignalFrac");
      //count = count + 10.0*(wspace->var(name7.Data())->getValV());
      RooRealVar twoLooseMuCount(twoLooseMuName+"_Count",twoLooseMuName+"_Count",count);
      twoLooseMuCount.setConstant();
      wspace->import(twoLooseMuCount,RecycleConflictNodes());
      wspace->extendSet("namesfordata",twoLooseMuCount.GetName());
      
      contentFile>>count;
      TString name8(twoTightMuLooseEName+"_SignalFrac");
      //count = count + 10.0*(wspace->var(name8.Data())->getValV());
      RooRealVar twoTightMuLooseECount(twoTightMuLooseEName+"_Count",twoTightMuLooseEName+"_Count",count);
      twoTightMuLooseECount.setConstant();
      wspace->import(twoTightMuLooseECount,RecycleConflictNodes());
      wspace->extendSet("namesfordata",twoTightMuLooseECount.GetName());
      
      contentFile>>count;
      TString name9(twoLooseMuLooseEName+"_SignalFrac");
      //count = count + 10.0*(wspace->var(name9.Data())->getValV());
      RooRealVar twoLooseMuLooseECount(twoLooseMuLooseEName+"_Count",twoLooseMuLooseEName+"_Count",count);
      twoLooseMuLooseECount.setConstant();
      wspace->import(twoLooseMuLooseECount,RecycleConflictNodes());
      wspace->extendSet("namesfordata",twoLooseMuLooseECount.GetName());
      
      contentFile>>count;
      TString name0(twoLooseEName+"_SignalFrac");
      //count = count + 10.0*(wspace->var(name0.Data())->getValV());
      RooRealVar twoLooseECount(twoLooseEName+"_Count",twoLooseEName+"_Count",count);
      twoLooseECount.setConstant();
      wspace->import(twoLooseECount,RecycleConflictNodes());
      wspace->extendSet("namesfordata",twoLooseECount.GetName());

      ///////////////////////////////////////////////////////////////////////////////////////////////////////

      contentFile.close();

          
    }// end of loop over bins setting counts


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////
  ////// READ IN SCALE FACTORS FOR ALL BINS, ADD TO WORKSPACE


  //  for(vector<TString>::iterator thisBin = binnames.begin(); thisBin != binnames.end() ; thisBin++)
  for( int i=0; i<binnames.size(); i++ )
    {
      
      TString thisBin = binnames.at(i);
      
      double scalefactor = 0.;
      double scalefactorerror = 0.;

      TString scaleFactorFileName = binPolarizationScaleFactorFileNames[thisBin];
      ifstream scaleFactorFile;
      scaleFactorFile.open(scaleFactorFileName.Data(),fstream::in);

      TString tauhadScaleFactorFileName = binTauHadScaleFactorFileNames[thisBin];
      ifstream tauhadScaleFactorFile; 
      tauhadScaleFactorFile.open(tauhadScaleFactorFileName.Data(),fstream::in);


      //      string fileLine;
 
      ////////////////////////////////////
      ////// READ IN POLARIZATION METHOD SCALE FACTORS, ADD TO WORKSPACE

      for( int j=1; j<6; j++ ){

	scaleFactorFile >> scalefactor >> scalefactorerror ;

	TString oneLeptonName("oneLepton_");
	oneLeptonName+=thisBin;
	oneLeptonName.Append("_Theta");
	oneLeptonName+=j;
	oneLeptonName.Append("_ScaleFactor");

	TString name1("blah1");
	name1+=j;
	TString name2("blah2");
	name2+=j;

	if( scalefactor<0.0001 ){
	  scalefactor=0.0001; scalefactorerror=0.0001; 
	}

	RooAbsArg* checksf0 = wspace->arg( oneLeptonName.Data() );
	if(checksf0 == NULL)
	  {
	    RooAbsArg* ScaleFactor =
	      getBetaPrimeConstraint( *wspace, oneLeptonName, "",
				      scalefactor, scalefactorerror,
				      name1, name2 );
	    wspace->import(*ScaleFactor,RecycleConflictNodes());
	  }

      }

      ////////////////////////////////////
      ////// READ IN DILEPTON METHOD SCALE FACTORS, ADD TO WORKSPACE

      scaleFactorFile>>scalefactor>>scalefactorerror;
      if( scalefactor<0.0001 ){
	scalefactor=0.0001; scalefactorerror=0.0001; 
      }
      
      TString mumuName("MuMu_");
      mumuName+=binnamesoutside[thisBin];
      mumuName.Append("_ScaleFactor");

      RooAbsArg* checksf1 = wspace->arg( mumuName.Data() );
      if(checksf1 == NULL){

	RooAbsArg* ScaleFactorMuMu =
	  getBetaPrimeConstraint( *wspace, mumuName, "",
				  scalefactor, scalefactorerror,
				  "blah1a", "blah2a" );
	wspace->import(*ScaleFactorMuMu,RecycleConflictNodes());
      }
      
      
      scaleFactorFile>>scalefactor>>scalefactorerror;
      if( scalefactor<0.0001 ){
	scalefactor=0.0001; scalefactorerror=0.0001; 
      }

      TString emuName("EMu_");
      emuName+=binnamesoutside[thisBin];
      emuName.Append("_ScaleFactor");

      RooAbsArg* checksf2 = wspace->arg( emuName.Data() );
      if(checksf2 == NULL){

	RooAbsArg *ScaleFactorEMu =
	  getBetaPrimeConstraint( *wspace, emuName, "",
				  scalefactor, scalefactorerror,
				  "blah1b", "blah2b" );
	wspace->import(*ScaleFactorEMu,RecycleConflictNodes());
      }

      scaleFactorFile>>scalefactor>>scalefactorerror;
      if( scalefactor<0.0001 ){
	scalefactor=0.0001; scalefactorerror=0.0001; 
      }

      TString eeName("EE_");
      eeName+=binnamesoutside[thisBin];
      eeName.Append("_ScaleFactor");

      RooAbsArg* checksf3 = wspace->arg( eeName.Data() );
      if(checksf3 == NULL){

	RooAbsArg *ScaleFactorEE =
	  getBetaPrimeConstraint( *wspace, eeName, "",
				  scalefactor, scalefactorerror,
				  "blah1c", "blah2c" );
	wspace->import(*ScaleFactorEE,RecycleConflictNodes());
      }
      


      //cout << " Done importing ee scale factors " << endl;

      ////////////////////////////////////
      ////// READ IN TAU/DITAU MC-BASED SCALE FACTORS, ADD TO WORKSPACE

      tauhadScaleFactorFile>>scalefactor>>scalefactorerror;
      if( scalefactor<0.0001 ){
	scalefactor=0.0001; scalefactorerror=0.0001; 
      }

      TString oneTauName("oneTau_");
      oneTauName+=binnamesoutside[thisBin];
      oneTauName.Append("_ScaleFactor");

      RooAbsArg* checksf4 = wspace->arg( oneTauName.Data() );
      if(checksf4 == NULL){

	RooAbsArg* ScaleFactorTau =
	  getBetaPrimeConstraint( *wspace, oneTauName, "",
				  scalefactor, scalefactorerror,
				  "blah1t", "blah2t" );
	wspace->import(*ScaleFactorTau,RecycleConflictNodes());
	
      }

      //////

      tauhadScaleFactorFile>>scalefactor>>scalefactorerror;
      if( scalefactor<0.0001 ){
	scalefactor=0.0001; scalefactorerror=0.0001; 
      }

      TString twoMuTauName("twoMuTau_");
      twoMuTauName+=binnamesoutside[thisBin];
      twoMuTauName.Append("_ScaleFactor");

      RooAbsArg* checksf5 = wspace->arg( twoMuTauName.Data() );
      if(checksf5 == NULL){

	RooAbsArg* ScaleFactorMuTau =
	  getBetaPrimeConstraint( *wspace, twoMuTauName, "",
				  scalefactor, scalefactorerror,
				  "blah1mt", "blah2mt" );
	wspace->import(*ScaleFactorMuTau,RecycleConflictNodes());
      }

      //////

      tauhadScaleFactorFile>>scalefactor>>scalefactorerror;
      if( scalefactor<0.0001 ){
	scalefactor=0.0001; scalefactorerror=0.0001; 
      }

      TString twoETauName("twoETau_");
      twoETauName+=binnamesoutside[thisBin];
      twoETauName.Append("_ScaleFactor");

      RooAbsArg* checksf6 = wspace->arg( twoETauName.Data() );
      if(checksf6 == NULL){

	RooAbsArg* ScaleFactorETau =
	  getBetaPrimeConstraint( *wspace, twoETauName, "",
				  scalefactor, scalefactorerror,
				  "blah1et", "blah2et" );
	wspace->import(*ScaleFactorETau,RecycleConflictNodes());
      }
      
      //////

      tauhadScaleFactorFile>>scalefactor>>scalefactorerror;
      if( scalefactor<0.0001 ){
	scalefactor=0.0001; scalefactorerror=0.0001; 
      }

      TString twoTauName("twoTau_");
      twoTauName+=binnamesoutside[thisBin];
      twoTauName.Append("_ScaleFactor");

      RooAbsArg* checksf7 = wspace->arg( twoTauName.Data() );
      if(checksf7 == NULL){
	
	RooAbsArg* ScaleFactorTauTau =
	  getBetaPrimeConstraint( *wspace, twoTauName, "",
				  scalefactor, scalefactorerror,
				  "blah1tt", "blah2tt" );
	wspace->import(*ScaleFactorTauTau,RecycleConflictNodes());
      }
      
      //////

      scaleFactorFile.close();
      tauhadScaleFactorFile.close();

      ////////////////////////////////////
      ////////////////////////////////////
      ////// PUTTING SOME CONSTRAINT-SETTING CODE HERE FOR NOW:


      makePolarizationConstraintsPredictions( *wspace, thisBin );
      //cout << " end of polarization constraints " << endl;
      makeDileptonConstraintsPredictions( *wspace, thisBin, binnamesoutside[thisBin] );
      //cout << " end of dilepton constraints " << endl;
      makeTauHadBinPrediction( *wspace, thisBin, binnamesoutside[thisBin] );
      //cout << " end of tauhad prediction " << endl;

      makePrediction( *wspace, thisBin );


    }// end of loop over bins


  
  RooArgSet Constraints(wspace->allPdfs());

  RooProdPdf model("model","model", Constraints );
  wspace->import(model);

  //cout << endl << endl << " FINAL SET OF CONSTRAINTS " << endl;
  //Constraints.Print("v");

  wspace->defineSet("poi","signalCrossSection");
  RooDataSet dataset("dataset","dataset",*wspace->set("namesfordata"));
  dataset.add(*wspace->set("namesfordata"));
  wspace->import(dataset);
  
  //cout << endl << endl << " FINAL SET OF DATA " << endl;
  //dataset.Print();

  /*
  
  ////////////////////////////////////

  cout << " MAKING MODEL CONFIG " << endl;

  ModelConfig * modelConfig = new ModelConfig("modelConfig");
  modelConfig->SetWorkspace(*wspace);
  modelConfig->SetPdf("model");
  modelConfig->SetParametersOfInterest(*wspace->set("poi"));
  modelConfig->SetNuisanceParameters(*wspace->set("nuisances"));

  cout << endl << endl << endl << " before model print " << endl;
  modelConfig->Print();

  cout << " after model print " << endl;


  // Declare parameter of interest
  RooArgSet parofinterest(*wspace->set("poi")) ;
  ProfileLikelihoodCalculator plc(dataset, *wspace->pdf("model"), parofinterest);
  //plc.SetTestSize(0.05);
  //ProfileLikelihoodCalculator plc(dataset, *modelConfig);
  plc.SetConfidenceLevel(0.95);

  cout << " before pl drawing" << endl;

  LikelihoodInterval* interval = plc.GetInterval() ;
  LikelihoodIntervalPlot*  lplot = new LikelihoodIntervalPlot(interval);

  //  interval->Print("v");

  RooRealVar* firstPOI = (RooRealVar*) parofinterest.first();

  double profileLikelihoodUpperLimit = interval->UpperLimit(*firstPOI);
  double profileLikelihoodLowerLimit = interval->LowerLimit(*firstPOI);

  cout << " ready to draw pl " << endl;
  
  firstPOI->setRange(profileLikelihoodLowerLimit,profileLikelihoodUpperLimit);
  
  TCanvas* canvas = new TCanvas("Result", "Result", 10, 10,500,500);
  lplot->Draw();
  //lplot->Save();
  canvas->Update();

  */

  cout << " before simple fit " << endl;
  
  ////////////////////////////////////
  model.fitTo(dataset);
 
  cout << " XSEC AFTER FIT " << wspace->var("signalCrossSection")->getValV() << endl;

  test << wspace->var("signalCrossSection")->getValV() << endl;
  test.close(); 

  ////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

  /* 

  cout << endl << endl << endl << endl;
  cout << " XSEC AFTER FIT " << wspace->var("signalCrossSection")->getValV() << endl;
  cout << " SIGFRAC AFTER FIT " << wspace->var("zeroLepton_3BhighHThighMET_SignalFrac")->getValV() << endl << endl;

  cout << " ONE LEPTON COUNTS 1 " << wspace->var("oneLepton_3BhighHThighMET_Theta1Count")->getValV() << endl;
  cout << " ONE LEPTON COUNTS 2 " << wspace->var("oneLepton_3BhighHThighMET_Theta2Count")->getValV() << endl;
  cout << " ONE LEPTON COUNTS 3 " << wspace->var("oneLepton_3BhighHThighMET_Theta3Count")->getValV() << endl;
  cout << " ONE LEPTON COUNTS 4 " << wspace->var("oneLepton_3BhighHThighMET_Theta4Count")->getValV() << endl;
  cout << " ONE LEPTON COUNTS 5 " << wspace->var("oneLepton_3BhighHThighMET_Theta5Count")->getValV() << endl << endl;
  
  cout << " ONE LEPTON SUM 1 " << wspace->function("oneLepton_3BhighHThighMET_YieldSumTheta1")->getVal() << endl;
  cout << " ONE LEPTON SUM 2 " << wspace->function("oneLepton_3BhighHThighMET_YieldSumTheta2")->getVal() << endl;
  cout << " ONE LEPTON SUM 3 " << wspace->function("oneLepton_3BhighHThighMET_YieldSumTheta3")->getVal() << endl;
  cout << " ONE LEPTON SUM 4 " << wspace->function("oneLepton_3BhighHThighMET_YieldSumTheta4")->getVal() << endl;
  cout << " ONE LEPTON SUM 5 " << wspace->function("oneLepton_3BhighHThighMET_YieldSumTheta5")->getVal() << endl << endl;

  cout << " ONE LEPTON ttbar 1 " << wspace->function("oneLepton_3BhighHThighMET_TopWJetsYieldTheta1")->getVal() << endl;
  cout << " ONE LEPTON ttbar 2 " << wspace->function("oneLepton_3BhighHThighMET_TopWJetsYieldTheta2")->getVal() << endl;
  cout << " ONE LEPTON ttbar 3 " << wspace->function("oneLepton_3BhighHThighMET_TopWJetsYieldTheta3")->getVal() << endl;
  cout << " ONE LEPTON ttbar 4 " << wspace->function("oneLepton_3BhighHThighMET_TopWJetsYieldTheta4")->getVal() << endl;
  cout << " ONE LEPTON ttbar 5 " << wspace->function("oneLepton_3BhighHThighMET_TopWJetsYieldTheta5")->getVal() << endl << endl;

  cout << " ONE LEPTON sig 1 " << wspace->function("oneLepton_3BhighHThighMET_SignalYieldTheta1")->getVal() << endl;
  cout << " ONE LEPTON sig 2 " << wspace->function("oneLepton_3BhighHThighMET_SignalYieldTheta2")->getVal() << endl;
  cout << " ONE LEPTON sig 3 " << wspace->function("oneLepton_3BhighHThighMET_SignalYieldTheta3")->getVal() << endl;
  cout << " ONE LEPTON sig 4 " << wspace->function("oneLepton_3BhighHThighMET_SignalYieldTheta4")->getVal() << endl;
  cout << " ONE LEPTON sig 5 " << wspace->function("oneLepton_3BhighHThighMET_SignalYieldTheta5")->getVal() << endl << endl;

  cout << " 0L SIGNAL AFTER FIT " << wspace->function("zeroLepton_3BhighHThighMET_SignalYield")->getVal() << endl;
  cout << " 0L TTBAR AFTER FIT  " << wspace->function("zeroLepton_3BhighHThighMET_TopWJetsYield")->getVal() << endl;
  cout << " ZERO LEPTON SUM " << wspace->function("zeroLepton_3BhighHThighMET_YieldSum")->getVal() << endl;
 
  cout << endl << endl << endl << endl << " zeroLeptonCounts " << wspace->var("zeroLepton_3BhighHThighMET_Count")->getValV() << endl << endl;
  
 cout << endl <<  endl;

 if( error==1 ) cout << " ERROR FOUND IN MATCHING BIN NAMES!!! " << endl;

*/
  



}// end of likelihood builder





