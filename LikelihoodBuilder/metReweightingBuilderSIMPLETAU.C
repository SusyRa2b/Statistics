
/**************************************************************

THIS LIKELIHOOD USES THE SIMPLIFIED TAU (MC SFs ONLY) PROCEDURE,
COMBINES LOOSE MUS WITH Es,
COMBINES ALL DILEPTONS (EE,MUMU,EMU), 
AND COMBINES ALL TAU-DILEPTONS (TAUTAU,MUTAU,ETAU)

USES STANDARD INPUT FILES
*** NEED TO FIX SFs!!! ***

**************************************************************/

#include "TROOT.h"


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

//#include "RooRatio.h"

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
#include "rooFitGaussianHelperFunctions.h"

#include "RooProdPdfLogSum.h"
#include "RooPoissonLogEval.h"
//#include "RooPosDefCorrGauss.h"

#include "RooRatio.h"
#include "RooBetaPdf.h"
#include "RooBetaPrimePdf.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooCorrelatedBetaGeneratorHelper.h"
#include "RooCorrelatedBetaPrimeGeneratorHelper.h"
#include "betaHelperFunctions.h"
#include "RooBetaInverseCDF.h"
#include "RooBetaPrimeInverseCDF.h"
#include "rooFitBetaHelperFunctions.h"
#include "RooNormalFromFlatPdf.h"

#include "TMath.h"

#include "metReweightingBuilderSIMPLETAU.h"

using namespace RooFit ;
using namespace RooStats ;



//////////////////////////////////////////////////////////////////////////////////////////////////////

void makePrediction( RooWorkspace& wspace, TString thisBin, bool standalone ){


      ////////////////////////////////////
      ////////////////////////////////////
      // PUTTING TOGETHER THE 0L PREDICTION FOR THIS BIN

      TString zeroLeptonName("zeroLepton_");
      zeroLeptonName+=thisBin;
      //////
      //cout << " BEFORE GETTING PREDICTION COMPONENTS " << endl;

      TString PolarizationName1(zeroLeptonName);
      PolarizationName1.Append("_Theta1_TopWJetsDataYield");
    
      TString PolarizationName2(zeroLeptonName);
      PolarizationName2.Append("_Theta2_TopWJetsDataYield");
        
      TString PolarizationName3(zeroLeptonName);
      PolarizationName3.Append("_Theta3_TopWJetsDataYield");
        
      TString PolarizationName4(zeroLeptonName);
      PolarizationName4.Append("_Theta4_TopWJetsDataYield");
        
      TString PolarizationName5(zeroLeptonName);
      PolarizationName5.Append("_Theta5_TopWJetsDataYield");
    
      //    cout << " GOT POLARIZATION PREDICTIONS" << endl;
      //////
      TString TauHadName1(zeroLeptonName);
      TauHadName1.Append("_1Tau_TopWJetsDataYield");

      TString TauHadName2(zeroLeptonName);
      TauHadName2.Append("_2Tau_TopWJetsDataYield");
      
      //    cout << " GOT TAU->HAD PREDICTIONS" << endl;
      //////
      TString DilepName1(zeroLeptonName);
      DilepName1.Append("_Dilep_TopWJetsDataYield");
      
      //    cout << " GOT DILEPTON PREDICTIONS" << endl;
      //////
      
      //cout << " GOT ALL COMPONENTS OF 0L PREDICTION, ADDING THEM... " << endl;
      
      RooAddition zeroLeptonTopWJetsPolarizationYield(zeroLeptonName+"_TopWJetsPolarizationYield",zeroLeptonName+"_TopWJetsPolarizationYield",
						      RooArgSet( *wspace.arg(PolarizationName1.Data()),
								 *wspace.arg(PolarizationName2.Data()),
								 *wspace.arg(PolarizationName3.Data()),
								 *wspace.arg(PolarizationName4.Data()),
								 *wspace.arg(PolarizationName5.Data())
								 ));
      

      cout << " pol pieces: ";
      (*wspace.function(PolarizationName1.Data())).Print();
      (*wspace.function(PolarizationName2.Data())).Print();
      (*wspace.function(PolarizationName3.Data())).Print();
      (*wspace.function(PolarizationName4.Data())).Print();
      (*wspace.function(PolarizationName5.Data())).Print();
      zeroLeptonTopWJetsPolarizationYield.Print();
      cout << endl;

      cout << " tau pieces: ";
      (*wspace.arg(TauHadName1.Data())).Print();
      (*wspace.arg(TauHadName2.Data())).Print();
      cout << endl;

      //cout << " dilep pieces: ";
      (*wspace.arg(DilepName1.Data())).Print(); 
      cout << endl;

      //      cout << " BACKGROUND PREDICTION PIECES:  " << zeroLeptonTopWJetsPolarizationYield.getValV() << "  " 
      //   << zeroLeptonTopWJetsTauHadYield.getValV() << "  " << zeroLeptonTopWJetsDileptonYield.getValV() << endl;
      //   << endl;
      
      //zeroLeptonTopWJetsYield for full likelihood 
      if( !standalone )
	{
	  RooAddition zeroLeptonTopWJetsYield(zeroLeptonName+"_TopWJetsYield",zeroLeptonName+"_TopWJetsYield",
					 RooArgSet( zeroLeptonTopWJetsPolarizationYield,
						    *wspace.arg(TauHadName1.Data()),
						    *wspace.arg(TauHadName2.Data()),
						    *wspace.arg(DilepName1.Data())
						    ));
	  wspace.import(zeroLeptonTopWJetsYield, RecycleConflictNodes());
	  cout << "Put " << zeroLeptonTopWJetsYield.GetName() << " into workspace" << endl;
	}



      /////////////////////////////////////////////////////////////////////////////////////
      // GET AND APPLY SIGNAL FRACTION FOR THIS 0L BIN, ADD TO ALL PREDICTIONS
      
      //cout << " GET 0L COUNTS AND SIG FRAC " << endl;

      if( standalone )
	{

	  TString zeroLeptonCountName(zeroLeptonName);
	  zeroLeptonCountName.Append("_Count");
	  
	  TString zeroLeptonSignalYieldName(zeroLeptonName);
	  zeroLeptonSignalYieldName.Append("_SignalYield");
	  
	  
	  
	  RooAddition zeroLeptonYieldSum(zeroLeptonName+"_YieldSum",zeroLeptonName+"_YieldSum",
					 RooArgSet( *wspace.arg(zeroLeptonSignalYieldName.Data()), 
						    zeroLeptonTopWJetsPolarizationYield,
						    *wspace.arg(TauHadName1.Data()),
						    *wspace.arg(TauHadName2.Data()),
						    *wspace.arg(DilepName1.Data())
						    ));
	  
	  
	  cout << " SIGNAL YIELD:  " << zeroLeptonSignalYieldName.Data() << "  " ;
	  (*wspace.arg(zeroLeptonSignalYieldName.Data())).Print();
	  
	  cout << " COUNTS:  " << zeroLeptonCountName.Data() << "  "
	       << (*wspace.var(zeroLeptonCountName.Data())).getValV() << "  " << endl;
	  
	  cout << " TOTAL PREDICTED IN BIN:  ";
	  zeroLeptonYieldSum.Print();
	  
	  
	  /////////////////////////////////////////////////////////////////////////////////////
	  // SET CONSTRAINTS FOR THIS 0L BIN
	  /****** SHOULD LEAVE THIS OUT IN COMPLETE LIKELIHOOD!!! ******/
	  //cout << " SETTING 0L CONSTRAINTS " << endl;
	  
	  //      RooPoisson zeroLeptonConstraint(zeroLeptonName+"_Constraint",zeroLeptonName+"_Constraint",
	  RooPoissonLogEval zeroLeptonConstraint(zeroLeptonName+"_Constraint",zeroLeptonName+"_Constraint",
						 *wspace.var(zeroLeptonCountName.Data()),zeroLeptonYieldSum);
	  
	  wspace.import( zeroLeptonConstraint,RecycleConflictNodes() );
	  
	  
	  cout << endl << endl;
	  zeroLeptonConstraint.Print();
	  cout << endl << endl;
	  
	}// END OF ZERO LEPTON STANDALONE PART     

}// END OF PREDICTIONS


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////



void makePolarizationConstraintsPredictions( RooWorkspace& wspace, TString binname, TString trigeffname ){

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
    
    TString oneLooseLepName("oneLooseLep_");
    oneLooseLepName+=binname;
    oneLooseLepName.Append(ending1);

    ////////////////////////////////////////////////////
    // GET CONTROL SAMPLE COUNTS FOR THIS DTHETA BIN 
    
    TString ending("_Count");

    TString oneTightMuCountName(oneTightMuName);
    oneTightMuCountName.Append(ending);

    TString oneLooseLepCountName(oneLooseLepName);
    oneLooseLepCountName.Append(ending);

    RooRealVar oneTightMuCount = *wspace.var(oneTightMuCountName.Data());
    RooRealVar oneLooseLepCount = *wspace.var(oneLooseLepCountName.Data());


    ////////////////////////////////////////////////////
    // GET SIGNAL YIELDS FOR THIS DTHETA BIN 
 
    //    cout << " getting signal yields for pol method " << endl;

    TString oneTightMuSignalYieldName(oneTightMuName);
    oneTightMuSignalYieldName.Append("_SignalYield");

    TString oneLooseLepSignalYieldName(oneLooseLepName);
    oneLooseLepSignalYieldName.Append("_SignalYield");

 
    ////////////////////////////////////////////////////
    // GET THE 1L->0L SCALE FACTOR FOR THIS DTHETA BIN 
    
    TString oneLeptonScaleFactorName(oneLeptonName);
    oneLeptonScaleFactorName.Append("_ScaleFactor");

    ////////////////////////////////////////////////////
    // INITIALIZE TTBAR COMPONENT OF THIS DTHETA BIN TO START WITH OBSERVED CONTENT
    // EACH YIELD IS A SEPARATE NUISANCE PARAMETER

    double initialvalue = oneTightMuCount.getVal();
    if( initialvalue==0 ) initialvalue = 1.;

    RooRealVar oneTightMuTopWJetsYield(oneTightMuName+"_TopWJetsYield",oneTightMuName+"_TopWJetsYield",
				       initialvalue,0.00001,20000.);
    wspace.import(oneTightMuTopWJetsYield,RecycleConflictNodes());
    wspace.extendSet("nuisances",oneTightMuTopWJetsYield.GetName());

    initialvalue = oneLooseLepCount.getVal();
    if( initialvalue==0 ) initialvalue = 1.;

    RooRealVar oneLooseLepTopWJetsYield(oneLooseLepName+"_TopWJetsYield",oneLooseLepName+"_TopWJetsYield",
				      initialvalue,0.00001,20000.);
    wspace.import(oneLooseLepTopWJetsYield,RecycleConflictNodes());
    wspace.extendSet("nuisances",oneLooseLepTopWJetsYield.GetName());


    ////////////////////////////////////////////////////
    // COMBINE SIGNAL AND TTWT COMPONENTS OF THIS DTHETA BIN AND SET CONTROL SAMPLE CONSTRAINTS
    
    /****** NOTE: INSERT SINGLE LEPTON EFFICIENCIES HERE??? ******/

    RooAddition oneTightMuYieldSum(oneTightMuName+"_YieldSum",oneTightMuName+"_YieldSum",
				   RooArgSet(*wspace.arg(oneTightMuSignalYieldName.Data()),oneTightMuTopWJetsYield));

    RooAddition oneLooseLepYieldSum(oneLooseLepName+"_YieldSum",oneLooseLepName+"_YieldSum",
				  RooArgSet(*wspace.arg(oneLooseLepSignalYieldName.Data()),oneLooseLepTopWJetsYield));
   

    RooPoissonLogEval oneTightMuConstraint(oneTightMuName+"_Constraint",oneTightMuName+"_Constraint",oneTightMuCount,oneTightMuYieldSum);
    wspace.import( oneTightMuConstraint,RecycleConflictNodes() );

    RooPoissonLogEval oneLooseLepConstraint(oneLooseLepName+"_Constraint",oneLooseLepName+"_Constraint",oneLooseLepCount,oneLooseLepYieldSum);
    wspace.import( oneLooseLepConstraint,RecycleConflictNodes() );


    ////////////////////////////////////////////////////
    // CONSTRUCT THE 0L BACKGROUND PREDICTION FROM THIS DTHETA BIN
    // ****** DONT FORGET ABOUT TRIGGER INEFFICIENCIES! ****** //
 
    RooAddition oneLeptonTopWJetsYield(oneLeptonName+"_TopWJetsYield",oneLeptonName+"_TopWJetsYield",
				       RooArgSet(oneTightMuTopWJetsYield,oneLooseLepTopWJetsYield));

    RooProduct  zeroLeptonTopWJetsYield(zeroLeptonName+"_TopWJetsYield",zeroLeptonName+"_TopWJetsYield",
					RooArgSet(*wspace.arg(oneLeptonScaleFactorName.Data()),oneLeptonTopWJetsYield));

    TString triggername = "oneLeptonTriggerEfficiency_";
    triggername.Append(trigeffname);
    RooAbsArg* triggerefficiency = wspace.arg(triggername.Data());

    RooProduct  zeroLeptonTopWJetsDataYield(zeroLeptonName+"_TopWJetsDataYield",zeroLeptonName+"_TopWJetsDataYield",
					    RooArgSet(zeroLeptonTopWJetsYield,*triggerefficiency));

    wspace.import(zeroLeptonTopWJetsDataYield,RecycleConflictNodes());
    wspace.extendSet("nuisances",zeroLeptonTopWJetsDataYield.GetName());


  }// END OF LOOP OVER DTHETA BINS
  
}// END OF BIG BIN IN NBS/HT/MET


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

void makeDileptonConstraintsPredictions( RooWorkspace& wspace, TString binname,  TString binname_outside, TString trigeffname ){


  TString twoTightMuName("twoTightMu_");
  twoTightMuName+=binname;

  TString twoLooseLepName("twoLooseLep_");
  twoLooseLepName+=binname;


  /////////////////////////////////////////////////////////////////////////////////////
  // GET DILEPTON COUNTS
  
  TString twoTightMuCountName(twoTightMuName);
  twoTightMuCountName.Append("_Count");

  TString twoLooseLepCountName(twoLooseLepName);
  twoLooseLepCountName.Append("_Count");

  RooRealVar twoTightMuCount = *wspace.var(twoTightMuCountName.Data());
  RooRealVar twoLooseLepCount = *wspace.var(twoLooseLepCountName.Data());

  //  cout << " END OF GETTING DILEPTON COUNTS " << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // GET DILEPTON SCALE FACTOR
  // SHARED ACROSS HT/MET BINS 

  TString DilepScaleFactorName("Dilep_");
  DilepScaleFactorName+=binname_outside;
  DilepScaleFactorName.Append("_ScaleFactor");

  /////////////////////////////////////////////////////////////////////////////////////
  // GET DILEPTON SIGNAL YIELDS

  TString twoTightMuSignalYieldName(twoTightMuName);
  twoTightMuSignalYieldName.Append("_SignalYield");

  TString twoLooseLepSignalYieldName(twoLooseLepName);
  twoLooseLepSignalYieldName.Append("_SignalYield");


  //  cout << " END OF GETTING DILEPTON SIGNAL YIELDS " << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // SET BACKGROUND CONTRIBUTIONS TO DILEPTON SAMPLES 
  // IMPORT BACKGROUND CONTRIBUTIONS TO DILEPTON SAMPLES TO WORKSPACE AS NUISANCE PARAMS

  // ****** TRIGGER INEFFICIENCIES MUST BE PUT IN SOMEWHERE FOR BACKGROUND AND SIGNAL ******//

  double initialvalue = twoTightMuCount.getVal();
  if( initialvalue==0 ) initialvalue = 1.;
  
  RooRealVar twoTightMuTopWJetsYield(twoTightMuName+"_TopWJetsYield",twoTightMuName+"_TopWJetsYield",
				     initialvalue,0.00001,1000.);
  wspace.import(twoTightMuTopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoTightMuTopWJetsYield.GetName());

  initialvalue = twoLooseLepCount.getVal();
  if( initialvalue==0 ) initialvalue = 1.;

  RooRealVar twoLooseLepTopWJetsYield(twoLooseLepName+"_TopWJetsYield",twoLooseLepName+"_TopWJetsYield",
				    initialvalue,0.00001,1000.);
  wspace.import(twoLooseLepTopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoLooseLepTopWJetsYield.GetName());


  /////////////////////////////////////////////////////////////////////////////////////
  // PRODUCE SIGNAL+BACKGROUND DILEPTON CONSTRAINTS (NOT EVER BINNED IN DTHETA) 
  // CONSTRAINTS ARE DONE SEPARATELY FOR EACH DILEPTON TYPE
 
  RooAddition twoTightMuYieldSum(twoTightMuName+"_YieldSum",twoTightMuName+"_YieldSum",
				RooArgSet( *wspace.arg(twoTightMuSignalYieldName.Data()), twoTightMuTopWJetsYield ));


  RooAddition twoLooseLepYieldSum(twoLooseLepName+"_YieldSum",twoLooseLepName+"_YieldSum",
				RooArgSet( *wspace.arg(twoLooseLepSignalYieldName.Data()), twoLooseLepTopWJetsYield ));
  

  RooPoissonLogEval twoTightMuConstraint(twoTightMuName+"_Constraint",twoTightMuName+"_Constraint",twoTightMuCount,twoTightMuYieldSum);
  wspace.import( twoTightMuConstraint,RecycleConflictNodes() );

  RooPoissonLogEval twoLooseLepConstraint(twoLooseLepName+"_Constraint",twoLooseLepName+"_Constraint",twoLooseLepCount,twoLooseLepYieldSum);
  wspace.import( twoLooseLepConstraint,RecycleConflictNodes() );


  cout << " END OF DILEPTON AND DITAU CONSTRAINTS " << endl;
  twoTightMuConstraint.Print();
  twoLooseLepConstraint.Print();
  cout << endl << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // PRODUCE 0L PREDICTIONS FOR MUMU, EMU, EE DILEPTONS (NOT EVER BINNED IN DTHETA) 
  // NOTE THAT DILEPTON TYPE BINS ARE COMBINED BEFORE SCALE FACTOR IS APPLIED

  TString DilepYieldName("Dilep_");
  DilepYieldName+=binname;



  RooAddition DilepTopWJetsYield(DilepYieldName+"_TopWJetsYield",DilepYieldName+"_TopWJetsYield",
				 RooArgSet(twoTightMuTopWJetsYield, twoLooseLepTopWJetsYield));

  RooProduct  zeroLeptonTopWJetsYieldDilep("zeroLepton_"+binname+"_Dilep_TopWJetsYield","zeroLepton_"+binname+"_Dilep_TopWJetsYield",
					   RooArgSet( *wspace.arg(DilepScaleFactorName.Data()), DilepTopWJetsYield ));

  TString triggername = "oneLeptonTriggerEfficiency_";
  triggername.Append(trigeffname);
  RooAbsArg* triggerefficiency = wspace.arg(triggername.Data());

  RooProduct  zeroLeptonTopWJetsDataYieldDilep("zeroLepton_"+binname+"_Dilep_TopWJetsDataYield","zeroLepton_"+binname+"_Dilep_TopWJetsDataYield",
					       RooArgSet( zeroLeptonTopWJetsYieldDilep, *triggerefficiency ));

  wspace.import(zeroLeptonTopWJetsDataYieldDilep,RecycleConflictNodes());
  wspace.extendSet("nuisances",zeroLeptonTopWJetsDataYieldDilep.GetName());


  //cout << " INITIAL GUESS FOR DILEPTON PREDICTION: " << endl;
  //zeroLeptonTopWJetsYieldDilep.Print();
  //cout << endl << endl;



}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
// MAKE TAU->HAD AND DITAU PREDICTIONS FOR THIS NB/MET/HT BIN 

void makeTauHadBinPrediction( RooWorkspace& wspace, TString binname, TString binname_outside, TString trigeffname ){

  cout << " IN TAU PREDICTION" << endl;
  
  TString oneTightMuName("oneTightMu_");
  oneTightMuName+=binname;
  TString twoTightMuName("twoTightMu_");
  twoTightMuName+=binname;
  
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

  cout << " GETTING DILEPTON YIELDS " << endl;

  TString twoTightMuTopWJetsYieldName(twoTightMuName+"_TopWJetsYield");	
  cout << (*wspace.var(twoTightMuTopWJetsYieldName.Data())).getValV() << endl;;
  

  /////////////////////////////////////////////////////////////////////////////////////
  // SUM MU SAMPLES

  RooAddition oneTauInput( "oneTauInput_"+binname,"oneTauInput_"+binname,
			   RooArgSet(*wspace.var(oneTightMuTopWJetsYieldName1.Data()),*wspace.var(oneTightMuTopWJetsYieldName2.Data()),
				     *wspace.var(oneTightMuTopWJetsYieldName4.Data()),*wspace.var(oneTightMuTopWJetsYieldName3.Data()),
				     *wspace.var(oneTightMuTopWJetsYieldName5.Data())) );
  
  
  /////////////////////////////////////////////////////////////////////////////////////
  // GET COMMON MU->TAU SCALE FACTORS
    
  TString oneTauSFName("oneTau_");
  oneTauSFName+=binname_outside;
  oneTauSFName.Append("_ScaleFactor");
  TString twoTauSFName("twoTau_");
  twoTauSFName+=binname_outside;
  twoTauSFName.Append("_ScaleFactor");
  
  
  cout << " tauhad MC Scale Factors: " << endl;
  (*wspace.arg(oneTauSFName.Data())).Print();

  cout << " di tauhad MC Scale Factors: " << endl;
  (*wspace.arg(twoTauSFName.Data())).Print();
  

  cout << " APPLYING SCALE FACTORS TO SUMS " << endl;

  TString triggername = "oneLeptonTriggerEfficiency_";
  triggername.Append(trigeffname);
  RooAbsArg* triggerefficiency = wspace.arg(triggername.Data());

  RooProduct  zeroLeptonTopWJetsYield1Tau("zeroLepton_"+binname+"_1Tau_TopWJetsYield","zeroLepton_"+binname+"_1Tau_TopWJetsYield",
					  RooArgSet( *wspace.arg(oneTauSFName.Data()), oneTauInput ));

  RooProduct  zeroLeptonTopWJetsDataYield1Tau("zeroLepton_"+binname+"_1Tau_TopWJetsDataYield","zeroLepton_"+binname+"_1Tau_TopWJetsDataYield",
					  RooArgSet( zeroLeptonTopWJetsYield1Tau, *triggerefficiency ));

  RooProduct  zeroLeptonTopWJetsYield2Tau("zeroLepton_"+binname+"_2Tau_TopWJetsYield","zeroLepton_"+binname+"_2Tau_TopWJetsYield",
					  RooArgSet( *wspace.arg(twoTauSFName.Data()), *wspace.var(twoTightMuTopWJetsYieldName.Data()) ));
  
  RooProduct  zeroLeptonTopWJetsDataYield2Tau("zeroLepton_"+binname+"_2Tau_TopWJetsDataYield","zeroLepton_"+binname+"_2Tau_TopWJetsDataYield",
					      RooArgSet( zeroLeptonTopWJetsYield2Tau, *triggerefficiency ));
  

  /////////////////////////////////////////////////////////////////////////////////////
  // ADD TAU->HAD AND DITAU CONTRIBUTIONS TO 0L PREDICTION TO WORKSPACE
  
  cout << " IMPORTING PREDICTIONS " << endl;
  /*
  zeroLeptonTopWJetsYield1Tau.Print();
  zeroLeptonTopWJetsYield2Tau.Print();
  */

  wspace.import(zeroLeptonTopWJetsDataYield1Tau,RecycleConflictNodes());
  wspace.extendSet("nuisances",zeroLeptonTopWJetsDataYield1Tau.GetName());

  wspace.import(zeroLeptonTopWJetsDataYield2Tau,RecycleConflictNodes());
  wspace.extendSet("nuisances",zeroLeptonTopWJetsDataYield2Tau.GetName());

  
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


//void buildMRLikelihood( TString outputFile, TString setupFileName ) 
void buildMRLikelihood( RooWorkspace& wspace, TString outputFile, TString setupFileName, bool standalone ) 
{

  ///////////////////////////////////////////////////

  //int error = 0;
  //ofstream test;
  //test.open(outputFile.Data(),ios::app);

  
  TFile *test;
  if( standalone ) 
    {
      test = new TFile(outputFile.Data(), "UPDATE" );
      //  RooWorkspace *wspace = (RooWorkspace*) test->Get("wspace"); 
      wspace.autoImportClassCode(true);
    }

  RooRealVar *signalCrossSection = wspace.var( "signalCrossSection" );

  cout << " EVERYTHING IN INITIAL WORKSPACE: " << endl;
  wspace.Print();

  ///////////////////////////////////////////////////


  vector<TString> allbinnames;
  allbinnames.clear();

  // DEFINE BINS AND LOCATION OF BIN INPUTS
  vector<TString> binnames;
  map<TString,TString> binnamesoutside;
  map<TString,TString> binContentFileNames;
  map<TString,TString> binPolarizationScaleFactorFileNames;
  map<TString,TString> binTauHadScaleFactorFileNames;
  map<TString,TString> binSignalFracFileNames;
  map<TString,TString> binTriggerEfficiencyNames;


  ////////////////////////////////////
  ////// READ IN FILE WITH LIST OF BIN NAMES AND THE CORRESPONDING INPUT FILES
  ////// FOR CONTENT AND SCALE FACTORS (repeated scale factors will be specified there)

  ifstream setupFile;

  string fileLine;
  TString index, outsidename, fileName1, fileName2, fileName3, fileName4, trigeffname;

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
      listOfFiles.NextToken();
      trigeffname = listOfFiles;
 
      //cout << index << " : " << fileName << endl;
      if(index == "") continue;
    
      binnames.push_back(index);
      binContentFileNames[index] = fileName1;
      binnamesoutside[index] = outsidename;
      binPolarizationScaleFactorFileNames[index] = fileName2;
      binTauHadScaleFactorFileNames[index] = fileName3;
      binSignalFracFileNames[index] = fileName4;
      binTriggerEfficiencyNames[index] = trigeffname;
    }
  
  setupFile.close();
  

  ////////////////////////////////////
  ////// READ IN SIGNAL YIELDS FROM OUTSIDE, ADD TO WORKSPACE
  ////// OR, READ FROM WORKSPACE
  ////// IN STANDALONE MODE
  if( standalone ){
    
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

      TString oneLeptonName("oneLepton_");
      oneLeptonName+=thisBin;

      TString oneTightMuName("oneTightMu_");
      oneTightMuName+=thisBin;

      TString oneLooseLepName("oneLooseLep_");
      oneLooseLepName+=thisBin;

      TString twoTightMuName("twoTightMu_");
      twoTightMuName+=thisBin;
 
      TString twoLooseLepName("twoLooseLep_");
      twoLooseLepName+=thisBin;


      ///////////////////////////////
      ////// THIS IS A PLACEHOLDER!!!

      TString triggername = "oneLeptonTriggerEfficiency_";
      triggername.Append(trigeffname);

      RooRealVar* TriggerEfficiency = (RooRealVar*)
	getBetaPrimeConstraint(wspace,triggername, "",
			       0.95,0.05,
			       "trig1","trig2");

      ////// 

      signalfracFile>>count;

      RooRealVar acount0(zeroLeptonName+"_SignalFrac",zeroLeptonName+"_SignalFrac",count); 
      RooProduct zeroLepton(zeroLeptonName+"_SignalYield",zeroLeptonName+"_SignalYield",RooArgSet(acount0,*signalCrossSection));
      wspace.import(zeroLepton,RecycleConflictNodes());
      
      TString outputthis(zeroLeptonName);
      outputthis.Append("_SignalYield");
      //cout << endl << endl << " ZERO LEPTON SIGNAL YIELD MADE HERE AS " << outputthis.Data() << "  " << zeroLepton.getValV() << endl << endl;

      //////
      
      for( int j=1; j<6; j++ ){
	
	TString oneTightMuThetaName(oneTightMuName+"_Theta");
	oneTightMuThetaName+=j;

	TString oneLooseLepThetaName(oneLooseLepName+"_Theta");
	oneLooseLepThetaName+=j;

	signalfracFile>>count;
	RooRealVar acount1(oneTightMuThetaName+"_SignalFrac",oneTightMuThetaName+"_SignalFrac",count); 
	RooProduct oneTightMuTheta(oneTightMuThetaName+"_SignalYield",oneTightMuThetaName+"_SignalYield",RooArgSet(acount1,*signalCrossSection));
	wspace.import(oneTightMuTheta,RecycleConflictNodes());

	// LOOSE MU AND E ALREADY COMBINED 

	signalfracFile>>count;
	RooRealVar acount3(oneLooseLepThetaName+"_SignalFrac",oneLooseLepThetaName+"_SignalFrac",count); 
	RooProduct oneLooseLepTheta(oneLooseLepThetaName+"_SignalYield",oneLooseLepThetaName+"_SignalYield",RooArgSet(acount3,*signalCrossSection));
	wspace.import(oneLooseLepTheta,RecycleConflictNodes());
		
      }

      signalfracFile>>count;
      RooRealVar acount4(twoTightMuName+"_SignalFrac",twoTightMuName+"_SignalFrac",count); 
      RooProduct twoTightMu(twoTightMuName+"_SignalYield",twoTightMuName+"_SignalYield",RooArgSet(acount4,*signalCrossSection));
      wspace.import(twoTightMu,RecycleConflictNodes());

      signalfracFile>>count;
      RooRealVar acount9(twoLooseLepName+"_SignalFrac",twoLooseLepName+"_SignalFrac",count); 
      RooProduct twoLooseLep(twoLooseLepName+"_SignalYield",twoLooseLepName+"_SignalYield",RooArgSet(acount9,*signalCrossSection));
      wspace.import(twoLooseLep,RecycleConflictNodes());

          
    }// end of loop over bins setting 

  }// end standalone mode

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
      
      TString oneLooseLepName("oneLooseLep_");
      oneLooseLepName+=thisBin;

      TString twoTightMuName("twoTightMu_");
      twoTightMuName+=thisBin;
 
      TString twoLooseLepName("twoLooseLep_");
      twoLooseLepName+=thisBin;

      ////// 
      // leaving placeholder in 0L count space
      contentFile>>count;

      // ONLY ADD ZERO LEPTON COUNT TO WORKSPACE IN STANDALONE MODE
      if( standalone ){
	RooRealVar zeroLeptonCount(zeroLeptonName+"_Count",zeroLeptonName+"_Count",count);
	zeroLeptonCount.setConstant();
	wspace.import(zeroLeptonCount,RecycleConflictNodes());
	wspace.extendSet("namesfordata",zeroLeptonCount.GetName());
      }
      
      //////
      
      for( int j=1; j<6; j++ ){
	
	TString oneTightMuThetaName(oneTightMuName+"_Theta");
	oneTightMuThetaName+=j;

	TString oneLooseLepThetaName(oneLooseLepName+"_Theta");
	oneLooseLepThetaName+=j;
	
	contentFile>>count;     
 	RooRealVar oneTightMuThetaCount(oneTightMuThetaName+"_Count",oneTightMuThetaName+"_Count",count);
	oneTightMuThetaCount.setConstant();
	wspace.import(oneTightMuThetaCount,RecycleConflictNodes());
	wspace.extendSet("namesfordata",oneTightMuThetaCount.GetName());
	
	//cout << endl << endl << oneTightMuThetaCount.getVal() << endl << endl;

	contentFile>>count;
	RooRealVar oneLooseLepThetaCount(oneLooseLepThetaName+"_Count",oneLooseLepThetaName+"_Count",count);
	oneLooseLepThetaCount.setConstant();
	wspace.import(oneLooseLepThetaCount,RecycleConflictNodes());
	wspace.extendSet("namesfordata",oneLooseLepThetaCount.GetName());
	
	//cout << endl << endl << oneLooseLepThetaCount.getVal() << endl << endl;

      }
      

      contentFile>>count;
      RooRealVar twoTightMuCount(twoTightMuName+"_Count",twoTightMuName+"_Count",count);
      twoTightMuCount.setConstant();
      wspace.import(twoTightMuCount,RecycleConflictNodes());
      wspace.extendSet("namesfordata",twoTightMuCount.GetName());

      //cout << endl << endl << twoTightMuCount.getVal() << endl << endl;

      contentFile>>count;
      RooRealVar twoLooseLepCount(twoLooseLepName+"_Count",twoLooseLepName+"_Count",count);
      twoLooseLepCount.setConstant();
      wspace.import(twoLooseLepCount,RecycleConflictNodes());
      wspace.extendSet("namesfordata",twoLooseLepCount.GetName());

      //cout << endl << endl << twoLooseLepCount.getVal() << endl << endl;

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

	
	//RooRealVar ScaleFactor(oneLeptonName,oneLeptonName,scalefactor);
	//ScaleFactor.setConstant();
	//wspace.import(ScaleFactor,RecycleConflictNodes());
	RooRealVar* ScaleFactor = (RooRealVar*)
	  getBetaPrimeConstraint(wspace,oneLeptonName, "",
				 scalefactor,scalefactorerror,
				 name1,name2);
	
      }

      ////////////////////////////////////
      ////// READ IN DILEPTON METHOD SCALE FACTORS, ADD TO WORKSPACE

      double scalefactor1, scalefactorerror1, scalefactor2, scalefactorerror2, scalefactor3, scalefactorerror3;

      scaleFactorFile>>scalefactor>>scalefactorerror;
    
      TString dilepName("Dilep_");
      dilepName+=binnamesoutside[thisBin];
      dilepName.Append("_ScaleFactor");

      TString name1("blah3");
      TString name2("blah4");

      //RooRealVar ScaleFactorDilep(dilepName,dilepName,scalefactor);
      //ScaleFactorDilep.setConstant();
      //wspace.import(ScaleFactorDilep,RecycleConflictNodes());
      RooRealVar* ScaleFactorDilep = (RooRealVar*)
	getBetaPrimeConstraint(wspace,dilepName, "",
			       scalefactor,scalefactorerror,
			       name1,name2);
	      
      ////////////////////////////////////
      ////// READ IN TAU/DITAU MC-BASED SCALE FACTORS, ADD TO WORKSPACE

      tauhadScaleFactorFile>>scalefactor>>scalefactorerror;

      TString oneTauName("oneTau_");
      oneTauName+=binnamesoutside[thisBin];
      oneTauName.Append("_ScaleFactor");

      TString name3("blah5");
      TString name4("blah6");

    //RooRealVar ScaleFactorTau(oneTauName,oneTauName,scalefactor);
      //ScaleFactorTau.setConstant();
      //wspace.import(ScaleFactorTau,RecycleConflictNodes());
      RooRealVar* ScaleFactorTau = (RooRealVar*)
	getBetaPrimeConstraint(wspace,oneTauName, "",
			       scalefactor,scalefactorerror,
			       name3,name4);

      //////

      tauhadScaleFactorFile>>scalefactor>>scalefactorerror;
      
      TString twoTauName("twoTau_");
      twoTauName+=binnamesoutside[thisBin];
      twoTauName.Append("_ScaleFactor");

      TString name5("blah7");
      TString name6("blah8");

     //RooRealVar ScaleFactorTauTau(twoTauName,twoTauName,scalefactor);
      //ScaleFactorTauTau.setConstant();
      //wspace.import(ScaleFactorTauTau,RecycleConflictNodes());
      RooRealVar* ScaleFactorTauTau = (RooRealVar*)
	getBetaPrimeConstraint(wspace,twoTauName, "",
			       scalefactor,scalefactorerror,
			       name5,name6);

      //////

      scaleFactorFile.close();
      tauhadScaleFactorFile.close();

      ////////////////////////////////////
      ////////////////////////////////////
      ////// PUTTING SOME CONSTRAINT-SETTING CODE HERE FOR NOW:


      makePolarizationConstraintsPredictions( wspace, thisBin, binTriggerEfficiencyNames[thisBin] );
      //cout << " end of polarization constraints " << endl;
      makeDileptonConstraintsPredictions( wspace, thisBin, binnamesoutside[thisBin], binTriggerEfficiencyNames[thisBin] );
      //cout << " end of dilepton constraints " << endl;
      makeTauHadBinPrediction( wspace, thisBin, binnamesoutside[thisBin], binTriggerEfficiencyNames[thisBin] );
      //cout << " end of tauhad prediction " << endl;

      makePrediction( wspace, thisBin , standalone);


    }// end of loop over bins


  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  if( standalone ){
  
  RooArgSet Constraints(wspace.allPdfs());

  //  RooProdPdf model("model","model", Constraints );
    RooProdPdfLogSum model("model","model", Constraints );
  //  wspace.import(model,RecycleConflictNodes());
  wspace.import(model);

  cout << endl << endl << " FINAL SET OF CONSTRAINTS " << endl;
  Constraints.Print("v");

  RooDataSet dataset( "dataset", "dataset", *wspace.set("namesfordata") );
  dataset.add(*wspace.set("namesfordata"));
  wspace.import(dataset);
  

  cout << endl << endl << " FINAL SET OF DATA " << endl;
  dataset.Print();

  wspace.defineSet("poi","signalCrossSection");



  ////////////////////////////////////

  RooFitResult *fitResult = model.fitTo(dataset, Minos(kTRUE), Save(true),"s");
  
  /*

  RooAbsReal* nll = model.createNLL(dataset);
  RooAbsReal* pll = nll->createProfile(*wspace.set("poi")) ;

  RooPlot* frame1 = signalCrossSection->frame();
  nll->plotOn(frame1) ;
  pll->plotOn(frame1,LineColor(kRed)) ;

  TCanvas *c1 =  new TCanvas();
  frame1->Draw();
  c1->SaveAs("nlltest_1.eps");
  

  ////////////////////////////////////

 /*

  cout << " MAKING MODEL CONFIG " << endl;

  ModelConfig * modelConfig = new ModelConfig("modelConfig");
  modelConfig->SetWorkspace(*wspace);
  modelConfig->SetPdf("model");
  modelConfig->SetParametersOfInterest(*wspace.set("poi"));
  modelConfig->SetNuisanceParameters(*wspace.set("nuisances"));

  cout << endl << endl << endl << " before model print " << endl;
  modelConfig->Print();

  cout << " after model print " << endl;


  // Declare parameter of interest
  RooArgSet parofinterest(*wspace.set("poi")) ;
  ProfileLikelihoodCalculator plc(dataset, *wspace.pdf("model"), parofinterest);
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
  //  lplot->Save();
  canvas->Update();
  canvas->SaveAs("testnll5.eps");

  /*

  cout << " before simple fit " << endl;


  RooAbsReal* nll = model.createNLL(dataset);
  cout << "*********************************************" << endl;
  cout << "*************** PRINT NLL *******************" << endl;
  cout << "*********************************************" << endl;
  nll->Print("T");

  cout << "*********************************************" << endl;
  cout << "*************** PRINT Server0 *******************" << endl;
  cout << "*********************************************" << endl;
  nll->findServer(0)->Print();

  cout << "*********************************************" << endl;
  cout << "*************** PRINT Sever1 *******************" << endl;
  cout << "*********************************************" << endl;
  nll->findServer(1)->Print("V");

  */
  
  ////////////////////////////////////

 
  cout << " XSEC AFTER FIT " << wspace.var("signalCrossSection")->getValV() << "  " 
     << wspace.var("signalCrossSection")->getError() << endl;


  cout << wspace.function("zeroLepton_bin301_Theta1_TopWJetsYield")->getTitle() << "  " <<  wspace.function("zeroLepton_bin301_Theta1_TopWJetsYield")->getVal() << "  " <<  wspace.function("zeroLepton_bin301_Theta1_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin301_Theta2_TopWJetsYield")->getTitle() << "  " <<  wspace.function("zeroLepton_bin301_Theta2_TopWJetsYield")->getVal() << "  " <<  wspace.function("zeroLepton_bin301_Theta2_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin301_Theta3_TopWJetsYield")->getTitle() << "  " <<  wspace.function("zeroLepton_bin301_Theta3_TopWJetsYield")->getVal() << "  " <<  wspace.function("zeroLepton_bin301_Theta3_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin301_Theta4_TopWJetsYield")->getTitle() << "  " <<  wspace.function("zeroLepton_bin301_Theta4_TopWJetsYield")->getVal() << "  " <<  wspace.function("zeroLepton_bin301_Theta4_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin301_Theta5_TopWJetsYield")->getTitle() << "  " <<  wspace.function("zeroLepton_bin301_Theta5_TopWJetsYield")->getVal() << "  " <<  wspace.function("zeroLepton_bin301_Theta5_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin301_1Tau_TopWJetsYield")->getTitle() << "  " <<  wspace.function("zeroLepton_bin301_1Tau_TopWJetsYield")->getVal() << "  " <<  wspace.function("zeroLepton_bin301_1Tau_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin301_Dilep_TopWJetsYield")->getTitle() << "  " <<  wspace.function("zeroLepton_bin301_Dilep_TopWJetsYield")->getVal() << "  " <<  wspace.function("zeroLepton_bin301_Dilep_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin301_2Tau_TopWJetsYield")->getTitle() << "  " <<  wspace.function("zeroLepton_bin301_2Tau_TopWJetsYield")->getVal() << "  " <<  wspace.function("zeroLepton_bin301_2Tau_TopWJetsYield")->getPropagatedError(*fitResult) << endl;

  cout << endl;
  cout << wspace.function("zeroLepton_bin333_Theta1_TopWJetsYield")->getTitle() << "  " << wspace.function("zeroLepton_bin333_Theta1_TopWJetsYield")->getVal() << "  " << wspace.function("zeroLepton_bin333_Theta1_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin333_Theta2_TopWJetsYield")->getTitle() << "  " << wspace.function("zeroLepton_bin333_Theta2_TopWJetsYield")->getVal() << "  " << wspace.function("zeroLepton_bin333_Theta2_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin333_Theta3_TopWJetsYield")->getTitle() << "  " << wspace.function("zeroLepton_bin333_Theta3_TopWJetsYield")->getVal() << "  " << wspace.function("zeroLepton_bin333_Theta3_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin333_Theta4_TopWJetsYield")->getTitle() << "  " << wspace.function("zeroLepton_bin333_Theta4_TopWJetsYield")->getVal() << "  " << wspace.function("zeroLepton_bin333_Theta4_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin333_Theta5_TopWJetsYield")->getTitle() << "  " << wspace.function("zeroLepton_bin333_Theta5_TopWJetsYield")->getVal() << "  " << wspace.function("zeroLepton_bin333_Theta5_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin333_1Tau_TopWJetsYield")->getTitle() << "  " << wspace.function("zeroLepton_bin333_1Tau_TopWJetsYield")->getVal() << "  " << wspace.function("zeroLepton_bin333_1Tau_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin333_Dilep_TopWJetsYield")->getTitle() << "  " << wspace.function("zeroLepton_bin333_Dilep_TopWJetsYield")->getVal() << "  " << wspace.function("zeroLepton_bin333_Dilep_TopWJetsYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("zeroLepton_bin333_2Tau_TopWJetsYield")->getTitle() << "  " << wspace.function("zeroLepton_bin333_2Tau_TopWJetsYield")->getVal() << "  " << wspace.function("zeroLepton_bin333_2Tau_TopWJetsYield")->getPropagatedError(*fitResult) << endl;

  cout << endl << endl;


  cout << wspace.function("zeroLepton_bin301_SignalYield")->getTitle() << "  " << wspace.function("zeroLepton_bin301_SignalYield")->getVal() << "  " << wspace.function("zeroLepton_bin301_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneLooseLep_bin301_Theta1_SignalYield")->getTitle() << "  " << wspace.function("oneLooseLep_bin301_Theta1_SignalYield")->getVal() << "  " << wspace.function("oneLooseLep_bin301_Theta1_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneLooseLep_bin301_Theta2_SignalYield")->getTitle() << "  " << wspace.function("oneLooseLep_bin301_Theta2_SignalYield")->getVal() << "  " << wspace.function("oneLooseLep_bin301_Theta2_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneLooseLep_bin301_Theta3_SignalYield")->getTitle() << "  " << wspace.function("oneLooseLep_bin301_Theta3_SignalYield")->getVal() << "  " << wspace.function("oneLooseLep_bin301_Theta3_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneLooseLep_bin301_Theta4_SignalYield")->getTitle() << "  " << wspace.function("oneLooseLep_bin301_Theta4_SignalYield")->getVal() << "  " << wspace.function("oneLooseLep_bin301_Theta4_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneLooseLep_bin301_Theta5_SignalYield")->getTitle() << "  " <<   wspace.function("oneLooseLep_bin301_Theta5_SignalYield")->getVal() << "  " <<  wspace.function("oneTightMu_bin301_Theta5_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneTightMu_bin301_Theta1_SignalYield")->getTitle() << "  " << wspace.function("oneTightMu_bin301_Theta1_SignalYield")->getVal() << "  " << wspace.function("oneTightMu_bin301_Theta1_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneTightMu_bin301_Theta2_SignalYield")->getTitle() << "  " << wspace.function("oneTightMu_bin301_Theta2_SignalYield")->getVal() << "  " << wspace.function("oneTightMu_bin301_Theta2_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneTightMu_bin301_Theta3_SignalYield")->getTitle() << "  " << wspace.function("oneTightMu_bin301_Theta3_SignalYield")->getVal() << "  " << wspace.function("oneTightMu_bin301_Theta3_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneTightMu_bin301_Theta4_SignalYield")->getTitle() << "  " << wspace.function("oneTightMu_bin301_Theta4_SignalYield")->getVal() << "  " << wspace.function("oneTightMu_bin301_Theta4_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneTightMu_bin301_Theta5_SignalYield")->getTitle() << "  " <<   wspace.function("oneTightMu_bin301_Theta5_SignalYield")->getVal() << "  " <<  wspace.function("oneTightMu_bin301_Theta5_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("twoLooseLep_bin301_SignalYield")->getTitle() << "  " << wspace.function("twoLooseLep_bin301_SignalYield")->getVal() << "  " << wspace.function("twoLooseLep_bin301_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("twoTightMu_bin301_SignalYield")->getTitle() << "  " << wspace.function("twoTightMu_bin301_SignalYield")->getVal() << "  " << wspace.function("twoTightMu_bin301_SignalYield")->getPropagatedError(*fitResult) << endl;

  cout << endl;


  cout << wspace.function("zeroLepton_bin333_SignalYield")->getTitle() << "  " << wspace.function("zeroLepton_bin333_SignalYield")->getVal() << "  " << wspace.function("zeroLepton_bin333_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneLooseLep_bin333_Theta1_SignalYield")->getTitle() << "  " << wspace.function("oneLooseLep_bin333_Theta1_SignalYield")->getVal() << "  " << wspace.function("oneLooseLep_bin333_Theta1_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneLooseLep_bin333_Theta2_SignalYield")->getTitle() << "  " << wspace.function("oneLooseLep_bin333_Theta2_SignalYield")->getVal() << "  " << wspace.function("oneLooseLep_bin333_Theta2_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneLooseLep_bin333_Theta3_SignalYield")->getTitle() << "  " << wspace.function("oneLooseLep_bin333_Theta3_SignalYield")->getVal() << "  " << wspace.function("oneLooseLep_bin333_Theta3_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneLooseLep_bin333_Theta4_SignalYield")->getTitle() << "  " << wspace.function("oneLooseLep_bin333_Theta4_SignalYield")->getVal() << "  " << wspace.function("oneLooseLep_bin333_Theta4_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneLooseLep_bin333_Theta5_SignalYield")->getTitle() << "  " << wspace.function("oneLooseLep_bin333_Theta5_SignalYield")->getVal() << "  " << wspace.function("oneLooseLep_bin333_Theta5_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneTightMu_bin333_Theta1_SignalYield")->getTitle() << "  " << wspace.function("oneTightMu_bin333_Theta1_SignalYield")->getVal() << "  " << wspace.function("oneTightMu_bin333_Theta1_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneTightMu_bin333_Theta2_SignalYield")->getTitle() << "  " << wspace.function("oneTightMu_bin333_Theta2_SignalYield")->getVal() << "  " << wspace.function("oneTightMu_bin333_Theta2_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneTightMu_bin333_Theta3_SignalYield")->getTitle() << "  " << wspace.function("oneTightMu_bin333_Theta3_SignalYield")->getVal() << "  " << wspace.function("oneTightMu_bin333_Theta3_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneTightMu_bin333_Theta4_SignalYield")->getTitle() << "  " << wspace.function("oneTightMu_bin333_Theta4_SignalYield")->getVal() << "  " << wspace.function("oneTightMu_bin333_Theta4_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("oneTightMu_bin333_Theta5_SignalYield")->getTitle() << "  " << wspace.function("oneTightMu_bin333_Theta5_SignalYield")->getVal() << "  " << wspace.function("oneTightMu_bin333_Theta5_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("twoLooseLep_bin333_SignalYield")->getTitle() << "  " << wspace.function("twoLooseLep_bin333_SignalYield")->getVal() << "  " << wspace.function("twoLooseLep_bin333_SignalYield")->getPropagatedError(*fitResult) << endl;
  cout << wspace.function("twoTightMu_bin333_SignalYield")->getTitle() << "  " << wspace.function("twoTightMu_bin333_SignalYield")->getVal() << "  " << wspace.function("twoTightMu_bin333_SignalYield")->getPropagatedError(*fitResult) << endl;    


  //test << wspace.var("signalCrossSection")->getValV() << endl;
  test->Write(); 
  test->Close(); 

  }

  //wspace.writeToFile( outputFile.Data(), false ); // DO NOT RECREATE ROOT FILE!!! 

  ////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	




}// end of likelihood builder





