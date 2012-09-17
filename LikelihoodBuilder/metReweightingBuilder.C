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

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "TMath.h"

using namespace RooFit ;
using namespace RooStats ;



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////



void makeBigBin( RooWorkspace& wspace, TString binname, ifstream &scalefactorsFile, ifstream &sigfracFile ){

  ////// GET STUFF YOU NEED FROM WORKSPACE

  RooRealVar* signalCrossSection = wspace.var("signalCrossSection");

  TString binname1;
  scalefactorsFile>>binname1;
  
  // cout << endl << " STARTING BIG BIN... " << endl << endl;

  //  cout << binname << "  " << binname1 << endl;

  if( binname!=binname1 ){
    cout << endl << endl << " ERROR IN BIN NAMES!!!! " << endl;
    cout << binname << "  " << binname1 << endl;

    //error = 1;
  }
  
  //  cout << " before loop " << endl;
  //////////////////////////////////////////////////////////////////
  // LOOP OVER STUFF IN DTHETA BINS!!!
  //////////////////////////////////////////////////////////////////

  for( int dtheta=1; dtheta<6; dtheta++ ){

    ////////////////////////////////////////////////////
    // GET CONTROL SAMPLE COUNTS FOR THIS DTHETA BIN 

    TString oneLeptonName("oneLepton_");
    oneLeptonName+=binname;
    TString zeroLeptonName("zeroLepton_");
    zeroLeptonName+=binname;
 
    TString oneTightMuName("oneTightMu_");
    oneTightMuName+=binname;
    TString oneLooseMuName("oneLooseMu_");
    oneLooseMuName+=binname;
    TString oneLooseEName("oneLooseE_");
    oneLooseEName+=binname;

    TString ending("_Theta");
    ending+=dtheta;
    ending.Append("_Count");

    TString oneTightMuCountName(oneTightMuName);
    oneTightMuCountName.Append(ending);
    TString oneLooseMuCountName(oneLooseMuName);
    oneLooseMuCountName.Append(ending);
    TString oneLooseECountName(oneLooseEName);
    oneLooseECountName.Append(ending);

    //    cout << oneTightMuCountName.Data() << "  " << oneLooseMuCountName.Data() << "  " << oneLooseECountName.Data() << endl;

    RooRealVar oneTightMuCount = *wspace.var(oneTightMuCountName.Data());
    RooRealVar oneLooseMuCount = *wspace.var(oneLooseMuCountName.Data());
    RooRealVar oneLooseECount = *wspace.var(oneLooseECountName.Data());

    //    cout << ending << "  " << oneTightMuCount.getVal() << "  " << oneLooseMuCount.getVal() << "  " << oneLooseECount.getVal() << "  " << endl;

    ////////////////////////////////////////////////////
    // GET SIGNAL FRACTIONS FOR THIS DTHETA BIN 
 
    double oneTightMusignalfrac, oneLooseMusignalfrac, oneLooseESsgnalfrac, scalefactor;

    sigfracFile>>oneTightMusignalfrac;
    sigfracFile>>oneLooseMusignalfrac;
    sigfracFile>>oneLooseESsgnalfrac;
    
    RooRealVar oneTightMuSignalFrac(oneTightMuName+"_SignalFrac",oneTightMuName+"_SignalFrac",oneTightMusignalfrac,0.000001,1000.);
    oneTightMuSignalFrac.setConstant();
    RooRealVar oneLooseMuSignalFrac(oneLooseMuName+"_SignalFrac",oneLooseMuName+"_SignalFrac",oneLooseMusignalfrac,0.000001,1000.);
    oneLooseMuSignalFrac.setConstant();
    RooRealVar oneLooseESignalFrac(oneLooseEName+"_SignalFrac",oneLooseEName+"_SignalFrac",oneLooseESsgnalfrac,0.000001,1000.);
    oneLooseESignalFrac.setConstant();
    
    ////////////////////////////////////////////////////
    // GET 1L->0L SCALE FACTORS FOR THIS DTHETA BIN 
    
    scalefactorsFile>>scalefactor;
    RooRealVar ScaleFactor("zeroLepton_"+binname+"_ScaleFactor","zeroLepton_"+binname+"_ScaleFactor",scalefactor,0.000001,100000.);
    ScaleFactor.setConstant();
    
    ////////////////////////////////////////////////////
    // INITIALIZE TTBAR COMPONENT OF THIS DTHETA BIN


    TString ThetaName("_TopWJetsYield_Theta");
    ThetaName+=dtheta;

    RooRealVar oneTightMuTopWJetsYield(oneTightMuName+ThetaName,oneTightMuName+ThetaName,oneTightMuCount.getVal(),0.000001,2000.);
    wspace.import(oneTightMuTopWJetsYield,RecycleConflictNodes());
    wspace.extendSet("nuisances",oneTightMuTopWJetsYield.GetName());

    RooRealVar oneLooseMuTopWJetsYield(oneLooseMuName+ThetaName,oneLooseMuName+ThetaName,oneLooseMuCount.getVal(),0.000001,2000.);
    wspace.import(oneLooseMuTopWJetsYield,RecycleConflictNodes());
    wspace.extendSet("nuisances",oneLooseMuTopWJetsYield.GetName());

    RooRealVar oneLooseETopWJetsYield(oneLooseEName+ThetaName,oneLooseEName+ThetaName,oneLooseECount.getVal(),0.000001,2000.);
    wspace.import(oneLooseETopWJetsYield,RecycleConflictNodes());
    wspace.extendSet("nuisances",oneLooseETopWJetsYield.GetName());


    ////////////////////////////////////////////////////
    // SET UP SIGNAL COMPONENT OF THIS DTHETA BIN AND CONTROL SAMPLE CONSTRAINTS
    
    /****** NOTE: INSERT SINGLE LEPTON EFFICIENCIES HERE??? ******/
    TString SignalYieldThetaName("_SignalYield_Theta");
    SignalYieldThetaName+=dtheta;
    TString YieldSumThetaName("_YieldSum_Theta");
    YieldSumThetaName+=dtheta;
    TString ConstraintThetaName("_Constraint_Theta");
    ConstraintThetaName+=dtheta;


    RooProduct oneTightMuSignalYield(oneTightMuName+SignalYieldThetaName,oneTightMuName+SignalYieldThetaName,RooArgSet(*signalCrossSection,oneTightMuSignalFrac));
    RooProduct oneLooseMuSignalYield(oneLooseMuName+SignalYieldThetaName,oneLooseMuName+SignalYieldThetaName,RooArgSet(*signalCrossSection,oneLooseMuSignalFrac));
    RooProduct oneLooseESignalYield(oneLooseEName+SignalYieldThetaName,oneLooseEName+SignalYieldThetaName,RooArgSet(*signalCrossSection,oneLooseESignalFrac));
    
    RooAddition oneTightMuYieldSum(oneTightMuName+YieldSumThetaName,oneTightMuName+YieldSumThetaName,RooArgSet(oneTightMuSignalYield,oneTightMuTopWJetsYield));
    RooAddition oneLooseMuYieldSum(oneLooseMuName+YieldSumThetaName,oneLooseMuName+YieldSumThetaName,RooArgSet(oneLooseMuSignalYield,oneLooseMuTopWJetsYield));
    RooAddition oneLooseEYieldSum(oneLooseEName+YieldSumThetaName,oneLooseEName+YieldSumThetaName,RooArgSet(oneLooseESignalYield,oneLooseETopWJetsYield));
   
    RooPoisson oneTightMuConstraint(oneTightMuName+ConstraintThetaName,oneTightMuName+ConstraintThetaName,oneTightMuCount,oneTightMuYieldSum);
    wspace.import( oneTightMuConstraint,RecycleConflictNodes() );
    RooPoisson oneLooseMuConstraint(oneLooseMuName+ConstraintThetaName,oneLooseMuName+ConstraintThetaName,oneLooseMuCount,oneLooseMuYieldSum);
    wspace.import( oneLooseMuConstraint,RecycleConflictNodes() );
    RooPoisson oneLooseEConstraint(oneLooseEName+ConstraintThetaName,oneLooseEName+ConstraintThetaName,oneLooseECount,oneLooseEYieldSum);
    wspace.import( oneLooseEConstraint,RecycleConflictNodes() );


    ////////////////////////////////////////////////////
    // CONSTRUCT SOME 0L BACKGROUND PREDICTION FROM THIS DTHETA BIN
 
    RooAddition oneLeptonTopWJetsYield(oneLeptonName+ThetaName,oneLeptonName+ThetaName,
				       RooArgSet(oneTightMuTopWJetsYield,oneLooseMuTopWJetsYield,oneLooseETopWJetsYield));
    RooProduct  zeroLeptonTopWJetsYield(zeroLeptonName+ThetaName,zeroLeptonName+ThetaName,
					RooArgSet(ScaleFactor,oneLeptonTopWJetsYield));
    wspace .import(zeroLeptonTopWJetsYield,RecycleConflictNodes());

  }// END OF LOOP OVER DTHETA BINS
  
}// END OF BIG BIN IN NBS/HT/MET


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

void makeDileptonConstraintsPrediction( RooWorkspace& wspace, TString binname, ifstream &scaleFactorsFile, ifstream &sigfracFile ){

  /////////////////////////////////////////////////////////////////////////////////////
  // GET DILEPTON COUNTS
  
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
 

  /////////////////////////////////////////////////////////////////////////////////////
  // GET DILEPTON SIGNAL FRACTIONS
  
   RooRealVar* signalCrossSection = wspace.var("signalCrossSection");

   double sfactorMuMu, sfactorEMu, sfactorEE, twoTightMusigfrac,  twoTightMuLooseMusigfrac, twoLooseMusigfrac, twoTightMuLooseEsigfrac, twoLooseMuLooseEsigfrac, twoLooseEsigfrac;

   scaleFactorsFile>>sfactorMuMu;
   scaleFactorsFile>>sfactorEMu;
   scaleFactorsFile>>sfactorEE;
     
   sigfracFile>>twoTightMusigfrac;
   sigfracFile>>twoTightMuLooseMusigfrac;
   sigfracFile>>twoLooseMusigfrac;
   sigfracFile>>twoTightMuLooseEsigfrac;
   sigfracFile>>twoLooseMuLooseEsigfrac;
   sigfracFile>>twoLooseEsigfrac;


   RooRealVar ScaleFactorMuMu(twoLooseMuName+"_ScaleFactor",twoLooseMuName+"_ScaleFactor",sfactorMuMu,0.000001,100000.);
   ScaleFactorMuMu.setConstant();
   RooRealVar ScaleFactorEMu(twoLooseMuLooseEName+"_ScaleFactor",twoLooseMuLooseEName+"_ScaleFactor",sfactorEMu,0.000001,100000.);
   ScaleFactorEMu.setConstant();
   RooRealVar ScaleFactorEE(twoLooseEName+"_ScaleFactor",twoLooseEName+"_ScaleFactor",sfactorEE,0.000001,100000.);
   ScaleFactorEE.setConstant();

    RooRealVar twoTightMuSignalFrac(twoTightMuName+"_SignalFraction",twoTightMuName+"_SignalFraction",twoTightMusigfrac,0.000001,100000.);
    twoTightMuSignalFrac.setConstant();
    RooRealVar twoTightMuLooseMuSignalFrac(twoTightMuLooseMuName+"_SignalFraction",twoTightMuLooseMuName+"_SignalFraction",twoTightMuLooseMusigfrac,0.000001,100000.);
    twoTightMuLooseMuSignalFrac.setConstant();
    RooRealVar twoLooseMuSignalFrac(twoLooseMuName+"_SignalFraction",twoLooseMuName+"_SignalFraction",twoLooseMusigfrac,0.000001,100000.);
    twoLooseMuSignalFrac.setConstant();
    RooRealVar twoTightMuLooseESignalFrac(twoTightMuLooseEName+"_SignalFraction",twoTightMuLooseEName+"_SignalFraction",twoTightMuLooseEsigfrac,0.000001,100000.);
    twoTightMuLooseESignalFrac.setConstant();
    RooRealVar twoLooseMuLooseESignalFrac(twoLooseMuLooseEName+"_SignalFraction",twoLooseMuLooseEName+"_SignalFraction",twoLooseMuLooseEsigfrac,0.000001,100000.);
    twoLooseMuLooseESignalFrac.setConstant();
    RooRealVar twoLooseESignalFrac(twoLooseEName+"_SignalFraction",twoLooseEName+"_SignalFraction",twoLooseEsigfrac,0.000001,100000.);
    twoLooseESignalFrac.setConstant();

  /////////////////////////////////////////////////////////////////////////////////////
  // IMPORT BACKGROUND CONTRIBUTIONS TO DILEPTON SAMPLES TO WORKSPACE
  
  RooRealVar twoTightMuTopWJetsYield(twoTightMuName+"_TopWJetsYield",twoTightMuName+"_TopWJetsYield",twoTightMuCount.getVal(),0.000001,2000.);
  wspace.import(twoTightMuTopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoTightMuTopWJetsYield.GetName());

  RooRealVar twoTightMuLooseMuTopWJetsYield(twoTightMuLooseMuName+"_TopWJetsYield",twoTightMuLooseMuName+"_TopWJetsYield",twoTightMuLooseMuCount.getVal(),0.000001,2000.);
  wspace.import(twoTightMuLooseMuTopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoTightMuLooseMuTopWJetsYield.GetName());

  RooRealVar twoLooseMuTopWJetsYield(twoLooseMuName+"_TopWJetsYield",twoLooseMuName+"_TopWJetsYield",twoLooseMuCount.getVal(),0.000001,2000.);
  wspace.import(twoLooseMuTopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoLooseMuTopWJetsYield.GetName());

  RooRealVar twoTightMuLooseETopWJetsYield(twoTightMuLooseEName+"_TopWJetsYield",twoTightMuLooseEName+"_TopWJetsYield",twoTightMuLooseECount.getVal(),0.000001,2000.);
  wspace.import(twoTightMuLooseETopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoTightMuLooseETopWJetsYield.GetName());

  RooRealVar twoLooseMuLooseETopWJetsYield(twoLooseMuLooseEName+"_TopWJetsYield",twoLooseMuLooseEName+"_TopWJetsYield",twoLooseMuLooseECount.getVal(),0.000001,2000.);
  wspace.import(twoLooseMuLooseETopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoLooseMuLooseETopWJetsYield.GetName());

  RooRealVar twoLooseETopWJetsYield(twoLooseEName+"_TopWJetsYield",twoLooseEName+"_TopWJetsYield",twoLooseECount.getVal(),0.000001,2000.);
  wspace.import(twoLooseETopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoLooseETopWJetsYield.GetName());


  /////////////////////////////////////////////////////////////////////////////////////
  // PRODUCE SIGNAL+BACKGROUND DILEPTON CONSTRAINTS (NOT EVER BINNED IN DTHETA) 
 
  RooProduct twoTightMuSignalYield(twoTightMuName+"_SignalYield",twoTightMuName+"_SignalYield",RooArgSet(*signalCrossSection,twoTightMuSignalFrac));
  RooProduct twoTightMuLooseMuSignalYield(twoTightMuLooseMuName+"_SignalYield",twoTightMuLooseMuName+"_SignalYield",RooArgSet(*signalCrossSection,twoTightMuLooseMuSignalFrac));
  RooProduct twoLooseMuSignalYield(twoLooseMuName+"_SignalYield",twoLooseMuName+"_SignalYield",RooArgSet(*signalCrossSection,twoLooseMuSignalFrac));
  RooProduct twoTightMuLooseESignalYield(twoTightMuLooseEName+"_SignalYield",twoTightMuLooseEName+"_SignalYield",RooArgSet(*signalCrossSection,twoTightMuLooseESignalFrac));
  RooProduct twoLooseMuLooseESignalYield(twoLooseMuLooseEName+"_SignalYield",twoLooseMuLooseEName+"_SignalYield",RooArgSet(*signalCrossSection,twoLooseMuLooseESignalFrac));
  RooProduct twoLooseESignalYield(twoLooseEName+"_SignalYield",twoLooseEName+"_SignalYield",RooArgSet(*signalCrossSection,twoLooseESignalFrac));
 
  RooAddition twoTightMuYieldSum(twoTightMuName+"_YieldSum",twoTightMuName+"_YieldSum",RooArgSet(twoTightMuSignalYield,twoTightMuTopWJetsYield));
  RooAddition twoTightMuLooseMuYieldSum(twoTightMuLooseMuName+"_YieldSum",twoTightMuLooseMuName+"_YieldSum",RooArgSet(twoTightMuLooseMuSignalYield,twoTightMuLooseMuTopWJetsYield));
  RooAddition twoLooseMuYieldSum(twoLooseMuName+"_YieldSum",twoLooseMuName+"_YieldSum",RooArgSet(twoLooseMuSignalYield,twoLooseMuTopWJetsYield));
  RooAddition twoTightMuLooseEYieldSum(twoTightMuLooseEName+"_YieldSum",twoTightMuLooseEName+"_YieldSum",RooArgSet(twoTightMuLooseESignalYield,twoTightMuLooseETopWJetsYield));
  RooAddition twoLooseMuLooseEYieldSum(twoLooseMuLooseEName+"_YieldSum",twoLooseMuLooseEName+"_YieldSum",RooArgSet(twoLooseMuLooseESignalYield,twoLooseMuLooseETopWJetsYield));
  RooAddition twoLooseEYieldSum(twoLooseEName+"_YieldSum",twoLooseEName+"_YieldSum",RooArgSet(twoLooseESignalYield,twoLooseETopWJetsYield));
  
  /*

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

  */

  /////////////////////////////////////////////////////////////////////////////////////
  // PRODUCE 0L PREDICTIONS FOR MUMU, EMU, EE DILEPTONS (NOT EVER BINNED IN DTHETA) 

  RooProduct  zeroLeptonTopWJetsYieldMuMu("zeroLepton_"+binname+"_TopWJetsYield_MuMu","zeroLepton_"+binname+"_TopWJetsYield_MuMu",
					  RooArgSet( ScaleFactorMuMu, twoLooseMuTopWJetsYield ));
  wspace.import(zeroLeptonTopWJetsYieldMuMu,RecycleConflictNodes());
 
  RooProduct  zeroLeptonTopWJetsYieldEMu("zeroLepton_"+binname+"_TopWJetsYield_EMu","zeroLepton_"+binname+"_TopWJetsYield_EMu",
					  RooArgSet( ScaleFactorEMu, twoLooseMuLooseETopWJetsYield ));
  wspace.import(zeroLeptonTopWJetsYieldEMu,RecycleConflictNodes());
 
  RooProduct  zeroLeptonTopWJetsYieldEE("zeroLepton_"+binname+"_TopWJetsYield_EE","zeroLepton_"+binname+"_TopWJetsYield_EE",
					  RooArgSet( ScaleFactorEE, twoLooseETopWJetsYield ));
  wspace.import(zeroLeptonTopWJetsYieldEE,RecycleConflictNodes());


}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

void makeOutsideUnbinnedConstraints( RooWorkspace& wspace, int nbs, double oneTightMusigfrac, double twoTightMusigfrac, double twoTightMuLooseMusigfrac, double twoTightMuLooseEsigfrac ){

  /////////////////////////////////////////////////////////////////////////////////////
  // GET TAU AND DITAU COUNTS ONLY, FOR EVENTS *OUTSIDE* MET/HT/ETC CUTS
  
  TString oneTightMuName("oneTightMu_");
  oneTightMuName+=nbs;
  oneTightMuName.Append("Outside");
  TString twoTightMuName("twoTightMu_");
  twoTightMuName+=nbs;
  twoTightMuName.Append("Outside");
  TString twoTightMuLooseMuName("twoTightMuLooseMu_");
  twoTightMuLooseMuName+=nbs;
  twoTightMuLooseMuName.Append("Outside");
  TString twoTightMuLooseEName("twoTightMuLooseE_");
  twoTightMuLooseEName+=nbs;
  twoTightMuLooseEName.Append("Outside");
  
  TString oneTightMuCountName(oneTightMuName);
  oneTightMuCountName.Append("_Count");
  TString twoTightMuCountName(twoTightMuName);
  twoTightMuCountName.Append("_Count");
  TString twoTightMuLooseMuCountName(twoTightMuLooseMuName);
  twoTightMuLooseMuCountName.Append("_Count");
  TString twoTightMuLooseECountName(twoTightMuLooseEName);
  twoTightMuLooseECountName.Append("_Count");

  RooRealVar oneTightMuCount = *wspace.var(oneTightMuCountName.Data());
  RooRealVar twoTightMuCount = *wspace.var(twoTightMuCountName.Data());
  RooRealVar twoTightMuLooseMuCount = *wspace.var(twoTightMuLooseMuCountName.Data());
  RooRealVar twoTightMuLooseECount = *wspace.var(twoTightMuLooseECountName.Data());
  
  /////////////////////////////////////////////////////////////////////////////////////
  // GET DILEPTON SIGNAL FRACTIONS
  
  RooRealVar* signalCrossSection = wspace.var("signalCrossSection");

  RooRealVar oneTightMuSignalFrac(oneTightMuName+"_SignalFraction",oneTightMuName+"_SignalFraction",oneTightMusigfrac,0.000001,100000.);
  oneTightMuSignalFrac.setConstant();

  RooRealVar twoTightMuSignalFrac(twoTightMuName+"_SignalFraction",twoTightMuName+"_SignalFraction",twoTightMusigfrac,0.000001,100000.);
  twoTightMuSignalFrac.setConstant();

  RooRealVar twoTightMuLooseMuSignalFrac(twoTightMuLooseMuName+"_SignalFraction",twoTightMuLooseMuName+"_SignalFraction",twoTightMuLooseMusigfrac,0.000001,100000.);
  twoTightMuLooseMuSignalFrac.setConstant();

  RooRealVar twoTightMuLooseESignalFrac(twoTightMuLooseEName+"_SignalFraction",twoTightMuLooseEName+"_SignalFraction",twoTightMuLooseEsigfrac,0.000001,100000.);
  twoTightMuLooseESignalFrac.setConstant();

  /////////////////////////////////////////////////////////////////////////////////////
  // IMPORT BACKGROUND CONTRIBUTIONS TO TAU AND DITAU SAMPLES TO WORKSPACE
 
  RooRealVar oneTightMuTopWJetsYield(oneTightMuName+"_TopWJetsYield",oneTightMuName+"_TopWJetsYield",oneTightMuCount.getVal(),0.000001,2000.);
  wspace.import(oneTightMuTopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",oneTightMuTopWJetsYield.GetName());
  
  RooRealVar twoTightMuTopWJetsYield(twoTightMuName+"_TopWJetsYield",twoTightMuName+"_TopWJetsYield",twoTightMuCount.getVal(),0.000001,2000.);
  wspace.import(twoTightMuTopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoTightMuTopWJetsYield.GetName());

  RooRealVar twoTightMuLooseMuTopWJetsYield(twoTightMuLooseMuName+"_TopWJetsYield",twoTightMuLooseMuName+"_TopWJetsYield",twoTightMuLooseMuCount.getVal(),0.000001,2000.);
  wspace.import(twoTightMuLooseMuTopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoTightMuLooseMuTopWJetsYield.GetName());

  RooRealVar twoTightMuLooseETopWJetsYield(twoTightMuLooseEName+"_TopWJetsYield",twoTightMuLooseEName+"_TopWJetsYield",twoTightMuLooseECount.getVal(),0.000001,2000.);
  wspace.import(twoTightMuLooseETopWJetsYield,RecycleConflictNodes());
  wspace.extendSet("nuisances",twoTightMuLooseETopWJetsYield.GetName());


  /////////////////////////////////////////////////////////////////////////////////////
  // PRODUCE SIGNAL+BACKGROUND DILEPTON CONSTRAINTS OUTSIDE HT/MET BINS 
  
  RooProduct oneTightMuSignalYield(oneTightMuName+"_SignalYield",oneTightMuName+"_SignalYield",RooArgSet(*signalCrossSection,oneTightMuSignalFrac));
  RooProduct twoTightMuSignalYield(twoTightMuName+"_SignalYield",twoTightMuName+"_SignalYield",RooArgSet(*signalCrossSection,twoTightMuSignalFrac));
  RooProduct twoTightMuLooseMuSignalYield(twoTightMuLooseMuName+"_SignalYield",twoTightMuLooseMuName+"_SignalYield",RooArgSet(*signalCrossSection,twoTightMuLooseMuSignalFrac));
  RooProduct twoTightMuLooseESignalYield(twoTightMuLooseEName+"_SignalYield",twoTightMuLooseEName+"_SignalYield",RooArgSet(*signalCrossSection,twoTightMuLooseESignalFrac));
  
  RooAddition oneTightMuYieldSum(oneTightMuName+"_YieldSum",oneTightMuName+"_YieldSum",RooArgSet(oneTightMuSignalYield,oneTightMuTopWJetsYield));
  RooAddition twoTightMuYieldSum(twoTightMuName+"_YieldSum",twoTightMuName+"_YieldSum",RooArgSet(twoTightMuSignalYield,twoTightMuTopWJetsYield));
  RooAddition twoTightMuLooseMuYieldSum(twoTightMuLooseMuName+"_YieldSum",twoTightMuLooseMuName+"_YieldSum",RooArgSet(twoTightMuLooseMuSignalYield,twoTightMuLooseMuTopWJetsYield));
  RooAddition twoTightMuLooseEYieldSum(twoTightMuLooseEName+"_YieldSum",twoTightMuLooseEName+"_YieldSum",RooArgSet(twoTightMuLooseESignalYield,twoTightMuLooseETopWJetsYield));

  /*
 
  RooPoisson oneTightMuConstraint(oneTightMuName+"_Constraint",oneTightMuName+"_Constraint",oneTightMuCount,oneTightMuYieldSum);
  wspace.import( oneTightMuConstraint,RecycleConflictNodes() );

  RooPoisson twoTightMuConstraint(twoTightMuName+"_Constraint",twoTightMuName+"_Constraint",twoTightMuCount,twoTightMuYieldSum);
  wspace.import( twoTightMuConstraint,RecycleConflictNodes() );

  RooPoisson twoTightMuLooseMuConstraint(twoTightMuLooseMuName+"_Constraint",twoTightMuLooseMuName+"_Constraint",twoTightMuLooseMuCount,twoTightMuLooseMuYieldSum);
  wspace.import( twoTightMuLooseMuConstraint,RecycleConflictNodes() );

  RooPoisson twoTightMuLooseEConstraint(twoTightMuLooseEName+"_Constraint",twoTightMuLooseEName+"_Constraint",twoTightMuLooseECount,twoTightMuLooseEYieldSum);
  wspace.import( twoTightMuLooseEConstraint,RecycleConflictNodes() );

  */

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
// MAKE TAU->HAD AND DITAU PREDICTIONS FOR THIS NB/MET/HT BIN 

void makeTauHadBinPrediction( RooWorkspace& wspace, TString binname, vector<TString> allbinnames, 
		    double sfactor1Tau, double sfactor2Tau, double sfactorMuTau, double sfactorETau ){

  TString want = binname(0,1);

  /////////////////////////////////////////////////////////////////////////////////////
  // COMBINE ALL YIELDS IN SAME NB BIN
  
  RooArgSet *Loose_oneTightMu = new RooArgSet("Loose_oneTightMu");
  RooArgSet *Loose_twoTightMu = new RooArgSet("Loose_twoTightMu");
  RooArgSet *Loose_twoTightMuLooseMu = new RooArgSet("Loose_twoTightMuLooseMu");
  RooArgSet *Loose_twoTightMuLooseE = new RooArgSet("Loose_twoTightMuLooseE");
  

  /////////////////////////////////////////////////////////////////////////////////////
  // MAKE TAU->HAD PREDICTION FOR THIS NB/MET/HT BIN BY LOOPING OVER ALL OTHER MET/HT BINS
  // INCLUDES DITAU COMPONENTS

  //   for( int tauhad=0; tauhad<27; tauhad++ ){// 27 NB/MET/HT BINS 
  for( int tauhad=0; tauhad<3; tauhad++ ){// 27 NB/MET/HT BINS 
    
    TString thisbin = allbinnames.at(tauhad);    
    //cout << endl << " LOOPING OVER ALL BIN CONTENT: " << thisbin << endl << endl;

    if( thisbin.BeginsWith(want) ){// LOOP OVER MET/HT BINS WITH CORRECT NB

      //cout << " ADDING BIN CONTENT IN: " << thisbin << endl << endl;
            
      TString oneTightMuName("oneTightMu_");
      oneTightMuName+=thisbin;
      TString twoTightMuName("twoTightMu_");
      twoTightMuName+=thisbin;
      TString twoTightMuLooseMuName("twoTightMuLooseMu_");
      twoTightMuLooseMuName+=thisbin;
      TString twoTightMuLooseEName("twoTightMuLooseE_");
      twoTightMuLooseEName+=thisbin;
      
      ///////////////////////////////////////////////////////////////////////////////////
      // THE 1 TIGHT MU PIECES ARE ALREADY INITIALIZED AND CONSTRAINED IN BINS OF DTHETA

      //cout << " GETTING DTHETA-BINNED TTBAR YIELDS " << endl;

      TString oneTightMuTopWJetsYieldName1(oneTightMuName+"_TopWJetsYield_Theta1");	
      RooRealVar oneTightMuTopWJetsYield1 = *wspace.var(oneTightMuTopWJetsYieldName1.Data());
	
      TString oneTightMuTopWJetsYieldName2(oneTightMuName+"_TopWJetsYield_Theta2");	
      RooRealVar oneTightMuTopWJetsYield2 = *wspace.var(oneTightMuTopWJetsYieldName2.Data());
	
      TString oneTightMuTopWJetsYieldName3(oneTightMuName+"_TopWJetsYield_Theta3");	
      RooRealVar oneTightMuTopWJetsYield3 = *wspace.var(oneTightMuTopWJetsYieldName3.Data());
	
      TString oneTightMuTopWJetsYieldName4(oneTightMuName+"_TopWJetsYield_Theta4");	
      RooRealVar oneTightMuTopWJetsYield4 = *wspace.var(oneTightMuTopWJetsYieldName4.Data());
	
      TString oneTightMuTopWJetsYieldName5(oneTightMuName+"_TopWJetsYield_Theta5");	
      RooRealVar oneTightMuTopWJetsYield5 = *wspace.var(oneTightMuTopWJetsYieldName5.Data());

      /////////////////////////////////////////////////////////////////////////////////////
      // PRODUCE DITAU PREDICTIONS FROM "THIS BIN" TO "BINNAME" (NOT EVER BINNED IN DTHETA) 

      //cout << " GETTING DILEPTON YIELDS " << endl;

      TString twoTightMuTopWJetsYieldName(twoTightMuName+"_TopWJetsYield");	
      RooRealVar twoTightMuTopWJetsYield = *wspace.var(twoTightMuTopWJetsYieldName.Data());

      TString twoTightMuLooseMuTopWJetsYieldName(twoTightMuLooseMuName+"_TopWJetsYield");	
      RooRealVar twoTightMuLooseMuTopWJetsYield = *wspace.var(twoTightMuLooseMuTopWJetsYieldName.Data());

      TString twoTightMuLooseETopWJetsYieldName(twoTightMuLooseEName+"_TopWJetsYield");	
      RooRealVar twoTightMuLooseETopWJetsYield = *wspace.var(twoTightMuLooseETopWJetsYieldName.Data());

      /////////////////////////////////////////////////////////////////////////////////////
      // SUM LOOSE TIGHT MU SAMPLES
  
      //cout << " ADDING HT/MET BINNED TAU PREDICTIONS TO COLLECTION " << endl;

      Loose_oneTightMu->addClone(oneTightMuTopWJetsYield1);
      Loose_oneTightMu->addClone(oneTightMuTopWJetsYield2);
      Loose_oneTightMu->addClone(oneTightMuTopWJetsYield3);
      Loose_oneTightMu->addClone(oneTightMuTopWJetsYield4);
      Loose_oneTightMu->addClone(oneTightMuTopWJetsYield5);
      Loose_twoTightMu->addClone(twoTightMuTopWJetsYield);
      Loose_twoTightMuLooseMu->addClone(twoTightMuLooseMuTopWJetsYield);
      Loose_twoTightMuLooseE->addClone(twoTightMuLooseETopWJetsYield);

      //Loose_oneTightMu->Print();

    }// end of other bins matched to correct NBS
  }// end of tau->had loop over other HT/MET bins
 
  /////////////////////////////////////////////////////////////////////////////////////
  // PRODUCE *OUTSIDE* (BY NB BIN) TAU->HAD AND DITAU CONTRIBUTIONS TO THIS 0L PREDICTION
  

   //cout << " GETTING TAU PREDICTIONS OUTSIDE HT/MET BINS " << endl;

  TString oneTightMuName("oneTightMu_"+want+"Outside");
  TString twoTightMuName("twoTightMu_"+want+"Outside");
  TString twoTightMuLooseMuName("twoTightMuLooseMu_"+want+"Outside");
  TString twoTightMuLooseEName("twoTightMuLooseE_"+want+"Outside");

  TString oneTightMuTopWJetsYieldName(oneTightMuName+"_TopWJetsYield");

  //cout << oneTightMuTopWJetsYieldName << endl;
  
  RooRealVar oneTightMuTopWJetsYield = *wspace.var(oneTightMuTopWJetsYieldName.Data());
  
  TString twoTightMuTopWJetsYieldName(twoTightMuName+"_TopWJetsYield");
  RooRealVar twoTightMuTopWJetsYield = *wspace.var(twoTightMuTopWJetsYieldName.Data());
  
  TString twoTightMuLooseMuTopWJetsYieldName(twoTightMuLooseMuName+"_TopWJetsYield");
  RooRealVar twoTightMuLooseMuTopWJetsYield = *wspace.var(twoTightMuLooseMuTopWJetsYieldName.Data());
  
  TString twoTightMuLooseETopWJetsYieldName(twoTightMuLooseEName+"_TopWJetsYield");
  RooRealVar twoTightMuLooseETopWJetsYield = *wspace.var(twoTightMuLooseETopWJetsYieldName.Data());
  
  ///////////////////////////////////////////////////////////////////////////////////
  // ADD "OUTSIDE" CONTRIBUTIONS TO LOOSE SELECTION

  Loose_oneTightMu->addClone(oneTightMuTopWJetsYield);
  Loose_twoTightMu->addClone(twoTightMuTopWJetsYield);
  Loose_twoTightMuLooseMu->addClone(twoTightMuLooseMuTopWJetsYield);
  Loose_twoTightMuLooseE->addClone(twoTightMuLooseETopWJetsYield);
  
  
  /////////////////////////////////////////////////////////////////////////////////////
  // PRODUCE DITAU PREDICTIONS FROM "LOOSE" TO "BINNAME" 
  
  
  TString oneTightMuName2("oneTightMu_");
  oneTightMuName2+=binname;
  TString twoTightMuName2("twoTightMu_");
  twoTightMuName2+=binname;
  TString twoTightMuLooseMuName2("twoTightMuLooseMu_");
  twoTightMuLooseMuName2+=binname;
  TString twoTightMuLooseEName2("twoTightMuLooseE_");
  twoTightMuLooseEName2+=binname;

  //  cout << " BEFORE SETTING SUMS " << endl;

  //Loose_oneTightMu->Print();

  RooAddition oneTightMuYieldSum(oneTightMuName2+"_TopWJetsYieldSum",oneTightMuName2+"_TopWJetsYieldSum",*Loose_oneTightMu);

  RooAddition twoTightMuYieldSum(twoTightMuName2+"_TopWJetsYieldSum",twoTightMuName2+"_TopWJetsYieldSum",*Loose_twoTightMu);

  RooAddition twoTightMuLooseMuYieldSum(twoTightMuLooseMuName2+"_TopWJetsYieldSum",twoTightMuLooseMuName2+"_TopWJetsYieldSum",*Loose_twoTightMuLooseMu);

  RooAddition twoTightMuLooseEYieldSum(twoTightMuLooseEName2+"_TopWJetsYieldSum",twoTightMuLooseEName2+"_TopWJetsYieldSum",*Loose_twoTightMuLooseE);

  RooRealVar ScaleFactor1Tau("oneTightMu_"+binname+"_ScaleFactorTau","oneTightMu_"+binname+"_ScaleFactorTau",sfactor1Tau,0.000001,100000.);
  ScaleFactor1Tau.setConstant();
  RooRealVar ScaleFactor2Tau("twoTightMu_"+binname+"_ScaleFactorTau","twoTightMu_"+binname+"_ScaleFactorTau",sfactor2Tau,0.000001,100000.);
  ScaleFactor2Tau.setConstant();
  RooRealVar ScaleFactorMuTau("twoTightMuLooseMu_"+binname+"_ScaleFactorTau","twoTightMuLooseMu_"+binname+"_ScaleFactorTau",sfactorMuTau,0.000001,100000.);
  ScaleFactorMuTau.setConstant();
  RooRealVar ScaleFactorETau("twoTightMuLooseE_"+binname+"_ScaleFactorTau","twoTightMuLooseE_"+binname+"_ScaleFactorTau",sfactorETau,0.000001,100000.);
  ScaleFactorETau.setConstant();

  //  cout << " APPLYING SCALE FACTORS TO SUMS " << endl;

  RooProduct  zeroLeptonTopWJetsYield1Tau("zeroLepton_"+binname+"_TopWJetsYield_1Tau","zeroLepton_"+binname+"_TopWJetsYield_1Tau",
					  RooArgSet(ScaleFactor1Tau,oneTightMuYieldSum));
  
  RooProduct  zeroLeptonTopWJetsYield2Tau("zeroLepton_"+binname+"_TopWJetsYield_2Tau","zeroLepton_"+binname+"_TopWJetsYield_2Tau",
					  RooArgSet(ScaleFactor2Tau,twoTightMuYieldSum));
  
  RooProduct  zeroLeptonTopWJetsYieldMuTau("zeroLepton_"+binname+"_TopWJetsYield_MuTau","zeroLepton_"+binname+"_TopWJetsYield_MuTau",
					   RooArgSet(ScaleFactorMuTau,twoTightMuLooseMuYieldSum));
  
  RooProduct  zeroLeptonTopWJetsYieldETau("zeroLepton_"+binname+"_TopWJetsYield_ETau","zeroLepton_"+binname+"_TopWJetsYield_ETau",
					  RooArgSet(ScaleFactorETau,twoTightMuLooseEYieldSum));
  
  
  /////////////////////////////////////////////////////////////////////////////////////
  // ADD TAU->HAD AND DITAU CONTRIBUTIONS TO 0L PREDICTION TO WORKSPACE
  
  //  cout << " IMPORTING PREDICTIONS " << endl;

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


void buildMRLikelihood( TString outputFile, TString bincontentFileName, TString outsidebincontentFileName, 
		      TString scalefactorsFileName, TString tauhadscalefactorsFileName,
		      TString signalfractionsFileName, TString outsidesignalfractionsFileName ) 
{

  ///////////////////////////////////////////////////

  int error = 0;
  TFile *test = new TFile(outputFile.Data(),"RECREATE");

  ///////////////////////////////////////////////////

  RooWorkspace* wspace = new RooWorkspace("wspace");

  RooRealVar signalCrossSection("signalCrossSection","signalCrossSection",0.,0.,10.);
  wspace->import(signalCrossSection);

  wspace->defineSet("namesfordata","");
  wspace->defineSet("nuisances","");
  vector<TString> allbinnames;
  allbinnames.clear();

  //  cout << " ACCESS INPUT FILES " << endl;

  ////////////////////////////////////
  ////// Get content from file

  ifstream contentFile;  
  contentFile.open(bincontentFileName.Data(),fstream::in);
 
  ////////////////////////////////////
  ////// Get "outside" content from file

  ifstream outsidecontentFile;  
  outsidecontentFile.open(outsidebincontentFileName.Data(),fstream::in);
  
  ////////////////////////////////////
  ////// Get SM dtheta and dilep scale factors from file

  ifstream scalefactorsFile;  
  scalefactorsFile.open(scalefactorsFileName.Data(),fstream::in);
 
  ////////////////////////////////////
  ////// Get SM tau->had and ditau scale factors from file

  ifstream tauhadscalefactorsFile;  
  tauhadscalefactorsFile.open(tauhadscalefactorsFileName.Data(),fstream::in);
  
  ////////////////////////////////////
  ////// Get SUSY model signal fractions from file

  ifstream sigfracFile;  
  sigfracFile.open(signalfractionsFileName.Data(),fstream::in);
    
  ////////////////////////////////////
  ////// Get "outside" SUSY model signal fractions from file

  ifstream sigfracoutsideFile;  
  sigfracoutsideFile.open(outsidesignalfractionsFileName.Data(),fstream::in);


  //  cout << " GET COUNTS AND SIG FRAC OUTSIDE POL METHOD " << endl;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////// 
  // GET AND SET oneTightMu, twoTightMu, twoTightMuLooseMu, and twoTightMuLooseE counts that are OUTSIDE STANDARD BINNING!!!!!!

  for( int i=1; i<4; i++ ){

    TString binname1, binname2;    
    double one_mutight, two_mutight, two_mutightmuloose, two_mutighteloose, sigone_mutight, sigtwo_mutight, sigtwo_mutightmuloose, sigtwo_mutighteloose;
    
    outsidecontentFile>>binname1;
    outsidecontentFile>>one_mutight;
    outsidecontentFile>>two_mutight;
    outsidecontentFile>>two_mutightmuloose;
    outsidecontentFile>>two_mutighteloose;
    
    sigfracoutsideFile>>binname2;
    sigfracoutsideFile>>sigone_mutight;
    sigfracoutsideFile>>sigtwo_mutight;
    sigfracoutsideFile>>sigtwo_mutightmuloose;
    sigfracoutsideFile>>sigtwo_mutighteloose;
    
    TString oneTightMuName("oneTightMu_");
    oneTightMuName+=i;
    oneTightMuName.Append("Outside_Count");
    RooRealVar oneTightMuCount(oneTightMuName,oneTightMuName,one_mutight);
    oneTightMuCount.setConstant();
    wspace->import(oneTightMuCount);
    wspace->extendSet("namesfordata",oneTightMuCount.GetName());
    
    TString twoTightMuName("twoTightMu_");
    twoTightMuName+=i;
    twoTightMuName.Append("Outside_Count");
    RooRealVar twoTightMuCount(twoTightMuName,twoTightMuName,two_mutight);
    twoTightMuCount.setConstant();
    wspace->import(twoTightMuCount);
    wspace->extendSet("namesfordata",twoTightMuCount.GetName());
    
    TString twoTightMuLooseMuName("twoTightMuLooseMu_");
    twoTightMuLooseMuName+=i;
    twoTightMuLooseMuName.Append("Outside_Count");
    RooRealVar twoTightMuLooseMuCount(twoTightMuLooseMuName,twoTightMuLooseMuName,two_mutightmuloose);
    twoTightMuLooseMuCount.setConstant();
    wspace->import(twoTightMuLooseMuCount);
    wspace->extendSet("namesfordata",twoTightMuLooseMuCount.GetName());
    
    TString twoTightMuLooseEName("twoTightMuLooseE_");
    twoTightMuLooseEName+=i;
    twoTightMuLooseEName.Append("Outside_Count");
    RooRealVar twoTightMuLooseECount(twoTightMuLooseEName,twoTightMuLooseEName,two_mutighteloose);
    twoTightMuLooseECount.setConstant();
    wspace->import(twoTightMuLooseECount);
    wspace->extendSet("namesfordata",twoTightMuLooseECount.GetName());
    
    ///////////////////////////////////////////////////////////////////////////////////////
    // CONSTRAIN "OUTSIDE" BINS BY NB ONLY
    
    makeOutsideUnbinnedConstraints( *wspace, i, sigone_mutight, sigtwo_mutight, sigtwo_mutightmuloose, sigtwo_mutighteloose );
    
  }

  //  cout << " GET COUNTS AND SIG FRAC FOR POL METHOD " << endl;

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  // GET AND SET ALL COUNTS DEPENDENT ON HT/MET BINS 


  // bins in NBs, HT, MET - these counts are found via the polarization code
  // 3 x 3 x 3
  //  for( int i=0; i<27; i++ ){
  for( int i=0; i<3; i++ ){
 
    TString binname, binname2;  
    double zerolepton, signal0, one_mutight_theta, one_muloose_theta, one_eloose_theta, two_mutight, two_mutightmuloose, two_muloose, two_mutighteloose, two_mulooseeloose, two_eloose;

    contentFile>>binname;
    contentFile>>zerolepton;

    sigfracFile>>binname2;
    sigfracFile>>signal0;
 
    allbinnames.push_back( binname );

    //    cout << " BEGINNING HT/MET LOOP " << binname << "  " << binname2 << endl;

    TString zeroLeptonName("zeroLepton_");
    zeroLeptonName+=binname;
    //////
    TString oneTightMuName("oneTightMu_");
    oneTightMuName+=binname;
    
    TString oneLooseMuName("oneLooseMu_");
    oneLooseMuName+=binname;

    TString oneLooseEName("oneLooseE_");
    oneLooseEName+=binname;
    //////
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
    ////// 
    RooRealVar zeroLeptonCount(zeroLeptonName+"_Count",zeroLeptonName+"_Count",zerolepton);
    zeroLeptonCount.setConstant();
    wspace->import(zeroLeptonCount);
    wspace->extendSet("namesfordata",zeroLeptonCount.GetName());
 
    RooRealVar zeroLeptonSignalFrac(zeroLeptonName+"_SignalFrac",zeroLeptonName+"_SignalFrac",signal0);
    zeroLeptonSignalFrac.setConstant();
    wspace->import(zeroLeptonSignalFrac);
 

    // input files must go: dtheta1(tmu,mu,e), dtheta2(tmu,mu,e), etc!, dileps

    for( int j=1; j<6; j++ ){
      
      contentFile>>one_mutight_theta;
      contentFile>>one_muloose_theta;
      contentFile>>one_eloose_theta;
 
      TString oneTightMuThetaName(oneTightMuName+"_Theta");
      oneTightMuThetaName+=j;
      TString oneLooseMuThetaName(oneLooseMuName+"_Theta");
      oneLooseMuThetaName+=j;
      TString oneLooseEThetaName(oneLooseEName+"_Theta");
      oneLooseEThetaName+=j;

      RooRealVar oneTightMuThetaCount(oneTightMuThetaName+"_Count",oneTightMuThetaName+"_Count",one_mutight_theta);
      oneTightMuThetaCount.setConstant();
      wspace->import(oneTightMuThetaCount);
      wspace->extendSet("namesfordata",oneTightMuThetaCount.GetName());
  
      RooRealVar oneLooseMuThetaCount(oneLooseMuThetaName+"_Count",oneLooseMuThetaName+"_Count",one_muloose_theta);
      oneLooseMuThetaCount.setConstant();
      wspace->import(oneLooseMuThetaCount);
      wspace->extendSet("namesfordata",oneLooseMuThetaCount.GetName());
 
      RooRealVar oneLooseEThetaCount(oneLooseEThetaName+"_Count",oneLooseEThetaName+"_Count",one_eloose_theta);
      oneLooseEThetaCount.setConstant();
      wspace->import(oneLooseEThetaCount);
      wspace->extendSet("namesfordata",oneLooseEThetaCount.GetName());
  
    }
    ////// end of loop over dtheta-binned SL content, now getting dilepton content
 
    contentFile>>two_mutight;
    contentFile>>two_mutightmuloose;
    contentFile>>two_muloose;
    contentFile>>two_mutighteloose;
    contentFile>>two_mulooseeloose;
    contentFile>>two_eloose;
    
    RooRealVar twoTightMuCount(twoTightMuName+"_Count",twoTightMuName+"_Count",two_mutight);
    twoTightMuCount.setConstant();
    wspace->import(twoTightMuCount);
    wspace->extendSet("namesfordata",twoTightMuCount.GetName());
    
    RooRealVar twoTightMuLooseMuCount(twoTightMuLooseMuName+"_Count",twoTightMuLooseMuName+"_Count",two_mutightmuloose);
    twoTightMuLooseMuCount.setConstant();
    wspace->import(twoTightMuLooseMuCount);
    wspace->extendSet("namesfordata",twoTightMuLooseMuCount.GetName());
    
    RooRealVar twoLooseMuCount(twoLooseMuName+"_Count",twoLooseMuName+"_Count",two_muloose);
    twoLooseMuCount.setConstant();
    wspace->import(twoLooseMuCount);
    wspace->extendSet("namesfordata",twoLooseMuCount.GetName());
    
    RooRealVar twoTightMuLooseECount(twoTightMuLooseEName+"_Count",twoTightMuLooseEName+"_Count",two_mutighteloose);
    twoTightMuLooseECount.setConstant();
    wspace->import(twoTightMuLooseECount);
    wspace->extendSet("namesfordata",twoTightMuLooseECount.GetName());
    
    RooRealVar twoLooseMuLooseECount(twoLooseMuLooseEName+"_Count",twoLooseMuLooseEName+"_Count",two_mulooseeloose);
    twoLooseMuLooseECount.setConstant();
    wspace->import(twoLooseMuLooseECount);
    wspace->extendSet("namesfordata",twoLooseMuLooseECount.GetName());
    
    RooRealVar twoLooseECount(twoLooseEName+"_Count",twoLooseEName+"_Count",two_eloose);
    twoLooseECount.setConstant();
    wspace->import(twoLooseECount);
    wspace->extendSet("namesfordata",twoLooseECount.GetName());

    //    cout << " last counts loaded:  " << twoTightMuLooseECount.getVal() << "  "
    // << twoLooseMuLooseECount.getVal() << "  " << twoLooseECount.getVal() << endl << endl;



    //    cout << " BEFORE MAIN CONSTRAINTS " << endl;

    ///////////////////////////////////////////////////////////////////////////////////////
    // CONSTRAIN 1L "INSIDE" BINS
    // MAKE POLARIZATION METHOD CONTRIBUTIONS TO 0L PREDICTION

    makeBigBin( *wspace, binname, scalefactorsFile, sigfracFile );

    //cout << endl << endl << " WORKSPACE AFTER MAIN CONSTRAINTS SET: " << endl; 
    //wspace->Print();
    //cout << endl << endl << endl;

    //    cout << " BEFORE DIPLEPTON CONSTRAINTS " << endl;

    ///////////////////////////////////////////////////////////////////////////////////////
    // CONSTRAIN DILEPTON "INSIDE" BINS

    makeDileptonConstraintsPrediction( *wspace, binname, scalefactorsFile, sigfracFile );

    //    cout << " AT END OF MAIN HT/MET LOOP " << endl;

  }// end of this loop over HT, MET bins
 

  //  cout << " GET COUNTS AND SIG FRAC FOR TAU METHOD AND FINAL PREDICTION" << endl;

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////


  // loop over other HT, MET regions here to get outside contributions from tau->had in this bin
  // MAKE BIN OF HT AND MET
  //  for( int i=0; i<27; i++ ){
  for( int i=0; i<3; i++ ){

    TString binname = allbinnames.at(i);

    //oneTightMu_1B0HT0MET_TopWJetsYieldYieldSum

    ///////////////////////////////////////////////////////////////////////////////////////
    // MAKE TAU->HAD AND DITAU CONTRIBUTIONS TO 0L PREDICTION (EVERYTHING ALREADY CONSTRAINED AT THIS POINT)

    TString binname1;
    double SFone_mutight, SFtwo_mutight, SFtwo_mutightmuloose, SFtwo_mutighteloose; // <------------------ for now these are just numbers, rather than beta prime functions with errors
    tauhadscalefactorsFile>>binname1;
    tauhadscalefactorsFile>>SFone_mutight;
    tauhadscalefactorsFile>>SFtwo_mutight;
    tauhadscalefactorsFile>>SFtwo_mutightmuloose;
    tauhadscalefactorsFile>>SFtwo_mutighteloose;

    if( binname!=binname1 ){
      cout << endl << endl << " ERROR IN BIN NAMES!!!! " << endl;
      //error = 1;
    }

    makeTauHadBinPrediction( *wspace, binname, allbinnames, SFone_mutight, SFtwo_mutight, SFtwo_mutightmuloose, SFtwo_mutighteloose );
    
    //    cout << " DONE WITH TAUHAD COMPONENT " << endl;

    /////////////////////////////////////////////////////////////////////////////////////
    // NOW COMBINE COMPONENTS OF 0L BACKGROUND PREDICTION IN THIS BIN
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    // GET PIECES OF TTBAR PREDICTION TO ADD TO 0L PREDICTION IN THIS HT/MET BIN
    
    TString zeroLeptonName("zeroLepton_");
    zeroLeptonName+=binname;
    //////
    // cout << " BEFORE GETTING PREDICTIONS " << endl;

    TString PolarizationName1(zeroLeptonName);
    PolarizationName1.Append("_TopWJetsYield_Theta1");
    //RooAbsReal Polarization1 = wspace->function(PolarizationName1.Data());

    TString PolarizationName2(zeroLeptonName);
    PolarizationName2.Append("_TopWJetsYield_Theta2");
    //RooAbsReal Polarization2 = wspace->function(PolarizationName2.Data());
    
    TString PolarizationName3(zeroLeptonName);
    PolarizationName3.Append("_TopWJetsYield_Theta3");
    //RooAbsReal Polarization3 = wspace->function(PolarizationName3.Data());
    
    TString PolarizationName4(zeroLeptonName);
    PolarizationName4.Append("_TopWJetsYield_Theta4");
    //RooAbsReal Polarization4 = wspace->function(PolarizationName4.Data());
    
    TString PolarizationName5(zeroLeptonName);
    PolarizationName5.Append("_TopWJetsYield_Theta5");
    //RooAbsReal Polarization5 = wspace->function(PolarizationName5.Data());

    //    cout << " GOT POLARIZATION PREDICTIONS" << endl;
    //////
    TString TauHadName1(zeroLeptonName);
    TauHadName1.Append("_TopWJetsYield_1Tau");
    //RooAbsReal TauHad = wspace->function(TauHadName1.Data());
    
    TString TauHadName2(zeroLeptonName);
    TauHadName2.Append("_TopWJetsYield_2Tau");
    //RooAbsReal TauTauHad = wspace->function(TauHadName2.Data());
    
    TString TauHadName3(zeroLeptonName);
    TauHadName3.Append("_TopWJetsYield_MuTau");
    //RooAbsReal MuTauHad = wspace->function(TauHadName3.Data());
    
    TString TauHadName4(zeroLeptonName);
    TauHadName4.Append("_TopWJetsYield_ETau");
    //RooAbsReal ETauHad = wspace->function(TauHadName4.Data());

    //    cout << " GOT TAU->HAD PREDICTIONS" << endl;
    //////
    TString DilepName1(zeroLeptonName);
    DilepName1.Append("_TopWJetsYield_MuMu");
    //RooAbsReal Dilep1 = wspace->function(DilepName1.Data());

    TString DilepName2(zeroLeptonName);
    DilepName2.Append("_TopWJetsYield_EMu");
    //RooAbsReal Dilep2 = wspace->function(DilepName2.Data());

    TString DilepName3(zeroLeptonName);
    DilepName3.Append("_TopWJetsYield_EE");
    //RooAbsReal Dilep3 = wspace->function(DilepName3.Data());

    //    cout << " GOT DILEPTON PREDICTIONS" << endl;
    //////
    
    //    cout << " GOT ALL COMPONENTS OF 0L PREDICTION, ADDING THEM... " << endl;

    RooAddition zeroLeptonTopWJetsPolarizationYield(zeroLeptonName+"_TopWJetsPolarizationYield",zeroLeptonName+"_TopWJetsPolarizationYield",
						    RooArgSet( *wspace->function(PolarizationName1.Data()),
							       *wspace->function(PolarizationName2.Data()),
							       *wspace->function(PolarizationName3.Data()),
							       *wspace->function(PolarizationName4.Data()),
							       *wspace->function(PolarizationName5.Data()) ));
    
    RooAddition zeroLeptonTopWJetsTauHadYield(zeroLeptonName+"_TopWJetsTauHadYield",zeroLeptonName+"_TopWJetsTauHadYield",
					      RooArgSet( *wspace->function(TauHadName1.Data()),
							 *wspace->function(TauHadName2.Data()),
							 *wspace->function(TauHadName3.Data()),
							 *wspace->function(TauHadName4.Data()) ));
    
    RooAddition zeroLeptonTopWJetsDileptonYield(zeroLeptonName+"_TopWJetsDileptonYield",zeroLeptonName+"_TopWJetsDileptonYield",
						RooArgSet( *wspace->function(DilepName1.Data()),
							   *wspace->function(DilepName2.Data()),
							   *wspace->function(DilepName3.Data()) )); 
    
       
    /////////////////////////////////////////////////////////////////////////////////////
    // GET AND APPLY SIGNAL FRACTION FOR THIS 0L BIN, ADD TO ALL PREDICTIONS
 
    //    cout << " GET 0L COUNTS AND SIG FRAC " << endl;

    TString zeroLeptonCountName(zeroLeptonName);
    zeroLeptonCountName.Append("_Count");
    //RooRealVar zeroLeptonCount = wspace->var(zeroLeptonCountName.Data());

    TString zeroLeptonSignalFracName(zeroLeptonName);
    zeroLeptonSignalFracName.Append("_SignalFrac");
    //RooRealVar zeroLeptonSignalFrac = wspace->var(zeroLeptonSignalFracName.Data());

 
    RooProduct zeroLeptonSignalYield(zeroLeptonName+"_SignalYield",zeroLeptonName+"_SignalYield",
				     RooArgSet(signalCrossSection,*wspace->var(zeroLeptonSignalFracName.Data())));

    RooAddition zeroLeptonYieldSum(zeroLeptonName+"_YieldSum",zeroLeptonName+"_YieldSum",
				   RooArgSet( zeroLeptonSignalYield, zeroLeptonTopWJetsPolarizationYield, zeroLeptonTopWJetsTauHadYield, zeroLeptonTopWJetsDileptonYield ));
		
	
    //    cout << " SET 0L CONSTRAINTS " << endl;

    /////////////////////////////////////////////////////////////////////////////////////
    // SET CONSTRAINTS FOR THIS 0L BIN
    
    RooPoisson zeroLeptonConstraint(zeroLeptonName+"_Constraint",zeroLeptonName+"_Constraint",
				    *wspace->var(zeroLeptonCountName.Data()),zeroLeptonYieldSum);
				    
    wspace->import( zeroLeptonConstraint,RecycleConflictNodes() );

    
  }// end of FINAL loop over HT, MET bins
  

  //  cout << endl << endl << " OUTSIDE HT MET BINS, SETTING ALL CONSTRAINTS " << endl;
 
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////
  contentFile.close();
  outsidecontentFile.close();
  scalefactorsFile.close();
  tauhadscalefactorsFile.close();
  sigfracFile.close();
  ////////////////////////////////////


  
  RooArgSet Constraints(wspace->allPdfs());
  
  RooProdPdf model("model","model", Constraints );
  wspace->import(model);

  //cout << "data, size: " << (*wspace->set("namesfordata")).getSize() << endl;
  //(*wspace->set("namesfordata")).Print("v");

  RooDataSet dataset("dataset","dataset",*wspace->set("namesfordata"));

  dataset.add(*wspace->set("namesfordata"));

  wspace->import(dataset);

 
  //  cout << endl << endl << " FINAL SET OF DATA " << endl;
  //dataset.Print("v");
  //  cout << endl << endl << " FINAL SET OF CONSTRAINTS " << endl;
  //Constraints.Print("v");
  
  wspace->defineSet("poi","signalCrossSection");
  
  ////////////////////////////////////

  //Hesse(false),Minos(false),
  //Save(),Constrain(Constraints)
  //	     );
    
 

  ModelConfig * modelConfig = new ModelConfig("modelConfig");
  modelConfig->SetWorkspace(*wspace);
  modelConfig->SetPdf("model");
  modelConfig->SetParametersOfInterest(*wspace->set("poi"));
  modelConfig->SetNuisanceParameters(*wspace->set("nuisances"));

  //cout << " before model print " << endl;
  //modelConfig->Print("v");
  //cout << " after model print " << endl;



  // Declare parameter of interest
  RooArgSet parofinterest(*wspace->set("poi")) ;
  //ProfileLikelihoodCalculator plc(dataset, *wspace->pdf("model"), parofinterest);
  //plc.SetTestSize(0.05);
  ProfileLikelihoodCalculator plc(dataset, *modelConfig);
  plc.SetConfidenceLevel(0.95);

  //  cout << " before pl " << endl;

  LikelihoodInterval* interval = plc.GetInterval() ;
  LikelihoodIntervalPlot*  lplot = new LikelihoodIntervalPlot(interval);

  //  interval->Print("v");

  RooRealVar* firstPOI = (RooRealVar*) parofinterest.first();

  double profileLikelihoodUpperLimit = interval->UpperLimit(*firstPOI);
  double profileLikelihoodLowerLimit = interval->LowerLimit(*firstPOI);
  
  firstPOI->setRange(profileLikelihoodLowerLimit,profileLikelihoodUpperLimit);
  
  TCanvas* canvas = new TCanvas("Result", "Result", 10, 10,500,500);
  lplot->Draw();
  //lplot->Save();
  canvas->Update();
 
  
  ////////////////////////////////////
  model.fitTo(dataset);
  ////////////////////////////////////


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
  
  test->Write();


}// end of likelihood builder





