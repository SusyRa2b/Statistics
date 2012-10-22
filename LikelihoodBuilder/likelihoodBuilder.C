#include <iostream>
#include <string.h>
#include <complex>
#include <map>
#include <cassert>

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

#include "RooCorrelatedBetaGeneratorHelper.h"
#include "RooCorrelatedBetaPrimeGeneratorHelper.h"
#include "betaHelperFunctions.h"
#include "RooBetaInverseCDF.h"
#include "RooBetaPrimeInverseCDF.h"
#include "rooFitBetaHelperFunctions.h"
#include "RooNormalFromFlatPdf.h"

//For truncated gaussians, uncomment this line and two lines in setup.C
//You must also check out RA2b/Statistics/3Dcode to get a couple of files
//#include "rooFitGaussianHelperFunctions.h"

#include "RooProdPdfLogSum.h"
#include "RooPoissonLogEval.h"

#include "metReweightingBuilder.h"

#include "RooStats/ModelConfig.h"

#include "TMath.h"

using namespace RooFit ;
using namespace RooStats ;

struct  likelihoodOptions
{
  bool skipTriggerEfficiency;
  TString qcdMethod;
  TString TopWJetsMethod;
};

struct channels
{
  double zeroLepton;
  double zeroLeptonLowDeltaPhiN;
  double oneLepton;
  double oneElectron;
  double oneMuon;
  TString diMuonName;
  double diMuon;
  TString diElectronName;
  double diElectron;
} ;

struct yields
{
  double zeroLepton;
  double zeroLeptonLowDeltaPhiN;
  double oneLepton;
} ;

struct abcdBinParameters
{
  double zeroLeptonTriggerEfficiency;
  double zeroLeptonTriggerEfficiencyError;
  double zeroLeptonLowDeltaPhiNTriggerEfficiency;
  double zeroLeptonLowDeltaPhiNTriggerEfficiencyError;
  double oneLeptonTriggerEfficiency;
  double oneLeptonTriggerEfficiencyError;
  double oneElectronTriggerEfficiency;
  double oneElectronTriggerEfficiencyError;
  double oneMuonTriggerEfficiency;
  double oneMuonTriggerEfficiencyError;
  double zeroLeptonLowDeltaPhiNMC;
  double topWJetsLowDeltaPhiNOverZeroLeptonRatioMC;
  double ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC;
  TString ZtoNuNubTagScalingName;
  double ZtoNuNubTagScaling;
  double ZtoNuNubTagScalingError;
  TString ZtomumuAcceptanceName;
  double ZtomumuAcceptance;
  double ZtomumuAcceptanceError;
  TString ZtoeeAcceptanceName;
  double ZtoeeAcceptance;
  double ZtoeeAcceptanceError;
  TString deltaPhiNRatioName;
  double deltaPhiNRatio;
  double deltaPhiNRatioError;
  TString lowDeltaPhiNScalingName;
  double qcdClosure;
  double qcdClosureError;
  double topWJetsClosure;
  double topWJetsClosureError;
  TString ZtoeeSystematicName;
  double ZtoeeSystematic;
  double ZtoeeSystematicError;
  TString ZtomumuSystematicName;
  double ZtomumuSystematic;
  double ZtomumuSystematicError;
} ;

struct allBinNames
{
  TString ZtoeeInvPurity;
  TString ZtomumuInvPurity;
  TString ZtoeeEfficiency;
  TString ZtomumuEfficiency;
  TString ZtollOverZtoNuNuRatio;
  TString singleLeptonScaling;
  TString MCUncertainty;
  TString signalCrossSection;
  TString observables;
  TString nuisances;
  TString signalUncertainty;
} ;

struct allBins
{
  double ZtollOverZtoNuNuRatio;
  double ZtoeePurity;
  double ZtoeePurityError;
  double ZtomumuPurity;
  double ZtomumuPurityError;
  double ZtoeeEfficiency;
  double ZtoeeEfficiencyError;
  double ZtomumuEfficiency;
  double ZtomumuEfficiencyError;
  double qcdClosure;
  double qcdClosureError;
  double topWJetsClosure;
  double topWJetsClosureError;
  double MCUncertainty;
} ;


struct sOutsideSignalMR
{
  TString oneTightMu_Outside_SignalRawCountName;
  int oneTightMu_Outside_SignalRawCount;
  TString twoTightMu_Outside_SignalRawCountName;
  int twoTightMu_Outside_SignalRawCount;
} ;

bool makeOneBin(const likelihoodOptions options, RooWorkspace& ws , TString& binName , allBinNames& names , const allBins& numbers, channels& observed , abcdBinParameters& abcd ,  double& oneLeptonTotal, double& zeroLeptonTopWJetsGuess)
{

  //////////////////////////
  // QCD
  //////////////////////////
  
  //Define unique names
  TString zeroLeptonName("zeroLepton_");
  zeroLeptonName+=binName;
  TString zeroLeptonLowDeltaPhiNName("zeroLeptonLowDeltaPhiN_");
  zeroLeptonLowDeltaPhiNName+=binName;
  
  //Define unique counts and add them to the workspace
  RooRealVar zeroLeptonCount(zeroLeptonName+"_Count",zeroLeptonName+"_Count",observed.zeroLepton);
  RooRealVar zeroLeptonLowDeltaPhiNCount(zeroLeptonLowDeltaPhiNName+"_Count",zeroLeptonLowDeltaPhiNName+"_Count",observed.zeroLeptonLowDeltaPhiN);
  zeroLeptonCount.setConstant();
  zeroLeptonLowDeltaPhiNCount.setConstant();

  ws.import(zeroLeptonCount);
  ws.import(zeroLeptonLowDeltaPhiNCount);
  ws.extendSet(names.observables,zeroLeptonCount.GetName());
  ws.extendSet(names.observables,zeroLeptonLowDeltaPhiNCount.GetName());

  TString scaling = "lowDeltaPhiNScaling_";   //for now, make this in option here. consider moving option to input file "scalingName"
  if(options.qcdMethod == "singleScaleWithCorrections") scaling +=  "all";
  else if (options.qcdMethod == "htDependent") scaling +=  abcd.lowDeltaPhiNScalingName;
  else assert(0);
  RooRealVar* lowDeltaPhiNScaling = ws.var(scaling);
  if(lowDeltaPhiNScaling == NULL) {
    RooRealVar lowDeltaPhiNScaling_temp(scaling, scaling, 0.0, 1e2);
    lowDeltaPhiNScaling_temp.setVal(0.2);
    ws.import(lowDeltaPhiNScaling_temp);
    ws.extendSet(names.nuisances,lowDeltaPhiNScaling_temp.GetName());
    lowDeltaPhiNScaling = ws.var(scaling);
  }
  double qcdGuess = observed.zeroLeptonLowDeltaPhiN - abcd.zeroLeptonLowDeltaPhiNMC;
  if(qcdGuess < 1e-5) qcdGuess = 1e-5;
  RooRealVar zeroLeptonLowDeltaPhiNQCDYield(zeroLeptonLowDeltaPhiNName+"_QCDYield",zeroLeptonLowDeltaPhiNName+"_QCDYield",qcdGuess,1e-5,1e6);
  ws.import(zeroLeptonLowDeltaPhiNQCDYield);
  ws.extendSet(names.nuisances,zeroLeptonLowDeltaPhiNQCDYield.GetName());
  RooRealVar* zeroLeptonQCDClosure = (RooRealVar*)
    getBetaPrimeConstraint(ws,"zeroLeptonQCDClosure_", binName,
			   abcd.qcdClosure,abcd.qcdClosureError,
			   names.observables,names.nuisances);
  RooProduct zeroLeptonQCDYield(zeroLeptonName+"_QCDYield",zeroLeptonName+"_QCDYield",RooArgSet(*lowDeltaPhiNScaling,*zeroLeptonQCDClosure,zeroLeptonLowDeltaPhiNQCDYield));
  
  
  //////////////////////////////
  // Top and W+jets 
  //////////////////////////////

  TString oneLeptonName("oneLepton_");
  oneLeptonName+=binName;

  RooRealVar oneLeptonCount(oneLeptonName+"_Count",oneLeptonName+"_Count",observed.oneLepton);

  oneLeptonCount.setConstant();
  ws.import(oneLeptonCount);
  ws.extendSet(names.observables,oneLeptonCount.GetName());
  
  RooRealVar* singleLeptonScaling = ws.var(names.singleLeptonScaling);
  RooRealVar  oneLeptonTopWJetsYield(oneLeptonName+"_TopWJetsYield",oneLeptonName+"_TopWJetsYield",observed.oneLepton,1e-5,1e5);
  ws.import(oneLeptonTopWJetsYield);
  ws.extendSet(names.nuisances,oneLeptonTopWJetsYield.GetName());
  oneLeptonTotal += oneLeptonTopWJetsYield.getVal();
  RooRealVar*  zeroLeptonTopWJetsClosure = (RooRealVar*)
    getBetaPrimeConstraint(ws,"zeroLeptonTopWJetsClosure_", binName,
			   abcd.topWJetsClosure,abcd.topWJetsClosureError,
			   names.observables,names.nuisances);
  RooProduct  zeroLeptonTopWJetsYield(zeroLeptonName+"_TopWJetsYield",zeroLeptonName+"_TopWJetsYield",RooArgSet(*singleLeptonScaling,*zeroLeptonTopWJetsClosure,oneLeptonTopWJetsYield));
  

  /////////////////////////////////
  // Z to invisible
  /////////////////////////////////
  
  //-----Define Z->ll observables
  
  RooRealVar* diElectronCount = ws.var("diElectron_"+observed.diElectronName+"_Count");
  if(diElectronCount == NULL) 
    {
      RooRealVar diElectronCount_temp("diElectron_"+observed.diElectronName+"_Count","diElectron_"+observed.diElectronName+"_Count",observed.diElectron); 
      diElectronCount_temp.setConstant();
      ws.import(diElectronCount_temp);
      diElectronCount = ws.var("diElectron_"+observed.diElectronName+"_Count");
      ws.extendSet(names.observables,diElectronCount->GetName());
    }
  RooRealVar* diMuonCount = ws.var("diMuon_"+observed.diMuonName+"_Count");
  if(diMuonCount == NULL) 
    {
      RooRealVar diMuonCount_temp("diMuon_"+observed.diMuonName+"_Count","diMuon_"+observed.diMuonName+"_Count",observed.diMuon);
      diMuonCount_temp.setConstant();
      ws.import(diMuonCount_temp);
      diMuonCount = ws.var("diMuon_"+observed.diMuonName+"_Count");
      ws.extendSet(names.observables,diMuonCount->GetName());
    }
  
  //-----Define Z->ll yields
  
  RooRealVar* ZtollOverZtoNuNuRatio = ws.var(names.ZtollOverZtoNuNuRatio);
  
  RooRealVar* ZtoNuNu = ws.var("ZtoNuNu_"+observed.diMuonName);//use MuonName, which should be the same as ElectronName
  if(ZtoNuNu == NULL) 
    {
      RooRealVar ZtoNuNu_temp("ZtoNuNu_"+observed.diMuonName,"ZtoNuNu_"+observed.diMuonName,(1./numbers.ZtollOverZtoNuNuRatio)*0.5*(observed.diMuon*numbers.ZtomumuPurity/(numbers.ZtomumuEfficiency*abcd.ZtomumuAcceptance) + observed.diElectron*numbers.ZtoeePurity/(numbers.ZtoeeEfficiency*abcd.ZtoeeAcceptance)),0.0,1e5);
      ws.import(ZtoNuNu_temp);
      ZtoNuNu = ws.var("ZtoNuNu_"+observed.diMuonName);
      ws.extendSet(names.nuisances,ZtoNuNu->GetName());
    }
  
  RooProduct* Ztoll = (RooProduct*)ws.arg("Ztoll_"+observed.diMuonName);//use MuonName, which should be the same as ElectronName
  if(Ztoll == NULL) 
    {
      RooProduct Ztoll_temp("Ztoll_"+observed.diMuonName,"Ztoll_"+observed.diMuonName,RooArgSet(*ZtoNuNu,*ZtollOverZtoNuNuRatio));
      ws.import(Ztoll_temp, RecycleConflictNodes());
      Ztoll = (RooProduct*)ws.arg("Ztoll_"+observed.diMuonName);//use MuonName, which should be the same as ElectronName
    }
  
  RooAbsArg* ZtomumuEfficiency = ws.arg(names.ZtomumuEfficiency);
  RooAbsArg* ZtomumuInvPurity = ws.arg(names.ZtomumuInvPurity);
  RooAbsArg* ZtoeeEfficiency = ws.arg(names.ZtoeeEfficiency);
  RooAbsArg* ZtoeeInvPurity = ws.arg(names.ZtoeeInvPurity);
  
  RooAbsArg* ZtoeeAcceptance = 
    getBetaConstraint(ws,"ZtoeeAcceptance_",abcd.ZtoeeAcceptanceName,
		      abcd.ZtoeeAcceptance,abcd.ZtoeeAcceptanceError,
		      names.observables,names.nuisances);

  RooAbsArg* ZtomumuAcceptance = 
	getBetaConstraint(ws,"ZtomumuAcceptance_",abcd.ZtomumuAcceptanceName,
			  abcd.ZtomumuAcceptance,abcd.ZtomumuAcceptanceError,
			  names.observables,names.nuisances);
  
  
  //BEN FIXME - this is not the correct implementation, but i'm doing it like this to copy OAK
  RooAbsArg*  ZtoeeSystematic = ws.arg("ZtoeeSystematic_"+abcd.ZtoeeSystematicName+"_Ratio");
  if(ZtoeeSystematic == NULL) 
    {
      ZtoeeSystematic = 
	getCorrelatedBetaPrimeConstraint(ws,"ZtoeeSystematic_",abcd.ZtoeeSystematicName,
					 abcd.ZtoeeSystematic,abcd.ZtoeeSystematicError,
					 names.observables,names.nuisances,
					 "ZtoNuNuSystematic_FIXME");
    }
  RooAbsArg*  ZtomumuSystematic = ws.arg("ZtomumuSystematic_"+abcd.ZtomumuSystematicName+"_Ratio");
  if(ZtomumuSystematic == NULL) 
    {
      ZtomumuSystematic = 
	getCorrelatedBetaPrimeConstraint(ws,"ZtomumuSystematic_",abcd.ZtomumuSystematicName,
					 abcd.ZtomumuSystematic,abcd.ZtomumuSystematicError,
					 names.observables,names.nuisances,
					 "ZtoNuNuSystematic_FIXME");
    }
  

  RooProduct* diMuonYield = (RooProduct*)ws.arg("diMuon_"+observed.diMuonName+"_Yield");//Assumes acceptance is only binned in zero or more dimensions of the count.
  if(diMuonYield == NULL) 
    {
      RooProduct diMuonYield_temp("diMuon_"+observed.diMuonName+"_Yield","diMuon_"+observed.diMuonName+"_Yield",RooArgSet(*Ztoll,*ZtomumuAcceptance,*ZtomumuEfficiency,*ZtomumuInvPurity,*ZtomumuSystematic));
      ws.import(diMuonYield_temp, RecycleConflictNodes());
      diMuonYield = (RooProduct*)ws.arg("diMuon_"+observed.diMuonName+"_Yield");//Assumes acceptance is only binned in zero or more dimensions of the count.
    }
  
  RooProduct* diElectronYield = (RooProduct*)ws.arg("diElectron_"+observed.diElectronName+"_Yield");//Assumes acceptance is only binned in zero or more dimensions of the count.
  if(diElectronYield == NULL) 
    {
      RooProduct diElectronYield_temp("diElectron_"+observed.diElectronName+"_Yield","diElectron_"+observed.diElectronName+"_Yield",RooArgSet(*Ztoll,*ZtoeeAcceptance,*ZtoeeEfficiency,*ZtoeeInvPurity,*ZtoeeSystematic));
      ws.import(diElectronYield_temp, RecycleConflictNodes());
      diElectronYield = (RooProduct*)ws.arg("diElectron_"+observed.diElectronName+"_Yield");//Assumes acceptance is only binned in zero or more dimensions of the count.
    }
  
  //-----Define Z->ll Poisson constraints
  
  RooPoissonLogEval* diMuonConstraint = (RooPoissonLogEval*)ws.arg("diMuon_"+observed.diMuonName+"_Constraint");
  if(diMuonConstraint == NULL) 
    {
      RooPoissonLogEval diMuonConstraint_temp("diMuon_"+observed.diMuonName+"_Constraint","diMuon_"+observed.diMuonName+"_Constraint",*diMuonCount,*diMuonYield);
      ws.import(diMuonConstraint_temp, RecycleConflictNodes());
      diMuonConstraint = (RooPoissonLogEval*)ws.arg("diMuon_"+observed.diMuonName+"_Constraint");
    }
  
  RooPoissonLogEval* diElectronConstraint = (RooPoissonLogEval*)ws.arg("diElectron_"+observed.diElectronName+"_Constraint");
  if(diElectronConstraint == NULL) 
    {
      RooPoissonLogEval diElectronConstraint_temp("diElectron_"+observed.diElectronName+"_Constraint","diElectron_"+observed.diElectronName+"_Constraint",*diElectronCount,*diElectronYield);
      ws.import(diElectronConstraint_temp, RecycleConflictNodes());
      diElectronConstraint = (RooPoissonLogEval*)ws.arg("diElectron_"+observed.diElectronName+"_Constraint");
    }
  
  //-----Define Z->nunu yield
  
  RooAbsArg* zeroLeptonZtoNuNubTagScaling = getBetaConstraint(ws,"zeroLeptonZtoNuNubTagScaling_",abcd.ZtoNuNubTagScalingName,
							      abcd.ZtoNuNubTagScaling,abcd.ZtoNuNubTagScalingError,
							      names.observables,names.nuisances);
  
  RooProduct zeroLeptonZtoNuNuYield(zeroLeptonName+"_ZtoNuNuYield",zeroLeptonName+"_ZtoNuNuYield",RooArgSet(*zeroLeptonZtoNuNubTagScaling,*ZtoNuNu));
  
  
  /////////////////////////////////////////////////////////////
  //  Non-QCD/Signal in Low Delta Phi_N region
  /////////////////////////////////////////////////////////////

  /*
  //Old non-QCD subtraction
  RooAbsArg* MCUncertainty = ws.arg(names.MCUncertainty);
  RooRealVar zeroLeptonLowDeltaPhiNMCCount(zeroLeptonLowDeltaPhiNName+"_MCCount",zeroLeptonLowDeltaPhiNName+"_MCCount",abcd.zeroLeptonLowDeltaPhiNMC);
  zeroLeptonLowDeltaPhiNMCCount.setConstant();
  RooProduct  zeroLeptonLowDeltaPhiNMCYield(zeroLeptonLowDeltaPhiNName+"_MCYield",zeroLeptonLowDeltaPhiNName+"_MCYield",RooArgSet(*MCUncertainty,zeroLeptonLowDeltaPhiNMCCount));
  */

  RooRealVar topWJetsLowDeltaPhiNOverZeroLeptonRatioMC(zeroLeptonLowDeltaPhiNName+"_TopWJetsLowDeltaPhiNOverZeroLeptonRatioMC",
						       zeroLeptonLowDeltaPhiNName+"_TopWJetsLowDeltaPhiNOverZeroLeptonRatioMC", 
						       abcd.topWJetsLowDeltaPhiNOverZeroLeptonRatioMC);
  topWJetsLowDeltaPhiNOverZeroLeptonRatioMC.setConstant();
  RooProduct zeroLeptonLowDeltaPhiNTopWJetsYield(zeroLeptonLowDeltaPhiNName+"_TopWJetsYield",zeroLeptonLowDeltaPhiNName+"_TopWJetsYield",RooArgSet(topWJetsLowDeltaPhiNOverZeroLeptonRatioMC,zeroLeptonTopWJetsYield));

  RooRealVar ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC(zeroLeptonLowDeltaPhiNName+"_ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC",
						       zeroLeptonLowDeltaPhiNName+"_ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC", 
						       abcd.ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC);
  ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC.setConstant();
  RooProduct zeroLeptonLowDeltaPhiNZtoNuNuYield(zeroLeptonLowDeltaPhiNName+"_ZtoNuNuYield",zeroLeptonLowDeltaPhiNName+"_ZtoNuNuYield",RooArgSet(ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC,zeroLeptonZtoNuNuYield));

  RooAddition zeroLeptonLowDeltaPhiNNonQCDYield(zeroLeptonLowDeltaPhiNName+"_NonQCDYield",zeroLeptonLowDeltaPhiNName+"_NonQCDYield",RooArgSet(zeroLeptonLowDeltaPhiNTopWJetsYield,zeroLeptonLowDeltaPhiNZtoNuNuYield));

 
  ////////////////////////////////////////////////
  // Constrain zero lepton yields
  ////////////////////////////////////////////////

  // Get signal yields from workspace
  
  RooAbsReal* zeroLeptonSignalYield = ws.function(zeroLeptonName+"_SignalYield");
  RooAbsReal* zeroLeptonLowDeltaPhiNSignalYield = ws.function(zeroLeptonLowDeltaPhiNName+"_SignalYield");
  RooAbsReal* oneLeptonSignalYield = ws.function(oneLeptonName+"_SignalYield");
  assert(zeroLeptonSignalYield && zeroLeptonLowDeltaPhiNSignalYield && oneLeptonSignalYield);
  
  // Setup yields in all bins
  
  RooAddition zeroLeptonYieldSum(zeroLeptonName+"_YieldSum",zeroLeptonName+"_YieldSum",RooArgSet(*zeroLeptonSignalYield,zeroLeptonZtoNuNuYield,zeroLeptonTopWJetsYield,zeroLeptonQCDYield));
  double topGuess = observed.zeroLepton - zeroLeptonZtoNuNuYield.getVal() - zeroLeptonQCDYield.getVal();
  if(topGuess > 0 ) zeroLeptonTopWJetsGuess += topGuess;
  //RooAddition zeroLeptonLowDeltaPhiNYieldSum(zeroLeptonLowDeltaPhiNName+"_YieldSum",zeroLeptonLowDeltaPhiNName+"_YieldSum",RooArgSet(*zeroLeptonLowDeltaPhiNSignalYield,zeroLeptonLowDeltaPhiNQCDYield,zeroLeptonLowDeltaPhiNMCYield));//Old non-QCD subtraction
  RooAddition zeroLeptonLowDeltaPhiNYieldSum(zeroLeptonLowDeltaPhiNName+"_YieldSum",zeroLeptonLowDeltaPhiNName+"_YieldSum",RooArgSet(*zeroLeptonLowDeltaPhiNSignalYield,zeroLeptonLowDeltaPhiNQCDYield,zeroLeptonLowDeltaPhiNNonQCDYield));
  RooAddition oneLeptonYieldSum(oneLeptonName+"_YieldSum",oneLeptonName+"_YieldSum",RooArgSet(*oneLeptonSignalYield,oneLeptonTopWJetsYield));
  
  
  // Setup trigger efficiencies

  RooRealVar* zeroLeptonTriggerEfficiency = 0;
  RooRealVar* zeroLeptonLowDeltaPhiNTriggerEfficiency = 0; 
  RooRealVar* oneLeptonTriggerEfficiency = 0;
  
  if (options.skipTriggerEfficiency == true)
    {
      RooRealVar zeroLeptonTriggerEfficiency_temp = RooRealVar("zeroLeptonTriggerEfficiency_"+binName, "zeroLeptonTriggerEfficiency_"+binName, 1.0);
      RooRealVar zeroLeptonLowDeltaPhiNTriggerEfficiency_temp = RooRealVar("zeroLeptonLowDeltaPhiNTriggerEfficiency_"+binName, "zeroLeptonLowDeltaPhiNTriggerEfficiency_"+binName, 1.0);
      RooRealVar oneLeptonTriggerEfficiency_temp = RooRealVar("oneLeptonTriggerEfficiency_"+binName, "oneLeptonTriggerEfficiency_"+binName, 1.0);
      
      zeroLeptonTriggerEfficiency_temp.setConstant();
      zeroLeptonLowDeltaPhiNTriggerEfficiency_temp.setConstant();
      oneLeptonTriggerEfficiency_temp.setConstant();
      
      ws.import( zeroLeptonTriggerEfficiency_temp );
      ws.import( zeroLeptonLowDeltaPhiNTriggerEfficiency_temp ); 
      ws.import( oneLeptonTriggerEfficiency_temp );
      
      zeroLeptonTriggerEfficiency = ws.var( zeroLeptonTriggerEfficiency_temp.GetName() );
      zeroLeptonLowDeltaPhiNTriggerEfficiency = ws.var(  zeroLeptonLowDeltaPhiNTriggerEfficiency_temp.GetName() );
      oneLeptonTriggerEfficiency = ws.var( oneLeptonTriggerEfficiency_temp.GetName() );
    }
  else 
    {
      zeroLeptonTriggerEfficiency = (RooRealVar*) 
	getBetaConstraint(ws,"zeroLeptonTriggerEfficiency_",binName,
			  abcd.zeroLeptonTriggerEfficiency,abcd.zeroLeptonTriggerEfficiencyError,
			  names.observables,names.nuisances);
      
      zeroLeptonLowDeltaPhiNTriggerEfficiency = (RooRealVar*) 
	getBetaConstraint(ws,"zeroLeptonLowDeltaPhiNTriggerEfficiency_",binName,
			  abcd.zeroLeptonLowDeltaPhiNTriggerEfficiency,abcd.zeroLeptonLowDeltaPhiNTriggerEfficiencyError,
			  names.observables,names.nuisances);
      
      oneLeptonTriggerEfficiency = (RooRealVar*) 
	getBetaConstraint(ws,"oneLeptonTriggerEfficiency_",binName,
			  abcd.oneLeptonTriggerEfficiency,abcd.oneLeptonTriggerEfficiencyError,
			  names.observables,names.nuisances);
    }
  
  // Total Yields in bins
  
  RooProduct zeroLeptonYield(zeroLeptonName+"_Yield",zeroLeptonName+"_Yield",RooArgSet(*zeroLeptonTriggerEfficiency,zeroLeptonYieldSum));
  RooProduct zeroLeptonLowDeltaPhiNYield(zeroLeptonLowDeltaPhiNName+"_Yield",zeroLeptonLowDeltaPhiNName+"_Yield",RooArgSet(*zeroLeptonLowDeltaPhiNTriggerEfficiency,zeroLeptonLowDeltaPhiNYieldSum));
  RooProduct oneLeptonYield(oneLeptonName+"_Yield",oneLeptonName+"_Yield",RooArgSet(*oneLeptonTriggerEfficiency,oneLeptonYieldSum));
  
  // Define poisson constraints
  
  RooPoissonLogEval zeroLeptonConstraint(zeroLeptonName+"_Constraint",zeroLeptonName+"_Constraint",zeroLeptonCount,zeroLeptonYield);
  RooPoissonLogEval zeroLeptonLowDeltaPhiNConstraint(zeroLeptonLowDeltaPhiNName+"_Constraint",zeroLeptonLowDeltaPhiNName+"_Constraint",zeroLeptonLowDeltaPhiNCount,zeroLeptonLowDeltaPhiNYield);
  RooPoissonLogEval oneLeptonConstraint(oneLeptonName+"_Constraint",oneLeptonName+"_Constraint",oneLeptonCount,oneLeptonYield);    
  
  // Cleanup
  
  ws.import(zeroLeptonConstraint, RecycleConflictNodes());
  ws.import(zeroLeptonLowDeltaPhiNConstraint, RecycleConflictNodes());
  ws.import(oneLeptonConstraint, RecycleConflictNodes());
    
  return true;
  
}

void makeGuess(const likelihoodOptions options, RooWorkspace& ws , TString& binName , channels& observed , abcdBinParameters& abcd)
{

  TString scaling = "lowDeltaPhiNScaling_";
  if(options.qcdMethod == "singleScaleWithCorrections") scaling +=  "all";
  else if (options.qcdMethod == "htDependent") scaling +=  abcd.lowDeltaPhiNScalingName;
  else assert(0);
  RooRealVar* lowDeltaPhiNScaling = ws.var(scaling);
  RooRealVar* zeroLeptonQCDClosure = ws.var("zeroLeptonQCDClosure_"+binName);
  RooRealVar* zeroLeptonLowDeltaPhiNCount = ws.var("zeroLeptonLowDeltaPhiN_"+binName+"_Count");
  RooAbsReal* zeroLeptonLowDeltaPhiNSignalYield = ws.function("zeroLeptonLowDeltaPhiN_"+binName+"_SignalYield");
  RooRealVar* zeroLeptonLowDeltaPhiNQCDYield = ws.var("zeroLeptonLowDeltaPhiN_"+binName+"_QCDYield");
  RooAbsReal* zeroLeptonLowDeltaPhiNNonQCDYield = ws.function("zeroLeptonLowDeltaPhiN_"+binName+"_NonQCDYield");
  RooAbsReal* zeroLeptonLowDeltaPhiNZtoNuNuYield = ws.function("zeroLeptonLowDeltaPhiN_"+binName+"_ZtoNuNuYield");
  RooAbsReal* zeroLeptonLowDeltaPhiNTopWJetsYield = ws.function("zeroLeptonLowDeltaPhiN_"+binName+"_TopWJetsYield");
  RooAbsReal* zeroLeptonLowDeltaPhiNYield = ws.function("zeroLeptonLowDeltaPhiN_"+binName+"_Yield");
  RooRealVar* zeroLeptonTopWJetsClosure = ws.var("zeroLeptonTopWJetsClosure_"+binName);
  RooAbsReal* zeroLeptonTopWJetsYield = ws.function("zeroLepton_"+binName+"_TopWJetsYield");
  RooRealVar* zeroLeptonCount = ws.var("zeroLepton_"+binName+"_Count");
  RooAbsReal* zeroLeptonSignalYield = ws.function("zeroLepton_"+binName+"_SignalYield");
  RooAbsReal* zeroLeptonZtoNuNuYield = ws.function("zeroLepton_"+binName+"_ZtoNuNuYield");
  RooAbsReal* zeroLeptonQCDYield = ws.function("zeroLepton_"+binName+"_QCDYield");
  RooRealVar* singleLeptonScaling = ws.var("singleLeptonScaling");
  RooAbsReal* oneLeptonTopWJetsYield = ws.function("oneLepton_"+binName+"_TopWJetsYield");
  RooAbsReal* zeroLeptonYield = ws.function("zeroLepton_"+binName+"_Yield");
  RooRealVar* topWJetsLowDeltaPhiNOverZeroLeptonRatioMC = ws.var("zeroLeptonLowDeltaPhiN_"+binName+"_TopWJetsLowDeltaPhiNOverZeroLeptonRatioMC");
  assert(lowDeltaPhiNScaling);
  assert(zeroLeptonQCDClosure);
  assert(zeroLeptonLowDeltaPhiNCount);
  assert(zeroLeptonLowDeltaPhiNSignalYield);
  assert(zeroLeptonLowDeltaPhiNNonQCDYield);
  assert(zeroLeptonLowDeltaPhiNQCDYield);
  assert(zeroLeptonLowDeltaPhiNZtoNuNuYield);
  assert(zeroLeptonLowDeltaPhiNTopWJetsYield);
  assert(zeroLeptonLowDeltaPhiNYield);
  assert(zeroLeptonTopWJetsClosure);
  assert(zeroLeptonTopWJetsYield);
  assert(zeroLeptonCount);
  assert(zeroLeptonSignalYield);
  assert(zeroLeptonZtoNuNuYield);
  assert(zeroLeptonQCDYield);
  assert(singleLeptonScaling);
  assert(oneLeptonTopWJetsYield);
  assert(zeroLeptonYield);
  
  /*
  //simple guesses
  //zeroLeptonLowDeltaPhiNQCDYield.setVal(min(1e-5,observed.zeroLeptonLowDeltaPhiN-zeroLeptonLowDeltaPhiNTopWJetsYield.getVal()-zeroLeptonLowDeltaPhiNZtoNuNuYield.getVal()));  // tried putting this in makeOneBin
  zeroLeptonQCDClosure->setVal( (zeroLeptonQCDYield->getVal()/lowDeltaPhiNScaling->getVal())*1.0/(zeroLeptonLowDeltaPhiNCount->getVal() - zeroLeptonLowDeltaPhiNSignalYield->getVal() - zeroLeptonLowDeltaPhiNNonQCDYield->getVal()) );
  zeroLeptonTopWJetsClosure->setVal( (zeroLeptonCount->getVal() - zeroLeptonSignalYield->getVal() - zeroLeptonZtoNuNuYield->getVal() - zeroLeptonQCDYield->getVal())/(singleLeptonScaling->getVal() * oneLeptonTopWJetsYield->getVal()) );
  // Printout
  cout << "*******************************************************************" << endl;
  cout << binName << endl;		
  cout << "zeroLeptonCount = " << zeroLeptonCount->getVal() << ", zeroLeptonYield = " << zeroLeptonYield->getVal() << endl;
  cout << "zeroLeptonLowDeltaPhiNCount = " << zeroLeptonLowDeltaPhiNCount->getVal() << ", zeroLeptonLowDeltaPhiNYield = " << zeroLeptonLowDeltaPhiNYield->getVal() << endl;
  cout << "*******************************************************************" << endl;
  */
  
  /*
  //Solving 2 eqns -- qcd closure and qcd ldp yield
  zeroLeptonQCDClosure->setVal( (-zeroLeptonCount->getVal()+zeroLeptonSignalYield->getVal()+zeroLeptonZtoNuNuYield->getVal()+zeroLeptonTopWJetsYield->getVal() ) / ( lowDeltaPhiNScaling->getVal() * (-zeroLeptonLowDeltaPhiNCount->getVal() + zeroLeptonLowDeltaPhiNSignalYield->getVal()+zeroLeptonLowDeltaPhiNNonQCDYield->getVal()  ) ) );	 
  zeroLeptonLowDeltaPhiNQCDYield->setVal( zeroLeptonLowDeltaPhiNCount->getVal() - zeroLeptonLowDeltaPhiNSignalYield->getVal() - zeroLeptonLowDeltaPhiNNonQCDYield->getVal()  );
  */
  
  //Solving 2 eqns -- TopWJets closure and qcd ldp yield
  double thisZeroLeptonTopWJetsClosure = 
    -( -zeroLeptonCount->getVal()+zeroLeptonSignalYield->getVal()+zeroLeptonZtoNuNuYield->getVal()+ lowDeltaPhiNScaling->getVal()*(zeroLeptonLowDeltaPhiNCount->getVal()
       -zeroLeptonLowDeltaPhiNSignalYield->getVal() -zeroLeptonLowDeltaPhiNZtoNuNuYield->getVal())*zeroLeptonQCDClosure->getVal()  ) 
    / ( singleLeptonScaling->getVal() * oneLeptonTopWJetsYield->getVal() - 
	lowDeltaPhiNScaling->getVal()*singleLeptonScaling->getVal() * topWJetsLowDeltaPhiNOverZeroLeptonRatioMC->getVal() * oneLeptonTopWJetsYield->getVal()*zeroLeptonQCDClosure->getVal() )  ;
  
  double thisZeroLeptonLowDeltaPhiNQCDYield =  
    (-zeroLeptonLowDeltaPhiNCount->getVal()+zeroLeptonLowDeltaPhiNSignalYield->getVal() + topWJetsLowDeltaPhiNOverZeroLeptonRatioMC->getVal() *(zeroLeptonCount->getVal() 
     -zeroLeptonSignalYield->getVal() - zeroLeptonZtoNuNuYield->getVal() ) + zeroLeptonLowDeltaPhiNZtoNuNuYield->getVal() )
    /(-1.0 + lowDeltaPhiNScaling->getVal()* topWJetsLowDeltaPhiNOverZeroLeptonRatioMC->getVal()*zeroLeptonQCDClosure->getVal()) ;
  
  if ((thisZeroLeptonTopWJetsClosure > zeroLeptonTopWJetsClosure->getMin()) && (thisZeroLeptonTopWJetsClosure < zeroLeptonTopWJetsClosure->getMax())
      && (thisZeroLeptonLowDeltaPhiNQCDYield > zeroLeptonLowDeltaPhiNQCDYield->getMin()) && (thisZeroLeptonLowDeltaPhiNQCDYield < zeroLeptonLowDeltaPhiNQCDYield->getMax())) 
    {
      zeroLeptonTopWJetsClosure->setVal( thisZeroLeptonTopWJetsClosure );
      zeroLeptonLowDeltaPhiNQCDYield->setVal( thisZeroLeptonLowDeltaPhiNQCDYield );
    }

}


void makeUnderlyingLikelihood(const likelihoodOptions options, RooWorkspace& ws ,allBinNames& names, allBins& numbers, double luminosityInput)
{
  ws.defineSet("observables","");
  names.observables = "observables";

  ws.defineSet("nuisances","");
  names.nuisances = "nuisances";

  //Universal parameters

  /*
  //Old non-QCD subtraction
  RooAbsArg* MCUncertainty = 
    getBetaPrimeConstraint(ws,"MCUncertainty","",
			   1.0,numbers.MCUncertainty,
			   names.observables,names.nuisances);
  names.MCUncertainty = MCUncertainty->GetName();
  */
  
  RooRealVar luminosity("luminosity","luminosity",luminosityInput,0.0,1e2);
  luminosity.setConstant();//should eventually have error
  ws.import(luminosity);
  
  RooRealVar signalCrossSection("signalCrossSection","signalCrossSection",0.0,0.0,1e5);
  ws.import(signalCrossSection);
  names.signalCrossSection = signalCrossSection.GetName();

  names.signalUncertainty = "signalUncertainty";

  RooRealVar singleLeptonScaling("singleLeptonScaling","singleLeptonScaling",0.0,1e3);
  ws.import(singleLeptonScaling);
  ws.extendSet(names.nuisances,singleLeptonScaling.GetName());
  names.singleLeptonScaling = singleLeptonScaling.GetName();


  //Objects for Z->invisible background:

  RooRealVar ZtollOverZtoNuNuRatio("ZtollOverZtoNuNuRatio","ZtollOverZtoNuNuRatio",numbers.ZtollOverZtoNuNuRatio);
  ZtollOverZtoNuNuRatio.setConstant();
  ws.import(ZtollOverZtoNuNuRatio);
  names.ZtollOverZtoNuNuRatio = ZtollOverZtoNuNuRatio.GetName();

  RooAbsArg* ZtomumuInvPurity = 
    getInverseBetaConstraint(ws,"ZtomumuInvPurity","",
			     numbers.ZtomumuPurity,numbers.ZtomumuPurityError,
			     names.observables,names.nuisances);
  names.ZtomumuInvPurity = ZtomumuInvPurity->GetName();
  
  RooAbsArg* ZtomumuEfficiency = 
    getBetaConstraint(ws,"ZtomumuEfficiency","",
		      numbers.ZtomumuEfficiency,numbers.ZtomumuEfficiencyError,
		      names.observables,names.nuisances);
  names.ZtomumuEfficiency = ZtomumuEfficiency->GetName();
  
  RooAbsArg* ZtoeeInvPurity = 
    getInverseBetaConstraint(ws,"ZtoeeInvPurity","",
			     numbers.ZtoeePurity,numbers.ZtoeePurityError,
			     names.observables,names.nuisances);
  names.ZtoeeInvPurity = ZtoeeInvPurity->GetName();

  RooAbsArg* ZtoeeEfficiency = 
    getBetaConstraint(ws,"ZtoeeEfficiency","",
		      numbers.ZtoeeEfficiency,numbers.ZtoeeEfficiencyError,
		      names.observables,names.nuisances);
  names.ZtoeeEfficiency = ZtoeeEfficiency->GetName();

}



void setupSignalModelMR(const TString binName, const TString binFileNameInsideSignalMR, map<TString,map<TString,int> >& insideSignalMR, 
			const TString binFileNameOutsideSignalMR, map<TString,map<TString,int> >& outsideSignalMR )
{
  
  cout << "Reading in MR signal inputs for " << binName << endl;
    
  ifstream insideFile;
  
  cout << "getting the file: " << binFileNameInsideSignalMR << endl;
  
  insideFile.open(binFileNameInsideSignalMR.Data(),fstream::in);
  assert(insideFile.is_open());
  
  string fileLine;
  
  TString index;
  int value;
  
  map<TString,int> thisInsideSignalMR;

  while(!insideFile.eof())
    {

      getline(insideFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken nameAndNumber(thisLine," ");
      nameAndNumber.NextToken();
      index = nameAndNumber;
      if(index == "") continue;
      nameAndNumber.NextToken();
      value = nameAndNumber.Atoi();
      cout << index << " : " << value << endl;
      
      if(index == "oneTightMu_Theta1_SignalRawCount"       ) { thisInsideSignalMR.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta1", value) ); }
      else if(index == "oneLooseLep_Theta1_SignalRawCount" ) { thisInsideSignalMR.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta1", value) ); }
      else if(index == "oneTightMu_Theta2_SignalRawCount"  ) { thisInsideSignalMR.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta2", value) ); } 
      else if(index == "oneLooseLep_Theta2_SignalRawCount" ) { thisInsideSignalMR.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta2", value) ); }
      else if(index == "oneTightMu_Theta3_SignalRawCount"  ) { thisInsideSignalMR.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta3", value) ); }
      else if(index == "oneLooseLep_Theta3_SignalRawCount" ) { thisInsideSignalMR.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta3", value) ); }
      else if(index == "oneTightMu_Theta4_SignalRawCount"  ) { thisInsideSignalMR.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta4", value) ); }
      else if(index == "oneLooseLep_Theta4_SignalRawCount" ) { thisInsideSignalMR.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta4", value) ); }
      else if(index == "oneTightMu_Theta5_SignalRawCount"  ) { thisInsideSignalMR.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta5", value) ); }
      else if(index == "oneLooseLep_Theta5_SignalRawCount" ) { thisInsideSignalMR.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta5", value) ); }
      else if(index == "twoTightMu_SignalRawCount"         ) { thisInsideSignalMR.insert( pair<TString,int>("twoTightMu_"+binName, value) ); }
      else if(index == "twoLooseLep_SignalRawCount"        ) { thisInsideSignalMR.insert( pair<TString,int>("twoLooseLep_"+binName, value) ); }
      else {assert(0);}
    }
  insideFile.close();
  
  //for ( map<TString,int>::iterator it=thisInsideSignalMR.begin() ; it != thisInsideSignalMR.end(); it++ ) cout << (*it).first << " => " << (*it).second << endl;
  insideSignalMR.insert( pair<TString, map<TString,int> >(binName, thisInsideSignalMR) );

  //Now the outside values
  
  ifstream outsideFile;
  
  cout << "getting the file: " << binFileNameOutsideSignalMR << endl;
  
  outsideFile.open(binFileNameOutsideSignalMR.Data(),fstream::in);
  assert(outsideFile.is_open());
  
  TString svalue;

  sOutsideSignalMR holdOutsideSignalMR;

  map<TString,int> thisOutsideSignalMR;

  while(!outsideFile.eof())
    {

      getline(outsideFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken nameAndNumber(thisLine," ");
      nameAndNumber.NextToken();
      index = nameAndNumber;
      if(index == "") continue;
      nameAndNumber.NextToken();
      svalue = nameAndNumber;
      if(svalue.IsDigit() == false) 
	{
	  //hack to get rid of potential whitespace in strings
	  string str = (string)svalue;
	  for(unsigned int i=0; i<str.length(); i++) 
	    {
	      if(str[i] == '\t') str.erase(i,1);
	      if(str[i] == ' ') str.erase(i,1);
	    }
	  svalue = (TString)str;
	}
      
      cout << index << " : " << svalue << endl;
      
      if(index == "oneTightMu_Outside_SignalRawCountName"      ) { holdOutsideSignalMR.oneTightMu_Outside_SignalRawCountName = svalue; }
      else if(index == "oneTightMu_Outside_SignalRawCount"     ) { holdOutsideSignalMR.oneTightMu_Outside_SignalRawCount = svalue.Atoi(); }
      else if(index == "twoTightMu_Outside_SignalRawCountName" ) { holdOutsideSignalMR.twoTightMu_Outside_SignalRawCountName = svalue; }
      else if(index == "twoTightMu_Outside_SignalRawCount"     ) { holdOutsideSignalMR.twoTightMu_Outside_SignalRawCount = svalue.Atoi(); }
      else{ assert(0); }
    }
  outsideFile.close();

  thisOutsideSignalMR.insert( pair<TString,int>("oneTightMu_"+holdOutsideSignalMR.oneTightMu_Outside_SignalRawCountName+"_Outside", holdOutsideSignalMR.oneTightMu_Outside_SignalRawCount) );
  thisOutsideSignalMR.insert( pair<TString,int>("twoTightMu_"+holdOutsideSignalMR.twoTightMu_Outside_SignalRawCountName+"_Outside", holdOutsideSignalMR.twoTightMu_Outside_SignalRawCount) );

  //for ( map<TString,int>::iterator it=thisOutsideSignalMR.begin() ; it != thisOutsideSignalMR.end(); it++ ) cout << (*it).first << " => " << (*it).second << endl;
  outsideSignalMR.insert( pair<TString, map<TString,int> >(binName, thisOutsideSignalMR) );
  
}



void setupSignalModel(const likelihoodOptions options, RooWorkspace& ws , vector<TString> binNames, allBinNames& names, TString signalModelFilesPath,
		      map<TString,map<TString,int> >& insideSignalMR, map<TString,map<TString,int> >& outsideSignalMR)
{
  
  int signalModelLineNumber = 0;      
  
  TString signalModelFileName = signalModelFilesPath;
  signalModelFileName += "/setupSignalNominal.dat";
  
  cout << "getting the file: " << signalModelFileName << endl;
  
  ifstream signalFile;
  
  signalFile.open(signalModelFileName.Data(),fstream::in);
  assert(signalFile.is_open());
  
  int thisLineNumber = 0;
  string fileLine;
  
  while( thisLineNumber < signalModelLineNumber && !signalFile.eof() )
    {
      thisLineNumber++;
      getline(signalFile,fileLine);
    }
  
  double m0,m12,susyGenerated;
  double susyInBin = 0;
  
  signalFile>>m0;
  cout << "m0 : " << m0 << endl;
  signalFile>>m12;
  cout << "m12 : " << m12 << endl;
  signalFile>>susyGenerated;
  cout << "SUSY Generated : " << susyGenerated << endl;
  
  RooArgSet signalRawYields("signalRawYields");
  
  //Loop over bins to get raw counts and make raw yields and constraints

  for(vector<TString>::iterator thisBin = binNames.begin() ; thisBin != binNames.end() ; thisBin++)
    {
      TString binName = *thisBin;
      
      yields thisSignal,thisSignalError;
      double valueHolder;
      signalFile>>valueHolder ; thisSignal.zeroLepton = valueHolder ;
      signalFile>>valueHolder ; thisSignalError.zeroLepton = valueHolder ;
      signalFile>>valueHolder ; thisSignal.zeroLeptonLowDeltaPhiN = valueHolder ;
      signalFile>>valueHolder ; thisSignalError.zeroLeptonLowDeltaPhiN = valueHolder ;
      signalFile>>valueHolder ; thisSignal.oneLepton = valueHolder ;
      signalFile>>valueHolder ; thisSignalError.oneLepton = valueHolder ;
      
      cout << "For selection " << *thisBin << " signal is " << endl;
      cout << "zero lepton bin               : " << thisSignal.zeroLepton             << endl ;
      cout << "zero lepton low delta phi bin : " << thisSignal.zeroLeptonLowDeltaPhiN << endl ;
      cout << "one lepton bin                : " << thisSignal.oneLepton              << endl ;
      //cout << "For selection " << *thisBin << " signal error is " << endl;
      //cout << "zero lepton bin               : " << thisSignalError.zeroLepton             << endl ;
      //cout << "zero lepton low delta phi bin : " << thisSignalError.zeroLeptonLowDeltaPhiN << endl ;
      //cout << "one lepton bin                : " << thisSignalError.oneLepton              << endl ;
      
      //Zero lepton and low delta phi n counts used in nominal and met reweighting modes

      TString zeroLeptonName("zeroLepton_");
      zeroLeptonName+=binName;
      TString zeroLeptonLowDeltaPhiNName("zeroLeptonLowDeltaPhiN_");
      zeroLeptonLowDeltaPhiNName+=binName;
      
      RooRealVar zeroLeptonSignalRawCount(zeroLeptonName+"_SignalRawCount", zeroLeptonName+"_SignalRawCount", thisSignal.zeroLepton); 
      RooRealVar zeroLeptonLowDeltaPhiNSignalRawCount(zeroLeptonLowDeltaPhiNName+"_SignalRawCount", zeroLeptonLowDeltaPhiNName+"_SignalRawCount", thisSignal.zeroLeptonLowDeltaPhiN); 
      
      zeroLeptonSignalRawCount.setConstant();
      zeroLeptonLowDeltaPhiNSignalRawCount.setConstant();
      
      RooRealVar zeroLeptonSignalRawYield(zeroLeptonName+"_SignalRawYield", zeroLeptonName+"_SignalRawYield", thisSignal.zeroLepton, 0, 1e5); 
      RooRealVar zeroLeptonLowDeltaPhiNSignalRawYield(zeroLeptonLowDeltaPhiNName+"_SignalRawYield", zeroLeptonLowDeltaPhiNName+"_SignalRawYield", thisSignal.zeroLeptonLowDeltaPhiN, 0, 1e5); 
      
      RooPoissonLogEval zeroLeptonSignalRawConstraint(zeroLeptonName+"_SignalRawConstraint", zeroLeptonName+"_SignalRawConstraint", zeroLeptonSignalRawCount, zeroLeptonSignalRawYield);
      ws.import(zeroLeptonSignalRawConstraint, RecycleConflictNodes());
      ws.extendSet(names.observables,zeroLeptonSignalRawCount.GetName());
      ws.extendSet(names.nuisances,zeroLeptonSignalRawYield.GetName());

      RooPoissonLogEval zeroLeptonLowDeltaPhiNSignalRawConstraint(zeroLeptonLowDeltaPhiNName+"_SignalRawConstraint", zeroLeptonLowDeltaPhiNName+"_SignalRawConstraint", zeroLeptonLowDeltaPhiNSignalRawCount, zeroLeptonLowDeltaPhiNSignalRawYield);
      ws.import(zeroLeptonLowDeltaPhiNSignalRawConstraint, RecycleConflictNodes());
      ws.extendSet(names.observables,zeroLeptonLowDeltaPhiNSignalRawCount.GetName());
      ws.extendSet(names.nuisances,zeroLeptonLowDeltaPhiNSignalRawYield.GetName());      

      signalRawYields.add(*(ws.var(zeroLeptonName+"_SignalRawYield")));
      signalRawYields.add(*(ws.var(zeroLeptonLowDeltaPhiNName+"_SignalRawYield")));

      susyInBin += thisSignal.zeroLepton;
      susyInBin += thisSignal.zeroLeptonLowDeltaPhiN;
  
      //one lepton counts 
      
      if(options.TopWJetsMethod == "ABCD")
	{
	  TString oneLeptonName("oneLepton_");
	  oneLeptonName+=binName;

	  RooRealVar oneLeptonSignalRawCount(oneLeptonName+"_SignalRawCount", oneLeptonName+"_SignalRawCount", thisSignal.oneLepton); 
	  oneLeptonSignalRawCount.setConstant();
	  
	  RooRealVar oneLeptonSignalRawYield(oneLeptonName+"_SignalRawYield", oneLeptonName+"_SignalRawYield", thisSignal.oneLepton, 0, 1e5); 
	  
	  RooPoissonLogEval oneLeptonSignalRawConstraint(oneLeptonName+"_SignalRawConstraint", oneLeptonName+"_SignalRawConstraint", oneLeptonSignalRawCount, oneLeptonSignalRawYield);
	  ws.import(oneLeptonSignalRawConstraint, RecycleConflictNodes());
	  ws.extendSet(names.observables,oneLeptonSignalRawCount.GetName());
	  ws.extendSet(names.nuisances,oneLeptonSignalRawYield.GetName());

	  signalRawYields.add(*(ws.var(oneLeptonName+"_SignalRawYield")));
	  susyInBin += thisSignal.oneLepton;
	}
      else if(options.TopWJetsMethod == "metReweighting") 
	{
	  //inside
	  map<TString,int> thisInsideSignalMR = insideSignalMR[binName];
	  for ( map<TString,int>::iterator it=thisInsideSignalMR.begin() ; it != thisInsideSignalMR.end(); it++ ) 
	    { 
	      TString signalName = (*it).first;
	      int signalRawCount = (*it).second;
	      
	      RooRealVar rooSignalRawCount(signalName+"_SignalRawCount", signalName+"_SignalRawCount", signalRawCount);
	      rooSignalRawCount.setConstant();

	      RooRealVar rooSignalRawYield(signalName+"_SignalRawYield", signalName+"_SignalRawYield", signalRawCount, 0, 1e5);

	      RooPoissonLogEval rooSignalRawConstraint(signalName+"_SignalRawConstraint", signalName+"_SignalRawConstraint", rooSignalRawCount, rooSignalRawYield);
	      ws.import(rooSignalRawConstraint, RecycleConflictNodes());
	      ws.extendSet(names.observables, rooSignalRawCount.GetName());
	      ws.extendSet(names.nuisances, rooSignalRawYield.GetName());

	      signalRawYields.add( *( ws.var(rooSignalRawYield.GetName()) ) );
	      susyInBin += signalRawCount;
	    }
	  
	  //outside
	  map<TString,int> thisOutsideSignalMR = outsideSignalMR[binName];
	  for ( map<TString,int>::iterator it=thisOutsideSignalMR.begin() ; it != thisOutsideSignalMR.end(); it++ ) 
	    { 
	      TString signalName = (*it).first;
	      int signalRawCount = (*it).second;
	      
	      //Continue if already in workspace 
	      if( ws.var(signalName+"_SignalRawCount") != NULL ) continue; 
	      
	      RooRealVar rooSignalRawCount(signalName+"_SignalRawCount", signalName+"_SignalRawCount", signalRawCount);
	      rooSignalRawCount.setConstant();
	      
	      RooRealVar rooSignalRawYield(signalName+"_SignalRawYield", signalName+"_SignalRawYield", signalRawCount, 0, 1e5);
	      
	      RooPoissonLogEval rooSignalRawConstraint(signalName+"_SignalRawConstraint", signalName+"_SignalRawConstraint", rooSignalRawCount, rooSignalRawYield);
	      ws.import(rooSignalRawConstraint, RecycleConflictNodes());
	      ws.extendSet(names.observables, rooSignalRawCount.GetName());
	      ws.extendSet(names.nuisances, rooSignalRawYield.GetName());

	      signalRawYields.add( *( ws.var(rooSignalRawYield.GetName()) ) );
	      susyInBin += signalRawCount;
	    }
	}
      else { assert(0); }
      
    }	  
  
  signalFile.close();
  
  //Make count/yield/constraint for events not in a bin and then define the total yield

  RooRealVar noBinSignalRawCount("noBinSignalRawCount", "noBinSignalRawCount", susyGenerated-susyInBin);
  noBinSignalRawCount.setConstant();
  
  RooRealVar noBinSignalRawYield("noBinSignalRawYield", "noBinSignalRawYield", susyGenerated-susyInBin, 0, 1e5);
  signalRawYields.add(noBinSignalRawYield);
  
  RooPoissonLogEval noBinSignalRawConstraint("noBinSignalRawConstraint", "noBinSignalRawConstraint", noBinSignalRawCount, noBinSignalRawYield);
  ws.import(noBinSignalRawConstraint, RecycleConflictNodes());
  ws.extendSet(names.observables,noBinSignalRawCount.GetName());
  ws.extendSet(names.nuisances, noBinSignalRawYield.GetName());
  
  cout <<"total raw signal yield" << endl;
  signalRawYields.Print();
  
  RooAddition totalSignalRawYield("totalSignalRawYield", "totalSignalRawYield", signalRawYields); 
    
  RooRealVar* luminosity = ws.var("luminosity");
  RooRealVar* signalCrossSection =  ws.var(names.signalCrossSection);
  
  //Loop over bins to make yields (rawYield/total*lumi*xsec)

  for(vector<TString>::iterator thisBin = binNames.begin() ; thisBin != binNames.end() ; thisBin++)
    {
      TString binName = *thisBin;
      
      //Make zero lepton and low delta phi n yields
      
      TString zeroLeptonName("zeroLepton_");
      zeroLeptonName+=binName;
      TString zeroLeptonLowDeltaPhiNName("zeroLeptonLowDeltaPhiN_");
      zeroLeptonLowDeltaPhiNName+=binName;
        
      RooRealVar* zeroLeptonSignalRawYield = ws.var(zeroLeptonName+"_SignalRawYield");
      RooRealVar* zeroLeptonLowDeltaPhiNSignalRawYield = ws.var(zeroLeptonLowDeltaPhiNName+"_SignalRawYield");
      
      RooRatio zeroLeptonSignalYieldFraction(zeroLeptonName+"_SignalYieldFraction", zeroLeptonName+"_SignalYieldFraction", *zeroLeptonSignalRawYield, totalSignalRawYield);
      RooRatio zeroLeptonLowDeltaPhiNSignalYieldFraction(zeroLeptonLowDeltaPhiNName+"_SignalYieldFraction", zeroLeptonLowDeltaPhiNName+"_SignalYieldFraction", *zeroLeptonLowDeltaPhiNSignalRawYield, totalSignalRawYield);
      
      RooProduct zeroLeptonSignalYield(zeroLeptonName+"_SignalYield", zeroLeptonName+"_SignalYield", RooArgSet(zeroLeptonSignalYieldFraction, *luminosity, *signalCrossSection) );
      RooProduct zeroLeptonLowDeltaPhiNSignalYield(zeroLeptonLowDeltaPhiNName+"_SignalYield", zeroLeptonLowDeltaPhiNName+"_SignalYield", RooArgSet(zeroLeptonLowDeltaPhiNSignalYieldFraction, *luminosity, *signalCrossSection) );
            
      ws.import(zeroLeptonSignalYield, RecycleConflictNodes());
      ws.import(zeroLeptonLowDeltaPhiNSignalYield, RecycleConflictNodes());

      //Make one lepton yields

      if(options.TopWJetsMethod == "ABCD")
	{
	  TString oneLeptonName("oneLepton_");
	  oneLeptonName+=binName;
	  
	  RooRealVar* oneLeptonSignalRawYield = ws.var(oneLeptonName+"_SignalRawYield");
	  RooRatio oneLeptonSignalYieldFraction(oneLeptonName+"_SignalYieldFraction", oneLeptonName+"_SignalYieldFraction", *oneLeptonSignalRawYield, totalSignalRawYield);
	  RooProduct oneLeptonSignalYield(oneLeptonName+"_SignalYield", oneLeptonName+"_SignalYield", RooArgSet(oneLeptonSignalYieldFraction, *luminosity, *signalCrossSection) );

	  ws.import(oneLeptonSignalYield, RecycleConflictNodes());
	  cout  << "zeroLeptonSignalYield=" << zeroLeptonSignalYield.getVal() << ", zeroLeptonLowDeltaPhiNSignalYield="  << zeroLeptonLowDeltaPhiNSignalYield.getVal() << ", oneLeptonSignalYield=" << oneLeptonSignalYield.getVal() << endl;
	}
      else if(options.TopWJetsMethod == "metReweighting")
	{
	  //inside
	  map<TString,int> thisInsideSignalMR = insideSignalMR[binName];
	  for ( map<TString,int>::iterator it=thisInsideSignalMR.begin() ; it != thisInsideSignalMR.end(); it++ ) 
	    { 
	      TString signalName = (*it).first;
	      RooRealVar* rooSignalRawYield = ws.var(signalName+"SignalRawYield");
	      RooRatio rooSignalYieldFraction(signalName+"_SignalYieldFraction", signalName+"_SignalYieldFraction", *rooSignalRawYield, totalSignalRawYield);
	      RooProduct rooSignalYield(signalName+"_SignalYield", signalName+"_SignalYield", RooArgSet(rooSignalYieldFraction, *luminosity, *signalCrossSection) );

	      ws.import(rooSignalYield, RecycleConflictNodes());
	    }
	  
	  //outside
	  map<TString,int> thisOutsideSignalMR = outsideSignalMR[binName];
	  for ( map<TString,int>::iterator it=thisOutsideSignalMR.begin() ; it != thisOutsideSignalMR.end(); it++ ) 
	    { 
	      TString signalName = (*it).first;
	      
	      //Continue if already in workspace 
	      if(ws.var(signalName+"_SignalYield") != NULL) continue;
	      
	      RooRealVar* rooSignalRawYield = ws.var(signalName+"SignalRawYield");
	      RooRatio rooSignalYieldFraction(signalName+"_SignalYieldFraction", signalName+"_SignalYieldFraction", *rooSignalRawYield, totalSignalRawYield);
	      RooProduct rooSignalYield(signalName+"_SignalYield", signalName+"_SignalYield", RooArgSet(rooSignalYieldFraction, *luminosity, *signalCrossSection) );

	      ws.import(rooSignalYield, RecycleConflictNodes());
	    }
	}
      else { assert(0); }
    }
  
}



void setupObservations(TString binName, TString binFileName, map<TString,abcdBinParameters>& bins, map<TString,channels>& observations)
{
  channels counts;
  abcdBinParameters abcd;
  
  ifstream binFile;
       
  cout << "getting the file: " << binFileName << " to fill bin " << binName << endl;
       
  binFile.open(binFileName.Data(),fstream::in);
  assert(binFile.is_open());
  
  string fileLine;
  
  TString index;
  TString value;
  
  while(!binFile.eof())
    {

      getline(binFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken nameAndNumber(thisLine," ");
      nameAndNumber.NextToken();
      index = nameAndNumber;
      nameAndNumber.NextToken();
      value = nameAndNumber;
      if(value.IsDigit() == false) 
	{
	  //hack to get rid of potential whitespace in strings
	  string str = (string)value;
	  for(unsigned int i=0; i<str.length(); i++) 
	    {
	    if(str[i] == '\t') str.erase(i,1);
	    if(str[i] == ' ') str.erase(i,1);
	    }
	  value = (TString)str;
	}
      cout << index << " : " << value << endl;
      if(index == "zeroLeptonCount"                                   ) counts.zeroLepton = value.Atof();		 
      else if(index == "zeroLeptonLowDeltaPhiNCount"                  ) counts.zeroLeptonLowDeltaPhiN = value.Atof();
      else if(index == "oneLeptonCount"                               ) counts.oneLepton = value.Atof();           
      else if(index == "oneElectronCount"                             ) counts.oneElectron = value.Atof();           
      else if(index == "oneMuonCount"                                 ) counts.oneMuon = value.Atof();		 
      else if(index == "diElectronCountName"                          ) counts.diElectronName = value;	     
      else if(index == "diElectronCount"                              ) counts.diElectron = value.Atof();	     
      else if(index == "diMuonCountName"	                      ) counts.diMuonName = value;	      
      else if(index == "diMuonCount"	                              ) counts.diMuon = value.Atof();	      
      else if(index == "zeroLeptonTriggerEfficiency"		      ) abcd.zeroLeptonTriggerEfficiency = value.Atof();			
      else if(index == "zeroLeptonTriggerEfficiencyError"	      ) abcd.zeroLeptonTriggerEfficiencyError = value.Atof();		
      else if(index == "zeroLeptonLowDeltaPhiNTriggerEfficiency"      ) abcd.zeroLeptonLowDeltaPhiNTriggerEfficiency = value.Atof();	
      else if(index == "zeroLeptonLowDeltaPhiNTriggerEfficiencyError" ) abcd.zeroLeptonLowDeltaPhiNTriggerEfficiencyError = value.Atof();	
      else if(index == "oneLeptonTriggerEfficiency"		      ) abcd.oneLeptonTriggerEfficiency = value.Atof();			
      else if(index == "oneLeptonTriggerEfficiencyError"    	      ) abcd.oneLeptonTriggerEfficiencyError = value.Atof();		
      else if(index == "oneElectronTriggerEfficiency"		      ) abcd.oneElectronTriggerEfficiency = value.Atof();			
      else if(index == "oneElectronTriggerEfficiencyError"    	      ) abcd.oneElectronTriggerEfficiencyError = value.Atof();		
      else if(index == "oneMuonTriggerEfficiency"	              ) abcd.oneMuonTriggerEfficiency = value.Atof();			
      else if(index == "oneMuonTriggerEfficiencyError"		      ) abcd.oneMuonTriggerEfficiencyError = value.Atof();			
      else if(index == "zeroLeptonLowDeltaPhiNMC"		      ) abcd.zeroLeptonLowDeltaPhiNMC = value.Atof();			
      else if(index == "topWJetsLowDeltaPhiNOverZeroLeptonRatioMC"    ) abcd.topWJetsLowDeltaPhiNOverZeroLeptonRatioMC = value.Atof();
      else if(index == "ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC"     ) abcd.ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC = value.Atof();
      else if(index == "ZtoNuNubTagScalingName"			      ) abcd.ZtoNuNubTagScalingName = value;				
      else if(index == "ZtoNuNubTagScaling"			      ) abcd.ZtoNuNubTagScaling = value.Atof();				
      else if(index == "ZtoNuNubTagScalingError"		      ) abcd.ZtoNuNubTagScalingError = value.Atof();			
      else if(index == "ZtomumuAcceptanceName"			      ) abcd.ZtomumuAcceptanceName = value;				
      else if(index == "ZtomumuAcceptance"			      ) abcd.ZtomumuAcceptance = value.Atof();				
      else if(index == "ZtomumuAcceptanceError"		              ) abcd.ZtomumuAcceptanceError = value.Atof();				
      else if(index == "ZtoeeAcceptanceName"			      ) abcd.ZtoeeAcceptanceName = value;				
      else if(index == "ZtoeeAcceptance"			      ) abcd.ZtoeeAcceptance = value.Atof();				
      else if(index == "ZtoeeAcceptanceError"		              ) abcd.ZtoeeAcceptanceError = value.Atof();				
      else if(index == "deltaPhiNRatioName"                           ) abcd.deltaPhiNRatioName = value;	     
      else if(index == "deltaPhiNRatio"	                              ) abcd.deltaPhiNRatio = value.Atof();	     
      else if(index == "deltaPhiNRatioError"                          ) abcd.deltaPhiNRatioError = value.Atof();  	
      else if(index == "lowDeltaPhiNScalingName"                      ) abcd.lowDeltaPhiNScalingName = value;	       
      else if(index == "qcdClosure"		        	      ) abcd.qcdClosure = value.Atof();					
      else if(index == "qcdClosureError"			      ) abcd.qcdClosureError = value.Atof();				
      else if(index == "topWJetsClosure"		              ) abcd.topWJetsClosure = value.Atof();				
      else if(index == "topWJetsClosureError"			      ) abcd.topWJetsClosureError = value.Atof();				
      else if(index == "ZtoeeSystematicName"		       	      ) abcd.ZtoeeSystematicName = value;					
      else if(index == "ZtoeeSystematic"			      ) abcd.ZtoeeSystematic = value.Atof();					
      else if(index == "ZtoeeSystematicError"                         ) abcd.ZtoeeSystematicError = value.Atof();                            
      else if(index == "ZtomumuSystematicName"		       	      ) abcd.ZtomumuSystematicName = value;					
      else if(index == "ZtomumuSystematic"			      ) abcd.ZtomumuSystematic = value.Atof();					
      else if(index == "ZtomumuSystematicError"                       ) abcd.ZtomumuSystematicError = value.Atof();                            
      else if(index != "") assert(0);
    }

  binFile.close();

  bins[binName] = abcd;
  observations[binName] = counts;
}

void setupUnderlyingModel(likelihoodOptions options, TString binFilesPath, map<TString,TString>& binFileNames,  
			  TString signalModelFilesPath, map<TString,TString>& binFileNamesInsideSignalMR,  map<TString,TString>& binFileNamesOutsideSignalMR,
			  vector<TString>& binNames, TString& modelFileName, TString& binFilesFileName, allBins& numbers, double& luminosityInput)
{
  ifstream setupFile;
       
  cout << "getting the file: " << modelFileName << " to setup likelihood " << endl;
  
  string fileLine;
  
  TString index;
  double value;
       
  setupFile.open(modelFileName.Data(),fstream::in);
  assert(setupFile.is_open());

  while(!setupFile.eof())
    {
      getline(setupFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken nameAndNumber(thisLine," ");
      nameAndNumber.NextToken();
      index = nameAndNumber;
      nameAndNumber.NextToken();
      value = nameAndNumber.Atof();
      cout << index << " : " << value << endl;
      if(index == "ZtollOverZtoNuNuRatio"  ) numbers.ZtollOverZtoNuNuRatio = value;
      else if(index == "ZtoeePurity"	        ) numbers.ZtoeePurity = value;	      	
      else if(index == "ZtoeePurityError"       ) numbers.ZtoeePurityError = value;     
      else if(index == "ZtomumuPurity"	        ) numbers.ZtomumuPurity = value;	      
      else if(index == "ZtomumuPurityError"     ) numbers.ZtomumuPurityError = value;   
      else if(index == "ZtoeeEfficiency"        ) numbers.ZtoeeEfficiency = value;	     
      else if(index == "ZtoeeEfficiencyError"   ) numbers.ZtoeeEfficiencyError = value; 
      else if(index == "ZtomumuEfficiency"      ) numbers.ZtomumuEfficiency = value;    
      else if(index == "ZtomumuEfficiencyError" ) numbers.ZtomumuEfficiencyError = value;	
      else if(index == "MCUncertainty" 	        ) numbers.MCUncertainty = value;	      	   
      else if(index == "Luminosity"   	        ) luminosityInput = value;
      else if(index != "") assert(0);
    }

  setupFile.close();

  ifstream binFilesFile;
  binFilesFile.open(binFilesFileName.Data(),fstream::in);

  TString fileName, fileNameInsideSignalMR, fileNameOutsideSignalMR;

  while(!binFilesFile.eof())
    {
      getline(binFilesFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken binString(thisLine," ");
      binString.NextToken();
      index = binString;
      if(index == "") continue;
      
      binString.NextToken();
      fileName = binString;
      binNames.push_back(index);
      binFileNames[index] = binFilesPath+fileName;
      cout << index << " : " << binFileNames[index] << endl;

      //MR
      if(options.TopWJetsMethod == "metReweighting") 
	{
	  binString.NextToken();
	  fileNameInsideSignalMR = binString;
	  binString.NextToken();
	  fileNameOutsideSignalMR = binString;
	  binFileNamesInsideSignalMR[index] = signalModelFilesPath+fileNameInsideSignalMR;
	  binFileNamesOutsideSignalMR[index] = signalModelFilesPath+fileNameOutsideSignalMR;
	  cout << index << " : " << binFileNamesInsideSignalMR[index] << ", " << binFileNamesOutsideSignalMR[index] << endl;
	}
   }
  binFilesFile.close();
}


void buildLikelihood( TString setupFileName, TString binFilesFileName, TString binFilesPath, TString signalModelFilesPath, TString workspaceName, TString outputFileName ) 
{
  
  double luminosityInput(1.);
  RooWorkspace ws (workspaceName) ;
  ws.autoImportClassCode(true);
  vector<TString> binNames;
  map<TString,TString> binFileNames;
  map<TString,TString> binFileNamesInsideSignalMR;
  map<TString,map<TString,int> > insideSignalMR;
  map<TString,TString> binFileNamesOutsideSignalMR;
  map<TString,map<TString,int> > outsideSignalMR;
  allBinNames names;
  allBins numbers;
  likelihoodOptions options;
  map<TString,abcdBinParameters> bins;
  map<TString,channels> observations;
  double oneLeptonTotal(0.);
  double zeroLeptonTopWJetsGuess(0.);

  
  options.skipTriggerEfficiency = true;
  options.qcdMethod = "singleScaleWithCorrections";//others: htDependent
  options.TopWJetsMethod = "ABCD"; //others: metReweighting
  
  //Read in setupFile and binFilesFile
  setupUnderlyingModel(options, binFilesPath, binFileNames, signalModelFilesPath, binFileNamesInsideSignalMR, binFileNamesOutsideSignalMR, binNames, setupFileName , binFilesFileName , numbers , luminosityInput);
  
  //Read in each bin file
  for(map<TString,TString>::iterator thisBin = binFileNames.begin(); thisBin != binFileNames.end() ; thisBin++)
    {
      setupObservations(thisBin->first , thisBin->second , bins, observations);
    }
  
  //Put "global" parameters into workspace
  makeUnderlyingLikelihood(options, ws , names, numbers, luminosityInput);

  //Read in MR signal
  if(options.TopWJetsMethod == "metReweighting")
    {
      for(vector<TString>::iterator thisBin = binNames.begin(); thisBin != binNames.end() ; thisBin++)
	{
	  setupSignalModelMR(*thisBin, binFileNamesInsideSignalMR[*thisBin], insideSignalMR, binFileNamesOutsideSignalMR[*thisBin], outsideSignalMR );
	}
    }
  
  //Read in signal and put into workspace
  setupSignalModel(options, ws , binNames , names, signalModelFilesPath, insideSignalMR, outsideSignalMR );
  
  //Do MET-reweighting method
  if(options.TopWJetsMethod == "metReweighting") buildMRLikelihood("","","","","","","",false);

  //Do nominal method, QCD, and ZtoNuNu
  for(vector<TString>::iterator thisBin = binNames.begin(); thisBin != binNames.end() ; thisBin++)
    {
      makeOneBin(options, ws , *thisBin , names , numbers, observations[*thisBin] , bins[*thisBin] , oneLeptonTotal, zeroLeptonTopWJetsGuess );
    }

  //Make guesses for parameters
  ws.var(names.singleLeptonScaling)->setVal(zeroLeptonTopWJetsGuess/oneLeptonTotal);
  for(vector<TString>::iterator thisBin = binNames.begin(); thisBin != binNames.end() ; thisBin++)
    {
      makeGuess(options, ws , *thisBin , observations[*thisBin] , bins[*thisBin] );
    }
  
  
  //Construct likelihood
  
  RooArgSet allpdfs = ws.allPdfs();
  
  cout << endl; cout << endl;
  cout << "allpdfs, size: " << allpdfs.getSize() << endl;
  allpdfs.Print("v");
  
  //RooProdPdf likelihood("likelihood","likelihood",allpdfs);
  RooProdPdfLogSum likelihood("likelihood","likelihood",allpdfs);
  
  ws.import(likelihood, RecycleConflictNodes());

  ws.defineSet("poi",names.signalCrossSection);  

  cout << endl; cout << endl;
  cout << "poi, size: " << (*ws.set("poi")).getSize() << endl;
  (*ws.set("poi")).Print("v");

  cout << endl; cout << endl;
  cout << "data, size: " << (*ws.set(names.observables)).getSize() << endl;
  (*ws.set(names.observables)).Print("v");

  RooDataSet data("data","data",*ws.set(names.observables));

  data.add(*ws.set(names.observables));

  ws.import(data);

  cout << "setting up models" << endl;

  ModelConfig sbModel("S+B_model",&ws);
  sbModel.SetPdf(*ws.pdf("likelihood"));
  sbModel.SetObservables(*ws.set(names.observables));
  sbModel.SetNuisanceParameters(*ws.set(names.nuisances));
  sbModel.SetParametersOfInterest(*ws.set("poi"));
  sbModel.SetProtoData(*ws.data("data"));

  ModelConfig bModel("B_model",&ws);
  bModel.SetPdf(*ws.pdf("likelihood"));
  bModel.SetObservables(*ws.set(names.observables));
  bModel.SetNuisanceParameters(*ws.set(names.nuisances));
  bModel.SetParametersOfInterest(*ws.set("poi"));
  bModel.SetProtoData(*ws.data("data"));
  ws.var(names.signalCrossSection)->setVal(0.0);
  bModel.SetSnapshot(*ws.set("poi"));

  ws.import (sbModel);
  ws.import (bModel);

  ws.Print() ;
  ws.writeToFile(outputFileName, true) ;

}


void likelihoodBuilder( TString setupFileName, TString binFilesFileName, TString binFilesPath, TString signalModelFilesPath, TString workspaceName, TString outputFileName ) {
  buildLikelihood( setupFileName, binFilesFileName, binFilesPath, signalModelFilesPath, workspaceName, outputFileName);
}
