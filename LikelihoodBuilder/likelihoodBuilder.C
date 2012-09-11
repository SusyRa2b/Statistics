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

#include "RooStats/ModelConfig.h"

#include "TMath.h"

using namespace RooFit ;
using namespace RooStats ;

struct  likelihoodOptions
{
  bool skipTriggerEfficiency;
  TString qcdMethod;
};

struct channels
{
  double zeroLepton;
  double zeroLeptonLowDeltaPhiN;
  double oneMuon;
  double oneElectron;
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
  double signal;
  double signalError;
} ;


bool makeOneBin(const likelihoodOptions options, RooWorkspace& ws , TString& binName , allBinNames& names , const allBins& numbers, channels& observed , abcdBinParameters& abcd , yields& signal , yields& signalError , double& oneLeptonTotal, double& zeroLeptonTopWJetsGuess)
{

  //Make sure that names are unique to this bin
  
  TString zeroLeptonName("zeroLepton_");
  zeroLeptonName+=binName;
  TString zeroLeptonLowDeltaPhiNName("zeroLeptonLowDeltaPhiN_");
  zeroLeptonLowDeltaPhiNName+=binName;
  TString oneLeptonName("oneLepton_");
  oneLeptonName+=binName;
  TString oneMuonName("oneMuon_");
  oneMuonName+=binName;
  TString oneElectronName("oneElectron_");
  oneElectronName+=binName;
  
  //Define unique counts and add them to the workspace

  RooRealVar zeroLeptonCount(zeroLeptonName+"_Count",zeroLeptonName+"_Count",observed.zeroLepton);
  RooRealVar zeroLeptonLowDeltaPhiNCount(zeroLeptonLowDeltaPhiNName+"_Count",zeroLeptonLowDeltaPhiNName+"_Count",observed.zeroLeptonLowDeltaPhiN);
  RooRealVar oneMuonCount(oneMuonName+"_Count",oneMuonName+"_Count",observed.oneMuon);
  RooRealVar oneElectronCount(oneElectronName+"_Count",oneElectronName+"_Count",observed.oneElectron);

  zeroLeptonCount.setConstant();
  zeroLeptonLowDeltaPhiNCount.setConstant();
  oneMuonCount.setConstant();
  oneElectronCount.setConstant();

  ws.import(zeroLeptonCount);
  ws.import(zeroLeptonLowDeltaPhiNCount);
  ws.import(oneMuonCount);
  ws.import(oneElectronCount);
  ws.extendSet(names.observables,zeroLeptonCount.GetName());
  ws.extendSet(names.observables,zeroLeptonLowDeltaPhiNCount.GetName());
  ws.extendSet(names.observables,oneMuonCount.GetName());
  ws.extendSet(names.observables,oneElectronCount.GetName());

  //Define QCD component
  /*
  //BEN FIXME -- old method commented out for now. should make a switch instead
  RooAbsArg* deltaPhiNRatio = ws.arg("deltaPhiNRatio_"+abcd.deltaPhiNRatioName+"_Ratio");
  if(deltaPhiNRatio == NULL) 
    {
      deltaPhiNRatio = 
	getBetaConstraint(ws,"deltaPhiNRatio_",abcd.deltaPhiNRatioName,
			  abcd.deltaPhiNRatio,abcd.deltaPhiNRatioError,
			  names.observables,names.nuisances);
    }
  double qcdGuess = observed.zeroLeptonLowDeltaPhiN - abcd.zeroLeptonLowDeltaPhiNMC;
  if(qcdGuess < 0) qcdGuess = 0;
  RooRealVar zeroLeptonLowDeltaPhiNQCDYield(zeroLeptonLowDeltaPhiNName+"_QCDYield",zeroLeptonLowDeltaPhiNName+"_QCDYield",qcdGuess,1e-5,10000);
  ws.import(zeroLeptonLowDeltaPhiNQCDYield);
  ws.extendSet(names.nuisances,zeroLeptonLowDeltaPhiNQCDYield.GetName());
  RooAbsArg* zeroLeptonQCDClosure = 
    getBetaPrimeConstraint(ws,"zeroLeptonQCDClosure_", binName,
			   abcd.qcdClosure,abcd.qcdClosureError,
			   names.observables,names.nuisances);
  RooProduct zeroLeptonQCDYield(zeroLeptonName+"_QCDYield",zeroLeptonName+"_QCDYield",RooArgSet(*deltaPhiNRatio,*zeroLeptonQCDClosure,zeroLeptonLowDeltaPhiNQCDYield));
  */

  //BEN FIXME -- for now, make this in option here.  consider moving option to input file "scalingName"
  TString scaling = "lowDeltaPhiNScaling_";
  if(options.qcdMethod == "singleScaleWithCorrections") scaling +=  "all";
  else if (options.qcdMethod == "htDependent") scaling +=  abcd.lowDeltaPhiNScalingName;
  else assert(0);
  RooRealVar* lowDeltaPhiNScaling = ws.var(scaling);
  if(lowDeltaPhiNScaling == NULL) {
    RooRealVar lowDeltaPhiNScaling_temp(scaling, scaling, 0.0, 1e3);
    lowDeltaPhiNScaling_temp.setVal(0.1);//BEN FIXME - automatic initial guess?
    ws.import(lowDeltaPhiNScaling_temp);
    ws.extendSet(names.nuisances,lowDeltaPhiNScaling_temp.GetName());
    lowDeltaPhiNScaling = ws.var(scaling);
  }
  double qcdGuess = observed.zeroLeptonLowDeltaPhiN - abcd.zeroLeptonLowDeltaPhiNMC;
  if(qcdGuess < 1e-5) qcdGuess = 1e-5;
  RooRealVar zeroLeptonLowDeltaPhiNQCDYield(zeroLeptonLowDeltaPhiNName+"_QCDYield",zeroLeptonLowDeltaPhiNName+"_QCDYield",qcdGuess,1e-5,1e5);
  ws.import(zeroLeptonLowDeltaPhiNQCDYield);
  ws.extendSet(names.nuisances,zeroLeptonLowDeltaPhiNQCDYield.GetName());
  RooAbsArg* zeroLeptonQCDClosure = 
    getBetaPrimeConstraint(ws,"zeroLeptonQCDClosure_", binName,
			   abcd.qcdClosure,abcd.qcdClosureError,
			   names.observables,names.nuisances);
  RooProduct zeroLeptonQCDYield(zeroLeptonName+"_QCDYield",zeroLeptonName+"_QCDYield",RooArgSet(*lowDeltaPhiNScaling,*zeroLeptonQCDClosure,zeroLeptonLowDeltaPhiNQCDYield));
  
  
  //Define top and W+jets component:
  
  RooRealVar* singleLeptonScaling = ws.var(names.singleLeptonScaling);
  RooRealVar  oneLeptonTopWJetsYield(oneLeptonName+"_TopWJetsYield",oneLeptonName+"_TopWJetsYield",0.5*(observed.oneMuon+observed.oneElectron),1e-5,1e5);
  ws.import(oneLeptonTopWJetsYield);
  ws.extendSet(names.nuisances,oneLeptonTopWJetsYield.GetName());
  oneLeptonTotal += oneLeptonTopWJetsYield.getVal();
  RooAbsArg*  zeroLeptonTopWJetsClosure = 
    getBetaPrimeConstraint(ws,"zeroLeptonTopWJetsClosure_", binName,
			   abcd.topWJetsClosure,abcd.topWJetsClosureError,
			   names.observables,names.nuisances);
  RooProduct  zeroLeptonTopWJetsYield(zeroLeptonName+"_TopWJetsYield",zeroLeptonName+"_TopWJetsYield",RooArgSet(*singleLeptonScaling,*zeroLeptonTopWJetsClosure,oneLeptonTopWJetsYield));
  

  //Define Z to invisible component:
  
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
      RooRealVar ZtoNuNu_temp("ZtoNuNu_"+observed.diMuonName,"ZtoNuNu_"+observed.diMuonName,(1./numbers.ZtollOverZtoNuNuRatio)*0.5*(observed.diMuon*numbers.ZtomumuPurity/(numbers.ZtomumuEfficiency*abcd.ZtomumuAcceptance) + observed.diElectron*numbers.ZtoeePurity/(numbers.ZtoeeEfficiency*abcd.ZtoeeAcceptance)),0.0,1e3);
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
  
  /*
  //BEN FIXME - systematics not implemented yet
  RooAbsArg*  ZtoeeSystematic = ws.arg("ZtoeeSystematic_"+abcd.ZtoeeSystematicName+"_Ratio");
  if(ZtoeeSystematic == NULL) 
    {
      ZtoeeSystematic = 
	getBetaPrimeConstraint(ws,"ZtoeeSystematic_",abcd.ZtoeeSystematicName,
			       abcd.ZtoeeSystematic,abcd.ZtoeeSystematicError,
			       names.observables,names.nuisances);
    }
  RooAbsArg*  ZtomumuSystematic = ws.arg("ZtomumuSystematic_"+abcd.ZtomumuSystematicName+"_Ratio");
  if(ZtomumuSystematic == NULL) 
    {
      ZtomumuSystematic = 
	getBetaPrimeConstraint(ws,"ZtomumuSystematic_",abcd.ZtomumuSystematicName,
			       abcd.ZtomumuSystematic,abcd.ZtomumuSystematicError,
			       names.observables,names.nuisances);
    }
  */

  RooProduct* diMuonYield = (RooProduct*)ws.arg("diMuon_"+observed.diMuonName+"_Yield");//Assumes acceptance is only binned in zero or more dimensions of the count.
  if(diMuonYield == NULL) 
    {
      RooProduct diMuonYield_temp("diMuon_"+observed.diMuonName+"_Yield","diMuon_"+observed.diMuonName+"_Yield",RooArgSet(*Ztoll,*ZtomumuAcceptance,*ZtomumuEfficiency,*ZtomumuInvPurity));
      ws.import(diMuonYield_temp, RecycleConflictNodes());
      diMuonYield = (RooProduct*)ws.arg("diMuon_"+observed.diMuonName+"_Yield");//Assumes acceptance is only binned in zero or more dimensions of the count.
    }
  
  RooProduct* diElectronYield = (RooProduct*)ws.arg("diElectron_"+observed.diElectronName+"_Yield");//Assumes acceptance is only binned in zero or more dimensions of the count.
  if(diElectronYield == NULL) 
    {
      RooProduct diElectronYield_temp("diElectron_"+observed.diElectronName+"_Yield","diElectron_"+observed.diElectronName+"_Yield",RooArgSet(*Ztoll,*ZtoeeAcceptance,*ZtoeeEfficiency,*ZtoeeInvPurity));
      ws.import(diElectronYield_temp, RecycleConflictNodes());
      diElectronYield = (RooProduct*)ws.arg("diElectron_"+observed.diElectronName+"_Yield");//Assumes acceptance is only binned in zero or more dimensions of the count.
    }
  
  //-----Define Z->ll Poisson constraints
  
  RooPoisson* diMuonConstraint = (RooPoisson*)ws.arg("diMuon_"+observed.diMuonName+"_Constraint");
  if(diMuonConstraint == NULL) 
    {
      RooPoisson diMuonConstraint_temp("diMuon_"+observed.diMuonName+"_Constraint","diMuon_"+observed.diMuonName+"_Constraint",*diMuonCount,*diMuonYield);
      ws.import(diMuonConstraint_temp, RecycleConflictNodes());
      diMuonConstraint = (RooPoisson*)ws.arg("diMuon_"+observed.diMuonName+"_Constraint");
    }
  
  RooPoisson* diElectronConstraint = (RooPoisson*)ws.arg("diElectron_"+observed.diElectronName+"_Constraint");
  if(diElectronConstraint == NULL) 
    {
      RooPoisson diElectronConstraint_temp("diElectron_"+observed.diElectronName+"_Constraint","diElectron_"+observed.diElectronName+"_Constraint",*diElectronCount,*diElectronYield);
      ws.import(diElectronConstraint_temp, RecycleConflictNodes());
      diElectronConstraint = (RooPoisson*)ws.arg("diElectron_"+observed.diElectronName+"_Constraint");
    }
  
  //-----Define Z->nunu yield
  
  RooAbsArg* zeroLeptonZtoNuNubTagScaling = getBetaConstraint(ws,"zeroLeptonZtoNuNubTagScaling_",abcd.ZtoNuNubTagScalingName,
							      abcd.ZtoNuNubTagScaling,abcd.ZtoNuNubTagScalingError,
							      names.observables,names.nuisances);
  
  RooProduct zeroLeptonZtoNuNuYield(zeroLeptonName+"_ZtoNuNuYield",zeroLeptonName+"_ZtoNuNuYield",RooArgSet(*zeroLeptonZtoNuNubTagScaling,*ZtoNuNu));
  
  //  Non-QCD/Signal in Low  Delta Phi_N  region

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


  // Setup signal yields

  RooRealVar* signalCrossSection = ws.var(names.signalCrossSection);
  
  RooAbsArg*  zeroLeptonSignalYieldFraction = 
    getCorrelatedBetaPrimeConstraint(ws,"zeroLeptonSignalYieldFraction_",binName,
				     signal.zeroLepton,signalError.zeroLepton,
				     names.observables,names.nuisances,
				     names.signalUncertainty);

  RooAbsArg*  zeroLeptonLowDeltaPhiNSignalYieldFraction = 
    getCorrelatedBetaPrimeConstraint(ws,"zeroLeptonLowDeltaPhiNSignalYieldFraction_",binName,
				     signal.zeroLeptonLowDeltaPhiN,signalError.zeroLeptonLowDeltaPhiN,
				     names.observables,names.nuisances,
				     names.signalUncertainty);
  
  RooAbsArg*  oneLeptonSignalYieldFraction = 
    getCorrelatedBetaPrimeConstraint(ws,"oneLeptonSignalYieldFraction_",binName,
				     signal.oneLepton,signalError.oneLepton,
				     names.observables,names.nuisances,
				     names.signalUncertainty);
  
  RooProduct zeroLeptonSignalYield(zeroLeptonName+"_SignalYield",zeroLeptonName+"_SignalYield",RooArgSet(*signalCrossSection,*zeroLeptonSignalYieldFraction));
  RooProduct zeroLeptonLowDeltaPhiNSignalYield(zeroLeptonLowDeltaPhiNName+"_SignalYield",zeroLeptonLowDeltaPhiNName+"_SignalYield",RooArgSet(*signalCrossSection,*zeroLeptonLowDeltaPhiNSignalYieldFraction));
  RooProduct oneLeptonSignalYield(oneLeptonName+"_SignalYield",oneLeptonName+"_SignalYield",RooArgSet(*signalCrossSection,*oneLeptonSignalYieldFraction));
  
  // Setup yields in all bins
  
  RooAddition zeroLeptonYieldSum(zeroLeptonName+"_YieldSum",zeroLeptonName+"_YieldSum",RooArgSet(zeroLeptonSignalYield,zeroLeptonZtoNuNuYield,zeroLeptonTopWJetsYield,zeroLeptonQCDYield));
  double topGuess = observed.zeroLepton - zeroLeptonZtoNuNuYield.getVal() - zeroLeptonQCDYield.getVal();
  if(topGuess > 0 ) zeroLeptonTopWJetsGuess += topGuess;
  //RooAddition zeroLeptonLowDeltaPhiNYieldSum(zeroLeptonLowDeltaPhiNName+"_YieldSum",zeroLeptonLowDeltaPhiNName+"_YieldSum",RooArgSet(zeroLeptonLowDeltaPhiNSignalYield,zeroLeptonLowDeltaPhiNQCDYield,zeroLeptonLowDeltaPhiNMCYield));//Old non-QCD subtraction
  RooAddition zeroLeptonLowDeltaPhiNYieldSum(zeroLeptonLowDeltaPhiNName+"_YieldSum",zeroLeptonLowDeltaPhiNName+"_YieldSum",RooArgSet(zeroLeptonLowDeltaPhiNSignalYield,zeroLeptonLowDeltaPhiNQCDYield,zeroLeptonLowDeltaPhiNNonQCDYield));
  RooAddition oneLeptonYieldSum(oneLeptonName+"_YieldSum",oneLeptonName+"_YieldSum",RooArgSet(oneLeptonSignalYield,oneLeptonTopWJetsYield));
  
  // Setup trigger efficiencies
  if (options.skipTriggerEfficiency == true)
    {
      
      // Total Yields in bins
      
      RooProduct zeroLeptonYield(zeroLeptonName+"_Yield",zeroLeptonName+"_Yield",RooArgSet(zeroLeptonYieldSum));
      RooProduct zeroLeptonLowDeltaPhiNYield(zeroLeptonLowDeltaPhiNName+"_Yield",zeroLeptonLowDeltaPhiNName+"_Yield",RooArgSet(zeroLeptonLowDeltaPhiNYieldSum));
      RooProduct oneMuonYield(oneMuonName+"_Yield",oneMuonName+"_Yield",RooArgSet(oneLeptonYieldSum));
      RooProduct oneElectronYield(oneElectronName+"_Yield",oneElectronName+"_Yield",RooArgSet(oneLeptonYieldSum));
      
      // Define poisson constraints
      
      RooPoisson zeroLeptonConstraint(zeroLeptonName+"_Constraint",zeroLeptonName+"_Constraint",zeroLeptonCount,zeroLeptonYield);
      RooPoisson zeroLeptonLowDeltaPhiNConstraint(zeroLeptonLowDeltaPhiNName+"_Constraint",zeroLeptonLowDeltaPhiNName+"_Constraint",zeroLeptonLowDeltaPhiNCount,zeroLeptonLowDeltaPhiNYield);
      RooPoisson oneMuonConstraint(oneMuonName+"_Constraint",oneMuonName+"_Constraint",oneMuonCount,oneMuonYield);
      RooPoisson oneElectronConstraint(oneElectronName+"_Constraint",oneElectronName+"_Constraint",oneElectronCount,oneElectronYield);
      
      // Cleanup
      
      ws.import(zeroLeptonConstraint, RecycleConflictNodes());
      ws.import(zeroLeptonLowDeltaPhiNConstraint, RecycleConflictNodes());
      ws.import(oneMuonConstraint, RecycleConflictNodes());
      ws.import(oneElectronConstraint, RecycleConflictNodes());

    }
  else
    {
      RooAbsArg* zeroLeptonTriggerEfficiency = 
	getBetaConstraint(ws,"zeroLeptonTriggerEfficiency_",binName,
			  abcd.zeroLeptonTriggerEfficiency,abcd.zeroLeptonTriggerEfficiencyError,
			  names.observables,names.nuisances);
      
      RooAbsArg* zeroLeptonLowDeltaPhiNTriggerEfficiency = 
	getBetaConstraint(ws,"zeroLeptonLowDeltaPhiNTriggerEfficiency_",binName,
			  abcd.zeroLeptonLowDeltaPhiNTriggerEfficiency,abcd.zeroLeptonLowDeltaPhiNTriggerEfficiencyError,
			  names.observables,names.nuisances);
      
      RooAbsArg* oneMuonTriggerEfficiency = 
	getBetaConstraint(ws,"oneMuonTriggerEfficiency_",binName,
			  abcd.oneMuonTriggerEfficiency,abcd.oneMuonTriggerEfficiencyError,
			  names.observables,names.nuisances);
      
      RooAbsArg* oneElectronTriggerEfficiency = 
	getBetaConstraint(ws,"oneElectronTriggerEfficiency_",binName,
			  abcd.oneElectronTriggerEfficiency,abcd.oneElectronTriggerEfficiencyError,
			  names.observables,names.nuisances);
      
      // Total Yields in bins
      
      RooProduct zeroLeptonYield(zeroLeptonName+"_Yield",zeroLeptonName+"_Yield",RooArgSet(*zeroLeptonTriggerEfficiency,zeroLeptonYieldSum));
      RooProduct zeroLeptonLowDeltaPhiNYield(zeroLeptonLowDeltaPhiNName+"_Yield",zeroLeptonLowDeltaPhiNName+"_Yield",RooArgSet(*zeroLeptonLowDeltaPhiNTriggerEfficiency,zeroLeptonLowDeltaPhiNYieldSum));
      RooProduct oneMuonYield(oneMuonName+"_Yield",oneMuonName+"_Yield",RooArgSet(*oneMuonTriggerEfficiency,oneLeptonYieldSum));
      RooProduct oneElectronYield(oneElectronName+"_Yield",oneElectronName+"_Yield",RooArgSet(*oneElectronTriggerEfficiency,oneLeptonYieldSum));
      
      // Define poisson constraints
      
      RooPoisson zeroLeptonConstraint(zeroLeptonName+"_Constraint",zeroLeptonName+"_Constraint",zeroLeptonYield,zeroLeptonCount);
      RooPoisson zeroLeptonLowDeltaPhiNConstraint(zeroLeptonLowDeltaPhiNName+"_Constraint",zeroLeptonLowDeltaPhiNName+"_Constraint",zeroLeptonLowDeltaPhiNYield,zeroLeptonLowDeltaPhiNCount);
      RooPoisson oneMuonConstraint(oneMuonName+"_Constraint",oneMuonName+"_Constraint",oneMuonYield,oneMuonCount);
      RooPoisson oneElectronConstraint(oneElectronName+"_Constraint",oneElectronName+"_Constraint",oneElectronYield,oneElectronCount);
      
      // Cleanup
      
      ws.import(zeroLeptonConstraint, RecycleConflictNodes());
      ws.import(zeroLeptonLowDeltaPhiNConstraint, RecycleConflictNodes());
      ws.import(oneMuonConstraint, RecycleConflictNodes());
      ws.import(oneElectronConstraint, RecycleConflictNodes());
    }

  return true;

}

void setupUnderlyingLikelihood(const likelihoodOptions options, RooWorkspace& ws ,allBinNames& names, allBins& numbers)
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

void setupSignalModel(vector<TString> binNames, TString signalModelFileName, int signalModelLineNumber, map<TString,yields>& signal,map<TString,yields>& signalError, const double& luminosity)
{
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

  signalFile>>m0;
  cout << "m0 : " << m0 << endl;
  signalFile>>m12;
  cout << "m12 : " << m12 << endl;
  signalFile>>susyGenerated;
  cout << "SUSY Generated : " << susyGenerated << endl;

  // Luminosity*sigma = number of events 

  double crossSectionScaling = luminosity/susyGenerated;

  for(vector<TString>::iterator thisBin = binNames.begin() ; thisBin != binNames.end() ; thisBin++)
    {
      yields thisSignal,thisSignalError;
      double valueHolder;
      signalFile>>valueHolder ; thisSignal.zeroLepton = valueHolder*crossSectionScaling ;
      signalFile>>valueHolder ; thisSignalError.zeroLepton = valueHolder*crossSectionScaling ;
      signalFile>>valueHolder ; thisSignal.zeroLeptonLowDeltaPhiN = valueHolder*crossSectionScaling ;
      signalFile>>valueHolder ; thisSignalError.zeroLeptonLowDeltaPhiN = valueHolder*crossSectionScaling ;
      signalFile>>valueHolder ; thisSignal.oneLepton = valueHolder*crossSectionScaling ;
      signalFile>>valueHolder ; thisSignalError.oneLepton = valueHolder*crossSectionScaling ;
      signal[*thisBin]=thisSignal;
      signalError[*thisBin]=thisSignalError;
      cout << "For selection " << *thisBin << " signal is " << endl;
      cout << "zero lepton bin               : " << thisSignal.zeroLepton             << endl ;
      cout << "zero lepton low delta phi bin : " << thisSignal.zeroLeptonLowDeltaPhiN << endl ;
      cout << "one lepton bin                : " << thisSignal.oneLepton              << endl ;
      cout << "For selection " << *thisBin << " signal error is " << endl;
      cout << "zero lepton bin               : " << thisSignalError.zeroLepton             << endl ;
      cout << "zero lepton low delta phi bin : " << thisSignalError.zeroLeptonLowDeltaPhiN << endl ;
      cout << "one lepton bin                : " << thisSignalError.oneLepton              << endl ;
    }

  signalFile.close();

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
      else if(index == "oneMuonCount"                                 ) counts.oneMuon = value.Atof();		 
      else if(index == "oneElectronCount"                             ) counts.oneElectron = value.Atof();           
      else if(index == "diElectronCountName"                          ) counts.diElectronName = value;	     
      else if(index == "diElectronCount"                              ) counts.diElectron = value.Atof();	     
      else if(index == "diMuonCountName"	                      ) counts.diMuonName = value;	      
      else if(index == "diMuonCount"	                              ) counts.diMuon = value.Atof();	      
      else if(index == "zeroLeptonTriggerEfficiency"		      ) abcd.zeroLeptonTriggerEfficiency = value.Atof();			
      else if(index == "zeroLeptonTriggerEfficiencyError"	      ) abcd.zeroLeptonTriggerEfficiencyError = value.Atof();		
      else if(index == "zeroLeptonLowDeltaPhiNTriggerEfficiency"      ) abcd.zeroLeptonLowDeltaPhiNTriggerEfficiency = value.Atof();	
      else if(index == "zeroLeptonLowDeltaPhiNTriggerEfficiencyError" ) abcd.zeroLeptonLowDeltaPhiNTriggerEfficiencyError = value.Atof();	
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

void setupUnderlyingModel(TString binFilesPath, map<TString,TString>& binFileNames, vector<TString>& binNames, TString& modelFileName , TString& binFilesFileName , allBins& numbers, double& luminosity)
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
      else if(index == "Luminosity"   	        ) luminosity = value;	      	   
      else if(index != "") assert(0);
    }

  setupFile.close();

  setupFile.open(binFilesFileName.Data(),fstream::in);

  TString fileName;

  while(!setupFile.eof())
    {
      getline(setupFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken nameAndNumber(thisLine," ");
      nameAndNumber.NextToken();
      index = nameAndNumber;
      nameAndNumber.NextToken();
      fileName = nameAndNumber;
      cout << index << " : " << fileName << endl;
      if(index == "") continue;
      binNames.push_back(index);
      binFileNames[index] = binFilesPath+fileName;
    }
  setupFile.close();
}


void buildLikelihood( TString binFilesPath, TString setupFileName, TString binFilesFileName, TString signalModelFileName, int signalModelFileLine, TString workspaceName, TString outputFileName ) 
{
  
  double luminosity(1.);
  RooWorkspace ws (workspaceName) ;
  ws.autoImportClassCode(true);
  vector<TString> binNames;
  map<TString,TString> binFileNames;
  allBinNames names;
  allBins numbers;
  likelihoodOptions options;
  map<TString,abcdBinParameters> bins;
  map<TString,channels> observations;
  map<TString,yields> signal;
  map<TString,yields> signalError;
  double oneLeptonTotal(0.);
  double zeroLeptonTopWJetsGuess(0.);

  options.skipTriggerEfficiency = true;
  options.qcdMethod = "singleScaleWithCorrections";//others: htDependent

  setupUnderlyingModel(binFilesPath, binFileNames, binNames, setupFileName , binFilesFileName , numbers , luminosity);
  for(map<TString,TString>::iterator thisBin = binFileNames.begin(); thisBin != binFileNames.end() ; thisBin++)
    {
      setupObservations(thisBin->first , thisBin->second , bins, observations);
    }

  setupSignalModel(binNames , signalModelFileName , signalModelFileLine , signal , signalError , luminosity);

  setupUnderlyingLikelihood(options, ws , names, numbers);

  for(vector<TString>::iterator thisBin = binNames.begin(); thisBin != binNames.end() ; thisBin++)
    {
      makeOneBin(options, ws , *thisBin , names , numbers, observations[*thisBin] , bins[*thisBin] , signal[*thisBin] , signalError[*thisBin] , oneLeptonTotal, zeroLeptonTopWJetsGuess );
    }

  ws.var(names.singleLeptonScaling)->setVal(zeroLeptonTopWJetsGuess/oneLeptonTotal);

  RooArgSet allpdfs = ws.allPdfs();
  
  cout << endl; cout << endl;
  cout << "allpdfs, size: " << allpdfs.getSize() << endl;
  allpdfs.Print("v");
  
  RooProdPdf likelihood("likelihood","likelihood",allpdfs);
  
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
  ws.writeToFile(outputFileName) ;

}


void likelihoodBuilder( TString binFilesPath, TString setupFileName, TString binFilesFileName, TString signalModelFileName, int signalModelFileLine, TString workspaceName, TString outputFileName ) {
  buildLikelihood(binFilesPath, setupFileName, binFilesFileName, signalModelFileName, signalModelFileLine, workspaceName, outputFileName);
}
