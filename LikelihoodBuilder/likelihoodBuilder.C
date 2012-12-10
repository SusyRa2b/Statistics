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
#include "RooGaussian.h"

#include "RooCorrelatedBetaGeneratorHelper.h"
#include "RooCorrelatedBetaPrimeGeneratorHelper.h"
#include "betaHelperFunctions.h"
#include "RooBetaInverseCDF.h"
#include "RooBetaPrimeInverseCDF.h"
#include "rooFitBetaHelperFunctions.h"
#include "RooNormalFromFlatPdf.h"

#include "rooFitLogNormalHelperFunctions.h"
//For truncated gaussians, uncomment this line and two lines in setup.C
//You must also check out RA2b/Statistics/3Dcode to get a couple of files
#include "rooFitGaussianHelperFunctions.h"


#include "RooProdPdfLogSum.h"
#include "RooPoissonLogEval.h"

//#include "metReweightingBuilder.h"
#include "metReweightingBuilderSIMPLETAU.h"

#include "RooStats/ModelConfig.h"

#include "TMath.h"

using namespace RooFit ;
using namespace RooStats ;

struct  likelihoodOptions
{
  bool skipTriggerEfficiency;
  TString qcdMethod;
  TString TopWJetsMethod;
  TString nuisanceOption;
};


struct mcCount
{
  double value;
  double error;
};


struct channels
{
  double zeroLepton;
  double zeroLeptonLowDeltaPhiN;
  double oneLepton;
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
  TString zeroLeptonTriggerEfficiencyName;
  double zeroLeptonTriggerEfficiency;
  double zeroLeptonTriggerEfficiencyError;
  TString oneLeptonTriggerEfficiencyName;
  double oneLeptonTriggerEfficiency;
  double oneLeptonTriggerEfficiencyError;
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
  TString lowDeltaPhiNScalingName;
  TString lowDeltaPhiNMETScaleFactorName;
  TString lowDeltaPhiNBTagScaleFactorName;
  double qcdClosure;
  double qcdClosureError;
  double topWJetsClosure;
  double topWJetsClosureError;
  double zeroLeptonDibosonMC;
  double zeroLeptonLowDeltaPhiNDibosonMC;
  double oneLeptonDibosonMC;
} ;

//Used to store strings to avoid typos from lots of hardcoding
struct allBinNames
{
  TString ZtoeeInvPurity;
  TString ZtomumuInvPurity;
  TString ZtoeeEfficiency;
  TString ZtomumuEfficiency;
  TString ZtoeeSystematic;
  TString ZtomumuSystematic;
  TString ZtollOverZtoNuNuRatio;
  TString singleLeptonScaling;
  TString MCUncertainty;
  TString dibosonMCUncertainty;
  TString signalCrossSection;
  TString signalGlobalUncertainty;
  TString observables;
  TString nuisances;
  TString globalObservables;
} ;

//Values for things used in all bins
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
  double ZtoeeSystematic;
  double ZtoeeSystematicError;
  double ZtomumuSystematic;
  double ZtomumuSystematicError;
  double MC;
  double MCUncertainty;
  double dibosonMC;
  double dibosonMCUncertainty;
  double Luminosity;
  double LuminosityError;
  double metCleaningError;
  double SFqcd_met3;
  double SFqcd_met3_err;
  double SFqcd_met4;
  double SFqcd_met4_err;
  double SFqcd_nb3;
  double SFqcd_nb3_err;
} ;


bool skipBin(TString binName) 
{

  return false; //no skipping right now

  if(binName=="bin37" || binName=="bin38" || binName=="bin39") return true;//FIXME  --  hardcoded!
  return false;
}


void makeTriggerEfficiencies(RooWorkspace& ws, allBinNames& names, abcdBinParameters& abcd )
{
  
  RooAbsArg* zeroLeptonTriggerEfficiency =
    getGaussianConstraint(ws, "zeroLeptonTriggerEfficiency_", abcd.zeroLeptonTriggerEfficiencyName, //FIXME change to beta after synch
			  abcd.zeroLeptonTriggerEfficiency, abcd.zeroLeptonTriggerEfficiencyError,
			  names.observables,names.nuisances, names.globalObservables);
  ws.import(*zeroLeptonTriggerEfficiency, RecycleConflictNodes());

  RooAbsArg* oneLeptonTriggerEfficiency =
    getGaussianConstraint(ws, "oneLeptonTriggerEfficiency_", abcd.oneLeptonTriggerEfficiencyName, //FIXME change to beta after synch
			  abcd.oneLeptonTriggerEfficiency, abcd.oneLeptonTriggerEfficiencyError,
			  names.observables,names.nuisances, names.globalObservables);
  ws.import(*oneLeptonTriggerEfficiency, RecycleConflictNodes());

}

bool makeOneBin(const likelihoodOptions options, RooWorkspace& ws , TString& binName , allBinNames& names , const allBins& numbers, channels& observed , abcdBinParameters& abcd ,  double& oneLeptonTotal, double& zeroLeptonTopWJetsGuess)
{
  
  //Get trigger efficiencies from workspace
  RooAbsArg* zeroLeptonTriggerEfficiency = ws.arg("zeroLeptonTriggerEfficiency_"+abcd.zeroLeptonTriggerEfficiencyName);
  RooAbsArg* oneLeptonTriggerEfficiency = ws.arg("oneLeptonTriggerEfficiency_"+abcd.oneLeptonTriggerEfficiencyName);
  assert(zeroLeptonTriggerEfficiency!=NULL); assert(oneLeptonTriggerEfficiency!=NULL);
  
  
  //zero lepton name and count
  TString zeroLeptonName("zeroLepton_");
  zeroLeptonName+=binName;

  RooRealVar zeroLeptonCount(zeroLeptonName+"_Count",zeroLeptonName+"_Count",observed.zeroLepton);
  zeroLeptonCount.setConstant();
  
  ws.import(zeroLeptonCount);
  ws.extendSet(names.observables,zeroLeptonCount.GetName());

  

  //////////////////////////
  // QCD
  //////////////////////////
  
  //LDP name and count
  TString zeroLeptonLowDeltaPhiNName("zeroLeptonLowDeltaPhiN_");
  zeroLeptonLowDeltaPhiNName+=binName;
  
  RooRealVar zeroLeptonLowDeltaPhiNCount(zeroLeptonLowDeltaPhiNName+"_Count",zeroLeptonLowDeltaPhiNName+"_Count",observed.zeroLeptonLowDeltaPhiN);
  zeroLeptonLowDeltaPhiNCount.setConstant();
  
  ws.import(zeroLeptonLowDeltaPhiNCount);
  ws.extendSet(names.observables,zeroLeptonLowDeltaPhiNCount.GetName());
  

  //zeroLepton over lowDeltaPhiN ratio
  TString lowDeltaPhiNScalingName = "lowDeltaPhiNScaling_";   //for now, make this in option here. consider moving option to input file "scalingName"
  if(options.qcdMethod == "singleScaleWithCorrections") lowDeltaPhiNScalingName +=  "all";
  else if (options.qcdMethod == "model4") lowDeltaPhiNScalingName +=  abcd.lowDeltaPhiNScalingName;
  else assert(0);
  
  RooRealVar* lowDeltaPhiNScaling = ws.var(lowDeltaPhiNScalingName);
  if(lowDeltaPhiNScaling == NULL) 
    {
      RooRealVar lowDeltaPhiNScaling_temp(lowDeltaPhiNScalingName, lowDeltaPhiNScalingName, 0.2, 0.0, 1e1);
      ws.import(lowDeltaPhiNScaling_temp);
      ws.extendSet(names.nuisances,lowDeltaPhiNScaling_temp.GetName());
      lowDeltaPhiNScaling = ws.var(lowDeltaPhiNScalingName);
    }

 
  //Correction scale factors (rather hardcoded right now)

  //MET
  TString lowDeltaPhiNMETScaleFactorName = "lowDeltaPhiNMETScaleFactor_";
  lowDeltaPhiNMETScaleFactorName += abcd.lowDeltaPhiNMETScaleFactorName;
  RooRealVar* lowDeltaPhiNMETScaleFactor = ws.var(lowDeltaPhiNMETScaleFactorName);
  if(lowDeltaPhiNMETScaleFactor == NULL)
    {
      RooRealVar lowDeltaPhiNMETScaleFactor_temp(lowDeltaPhiNMETScaleFactorName, lowDeltaPhiNMETScaleFactorName, 1.0, 0.0, 3.0);
      lowDeltaPhiNMETScaleFactor_temp.setConstant();
      if( (options.qcdMethod == "model4") && (abcd.lowDeltaPhiNMETScaleFactorName != "M1") )
	{
	  lowDeltaPhiNMETScaleFactor_temp.setConstant(kFALSE);
	  if(abcd.lowDeltaPhiNMETScaleFactorName == "M3")
	    {
	      getGaussianConstraint(ws, lowDeltaPhiNMETScaleFactorName, "",
				    numbers.SFqcd_met3, numbers.SFqcd_met3_err,
				    names.observables,names.nuisances, names.globalObservables);		
	    }
	  else if(abcd.lowDeltaPhiNMETScaleFactorName == "M4")
	    {
	      getGaussianConstraint(ws, lowDeltaPhiNMETScaleFactorName, "",
				    numbers.SFqcd_met4, numbers.SFqcd_met4_err,
				    names.observables,names.nuisances, names.globalObservables);		
	    }
	} 
      ws.import(lowDeltaPhiNMETScaleFactor_temp, RecycleConflictNodes());
      ws.extendSet(names.nuisances, lowDeltaPhiNMETScaleFactor_temp.GetName());
      lowDeltaPhiNMETScaleFactor = ws.var(lowDeltaPhiNMETScaleFactorName);
    }
  
  
  //BTag
  TString lowDeltaPhiNBTagScaleFactorName = "lowDeltaPhiNBTagScaleFactor_";
  lowDeltaPhiNBTagScaleFactorName += abcd.lowDeltaPhiNBTagScaleFactorName;
  RooRealVar* lowDeltaPhiNBTagScaleFactor = ws.var(lowDeltaPhiNBTagScaleFactorName);
  if(lowDeltaPhiNBTagScaleFactor == NULL)
    {
      RooRealVar lowDeltaPhiNBTagScaleFactor_temp(lowDeltaPhiNBTagScaleFactorName, lowDeltaPhiNBTagScaleFactorName, 1.0, 0.0, 3.0);
      lowDeltaPhiNBTagScaleFactor_temp.setConstant();
      if( (options.qcdMethod == "model4") && (abcd.lowDeltaPhiNBTagScaleFactorName != "1b") )
	{
	  lowDeltaPhiNBTagScaleFactor_temp.setConstant(kFALSE);
	  if(abcd.lowDeltaPhiNBTagScaleFactorName == "3b")
	    {
	      getGaussianConstraint(ws, lowDeltaPhiNBTagScaleFactorName, "",
				    numbers.SFqcd_nb3, numbers.SFqcd_nb3_err,
				    names.observables,names.nuisances, names.globalObservables);		
	    }
	} 
      ws.import(lowDeltaPhiNBTagScaleFactor_temp, RecycleConflictNodes());
      ws.extendSet(names.nuisances, lowDeltaPhiNBTagScaleFactor_temp.GetName());
      lowDeltaPhiNBTagScaleFactor = ws.var(lowDeltaPhiNBTagScaleFactorName);
    }
  
  
  double qcdGuess = observed.zeroLeptonLowDeltaPhiN - abcd.zeroLeptonLowDeltaPhiNMC;
  if(qcdGuess < 1e-5) qcdGuess = 1e-5;
  
  RooRealVar zeroLeptonLowDeltaPhiNQCDYield(zeroLeptonLowDeltaPhiNName+"_QCDYield",zeroLeptonLowDeltaPhiNName+"_QCDYield",qcdGuess,0,1e5);
  ws.import(zeroLeptonLowDeltaPhiNQCDYield);
  ws.extendSet(names.nuisances,zeroLeptonLowDeltaPhiNQCDYield.GetName());
  
  RooRealVar* zeroLeptonQCDClosure = (RooRealVar*)
    //getBetaPrimeConstraint(ws,"zeroLeptonQCDClosure_", binName,
    getGaussianConstraint(ws,"zeroLeptonQCDClosure_", binName,
			   abcd.qcdClosure,abcd.qcdClosureError,
			   names.observables,names.nuisances, names.globalObservables);
  
  RooProduct zeroLeptonQCDYield(zeroLeptonName+"_QCDYield",zeroLeptonName+"_QCDYield",RooArgSet(*lowDeltaPhiNScaling,*zeroLeptonQCDClosure,*lowDeltaPhiNMETScaleFactor,*lowDeltaPhiNBTagScaleFactor,zeroLeptonLowDeltaPhiNQCDYield));
  
  

  //////////////////////////////
  // Top and W+jets 
  //////////////////////////////
  TString oneLeptonName("oneLepton_");
  oneLeptonName+=binName;
  
  if(options.TopWJetsMethod == "ABCD")
    {
      RooRealVar* singleLeptonScaling = ws.var(names.singleLeptonScaling);
      RooRealVar  oneLeptonTopWJetsYield(oneLeptonName+"_TopWJetsYield",oneLeptonName+"_TopWJetsYield",observed.oneLepton+1.0,0,1e5);
      ws.import(oneLeptonTopWJetsYield);
      ws.extendSet(names.nuisances,oneLeptonTopWJetsYield.GetName());
      oneLeptonTotal += oneLeptonTopWJetsYield.getVal();
      RooRealVar*  zeroLeptonTopWJetsClosure = (RooRealVar*)
	//getBetaPrimeConstraint(ws,"zeroLeptonTopWJetsClosure_", binName,
	getGaussianConstraint(ws,"zeroLeptonTopWJetsClosure_", binName,
			      abcd.topWJetsClosure,abcd.topWJetsClosureError,
			      names.observables,names.nuisances, names.globalObservables);
      RooProduct  zeroLeptonTopWJetsYield(zeroLeptonName+"_TopWJetsYield",zeroLeptonName+"_TopWJetsYield",RooArgSet(*singleLeptonScaling,*zeroLeptonTopWJetsClosure,oneLeptonTopWJetsYield));
      ws.import(zeroLeptonTopWJetsYield, RecycleConflictNodes() );
    }
  
  RooAbsReal* zeroLeptonTopWJetsYield = ws.function(zeroLeptonName+"_TopWJetsYield");
  assert(zeroLeptonTopWJetsYield != NULL);

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
  
  RooRealVar* ZtoNuNu = ws.var("ZtoNuNu_VeryLooseBtagYield_"+observed.diMuonName);//use MuonName, which should be the same as ElectronName
  if(ZtoNuNu == NULL) 
    {
      RooRealVar ZtoNuNu_temp("ZtoNuNu_VeryLooseBtagYield_"+observed.diMuonName,"ZtoNuNu_VeryLooseBtagYield_"+observed.diMuonName,(1./numbers.ZtollOverZtoNuNuRatio)*0.5*(observed.diMuon*numbers.ZtomumuPurity/(numbers.ZtomumuEfficiency*abcd.ZtomumuAcceptance) + observed.diElectron*numbers.ZtoeePurity/(numbers.ZtoeeEfficiency*abcd.ZtoeeAcceptance))+1.0,0.0,1e5);
      ws.import(ZtoNuNu_temp);
      ZtoNuNu = ws.var(ZtoNuNu_temp.GetName());
      ws.extendSet(names.nuisances,ZtoNuNu_temp.GetName());
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
  RooAbsArg* ZtomumuSystematic = ws.arg(names.ZtomumuSystematic);
    
  RooAbsArg* ZtoeeEfficiency = ws.arg(names.ZtoeeEfficiency);
  RooAbsArg* ZtoeeInvPurity = ws.arg(names.ZtoeeInvPurity);
  RooAbsArg* ZtoeeSystematic = ws.arg(names.ZtoeeSystematic);  

  RooAbsArg* ZtoeeAcceptance = 
    //getBetaConstraint(ws,"ZtoeeAcceptance_",abcd.ZtoeeAcceptanceName,
    getGaussianConstraint(ws,"ZtoeeAcceptance_",abcd.ZtoeeAcceptanceName,
		      abcd.ZtoeeAcceptance,abcd.ZtoeeAcceptanceError,
		      names.observables,names.nuisances, names.globalObservables);
  
  RooAbsArg* ZtomumuAcceptance = 
    //getBetaConstraint(ws,"ZtomumuAcceptance_",abcd.ZtomumuAcceptanceName,
    getGaussianConstraint(ws,"ZtomumuAcceptance_",abcd.ZtomumuAcceptanceName,
		      abcd.ZtomumuAcceptance,abcd.ZtomumuAcceptanceError,
		      names.observables,names.nuisances, names.globalObservables);
  
  
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
  RooAbsArg* zeroLeptonZtoNuNubTagScaling = getGaussianConstraint(ws,"zeroLeptonZtoNuNubTagScaling_",abcd.ZtoNuNubTagScalingName,
							      abcd.ZtoNuNubTagScaling,abcd.ZtoNuNubTagScalingError,
							      names.observables,names.nuisances, names.globalObservables);
  
  RooProduct zeroLeptonZtoNuNuYield(zeroLeptonName+"_ZtoNuNuYield",zeroLeptonName+"_ZtoNuNuYield",RooArgSet(*zeroLeptonZtoNuNubTagScaling,*ZtoNuNu));
  
  
  //////////////////////
  //  Diboson
  //////////////////////
  RooAbsArg* dibosonMCUncertainty = ws.arg(names.dibosonMCUncertainty);
  
  RooRealVar zeroLeptonDibosonMCYield(zeroLeptonName+"_DibosonMCYield", zeroLeptonName+"_DibosonMCYield", abcd.zeroLeptonDibosonMC);
  RooRealVar zeroLeptonLowDeltaPhiNDibosonMCYield(zeroLeptonLowDeltaPhiNName+"_DibosonMCYield", zeroLeptonLowDeltaPhiNName+"_DibosonMCYield", abcd.zeroLeptonLowDeltaPhiNDibosonMC);
  RooRealVar oneLeptonDibosonMCYield(oneLeptonName+"_DibosonMCYield", oneLeptonName+"_DibosonMCYield", abcd.oneLeptonDibosonMC);
  zeroLeptonDibosonMCYield.setConstant();
  zeroLeptonLowDeltaPhiNDibosonMCYield.setConstant();
  oneLeptonDibosonMCYield.setConstant();
  
  RooProduct zeroLeptonDibosonYield(zeroLeptonName+"_DibosonYield", zeroLeptonName+"_DibosonYield", RooArgSet(*dibosonMCUncertainty,zeroLeptonDibosonMCYield));
  RooProduct zeroLeptonLowDeltaPhiNDibosonYield(zeroLeptonLowDeltaPhiNName+"_DibosonYield", zeroLeptonLowDeltaPhiNName+"_DibosonYield", RooArgSet(*dibosonMCUncertainty,zeroLeptonLowDeltaPhiNDibosonMCYield));
  RooProduct oneLeptonDibosonYield(oneLeptonName+"_DibosonYield", oneLeptonName+"_DibosonYield", RooArgSet(*dibosonMCUncertainty,oneLeptonDibosonMCYield));
  
  
  ///////////////////////////////////////////////////////////////////////////
  //  Non-QCD SM in Low Delta Phi_N region (TopWJets, ZtoNuNu, Diboson)
  ///////////////////////////////////////////////////////////////////////////
  RooRealVar topWJetsLowDeltaPhiNOverZeroLeptonRatioMC(zeroLeptonLowDeltaPhiNName+"_TopWJetsLowDeltaPhiNOverZeroLeptonRatioMC",
						       zeroLeptonLowDeltaPhiNName+"_TopWJetsLowDeltaPhiNOverZeroLeptonRatioMC", 
						       abcd.topWJetsLowDeltaPhiNOverZeroLeptonRatioMC);
  topWJetsLowDeltaPhiNOverZeroLeptonRatioMC.setConstant();
  RooProduct zeroLeptonLowDeltaPhiNTopWJetsYield(zeroLeptonLowDeltaPhiNName+"_TopWJetsYield",zeroLeptonLowDeltaPhiNName+"_TopWJetsYield",RooArgSet(topWJetsLowDeltaPhiNOverZeroLeptonRatioMC,*zeroLeptonTopWJetsYield));
  
  RooRealVar ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC(zeroLeptonLowDeltaPhiNName+"_ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC",
						      zeroLeptonLowDeltaPhiNName+"_ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC", 
						      abcd.ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC);
  ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC.setConstant();
  RooProduct zeroLeptonLowDeltaPhiNZtoNuNuYield(zeroLeptonLowDeltaPhiNName+"_ZtoNuNuYield",zeroLeptonLowDeltaPhiNName+"_ZtoNuNuYield",RooArgSet(ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC,zeroLeptonZtoNuNuYield));
  
  RooAddition zeroLeptonLowDeltaPhiNNonQCDYield_NoSystematic(zeroLeptonLowDeltaPhiNName+"_NonQCDYield_NoSystematic",zeroLeptonLowDeltaPhiNName+"_NonQCDYield_NoSystematic",RooArgSet(zeroLeptonLowDeltaPhiNTopWJetsYield, zeroLeptonLowDeltaPhiNZtoNuNuYield, zeroLeptonLowDeltaPhiNDibosonYield));
  
  RooAbsArg* MCUncertainty = ws.arg(names.MCUncertainty);
  RooProduct zeroLeptonLowDeltaPhiNNonQCDYield(zeroLeptonLowDeltaPhiNName+"_NonQCDYield", zeroLeptonLowDeltaPhiNName+"_NonQCDYield", RooArgSet(zeroLeptonLowDeltaPhiNNonQCDYield_NoSystematic, *MCUncertainty));
  
  ////////////////////////////////////////////////
  // Constrain yields
  ////////////////////////////////////////////////
  
  // Get signal yields from workspace
  RooAbsReal* zeroLeptonSignalYield = ws.function(zeroLeptonName+"_SignalYield");
  RooAbsReal* zeroLeptonLowDeltaPhiNSignalYield = ws.function(zeroLeptonLowDeltaPhiNName+"_SignalYield");
  assert(zeroLeptonSignalYield);
  assert(zeroLeptonLowDeltaPhiNSignalYield);
  
  //Multiply yields by trigger efficiency to make "data yield"
  RooProduct zeroLeptonTopWJetsDataYield(zeroLeptonName+"_TopWJetsDataYield", zeroLeptonName+"_TopWJetsDataYield", RooArgSet(*zeroLeptonTopWJetsYield, *oneLeptonTriggerEfficiency));
  RooProduct zeroLeptonQCDDataYield(zeroLeptonName+"_QCDDataYield", zeroLeptonName+"_QCDDataYield", RooArgSet(zeroLeptonQCDYield, *zeroLeptonTriggerEfficiency));
  RooProduct zeroLeptonZtoNuNuDataYield(zeroLeptonName+"_ZtoNuNuDataYield", zeroLeptonName+"_ZtoNuNuDataYield", RooArgSet(zeroLeptonZtoNuNuYield, *oneLeptonTriggerEfficiency));
  RooProduct zeroLeptonDibosonDataYield(zeroLeptonName+"_DibosonDataYield", zeroLeptonName+"_DibosonDataYield", RooArgSet(zeroLeptonDibosonYield, *oneLeptonTriggerEfficiency));
  RooProduct zeroLeptonSignalDataYield(zeroLeptonName+"_SignalDataYield", zeroLeptonName+"_SignalDataYield", RooArgSet(*zeroLeptonSignalYield, *oneLeptonTriggerEfficiency));

  RooProduct zeroLeptonLowDeltaPhiNQCDDataYield(zeroLeptonLowDeltaPhiNName+"_QCDDataYield", zeroLeptonLowDeltaPhiNName+"_QCDDataYield", RooArgSet(zeroLeptonLowDeltaPhiNQCDYield, *zeroLeptonTriggerEfficiency));
  RooProduct zeroLeptonLowDeltaPhiNNonQCDDataYield(zeroLeptonLowDeltaPhiNName+"_NonQCDDataYield", zeroLeptonLowDeltaPhiNName+"_NonQCDDataYield", RooArgSet(zeroLeptonLowDeltaPhiNNonQCDYield, *oneLeptonTriggerEfficiency));
  RooProduct zeroLeptonLowDeltaPhiNSignalDataYield(zeroLeptonLowDeltaPhiNName+"_SignalDataYield", zeroLeptonLowDeltaPhiNName+"_SignalDataYield", RooArgSet(*zeroLeptonLowDeltaPhiNSignalYield, *oneLeptonTriggerEfficiency));

  //Some additional parameters that are used to make plots later.  Consider rewriting  LDP components to use these.
  RooProduct zeroLeptonLowDeltaPhiNTopWJetsDataYield(zeroLeptonLowDeltaPhiNName+"_TopWJetsDataYield", zeroLeptonLowDeltaPhiNName+"_TopWJetsDataYield", RooArgSet(zeroLeptonLowDeltaPhiNTopWJetsYield, *MCUncertainty, *oneLeptonTriggerEfficiency));
  RooProduct zeroLeptonLowDeltaPhiNZtoNuNuDataYield(zeroLeptonLowDeltaPhiNName+"_ZtoNuNuDataYield", zeroLeptonLowDeltaPhiNName+"_ZtoNuNuDataYield", RooArgSet(zeroLeptonLowDeltaPhiNZtoNuNuYield, *MCUncertainty, *oneLeptonTriggerEfficiency));
  RooProduct zeroLeptonLowDeltaPhiNDibosonDataYield(zeroLeptonLowDeltaPhiNName+"_DibosonDataYield", zeroLeptonLowDeltaPhiNName+"_DibosonDataYield", RooArgSet(zeroLeptonLowDeltaPhiNDibosonYield, *MCUncertainty, *oneLeptonTriggerEfficiency));
  ws.import(zeroLeptonLowDeltaPhiNTopWJetsDataYield, RecycleConflictNodes());
  ws.import(zeroLeptonLowDeltaPhiNZtoNuNuDataYield, RecycleConflictNodes());
  ws.import(zeroLeptonLowDeltaPhiNDibosonDataYield, RecycleConflictNodes());

  //sum data yields
  RooAddition zeroLeptonYield(zeroLeptonName+"_Yield", zeroLeptonName+"_Yield", RooArgSet(zeroLeptonTopWJetsDataYield, zeroLeptonQCDDataYield, zeroLeptonZtoNuNuDataYield, zeroLeptonDibosonDataYield, zeroLeptonSignalDataYield));
  RooAddition zeroLeptonLowDeltaPhiNYield(zeroLeptonLowDeltaPhiNName+"_Yield", zeroLeptonLowDeltaPhiNName+"_Yield", RooArgSet(zeroLeptonLowDeltaPhiNQCDDataYield, zeroLeptonLowDeltaPhiNNonQCDDataYield, zeroLeptonLowDeltaPhiNSignalDataYield));
      
  // Define poisson constraints
  RooPoissonLogEval zeroLeptonConstraint(zeroLeptonName+"_Constraint",zeroLeptonName+"_Constraint",zeroLeptonCount,zeroLeptonYield);
  RooPoissonLogEval zeroLeptonLowDeltaPhiNConstraint(zeroLeptonLowDeltaPhiNName+"_Constraint",zeroLeptonLowDeltaPhiNName+"_Constraint",zeroLeptonLowDeltaPhiNCount,zeroLeptonLowDeltaPhiNYield);
  ws.import(zeroLeptonConstraint, RecycleConflictNodes());
  ws.import(zeroLeptonLowDeltaPhiNConstraint, RecycleConflictNodes());

  //Same for one lepton sample if using nominal method
  if(options.TopWJetsMethod == "ABCD")
    {
      RooAbsReal* oneLeptonSignalYield = ws.function(oneLeptonName+"_SignalYield");
      RooRealVar* oneLeptonTopWJetsYield = ws.var(oneLeptonName+"_TopWJetsYield");
      assert(oneLeptonSignalYield != NULL);
      assert(oneLeptonTopWJetsYield != NULL);

      RooRealVar oneLeptonCount(oneLeptonName+"_Count",oneLeptonName+"_Count",observed.oneLepton);
      oneLeptonCount.setConstant();
      ws.import(oneLeptonCount);
      ws.extendSet(names.observables,oneLeptonCount.GetName());
      
      RooProduct oneLeptonTopWJetsDataYield(oneLeptonName+"_TopWJetsDataYield", oneLeptonName+"_TopWJetsDataYield", RooArgSet(*oneLeptonTopWJetsYield, *oneLeptonTriggerEfficiency));
      RooProduct oneLeptonDibosonDataYield(oneLeptonName+"_DibosonDataYield", oneLeptonName+"_DibosonDataYield", RooArgSet(oneLeptonDibosonYield, *oneLeptonTriggerEfficiency));
      RooProduct oneLeptonSignalDataYield(oneLeptonName+"_SignalDataYield", oneLeptonName+"_SignalDataYield", RooArgSet(*oneLeptonSignalYield, *oneLeptonTriggerEfficiency));

      RooAddition oneLeptonYield(oneLeptonName+"_Yield", oneLeptonName+"_Yield", RooArgSet(oneLeptonTopWJetsDataYield, oneLeptonDibosonDataYield, oneLeptonSignalDataYield));  
      
      RooPoissonLogEval oneLeptonConstraint(oneLeptonName+"_Constraint",oneLeptonName+"_Constraint",oneLeptonCount,oneLeptonYield);    
      
      ws.import(oneLeptonConstraint, RecycleConflictNodes());
  }  
  
  return true;
}


void makeUnderlyingLikelihood(const likelihoodOptions options, RooWorkspace& ws ,allBinNames& names, allBins& numbers)
{
  names.observables = "observables";
  ws.defineSet(names.observables,"");

  names.nuisances = "nuisances";
  ws.defineSet(names.nuisances,"");

  names.globalObservables = "globalObservables";
  ws.defineSet(names.globalObservables,"");



  //Universal parameters
  
  RooAbsArg* dibosonMCUncertainty = 
    getGaussianConstraint(ws,"dibosonMCUncertainty","",
			  numbers.dibosonMC, numbers.dibosonMCUncertainty,
			  names.observables,names.nuisances, names.globalObservables);
  names.dibosonMCUncertainty = dibosonMCUncertainty->GetName();
  
  RooAbsArg* MCUncertainty = 
    //getBetaPrimeConstraint(ws,"MCUncertainty","",
    getGaussianConstraint(ws,"MCUncertainty","",
			   numbers.MC,numbers.MCUncertainty,
			   names.observables,names.nuisances, names.globalObservables);
  names.MCUncertainty = MCUncertainty->GetName();
    
  RooRealVar luminosity("luminosity","luminosity",numbers.Luminosity,0.0,1e2);
  luminosity.setConstant();//should eventually have error
  ws.import(luminosity);
  
  RooRealVar signalCrossSection("signalCrossSection","signalCrossSection",0.0,0.0,1e3);
  ws.import(signalCrossSection);
  names.signalCrossSection = signalCrossSection.GetName();

  double globalUncertainty = sqrt(numbers.LuminosityError*numbers.LuminosityError + numbers.metCleaningError*numbers.metCleaningError);
  RooAbsArg* signalGlobalUncertainty = 
    getGaussianConstraint(ws,"signalGlobalUncertainty","",
			  1.0, globalUncertainty,
			  names.observables,names.nuisances, names.globalObservables);
  names.signalGlobalUncertainty = signalGlobalUncertainty->GetName();


  if(options.TopWJetsMethod == "ABCD")
    {
      RooRealVar singleLeptonScaling("singleLeptonScaling","singleLeptonScaling",0.0,5);
      ws.import(singleLeptonScaling);
      ws.extendSet(names.nuisances,singleLeptonScaling.GetName());
      names.singleLeptonScaling = singleLeptonScaling.GetName();
    }

  //Objects for Z->invisible background:

  RooRealVar ZtollOverZtoNuNuRatio("ZtollOverZtoNuNuRatio","ZtollOverZtoNuNuRatio",numbers.ZtollOverZtoNuNuRatio);
  ZtollOverZtoNuNuRatio.setConstant();
  ws.import(ZtollOverZtoNuNuRatio);
  names.ZtollOverZtoNuNuRatio = ZtollOverZtoNuNuRatio.GetName();

  // RooAbsArg* ZtomumuInvPurity = 
  //   getInverseBetaConstraint(ws,"ZtomumuInvPurity","",
  // 			     numbers.ZtomumuPurity,numbers.ZtomumuPurityError,
  // 			     names.observables,names.nuisances, names.globalObservables);
  // names.ZtomumuInvPurity = ZtomumuInvPurity->GetName();
  
  // RooAbsArg* ZtoeeInvPurity = 
  //   getInverseBetaConstraint(ws,"ZtoeeInvPurity","",
  // 			     numbers.ZtoeePurity,numbers.ZtoeePurityError,
  // 			     names.observables,names.nuisances, names.globalObservables);
  // names.ZtoeeInvPurity = ZtoeeInvPurity->GetName();

  RooRealVar* ZtomumuPurity = (RooRealVar*)
    getGaussianConstraint(ws,"ZtomumuPurity","",
			  numbers.ZtomumuPurity,numbers.ZtomumuPurityError,
			  names.observables,names.nuisances, names.globalObservables);
  RooRatio ZtomumuInvPurity("ZtomumuInvPurity", "ZtomumuInvPurity", RooConst(1.), *ZtomumuPurity);
  ws.import(ZtomumuInvPurity);
  names.ZtomumuInvPurity = ZtomumuInvPurity.GetName();

  RooRealVar* ZtoeePurity = (RooRealVar*)
    getGaussianConstraint(ws,"ZtoeePurity","",
			  numbers.ZtoeePurity,numbers.ZtoeePurityError,
			  names.observables,names.nuisances, names.globalObservables);
  RooRatio ZtoeeInvPurity("ZtoeeInvPurity", "ZtoeeInvPurity", RooConst(1.), *ZtoeePurity);
  ws.import(ZtoeeInvPurity);
  names.ZtoeeInvPurity = ZtoeeInvPurity.GetName();
  

  RooAbsArg* ZtomumuEfficiency = 
    //getBetaConstraint(ws,"ZtomumuEfficiency","",
    getGaussianConstraint(ws,"ZtomumuEfficiency","",
		      numbers.ZtomumuEfficiency,numbers.ZtomumuEfficiencyError,
		      names.observables,names.nuisances, names.globalObservables);
  names.ZtomumuEfficiency = ZtomumuEfficiency->GetName();

  RooAbsArg* ZtoeeEfficiency = 
    //getBetaConstraint(ws,"ZtoeeEfficiency","",	
    getGaussianConstraint(ws,"ZtoeeEfficiency","",
		      numbers.ZtoeeEfficiency,numbers.ZtoeeEfficiencyError,
		      names.observables,names.nuisances, names.globalObservables);
  names.ZtoeeEfficiency = ZtoeeEfficiency->GetName();

  RooAbsArg* ZtoeeSystematic = 
    getCorrelatedGaussianConstraint(ws,"ZtoeeSystematic","",
				    numbers.ZtoeeSystematic,numbers.ZtoeeSystematicError,
				    names.observables,names.nuisances, names.globalObservables,
				    "ZtollSystematic");
  names.ZtoeeSystematic = ZtoeeSystematic->GetName();
  
  RooAbsArg* ZtomumuSystematic = 
    getCorrelatedGaussianConstraint(ws,"ZtomumuSystematic","",
				    numbers.ZtomumuSystematic,numbers.ZtomumuSystematicError,
				    names.observables,names.nuisances, names.globalObservables,
				    "ZtollSystematic");
  names.ZtomumuSystematic = ZtomumuSystematic->GetName();
  
}



void setupSignalModelMR(const TString binName, const TString binFileNameInsideSignalMR, map<TString,map<TString,mcCount> >& insideSignalMR, 
			const TString binFileNameOutsideSignalMR, map<TString,map<TString,mcCount> >& outsideSignalMR )
{
  
  cout << "Reading in MR signal inputs for " << binName << endl;
    
  ifstream insideFile;
  
  cout << "getting the file: " << binFileNameInsideSignalMR << endl;
  
  insideFile.open(binFileNameInsideSignalMR.Data(),fstream::in);
  assert(insideFile.is_open());
  
  string fileLine;
  
  TString index;
  double value;
  
  map<TString,double> thisInsideSignalMRValue;
  map<TString,double> thisInsideSignalMRError;

  while(!insideFile.eof())
    {

      getline(insideFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken nameAndNumber(thisLine," ");
      nameAndNumber.NextToken();
      index = nameAndNumber;
      if(index == "") continue;
      nameAndNumber.NextToken();
      value = nameAndNumber.Atof();
      cout << index << " : " << value << endl;
      
      if(index == "oneTightMu_Theta1_SignalCount"       ) { thisInsideSignalMRValue.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta1", value) ); }
      else if(index == "oneLooseLep_Theta1_SignalCount" ) { thisInsideSignalMRValue.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta1", value) ); }
      else if(index == "oneTightMu_Theta2_SignalCount"  ) { thisInsideSignalMRValue.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta2", value) ); } 
      else if(index == "oneLooseLep_Theta2_SignalCount" ) { thisInsideSignalMRValue.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta2", value) ); }
      else if(index == "oneTightMu_Theta3_SignalCount"  ) { thisInsideSignalMRValue.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta3", value) ); }
      else if(index == "oneLooseLep_Theta3_SignalCount" ) { thisInsideSignalMRValue.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta3", value) ); }
      else if(index == "oneTightMu_Theta4_SignalCount"  ) { thisInsideSignalMRValue.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta4", value) ); }
      else if(index == "oneLooseLep_Theta4_SignalCount" ) { thisInsideSignalMRValue.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta4", value) ); }
      else if(index == "oneTightMu_Theta5_SignalCount"  ) { thisInsideSignalMRValue.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta5", value) ); }
      else if(index == "oneLooseLep_Theta5_SignalCount" ) { thisInsideSignalMRValue.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta5", value) ); }
      else if(index == "twoTightMu_SignalCount"         ) { thisInsideSignalMRValue.insert( pair<TString,int>("twoTightMu_"+binName, value) ); }
      else if(index == "twoLooseLep_SignalCount"        ) { thisInsideSignalMRValue.insert( pair<TString,int>("twoLooseLep_"+binName, value) ); }

      else if(index == "oneTightMu_Theta1_SignalCountError"  ) { thisInsideSignalMRError.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta1", value) ); }
      else if(index == "oneLooseLep_Theta1_SignalCountError" ) { thisInsideSignalMRError.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta1", value) ); }
      else if(index == "oneTightMu_Theta2_SignalCountError"  ) { thisInsideSignalMRError.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta2", value) ); } 
      else if(index == "oneLooseLep_Theta2_SignalCountError" ) { thisInsideSignalMRError.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta2", value) ); }
      else if(index == "oneTightMu_Theta3_SignalCountError"  ) { thisInsideSignalMRError.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta3", value) ); }
      else if(index == "oneLooseLep_Theta3_SignalCountError" ) { thisInsideSignalMRError.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta3", value) ); }
      else if(index == "oneTightMu_Theta4_SignalCountError"  ) { thisInsideSignalMRError.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta4", value) ); }
      else if(index == "oneLooseLep_Theta4_SignalCountError" ) { thisInsideSignalMRError.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta4", value) ); }
      else if(index == "oneTightMu_Theta5_SignalCountError"  ) { thisInsideSignalMRError.insert( pair<TString,int>("oneTightMu_"+binName+"_Theta5", value) ); }
      else if(index == "oneLooseLep_Theta5_SignalCountError" ) { thisInsideSignalMRError.insert( pair<TString,int>("oneLooseLep_"+binName+"_Theta5", value) ); }
      else if(index == "twoTightMu_SignalCountError"         ) { thisInsideSignalMRError.insert( pair<TString,int>("twoTightMu_"+binName, value) ); }
      else if(index == "twoLooseLep_SignalCountError"        ) { thisInsideSignalMRError.insert( pair<TString,int>("twoLooseLep_"+binName, value) ); }

      else {assert(0);}
    }
  insideFile.close();
  
  map<TString,mcCount> thisInsideSignalMR;
  for( map<TString,double>::iterator it=thisInsideSignalMRValue.begin(); it!=thisInsideSignalMRValue.end(); it++ )
    {
      TString name = (*it).first;
      mcCount thisCount;
      thisCount.value = (*it).second;
      thisCount.error = thisInsideSignalMRError[name];
      thisInsideSignalMR.insert( pair<TString,mcCount>(name, thisCount) );
    }
  //for( map<TString,mcCount>::iterator it=thisInsideSignalMR.begin(); it!=thisInsideSignalMR.end(); it++ ) cout << "Map: " << (*it).first << " = " << ((*it).second).value << " +- " << ((*it).second).error << endl;
  insideSignalMR.insert( pair<TString, map<TString,mcCount> >(binName, thisInsideSignalMR) );


  //Now the outside values
  
  ifstream outsideFile;
  
  cout << "getting the file: " << binFileNameOutsideSignalMR << endl;
  
  outsideFile.open(binFileNameOutsideSignalMR.Data(),fstream::in);
  assert(outsideFile.is_open());
  
  TString oneMuName = "", twoMuName = "";
  mcCount oneMu, twoMu;
  oneMu.value=0; oneMu.error=0; twoMu.value=0; twoMu.error=0;
  
  TString svalue;
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
      
      if(index == "oneTightMu_Outside_SignalCountName"        ) { oneMuName = svalue; }
      else if(index == "oneTightMu_Outside_SignalCount"       ) { oneMu.value = svalue.Atof(); }
      else if(index == "oneTightMu_Outside_SignalCountError"  ) { oneMu.error = svalue.Atof(); }
      else if(index == "twoTightMu_Outside_SignalCountName"   ) { twoMuName = svalue; }
      else if(index == "twoTightMu_Outside_SignalCount"       ) { twoMu.value = svalue.Atof(); }
      else if(index == "twoTightMu_Outside_SignalCountError"  ) { twoMu.error = svalue.Atof(); }
      else{ assert(0); }
    }
  outsideFile.close();

  map<TString,mcCount> thisOutsideSignalMR;
  thisOutsideSignalMR.insert( pair<TString,mcCount>("oneTightMu_"+oneMuName+"_Outside",oneMu) );
  thisOutsideSignalMR.insert( pair<TString,mcCount>("twoTightMu_"+twoMuName+"_Outside",twoMu) );
  
  //for( map<TString,mcCount>::iterator it=thisOutsideSignalMR.begin(); it!=thisOutsideSignalMR.end(); it++ ) cout << "Map: " << (*it).first << " = " << ((*it).second).value << " +- " << ((*it).second).error << endl;
  
  outsideSignalMR.insert( pair<TString, map<TString,mcCount> >(binName, thisOutsideSignalMR) );
  
}


//To copy OAK's signal modeling 
void setupSignalModelOAK( vector<TString> binNames, TString signalModelFilesPath, double &nGenerated,
			  map<TString,yields> &signalFractionsOAK, map<TString,yields> &signalStatisticalErrorOAK, map<TString,yields> &signalBTagEfficiencyErrorOAK, map<TString,yields> &signalJesErrorOAK)
{

  //fractions and statistical error

  TString signalModelFileName = signalModelFilesPath;
  signalModelFileName += "/subset_sigcounts.txt";
  
  ifstream signalFile;
  
  signalFile.open(signalModelFileName.Data(),fstream::in);
  assert(signalFile.is_open());

  double ArrayContent[3+4*binNames.size()];//no SL
  for(unsigned int i = 0; signalFile && i< (3+4*binNames.size()); i++)
    {
      signalFile >> ArrayContent[i];
    }

  cout << "OAK signal model, m0 = " << ArrayContent[0] << ", m12 = " << ArrayContent[1] << ", nGenerated = " << ArrayContent[2] << endl;
  
  nGenerated = ArrayContent[2];

  double zeroLeptonTotal = 0;
  cout << endl;
  cout << "Fractions and Statistical error" << endl;
  for(unsigned int i = 0; i < binNames.size(); i++)
    {
      TString binName = binNames[i];
      
      yields thisSignalFractionsOAK;
      yields thisSignalStatisticalErrorOAK;
      
      zeroLeptonTotal += ArrayContent[3 + i];
      
      thisSignalFractionsOAK.zeroLepton = ArrayContent[3 + i] / nGenerated;
      thisSignalFractionsOAK.oneLepton = 0.0;
      thisSignalFractionsOAK.zeroLeptonLowDeltaPhiN = ArrayContent[3 + binNames.size() + i] / nGenerated;
      
      thisSignalStatisticalErrorOAK.zeroLepton = (ArrayContent[3 + i] > 1e-5) ? ArrayContent[3 + 2*binNames.size() + i]/ArrayContent[3 + i] : 0.10/100.0;
      thisSignalStatisticalErrorOAK.oneLepton = 0.10/100.0;
      thisSignalStatisticalErrorOAK.zeroLeptonLowDeltaPhiN = (ArrayContent[3 + binNames.size() + i] > 1e-5) ? ArrayContent[3 + 3*binNames.size() + i]/ArrayContent[3 + binNames.size() + i] : 0.10/100.0;
      cout << binName << " " << thisSignalFractionsOAK.zeroLepton << " +- " << thisSignalStatisticalErrorOAK.zeroLepton << ", " 
	   << thisSignalFractionsOAK.oneLepton << " +- " << thisSignalStatisticalErrorOAK.oneLepton << ", " 
	   << thisSignalFractionsOAK.zeroLeptonLowDeltaPhiN << " +- " << thisSignalStatisticalErrorOAK.zeroLeptonLowDeltaPhiN << endl;
      
      signalFractionsOAK.insert( pair<TString,yields>(binName, thisSignalFractionsOAK) );
      signalStatisticalErrorOAK.insert( pair<TString,yields>(binName, thisSignalStatisticalErrorOAK) );
    }

  signalFile.close();

  cout << "Total zero lepton = " << zeroLeptonTotal;

  
  //btag eff
  cout << endl;
  cout << "bTag Efficiency" << endl; 
  TString bTagEfficiencyFileName = signalModelFilesPath;
  bTagEfficiencyFileName += "/subset_btagEff.txt";
  
  ifstream bTagEfficiencyFile;

  bTagEfficiencyFile.open(bTagEfficiencyFileName.Data(),fstream::in);
  assert(bTagEfficiencyFile.is_open());

  double bTagEfficiencyArrayContent[2+3*binNames.size()];
  for(unsigned int i = 0; bTagEfficiencyFile && i< (2+3*binNames.size()); i++)
    {
      bTagEfficiencyFile >> bTagEfficiencyArrayContent[i];
    }

  for(unsigned int i = 0; i < binNames.size(); i++)
    {
      TString binName = binNames[i];
      yields thisSignalBTagEfficiencyErrorOAK;
      
      thisSignalBTagEfficiencyErrorOAK.zeroLepton = bTagEfficiencyArrayContent[2 + i]; 
      thisSignalBTagEfficiencyErrorOAK.oneLepton = bTagEfficiencyArrayContent[2 + binNames.size() + i];
      thisSignalBTagEfficiencyErrorOAK.zeroLeptonLowDeltaPhiN = bTagEfficiencyArrayContent[2 + 2*binNames.size() + i]; 
      cout << binName << " " << thisSignalBTagEfficiencyErrorOAK.zeroLepton << " " << thisSignalBTagEfficiencyErrorOAK.oneLepton << " " << thisSignalBTagEfficiencyErrorOAK.zeroLeptonLowDeltaPhiN << endl;

      signalBTagEfficiencyErrorOAK.insert( pair<TString,yields>(binName, thisSignalBTagEfficiencyErrorOAK) );
    }
  bTagEfficiencyFile.close();
  

  //JES
  cout << endl;
  cout << "JES" << endl;
  TString jesFileName = signalModelFilesPath;
  jesFileName += "/subset_JES.txt";

  ifstream jesFile;

  jesFile.open(jesFileName.Data(),fstream::in);
  assert(jesFile.is_open());

  double jesArrayContent[2+2*binNames.size()];
  for(unsigned int i = 0; jesFile && i < (2+2*binNames.size()); i++)
    {
      jesFile >> jesArrayContent[i];
    }

  for(unsigned int i = 0; i < binNames.size(); i++)
    {
      TString binName = binNames[i];
      yields thisJesErrorOAK;
      
      thisJesErrorOAK.zeroLepton = jesArrayContent[2 + i];
      thisJesErrorOAK.oneLepton = 0.0;
      thisJesErrorOAK.zeroLeptonLowDeltaPhiN = jesArrayContent[2 + binNames.size() + i];
      cout << binName << " " << thisJesErrorOAK.zeroLepton << " " << thisJesErrorOAK.oneLepton << " " << thisJesErrorOAK.zeroLeptonLowDeltaPhiN << endl; 
      
      signalJesErrorOAK.insert( pair<TString,yields>(binName, thisJesErrorOAK) );
    }
}


void makeSignalModel(const likelihoodOptions options, RooWorkspace& ws , vector<TString> binNames, allBinNames& names,  
		     map<TString,abcdBinParameters> bins, double nGenerated,
		     map<TString,yields> signalFractionsOAK, map<TString,yields> signalStatisticalErrorOAK,
		     map<TString,yields> signalBTagEfficiencyErrorOAK, map<TString,yields> signalJesErrorOAK,
		     map<TString,map<TString,mcCount> > insideSignalMR, map<TString,map<TString,mcCount> > outsideSignalMR )
{
  RooRealVar* luminosity = ws.var("luminosity");
  RooRealVar* signalCrossSection =  ws.var(names.signalCrossSection);
  RooAbsArg* signalGlobalUncertainty = ws.arg(names.signalGlobalUncertainty);

  for(vector<TString>::iterator thisBin = binNames.begin() ; thisBin != binNames.end() ; thisBin++)
    {
      TString binName = *thisBin;
      if( skipBin(binName) ) continue;
      
      if(options.TopWJetsMethod == "metReweighting")
	{
	  //Zero lepton and zero lepton low deltaPhiN
	  
	  TString zeroLeptonName("zeroLepton_");
	  zeroLeptonName+=binName;
	  TString zeroLeptonLowDeltaPhiNName("zeroLeptonLowDeltaPhiN_");
	  zeroLeptonLowDeltaPhiNName+=binName;
	  	  
	  //Fractions
	  yields thisSignalFractionsOAK = signalFractionsOAK[binName];
	  RooRealVar zeroLeptonFraction(zeroLeptonName+"_SignalFractionOAK", zeroLeptonName+"_SignalFractionOAK", thisSignalFractionsOAK.zeroLepton);
	  RooRealVar zeroLeptonLowDeltaPhiNFraction(zeroLeptonLowDeltaPhiNName+"_SignalFractionOAK", zeroLeptonLowDeltaPhiNName+"_SignalFractionOAK", thisSignalFractionsOAK.zeroLeptonLowDeltaPhiN);
	  zeroLeptonFraction.setConstant();
	  zeroLeptonLowDeltaPhiNFraction.setConstant();
	  
	  //Statistical Error
	  yields thisSignalStatisticalErrorOAK = signalStatisticalErrorOAK[binName]; 
	  if(options.nuisanceOption == "noWidths") 
	    {
	      thisSignalStatisticalErrorOAK.zeroLepton = 0; 
	      thisSignalStatisticalErrorOAK.zeroLeptonLowDeltaPhiN = 0;
	    }

	  RooAbsArg* zeroLeptonError = getGaussianConstraint(ws,"zeroLeptonSignalError_", binName,
							     1.0, thisSignalStatisticalErrorOAK.zeroLepton,
							     names.observables,names.nuisances, names.globalObservables);
	  
	  RooAbsArg* zeroLeptonLowDeltaPhiNError = getGaussianConstraint(ws,"zeroLeptonLowDeltaPhiNSignalError_", binName,
									 1.0, thisSignalStatisticalErrorOAK.zeroLeptonLowDeltaPhiN,
									 names.observables,names.nuisances, names.globalObservables);
	  
	  //B-tag efficiency systematic
	  yields thisSignalBTagEfficiencyErrorOAK = signalBTagEfficiencyErrorOAK[binName];
	  TString signalBTagEfficiencyName = "signalBTagEfficiencyCorrelated";
	  if(options.nuisanceOption == "noWidths") 
	    {
	      thisSignalBTagEfficiencyErrorOAK.zeroLepton = 0; 
	      thisSignalBTagEfficiencyErrorOAK.zeroLeptonLowDeltaPhiN = 0;
	    }
	  
	  bool changeSign = false;
	  if(thisSignalBTagEfficiencyErrorOAK.zeroLepton < 0.0) changeSign=true;
	  //RooAbsArg* zeroLeptonBTagEfficiencyError = getCorrelatedGaussianConstraint(ws,"zeroLeptonBTagEfficiencyError_", binName,
	  RooAbsArg* zeroLeptonBTagEfficiencyError = getCorrelatedLogNormalConstraint(ws,"zeroLeptonBTagEfficiencyError_", binName,
										      1.0, fabs(thisSignalBTagEfficiencyErrorOAK.zeroLepton),
										      names.observables,names.nuisances, names.globalObservables,
										      signalBTagEfficiencyName, changeSign);
	  
	  changeSign = false;
	  if(thisSignalBTagEfficiencyErrorOAK.zeroLeptonLowDeltaPhiN < 0.0) changeSign=true;
	  //RooAbsArg* zeroLeptonLowDeltaPhiNBTagEfficiencyError = getCorrelatedGaussianConstraint(ws,"zeroLeptonLowDeltaPhiNBTagEfficiencyError_", binName,
	  RooAbsArg* zeroLeptonLowDeltaPhiNBTagEfficiencyError = getCorrelatedLogNormalConstraint(ws,"zeroLeptonLowDeltaPhiNBTagEfficiencyError_", binName,
												  1.0, fabs(thisSignalBTagEfficiencyErrorOAK.zeroLeptonLowDeltaPhiN),
												  names.observables,names.nuisances, names.globalObservables,
												  signalBTagEfficiencyName, changeSign);
	  
	  //JES systematic
	  yields thisSignalJesErrorOAK = signalJesErrorOAK[binName];
	  TString signalJesErrorName = "signalJesErrorCorrelated";
	  if(options.nuisanceOption == "noWidths") 
	    {
	      thisSignalJesErrorOAK.zeroLepton = 0; 
	      thisSignalJesErrorOAK.zeroLeptonLowDeltaPhiN = 0;
	    }

	  changeSign = false;
	  if(thisSignalJesErrorOAK.zeroLepton < 0.0) changeSign=true;
	  //RooAbsArg* zeroLeptonJesError = getCorrelatedGaussianConstraint(ws,"zeroLeptonJesError_", binName,
	  RooAbsArg* zeroLeptonJesError = getCorrelatedLogNormalConstraint(ws,"zeroLeptonJesError_", binName,
									   1.0, fabs(thisSignalJesErrorOAK.zeroLepton),
									   names.observables,names.nuisances, names.globalObservables,
									   signalJesErrorName, changeSign);
	  
	  changeSign = false;
	  if(thisSignalJesErrorOAK.zeroLeptonLowDeltaPhiN < 0.0) changeSign=true;
	  //RooAbsArg* zeroLeptonLowDeltaPhiNJesError = getCorrelatedGaussianConstraint(ws,"zeroLeptonLowDeltaPhiNJesError_", binName,
	  RooAbsArg* zeroLeptonLowDeltaPhiNJesError = getCorrelatedLogNormalConstraint(ws,"zeroLeptonLowDeltaPhiNJesError_", binName,
										       1.0, fabs(thisSignalJesErrorOAK.zeroLeptonLowDeltaPhiN),
										       names.observables,names.nuisances, names.globalObservables,
										       signalJesErrorName, changeSign);	  
	  	  
	  //Setup yields
	  RooProduct zeroLeptonSignalYieldOAK(zeroLeptonName+"_SignalYield", zeroLeptonName+"_SignalYield", RooArgSet(*luminosity, *signalCrossSection, zeroLeptonFraction, *signalGlobalUncertainty, *zeroLeptonError, *zeroLeptonBTagEfficiencyError, *zeroLeptonJesError) );
	  RooProduct zeroLeptonLowDeltaPhiNSignalYieldOAK(zeroLeptonLowDeltaPhiNName+"_SignalYield", zeroLeptonLowDeltaPhiNName+"_SignalYield", RooArgSet(*luminosity, *signalCrossSection, zeroLeptonLowDeltaPhiNFraction, *signalGlobalUncertainty, *zeroLeptonLowDeltaPhiNError, *zeroLeptonLowDeltaPhiNBTagEfficiencyError, *zeroLeptonLowDeltaPhiNJesError) );
	  	  
	  ws.import(zeroLeptonSignalYieldOAK, RecycleConflictNodes());
	  ws.import(zeroLeptonLowDeltaPhiNSignalYieldOAK, RecycleConflictNodes());
	  
	  
	  //Single lepton

	  //inside
	  map<TString,mcCount> thisInsideSignalMR = insideSignalMR[binName];
	  for ( map<TString,mcCount>::iterator it=thisInsideSignalMR.begin() ; it != thisInsideSignalMR.end(); it++ ) 
	    { 
	      TString signalName = (*it).first;
	      double signalValue = ((*it).second).value;
	      double signalError = ((*it).second).error;
	      double percentError = (signalValue > 0) ? signalError/signalValue : 0.10/100.0;
	      if(options.nuisanceOption == "noWidths") percentError = 0;

	      RooRealVar signalFraction(signalName+"_SignalFraction", signalName+"_SignalFraction", signalValue/nGenerated);
	      signalFraction.setConstant();
	      
	      RooAbsArg* signalStatisticalError = getGaussianConstraint(ws, signalName+"_SignalError", "",
							     1.0, percentError,
							     names.observables, names.nuisances, names.globalObservables);

	      RooProduct signalYield(signalName+"_SignalYield", signalName+"_SignalYield", RooArgSet(*luminosity, *signalCrossSection, signalFraction, *signalGlobalUncertainty, *signalStatisticalError) );
	      ws.import(signalYield, RecycleConflictNodes());
	    }
	  
	  //outside
	  /*
	  map<TString,mcCount> thisOutsideSignalMR = outsideSignalMR[binName];
	  for ( map<TString,mcCount>::iterator it=thisOutsideSignalMR.begin() ; it != thisOutsideSignalMR.end(); it++ ) 
	    { 
	      TString signalName = (*it).first;
	      double signalValue = ((*it).second).value;
	      double signalError = ((*it).second).error;
	      double percentError = (signalValue > 0) ? signalError/signalValue : 0.10/100.0;
	      if(options.nuisanceOption == "noWidths") percentError = 0;

	      //Continue if already in workspace 
	      if( ws.var(signalName+"_SignalYield") != NULL ) continue; 
	      
	      RooRealVar signalFraction(signalName+"_SignalFraction", signalName+"_SignalFraction", signalValue/nGenerated);
	      signalFraction.setConstant();
	      
	      RooAbsArg* signalStatisticalError = getGaussianConstraint(ws, signalName+"_SignalError", "",
							     1.0, percentError,
							     names.observables, names.nuisances, names.globalObservables);
	      
	      RooProduct signalYield(signalName+"_SignalYield", signalName+"_SignalYield", RooArgSet(*luminosity, *signalCrossSection, signalFraction, *signalGlobalUncertainty, *signalStatisticalError) );
	      ws.import(signalYield, RecycleConflictNodes());
	    }
	  */
	}
      else if(options.TopWJetsMethod == "ABCD")
	{
	  TString zeroLeptonName("zeroLepton_");
	  zeroLeptonName+=binName;
	  TString zeroLeptonLowDeltaPhiNName("zeroLeptonLowDeltaPhiN_");
	  zeroLeptonLowDeltaPhiNName+=binName;
	  TString oneLeptonNameOAK("oneLepton_");
	  oneLeptonNameOAK+=binName;
	  
	  //Fractions
	  yields thisSignalFractionsOAK = signalFractionsOAK[binName];
	  RooRealVar zeroLeptonFraction(zeroLeptonName+"_SignalFractionOAK", zeroLeptonName+"_SignalFractionOAK", thisSignalFractionsOAK.zeroLepton);
	  RooRealVar zeroLeptonLowDeltaPhiNFraction(zeroLeptonLowDeltaPhiNName+"_SignalFractionOAK", zeroLeptonLowDeltaPhiNName+"_SignalFractionOAK", thisSignalFractionsOAK.zeroLeptonLowDeltaPhiN);
	  RooRealVar oneLeptonFraction(oneLeptonNameOAK+"_SignalFractionOAK", oneLeptonNameOAK+"_SignalFractionOAK", thisSignalFractionsOAK.oneLepton);
	  zeroLeptonFraction.setConstant();
	  zeroLeptonLowDeltaPhiNFraction.setConstant();
	  oneLeptonFraction.setConstant();
	  
	  //Statistical Error
	  yields thisSignalStatisticalErrorOAK = signalStatisticalErrorOAK[binName];
	  if(options.nuisanceOption == "noWidths") 
	    {
	      thisSignalStatisticalErrorOAK.zeroLepton = 0; 
	      thisSignalStatisticalErrorOAK.zeroLeptonLowDeltaPhiN = 0;
	      thisSignalStatisticalErrorOAK.oneLepton = 0;
	    }

	  RooAbsArg* zeroLeptonError = getGaussianConstraint(ws,"zeroLeptonSignalError_", binName,
							     1.0, thisSignalStatisticalErrorOAK.zeroLepton,
							     names.observables,names.nuisances, names.globalObservables);
	  
	  RooAbsArg* zeroLeptonLowDeltaPhiNError = getGaussianConstraint(ws,"zeroLeptonLowDeltaPhiNSignalError_", binName,
									 1.0, thisSignalStatisticalErrorOAK.zeroLeptonLowDeltaPhiN,
									 names.observables,names.nuisances, names.globalObservables);
	  
	  RooAbsArg* oneLeptonError = getGaussianConstraint(ws,"oneLeptonSignalError_", binName,
							    1.0, thisSignalStatisticalErrorOAK.oneLepton,
							    names.observables,names.nuisances, names.globalObservables);
	  
	  //B-tag efficiency systematic
	  yields thisSignalBTagEfficiencyErrorOAK = signalBTagEfficiencyErrorOAK[binName];
	  TString signalBTagEfficiencyName = "signalBTagEfficiencyCorrelated";
	  if(options.nuisanceOption == "noWidths") 
	    {
	      thisSignalBTagEfficiencyErrorOAK.zeroLepton = 0; 
	      thisSignalBTagEfficiencyErrorOAK.zeroLeptonLowDeltaPhiN = 0;
	      thisSignalBTagEfficiencyErrorOAK.oneLepton = 0;
	    }	  

	  bool changeSign = false;
	  if(thisSignalBTagEfficiencyErrorOAK.zeroLepton < 0.0) changeSign=true;
	  //RooAbsArg* zeroLeptonBTagEfficiencyError = getCorrelatedGaussianConstraint(ws,"zeroLeptonBTagEfficiencyError_", binName,
	  RooAbsArg* zeroLeptonBTagEfficiencyError = getCorrelatedLogNormalConstraint(ws,"zeroLeptonBTagEfficiencyError_", binName,
										      1.0, fabs(thisSignalBTagEfficiencyErrorOAK.zeroLepton),
										      names.observables,names.nuisances, names.globalObservables,
										      signalBTagEfficiencyName, changeSign);
	  
	  changeSign = false;
	  if(thisSignalBTagEfficiencyErrorOAK.zeroLeptonLowDeltaPhiN < 0.0) changeSign=true;
	  //RooAbsArg* zeroLeptonLowDeltaPhiNBTagEfficiencyError = getCorrelatedGaussianConstraint(ws,"zeroLeptonLowDeltaPhiNBTagEfficiencyError_", binName,
	  RooAbsArg* zeroLeptonLowDeltaPhiNBTagEfficiencyError = getCorrelatedLogNormalConstraint(ws,"zeroLeptonLowDeltaPhiNBTagEfficiencyError_", binName,
												  1.0, fabs(thisSignalBTagEfficiencyErrorOAK.zeroLeptonLowDeltaPhiN),
												  names.observables,names.nuisances, names.globalObservables,
												  signalBTagEfficiencyName, changeSign);
	  
	  changeSign = false;
	  if(thisSignalBTagEfficiencyErrorOAK.oneLepton < 0.0) changeSign=true;
	  //RooAbsArg* oneLeptonBTagEfficiencyError = getCorrelatedGaussianConstraint(ws,"oneLeptonBTagEfficiencyError_", binName,
	  RooAbsArg* oneLeptonBTagEfficiencyError = getCorrelatedLogNormalConstraint(ws,"oneLeptonBTagEfficiencyError_", binName,
										     1.0, fabs(thisSignalBTagEfficiencyErrorOAK.oneLepton),
										     names.observables,names.nuisances, names.globalObservables,
										     signalBTagEfficiencyName, changeSign);
	  
	  //JES systematic
	  yields thisSignalJesErrorOAK = signalJesErrorOAK[binName];
	  TString signalJesErrorName = "signalJesErrorCorrelated";
	  if(options.nuisanceOption == "noWidths") 
	    {
	      thisSignalJesErrorOAK.zeroLepton = 0; 
	      thisSignalJesErrorOAK.zeroLeptonLowDeltaPhiN = 0;
	      thisSignalJesErrorOAK.oneLepton = 0;
	    }	 

	  changeSign = false;
	  if(thisSignalJesErrorOAK.zeroLepton < 0.0) changeSign=true;
	  //RooAbsArg* zeroLeptonJesError = getCorrelatedGaussianConstraint(ws,"zeroLeptonJesError_", binName,
	  RooAbsArg* zeroLeptonJesError = getCorrelatedLogNormalConstraint(ws,"zeroLeptonJesError_", binName,
									   1.0, fabs(thisSignalJesErrorOAK.zeroLepton),
									   names.observables,names.nuisances, names.globalObservables,
									   signalJesErrorName, changeSign);
	  
	  changeSign = false;
	  if(thisSignalJesErrorOAK.zeroLeptonLowDeltaPhiN < 0.0) changeSign=true;
	  //RooAbsArg* zeroLeptonLowDeltaPhiNJesError = getCorrelatedGaussianConstraint(ws,"zeroLeptonLowDeltaPhiNJesError_", binName,
	  RooAbsArg* zeroLeptonLowDeltaPhiNJesError = getCorrelatedLogNormalConstraint(ws,"zeroLeptonLowDeltaPhiNJesError_", binName,
										       1.0, fabs(thisSignalJesErrorOAK.zeroLeptonLowDeltaPhiN),
										       names.observables,names.nuisances, names.globalObservables,
										       signalJesErrorName, changeSign);
	  
	  changeSign = false;
	  if(thisSignalJesErrorOAK.oneLepton < 0.0) changeSign=true;
	  //RooAbsArg* oneLeptonJesError = getCorrelatedGaussianConstraint(ws,"oneLeptonJesError_", binName,
	  RooAbsArg* oneLeptonJesError = getCorrelatedLogNormalConstraint(ws,"oneLeptonJesError_", binName,
									  1.0, fabs(thisSignalJesErrorOAK.oneLepton),
									  names.observables,names.nuisances, names.globalObservables,
									  signalJesErrorName, changeSign);
	  
	  
	  //Setup yields
	  RooProduct zeroLeptonSignalYieldOAK(zeroLeptonName+"_SignalYield", zeroLeptonName+"_SignalYield", RooArgSet(*luminosity, *signalCrossSection, zeroLeptonFraction, *signalGlobalUncertainty, *zeroLeptonError, *zeroLeptonBTagEfficiencyError, *zeroLeptonJesError) );
	  RooProduct zeroLeptonLowDeltaPhiNSignalYieldOAK(zeroLeptonLowDeltaPhiNName+"_SignalYield", zeroLeptonLowDeltaPhiNName+"_SignalYield", RooArgSet(*luminosity, *signalCrossSection, zeroLeptonLowDeltaPhiNFraction, *signalGlobalUncertainty, *zeroLeptonLowDeltaPhiNError, *zeroLeptonLowDeltaPhiNBTagEfficiencyError, *zeroLeptonLowDeltaPhiNJesError) );
	  RooProduct oneLeptonSignalYieldOAK(oneLeptonNameOAK+"_SignalYield", oneLeptonNameOAK+"_SignalYield", RooArgSet(*luminosity, *signalCrossSection, oneLeptonFraction, *signalGlobalUncertainty, *oneLeptonError, *oneLeptonBTagEfficiencyError, *oneLeptonJesError) );
	  
	  ws.import(zeroLeptonSignalYieldOAK, RecycleConflictNodes());
	  ws.import(zeroLeptonLowDeltaPhiNSignalYieldOAK, RecycleConflictNodes());
	  ws.import(oneLeptonSignalYieldOAK, RecycleConflictNodes());
	  
	}
      else { assert(0); }
      
    }//end loop over bins
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
      else if(index == "diElectronCountName"                          ) counts.diElectronName = value;	     
      else if(index == "diElectronCount"                              ) counts.diElectron = value.Atof();	     
      else if(index == "diMuonCountName"	                      ) counts.diMuonName = value;	      
      else if(index == "diMuonCount"	                              ) counts.diMuon = value.Atof();	      
      else if(index == "zeroLeptonTriggerEfficiencyName"	      ) abcd.zeroLeptonTriggerEfficiencyName = value;			
      else if(index == "zeroLeptonTriggerEfficiency"		      ) abcd.zeroLeptonTriggerEfficiency = value.Atof();			
      else if(index == "zeroLeptonTriggerEfficiencyError"	      ) abcd.zeroLeptonTriggerEfficiencyError = value.Atof();		
      else if(index == "oneLeptonTriggerEfficiencyName"		      ) abcd.oneLeptonTriggerEfficiencyName = value;			
      else if(index == "oneLeptonTriggerEfficiency"		      ) abcd.oneLeptonTriggerEfficiency = value.Atof();			
      else if(index == "oneLeptonTriggerEfficiencyError"    	      ) abcd.oneLeptonTriggerEfficiencyError = value.Atof();		
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
      else if(index == "lowDeltaPhiNScalingName"                      ) abcd.lowDeltaPhiNScalingName = value;	       
      else if(index == "lowDeltaPhiNMETScaleFactorName"               ) abcd.lowDeltaPhiNMETScaleFactorName = value;	       
      else if(index == "lowDeltaPhiNBTagScaleFactorName"              ) abcd.lowDeltaPhiNBTagScaleFactorName = value;	       
      else if(index == "qcdClosure"		        	      ) abcd.qcdClosure = value.Atof();					
      else if(index == "qcdClosureError"			      ) abcd.qcdClosureError = value.Atof();				
      else if(index == "topWJetsClosure"		              ) abcd.topWJetsClosure = value.Atof();				
      else if(index == "topWJetsClosureError"			      ) abcd.topWJetsClosureError = value.Atof();				
      else if(index == "zeroLeptonDibosonMC"                          ) abcd.zeroLeptonDibosonMC = value.Atof();
      else if(index == "zeroLeptonLowDeltaPhiNDibosonMC"              ) abcd.zeroLeptonLowDeltaPhiNDibosonMC = value.Atof();
      else if(index == "oneLeptonDibosonMC"                           ) abcd.oneLeptonDibosonMC = value.Atof();
      else if(index != "") assert(0);
    }

  binFile.close();

  bins[binName] = abcd;
  observations[binName] = counts;
}

void setupUnderlyingModel(likelihoodOptions options, TString binFilesPath, map<TString,TString>& binFileNames,  
			  TString signalModelFilesPath, map<TString,TString>& binFileNamesInsideSignalMR,  map<TString,TString>& binFileNamesOutsideSignalMR,
			  vector<TString>& binNames, TString& modelFileName, TString& binFilesFileName, allBins& numbers)
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
      if(index == "ZtollOverZtoNuNuRatio"       ) numbers.ZtollOverZtoNuNuRatio = value;
      else if(index == "ZtoeePurity"	        ) numbers.ZtoeePurity = value;	      	
      else if(index == "ZtoeePurityError"       ) numbers.ZtoeePurityError = value;     
      else if(index == "ZtomumuPurity"	        ) numbers.ZtomumuPurity = value;	      
      else if(index == "ZtomumuPurityError"     ) numbers.ZtomumuPurityError = value;   
      else if(index == "ZtoeeEfficiency"        ) numbers.ZtoeeEfficiency = value;	     
      else if(index == "ZtoeeEfficiencyError"   ) numbers.ZtoeeEfficiencyError = value; 
      else if(index == "ZtomumuEfficiency"      ) numbers.ZtomumuEfficiency = value;    
      else if(index == "ZtomumuEfficiencyError" ) numbers.ZtomumuEfficiencyError = value;	
      else if(index == "ZtoeeSystematic"	) numbers.ZtoeeSystematic = value;					
      else if(index == "ZtoeeSystematicError"   ) numbers.ZtoeeSystematicError = value;                            
      else if(index == "ZtomumuSystematic"      ) numbers.ZtomumuSystematic = value;					
      else if(index == "ZtomumuSystematicError" ) numbers.ZtomumuSystematicError = value;                            
      else if(index == "MC" 	                ) numbers.MC = value;	      	   
      else if(index == "MCUncertainty" 	        ) numbers.MCUncertainty = value;	      	   
      else if(index == "dibosonMC" 	        ) numbers.dibosonMC = value;	      	   
      else if(index == "dibosonMCUncertainty"   ) numbers.dibosonMCUncertainty = value;	      	   
      else if(index == "Luminosity"   	        ) numbers.Luminosity = value;
      else if(index == "LuminosityError"        ) numbers.LuminosityError = value;
      else if(index == "metCleaningError"       ) numbers.metCleaningError = value;
      else if(index == "SFqcd_met3"             ) numbers.SFqcd_met3 = value;
      else if(index == "SFqcd_met3_err"         ) numbers.SFqcd_met3_err = value;
      else if(index == "SFqcd_met4"             ) numbers.SFqcd_met4 = value;
      else if(index == "SFqcd_met4_err"         ) numbers.SFqcd_met4_err = value;
      else if(index == "SFqcd_nb3"              ) numbers.SFqcd_nb3 = value;
      else if(index == "SFqcd_nb3_err"          ) numbers.SFqcd_nb3_err = value;
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


void buildLikelihood( TString setupFileName, TString binFilesFileName, TString binFilesPath, TString signalModelFilesPath, TString workspaceName, TString outputFileName, TString binFilesFileNameMR ) 
{
  
  RooWorkspace ws (workspaceName) ;
  ws.autoImportClassCode(true);
  
  vector<TString> binNames;
  map<TString,TString> binFileNames;

  allBinNames names;
  allBins numbers;
  likelihoodOptions options;
  map<TString,abcdBinParameters> bins;
  map<TString,channels> observations;

  //Objets for signal model
  double nGenerated = 0;
  map<TString,yields> signalFractionsOAK;
  map<TString,yields> signalStatisticalErrorOAK;
  map<TString,yields> signalBTagEfficiencyErrorOAK;
  map<TString,yields> signalJesErrorOAK;
  map<TString,TString> binFileNamesInsideSignalMR;
  map<TString,map<TString,mcCount> > insideSignalMR;
  map<TString,TString> binFileNamesOutsideSignalMR;
  map<TString,map<TString,mcCount> > outsideSignalMR;
  
  double oneLeptonTotal(0.);
  double zeroLeptonTopWJetsGuess(0.);

  //Set options 
  options.skipTriggerEfficiency = false;//option no longer works
  options.qcdMethod = "model4";//choices: htDependent, singleScaleWithCorrections
  //options.TopWJetsMethod = "ABCD"; //choices: ABCD, metReweighting
  options.TopWJetsMethod = "metReweighting";
  options.nuisanceOption = "noWidths"; //choices: allWidths, noWidths

  //Read in setupFile and binFilesFile
  setupUnderlyingModel(options, binFilesPath, binFileNames, signalModelFilesPath, binFileNamesInsideSignalMR, binFileNamesOutsideSignalMR, binNames, setupFileName , binFilesFileName , numbers);
  
  //Read in each bin file
  for(map<TString,TString>::iterator thisBin = binFileNames.begin(); thisBin != binFileNames.end() ; thisBin++)
    {
      setupObservations(thisBin->first , thisBin->second , bins, observations);
    }
  
  //Put "global" parameters into workspace
  makeUnderlyingLikelihood(options, ws , names, numbers);

  //Put trigger efficiencies into workspace
  //To get later, do e.g. RooAbsArg* zeroLeptonTriggerEfficiency = ws.arg("zeroLeptonTriggerEfficiency_"+bins["bin12"].zeroLeptonTriggerEfficiencyName); 
  for(vector<TString>::iterator thisBin = binNames.begin(); thisBin != binNames.end() ; thisBin++)
    {
      if( skipBin(*thisBin) ) continue;
      makeTriggerEfficiencies(ws, names, bins[*thisBin] );
    }
  
  //Read in MR signal
  if(options.TopWJetsMethod == "metReweighting")
    {
      for(vector<TString>::iterator thisBin = binNames.begin(); thisBin != binNames.end() ; thisBin++)
	{
	  setupSignalModelMR(*thisBin, binFileNamesInsideSignalMR[*thisBin], insideSignalMR, binFileNamesOutsideSignalMR[*thisBin], outsideSignalMR );
	}
    }
  
  //Read in nominal signal
  setupSignalModelOAK(binNames, signalModelFilesPath, nGenerated, signalFractionsOAK, signalStatisticalErrorOAK, signalBTagEfficiencyErrorOAK, signalJesErrorOAK);

  //Put signal into workspace
  makeSignalModel(options, ws, binNames, names, bins, nGenerated, signalFractionsOAK, signalStatisticalErrorOAK, signalBTagEfficiencyErrorOAK, signalJesErrorOAK, insideSignalMR, outsideSignalMR );

  //Do MET-reweighting method
  if(options.TopWJetsMethod == "metReweighting") buildMRLikelihood(ws, "", binFilesFileNameMR, false, options.nuisanceOption);

  //Do nominal method, QCD, and ZtoNuNu
  for(vector<TString>::iterator thisBin = binNames.begin(); thisBin != binNames.end() ; thisBin++)
    {
      if( skipBin(*thisBin) ) continue;
      makeOneBin(options, ws , *thisBin , names , numbers, observations[*thisBin] , bins[*thisBin] , oneLeptonTotal, zeroLeptonTopWJetsGuess );
    }
  
  //Construct likelihood
  RooArgSet allpdfs = ws.allPdfs();
  cout << endl; cout << endl;
  cout << "allpdfs, size: " << allpdfs.getSize() << endl;
  allpdfs.Print("v");
  
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

  cout << endl; cout << endl;
  cout << "globalObservables, size: " <<   (*ws.set(names.globalObservables)).getSize() << endl;
  (*ws.set(names.globalObservables)).Print("v");

  cout << endl; cout << endl;
  cout << "setting up models" << endl;

  ModelConfig sbModel("S+B_model",&ws);
  sbModel.SetPdf(*ws.pdf("likelihood"));
  sbModel.SetObservables(*ws.set(names.observables));
  sbModel.SetNuisanceParameters(*ws.set(names.nuisances));
  sbModel.SetGlobalObservables(*ws.set(names.globalObservables));
  sbModel.SetParametersOfInterest(*ws.set("poi"));
  sbModel.SetProtoData(*ws.data("data"));

  ModelConfig bModel("B_model",&ws);
  bModel.SetPdf(*ws.pdf("likelihood"));
  bModel.SetObservables(*ws.set(names.observables));
  bModel.SetNuisanceParameters(*ws.set(names.nuisances));
  bModel.SetGlobalObservables(*ws.set(names.globalObservables));
  bModel.SetParametersOfInterest(*ws.set("poi"));
  bModel.SetProtoData(*ws.data("data"));
  ws.var(names.signalCrossSection)->setVal(0.0);
  bModel.SetSnapshot(*ws.set("poi"));

  ws.import (sbModel);
  ws.import (bModel);

  ws.Print() ;
  ws.writeToFile(outputFileName, true) ;

}


void likelihoodBuilder( TString setupFileName, TString binFilesFileName, TString binFilesPath, TString signalModelFilesPath, TString workspaceName, TString outputFileName, TString binFilesFileNameMR ) {
  buildLikelihood( setupFileName, binFilesFileName, binFilesPath, signalModelFilesPath, workspaceName, outputFileName, binFilesFileNameMR );
}
