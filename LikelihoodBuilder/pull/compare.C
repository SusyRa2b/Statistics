#include <iostream>
#include <cassert>

#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPRegexp.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TStyle.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooProduct.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooAddition.h"
#include "RooAbsPdf.h"
#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

using namespace std;

using namespace RooStats;



TString binTranslate(TString name) 
{


  //OAK->LB
  if(name == "M1_H1_1b") return "bin1";
  if(name == "M1_H1_2b") return "bin2";
  if(name == "M1_H1_3b") return "bin3";
  if(name == "M1_H2_1b") return "bin4";
  if(name == "M1_H2_2b") return "bin5";
  if(name == "M1_H2_3b") return "bin6";
  if(name == "M1_H3_1b") return "bin7";
  if(name == "M1_H3_2b") return "bin8";
  if(name == "M1_H3_3b") return "bin9";
  if(name == "M1_H4_1b") return "bin10";
  if(name == "M1_H4_2b") return "bin11";
  if(name == "M1_H4_3b") return "bin12";
  if(name == "M2_H1_1b") return "bin13";
  if(name == "M2_H1_2b") return "bin14";
  if(name == "M2_H1_3b") return "bin15";
  if(name == "M2_H2_1b") return "bin16";
  if(name == "M2_H2_2b") return "bin17";
  if(name == "M2_H2_3b") return "bin18";
  if(name == "M2_H3_1b") return "bin19";
  if(name == "M2_H3_2b") return "bin20";
  if(name == "M2_H3_3b") return "bin21";
  if(name == "M2_H4_1b") return "bin22";
  if(name == "M2_H4_2b") return "bin23";
  if(name == "M2_H4_3b") return "bin24";
  if(name == "M3_H1_1b") return "bin25";
  if(name == "M3_H1_2b") return "bin26";
  if(name == "M3_H1_3b") return "bin27";
  if(name == "M3_H2_1b") return "bin28";
  if(name == "M3_H2_2b") return "bin29";
  if(name == "M3_H2_3b") return "bin30";
  if(name == "M3_H3_1b") return "bin31";
  if(name == "M3_H3_2b") return "bin32";
  if(name == "M3_H3_3b") return "bin33";
  if(name == "M3_H4_1b") return "bin34";
  if(name == "M3_H4_2b") return "bin35";
  if(name == "M3_H4_3b") return "bin36";
  if(name == "M4_H1_1b") return "bin37";
  if(name == "M4_H1_2b") return "bin38";
  if(name == "M4_H1_3b") return "bin39";
  if(name == "M4_H2_1b") return "bin40";
  if(name == "M4_H2_2b") return "bin41";
  if(name == "M4_H2_3b") return "bin42";
  if(name == "M4_H3_1b") return "bin43";
  if(name == "M4_H3_2b") return "bin44";
  if(name == "M4_H3_3b") return "bin45";
  if(name == "M4_H4_1b") return "bin46";
  if(name == "M4_H4_2b") return "bin47";
  if(name == "M4_H4_3b") return "bin48";

  //LB->OAK
  if(name == "bin1") return "M1_H1_1b";
  if(name == "bin2") return "M1_H1_2b";
  if(name == "bin3") return "M1_H1_3b";
  if(name == "bin4") return "M1_H2_1b";
  if(name == "bin5") return "M1_H2_2b";
  if(name == "bin6") return "M1_H2_3b";
  if(name == "bin7") return "M1_H3_1b";
  if(name == "bin8") return "M1_H3_2b";
  if(name == "bin9") return "M1_H3_3b";
  if(name == "bin10") return "M1_H4_1b";
  if(name == "bin11") return "M1_H4_2b";
  if(name == "bin12") return "M1_H4_3b";
  if(name == "bin13") return "M2_H1_1b";
  if(name == "bin14") return "M2_H1_2b";
  if(name == "bin15") return "M2_H1_3b";
  if(name == "bin16") return "M2_H2_1b";
  if(name == "bin17") return "M2_H2_2b";
  if(name == "bin18") return "M2_H2_3b";
  if(name == "bin19") return "M2_H3_1b";
  if(name == "bin20") return "M2_H3_2b";
  if(name == "bin21") return "M2_H3_3b";
  if(name == "bin22") return "M2_H4_1b";
  if(name == "bin23") return "M2_H4_2b";
  if(name == "bin24") return "M2_H4_3b";
  if(name == "bin25") return "M3_H1_1b";
  if(name == "bin26") return "M3_H1_2b";
  if(name == "bin27") return "M3_H1_3b";
  if(name == "bin28") return "M3_H2_1b";
  if(name == "bin29") return "M3_H2_2b";
  if(name == "bin30") return "M3_H2_3b";
  if(name == "bin31") return "M3_H3_1b";
  if(name == "bin32") return "M3_H3_2b";
  if(name == "bin33") return "M3_H3_3b";
  if(name == "bin34") return "M3_H4_1b";
  if(name == "bin35") return "M3_H4_2b";
  if(name == "bin36") return "M3_H4_3b";
  if(name == "bin37") return "M4_H1_1b";
  if(name == "bin38") return "M4_H1_2b";
  if(name == "bin39") return "M4_H1_3b";
  if(name == "bin40") return "M4_H2_1b";
  if(name == "bin41") return "M4_H2_2b";
  if(name == "bin42") return "M4_H2_3b";
  if(name == "bin43") return "M4_H3_1b";
  if(name == "bin44") return "M4_H3_2b";
  if(name == "bin45") return "M4_H3_3b";
  if(name == "bin46") return "M4_H4_1b";
  if(name == "bin47") return "M4_H4_2b";
  if(name == "bin48") return "M4_H4_3b";
  
  else 
    {
      cout << "translation not found!" << endl;
      assert(0);
    }
  return name;

}


TString translate(TString name)
{

  cout << "translating " << name << endl;
  //OAK->LB

  //bin dependent -- eff_sf_M1_H1_1b, eff_sf_ldp_M1_H1_1b, eff_sf_sl_M1_H1_1b, mu_qcd_ldp_M1_H1_1b, mu_ttwj_sl_M1_H1_1b, sf_qcd_M1_H1_1b, sf_ttwj_M1_H1_1b
  if(name.Contains("prim_eff_sf_") || name.Contains("prim_eff_sf_ldp_") || name.Contains("prim_eff_sf_sl_") || 
     name.Contains("mu_qcd_ldp_") || name.Contains("mu_ttwj_sl_") ||
     name.Contains("prim_sf_qcd_") || name.Contains("prim_sf_ttwj_") ) 
    {
      TObjArray *subStrL = TPRegexp("(prim_eff_sf_M|prim_eff_sf_ldp_M|prim_eff_sf_sl_M|mu_qcd_ldp_M|mu_ttwj_sl_M|prim_sf_qcd_M|prim_sf_ttwj_M)(\\S+)").MatchS(name);
      TString prefix = ( ((TObjString *)subStrL->At(1))->GetString() );
      TString bin = "M";
      bin+=( ((TObjString *)subStrL->At(2))->GetString() );
      cout << "prefix: " << prefix << ", bin: " << bin << endl;
      
      TString translatedBin = binTranslate(bin);
      
      if(prefix == "prim_eff_sf_M") return "prim_zeroLeptonSignalError_"+translatedBin;
      if(prefix == "prim_eff_sf_ldp_M") return "prim_zeroLeptonLowDeltaPhiNSignalError_"+translatedBin;
      if(prefix == "prim_eff_sf_sl_M") return "prim_oneLeptonSignalError_"+translatedBin;
      if(prefix == "mu_qcd_ldp_M") return "zeroLeptonLowDeltaPhiN_"+translatedBin+"_QCDYield";
      if(prefix == "mu_ttwj_sl_M") return "oneLepton_"+translatedBin+"_TopWJetsYield";
      if(prefix == "prim_sf_qcd_M") return "prim_zeroLeptonQCDClosure_"+translatedBin;
      if(prefix == "prim_sf_ttwj_M") return "prim_zeroLeptonTopWJetsClosure_"+translatedBin;
    }

  if(name == "prim_JES_sf") return "signalJesErrorCorrelated";

  if(name == "SFqcd_met2" ) return "lowDeltaPhiNMETScaleFactor_M2";
  if(name == "prim_SFqcd_met3" ) return "prim_lowDeltaPhiNMETScaleFactor_M3";
  if(name == "prim_SFqcd_met4" ) return "prim_lowDeltaPhiNMETScaleFactor_M4";
  if(name == "SFqcd_nb2"  ) return "lowDeltaPhiNBTagScaleFactor_2b";
  if(name == "prim_SFqcd_nb3"  ) return "prim_lowDeltaPhiNBTagScaleFactor_3b";
  
  if(name == "acc_Zee_M1" ) return "ZtoeeAcceptance_M1";
  if(name == "acc_Zee_M2" ) return "ZtoeeAcceptance_M2";
  if(name == "acc_Zee_M3" ) return "ZtoeeAcceptance_M3";
  if(name == "acc_Zee_M4" ) return "ZtoeeAcceptance_M4";
  if(name == "acc_Zmm_M1" ) return "ZtomumuAcceptance_M1";
  if(name == "acc_Zmm_M2" ) return "ZtomumuAcceptance_M2";
  if(name == "acc_Zmm_M3" ) return "ZtomumuAcceptance_M3";
  if(name == "acc_Zmm_M4" ) return "ZtomumuAcceptance_M4";

  if(name == "prim_all_gu" ) return "prim_signalGlobalUncertainty";

  if(name == "prim_btageff_sf" ) return "signalBTagEfficiencyCorrelated";
  if(name == "prim_btageff_lf_sf" ) return "signalLTagEfficiencyCorrelated";
  
  if(name == "prim_singletop_xsec" ) return "singleTopShapeSystematic";
  if(name == "prim_wjets_xsec" ) return "wJetsShapeSystematic";
  if(name == "prim_pdfsyst_sf" ) return "signalPdfErrorCorrelated";

  if(name == "eff_Zee" ) return "ZtoeeEfficiency";	
  if(name == "eff_Zmm" ) return "ZtomumuEfficiency";	
  
  if(name == "prim_knn_1b_M1" ) return "prim_zeroLeptonZtoNuNubTagScaling_M1_1b";
  if(name == "prim_knn_1b_M2" ) return "prim_zeroLeptonZtoNuNubTagScaling_M2_1b";
  if(name == "prim_knn_1b_M3" ) return "prim_zeroLeptonZtoNuNubTagScaling_M3_1b";
  if(name == "prim_knn_1b_M4" ) return "prim_zeroLeptonZtoNuNubTagScaling_M4_1b";
  if(name == "prim_knn_2b" ) return "prim_zeroLeptonZtoNuNubTagScaling_2b";
  if(name == "prim_knn_3b" ) return "prim_zeroLeptonZtoNuNubTagScaling_3b";

  if(name == "pur_Zee" ) return "ZtoeePurity";
  if(name == "pur_Zmm" ) return "ZtomumuPurity";

  if(name == "qcd_0lepLDP_ratio_H1" ) return "lowDeltaPhiNScaling_H1";
  if(name == "qcd_0lepLDP_ratio_H2" ) return "lowDeltaPhiNScaling_H2";
  if(name == "qcd_0lepLDP_ratio_H3" ) return "lowDeltaPhiNScaling_H3";
  if(name == "qcd_0lepLDP_ratio_H4" ) return "lowDeltaPhiNScaling_H4";
  
  if(name == "prim_rar_vv_sf" ) return "prim_dibosonMCUncertainty";
  
  if(name == "prim_sf_ll" ) return "ZtollSystematic";
  
  if(name == "prim_sf_mc" ) return "prim_MCUncertainty";
  
  if(name == "trigeff_M1_H1" ) return "zeroLeptonTriggerEfficiency_M1_H1";
  if(name == "trigeff_M1_H2" ) return "zeroLeptonTriggerEfficiency_M1_H2";
  if(name == "trigeff_M1_H3" ) return "zeroLeptonTriggerEfficiency_M1_H3";
  if(name == "trigeff_M1_H4" ) return "zeroLeptonTriggerEfficiency_M1_H4";
  if(name == "trigeff_M2_H1" ) return "zeroLeptonTriggerEfficiency_M2_H1";
  if(name == "trigeff_M2_H2" ) return "zeroLeptonTriggerEfficiency_M2_H2";
  if(name == "trigeff_M2_H3" ) return "zeroLeptonTriggerEfficiency_M2_H3";
  if(name == "trigeff_M2_H4" ) return "zeroLeptonTriggerEfficiency_M2_H4";
  if(name == "trigeff_M3_H1" ) return "zeroLeptonTriggerEfficiency_M3_H1";
  if(name == "trigeff_M3_H2" ) return "zeroLeptonTriggerEfficiency_M3_H2";
  if(name == "trigeff_M3_H3" ) return "zeroLeptonTriggerEfficiency_M3_H3";
  if(name == "trigeff_M3_H4" ) return "zeroLeptonTriggerEfficiency_M3_H4";
  if(name == "trigeff_M4_H2" ) return "zeroLeptonTriggerEfficiency_M4_H2";
  if(name == "trigeff_M4_H3" ) return "zeroLeptonTriggerEfficiency_M4_H3";
  if(name == "trigeff_M4_H4" ) return "zeroLeptonTriggerEfficiency_M4_H4";

  if(name == "trigeff_sl_M1_H1" ) return "oneLeptonTriggerEfficiency_M1_H1";
  if(name == "trigeff_sl_M1_H2" ) return "oneLeptonTriggerEfficiency_M1_H2";
  if(name == "trigeff_sl_M1_H3" ) return "oneLeptonTriggerEfficiency_M1_H3";
  if(name == "trigeff_sl_M1_H4" ) return "oneLeptonTriggerEfficiency_M1_H4";
  if(name == "trigeff_sl_M2_H1" ) return "oneLeptonTriggerEfficiency_M2_H1";
  if(name == "trigeff_sl_M2_H2" ) return "oneLeptonTriggerEfficiency_M2_H2";
  if(name == "trigeff_sl_M2_H3" ) return "oneLeptonTriggerEfficiency_M2_H3";
  if(name == "trigeff_sl_M2_H4" ) return "oneLeptonTriggerEfficiency_M2_H4";
  if(name == "trigeff_sl_M3_H1" ) return "oneLeptonTriggerEfficiency_M3_H1";
  if(name == "trigeff_sl_M3_H2" ) return "oneLeptonTriggerEfficiency_M3_H2";
  if(name == "trigeff_sl_M3_H3" ) return "oneLeptonTriggerEfficiency_M3_H3";
  if(name == "trigeff_sl_M3_H4" ) return "oneLeptonTriggerEfficiency_M3_H4";
  if(name == "trigeff_sl_M4_H2" ) return "oneLeptonTriggerEfficiency_M4_H2";
  if(name == "trigeff_sl_M4_H3" ) return "oneLeptonTriggerEfficiency_M4_H3";
  if(name == "trigeff_sl_M4_H4" ) return "oneLeptonTriggerEfficiency_M4_H4";

  if(name == "ttwj_0lep1lep_ratio" ) return "singleLeptonScaling";
  
  //mu_susy_all0lep
  
  /*
    if(name == "mu_znn_M1_H1_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M1_H1";
    if(name == "mu_znn_M1_H2_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M1_H2";
    if(name == "mu_znn_M1_H3_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M1_H3";
    if(name == "mu_znn_M1_H4_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M1_H4";
    if(name == "mu_znn_M2_H1_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M2_H1";
    if(name == "mu_znn_M2_H2_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M2_H2";
    if(name == "mu_znn_M2_H3_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M2_H3";
    if(name == "mu_znn_M2_H4_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M2_H4";
    if(name == "mu_znn_M3_H1_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M3_H1";
    if(name == "mu_znn_M3_H2_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M3_H2";
    if(name == "mu_znn_M3_H3_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M3_H3";
    if(name == "mu_znn_M3_H4_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M3_H4";
    if(name == "mu_znn_M4_H2_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M4_H2";
    if(name == "mu_znn_M4_H3_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M4_H3";
    if(name == "mu_znn_M4_H4_1b" ) return "ZtoNuNu_VeryLooseBtagYield_M4_H4";
  */
  

  
  //LB->OAK

  
  else 
    {
      cout << "translation not found!" << endl;
      assert(0);
    }
  return name;
}





TString translatePDF(TString name)
{

  //cout << "translating pdf " << name << endl;
  //OAK->LB

  //bin dependent -- eff_sf_M1_H1_1b, eff_sf_ldp_M1_H1_1b, eff_sf_sl_M1_H1_1b, mu_qcd_ldp_M1_H1_1b, mu_ttwj_sl_M1_H1_1b, sf_qcd_M1_H1_1b, sf_ttwj_M1_H1_1b
  if(name.Contains("pdf_prim_eff_sf_") || name.Contains("pdf_prim_eff_sf_ldp_") || name.Contains("pdf_prim_eff_sf_sl_") || 
     name.Contains("pdf_N_0lep_") || name.Contains("pdf_N_1lep_") || name.Contains("pdf_N_ldp_") ||
     name.Contains("pdf_prim_sf_qcd_") || name.Contains("pdf_prim_sf_ttwj_") ) 
    {
      TObjArray *subStrL = TPRegexp("(pdf_prim_eff_sf_M|pdf_prim_eff_sf_ldp_M|pdf_prim_eff_sf_sl_M|pdf_N_0lep_M|pdf_N_1lep_M|pdf_N_ldp_M|pdf_prim_sf_qcd_M|pdf_prim_sf_ttwj_M)(\\S+)").MatchS(name);
      TString prefix = ( ((TObjString *)subStrL->At(1))->GetString() );
      TString bin = "M";
      bin+=( ((TObjString *)subStrL->At(2))->GetString() );
      //cout << "prefix: " << prefix << ", bin: " << bin << endl;
      
      TString translatedBin = binTranslate(bin);
      
      if(prefix == "pdf_prim_eff_sf_M") return "pdf_prim_zeroLeptonSignalError_"+translatedBin;
      if(prefix == "pdf_prim_eff_sf_ldp_M") return "pdf_prim_zeroLeptonLowDeltaPhiNSignalError_"+translatedBin;
      if(prefix == "pdf_prim_eff_sf_sl_M") return "pdf_prim_oneLeptonSignalError_"+translatedBin;
      if(prefix == "pdf_N_0lep_M") return "zeroLepton_"+translatedBin+"_Constraint";
      if(prefix == "pdf_N_1lep_M") return "oneLepton_"+translatedBin+"_Constraint";
      if(prefix == "pdf_N_ldp_M") return "zeroLeptonLowDeltaPhiN_"+translatedBin+"_Constraint";
      if(prefix == "pdf_prim_sf_qcd_M") return "pdf_prim_zeroLeptonQCDClosure_"+translatedBin;
      if(prefix == "pdf_prim_sf_ttwj_M") return "pdf_prim_zeroLeptonTopWJetsClosure_"+translatedBin;
    }

  if(name == "pdf_prim_JES_sf") return "pdf_signalJesErrorCorrelated";

  if(name == "pdf_prim_SFqcd_met3" ) return "pdf_prim_lowDeltaPhiNMETScaleFactor_M3";
  if(name == "pdf_prim_SFqcd_met4" ) return "pdf_prim_lowDeltaPhiNMETScaleFactor_M4";
  if(name == "pdf_prim_SFqcd_nb3"  ) return "pdf_prim_lowDeltaPhiNBTagScaleFactor_3b";
  
  if(name == "betapdf_acc_Zee_M1" ) return "ZtoeeAcceptance_M1_Constraint";
  if(name == "betapdf_acc_Zee_M2" ) return "ZtoeeAcceptance_M2_Constraint";
  if(name == "betapdf_acc_Zee_M3" ) return "ZtoeeAcceptance_M3_Constraint";
  if(name == "betapdf_acc_Zee_M4" ) return "ZtoeeAcceptance_M4_Constraint";
  if(name == "betapdf_acc_Zmm_M1" ) return "ZtomumuAcceptance_M1_Constraint";
  if(name == "betapdf_acc_Zmm_M2" ) return "ZtomumuAcceptance_M2_Constraint";
  if(name == "betapdf_acc_Zmm_M3" ) return "ZtomumuAcceptance_M3_Constraint";
  if(name == "betapdf_acc_Zmm_M4" ) return "ZtomumuAcceptance_M4_Constraint";

  if(name == "pdf_prim_all_gu" ) return "pdf_prim_signalGlobalUncertainty";

  if(name == "pdf_prim_btageff_sf" ) return "pdf_signalBTagEfficiencyCorrelated";
  if(name == "pdf_prim_btageff_lf_sf" ) return "pdf_signalLTagEfficiencyCorrelated";
  
  if(name == "pdf_prim_singletop_xsec" ) return "pdf_singleTopShapeSystematic";
  if(name == "pdf_prim_wjets_xsec" ) return "pdf_wJetsShapeSystematic";

  if(name == "pdf_prim_pdfsyst_sf" ) return "pdf_signalPdfErrorCorrelated";

  if(name == "betapdf_eff_Zee" ) return "ZtoeeEfficiency_Constraint";	
  if(name == "betapdf_eff_Zmm" ) return "ZtomumuEfficiency_Constraint";	
  
  if(name == "pdf_prim_knn_1b_M1" ) return "pdf_prim_zeroLeptonZtoNuNubTagScaling_M1_1b";
  if(name == "pdf_prim_knn_1b_M2" ) return "pdf_prim_zeroLeptonZtoNuNubTagScaling_M2_1b";
  if(name == "pdf_prim_knn_1b_M3" ) return "pdf_prim_zeroLeptonZtoNuNubTagScaling_M3_1b";
  if(name == "pdf_prim_knn_1b_M4" ) return "pdf_prim_zeroLeptonZtoNuNubTagScaling_M4_1b";
  if(name == "pdf_prim_knn_2b" ) return "pdf_prim_zeroLeptonZtoNuNubTagScaling_2b";
  if(name == "pdf_prim_knn_3b" ) return "pdf_prim_zeroLeptonZtoNuNubTagScaling_3b";

  if(name == "betapdf_pur_Zee" ) return "ZtoeePurity_Constraint";
  if(name == "betapdf_pur_Zmm" ) return "ZtomumuPurity_Constraint";

  if(name == "pdf_prim_rar_vv_sf" ) return "pdf_prim_dibosonMCUncertainty";
  
  if(name == "pdf_prim_sf_ll" ) return "pdf_ZtollSystematic";
  
  if(name == "pdf_prim_sf_mc" ) return "pdf_prim_MCUncertainty";
  
  if(name == "betapdf_trigeff_M1_H1" ) return "zeroLeptonTriggerEfficiency_M1_H1_Constraint";
  if(name == "betapdf_trigeff_M1_H2" ) return "zeroLeptonTriggerEfficiency_M1_H2_Constraint";
  if(name == "betapdf_trigeff_M1_H3" ) return "zeroLeptonTriggerEfficiency_M1_H3_Constraint";
  if(name == "betapdf_trigeff_M1_H4" ) return "zeroLeptonTriggerEfficiency_M1_H4_Constraint";
  if(name == "betapdf_trigeff_M2_H1" ) return "zeroLeptonTriggerEfficiency_M2_H1_Constraint";
  if(name == "betapdf_trigeff_M2_H2" ) return "zeroLeptonTriggerEfficiency_M2_H2_Constraint";
  if(name == "betapdf_trigeff_M2_H3" ) return "zeroLeptonTriggerEfficiency_M2_H3_Constraint";
  if(name == "betapdf_trigeff_M2_H4" ) return "zeroLeptonTriggerEfficiency_M2_H4_Constraint";
  if(name == "betapdf_trigeff_M3_H1" ) return "zeroLeptonTriggerEfficiency_M3_H1_Constraint";
  if(name == "betapdf_trigeff_M3_H2" ) return "zeroLeptonTriggerEfficiency_M3_H2_Constraint";
  if(name == "betapdf_trigeff_M3_H3" ) return "zeroLeptonTriggerEfficiency_M3_H3_Constraint";
  if(name == "betapdf_trigeff_M3_H4" ) return "zeroLeptonTriggerEfficiency_M3_H4_Constraint";
  if(name == "betapdf_trigeff_M4_H2" ) return "zeroLeptonTriggerEfficiency_M4_H2_Constraint";
  if(name == "betapdf_trigeff_M4_H3" ) return "zeroLeptonTriggerEfficiency_M4_H3_Constraint";
  if(name == "betapdf_trigeff_M4_H4" ) return "zeroLeptonTriggerEfficiency_M4_H4_Constraint";

  if(name == "betapdf_trigeff_sl_M1_H1" ) return "oneLeptonTriggerEfficiency_M1_H1_Constraint";
  if(name == "betapdf_trigeff_sl_M1_H2" ) return "oneLeptonTriggerEfficiency_M1_H2_Constraint";
  if(name == "betapdf_trigeff_sl_M1_H3" ) return "oneLeptonTriggerEfficiency_M1_H3_Constraint";
  if(name == "betapdf_trigeff_sl_M1_H4" ) return "oneLeptonTriggerEfficiency_M1_H4_Constraint";
  if(name == "betapdf_trigeff_sl_M2_H1" ) return "oneLeptonTriggerEfficiency_M2_H1_Constraint";
  if(name == "betapdf_trigeff_sl_M2_H2" ) return "oneLeptonTriggerEfficiency_M2_H2_Constraint";
  if(name == "betapdf_trigeff_sl_M2_H3" ) return "oneLeptonTriggerEfficiency_M2_H3_Constraint";
  if(name == "betapdf_trigeff_sl_M2_H4" ) return "oneLeptonTriggerEfficiency_M2_H4_Constraint";
  if(name == "betapdf_trigeff_sl_M3_H1" ) return "oneLeptonTriggerEfficiency_M3_H1_Constraint";
  if(name == "betapdf_trigeff_sl_M3_H2" ) return "oneLeptonTriggerEfficiency_M3_H2_Constraint";
  if(name == "betapdf_trigeff_sl_M3_H3" ) return "oneLeptonTriggerEfficiency_M3_H3_Constraint";
  if(name == "betapdf_trigeff_sl_M3_H4" ) return "oneLeptonTriggerEfficiency_M3_H4_Constraint";
  if(name == "betapdf_trigeff_sl_M4_H2" ) return "oneLeptonTriggerEfficiency_M4_H2_Constraint";
  if(name == "betapdf_trigeff_sl_M4_H3" ) return "oneLeptonTriggerEfficiency_M4_H3_Constraint";
  if(name == "betapdf_trigeff_sl_M4_H4" ) return "oneLeptonTriggerEfficiency_M4_H4_Constraint";

  
  if(name == "pdf_N_Zee_M1_H1" ) return "diElectron_M1_H1_Constraint";
  if(name == "pdf_N_Zee_M1_H2" ) return "diElectron_M1_H2_Constraint";
  if(name == "pdf_N_Zee_M1_H3" ) return "diElectron_M1_H3_Constraint";
  if(name == "pdf_N_Zee_M1_H4" ) return "diElectron_M1_H4_Constraint";
  if(name == "pdf_N_Zee_M2_H1" ) return "diElectron_M2_H1_Constraint";
  if(name == "pdf_N_Zee_M2_H2" ) return "diElectron_M2_H2_Constraint";
  if(name == "pdf_N_Zee_M2_H3" ) return "diElectron_M2_H3_Constraint";
  if(name == "pdf_N_Zee_M2_H4" ) return "diElectron_M2_H4_Constraint";
  if(name == "pdf_N_Zee_M3_H1" ) return "diElectron_M3_H1_Constraint";
  if(name == "pdf_N_Zee_M3_H2" ) return "diElectron_M3_H2_Constraint";
  if(name == "pdf_N_Zee_M3_H3" ) return "diElectron_M3_H3_Constraint";
  if(name == "pdf_N_Zee_M3_H4" ) return "diElectron_M3_H4_Constraint";
  if(name == "pdf_N_Zee_M4_H2" ) return "diElectron_M4_H2_Constraint";
  if(name == "pdf_N_Zee_M4_H3" ) return "diElectron_M4_H3_Constraint";
  if(name == "pdf_N_Zee_M4_H4" ) return "diElectron_M4_H4_Constraint";

  if(name == "pdf_N_Zmm_M1_H1" ) return "diMuon_M1_H1_Constraint";
  if(name == "pdf_N_Zmm_M1_H2" ) return "diMuon_M1_H2_Constraint";
  if(name == "pdf_N_Zmm_M1_H3" ) return "diMuon_M1_H3_Constraint";
  if(name == "pdf_N_Zmm_M1_H4" ) return "diMuon_M1_H4_Constraint";
  if(name == "pdf_N_Zmm_M2_H1" ) return "diMuon_M2_H1_Constraint";
  if(name == "pdf_N_Zmm_M2_H2" ) return "diMuon_M2_H2_Constraint";
  if(name == "pdf_N_Zmm_M2_H3" ) return "diMuon_M2_H3_Constraint";
  if(name == "pdf_N_Zmm_M2_H4" ) return "diMuon_M2_H4_Constraint";
  if(name == "pdf_N_Zmm_M3_H1" ) return "diMuon_M3_H1_Constraint";
  if(name == "pdf_N_Zmm_M3_H2" ) return "diMuon_M3_H2_Constraint";
  if(name == "pdf_N_Zmm_M3_H3" ) return "diMuon_M3_H3_Constraint";
  if(name == "pdf_N_Zmm_M3_H4" ) return "diMuon_M3_H4_Constraint";
  if(name == "pdf_N_Zmm_M4_H2" ) return "diMuon_M4_H2_Constraint";
  if(name == "pdf_N_Zmm_M4_H3" ) return "diMuon_M4_H3_Constraint";
  if(name == "pdf_N_Zmm_M4_H4" ) return "diMuon_M4_H4_Constraint";
    

  
  //LB->OAK

  
  else 
    {
      cout << "translation not found!" << endl;
      assert(0);
    }
  return name;
}






//void makeEqual() { }

double percentDifference(double a, double b) 
{
  return 100.0*(a-b)/(0.5*(a+b));
}


void printDiff_RRV_RRV(RooWorkspace* ws1, TString name1, RooWorkspace* ws2, TString name2)
{
  cout.precision(10);
  double v1 = (ws1->var(name1))->getVal();
  double v2 = (ws2->var(name2))->getVal();
  cout << name1 << " " << name2 << " " << v1 << " " << v2 << " " << v1-v2 << " " << percentDifference(v1, v2) << endl;
}

void printDiff_RAR_RRV(RooWorkspace* ws1, TString name1, RooWorkspace* ws2, TString name2)
{
  double v1 = (ws1->function(name1))->getVal();
  double v2 = (ws2->var(name2))->getVal();
  cout << name1 << " " << name2 << " " << v1 << " " << v2 << " " << v1-v2 << " " << percentDifference(v1, v2) << endl;
}

void printDiff_RAR_RAR(RooWorkspace* ws1, TString name1, RooWorkspace* ws2, TString name2)
{
  double v1 = (ws1->function(name1))->getVal();
  double v2 = (ws2->function(name2))->getVal();
  cout << name1 << " " << name2 << " " << v1 << " " << v2 << " " << v1-v2 << " " << percentDifference(v1, v2) << endl;
}

void printDiff_PDF_PDF(RooWorkspace* ws1, TString name1, RooWorkspace* ws2, TString name2)
{

  cout.precision(20);
  //  cout << "Printing " << name1 << endl;
  double v1 = (ws1->pdf(name1))->getVal();
  //  cout << "Printing " << name2 << endl;
  double v2 = (ws2->pdf(name2))->getVal();
  cout << name1 << " " << name2 << " " << v1 << " " << v2 << " " << v1-v2 << " " << percentDifference(v1, v2) << endl;
}


void printZtoNuNu(RooWorkspace* ws1, TString name1, RooWorkspace* ws2, TString bin2, TString scale2)
{
  double znnOAK = (ws1->function(name1))->getVal();
  RooAbsReal* ZtoNuNuVeryLoose = ws2->function("ZtoNuNu_VeryLooseBtagYield_"+bin2);
  RooAbsReal* knn = ws2->function("zeroLeptonZtoNuNubTagScaling_"+scale2);
  RooProduct ZtoNuNu("ZtoNuNu", "ZtoNuNu", RooArgSet(*ZtoNuNuVeryLoose, *knn));
  cout << name1 << " " << "xxx" << " " << znnOAK << " " << ZtoNuNu.getVal() << " " << znnOAK - ZtoNuNu.getVal() << " " << percentDifference(znnOAK, ZtoNuNu.getVal()) << endl;

}

double getSignalLB(RooWorkspace* wsLB)
{
  
  RooArgSet* signalYieldSet = new RooArgSet("signalYields");
  for(unsigned int i =1; i<=48; i++) {
    TString binName = "bin";
    binName+=i;
    if(binName=="bin37" || binName=="bin38" || binName=="bin39") continue;
    signalYieldSet->add( *(wsLB->arg("zeroLepton_"+binName+"_SignalYield")) );
    //cout << (wsLB->function("zeroLepton_"+binName+"_SignalYield"))->getVal() << endl;
    //cout << binName << " signal yield pointer = " << (wsLB->arg("zeroLepton_"+binName+"_SignalYield"))  << endl;
  }
  RooAddition* signalYield = new RooAddition("signalYield", "signalYield", *signalYieldSet);

  return signalYield->getVal();
}

void compareFit()
{
  
  //Open workspaces
  TFile fOAK("pull/OAK.root", "READ");
  RooWorkspace* wsOAK = (RooWorkspace*)fOAK.Get("ws");

  //TFile fLB("likelihood_toy100_Input-met4-ht4-v15-newqcdsyst-model4-exp0lep-ttwjave-wtrig-toy0000.root", "READ");
  TFile fLB("pull/LB.root", "READ");
  RooWorkspace *wsLB = (RooWorkspace*)fLB.Get("workspace");
  

  //First the things that don't need special treatment
  //eff_sf_M1_H1_1b, eff_sf_ldp_M1_H1_1b, eff_sf_sl_M1_H1_1b, mu_qcd_ldp_M1_H1_1b, mu_ttwj_sl_M1_H1_1b, sf_qcd_M1_H1_1b, sf_ttwj_M1_H1_1b
  vector<TString> names;
  names.push_back("eff_sf_");
  names.push_back("eff_sf_ldp_");
  names.push_back("eff_sf_sl_");
  names.push_back("mu_qcd_ldp_");
  names.push_back("mu_ttwj_sl_");
  names.push_back("sf_qcd_");
  names.push_back("sf_ttwj_");
  for(vector<TString>::iterator it = names.begin(); it != names.end(); ++it)
    {
      TString name = *it;
      //cout << name << endl;

      for(int m = 1; m<=4; m++)
	{
	  for(int h = 1; h<=4; h++)
	    {
	      if(m==4 && h==1) continue;
	      for(int b = 1; b<=3; b++)
		{
		  TString bin = "M"; bin+=m; bin+="_H"; bin+=h; bin+="_"; bin+=b; bin+="b";
		  
		  TString varName = name;
		  varName += bin;
		  //cout << varName << " " << translate(varName) << endl;
		  
		  if(name.Contains("eff_sf_") || name.Contains("sf_qcd_") || name.Contains("sf_ttwj_") )
		    {
		      printDiff_RAR_RRV(wsOAK, varName, wsLB, translate(varName));
		    }

		  if(name.Contains("mu_"))
		    {
		     printDiff_RRV_RRV(wsOAK, varName, wsLB, translate(varName));
		    }

		}//b
	    }//h
	}//m

    }//vector of names
  
  
  //Others
  printDiff_RAR_RRV(wsOAK, "JES_sf", wsLB, translate("JES_sf"));
  
  printDiff_RAR_RRV(wsOAK, "SFqcd_met2", wsLB, translate("SFqcd_met2"));
  printDiff_RAR_RRV(wsOAK, "SFqcd_met3", wsLB, translate("SFqcd_met3"));
  printDiff_RAR_RRV(wsOAK, "SFqcd_met4", wsLB, translate("SFqcd_met4"));
  printDiff_RAR_RRV(wsOAK, "SFqcd_nb2", wsLB, translate("SFqcd_nb2"));
  printDiff_RAR_RRV(wsOAK, "SFqcd_nb3", wsLB, translate("SFqcd_nb3"));

  printDiff_RAR_RRV(wsOAK, "acc_Zee_M1", wsLB, translate("acc_Zee_M1"));
  printDiff_RAR_RRV(wsOAK, "acc_Zee_M2", wsLB, translate("acc_Zee_M2"));
  printDiff_RAR_RRV(wsOAK, "acc_Zee_M3", wsLB, translate("acc_Zee_M3"));
  printDiff_RAR_RRV(wsOAK, "acc_Zee_M4", wsLB, translate("acc_Zee_M4"));
  printDiff_RAR_RRV(wsOAK, "acc_Zmm_M1", wsLB, translate("acc_Zmm_M1"));
  printDiff_RAR_RRV(wsOAK, "acc_Zmm_M2", wsLB, translate("acc_Zmm_M2"));
  printDiff_RAR_RRV(wsOAK, "acc_Zmm_M3", wsLB, translate("acc_Zmm_M3"));
  printDiff_RAR_RRV(wsOAK, "acc_Zmm_M4", wsLB, translate("acc_Zmm_M4"));

  printDiff_RAR_RRV(wsOAK, "all_gu", wsLB, translate("all_gu"));
  
  printDiff_RAR_RRV(wsOAK, "btageff_sf", wsLB, translate("btageff_sf"));

  printDiff_RAR_RRV(wsOAK, "eff_Zee", wsLB, translate("eff_Zee"));
  printDiff_RAR_RRV(wsOAK, "eff_Zmm", wsLB, translate("eff_Zmm"));
  
  printDiff_RAR_RRV(wsOAK, "knn_1b_M1", wsLB, translate("knn_1b_M1"));
  printDiff_RAR_RRV(wsOAK, "knn_1b_M2", wsLB, translate("knn_1b_M2"));
  printDiff_RAR_RRV(wsOAK, "knn_1b_M3", wsLB, translate("knn_1b_M3"));
  printDiff_RAR_RRV(wsOAK, "knn_1b_M4", wsLB, translate("knn_1b_M4"));
  printDiff_RAR_RRV(wsOAK, "knn_2b", wsLB, translate("knn_2b"));
  printDiff_RAR_RRV(wsOAK, "knn_3b", wsLB, translate("knn_3b"));

  printDiff_RAR_RRV(wsOAK, "pur_Zee", wsLB, translate("pur_Zee"));
  printDiff_RAR_RRV(wsOAK, "pur_Zmm", wsLB, translate("pur_Zmm"));

  printDiff_RRV_RRV(wsOAK, "qcd_0lepLDP_ratio_H1", wsLB, translate("qcd_0lepLDP_ratio_H1"));
  printDiff_RRV_RRV(wsOAK, "qcd_0lepLDP_ratio_H2", wsLB, translate("qcd_0lepLDP_ratio_H2"));
  printDiff_RRV_RRV(wsOAK, "qcd_0lepLDP_ratio_H3", wsLB, translate("qcd_0lepLDP_ratio_H3"));
  printDiff_RRV_RRV(wsOAK, "qcd_0lepLDP_ratio_H4", wsLB, translate("qcd_0lepLDP_ratio_H4"));

  printDiff_RAR_RRV(wsOAK, "rar_vv_sf", wsLB, translate("rar_vv_sf"));

  printDiff_RAR_RRV(wsOAK, "sf_ll", wsLB, translate("sf_ll"));
  
  printDiff_RAR_RRV(wsOAK, "sf_mc", wsLB, translate("sf_mc"));
  
  printDiff_RAR_RRV(wsOAK, "trigeff_M1_H1", wsLB, translate("trigeff_M1_H1"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M1_H2", wsLB, translate("trigeff_M1_H2"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M1_H3", wsLB, translate("trigeff_M1_H3"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M1_H4", wsLB, translate("trigeff_M1_H4"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M2_H1", wsLB, translate("trigeff_M2_H1"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M2_H2", wsLB, translate("trigeff_M2_H2"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M2_H3", wsLB, translate("trigeff_M2_H3"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M2_H4", wsLB, translate("trigeff_M2_H4"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M3_H1", wsLB, translate("trigeff_M3_H1"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M3_H2", wsLB, translate("trigeff_M3_H2"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M3_H3", wsLB, translate("trigeff_M3_H3"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M3_H4", wsLB, translate("trigeff_M3_H4"));
  //printDiff_RAR_RRV(wsOAK, "trigeff_M4_H1", wsLB, translate("trigeff_M4_H1"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M4_H2", wsLB, translate("trigeff_M4_H2"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M4_H3", wsLB, translate("trigeff_M4_H3"));
  printDiff_RAR_RRV(wsOAK, "trigeff_M4_H4", wsLB, translate("trigeff_M4_H4"));

  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M1_H1", wsLB, translate("trigeff_sl_M1_H1"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M1_H2", wsLB, translate("trigeff_sl_M1_H2"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M1_H3", wsLB, translate("trigeff_sl_M1_H3"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M1_H4", wsLB, translate("trigeff_sl_M1_H4"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M2_H1", wsLB, translate("trigeff_sl_M2_H1"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M2_H2", wsLB, translate("trigeff_sl_M2_H2"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M2_H3", wsLB, translate("trigeff_sl_M2_H3"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M2_H4", wsLB, translate("trigeff_sl_M2_H4"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M3_H1", wsLB, translate("trigeff_sl_M3_H1"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M3_H2", wsLB, translate("trigeff_sl_M3_H2"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M3_H3", wsLB, translate("trigeff_sl_M3_H3"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M3_H4", wsLB, translate("trigeff_sl_M3_H4"));
  //printDiff_RAR_RRV(wsOAK, "trigeff_sl_M4_H1", wsLB, translate("trigeff_sl_M4_H1"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M4_H2", wsLB, translate("trigeff_sl_M4_H2"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M4_H3", wsLB, translate("trigeff_sl_M4_H3"));
  printDiff_RAR_RRV(wsOAK, "trigeff_sl_M4_H4", wsLB, translate("trigeff_sl_M4_H4"));
  
  
  printZtoNuNu(wsOAK, "mu_znn_M1_H1_1b", wsLB, "M1_H1", "M1_1b");
  printZtoNuNu(wsOAK, "mu_znn_M1_H2_1b", wsLB, "M1_H2", "M1_1b");
  printZtoNuNu(wsOAK, "mu_znn_M1_H3_1b", wsLB, "M1_H3", "M1_1b");
  printZtoNuNu(wsOAK, "mu_znn_M1_H4_1b", wsLB, "M1_H4", "M1_1b");
  printZtoNuNu(wsOAK, "mu_znn_M2_H1_1b", wsLB, "M2_H1", "M2_1b");
  printZtoNuNu(wsOAK, "mu_znn_M2_H2_1b", wsLB, "M2_H2", "M2_1b");
  printZtoNuNu(wsOAK, "mu_znn_M2_H3_1b", wsLB, "M2_H3", "M2_1b");
  printZtoNuNu(wsOAK, "mu_znn_M2_H4_1b", wsLB, "M2_H4", "M2_1b");
  printZtoNuNu(wsOAK, "mu_znn_M3_H1_1b", wsLB, "M3_H1", "M3_1b");
  printZtoNuNu(wsOAK, "mu_znn_M3_H2_1b", wsLB, "M3_H2", "M3_1b");
  printZtoNuNu(wsOAK, "mu_znn_M3_H3_1b", wsLB, "M3_H3", "M3_1b");
  printZtoNuNu(wsOAK, "mu_znn_M3_H4_1b", wsLB, "M3_H4", "M3_1b");
  printZtoNuNu(wsOAK, "mu_znn_M4_H2_1b", wsLB, "M4_H2", "M4_1b");
  printZtoNuNu(wsOAK, "mu_znn_M4_H3_1b", wsLB, "M4_H3", "M4_1b");
  printZtoNuNu(wsOAK, "mu_znn_M4_H4_1b", wsLB, "M4_H4", "M4_1b");

  printDiff_RRV_RRV(wsOAK, "ttwj_0lep1lep_ratio", wsLB, translate("ttwj_0lep1lep_ratio"));

  double susyOAK = (wsOAK->var("mu_susy_all0lep"))->getVal();
  double susyLB = getSignalLB(wsLB);
  cout << "mu_susy_all0lep" << " " << "xxx" << " " << susyOAK << " " << susyLB << " " << susyOAK-susyLB << " " << percentDifference(susyOAK, susyLB) << endl;


  fOAK.Close();
  fLB.Close();

  return;
}



void comparePDFs(RooWorkspace* wsLB, RooWorkspace* wsOAK)
{

  /*
   //Open workspaces
  TFile fOAK("ws3b_Input-met4-ht4-v15-newqcdsyst-model4-exp0lep-ttwjave-wtrig-toy0000.root", "READ");
  RooWorkspace* wsOAK = (RooWorkspace*)fOAK.Get("ws");

  TFile fLB("likelihood_toy100_Input-met4-ht4-v15-newqcdsyst-model4-exp0lep-ttwjave-wtrig-toy0000.root", "READ");
  RooWorkspace *wsLB = (RooWorkspace*)fLB.Get("workspace");
  */

  //First the things that don't need special treatment
  //eff_sf_M1_H1_1b, eff_sf_ldp_M1_H1_1b, eff_sf_sl_M1_H1_1b, mu_qcd_ldp_M1_H1_1b, mu_ttwj_sl_M1_H1_1b, sf_qcd_M1_H1_1b, sf_ttwj_M1_H1_1b
  vector<TString> names;
  names.push_back("pdf_prim_eff_sf_");
  names.push_back("pdf_prim_eff_sf_ldp_");
  names.push_back("pdf_prim_eff_sf_sl_");
  names.push_back("pdf_N_0lep_");
  names.push_back("pdf_N_1lep_");
  names.push_back("pdf_N_ldp_");
  names.push_back("pdf_prim_sf_qcd_");
  names.push_back("pdf_prim_sf_ttwj_");
  for(vector<TString>::iterator it = names.begin(); it != names.end(); ++it)
    {
      TString name = *it;
      //cout << name << endl;

      for(int m = 1; m<=4; m++)
	{
	  for(int h = 1; h<=4; h++)
	    {
	      if(m==4 && h==1) continue;
	      for(int b = 1; b<=3; b++)
		{
		  TString bin = "M"; bin+=m; bin+="_H"; bin+=h; bin+="_"; bin+=b; bin+="b";
		  
		  TString varName = name;
		  varName += bin;
		  
		  printDiff_PDF_PDF(wsOAK, varName, wsLB, translatePDF(varName));
		  
		}//b
	    }//h
	}//m

    }//vector of names
  
  
  printDiff_PDF_PDF(wsOAK, "pdf_prim_JES_sf", wsLB, translatePDF("pdf_prim_JES_sf"));
  
  printDiff_PDF_PDF(wsOAK, "pdf_prim_SFqcd_met3", wsLB, translatePDF("pdf_prim_SFqcd_met3"));
  printDiff_PDF_PDF(wsOAK, "pdf_prim_SFqcd_met4", wsLB, translatePDF("pdf_prim_SFqcd_met4"));
  printDiff_PDF_PDF(wsOAK, "pdf_prim_SFqcd_nb3", wsLB, translatePDF("pdf_prim_SFqcd_nb3"));

  printDiff_PDF_PDF(wsOAK, "betapdf_acc_Zee_M1", wsLB, translatePDF("betapdf_acc_Zee_M1"));
  printDiff_PDF_PDF(wsOAK, "betapdf_acc_Zee_M2", wsLB, translatePDF("betapdf_acc_Zee_M2"));
  printDiff_PDF_PDF(wsOAK, "betapdf_acc_Zee_M3", wsLB, translatePDF("betapdf_acc_Zee_M3"));
  printDiff_PDF_PDF(wsOAK, "betapdf_acc_Zee_M4", wsLB, translatePDF("betapdf_acc_Zee_M4"));
  printDiff_PDF_PDF(wsOAK, "betapdf_acc_Zmm_M1", wsLB, translatePDF("betapdf_acc_Zmm_M1"));
  printDiff_PDF_PDF(wsOAK, "betapdf_acc_Zmm_M2", wsLB, translatePDF("betapdf_acc_Zmm_M2"));
  printDiff_PDF_PDF(wsOAK, "betapdf_acc_Zmm_M3", wsLB, translatePDF("betapdf_acc_Zmm_M3"));
  printDiff_PDF_PDF(wsOAK, "betapdf_acc_Zmm_M4", wsLB, translatePDF("betapdf_acc_Zmm_M4"));

  printDiff_PDF_PDF(wsOAK, "pdf_prim_all_gu", wsLB, translatePDF("pdf_prim_all_gu"));
  
  printDiff_PDF_PDF(wsOAK, "pdf_prim_btageff_sf", wsLB, translatePDF("pdf_prim_btageff_sf"));
  printDiff_PDF_PDF(wsOAK, "pdf_prim_btageff_lf_sf", wsLB, translatePDF("pdf_prim_btageff_lf_sf"));

  printDiff_PDF_PDF(wsOAK, "pdf_prim_singletop_xsec", wsLB, translatePDF("pdf_prim_singletop_xsec"));
  printDiff_PDF_PDF(wsOAK, "pdf_prim_wjets_xsec", wsLB, translatePDF("pdf_prim_wjets_xsec"));

  printDiff_PDF_PDF(wsOAK, "pdf_prim_pdfsyst_sf", wsLB, translatePDF("pdf_prim_pdfsyst_sf"));

  printDiff_PDF_PDF(wsOAK, "betapdf_eff_Zee", wsLB, translatePDF("betapdf_eff_Zee"));
  printDiff_PDF_PDF(wsOAK, "betapdf_eff_Zmm", wsLB, translatePDF("betapdf_eff_Zmm"));
  
  printDiff_PDF_PDF(wsOAK, "pdf_prim_knn_1b_M1", wsLB, translatePDF("pdf_prim_knn_1b_M1"));
  printDiff_PDF_PDF(wsOAK, "pdf_prim_knn_1b_M2", wsLB, translatePDF("pdf_prim_knn_1b_M2"));
  printDiff_PDF_PDF(wsOAK, "pdf_prim_knn_1b_M3", wsLB, translatePDF("pdf_prim_knn_1b_M3"));
  printDiff_PDF_PDF(wsOAK, "pdf_prim_knn_1b_M4", wsLB, translatePDF("pdf_prim_knn_1b_M4"));
  printDiff_PDF_PDF(wsOAK, "pdf_prim_knn_2b", wsLB, translatePDF("pdf_prim_knn_2b"));
  printDiff_PDF_PDF(wsOAK, "pdf_prim_knn_3b", wsLB, translatePDF("pdf_prim_knn_3b"));

  printDiff_PDF_PDF(wsOAK, "betapdf_pur_Zee", wsLB, translatePDF("betapdf_pur_Zee"));
  printDiff_PDF_PDF(wsOAK, "betapdf_pur_Zmm", wsLB, translatePDF("betapdf_pur_Zmm"));

  printDiff_PDF_PDF(wsOAK, "pdf_prim_rar_vv_sf", wsLB, translatePDF("pdf_prim_rar_vv_sf"));

  printDiff_PDF_PDF(wsOAK, "pdf_prim_sf_ll", wsLB, translatePDF("pdf_prim_sf_ll"));
  
  printDiff_PDF_PDF(wsOAK, "pdf_prim_sf_mc", wsLB, translatePDF("pdf_prim_sf_mc"));
  
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M1_H1", wsLB, translatePDF("betapdf_trigeff_M1_H1"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M1_H2", wsLB, translatePDF("betapdf_trigeff_M1_H2"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M1_H3", wsLB, translatePDF("betapdf_trigeff_M1_H3"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M1_H4", wsLB, translatePDF("betapdf_trigeff_M1_H4"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M2_H1", wsLB, translatePDF("betapdf_trigeff_M2_H1"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M2_H2", wsLB, translatePDF("betapdf_trigeff_M2_H2"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M2_H3", wsLB, translatePDF("betapdf_trigeff_M2_H3"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M2_H4", wsLB, translatePDF("betapdf_trigeff_M2_H4"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M3_H1", wsLB, translatePDF("betapdf_trigeff_M3_H1"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M3_H2", wsLB, translatePDF("betapdf_trigeff_M3_H2"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M3_H3", wsLB, translatePDF("betapdf_trigeff_M3_H3"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M3_H4", wsLB, translatePDF("betapdf_trigeff_M3_H4"));
  //printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M4_H1", wsLB, translatePDF("betapdf_trigeff_M4_H1"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M4_H2", wsLB, translatePDF("betapdf_trigeff_M4_H2"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M4_H3", wsLB, translatePDF("betapdf_trigeff_M4_H3"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_M4_H4", wsLB, translatePDF("betapdf_trigeff_M4_H4"));

  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M1_H1", wsLB, translatePDF("betapdf_trigeff_sl_M1_H1"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M1_H2", wsLB, translatePDF("betapdf_trigeff_sl_M1_H2"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M1_H3", wsLB, translatePDF("betapdf_trigeff_sl_M1_H3"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M1_H4", wsLB, translatePDF("betapdf_trigeff_sl_M1_H4"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M2_H1", wsLB, translatePDF("betapdf_trigeff_sl_M2_H1"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M2_H2", wsLB, translatePDF("betapdf_trigeff_sl_M2_H2"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M2_H3", wsLB, translatePDF("betapdf_trigeff_sl_M2_H3"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M2_H4", wsLB, translatePDF("betapdf_trigeff_sl_M2_H4"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M3_H1", wsLB, translatePDF("betapdf_trigeff_sl_M3_H1"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M3_H2", wsLB, translatePDF("betapdf_trigeff_sl_M3_H2"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M3_H3", wsLB, translatePDF("betapdf_trigeff_sl_M3_H3"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M3_H4", wsLB, translatePDF("betapdf_trigeff_sl_M3_H4"));
  //printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M4_H1", wsLB, translatePDF("betapdf_trigeff_sl_M4_H1"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M4_H2", wsLB, translatePDF("betapdf_trigeff_sl_M4_H2"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M4_H3", wsLB, translatePDF("betapdf_trigeff_sl_M4_H3"));
  printDiff_PDF_PDF(wsOAK, "betapdf_trigeff_sl_M4_H4", wsLB, translatePDF("betapdf_trigeff_sl_M4_H4"));
  


  //for anders
  //cout << "anders" << endl;
  //printDiff_RRV_RRV(wsOAK, "mean_trigeff_sl_M2_H3", wsLB, "mean_oneLeptonTriggerEfficiency_M2_H3");
  //printDiff_RRV_RRV(wsOAK, "sigma_trigeff_sl_M2_H3", wsLB, "sigma_oneLeptonTriggerEfficiency_M2_H3");

  
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M1_H1", wsLB, translatePDF("pdf_N_Zee_M1_H1"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M1_H2", wsLB, translatePDF("pdf_N_Zee_M1_H2"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M1_H3", wsLB, translatePDF("pdf_N_Zee_M1_H3"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M1_H4", wsLB, translatePDF("pdf_N_Zee_M1_H4"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M2_H1", wsLB, translatePDF("pdf_N_Zee_M2_H1"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M2_H2", wsLB, translatePDF("pdf_N_Zee_M2_H2"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M2_H3", wsLB, translatePDF("pdf_N_Zee_M2_H3"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M2_H4", wsLB, translatePDF("pdf_N_Zee_M2_H4"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M3_H1", wsLB, translatePDF("pdf_N_Zee_M3_H1"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M3_H2", wsLB, translatePDF("pdf_N_Zee_M3_H2"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M3_H3", wsLB, translatePDF("pdf_N_Zee_M3_H3"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M3_H4", wsLB, translatePDF("pdf_N_Zee_M3_H4"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M4_H2", wsLB, translatePDF("pdf_N_Zee_M4_H2"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M4_H3", wsLB, translatePDF("pdf_N_Zee_M4_H3"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zee_M4_H4", wsLB, translatePDF("pdf_N_Zee_M4_H4"));

  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M1_H1", wsLB, translatePDF("pdf_N_Zmm_M1_H1"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M1_H2", wsLB, translatePDF("pdf_N_Zmm_M1_H2"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M1_H3", wsLB, translatePDF("pdf_N_Zmm_M1_H3"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M1_H4", wsLB, translatePDF("pdf_N_Zmm_M1_H4"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M2_H1", wsLB, translatePDF("pdf_N_Zmm_M2_H1"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M2_H2", wsLB, translatePDF("pdf_N_Zmm_M2_H2"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M2_H3", wsLB, translatePDF("pdf_N_Zmm_M2_H3"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M2_H4", wsLB, translatePDF("pdf_N_Zmm_M2_H4"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M3_H1", wsLB, translatePDF("pdf_N_Zmm_M3_H1"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M3_H2", wsLB, translatePDF("pdf_N_Zmm_M3_H2"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M3_H3", wsLB, translatePDF("pdf_N_Zmm_M3_H3"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M3_H4", wsLB, translatePDF("pdf_N_Zmm_M3_H4"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M4_H2", wsLB, translatePDF("pdf_N_Zmm_M4_H2"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M4_H3", wsLB, translatePDF("pdf_N_Zmm_M4_H3"));
  printDiff_PDF_PDF(wsOAK, "pdf_N_Zmm_M4_H4", wsLB, translatePDF("pdf_N_Zmm_M4_H4"));


  //special -- can remove later -- not pdfs -- just testing
  /*
  printZtoNuNu(wsOAK, "mu_znn_M1_H1_1b", wsLB, "M1_H1", "M1_1b");
  printZtoNuNu(wsOAK, "mu_znn_M1_H2_1b", wsLB, "M1_H2", "M1_1b");
  printZtoNuNu(wsOAK, "mu_znn_M1_H3_1b", wsLB, "M1_H3", "M1_1b");
  printZtoNuNu(wsOAK, "mu_znn_M1_H4_1b", wsLB, "M1_H4", "M1_1b");
  printZtoNuNu(wsOAK, "mu_znn_M2_H1_1b", wsLB, "M2_H1", "M2_1b");
  printZtoNuNu(wsOAK, "mu_znn_M2_H2_1b", wsLB, "M2_H2", "M2_1b");
  printZtoNuNu(wsOAK, "mu_znn_M2_H3_1b", wsLB, "M2_H3", "M2_1b");
  printZtoNuNu(wsOAK, "mu_znn_M2_H4_1b", wsLB, "M2_H4", "M2_1b");
  printZtoNuNu(wsOAK, "mu_znn_M3_H1_1b", wsLB, "M3_H1", "M3_1b");
  printZtoNuNu(wsOAK, "mu_znn_M3_H2_1b", wsLB, "M3_H2", "M3_1b");
  printZtoNuNu(wsOAK, "mu_znn_M3_H3_1b", wsLB, "M3_H3", "M3_1b");
  printZtoNuNu(wsOAK, "mu_znn_M3_H4_1b", wsLB, "M3_H4", "M3_1b");
  printZtoNuNu(wsOAK, "mu_znn_M4_H2_1b", wsLB, "M4_H2", "M4_1b");
  printZtoNuNu(wsOAK, "mu_znn_M4_H3_1b", wsLB, "M4_H3", "M4_1b");
  printZtoNuNu(wsOAK, "mu_znn_M4_H4_1b", wsLB, "M4_H4", "M4_1b");
  printDiff_RAR_RAR(wsOAK, "n_ee_M1_H1", wsLB, "diElectron_M1_H1_Yield");
  //printDiff_RAR_RAR(wsOAK, "n_ee_M1_H1", wsLB, "Ztoll_M1_H1");
  printDiff_RAR_RAR(wsOAK, "acc_Zee_M1", wsLB, "ZtoeeAcceptance_M1");
  printDiff_RAR_RAR(wsOAK, "eff_Zee", wsLB, "ZtoeeEfficiency");
  printDiff_RAR_RAR(wsOAK, "pur_Zee", wsLB, "ZtoeeInvPurity");
  printDiff_RAR_RAR(wsOAK, "sf_ee", wsLB, "ZtoeeSystematic");
  printDiff_RAR_RAR(wsOAK, "mu_znn_M1_H1_1b", wsLB, "ZtoNuNu_VeryLooseBtagYield_M1_H1");
  printDiff_RAR_RAR(wsOAK, "knn_1b_M1", wsLB, "ZtoNuNu_VeryLooseBtagYield_M1_H1");
  printDiff_RAR_RAR(wsOAK, "znnoverll_bfratio", wsLB, "ZtollOverZtoNuNuRatio");
  printDiff_RAR_RAR(wsOAK, "dataoverll_lumiratio", wsLB, "ZtollOverZtoNuNuRatio");
  */


  //more tests
  /*
  printDiff_RAR_RAR(wsOAK, "mu_ttwj_M4_H2_3b", wsLB, "zeroLepton_bin42_TopWJetsYield");
  printDiff_RAR_RAR(wsOAK, "mu_vv_M4_H2_3b", wsLB, "zeroLepton_bin42_DibosonYield");
  printDiff_RAR_RAR(wsOAK, "mu_znn_M4_H2_3b", wsLB, "zeroLepton_bin42_ZtoNuNuYield");
  printDiff_RAR_RAR(wsOAK, "mu_qcd_M4_H2_3b", wsLB, "zeroLepton_bin42_QCDYield");
  printDiff_RAR_RAR(wsOAK, "mu_qcd_ldp_M4_H2_3b", wsLB, "zeroLeptonLowDeltaPhiN_bin42_QCDYield");
  printDiff_RAR_RAR(wsOAK, "sf_qcd_M4_H2_3b", wsLB, "zeroLeptonQCDClosure_bin42");
  
  //qcd test
  printDiff_RAR_RAR(wsOAK, "qcd_0lepLDP_ratio_H2", wsLB, "lowDeltaPhiNScaling_H2");
  printDiff_RAR_RAR(wsOAK, "SFqcd_met4", wsLB, "lowDeltaPhiNMETScaleFactor_M4");
  printDiff_RAR_RAR(wsOAK, "SFqcd_nb3", wsLB, "lowDeltaPhiNBTagScaleFactor_3b");
  printDiff_RAR_RAR(wsOAK, "prim_SFqcd_met4", wsLB, "prim_lowDeltaPhiNMETScaleFactor_M4");
  printDiff_RAR_RAR(wsOAK, "mean_SFqcd_met4", wsLB, "mean_lowDeltaPhiNMETScaleFactor_M4");
  printDiff_RAR_RAR(wsOAK, "sigma_SFqcd_met4", wsLB, "sigma_lowDeltaPhiNMETScaleFactor_M4");
  */

  
  //susy tests
  /*
  printDiff_RAR_RAR(wsOAK, "all_gu", wsLB, "signalGlobalUncertainty");
  printDiff_RAR_RAR(wsOAK, "eff_sf_M4_H4_3b", wsLB, "zeroLeptonSignalError_bin48");
  printDiff_RAR_RAR(wsOAK, "btageff_sf_M4_H4_3b", wsLB, "zeroLeptonBTagEfficiencyError_bin48");
  printDiff_RAR_RAR(wsOAK, "sigma_btageff_sf_M4_H4_3b", wsLB, "sigma_zeroLeptonBTagEfficiencyError_bin48");
  printDiff_RAR_RAR(wsOAK, "btageff_lf_sf_M4_H4_3b", wsLB, "zeroLeptonLTagEfficiencyError_bin48");
  printDiff_RAR_RAR(wsOAK, "JES_sf_M4_H4_3b", wsLB, "zeroLeptonJesError_bin48");
  printDiff_RAR_RAR(wsOAK, "pdfsyst_sf_M4_H4_3b", wsLB, "zeroLeptonPdfError_bin48");
  
  printDiff_RAR_RAR(wsOAK, "mu_susy_M4_H4_3b", wsLB, "zeroLepton_bin48_SignalYield");
  printDiff_RAR_RAR(wsOAK, "mu_susymc_M4_H4_3b", wsLB, "zeroLepton_bin48_SignalFractionOAK");
  printDiff_RAR_RAR(wsOAK, "mu_susy_all0lep", wsLB, "zeroLepton_bin48_SignalYield");
  printDiff_RAR_RAR(wsOAK, "mu_susymc_all0lep", wsLB, "zeroLepton_bin48_SignalYield");
  */

  printDiff_RAR_RAR(wsOAK, "all_gu", wsLB, "signalGlobalUncertainty");
  printDiff_RAR_RAR(wsOAK, "eff_sf_sl_M4_H4_3b", wsLB, "oneLeptonSignalError_bin48");
  printDiff_RAR_RAR(wsOAK, "btageff_sf_sl_M4_H4_3b", wsLB, "oneLeptonBTagEfficiencyError_bin48");
  printDiff_RAR_RAR(wsOAK, "btageff_lf_sf_sl_M4_H4_3b", wsLB, "oneLeptonLTagEfficiencyError_bin48");
  printDiff_RAR_RAR(wsOAK, "JES_sf_sl_M4_H4_3b", wsLB, "oneLeptonJesError_bin48");
  printDiff_RAR_RAR(wsOAK, "pdfsyst_sf_sl_M4_H4_3b", wsLB, "oneLeptonPdfError_bin48");
  
  printDiff_RAR_RAR(wsOAK, "mu_susy_sl_M4_H4_3b", wsLB, "oneLepton_bin48_SignalYield");
  printDiff_RAR_RAR(wsOAK, "mu_susymc_sl_M4_H4_3b", wsLB, "oneLepton_bin48_SignalFractionOAK");
  printDiff_RAR_RAR(wsOAK, "mu_susy_all0lep", wsLB, "oneLepton_bin48_SignalYield");
  printDiff_RAR_RAR(wsOAK, "mu_susymc_all0lep", wsLB, "oneLepton_bin48_SignalYield");

}



void setVal_RRV_RRV(double val, RooWorkspace* ws1, TString name1, RooWorkspace* ws2, TString name2)
{
  cout << "Setting " << name1 << " " << name2 << endl;
  (ws1->var(name1))->setVal(val);
  (ws2->var(name2))->setVal(val);
  //cout << "Set" << endl;
}



void setVal_RAR_RRV(double val, RooWorkspace* ws1, TString name1, RooWorkspace* ws2, TString name2)
{
  cout << "Setting " << name1 << " " << name2 << endl;
  ((RooRealVar*)(ws1->function(name1)))->setVal(val);
  (ws2->var(name2))->setVal(val);
  //cout << "Set" << endl;
}

void setVal_RAR_RAR(double val, RooWorkspace* ws1, TString name1, RooWorkspace* ws2, TString name2)
{
  cout << "Setting " << name1 << " " << ws1->function(name1) << endl;
  ((RooRealVar*)(ws1->function(name1)))->setVal(val);
  cout << "Setting " << name2 << " " << ws2->function(name2) << endl;
  ((RooRealVar*)(ws2->function(name2)))->setVal(val);
  //  (ws2->var(name2))->setVal(val);
  //cout << "Set" << endl;
}


/*

void setZtoNuNu(RooWorkspace* ws1, TString name1, RooWorkspace* ws2, TString bin2, TString scale2)
{
  //OAK Znn
  ((RooRealVar*)(ws1->function(name1)))->setVal(1);
  
  //LB Znn = ZtoNuNuVeryLoose*knn
  RooRealVar* ZtoNuNuVeryLoose = ws2->var("ZtoNuNu_VeryLooseBtagYield_"+bin2);
  ZtoNuNuVeryLoose->setVal(1);
  RooRealVar* knn = ws2->var("zeroLeptonZtoNuNubTagScaling_"+scale2);
  knn->setVal(1);
  
  cout << "Set" << endl;
}

*/

 /*
void setZtoNuNu(RooWorkspace* ws1, TString name1, double val1, RooWorkspace* ws2, TString bin2, double val2, TString scale2, double val3)
{
  //  setZtoNuNu(wsOAK, "mu_znn_M1_H1_1b", 128, wsLB, "M1_H1", 128/0.4, "M1_1b", 0.4);

  //OAK Znn
  ((RooRealVar*)(ws1->function(name1)))->setVal(val1);
  
  //LB Znn = ZtoNuNuVeryLoose*knn
  RooRealVar* ZtoNuNuVeryLoose = ws2->var("ZtoNuNu_VeryLooseBtagYield_"+bin2);
  ZtoNuNuVeryLoose->setVal(val2);
  RooRealVar* knn = (RooRealVar*)ws2->function("prim_zeroLeptonZtoNuNubTagScaling_"+scale2);
  knn->setVal(val3);
  
  //cout << "Set" << endl;
}
 */

 void setZtoNuNu(RooWorkspace* ws1, TString name1, double val1, RooWorkspace* ws2, TString bin2, double val2, TString scale2, double val3)
 {
   //OAK Znn
   ((RooRealVar*)(ws1->function(name1)))->setVal(val1);
   
   //LB
   RooRealVar* ZtoNuNuVeryLoose = ws2->var("ZtoNuNu_VeryLooseBtagYield_"+bin2);
   RooAbsReal* knn = ws2->function("zeroLeptonZtoNuNubTagScaling_"+scale2);
   ZtoNuNuVeryLoose->setVal(val1/(knn->getVal()));
 }


void setSignal(RooWorkspace* ws1, RooWorkspace* ws2)
{

  RooRealVar* rv_mu_susy_all0lep = ws1->var("mu_susy_all0lep");
  rv_mu_susy_all0lep->setVal(200);
  //rv_mu_susy_all0lep->setVal(0);
  //cout << "mu_susy_all0lep " << rv_mu_susy_all0lep->getVal() << endl;
  
  RooRealVar* signalCrossSection = ws2->var("signalCrossSection");
  //signalCrossSection->setVal( 1000.0 * 200.0 / 1.9399 / 2740.09 );
  //signalCrossSection->setVal( 1000.0 * 200.0 / 1.9399 / 2651.3400760000008631 );
  signalCrossSection->setVal( 200.0 * 49996.0  / (19.399 * 13496.559610000000248) ) ;
  //signalCrossSection->setVal( 200.0 * 49996.0  / (19.399 * 13496.559619545936584) ) ;
  //signalCrossSection->setVal(0);
  //cout << "signalCrossSection " << signalCrossSection->getVal() << endl;

  //RooRealVar* rv_mu_susymc_all0lep = ws1->var("mu_susymc_all0lep");
  //rv_mu_susymc_all0lep->setVal(100);
  //cout << "mu_susymc_all0lep " << rv_mu_susymc_all0lep->getVal() << endl;

  //RooRealVar* luminosity = ws2->var("luminosity");
  //luminosity->setVal(1);
  //cout << "luminosity " << luminosity->getVal() << endl;
  
  

  //just the cross section
  //signalCrossSection->setVal();


  
  for(int m = 1; m<=4; m++)
    {
      for(int h = 1; h<=4; h++)
	{
	  if(m==4 && h==1) continue;
	  for(int b = 1; b<=3; b++)
	    {
	      TString binName = "M"; binName+=m; binName+="_H"; binName+=h; binName+="_"; binName+=b; binName+="b";
	      
	      RooRealVar* rv_mu_susymc = ws1->var("mu_susymc_"+binName);
	      cout <<  rv_mu_susymc->GetName() << " " << rv_mu_susymc->getVal() << endl;


	    }//b
	}//h
    }//m
    
  //cout << "Set" << endl;
}


void setSame()
{

  //Open workspaces
  TFile fOAK("pull/OAK.root", "READ");
  RooWorkspace* wsOAK = (RooWorkspace*)fOAK.Get("ws");
  
  TFile fLB("pull/LB.root", "READ");
  RooWorkspace *wsLB = (RooWorkspace*)fLB.Get("workspace");
  
  cout << "workspace pointers " << wsOAK << " " << wsLB << endl;

  //First the things that don't need special treatment
  //eff_sf_M1_H1_1b, eff_sf_ldp_M1_H1_1b, eff_sf_sl_M1_H1_1b, mu_qcd_ldp_M1_H1_1b, mu_ttwj_sl_M1_H1_1b, sf_qcd_M1_H1_1b, sf_ttwj_M1_H1_1b
  vector<TString> names;
  names.push_back("prim_eff_sf_");
  names.push_back("prim_eff_sf_ldp_");
  names.push_back("prim_eff_sf_sl_");
  names.push_back("mu_qcd_ldp_");
  names.push_back("mu_ttwj_sl_");
  names.push_back("prim_sf_qcd_");
  names.push_back("prim_sf_ttwj_");
  for(vector<TString>::iterator it = names.begin(); it != names.end(); ++it)
    {
      TString name = *it;
      //cout << name << endl;
      
      for(int m = 1; m<=4; m++)
	{
	  for(int h = 1; h<=4; h++)
	    {
	      if(m==4 && h==1) continue;
	      for(int b = 1; b<=3; b++)
		{
		  TString bin = "M"; bin+=m; bin+="_H"; bin+=h; bin+="_"; bin+=b; bin+="b";
		  
		  TString varName = name;
		  varName += bin;
		  //cout << varName << " " << translate(varName) << endl;
		  
		  if(name.Contains("eff_sf_") || name.Contains("sf_qcd_") || name.Contains("sf_ttwj_") )
		    {
		      //setVal_RAR_RRV(0.9, wsOAK, varName, wsLB, translate(varName));
		      setVal_RAR_RAR(0.9, wsOAK, varName, wsLB, translate(varName));
		    }
		  
		  //set these later
		  //if(name.Contains("mu_"))
		  // {
		  //    setVal_RAR_RAR(1, wsOAK, varName, wsLB, translate(varName));
		  // }
		  
		}//b
	    }//h
	}//m
      
    }//vector of names
  
  
  setVal_RRV_RRV(15064, wsOAK, "mu_qcd_ldp_M1_H1_1b", wsLB, translate("mu_qcd_ldp_M1_H1_1b"));
  setVal_RRV_RRV(4643, wsOAK, "mu_qcd_ldp_M1_H1_2b", wsLB, translate("mu_qcd_ldp_M1_H1_2b"));
  setVal_RRV_RRV(18, wsOAK, "mu_qcd_ldp_M1_H1_3b", wsLB, translate("mu_qcd_ldp_M1_H1_3b"));
  setVal_RRV_RRV(23386, wsOAK, "mu_qcd_ldp_M1_H2_1b", wsLB, translate("mu_qcd_ldp_M1_H2_1b"));
  setVal_RRV_RRV(5495, wsOAK, "mu_qcd_ldp_M1_H2_2b", wsLB, translate("mu_qcd_ldp_M1_H2_2b"));
  setVal_RRV_RRV(72, wsOAK, "mu_qcd_ldp_M1_H2_3b", wsLB, translate("mu_qcd_ldp_M1_H2_3b"));
  setVal_RRV_RRV(4847, wsOAK, "mu_qcd_ldp_M1_H3_1b", wsLB, translate("mu_qcd_ldp_M1_H3_1b"));
  setVal_RRV_RRV(1174, wsOAK, "mu_qcd_ldp_M1_H3_2b", wsLB, translate("mu_qcd_ldp_M1_H3_2b"));
  setVal_RRV_RRV(47, wsOAK, "mu_qcd_ldp_M1_H3_3b", wsLB, translate("mu_qcd_ldp_M1_H3_3b"));
  setVal_RRV_RRV(4366, wsOAK, "mu_qcd_ldp_M1_H4_1b", wsLB, translate("mu_qcd_ldp_M1_H4_1b"));
  setVal_RRV_RRV(629, wsOAK, "mu_qcd_ldp_M1_H4_2b", wsLB, translate("mu_qcd_ldp_M1_H4_2b"));
  setVal_RRV_RRV(64, wsOAK, "mu_qcd_ldp_M1_H4_3b", wsLB, translate("mu_qcd_ldp_M1_H4_3b"));
  setVal_RRV_RRV(5974, wsOAK, "mu_qcd_ldp_M2_H1_1b", wsLB, translate("mu_qcd_ldp_M2_H1_1b"));
  setVal_RRV_RRV(1640, wsOAK, "mu_qcd_ldp_M2_H1_2b", wsLB, translate("mu_qcd_ldp_M2_H1_2b"));
  setVal_RRV_RRV(164, wsOAK, "mu_qcd_ldp_M2_H1_3b", wsLB, translate("mu_qcd_ldp_M2_H1_3b"));
  setVal_RRV_RRV(15187, wsOAK, "mu_qcd_ldp_M2_H2_1b", wsLB, translate("mu_qcd_ldp_M2_H2_1b"));
  setVal_RRV_RRV(4385, wsOAK, "mu_qcd_ldp_M2_H2_2b", wsLB, translate("mu_qcd_ldp_M2_H2_2b"));
  setVal_RRV_RRV(247, wsOAK, "mu_qcd_ldp_M2_H2_3b", wsLB, translate("mu_qcd_ldp_M2_H2_3b"));
  setVal_RRV_RRV(4135, wsOAK, "mu_qcd_ldp_M2_H3_1b", wsLB, translate("mu_qcd_ldp_M2_H3_1b"));
  setVal_RRV_RRV(1014, wsOAK, "mu_qcd_ldp_M2_H3_2b", wsLB, translate("mu_qcd_ldp_M2_H3_2b"));
  setVal_RRV_RRV(56, wsOAK, "mu_qcd_ldp_M2_H3_3b", wsLB, translate("mu_qcd_ldp_M2_H3_3b"));
  setVal_RRV_RRV(3726, wsOAK, "mu_qcd_ldp_M2_H4_1b", wsLB, translate("mu_qcd_ldp_M2_H4_1b"));
  setVal_RRV_RRV(670, wsOAK, "mu_qcd_ldp_M2_H4_2b", wsLB, translate("mu_qcd_ldp_M2_H4_2b"));
  setVal_RRV_RRV(47, wsOAK, "mu_qcd_ldp_M2_H4_3b", wsLB, translate("mu_qcd_ldp_M2_H4_3b"));
  setVal_RRV_RRV(149, wsOAK, "mu_qcd_ldp_M3_H1_1b", wsLB, translate("mu_qcd_ldp_M3_H1_1b"));
  setVal_RRV_RRV(0.6, wsOAK, "mu_qcd_ldp_M3_H1_2b", wsLB, translate("mu_qcd_ldp_M3_H1_2b"));
  setVal_RRV_RRV(0.3, wsOAK, "mu_qcd_ldp_M3_H1_3b", wsLB, translate("mu_qcd_ldp_M3_H1_3b"));
  setVal_RRV_RRV(457, wsOAK, "mu_qcd_ldp_M3_H2_1b", wsLB, translate("mu_qcd_ldp_M3_H2_1b"));
  setVal_RRV_RRV(124, wsOAK, "mu_qcd_ldp_M3_H2_2b", wsLB, translate("mu_qcd_ldp_M3_H2_2b"));
  setVal_RRV_RRV(0.5, wsOAK, "mu_qcd_ldp_M3_H2_3b", wsLB, translate("mu_qcd_ldp_M3_H2_3b"));
  setVal_RRV_RRV(347, wsOAK, "mu_qcd_ldp_M3_H3_1b", wsLB, translate("mu_qcd_ldp_M3_H3_1b"));
  setVal_RRV_RRV(54, wsOAK, "mu_qcd_ldp_M3_H3_2b", wsLB, translate("mu_qcd_ldp_M3_H3_2b"));
  setVal_RRV_RRV(18, wsOAK, "mu_qcd_ldp_M3_H3_3b", wsLB, translate("mu_qcd_ldp_M3_H3_3b"));
  setVal_RRV_RRV(527, wsOAK, "mu_qcd_ldp_M3_H4_1b", wsLB, translate("mu_qcd_ldp_M3_H4_1b"));
  setVal_RRV_RRV(102, wsOAK, "mu_qcd_ldp_M3_H4_2b", wsLB, translate("mu_qcd_ldp_M3_H4_2b"));
  setVal_RRV_RRV(12, wsOAK, "mu_qcd_ldp_M3_H4_3b", wsLB, translate("mu_qcd_ldp_M3_H4_3b"));
  setVal_RRV_RRV(17, wsOAK, "mu_qcd_ldp_M4_H2_1b", wsLB, translate("mu_qcd_ldp_M4_H2_1b"));
  setVal_RRV_RRV(7, wsOAK, "mu_qcd_ldp_M4_H2_2b", wsLB, translate("mu_qcd_ldp_M4_H2_2b"));
  setVal_RRV_RRV(2, wsOAK, "mu_qcd_ldp_M4_H2_3b", wsLB, translate("mu_qcd_ldp_M4_H2_3b"));
  setVal_RRV_RRV(37, wsOAK, "mu_qcd_ldp_M4_H3_1b", wsLB, translate("mu_qcd_ldp_M4_H3_1b"));
  setVal_RRV_RRV(4, wsOAK, "mu_qcd_ldp_M4_H3_2b", wsLB, translate("mu_qcd_ldp_M4_H3_2b"));
  setVal_RRV_RRV(0.00002, wsOAK, "mu_qcd_ldp_M4_H3_3b", wsLB, translate("mu_qcd_ldp_M4_H3_3b"));
  setVal_RRV_RRV(103, wsOAK, "mu_qcd_ldp_M4_H4_1b", wsLB, translate("mu_qcd_ldp_M4_H4_1b"));
  setVal_RRV_RRV(7, wsOAK, "mu_qcd_ldp_M4_H4_2b", wsLB, translate("mu_qcd_ldp_M4_H4_2b"));
  setVal_RRV_RRV(0.5, wsOAK, "mu_qcd_ldp_M4_H4_3b", wsLB, translate("mu_qcd_ldp_M4_H4_3b"));
		 
  setVal_RRV_RRV(978, wsOAK, "mu_ttwj_sl_M1_H1_1b", wsLB, translate("mu_ttwj_sl_M1_H1_1b"));
  setVal_RRV_RRV(609, wsOAK, "mu_ttwj_sl_M1_H1_2b", wsLB, translate("mu_ttwj_sl_M1_H1_2b"));
  setVal_RRV_RRV(49, wsOAK, "mu_ttwj_sl_M1_H1_3b", wsLB, translate("mu_ttwj_sl_M1_H1_3b"));
  setVal_RRV_RRV(778, wsOAK, "mu_ttwj_sl_M1_H2_1b", wsLB, translate("mu_ttwj_sl_M1_H2_1b"));
  setVal_RRV_RRV(436, wsOAK, "mu_ttwj_sl_M1_H2_2b", wsLB, translate("mu_ttwj_sl_M1_H2_2b"));
  setVal_RRV_RRV(46, wsOAK, "mu_ttwj_sl_M1_H2_3b", wsLB, translate("mu_ttwj_sl_M1_H2_3b"));
  setVal_RRV_RRV(77, wsOAK, "mu_ttwj_sl_M1_H3_1b", wsLB, translate("mu_ttwj_sl_M1_H3_1b"));
  setVal_RRV_RRV(39, wsOAK, "mu_ttwj_sl_M1_H3_2b", wsLB, translate("mu_ttwj_sl_M1_H3_2b"));
  setVal_RRV_RRV(3, wsOAK, "mu_ttwj_sl_M1_H3_3b", wsLB, translate("mu_ttwj_sl_M1_H3_3b"));
  setVal_RRV_RRV(22, wsOAK, "mu_ttwj_sl_M1_H4_1b", wsLB, translate("mu_ttwj_sl_M1_H4_1b"));
  setVal_RRV_RRV(13, wsOAK, "mu_ttwj_sl_M1_H4_2b", wsLB, translate("mu_ttwj_sl_M1_H4_2b"));
  setVal_RRV_RRV(0.00001, wsOAK, "mu_ttwj_sl_M1_H4_3b", wsLB, translate("mu_ttwj_sl_M1_H4_3b"));
  setVal_RRV_RRV(1650, wsOAK, "mu_ttwj_sl_M2_H1_1b", wsLB, translate("mu_ttwj_sl_M2_H1_1b"));
  setVal_RRV_RRV(900, wsOAK, "mu_ttwj_sl_M2_H1_2b", wsLB, translate("mu_ttwj_sl_M2_H1_2b"));
  setVal_RRV_RRV(86, wsOAK, "mu_ttwj_sl_M2_H1_3b", wsLB, translate("mu_ttwj_sl_M2_H1_3b"));
  setVal_RRV_RRV(1455, wsOAK, "mu_ttwj_sl_M2_H2_1b", wsLB, translate("mu_ttwj_sl_M2_H2_1b"));
  setVal_RRV_RRV(897, wsOAK, "mu_ttwj_sl_M2_H2_2b", wsLB, translate("mu_ttwj_sl_M2_H2_2b"));
  setVal_RRV_RRV(70, wsOAK, "mu_ttwj_sl_M2_H2_3b", wsLB, translate("mu_ttwj_sl_M2_H2_3b"));
  setVal_RRV_RRV(154, wsOAK, "mu_ttwj_sl_M2_H3_1b", wsLB, translate("mu_ttwj_sl_M2_H3_1b"));
  setVal_RRV_RRV(91, wsOAK, "mu_ttwj_sl_M2_H3_2b", wsLB, translate("mu_ttwj_sl_M2_H3_2b"));
  setVal_RRV_RRV(13, wsOAK, "mu_ttwj_sl_M2_H3_3b", wsLB, translate("mu_ttwj_sl_M2_H3_3b"));
  setVal_RRV_RRV(70, wsOAK, "mu_ttwj_sl_M2_H4_1b", wsLB, translate("mu_ttwj_sl_M2_H4_1b"));
  setVal_RRV_RRV(36, wsOAK, "mu_ttwj_sl_M2_H4_2b", wsLB, translate("mu_ttwj_sl_M2_H4_2b"));
  setVal_RRV_RRV(5, wsOAK, "mu_ttwj_sl_M2_H4_3b", wsLB, translate("mu_ttwj_sl_M2_H4_3b"));
  setVal_RRV_RRV(262, wsOAK, "mu_ttwj_sl_M3_H1_1b", wsLB, translate("mu_ttwj_sl_M3_H1_1b"));
  setVal_RRV_RRV(111, wsOAK, "mu_ttwj_sl_M3_H1_2b", wsLB, translate("mu_ttwj_sl_M3_H1_2b"));
  setVal_RRV_RRV(9, wsOAK, "mu_ttwj_sl_M3_H1_3b", wsLB, translate("mu_ttwj_sl_M3_H1_3b"));
  setVal_RRV_RRV(438, wsOAK, "mu_ttwj_sl_M3_H2_1b", wsLB, translate("mu_ttwj_sl_M3_H2_1b"));
  setVal_RRV_RRV(193, wsOAK, "mu_ttwj_sl_M3_H2_2b", wsLB, translate("mu_ttwj_sl_M3_H2_2b"));
  setVal_RRV_RRV(119, wsOAK, "mu_ttwj_sl_M3_H2_3b", wsLB, translate("mu_ttwj_sl_M3_H2_3b"));
  setVal_RRV_RRV(59, wsOAK, "mu_ttwj_sl_M3_H3_1b", wsLB, translate("mu_ttwj_sl_M3_H3_1b"));
  setVal_RRV_RRV(32, wsOAK, "mu_ttwj_sl_M3_H3_2b", wsLB, translate("mu_ttwj_sl_M3_H3_2b"));
  setVal_RRV_RRV(2, wsOAK, "mu_ttwj_sl_M3_H3_3b", wsLB, translate("mu_ttwj_sl_M3_H3_3b"));
  setVal_RRV_RRV(25, wsOAK, "mu_ttwj_sl_M3_H4_1b", wsLB, translate("mu_ttwj_sl_M3_H4_1b"));
  setVal_RRV_RRV(12, wsOAK, "mu_ttwj_sl_M3_H4_2b", wsLB, translate("mu_ttwj_sl_M3_H4_2b"));
  setVal_RRV_RRV(1, wsOAK, "mu_ttwj_sl_M3_H4_3b", wsLB, translate("mu_ttwj_sl_M3_H4_3b"));
  setVal_RRV_RRV(134, wsOAK, "mu_ttwj_sl_M4_H2_1b", wsLB, translate("mu_ttwj_sl_M4_H2_1b"));
  setVal_RRV_RRV(51, wsOAK, "mu_ttwj_sl_M4_H2_2b", wsLB, translate("mu_ttwj_sl_M4_H2_2b"));
  setVal_RRV_RRV(5, wsOAK, "mu_ttwj_sl_M4_H2_3b", wsLB, translate("mu_ttwj_sl_M4_H2_3b"));
  setVal_RRV_RRV(42, wsOAK, "mu_ttwj_sl_M4_H3_1b", wsLB, translate("mu_ttwj_sl_M4_H3_1b"));
  setVal_RRV_RRV(14, wsOAK, "mu_ttwj_sl_M4_H3_2b", wsLB, translate("mu_ttwj_sl_M4_H3_2b"));
  setVal_RRV_RRV(2, wsOAK, "mu_ttwj_sl_M4_H3_3b", wsLB, translate("mu_ttwj_sl_M4_H3_3b"));
  setVal_RRV_RRV(27, wsOAK, "mu_ttwj_sl_M4_H4_1b", wsLB, translate("mu_ttwj_sl_M4_H4_1b"));
  setVal_RRV_RRV(11, wsOAK, "mu_ttwj_sl_M4_H4_2b", wsLB, translate("mu_ttwj_sl_M4_H4_2b"));
  setVal_RRV_RRV(1, wsOAK, "mu_ttwj_sl_M4_H4_3b", wsLB, translate("mu_ttwj_sl_M4_H4_3b"));

  setVal_RAR_RAR(-0.05, wsOAK, "prim_JES_sf", wsLB, translate("prim_JES_sf"));
  
  setVal_RAR_RAR(1.15, wsOAK, "SFqcd_met2", wsLB, translate("SFqcd_met2"));
  setVal_RAR_RAR(1.38, wsOAK, "prim_SFqcd_met3", wsLB, translate("prim_SFqcd_met3"));
  setVal_RAR_RAR(1.89, wsOAK, "prim_SFqcd_met4", wsLB, translate("prim_SFqcd_met4"));
  setVal_RAR_RAR(0.7, wsOAK, "SFqcd_nb2", wsLB, translate("SFqcd_nb2"));
  setVal_RAR_RAR(0.7, wsOAK, "prim_SFqcd_nb3", wsLB, translate("prim_SFqcd_nb3"));

  setVal_RAR_RAR(0.67, wsOAK, "acc_Zee_M1", wsLB, translate("acc_Zee_M1"));
  setVal_RAR_RAR(0.69, wsOAK, "acc_Zee_M2", wsLB, translate("acc_Zee_M2"));
  setVal_RAR_RAR(0.74, wsOAK, "acc_Zee_M3", wsLB, translate("acc_Zee_M3"));
  setVal_RAR_RAR(0.79, wsOAK, "acc_Zee_M4", wsLB, translate("acc_Zee_M4"));
  setVal_RAR_RAR(0.67, wsOAK, "acc_Zmm_M1", wsLB, translate("acc_Zmm_M1"));
  setVal_RAR_RAR(0.73, wsOAK, "acc_Zmm_M2", wsLB, translate("acc_Zmm_M2"));
  setVal_RAR_RAR(0.80, wsOAK, "acc_Zmm_M3", wsLB, translate("acc_Zmm_M3"));
  setVal_RAR_RAR(0.85, wsOAK, "acc_Zmm_M4", wsLB, translate("acc_Zmm_M4"));

  setVal_RAR_RAR(0.9, wsOAK, "prim_all_gu", wsLB, translate("prim_all_gu"));
  
  setVal_RAR_RAR(-0.09, wsOAK, "prim_btageff_sf", wsLB, translate("prim_btageff_sf"));
  setVal_RAR_RAR(-0.09, wsOAK, "prim_btageff_lf_sf", wsLB, translate("prim_btageff_lf_sf"));

  setVal_RAR_RAR(-0.05, wsOAK, "prim_singletop_xsec", wsLB, translate("prim_singletop_xsec"));
  setVal_RAR_RAR(0.05, wsOAK, "prim_wjets_xsec", wsLB, translate("prim_wjets_xsec"));
  setVal_RAR_RAR(0.05, wsOAK, "prim_pdfsyst_sf", wsLB, translate("prim_pdfsyst_sf"));

  setVal_RAR_RRV(0.85, wsOAK, "eff_Zee", wsLB, translate("eff_Zee"));
  setVal_RAR_RRV(0.88, wsOAK, "eff_Zmm", wsLB, translate("eff_Zmm"));



  setVal_RAR_RRV(0.9, wsOAK, "pur_Zee", wsLB, translate("pur_Zee"));
  setVal_RAR_RRV(0.77, wsOAK, "pur_Zmm", wsLB, translate("pur_Zmm"));

  setVal_RRV_RRV(0.27, wsOAK, "qcd_0lepLDP_ratio_H1", wsLB, translate("qcd_0lepLDP_ratio_H1"));
  setVal_RRV_RRV(0.17, wsOAK, "qcd_0lepLDP_ratio_H2", wsLB, translate("qcd_0lepLDP_ratio_H2"));
  setVal_RRV_RRV(0.11, wsOAK, "qcd_0lepLDP_ratio_H3", wsLB, translate("qcd_0lepLDP_ratio_H3"));
  setVal_RRV_RRV(0.05, wsOAK, "qcd_0lepLDP_ratio_H4", wsLB, translate("qcd_0lepLDP_ratio_H4"));

  setVal_RAR_RAR(0.9, wsOAK, "prim_rar_vv_sf", wsLB, translate("prim_rar_vv_sf"));

  setVal_RAR_RAR(-0.23, wsOAK, "prim_sf_ll", wsLB, translate("prim_sf_ll"));
  
  setVal_RAR_RAR(0.98, wsOAK, "prim_sf_mc", wsLB, translate("prim_sf_mc"));
  
  setVal_RAR_RRV(0.61, wsOAK, "trigeff_M1_H1", wsLB, translate("trigeff_M1_H1"));
  setVal_RAR_RRV(0.72, wsOAK, "trigeff_M1_H2", wsLB, translate("trigeff_M1_H2"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_M1_H3", wsLB, translate("trigeff_M1_H3"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_M1_H4", wsLB, translate("trigeff_M1_H4"));
  setVal_RAR_RRV(0.76, wsOAK, "trigeff_M2_H1", wsLB, translate("trigeff_M2_H1"));
  setVal_RAR_RRV(0.83, wsOAK, "trigeff_M2_H2", wsLB, translate("trigeff_M2_H2"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_M2_H3", wsLB, translate("trigeff_M2_H3"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_M2_H4", wsLB, translate("trigeff_M2_H4"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_M3_H1", wsLB, translate("trigeff_M3_H1"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_M3_H2", wsLB, translate("trigeff_M3_H2"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_M3_H3", wsLB, translate("trigeff_M3_H3"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_M3_H4", wsLB, translate("trigeff_M3_H4"));
  //setVal_RAR_RRV(1, wsOAK, "trigeff_M4_H1", wsLB, translate("trigeff_M4_H1"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_M4_H2", wsLB, translate("trigeff_M4_H2"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_M4_H3", wsLB, translate("trigeff_M4_H3"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_M4_H4", wsLB, translate("trigeff_M4_H4"));

  setVal_RAR_RRV(0.9, wsOAK, "trigeff_sl_M1_H1", wsLB, translate("trigeff_sl_M1_H1"));
  setVal_RAR_RRV(0.987, wsOAK, "trigeff_sl_M1_H2", wsLB, translate("trigeff_sl_M1_H2"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M1_H3", wsLB, translate("trigeff_sl_M1_H3"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M1_H4", wsLB, translate("trigeff_sl_M1_H4"));
  setVal_RAR_RRV(0.95, wsOAK, "trigeff_sl_M2_H1", wsLB, translate("trigeff_sl_M2_H1"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M2_H2", wsLB, translate("trigeff_sl_M2_H2"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M2_H3", wsLB, translate("trigeff_sl_M2_H3"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M2_H4", wsLB, translate("trigeff_sl_M2_H4"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M3_H1", wsLB, translate("trigeff_sl_M3_H1"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M3_H2", wsLB, translate("trigeff_sl_M3_H2"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M3_H3", wsLB, translate("trigeff_sl_M3_H3"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M3_H4", wsLB, translate("trigeff_sl_M3_H4"));
  //setVal_RAR_RRV(1, wsOAK, "trigeff_sl_M4_H1", wsLB, translate("trigeff_sl_M4_H1"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M4_H2", wsLB, translate("trigeff_sl_M4_H2"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M4_H3", wsLB, translate("trigeff_sl_M4_H3"));
  setVal_RAR_RRV(0.99, wsOAK, "trigeff_sl_M4_H4", wsLB, translate("trigeff_sl_M4_H4"));
  
  setVal_RAR_RAR(0.4, wsOAK, "prim_knn_1b_M1", wsLB, translate("prim_knn_1b_M1"));
  setVal_RAR_RAR(0.4, wsOAK, "prim_knn_1b_M2", wsLB, translate("prim_knn_1b_M2"));
  setVal_RAR_RAR(0.4, wsOAK, "prim_knn_1b_M3", wsLB, translate("prim_knn_1b_M3"));
  setVal_RAR_RAR(0.4, wsOAK, "prim_knn_1b_M4", wsLB, translate("prim_knn_1b_M4"));
  setVal_RAR_RAR(0.1, wsOAK, "prim_knn_2b", wsLB, translate("prim_knn_2b"));
  setVal_RAR_RAR(0.01, wsOAK, "prim_knn_3b", wsLB, translate("prim_knn_3b"));
  //NOTE -- LB Knn will be reset by setZtoNuNu function
  
  setZtoNuNu(wsOAK, "mu_znn_M1_H1_1b", 128, wsLB, "M1_H1", 128/0.4, "M1_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M1_H2_1b", 89, wsLB, "M1_H2", 89/0.4, "M1_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M1_H3_1b", 7, wsLB, "M1_H3", 7/0.4, "M1_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M1_H4_1b", 4, wsLB, "M1_H4", 4/0.4, "M1_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M2_H1_1b", 268, wsLB, "M2_H1", 268/0.4, "M2_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M2_H2_1b", 287, wsLB, "M2_H2", 287/0.4, "M2_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M2_H3_1b", 36, wsLB, "M2_H3", 36/0.4, "M2_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M2_H4_1b", 15, wsLB, "M2_H4", 15/0.4, "M2_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M3_H1_1b", 121, wsLB, "M3_H1", 121/0.4, "M3_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M3_H2_1b", 105, wsLB, "M3_H2", 105/0.4, "M3_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M3_H3_1b", 11, wsLB, "M3_H3", 11/0.4, "M3_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M3_H4_1b", 9, wsLB, "M3_H4", 9/0.4, "M3_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M4_H2_1b", 49, wsLB, "M4_H2", 49/0.4, "M4_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M4_H3_1b", 8, wsLB, "M4_H3", 8/0.4, "M4_1b", 0.4);
  setZtoNuNu(wsOAK, "mu_znn_M4_H4_1b", 9, wsLB, "M4_H4", 9/0.4, "M4_1b", 0.4);

  
  setSignal(wsOAK, wsLB);

  cout << "signal" << endl;
  RooRealVar* rv_mu_susy_all0lep = wsOAK->var("mu_susy_all0lep");
  cout << "rv_mu_susy_all0lep " << rv_mu_susy_all0lep->getVal() << endl;
  RooRealVar* rv_mu_susymc_all0lep = wsOAK->var("mu_susymc_all0lep");
  cout << "rv_mu_susymc_all0lep " << rv_mu_susymc_all0lep->getVal() << endl;

  
  setVal_RAR_RAR(1.21, wsOAK, "ttwj_0lep1lep_ratio", wsLB, translate("ttwj_0lep1lep_ratio"));



  if(0)
    {
      RooAbsPdf* likelihoodLB = wsLB->pdf( "likelihood" );
      RooDataSet* rdsLB = (RooDataSet*) wsLB->obj( "data" );
      //likelihoodLB->Print("v");
      //RooAbsReal* nllLB = likelihoodLB->createNLL(*rdsLB);
      //cout << nllLB->getVal() << endl;
      wsLB->Print();
      
      
      ModelConfig* modelConfig = (ModelConfig*) wsOAK->obj( "SbModel" ) ;
      RooAbsPdf* likelihoodOAK = modelConfig->GetPdf() ;
      RooDataSet* rdsOAK = (RooDataSet*) wsOAK->obj( "ra2b_observed_rds" ) ;
      //  RooAbsReal* nllOAK = likelihoodOAK->createNLL(*rdsOAK);
      //cout << nllOAK->getVal() << endl;
      wsOAK->Print();
    }


  comparePDFs(wsLB, wsOAK);

  fOAK.Close();
  fLB.Close();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compare Zero Lepton Yields


void getValueAndError(RooWorkspace* ws, RooFitResult* fitResult, TString varName, double &value, double &error)
{

  RooAbsReal* var = ws->function(varName);
  assert(var);
  
  value = var->getVal();
  error = var->getPropagatedError(*fitResult);

  //cout << varName << " " << value << " +- " << error << endl;

}



void compareZL() 
{
  gStyle->SetOptStat(0);

  //likelihood models
  TString model1 = "LB"; //options are "OAK" and "LB"
  TString wsname1 = "workspace";
  TString f1path = "pull/ws1.root";
  TString fitresultname1 = "fitresult_likelihood_data";

  TString model2 = "LB"; //options are "OAK" and "LB"
  TString wsname2 = "workspace";
  TString f2path = "pull/ws2.root";
  TString fitresultname2 = "fitresult_likelihood_data";

  //Open workspaces
  TFile f1(f1path, "READ");
  RooWorkspace* ws1 = (RooWorkspace*)f1.Get(wsname1);
  RooFitResult *fitResult1 = (RooFitResult*)f1.Get(fitresultname1);

  TFile f2(f2path, "READ");
  RooWorkspace *ws2 = (RooWorkspace*)f2.Get(wsname2);
  RooFitResult *fitResult2 = (RooFitResult*)f2.Get(fitresultname2);
  
  vector<TString> names;
  names.push_back("TopWJets");
  
  names.push_back("Diboson");
  names.push_back("ZtoNuNu");
  names.push_back("QCD");
  names.push_back("Signal");
  for(vector<TString>::iterator it = names.begin(); it != names.end(); ++it)
    {
      TString componentName = *it;
      cout << componentName << endl;
      
      TH1D h1("h1", componentName+" offset yield = yield - <yield>", 48, 0.5, 48.5);
      TH1D h2("h2", componentName+" offset yield = yield - <yield>", 48, 0.5, 48.5);

      TH1D hpercentdiff("hpercentdiff", componentName+" percent difference = 100*(y1-y2)/<y>", 48, 0.5, 48.5);
      TH1D hsignificance("hsignificance", componentName+" difference significance = (y1-y2)/sqrt(e1^2+e2^2)", 48, 0.5, 48.5);

      int binNum = 0;
      for(int m = 1; m<=4; m++)
	{
	  for(int h = 1; h<=4; h++)
	    {
	      for(int b = 1; b<=3; b++)
		{
		  binNum++;
		  TString binCheck = "bin"; binCheck += binNum;
		  if(m==4 && h==1) continue;

		  TString bin = "M"; bin+=m; bin+="_H"; bin+=h; bin+="_"; bin+=b; bin+="b";
		  TString binLB = binTranslate(bin);
		  assert(binLB == binCheck);//sanity check
		  
		  TString varName = "zeroLepton_";
		  varName += binLB;
		  varName += "_";
		  varName += componentName;
		  varName += "Yield";
		  //cout << varName << endl;
	       
		  double value1, error1;
		  getValueAndError(ws1, fitResult1, varName, value1, error1);

		  double value2, error2;
		  getValueAndError(ws2, fitResult2, varName, value2, error2);

		  //cout << value1 << " " << error1 << " " << value2 << " " << error2 << endl;

		  double averageValue = (value1+value2)/2.0;

		  h1.SetBinContent(binNum, value1-averageValue);
		  h1.SetBinError(binNum, error1);
		  h2.SetBinContent(binNum, value2-averageValue);
		  h2.SetBinError(binNum, error2);
		  
		  if( (value1!=0) && (value2!=0) ) hpercentdiff.SetBinContent(binNum, 100.0*(value1-value2)/((value1+value2)/2.0) );
		  if( (error1!=0) && (error2!=0) ) hsignificance.SetBinContent(binNum, (value1-value2)/(sqrt(error1*error1+error2*error2)) );

		  h1.GetXaxis()->SetBinLabel(binNum, bin);
		  h2.GetXaxis()->SetBinLabel(binNum, bin);
		  hpercentdiff.GetXaxis()->SetBinLabel(binNum, bin);
		  hsignificance.GetXaxis()->SetBinLabel(binNum, bin);

		}//b
	    }//h
	}//m


      h1.SetLineColor(kRed);
      h2.SetLineColor(kBlue);

      h1.SetMarkerColor(kRed);
      h2.SetMarkerColor(kBlue);
      
      h1.SetMarkerStyle(20);
      h2.SetMarkerStyle(20);
      hpercentdiff.SetMarkerStyle(20);
      hsignificance.SetMarkerStyle(20);

      h1.SetMarkerSize(1);
      h2.SetMarkerSize(1);
      hpercentdiff.SetMarkerSize(1);
      hsignificance.SetMarkerSize(1);

      //h1.SetMinimum(-5);
      //h1.SetMaximum(5);

      h1.GetXaxis()->LabelsOption("v");
      h2.GetXaxis()->LabelsOption("v");
      hpercentdiff.GetXaxis()->LabelsOption("v");
      hsignificance.GetXaxis()->LabelsOption("v");

      TCanvas c("c", "c", 800, 400);
      c.cd();
      h1.Draw("E1 P");
      h2.Draw("SAME E1 P");
      gPad->SetBottomMargin(0.15);
      c.Print("plots/"+componentName+"_yields.pdf");

      TCanvas cpercentdiff("cpercentdiff", "cpercentdiff", 800, 400);
      cpercentdiff.cd();
      hpercentdiff.Draw("E1 P");
      gPad->SetBottomMargin(0.15);
      cpercentdiff.Print("plots/"+componentName+"_percentdiff.pdf");

      TCanvas csignificance("csignificance", "csignificance", 800, 400);
      csignificance.cd();
      hsignificance.Draw("E1 P");
      gPad->SetBottomMargin(0.15);
      csignificance.Print("plots/"+componentName+"_significance.pdf");

    }//vector of names
  

  f1.Close();
  f2.Close();
}

















void compare()
{

  setSame();

  return;
}
