#include <iostream>
#include <fstream>

#include "TFile.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"

#include "fixParameters.h"

using namespace RooFit;
using namespace std; 


//Note: this code assumes a specific bin naming convention and relies on files created in a painful process :(


void test(TString fileName)
{
  TFile fws(fileName, "READ");
  RooWorkspace* ws = (RooWorkspace*)fws.Get("workspace");
  fixParameters(*ws);
  fws.Close();
}


TString translateBin(TString binName)
{
  //LB->OAK
  if(binName=="bin1") return "M1_H1_1b";
  if(binName=="bin2") return "M1_H1_2b";
  if(binName=="bin3") return "M1_H1_3b";
  if(binName=="bin4") return "M1_H2_1b";
  if(binName=="bin5") return "M1_H2_2b";
  if(binName=="bin6") return "M1_H2_3b";
  if(binName=="bin7") return "M1_H3_1b";
  if(binName=="bin8") return "M1_H3_2b";
  if(binName=="bin9") return "M1_H3_3b";
  if(binName=="bin10") return "M1_H4_1b";
  if(binName=="bin11") return "M1_H4_2b";
  if(binName=="bin12") return "M1_H4_3b";
  if(binName=="bin13") return "M2_H1_1b";
  if(binName=="bin14") return "M2_H1_2b";
  if(binName=="bin15") return "M2_H1_3b";
  if(binName=="bin16") return "M2_H2_1b";
  if(binName=="bin17") return "M2_H2_2b";
  if(binName=="bin18") return "M2_H2_3b";
  if(binName=="bin19") return "M2_H3_1b";
  if(binName=="bin20") return "M2_H3_2b";
  if(binName=="bin21") return "M2_H3_3b";
  if(binName=="bin22") return "M2_H4_1b";
  if(binName=="bin23") return "M2_H4_2b";
  if(binName=="bin24") return "M2_H4_3b";
  if(binName=="bin25") return "M3_H1_1b";
  if(binName=="bin26") return "M3_H1_2b";
  if(binName=="bin27") return "M3_H1_3b";
  if(binName=="bin28") return "M3_H2_1b";
  if(binName=="bin29") return "M3_H2_2b";
  if(binName=="bin30") return "M3_H2_3b";
  if(binName=="bin31") return "M3_H3_1b";
  if(binName=="bin32") return "M3_H3_2b";
  if(binName=="bin33") return "M3_H3_3b";
  if(binName=="bin34") return "M3_H4_1b";
  if(binName=="bin35") return "M3_H4_2b";
  if(binName=="bin36") return "M3_H4_3b";
  if(binName=="bin37") return "M4_H1_1b";
  if(binName=="bin38") return "M4_H1_2b";
  if(binName=="bin39") return "M4_H1_3b";
  if(binName=="bin40") return "M4_H2_1b";
  if(binName=="bin41") return "M4_H2_2b";
  if(binName=="bin42") return "M4_H2_3b";
  if(binName=="bin43") return "M4_H3_1b";
  if(binName=="bin44") return "M4_H3_2b";
  if(binName=="bin45") return "M4_H3_3b";
  if(binName=="bin46") return "M4_H4_1b";
  if(binName=="bin47") return "M4_H4_2b";
  if(binName=="bin48") return "M4_H4_3b";
  assert(0);
  return "";
}


void fixRRV(RooWorkspace &ws, TString name, double value)
{
  //cout << "Fixing " << name << " to " << value << endl;
  RooRealVar* myRRV = ws.var(name);
  //cout << " -- was " << myRRV->getVal() << endl;

  if(value>=myRRV->getMin() && value<=myRRV->getMax())
    {
      myRRV->setVal(value);
    }
  else if(value<myRRV->getMin())
    {
      cout << "Warning: Hit min.  Fixing " << name << " to " << myRRV->getMin() << endl;
      myRRV->setVal(myRRV->getMin());
    }
  else if(value>myRRV->getMax())
    {
      cout << "Warning: Hit max.  Fixing " << name << " to " << myRRV->getMax() << endl;
      myRRV->setVal(myRRV->getMax());
    }

  myRRV->setConstant();
}


void fixParameters(RooWorkspace &ws) 
{

  //RQCD, VeryLooseZYields, TopWMRYields, QCDLDPYields
  TString whatToFix = "RQCD-VeryLooseZYields-TopWMRYields-QCDLDPYields";
  cout << "whatToFix " << whatToFix << endl;

  //Get true values 
  //Have to run combineTrue function in pull/pull.C to generate this file, which uses stuff from running 3DCode
  map<TString,double> genMap;
  
  fstream fin("pull/combinedTrue.dat", ios::in);
  assert(fin.is_open());
  
  while(fin.good())
    {
      string sname;
      float value;
      fin >> sname >> value;
      TString name = sname;
      if(name=="") continue; 

      //cout << name << " " << value << endl;
      genMap.insert( pair<TString,float>(name,value) );
    }

  if(whatToFix.Contains("RQCD"))
    {

      fixRRV(ws, "lowDeltaPhiNScaling_H1", genMap["R_0lepLDP_HT_1"]);
      fixRRV(ws, "lowDeltaPhiNScaling_H2", genMap["R_0lepLDP_HT_2"]);
      fixRRV(ws, "lowDeltaPhiNScaling_H3", genMap["R_0lepLDP_HT_3"]);
      fixRRV(ws, "lowDeltaPhiNScaling_H4", genMap["R_0lepLDP_HT_4"]);
      
      fixRRV(ws, "lowDeltaPhiNMETScaleFactor_M2", genMap["SF_0lepLDP_MET_2"]);
      fixRRV(ws, "lowDeltaPhiNMETScaleFactor_M3", genMap["SF_0lepLDP_MET_3"]);
      fixRRV(ws, "lowDeltaPhiNMETScaleFactor_M4", genMap["SF_0lepLDP_MET_4"]);

      fixRRV(ws, "lowDeltaPhiNBTagScaleFactor_2b", genMap["SF_0lepLDP_b_2"]);
      fixRRV(ws, "lowDeltaPhiNBTagScaleFactor_3b", genMap["SF_0lepLDP_b_3"]);

    }

  if(whatToFix.Contains("VeryLooseZYields"))
    {

      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M1_H1", genMap["N_Znn_M1_H1"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M1_H2", genMap["N_Znn_M1_H2"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M1_H3", genMap["N_Znn_M1_H3"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M1_H4", genMap["N_Znn_M1_H4"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M2_H1", genMap["N_Znn_M2_H1"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M2_H2", genMap["N_Znn_M2_H2"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M2_H3", genMap["N_Znn_M2_H3"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M2_H4", genMap["N_Znn_M2_H4"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M3_H1", genMap["N_Znn_M3_H1"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M3_H2", genMap["N_Znn_M3_H2"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M3_H3", genMap["N_Znn_M3_H3"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M3_H4", genMap["N_Znn_M3_H4"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M4_H2", genMap["N_Znn_M4_H2"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M4_H3", genMap["N_Znn_M4_H3"]);
      fixRRV(ws, "ZtoNuNu_VeryLooseBtagYield_M4_H4", genMap["N_Znn_M4_H4"]);

    }
  
  for(int i = 1; i<=48; i++)
    {
      if(i==37 || i==38 || i==39) continue;

      TString binName = "bin";
      binName += i;

      if(whatToFix.Contains("QCDLDPYields"))
	{
	  	  
	  fixRRV(ws, "zeroLeptonLowDeltaPhiN_"+binName+"_QCDYield", genMap["N_ldp_qcd_"+translateBin(binName)]);

	}//qcd
      
      if(whatToFix.Contains("TopWMRYields"))
	{
	  
	  for(int j=1; j<=5; j++)
	    {
	      TString thetaName = "Theta";
	      thetaName += j;
	      
	      fixRRV(ws, "oneLooseLep_"+binName+"_"+thetaName+"_TopWJetsYield", genMap["N_oneLooseLep_"+thetaName+"_TopWJets_"+translateBin(binName)]);
	      fixRRV(ws, "oneTightMu_"+binName+"_"+thetaName+"_TopWJetsYield", genMap["N_oneTightMu_"+thetaName+"_TopWJets_"+translateBin(binName)]);
	    }//theta

	  fixRRV(ws, "twoLooseLep_"+binName+"_TopWJetsYield", genMap["N_twoLooseLep_TopWJets_"+translateBin(binName)]);
	  fixRRV(ws, "twoTightMu_"+binName+"_TopWJetsYield", genMap["N_twoTightMu_TopWJets_"+translateBin(binName)]);

	}//top+W

    }//bin
  


}

