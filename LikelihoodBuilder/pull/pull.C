#include <iostream>
#include <fstream>
#include <cassert>
#include <map>

#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFitResult.h"

using namespace std;

bool useAdjusted = true;

void styleHistAndPrint(TH1D* hist, TString extra = "")
{
  hist->SetFillColor(kGray+1);
  hist->GetXaxis()->SetTitle("pull");

  TString name = hist->GetName();
  TCanvas canvas("c_"+name, "c_"+name, 640, 480);
  canvas.cd();
  hist->Draw();
  if(extra != "") { name += "_"; name += extra; }
  name += ".pdf";
  canvas.SaveAs(name);
}

void readNominal(map<TString,float> &genMap) 
{

  fstream fin("generate.dat", ios::in);
  assert(fin.is_open());

  while(fin.good())
    {
      string sname;
      float value;
      fin >> sname >> value;
      TString name = sname;

      genMap.insert( pair<TString,float>(name,value) );
    }

}


void readAdjusted(map<TString,float> &genMap) 
{

  fstream fin("adjusted.dat", ios::in);
  assert(fin.is_open());

  while(fin.good())
    {
      string dummy;
      string sname;
      float value;
      fin >> dummy >> sname >> value;
      TString name = sname;
      //cout << name << " " << value << endl;

      genMap.insert( pair<TString,float>(name,value) );
    }

}


void combineTrue()
{

  map<TString,float> genMap;
  
  fstream fin("generate.dat", ios::in);
  assert(fin.is_open());
  
  while(fin.good())
    {
      string sname;
      float value;
      fin >> sname >> value;
      TString name = sname;
      if(name=="") continue;
      cout << "name = " << name << endl;
      genMap.insert( pair<TString,float>(name,value) );
    }


  fstream fin2("adjusted17_susy200.dat", ios::in);
  assert(fin2.is_open());
  
  while(fin2.good())
    {
      string dummy;
      string sname;
      float value;
      fin2 >> dummy >> sname >> value;
      TString name = sname;
      if(name=="") continue;

      if( genMap.find(name) == genMap.end()) { //not found
	//cout << "Adding " << name << endl;
	cout << "name = " << name << endl;
	genMap.insert( pair<TString,float>(name,value) ); //simply add to map
      }
      else { //found
	//cout << "Replacing " << name << endl;
	genMap[name] = value; //replace existing value
      }
      
    }
  
  fstream fout("combinedTrue.dat", ios::out | ios::trunc);
  assert(fout.is_open());

  for(map<TString,float>::const_iterator i = genMap.begin(); i != genMap.end(); ++i)
    {
      TString k = i->first;
      float v = i->second;
      fout << k << " " << v << endl;
    }
  
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

TString translateRegion(TString region) 
{
  if(region=="zeroLepton") return "0lep";
  if(region=="zeroLeptonLowDeltaPhiN") return "ldp";
  if(region=="twoLooseLep") return "twoLooseLep";
  if(region=="twoTightMu") return "twoTightMu";
  if(region=="oneLooseLep") return "oneLooseLep";
  if(region=="oneTightMu") return "oneTightMu";
  assert(0);
  return "";
}

TString translateContribution(TString contribution)
{
  if(contribution=="Signal") return "susy";
  if(contribution=="QCD") return "qcd";
  if(contribution=="ZtoNuNu") return "znn";
  if(contribution=="Diboson") return "vv";
  if(contribution=="tt") return "tt";
  if(contribution=="wjets") return "wjets";
  if(contribution=="TopWJets") return "TopWJets";
  assert(0);
  return "";
}

bool skipBin(int i)
{
  
  if(i==37 || i==38 || i==39) { return true; }
  else { return false; }

}



void wstest(TString fileName)
{
  cout << "wstest " << fileName << endl;
  //TFile fws(fileName, "READ");
   //TFile *fws = new TFile(fileName, "READ");
  
  TFile fws(fileName, "READ");
  RooWorkspace *ws;
  fws.GetObject("workspace", ws);
  ws->Delete();
  //delete ws;  
  fws.Close();
  
   //RooWorkspace* ws = (RooWorkspace*)fws.Get("workspace");
   //fws.Close();
}

void fillMaps(map<TString,float> genMap, map<TString,float> adjustedMap, TString fileName, map<TString,float> &toyValues, map<TString,float> &toyErrors, map<TString,float> &trueValues, bool skipTrue=false)
{
  
  TFile fws(fileName, "READ");
  RooWorkspace* ws = (RooWorkspace*)fws.Get("workspace");

  RooFitResult* fitResult = (RooFitResult*)fws.Get("fitresult_likelihood_data");
  cout << "fitResult " << fitResult << endl;

  std::vector<TString> region;
  region.push_back("zeroLepton");
  region.push_back("zeroLeptonLowDeltaPhiN");
  region.push_back("twoLooseLep");
  region.push_back("twoTightMu");
  region.push_back("oneLooseLep");
  region.push_back("oneTightMu");  
  
  std::vector<TString> contribution;
  //contribution.push_back("Signal");
  contribution.push_back("TopWJets");
  contribution.push_back("QCD");
  contribution.push_back("ZtoNuNu");
  contribution.push_back("Diboson");
  
  int numParameters = 0;
  for(int i=1; i<=48; i++)
    {
      if(skipBin(i)) continue;
      
      TString binName = "bin";
      binName+=i;
      //cout << binName << endl;
      
      for(unsigned int j=0; j<region.size(); j++)
	{

	  for(unsigned int k=0; k<contribution.size(); k++)
	    {
	      
	      //In MR SL/DL, only Signal and TopWJets
	      if( (region[j].Contains("LooseLep") || region[j].Contains("TightMu")) && !(contribution[k].Contains("TopWJets") || contribution[k].Contains("Signal"))) continue;
	      
	      std::vector<TString> theta;
	      if( region[j].Contains("oneLooseLep") || region[j].Contains("oneTightMu") )
		{
		  theta.push_back("Theta1_");
		  theta.push_back("Theta2_");
		  theta.push_back("Theta3_");
		  theta.push_back("Theta4_");
		  theta.push_back("Theta5_");
		}
	      else { theta.push_back(""); }

	      for(unsigned int l=0; l<theta.size(); l++)
		{
	 
		  // WORKSPACE
		  //////////////////////////////////////
		  TString varName = region[j];
		  varName+="_";
		  varName+=binName;
		  varName+="_";
		  varName+=theta[l];
		  varName+=contribution[k];
		  varName+="Yield";
		  //cout << varName << endl;
	
		  //TEST
		  //float value1 = 1;
		  //float error1 = 1;
		  //float mapValue1 = 1;
		  //toyValues.insert(  pair<TString,float>(varName,value1) );
		  //toyErrors.insert(  pair<TString,float>(varName,error1) );
		  //trueValues.insert( pair<TString,float>(varName,mapValue1) );
		  //continue;
		  
		  RooAbsReal* myRar = ws->function(varName);
		  float value = myRar->getVal();
		  float error = myRar->getPropagatedError(*fitResult);
		  ////float value = (ws->function(varName))->getVal();
		  ////float error = (ws->function(varName))->getPropagatedError(*fitResult);
		  //cout << "ws " << varName << " " << value << " " << error << endl;
		  
		  numParameters += 1;
		  toyValues.insert(  pair<TString,float>(varName,value) );
		  toyErrors.insert(  pair<TString,float>(varName,error) );

		  
		  if(skipTrue) continue;
		  // MAP
		  ////////////////////////////////////////
		  
		  //combine tt and wjets if TopWJets in zeroLepton or zeroLeptonLowDeltaPhiN
		  std::vector<TString> subContribution;
		  if(contribution[k].Contains("TopWJets") && (region[j].Contains("zeroLepton")) )
		     {
		       subContribution.push_back("tt");
		       subContribution.push_back("wjets");
		     }
		  else{ subContribution.push_back(contribution[k]); }
		  
		  float mapValue=0;
		  for(unsigned int m=0; m<subContribution.size(); m++)
		    {
		      TString mapName = "N_";
		      mapName+=translateRegion(region[j]);
		      mapName+="_";
		      mapName+=theta[l];
		      mapName+=translateContribution(subContribution[m]) ;
		      mapName+="_";
		      mapName+=translateBin(binName);
		      //cout << "map " << mapName << " " << genMap[mapName] << endl;
		      
		      double thisMapValue = genMap[mapName];
		      mapValue+= thisMapValue;
		      
		    }//subContribution loop
		  //cout << varName << " " << mapValue << endl;

		  if(useAdjusted)
		    {
		      TString mapName = "N_";
		      mapName+=translateRegion(region[j]);
		      mapName+="_";
		      mapName+=theta[l];
		      mapName+= (contribution[k]=="TopWJets" && region[j].Contains("zeroLepton")) ? "ttwj" : translateContribution(contribution[k]) ;
		      mapName+="_";
		      mapName+=translateBin(binName);
		      
		      if( adjustedMap.find(mapName.Data()) == adjustedMap.end()) {
			// not found
		      }
		      else {
			mapValue = adjustedMap[mapName];
			if(mapName.Contains("susy")) cout << "adjusted " << mapName << " " << mapValue << ", fit " << value << " +- " << error << endl;
		      }

		    }
		  
		  trueValues.insert( pair<TString,float>(varName,mapValue) );
		  
		}//theta loop
	      
	    }//contribution loop
	  
	}//region loop
      
    }//bin loop
  //delete ws;
  fws.Close();
  //delete fitResult;
  cout << "number of parameters: " << numParameters << endl;
}




//creates histograms of pulls of various categories from the parameters of one fit. the parameters this considers are defined in fillMaps.
void pull1(TString fileName, TString plotName = "")
{
  gStyle->SetStatY(0.98);
  gStyle->SetStatX(0.98);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3); 
  
  map<TString, float> genMap;
  readNominal(genMap);

  map<TString, float> adjustedMap;
  readAdjusted(adjustedMap);

  map<TString,float> trueValues;
  map<TString,float> toyValues;
  map<TString,float> toyErrors;
  fillMaps(genMap, adjustedMap, fileName, toyValues, toyErrors, trueValues, false);
  
  //Book histograms
  TFile fsave("pull1Histograms.root", "RECREATE");
  fsave.cd();

  //PULL
  float pullmin = -10;
  float pullmax = 10;

  TH1D* hPull_LDP_Signal = new TH1D("hPull_LDP_Signal", "hPull_LDP_Signal", 50, pullmin, pullmax);
  TH1D* hPull_LDP_TopWJets = new TH1D("hPull_LDP_TopWJets", "hPull_LDP_TopWJets", 50, pullmin, pullmax);
  TH1D* hPull_LDP_QCD = new TH1D("hPull_LDP_QCD", "hPull_LDP_QCD", 50, pullmin, pullmax);
  TH1D* hPull_LDP_ZtoNuNu = new TH1D("hPull_LDP_ZtoNuNu", "hPull_LDP_ZtoNuNu", 50, pullmin, pullmax);
  TH1D* hPull_LDP_Diboson = new TH1D("hPull_LDP_Diboson", "hPull_LDP_Diboson", 50, pullmin, pullmax);

  TH1D* hPull_ZL_Signal = new TH1D("hPull_ZL_Signal", "hPull_ZL_Signal", 50, pullmin, pullmax);
  TH1D* hPull_ZL_TopWJets = new TH1D("hPull_ZL_TopWJets", "hPull_ZL_TopWJets", 50, pullmin, pullmax);
  TH1D* hPull_ZL_QCD = new TH1D("hPull_ZL_QCD", "hPull_ZL_QCD", 50, pullmin, pullmax);
  TH1D* hPull_ZL_ZtoNuNu = new TH1D("hPull_ZL_ZtoNuNu", "hPull_ZL_ZtoNuNu", 50, pullmin, pullmax);
  TH1D* hPull_ZL_Diboson = new TH1D("hPull_ZL_Diboson", "hPull_ZL_Diboson", 50, pullmin, pullmax);
  
  TH1D* hPull_SL_Signal = new TH1D("hPull_SL_Signal", "hPull_SL_Signal", 50, pullmin, pullmax);
  TH1D* hPull_SL_TopWJets = new TH1D("hPull_SL_TopWJets", "hPull_SL_TopWJets", 50, pullmin, pullmax);

  TH1D* hPull_MRDL_Signal = new TH1D("hPull_MRDL_Signal", "hPull_MRDL_Signal", 50, pullmin, pullmax);
  TH1D* hPull_MRDL_TopWJets = new TH1D("hPull_MRDL_TopWJets", "hPull_MRDL_TopWJets", 50, pullmin, pullmax);  
  
  //create vector of number names
  std::vector<TString> valueNames;
  for(map<TString,float>::iterator it=trueValues.begin(); it!=trueValues.end(); it++)
    {
      valueNames.push_back((*it).first);
    }
  cout << "Number of variables = " << valueNames.size() << endl;
  
  
  //loop over all variables
  for(unsigned int i=0; i<valueNames.size(); i++)
    {
      TString thisValueName = valueNames[i];

      float trueValue = trueValues[thisValueName];
      float toyValue = toyValues[thisValueName];
      float toyError = toyErrors[thisValueName];
      
      float pull = (toyValue - trueValue)/toyError;
      //if(fabs(pull)>4 && toyError>0) cout << thisValueName << ", toyValue " << toyValue << ", trueValue " << trueValue << ", toyError " << toyError << ", pull " << pull << endl;
      if(fabs(pull)>4 && toyError>0) cout << thisValueName << ", true value = " << trueValue << ", fitted value = " << toyValue << " +/- " << toyError << ", pull = " << pull << endl; 

      if(thisValueName.Contains("oneLooseLep_bin20_Theta1_TopWJetsYield")) cout << "handpicked " << thisValueName << ", toyValue " << toyValue << ", trueValue " << trueValue << ", toyError " << toyError << ", pull " << pull << endl;

      //Fill Histograms
      if(thisValueName.Contains("zeroLeptonLowDeltaPhiN") && thisValueName.Contains("Signal")) hPull_LDP_Signal->Fill(pull);
      if(thisValueName.Contains("zeroLeptonLowDeltaPhiN") && thisValueName.Contains("TopWJets")) hPull_LDP_TopWJets->Fill(pull);
      if(thisValueName.Contains("zeroLeptonLowDeltaPhiN") && thisValueName.Contains("QCD")) hPull_LDP_QCD->Fill(pull);
      if(thisValueName.Contains("zeroLeptonLowDeltaPhiN") && thisValueName.Contains("ZtoNuNu")) hPull_LDP_ZtoNuNu->Fill(pull);
      if(thisValueName.Contains("zeroLeptonLowDeltaPhiN") && thisValueName.Contains("Diboson")) hPull_LDP_Diboson->Fill(pull);
      
      if(thisValueName.Contains("zeroLepton_") && thisValueName.Contains("Signal")) hPull_ZL_Signal->Fill(pull);
      if(thisValueName.Contains("zeroLepton_") && thisValueName.Contains("TopWJets")) hPull_ZL_TopWJets->Fill(pull);
      if(thisValueName.Contains("zeroLepton_") && thisValueName.Contains("QCD")) hPull_ZL_QCD->Fill(pull);
      if(thisValueName.Contains("zeroLepton_") && thisValueName.Contains("ZtoNuNu")) hPull_ZL_ZtoNuNu->Fill(pull);
      if(thisValueName.Contains("zeroLepton_") && thisValueName.Contains("Diboson")) hPull_ZL_Diboson->Fill(pull);
      
      if((thisValueName.Contains("oneLooseLep") || thisValueName.Contains("oneTightMu")) && thisValueName.Contains("Signal")) hPull_SL_Signal->Fill(pull);
      if((thisValueName.Contains("oneLooseLep") || thisValueName.Contains("oneTightMu")) && thisValueName.Contains("TopWJets")) hPull_SL_TopWJets->Fill(pull);

      if(thisValueName.Contains("two") && thisValueName.Contains("Signal")) hPull_MRDL_Signal->Fill(pull);
      if(thisValueName.Contains("two") && thisValueName.Contains("TopWJets")) hPull_MRDL_TopWJets->Fill(pull);
      
    }//loop over variables
  
  //Style and print to pdf
  styleHistAndPrint(hPull_LDP_Signal, plotName);
  styleHistAndPrint(hPull_LDP_TopWJets, plotName);
  styleHistAndPrint(hPull_LDP_QCD, plotName);
  styleHistAndPrint(hPull_LDP_ZtoNuNu, plotName);
  styleHistAndPrint(hPull_LDP_Diboson, plotName);
  styleHistAndPrint(hPull_ZL_Signal, plotName);
  styleHistAndPrint(hPull_ZL_TopWJets, plotName);
  styleHistAndPrint(hPull_ZL_QCD, plotName);
  styleHistAndPrint(hPull_ZL_ZtoNuNu, plotName);
  styleHistAndPrint(hPull_ZL_Diboson, plotName);
  styleHistAndPrint(hPull_SL_TopWJets, plotName);
  styleHistAndPrint(hPull_MRDL_TopWJets, plotName);
  styleHistAndPrint(hPull_SL_Signal, plotName);
  styleHistAndPrint(hPull_MRDL_Signal, plotName);

  //Save Histograms
  hPull_LDP_TopWJets->Write();
  hPull_LDP_QCD->Write();
  hPull_LDP_ZtoNuNu->Write();
  hPull_LDP_Diboson->Write();

  hPull_ZL_TopWJets->Write();
  hPull_ZL_QCD->Write();
  hPull_ZL_ZtoNuNu->Write();
  hPull_ZL_Diboson->Write();

  hPull_SL_TopWJets->Write();
  hPull_SL_Signal->Write();

  hPull_MRDL_TopWJets->Write();
  hPull_MRDL_Signal->Write();


  fsave.Close();
  
}


//creates a histogram of pulls for each parameter defined in fillMaps.  each entry is one toy.
void pull2(TString fileList="toylist.dat")
{
  
  //get map of true values
  map<TString, float> genMap;
  readNominal(genMap);
  map<TString, float> adjustedMap;
  readAdjusted(adjustedMap);
  
  //loop through list of root files
  // -- create vector of file names
  std::vector<TString> fileNames;

  fstream fin(fileList.Data(), ios::in);
  assert(fin.is_open());

  while(fin.good())
    {
      string sname;
      fin >> sname;
      TString name = sname;
      if(name != "") fileNames.push_back(name);
    }
  cout << "Number of toys = " << fileNames.size() << endl;
  fin.close();

  
  //loop through vector  of file names
  // -- call function to open file, put all numbers into a map, put map into a vector
  bool skipTrue = false;
  map<TString,float> trueValues;
  std::vector< map<TString,float> > allToyValues;
  std::vector< map<TString,float> > allToyErrors;
  for(unsigned int i=0; i<fileNames.size(); i++)
    {
      cout << "Filling map for toy " << i << endl;
      map<TString,float> toyValues;
      map<TString,float> toyErrors;
      if(i != 0) skipTrue=true;
      fillMaps(genMap, adjustedMap, fileNames[i], toyValues, toyErrors, trueValues, skipTrue);
      allToyValues.push_back(toyValues);
      allToyErrors.push_back(toyErrors);
    }

  
  //create vector of number names
  std::vector<TString> valueNames;
  for(map<TString,float>::iterator it=trueValues.begin(); it!=trueValues.end(); it++)
    {
      valueNames.push_back((*it).first);
    }
  cout << "Number of variables = " << valueNames.size() << endl;


  //draw histogram of all pulls
  // -- declare histograms in loop over vector of number names
  // -- fill by looping over maps
  // -- draw true value
  // -- save histogram in root file

  TFile fpull2("pull2Histograms_"+fileList+".root", "RECREATE");
  fpull2.cd();

  //loop over all variables
  for(unsigned int i=0; i<valueNames.size(); i++)
    {
      TString thisValueName = valueNames[i];

      float trueValue = trueValues[thisValueName];

      float pullmin = -10;
      float pullmax = 10;
      //TH1D hPull("hPull_"+thisValueName, "hPull_"+thisValueName, 50, pullmin, pullmax);
      TH1D hPull("h_"+thisValueName+"_Pull", "h_"+thisValueName+"_Pull", 50, pullmin, pullmax);
      TH1D hDiff("h_"+thisValueName+"_Diff", "h_"+thisValueName+"_Diff", 200, -100, 100);

      //loop over all toys
      for(unsigned int j=0; j<allToyValues.size(); j++)
	{
	  map<TString,float> toyValues = allToyValues[j];
	  map<TString,float> toyErrors = allToyErrors[j];
	  
	  float toyValue = toyValues[thisValueName];
	  float toyError = toyErrors[thisValueName];
	  
	  float pull = (toyValue - trueValue)/toyError;
	  float diff = toyValue - trueValue;

	  hPull.Fill(pull);
	  hDiff.Fill(diff);
	  
	}//loop over toys
      
      hPull.Write();
      hDiff.Write();

    }//loop over variables

  fpull2.Close();

}//end of pull2

void pull()
{
  
}
