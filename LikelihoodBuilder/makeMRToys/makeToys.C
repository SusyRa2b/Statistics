#include <iostream>
#include <fstream>
#include <map>
#include <cassert>

#include <sys/stat.h>

#include "TString.h"
#include "TRandom3.h"
using namespace std;

/*
  to run, do
  >root -l
  >.L toysMR.C+
  >toysMR("yourListFileDescribedBelow.dat")
  (or "root -l -b -q toysMR.C+" if your file is called listMR.txt)
  
  NOTE THIS ASSUMES THE SIMPLE TAU METHOD IS BEING USED
  -- Does not handle the "outside" bins

  inputList is a file with the MR input files that looks like e.g. 
  topwfile001.txt signalfile001.txt
  topwfile002.txt signalfile002.txt
  ...
*/



struct MRmeans {
  double oneTightMu_Theta1;
  double oneLooseLep_Theta1;
  double oneTightMu_Theta2;
  double oneLooseLep_Theta2;
  double oneTightMu_Theta3;
  double oneLooseLep_Theta3;
  double oneTightMu_Theta4;
  double oneLooseLep_Theta4;
  double oneTightMu_Theta5;
  double oneLooseLep_Theta5;
  double twoTightMu;
  double twoLooseLep;
};


double getMean(double topW, double signal)
{
  //double lumiScale = 19.399/12.0; // applied to topW only
  double lumiScale = 1.0;

  //T1tttt, 1175 400
  //double signalCrossSection = 200 * 49996.0  / (19.399 * 13496.559610000000248);
  double signalCrossSection = 0.0;
  double lumi = 19.399;
  double Ntot = 49996.0;

  //topW is in events weighted for lumi
  double topMean = topW*lumiScale;

  //signal is in events, not weighted for lumi. Apply lumi weighting here.
  double signalMean = signal * signalCrossSection * lumi / Ntot;
  
  //cout << "topMean = " << topMean << ", signalMean = " << signalMean << endl;
  double mean = topMean + signalMean;
  return mean;
}

void toysMR(TString inputList)
{

  TRandom3 rand(54321);
  
  int ntoys = 400;
  
  //READ INPUTS INTO MEMORY

  map<TString, MRmeans> allTopWMeans;
  map<TString, MRmeans> allSignalMeans;

  //Setup up to read list of files to make toys from
  fstream fileListFile(inputList.Data(), ios::in);
  assert(fileListFile.is_open());
  
  while(!fileListFile.eof())
    {
      
      //GET FILE NAMES FOR THIS M/H/b BIN
      string fileName;
      string signalFileName;
      fileListFile >> fileName >> signalFileName;
      if(fileName=="") continue;
      cout << "-" << fileName << "- -" << signalFileName << "-" << endl;
  
      
      //Top+W INPUTS
      
      MRmeans topWMeans;

      //Setup to read this M/H/b bin file
      fstream contentFile(fileName.c_str(), ios::in);
      assert(contentFile.is_open());

      unsigned int lineCount = 0;
      while(!contentFile.eof())
	{
	  string value_s;
	  contentFile>>value_s;
	  if(value_s=="") continue;
	  lineCount++;

	  TString value_ts = value_s;
	  double value = value_ts.Atof();
	  
	  //cout << value << endl;
	  if(lineCount==1){ }
	  else if(lineCount==2){ topWMeans.oneTightMu_Theta1 = value ; }
	  else if(lineCount==3){ topWMeans.oneLooseLep_Theta1 = value ; }
	  else if(lineCount==4){ topWMeans.oneTightMu_Theta2 = value ; }
	  else if(lineCount==5){ topWMeans.oneLooseLep_Theta2 = value ; }
	  else if(lineCount==6){ topWMeans.oneTightMu_Theta3 = value ; }
	  else if(lineCount==7){ topWMeans.oneLooseLep_Theta3 = value ; }
	  else if(lineCount==8){ topWMeans.oneTightMu_Theta4 = value ; }
	  else if(lineCount==9){ topWMeans.oneLooseLep_Theta4 = value ; }
	  else if(lineCount==10){ topWMeans.oneTightMu_Theta5 = value ; }
	  else if(lineCount==11){ topWMeans.oneLooseLep_Theta5 = value ; }
	  else if(lineCount==12){ topWMeans.twoTightMu = value ; }
	  else if(lineCount==13){ topWMeans.twoLooseLep = value ; }
	  else {assert(0);}
	  
	}
      contentFile.close();
      allTopWMeans.insert(pair<TString,MRmeans>(fileName, topWMeans));
      
      
      //SIGNAL INPUTS
      
      MRmeans signalMeans;
      
      //Setup to read this M/H/b bin file
      fstream signalContentFile(signalFileName.c_str(), ios::in);
      assert(signalContentFile.is_open());
      
      while(!signalContentFile.eof())
	{

	  string binName_s, value_s;
	  signalContentFile>>binName_s>>value_s;
	  if(binName_s=="") continue;
	  
	  TString binName = binName_s;
	  TString value_ts = value_s;
	  double value = value_ts.Atof();
	  
	  //cout << binName << " " << value << endl;
	  if(binName.Contains("SignalCountError")) continue;
	  else if(binName.Contains("oneTightMu_Theta1")) signalMeans.oneTightMu_Theta1 = value;
	  else if(binName.Contains("oneTightMu_Theta2")) signalMeans.oneTightMu_Theta2 = value;
	  else if(binName.Contains("oneTightMu_Theta3")) signalMeans.oneTightMu_Theta3 = value;
	  else if(binName.Contains("oneTightMu_Theta4")) signalMeans.oneTightMu_Theta4 = value;
	  else if(binName.Contains("oneTightMu_Theta5")) signalMeans.oneTightMu_Theta5 = value;
	  else if(binName.Contains("oneLooseLep_Theta1")) signalMeans.oneLooseLep_Theta1 = value;
	  else if(binName.Contains("oneLooseLep_Theta2")) signalMeans.oneLooseLep_Theta2 = value;
	  else if(binName.Contains("oneLooseLep_Theta3")) signalMeans.oneLooseLep_Theta3 = value;
	  else if(binName.Contains("oneLooseLep_Theta4")) signalMeans.oneLooseLep_Theta4 = value;
	  else if(binName.Contains("oneLooseLep_Theta5")) signalMeans.oneLooseLep_Theta5 = value;
	  else if(binName.Contains("twoTightMu")) signalMeans.twoTightMu = value;
	  else if(binName.Contains("twoLooseLep")) signalMeans.twoLooseLep = value;
	  else{assert(0);}
	  
	}
      signalContentFile.close();
      allSignalMeans.insert(pair<TString,MRmeans>(fileName, signalMeans));

      
    }
  fileListFile.close();
  cout << "allTopWMeans size = " << allTopWMeans.size() << endl;
  cout << "allSignalMeans size = " << allSignalMeans.size() << endl;
  

  //MAKE TOYS

  for(int i=0; i<ntoys; i++)
    {
      //Make directory for toy input files
      TString dirName = "toyMR_";
      dirName+=i;
      mkdir(dirName.Data(), S_IRGRP);

      //loop over files
      for ( map<TString,MRmeans>::iterator it=allTopWMeans.begin() ; it != allTopWMeans.end(); it++ ) 
	{

	  TString fileName = (*it).first;
	  MRmeans topWMeans = (*it).second;
	  MRmeans signalMeans = allSignalMeans[fileName];

	  //File with toy counts
	  TString outName = dirName;
	  outName += "/";
	  outName += fileName;
	  fstream contentFileOut(outName.Data(), ios::out);
	  assert(contentFileOut.is_open());
	  

	  contentFileOut << "0" << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.oneTightMu_Theta1,  signalMeans.oneTightMu_Theta1)  ) << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.oneLooseLep_Theta1, signalMeans.oneLooseLep_Theta1) ) << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.oneTightMu_Theta2,  signalMeans.oneTightMu_Theta2)  ) << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.oneLooseLep_Theta2, signalMeans.oneLooseLep_Theta2) ) << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.oneTightMu_Theta3,  signalMeans.oneTightMu_Theta3)  ) << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.oneLooseLep_Theta3, signalMeans.oneLooseLep_Theta3) ) << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.oneTightMu_Theta4,  signalMeans.oneTightMu_Theta4)  ) << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.oneLooseLep_Theta4, signalMeans.oneLooseLep_Theta4) ) << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.oneTightMu_Theta5,  signalMeans.oneTightMu_Theta5)  ) << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.oneLooseLep_Theta5, signalMeans.oneLooseLep_Theta5) ) << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.twoTightMu,  signalMeans.twoTightMu)  ) << endl;
	  contentFileOut << rand.Poisson( getMean(topWMeans.twoLooseLep, signalMeans.twoLooseLep) ) << endl;
	  
	  contentFileOut.close();
	}

    }//toys
  
}


void makeToys()
{

  toysMR("listMR.txt");

}
