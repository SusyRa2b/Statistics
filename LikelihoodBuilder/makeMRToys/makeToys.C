#include <iostream>
#include <fstream>
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
(or "root -l -b -q toysMR.C+" if your file is called list.txt)

inputList is a file with the MR input files that looks like e.g. 
smmctest001.txt
smmctest002.txt
smmctest003.txt
*/

void toysMR(TString inputList)
{

  TRandom3 rand(54321);
  
  //poorly written that re-reads central values each time -- should change if speed becomes issue
  int ntoys = 400;
  for(int i=0; i<ntoys; i++)
    {
      //Make directory for toy input files
      TString dirName = "toyMR_";
      dirName+=i;
      mkdir(dirName.Data(), S_IRGRP);
      
      //Setup up to read list of files to make toys from
      fstream fileListFile(inputList.Data(), ios::in);
      assert(fileListFile.is_open());
      
      while(!fileListFile.eof())
	{
	  
	  //Get name of file to make toys from
	  string fileName;
	  fileListFile >> fileName;
	  if(fileName=="") continue;
	  cout << fileName << endl;
	  
	  //Setup to read this file
	  fstream contentFile(fileName.c_str(), ios::in);
	  assert(contentFile.is_open());

	  //File with toy counts
	  TString outName = dirName;
	  outName += "/";
	  outName += fileName;
	  fstream contentFileOut(outName.Data(), ios::out);
	  assert(contentFileOut.is_open());
	  
	  while(!contentFile.eof())
	    {
	      string value_s;
	      contentFile>>value_s;
	      if(value_s=="") continue;

	      TString value_ts = value_s;
	      double value = value_ts.Atof();
	      int fluctuatedValue = rand.Poisson(value);
	      //cout << value << " " << fluctuatedValue << endl;
	      contentFileOut<<fluctuatedValue<<endl;
	    }
	  contentFile.close();
	  contentFileOut.close();

	}
      fileListFile.close();

    }//toys
}


void makeToys()
{

  toysMR("listMR.txt");

}
