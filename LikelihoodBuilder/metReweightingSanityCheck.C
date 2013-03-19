#include <iostream>
#include <fstream>
#include <cassert>

#include "TString.h"

using namespace std;

float get_exp_0lep_ttwj_mr(int bm, int bh, int bb, float sltrigeff) {
  bb=bb+1;//Kristen's convention
  
  TString binLabel = "";
  binLabel += bh; binLabel+=bm; binLabel+=bb;
  
  TString binLabelBen = "M";
  binLabelBen += bm+1;
  binLabelBen += "_H";
  binLabelBen += bh+1;
  binLabelBen += "_";
  binLabelBen += bb;
  binLabelBen += "b";

  //TString tauhadLabel = "00";
  //tauhadLabel += bb;
  TString tauhadLabel = binLabel;

  //content file
  TString contentFileName = "goodinputsMR/smmcttwt";
  contentFileName += binLabel;
  contentFileName += ".txt";
  
  ifstream contentFile(contentFileName.Data(),fstream::in);
  assert(contentFile.is_open());
  
  double zlcount, tmcounts[5], llcounts[5], tm2count, ll2count;
  contentFile>>zlcount;
  for(int i=0; i<5; i++) {
    contentFile>>tmcounts[i];
    contentFile>>llcounts[i];
  }
  contentFile>>tm2count;
  contentFile>>ll2count;
  
  
  //scale factor file
  TString scaleFactorFileName = "goodinputsMR/polarizationScaleFactorsMT";
  scaleFactorFileName += binLabel;
  scaleFactorFileName += ".txt";
  
  ifstream scaleFactorFile(scaleFactorFileName.Data(),fstream::in);
  assert(scaleFactorFile.is_open());
  
  double polscalefactors[5], dummy, dilepscalefactor;
  for(int i=0; i<5; i++) {
    scaleFactorFile>>polscalefactors[i]>>dummy;
  }
  scaleFactorFile>>dilepscalefactor>>dummy;
  
  
  //tau->had scale factors
  TString tauhadScaleFactorFileName = "goodinputsMR/tauhadScaleFactorsMT";
  tauhadScaleFactorFileName += tauhadLabel;
  tauhadScaleFactorFileName += ".txt";
  
  ifstream tauhadScaleFactorFile(tauhadScaleFactorFileName.Data(),fstream::in);
  assert(tauhadScaleFactorFile.is_open());
  
  double onetauscalefactor, twotauscalefactor;
  tauhadScaleFactorFile>>onetauscalefactor>>dummy;
  tauhadScaleFactorFile>>twotauscalefactor>>dummy;
  
  
  //zero lepton prediction
  double exp_0lep_ttwj = 0;
  
  for(int i=0; i<5; i++) {
    exp_0lep_ttwj += polscalefactors[i]*(tmcounts[i] + llcounts[i]);
    exp_0lep_ttwj += onetauscalefactor*tmcounts[i];
    
    /*
    cout << binLabel << " " << "DEBUG i=" << i << " polscalefactor " << polscalefactors[i] << endl;
    cout << binLabel << " " << "DEBUG i=" << i << " tmcount " << tmcounts[i] << endl;
    cout << binLabel << " " << "DEBUG i=" << i << " llcount " << llcounts[i] << endl;
    cout << binLabel << " " << "DEBUG i=" << i << " pol " << polscalefactors[i]*(tmcounts[i] + llcounts[i]) << endl;
    cout << binLabel << " " << "DEBUG i=" << i << " onetauscalefactor " << onetauscalefactor << endl;
    cout << binLabel << " " << "DEBUG i=" << i << " onetau " << onetauscalefactor*tmcounts[i] << endl;
    */
  }
  exp_0lep_ttwj += twotauscalefactor*tm2count;
  exp_0lep_ttwj += dilepscalefactor*(tm2count + ll2count);
  
  /*
    cout << binLabel << " " << "DEBUG twotauscalefactor " << twotauscalefactor << endl;
    cout << binLabel << " " << "DEBUG tm2count " << tm2count << endl;
    cout << binLabel << " " << "DEBUG twotau " << twotauscalefactor*tm2count << endl;
    cout << binLabel << " " << "DEBUG ll2count " << ll2count << endl;
    cout << binLabel << " " << "DEBUG dilepscalefactor " << dilepscalefactor << endl;
    cout << binLabel << " " << "DEBUG dilep " << dilepscalefactor*(tm2count + ll2count) << endl;
  */

  //Trigger efficiency
  //cout << binLabel << " ZL truth: " << zlcount << endl;
  //cout << binLabel << " total prediction: " << exp_0lep_ttwj  << endl;
  cout << "bin " << binLabelBen << " (" << binLabel << ") true = " << zlcount << ", prediction = " << exp_0lep_ttwj << endl;
  exp_0lep_ttwj = exp_0lep_ttwj/sltrigeff;
  //cout << "DEBUG total post trig " << exp_0lep_ttwj  << endl;
  
  return exp_0lep_ttwj;
}//get_exp_0lep_ttwj_mr



void metReweightingSanityCheck()
{
  for(int i=0; i<=3; i++)
    { 
      for(int j=0; j<=3; j++)
	{
	  for(int k=0; k<=2; k++)
	    {
	      get_exp_0lep_ttwj_mr(i, j, k, 1);
	    }
	}
    }
}
