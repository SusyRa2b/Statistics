#include <iostream>
#include <fstream>

#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TFile.h"

using namespace std;

double tZeroLeptonTopWJetsYieldTotal=0;
double tZeroLeptonQCDYieldTotal=0;
double tZeroLeptonZtoNuNuYieldTotal=0;

int trueN = 0; 
TString trueXsec_string = "0.0";
double trueXsec=0;


TTree* treeToys = new TTree("treeToys","treeToys");

bool drawTrue = true;




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







void makeTree(int trueN_in, TString inFile, TString fitter = "LB") {

  
  TString branchDescriptor; 
  if(fitter == "OAK") {
    branchDescriptor = "zeroLeptonSignalYieldTotal:zeroLeptonTopWJetsYieldTotal:zeroLeptonQCDYieldTotal:zeroLeptonZtoNuNuYieldTotal";
  }
  else if(fitter == "test") {
    //branchDescriptor = "zeroLeptonCountTotal";
    branchDescriptor = "XXsignalCrossSection:signalCrossSectionError:signalCrossSection";
  }
  else if(fitter == "sensitivity") {
    branchDescriptor = "signalCrossSection/D:signalCrossSectionError";
    branchDescriptor += ":lowerLimit95:upperLimit95";
    branchDescriptor += ":susyFloatingNLL:susyFixedToZeroNLL";
  }
  else if(fitter == "bins") {
    branchDescriptor = "bin1/D:bin2:bin3:bin4:bin5:bin6:bin7:bin8:bin9:bin10:bin11:bin12:bin13:bin14:bin15:bin16";
    branchDescriptor += ":bin17:bin18:bin19:bin20:bin21:bin22:bin23:bin24:bin25:bin26:bin27:bin28:bin29:bin30:bin31:bin32";
    branchDescriptor += ":bin33:bin34:bin35:bin36:bin40:bin41:bin42:bin43:bin44:bin45:bin46:bin47:bin48";
  }
  else {
    branchDescriptor = "signalCrossSection/D:signalCrossSectionError";
    branchDescriptor += ":lowerLimit68:upperLimit68:lowerLimit95:upperLimit95";
    branchDescriptor += ":signalUncertainty";
    branchDescriptor += ":zeroLeptonSignalYieldTotal:zeroLeptonSignalYieldTotalError";
    branchDescriptor += ":zeroLeptonTopWJetsYieldTotal:zeroLeptonQCDYieldTotal:zeroLeptonZtoNuNuYieldTotal:zeroLeptonDibosonYieldTotal:zeroLeptonCountTotal";
  }

  treeToys->ReadFile(inFile,branchDescriptor);
  
  
  //TString foutName = ""; foutName += "tree_"; foutName += trueN; foutName += "_"; foutName += fitter; foutName += ".root";
  //TFile fout(foutName, "RECREATE");
  //treeToys->Write();
  //fout.Close();
  

}


void readNominal(map<TString,float> &genMap)
{

  fstream fin("saveGenerate.dat", ios::in);
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


void binTruth(TString componentName = "QCD") {
  
  TString trueN_string = ""; trueN_string+=trueN;
  TFile f1("fbins_"+trueN_string+".root", "UPDATE");
  f1.cd();

  map<TString, float> genMap;
  readNominal(genMap);
  
  TH1D h1("h1", componentName+" offset yield = fitted yield - true yield", 48, 0.5, 48.5);
  TH1D h1_percent("h1_percent", componentName+" percent difference from truth = 100*(fitted yield - true yield)/(true yield)", 48, 0.5, 48.5);
  //TH1D hpercentdiff("hpercentdiff", componentName+" percent difference = 100*(y1-y2)/<y>", 48, 0.5, 48.5);
  //TH1D hsignificance("hsignificance", componentName+" difference significance = (y1-y2)/sqrt(e1^2+e2^2)", 48, 0.5, 48.5);
  
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
	      
	      TString translatedName = "";
	      //if(componentName == "Signal")   translatedName = "susy";
	      if(componentName == "TopWJets")     translatedName = "tt";
	      else if(componentName == "QCD")      translatedName = "qcd";
	      else if(componentName == "ZtoNuNu")  translatedName = "znn";
	      else if(componentName == "Diboson")  translatedName = "vv";
	      else{ assert(0); }

	      TString mapName = "N_0lep_";
	      mapName += translatedName;
	      mapName += "_";
	      mapName += bin;
	      
	      //get true value
	      double trueValue = genMap[mapName];
	      cout << mapName << " " << trueValue << endl;
	      if(componentName == "TopWJets")
		{
		  mapName = "N_0lep_wjets_";
		  mapName += bin;
		  trueValue += genMap[mapName];
		  cout << mapName << " add " << genMap[mapName] << " to get " << trueValue << endl;
		}
	      
	      //get fitted value
	      TH1D htemp("htemp", "htemp", 1, 0, 1e9);
	      treeToys->Project("htemp",binLB);
	      double value1 = htemp.GetMean();
	      double error1 = htemp.GetRMS();
	      cout << "fit: " << value1 << " +- " << error1 << endl;
	      
	      
	      h1.SetBinContent(binNum, value1-trueValue);
	      h1.SetBinError(binNum, error1);
	      
	      if( trueValue>0 ) h1_percent.SetBinContent(binNum, 100.0*(value1-trueValue)/(trueValue));
							   
	      h1.GetXaxis()->SetBinLabel(binNum, bin);
	      h1_percent.GetXaxis()->SetBinLabel(binNum, bin);
	      //hsignificance.GetXaxis()->SetBinLabel(binNum, bin);
	      
	    }//b
	}//h
    }//m
  
  
  h1.SetMarkerStyle(20);
  h1_percent.SetMarkerStyle(20);
  //hsignificance.SetMarkerStyle(20);
  
  h1.SetMarkerSize(1);
  h1_percent.SetMarkerSize(1);
  //hsignificance.SetMarkerSize(1);
  
  h1.GetXaxis()->LabelsOption("v");
  h1_percent.GetXaxis()->LabelsOption("v");
  //hsignificance.GetXaxis()->LabelsOption("v");
    
  h1.Write();
  h1_percent.Write();
  f1.Close();

}


void drawTest() {
  treeToys->Draw("zeroLeptonTopWJetsYieldTotal");
}


void makeHists(TString fitter = "LB") {
  
  gStyle->SetStatY(0.98);
  gStyle->SetStatX(0.98);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);

  /*
  TH1D* hSignalCrossSection = new TH1D("hSignalCrossSection", "Signal cross section", 50, 0, 200);
  TH1D* hPull = new TH1D("hPull", "Pull using fit errors", 50, -5, 5);
  TH1D* hPullPL = new TH1D("hPullPL", "Pull using PL errors", 50, -5, 5);
  TH1D* hZeroLeptonSignalYieldTotal = new TH1D("hZeroLeptonSignalYieldTotal", "Zero lepton signal yield", 50, 0, 200);
  TH1D* hZeroLeptonTopWJetsYieldTotal = new TH1D("hZeroLeptonTopWJetsYieldTotal", "Zero lepton top W-jets yield", 50, 1500, 2200);
  TH1D* hZeroLeptonQCDYieldTotal = new TH1D("hZeroLeptonQCDYieldTotal", "Zero lepton QCD yield", 50, 0, 500);
  TH1D* hZeroLeptonZtoNuNuYieldTotal = new TH1D("hZeroLeptonZtoNuNuYieldTotal", "Zero lepton Z-invisible yield", 50, 150, 450);
  */
  TH1D* hSignalCrossSection = new TH1D("hSignalCrossSection", "Signal cross section", 50, 0, 100);
  TH1D* hPull = new TH1D("hPull", "Pull using fit errors", 50, -5, 5);
  TH1D* hPullPL = new TH1D("hPullPL", "Pull using PL errors", 50, -5, 5);
  TH1D* hZeroLeptonSignalYieldTotal = new TH1D("hZeroLeptonSignalYieldTotal", "Zero lepton signal yield", 50, 0, 500);
  TH1D* hZeroLeptonTopWJetsYieldTotal = new TH1D("hZeroLeptonTopWJetsYieldTotal", "Zero lepton top W-jets yield", 50, 9000, 16000);
  TH1D* hZeroLeptonQCDYieldTotal = new TH1D("hZeroLeptonQCDYieldTotal", "Zero lepton QCD yield", 50, 15000, 22000);
  TH1D* hZeroLeptonZtoNuNuYieldTotal = new TH1D("hZeroLeptonZtoNuNuYieldTotal", "Zero lepton Z-invisible yield", 50, 1000, 2000);
  TH1D* hZeroLeptonYieldTotal = new TH1D("hZeroLeptonYieldTotal", "Zero lepton yield", 50, 25000, 28000);
  TH1D* hZeroLeptonCountTotal = new TH1D("hZeroLeptonCountTotal", "Zero lepton count", 50, 25000, 28000);
  TH1D* hUpperLimit95 = new TH1D("hUpperLimit95", "95% Confidence Level Upper Limit", 50, 0, 200);
  TH1D* hSignificance = new TH1D("hSignificance", "Significance", 50, 0, 5);


  hSignalCrossSection->SetFillColor(6);
  hZeroLeptonSignalYieldTotal->SetFillColor(6);  
  hZeroLeptonTopWJetsYieldTotal->SetFillColor(kBlue-9);
  hZeroLeptonQCDYieldTotal->SetFillColor(2);
  hZeroLeptonZtoNuNuYieldTotal->SetFillColor(kGreen-3);
  hZeroLeptonYieldTotal->SetFillColor(kGray+1);
  hZeroLeptonCountTotal->SetFillColor(kGray+1);
  hUpperLimit95->SetFillColor(kGray+1);
  hSignificance->SetFillColor(kGray+1);

  treeToys->Project("hSignalCrossSection","signalCrossSection");
  treeToys->Project("hZeroLeptonSignalYieldTotal", "zeroLeptonSignalYieldTotal");
  treeToys->Project("hZeroLeptonTopWJetsYieldTotal", "zeroLeptonTopWJetsYieldTotal");
  treeToys->Project("hZeroLeptonQCDYieldTotal", "zeroLeptonQCDYieldTotal");
  treeToys->Project("hZeroLeptonZtoNuNuYieldTotal", "zeroLeptonZtoNuNuYieldTotal");
  treeToys->Project("hZeroLeptonYieldTotal", "zeroLeptonSignalYieldTotal+zeroLeptonTopWJetsYieldTotal+zeroLeptonQCDYieldTotal+zeroLeptonZtoNuNuYieldTotal+zeroLeptonDibosonYieldTotal");
  treeToys->Project("hZeroLeptonCountTotal", "zeroLeptonCountTotal");

  treeToys->Project("hPullPL", "(signalCrossSection-"+trueXsec_string+")/( (upperLimit68-lowerLimit68)/2.0 )");
  treeToys->Project("hPull", "(signalCrossSection-"+trueXsec_string+")/( signalCrossSectionError )");

  treeToys->Project("hUpperLimit95", "upperLimit95");
  treeToys->Project("hSignificance", "sqrt(2.0*(susyFixedToZeroNLL - susyFloatingNLL))");

  cout << "CrossSection = " << hSignalCrossSection->GetMean() << " +- " << hSignalCrossSection->GetRMS() << endl;
  cout << "Yield Total = " << hZeroLeptonYieldTotal->GetMean() << " +- " << hZeroLeptonYieldTotal->GetRMS() << endl;
  cout << "Count Total = " << hZeroLeptonCountTotal->GetMean() << " +- " << hZeroLeptonCountTotal->GetRMS() << endl;

  //outfile
  TString foutName = ""; foutName += "output_"; foutName += trueN; foutName += "_"; foutName += fitter; foutName += ".dat";
  ofstream fout;
  fout.open(foutName.Data(), ios::out | ios::trunc);
  assert(fout.is_open());
  fout << trueN << " " << hZeroLeptonSignalYieldTotal->GetMean() << " " << hZeroLeptonSignalYieldTotal->GetRMS() << " ";
  fout << hZeroLeptonTopWJetsYieldTotal->GetMean() << " " << hZeroLeptonTopWJetsYieldTotal->GetRMS() << " " ;
  fout << hZeroLeptonQCDYieldTotal->GetMean() << " " << hZeroLeptonQCDYieldTotal->GetRMS() << " " ;
  fout << hZeroLeptonZtoNuNuYieldTotal->GetMean() << " " << hZeroLeptonZtoNuNuYieldTotal->GetRMS() << endl;
  fout.close();


  //Plots
  TString trueN_string = ""; trueN_string+=trueN;
  TCanvas* cBackground = new TCanvas("cBackground", "Backgrounds", 1200, 330);
  cBackground->Divide(3,1);
  cBackground->cd(1);
  hZeroLeptonTopWJetsYieldTotal->Draw();
  TLine * lttwj = new TLine(tZeroLeptonTopWJetsYieldTotal, 0, tZeroLeptonTopWJetsYieldTotal, hZeroLeptonTopWJetsYieldTotal->GetMaximum());
  lttwj->SetLineWidth(3);
  lttwj->SetLineColor(kOrange-3);
  lttwj->SetLineStyle(2);
  if(drawTrue) lttwj->Draw();

  cBackground->cd(2);
  hZeroLeptonQCDYieldTotal->Draw();
  TLine * lqcd = new TLine(tZeroLeptonQCDYieldTotal, 0, tZeroLeptonQCDYieldTotal, hZeroLeptonQCDYieldTotal->GetMaximum());
  lqcd->SetLineWidth(3);
  lqcd->SetLineColor(kOrange-3);
  lqcd->SetLineStyle(2);
  if(drawTrue) lqcd->Draw();

  cBackground->cd(3);
  hZeroLeptonZtoNuNuYieldTotal->Draw();
  TLine * lznn = new TLine(tZeroLeptonZtoNuNuYieldTotal, 0, tZeroLeptonZtoNuNuYieldTotal, hZeroLeptonZtoNuNuYieldTotal->GetMaximum());
  lznn->SetLineWidth(3);
  lznn->SetLineColor(kOrange-3);
  lznn->SetLineStyle(2);
  if(drawTrue) lznn->Draw();

  cBackground->Print("cBackground_"+trueN_string+"_"+fitter+".pdf");

  TCanvas* cSignalCrossSection = new TCanvas("cSignalCrossSection", "Signal cross section", 800, 330);
  cSignalCrossSection->Divide(2,1);
  cSignalCrossSection->cd(1);
  hSignalCrossSection->Draw();
  TLine * lsig = new TLine(trueXsec, 0, trueXsec, hSignalCrossSection->GetMaximum());
  lsig->SetLineWidth(3);
  lsig->SetLineColor(kOrange-3);
  lsig->SetLineStyle(2);
  if(drawTrue) lsig->Draw();

  cSignalCrossSection->cd(2);
  //hPullPL->Draw();
  hPull->Draw();
  cSignalCrossSection->Print("cSignalCrossSection_"+trueN_string+"_"+fitter+".pdf");

  TCanvas * cZeroLeptonSignalYieldTotal = new TCanvas("cZeroLeptonSignalYieldTotal", "Zero lepton signal yield", 400, 330);
  cZeroLeptonSignalYieldTotal->cd();
  hZeroLeptonSignalYieldTotal->Draw();
  TLine * lsigy = new TLine(trueN, 0, trueN, hZeroLeptonSignalYieldTotal->GetMaximum());
  lsigy->SetLineWidth(3);
  lsigy->SetLineColor(kOrange-3);
  lsigy->SetLineStyle(2);
  if(drawTrue) lsigy->Draw();

  cZeroLeptonSignalYieldTotal->Print("cZeroLeptonSignalYieldTotal_"+trueN_string+"_"+fitter+".pdf");
  
  TCanvas* cZeroLeptonYieldTotal = new TCanvas("cZeroLeptonYieldTotal", "Zero lepton yield", 2*400, 330);
  cZeroLeptonYieldTotal->Divide(2,1);
  cZeroLeptonYieldTotal->cd(1);
  hZeroLeptonYieldTotal->Draw();
  cZeroLeptonYieldTotal->cd(2);
  hZeroLeptonCountTotal->Draw();
  cZeroLeptonYieldTotal->Print("cZeroLeptonYield_"+trueN_string+"_"+fitter+".pdf");

  TCanvas* cSensitivity = new TCanvas("cSensitivity", "Sensitivity", 800, 330);
  cSensitivity->Divide(2,1);
  cSensitivity->cd(1);
  hUpperLimit95->Draw();
  cSensitivity->cd(2);
  hSignificance->Draw();
  cSensitivity->Print("cSensitivity_"+trueN_string+"_"+fitter+".pdf");


  return;
}


void initialize(int trueN_in, TString inFile, TString fitter = "LB") {

  trueN = trueN_in;
  makeTree(trueN, inFile, fitter);

  double trueN_double = trueN;

  //last number should not included skipped bins

  //signal point t1bbbb 850 600, lumi = 19.399/fb
  trueXsec = 1000.0 * trueN_double / 1.9399 / 2651.34 ; 
  
  //signal point t1tttt 1175 400, lumi = 19.399/fb
  //trueXsec = trueN_double * 49996.0  / (19.399 * 13496.559610000000248) ;

  cout << "true xsec = " << trueXsec << endl;

  trueXsec_string = "";
  trueXsec_string += trueXsec;

  //true ttbar, qcd, znn, lumi = 15/fb, useExpected0lep
  tZeroLeptonTopWJetsYieldTotal = 18825.470703;
  tZeroLeptonQCDYieldTotal = 18169.785156;
  tZeroLeptonZtoNuNuYieldTotal = 1451.139893;

}

