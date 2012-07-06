#include <iostream>

void GenerateDummyInputFile() {

  gROOT->Reset();

  int nBinsMET   = 4 ;
  int nBinsHT    = 4 ;
  int nBinsBjets = 3 ;   // this must always be 3

  TString sMbins[nBinsMET];
  TString sHbins[nBinsHT];
  TString sBbins[3] = {"_1b","_2b","_3b"};

  for (int i = 0 ; i < nBinsMET ; i++) {
    TString base = "_M";
    stringstream sbin;
    sbin << i+1;
    base += sbin.str();
    sMbins[i] = base;
  }

  for (int j = 0 ; j < nBinsHT ; j++) {
    TString base = "_H";
    stringstream sbin;
    sbin << j+1;
    base += sbin.str();
    sHbins[j] = base;
  }


  int dummyInt = 99;
  float dummyFloat = 9.999;
  float dummyZero = 0.;
  float dummyOne = 1.;
  float dummyErr = .1;

  ofstream inFile;
  inFile.open("dummy_Input.dat");

  inFile << "Dummy inputs:" << endl;


  // 0lep observables

  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString obs_0lep = "N_0lep" ;
	obs_0lep = obs_0lep+sMbins[i]+sHbins[j]+sBbins[k] ;
	
	inFile << obs_0lep << "  \t" << dummyInt << endl;

      }
    }
  }


  // single lepton observables

  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString obs_1lep = "N_1lep" ;
	obs_1lep = obs_1lep+sMbins[i]+sHbins[j]+sBbins[k] ;
	
	inFile << obs_1lep << "  \t" << dummyInt << endl;

      }
    }
  }


  // ldp observables

  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString obs_ldp = "N_ldp" ;
	obs_ldp = obs_ldp+sMbins[i]+sHbins[j]+sBbins[k] ;
	
	inFile << obs_ldp << "  \t" << dummyInt << endl;

      }
    }
  }


  // R_lsb

  for (int j = 0 ; j < nBinsHT ; j++) {
    for (int k = 0 ; k < nBinsBjets ; k++) {

      TString Rlsb = "R_lsb" ;
      Rlsb = Rlsb+sHbins[j]+sBbins[k] ;
  
      inFile << Rlsb << "  \t" << dummyOne << endl;

      Rlsb = Rlsb+"_err" ;
      inFile << Rlsb << "  \t" << dummyErr << endl;
  
    }
  }

  
  // Z -> ee observables

  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {

      TString obs_Zee = "N_Zee" ;
      obs_Zee = obs_Zee+sMbins[i]+sHbins[j] ;
      
      inFile << obs_Zee << "  \t" << dummyInt << endl;
      
    }
  }

  
  // Z -> mm observables

  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {

      TString obs_Zmm = "N_Zmm" ;
      obs_Zmm = obs_Zmm+sMbins[i]+sHbins[j] ;
      
      inFile << obs_Zmm << "  \t" << dummyInt << endl;
      
    }
  }


  // Nttbarsingletopzjetsmc_ldp

  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString obs_ttbarsingletopzjetsmc_ldp = "N_ttbarsingletopzjetsmc_ldp" ;
	obs_ttbarsingletopzjetsmc_ldp = obs_ttbarsingletopzjetsmc_ldp+sMbins[i]+sHbins[j]+sBbins[k] ;
	
	inFile << obs_ttbarsingletopzjetsmc_ldp << "  \t" << dummyZero << endl;

      }
    }
  }


  // NWJmc_ldp

  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString obs_WJmc_ldp = "N_WJmc_ldp" ;
	obs_WJmc_ldp = obs_WJmc_ldp+sMbins[i]+sHbins[j]+sBbins[k] ;
	
	inFile << obs_WJmc_ldp << "  \t" << dummyZero << endl;

      }
    }
  }


  // NZnnmc_ldp

  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString obs_Znnmc_ldp = "N_Znnmc_ldp" ;
	obs_Znnmc_ldp = obs_Znnmc_ldp+sMbins[i]+sHbins[j]+sBbins[k] ;
	
	inFile << obs_Znnmc_ldp << "  \t" << dummyZero << endl;

      }
    }
  }


  // various parameters needed for Z -> invis.

  // Z -> ee acceptance

  for (int i = 0 ; i < nBinsMET ; i++) {

    TString Zee_acc = "acc_Zee";
    Zee_acc = Zee_acc+sMbins[i];

    inFile << Zee_acc << "  \t" << dummyOne << endl;

    Zee_acc = Zee_acc+"_err";
    inFile << Zee_acc << "  \t" << dummyErr << endl;

  }


  // Z -> mm acceptance

  for (int i = 0 ; i < nBinsMET ; i++) {

    TString Zmm_acc = "acc_Zmm";
    Zmm_acc = Zmm_acc+sMbins[i];

    inFile << Zmm_acc << "  \t" << dummyOne << endl;

    Zmm_acc = Zmm_acc+"_err";
    inFile << Zmm_acc << "  \t" << dummyErr << endl;

  }


  // Z -> ll efficiencies

  inFile << "Z_ee_eff  \t" << dummyOne << endl;
  inFile << "Z_ee_eff_err  \t" << dummyErr << endl;
  inFile << "Z_mm_eff  \t" << dummyOne << endl;
  inFile << "Z_mm_eff_err  \t" << dummyErr << endl;


  // Z -> ll VL to nominal scale factors

  for (int k = 0 ; k < nBinsBjets ; k++) {

    TString knn = "knn" ;
    knn = knn+sBbins[k] ;

    inFile << knn << "  \t" << dummyOne << endl;
    
    knn = knn+"_err" ;
    inFile << knn << "  \t" << dummyErr << endl;

  }


  // Z -> ll purity

  inFile << "Z_ee_pur  \t" << dummyOne << endl;
  inFile << "Z_ee_pur_err  \t" << dummyErr << endl;
  inFile << "Z_mm_pur  \t" << dummyOne << endl;
  inFile << "Z_mm_pur_err  \t" << dummyErr << endl;


  // scale factors:

  inFile << "sf_mc  \t" << dummyOne << endl ; 
  inFile << "sf_mc_err  \t" << dummyErr << endl; 


  // sf_qcd

  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString sf_qcd = "sf_qcd" ;
	sf_qcd = sf_qcd+sMbins[i]+sHbins[j]+sBbins[k] ;
	
	inFile << sf_qcd << "  \t" << dummyOne << endl;	

	sf_qcd = sf_qcd+"_err" ;
	inFile << sf_qcd << "  \t" << dummyErr << endl;	

      }
    }
  }


  // sf_ttwj

  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString sf_ttwj = "sf_ttwj" ;
	sf_ttwj = sf_ttwj+sMbins[i]+sHbins[j]+sBbins[k] ;
	
	inFile << sf_ttwj << "  \t" << dummyOne << endl;	

	sf_ttwj = sf_ttwj+"_err" ;
	inFile << sf_ttwj << "  \t" << dummyErr << endl;	

      }
    }
  }


  // sf_ee

  for (int k = 0 ; k < nBinsBjets ; k++) {

    TString sf_ee = "sf_ee" ;
    sf_ee = sf_ee+sBbins[k] ;

    inFile << sf_ee << "  \t" << dummyOne << endl;
    
    sf_ee = sf_ee+"_err" ;
    inFile << sf_ee << "  \t" << dummyErr << endl;

  }


  // sf_mm

  for (int k = 0 ; k < nBinsBjets ; k++) {

    TString sf_mm = "sf_mm" ;
    sf_mm = sf_mm+sBbins[k] ;

    inFile << sf_mm << "  \t" << dummyOne << endl;
    
    sf_mm = sf_mm+"_err" ;
    inFile << sf_mm << "  \t" << dummyErr << endl;

  }


  inFile.close();
  return;

}
