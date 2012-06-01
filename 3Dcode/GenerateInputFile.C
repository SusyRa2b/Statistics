#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TROOT.h"
#include "TH1F.h"

  using std::stringstream ;
  using std::ofstream ;
  using std::endl ;

// to add in: nMu, nEl, minDelPhi

void GenerateInputFile() {

  TChain* dyTree = new TChain("treeZ") ;
  int nAdded = dyTree->Add("files5fb/DY.root") ;
  if ( nAdded <= 0 ) {
     printf("\n\n\n *** No treeZ in files5fb/DY.root\n\n\n") ;
     return ;
  }

  TChain chain("tree");
  chain.Add("files5fb/QCD-1800.root");
  chain.Add("files5fb/QCD-1400to1800.root");
  chain.Add("files5fb/QCD-1000to1400.root");
  chain.Add("files5fb/QCD-800to1000.root");
  chain.Add("files5fb/QCD-600to800.root");
  chain.Add("files5fb/QCD-470to600.root");
  chain.Add("files5fb/QCD-300to470.root");
  chain.Add("files5fb/QCD-170to300.root");
  chain.Add("files5fb/QCD-120to170.root");
  chain.Add("files5fb/QCD-80to120.root");
  chain.Add("files5fb/QCD-50to80.root");
  chain.Add("files5fb/TT.root");
  chain.Add("files5fb/Zinv.root");

  TChain chainTZ("tree");
  chainTZ.Add("files5fb/TT.root");
  chainTZ.Add("files5fb/Zinv.root");

  gROOT->Reset();

  const int nBinsBjets = 3 ;   // this must always be 3

  //-- met2-ht1-v1
//const int nBinsMET   = 2 ;
//const int nBinsHT    = 1 ;
//float Mbins[nBinsMET+1] = {150.,250.,99999.};
//float Hbins[nBinsHT+1] = {400.,99999.};

  //-- met3-ht2-v1
//const int nBinsMET   = 3 ;
//const int nBinsHT    = 2 ;
//float Mbins[nBinsMET+1] = {150.,250.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,800.,99999.};

  //-- met3-ht3-v1
//const int nBinsMET   = 3 ;
//const int nBinsHT    = 3 ;
//float Mbins[nBinsMET+1] = {150.,250.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,1000.,99999.};

////-- met4-ht3-v1
//const int nBinsMET   = 4 ;
//const int nBinsHT    = 3 ;
//float Mbins[nBinsMET+1] = {150.,250.,350.,500.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,1000.,99999.};

  //-- met5-ht4-v1
//const int nBinsMET   = 5 ;
//const int nBinsHT    = 4 ;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,600.,800.,1000.,99999.};

  //-- met5-ht5-v1
  const int nBinsMET   = 4 ;
  const int nBinsHT    = 4 ;
  float Mbins[nBinsMET+1] = {150.,200.,250.,300.,99999.};
  float Hbins[nBinsHT+1] = {400.,500.,600.,800.,99999.};

  //-- met5-ht5-v1
//const int nBinsMET   = 5 ;
//const int nBinsHT    = 5 ;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,99999.};
//float Hbins[nBinsHT+1] = {400.,500.,600.,800.,1000.,99999.};

  //-- met8-ht8-v1
//const int nBinsMET   = 6 ;
//const int nBinsHT    = 6 ;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,400.,99999.};
//float Hbins[nBinsHT+1] = {400.,500.,600.,700.,800.,900.,99999.};

  //-- met8-ht8-v1
//const int nBinsMET   = 7 ;
//const int nBinsHT    = 7 ;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,400.,500.,99999.};
//float Hbins[nBinsHT+1] = {400.,500.,600.,700.,800.,900.,1000.,99999.};

  //-- met8-ht8-v1
//const int nBinsMET   = 8 ;
//const int nBinsHT    = 8 ;
//float Mbins[nBinsMET+1] = {150.,200.,250.,300.,350.,400.,450.,600.,99999.};
//float Hbins[nBinsHT+1] = {400.,500.,600.,700.,800.,900.,1000.,1200.,99999.};


  TString sMbins[nBinsMET];
  TString sHbins[nBinsHT];
  TString sBbins[3] = {"_1b","_2b","_3b"};

  TString cMbins[nBinsMET];
  TString cHbins[nBinsHT];
//TString cBbins[3] = {"&&nB==1","&&nB==2","&&nB>=3"};

  for (int i = 0 ; i < nBinsMET ; i++) {
    TString base = "_M";
    stringstream sbin;
    sbin << i+1;
    base += sbin.str();
    sMbins[i] = base;
    base = "&&MET>";
    stringstream cbin;
    cbin << Mbins[i] << "&&MET<" << Mbins[i+1];
    base += cbin.str();
    cMbins[i] = base;
  }

  for (int j = 0 ; j < nBinsHT ; j++) {
    TString base = "_H";
    stringstream sbin;
    sbin << j+1;
    base += sbin.str();
    sHbins[j] = base;
    base = "&&HT>";
    stringstream cbin;
    cbin << Hbins[j] << "&&HT<" << Hbins[j+1];
    base += cbin.str();
    cHbins[j] = base;
  }

//int dummyInt = 99;
//float dummyFloat = 9.999;
  float dummyZero = 0.;
  float dummyOne = 1.;
  float dummyErr = .1;

  ofstream inFile;
  char outfile[10000] ;
  sprintf( outfile, "Input-met%d-ht%d.dat", nBinsMET, nBinsHT ) ;
  inFile.open( outfile );

  // print out header line:

  inFile << "Using HT bins:  " ;
  for (int j = 0 ; j <= nBinsHT ; j++ ) {
    inFile << Hbins[j] ;
    if ( j < nBinsHT ) inFile << "-" ;
  }

  inFile << "\t Using MET bins: " ;
  for (int i = 0 ; i <= nBinsMET ; i++ ) {
    inFile << Mbins[i] ;
    if ( i < nBinsMET ) inFile << "-" ;
  }
  
  inFile << endl ;


  // 0lep observables

  printf("\n\n-----------------------------------------------------------------\n\n") ;

  TH1F* ht = new TH1F("ht","ht",10,0,10000);
  TString cuts0lep = "minDelPhiN>4&&nMu==0&&nEl==0&&";
  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString obs_0lep = "N_0lep" ;
	obs_0lep = obs_0lep+sMbins[i]+sHbins[j]+sBbins[k] ;
	
	TString cut = "HT>";
	cut += Hbins[j];
	cut += "&&HT<";
	cut += Hbins[j+1];
	cut += "&&MET>";
	cut += Mbins[i];
	cut += "&&MET<";
	cut += Mbins[i+1];
	cut += "&&nB==";
	cut += k+1;

	chain.Project("ht","HT",cuts0lep+cut);
	//inFile << obs_0lep << "  \t" << chain.GetEntries(cuts0lep+cut) << endl;
        inFile << obs_0lep << "  \t" << (int)ht->GetSumOfWeights() << endl;
        TString allcuts = cuts0lep+cut ;
        printf(" N_0lep -- HT,MET,nbjet bins (%d,%d,%d): npass=%7.1f, cuts=%s\n", j,i,k,ht->GetSumOfWeights(), allcuts.Data()) ;
        ht->Reset() ;

        // signal selection, so MET>250, HT>400, >=1 b, mindelphi>4, 0L, nJets >= 3

      }
    }
  }


  printf("\n\n-----------------------------------------------------------------\n\n") ;

  // single lepton observables
  TString cuts1lep = "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&";
  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString obs_1lep = "N_1lep" ;
	obs_1lep = obs_1lep+sMbins[i]+sHbins[j]+sBbins[k] ;

	TString cut = "HT>";
	cut += Hbins[j];
	cut += "&&HT<";
	cut += Hbins[j+1];
	cut += "&&MET>";
	cut += Mbins[i];
	cut += "&&MET<";
	cut += Mbins[i+1];
	cut += "&&nB==";
	cut += k+1;
	
	chain.Project("ht","HT",cuts1lep+cut);
	//inFile << obs_1lep << "  \t" << chain.GetEntries() << endl;
        inFile << obs_1lep << "  \t" << (int)ht->GetSumOfWeights() << endl;
        TString allcuts = cuts1lep+cut ;
        printf(" N_1lep -- HT,MET,nbjet bins (%d,%d,%d): npass=%7.1f, cuts=%s\n", j,i,k,ht->GetSumOfWeights(),allcuts.Data()) ;
        ht->Reset() ;

        // signal selection but with 1L, so MET>250, HT>400, >=1 b, mindelphi>4, 1L, nJets >= 3

      }
    }
  }

  printf("\n\n-----------------------------------------------------------------\n\n") ;

  // ldp observables
  TString cutsldp = "minDelPhiN<4&&nMu==0&&nEl==0&&";
  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString obs_ldp = "N_ldp" ;
	obs_ldp = obs_ldp+sMbins[i]+sHbins[j]+sBbins[k] ;
	
	TString cut = "HT>";
	cut += Hbins[j];
	cut += "&&HT<";
	cut += Hbins[j+1];
	cut += "&&MET>";
	cut += Mbins[i];
	cut += "&&MET<";
	cut += Mbins[i+1];
	cut += "&&nB==";
	cut += k+1;
	
	chain.Project("ht","HT",cutsldp+cut);
	inFile << obs_ldp << "  \t" << (int)ht->GetSumOfWeights() << endl;
        TString allcuts = cutsldp+cut ;
        printf(" N_ldp -- HT,MET,nbjet bins (%d,%d,%d): npass=%7.1f, cuts=%s\n", j,i,k,ht->GetSumOfWeights(),allcuts.Data()) ;
        ht->Reset() ;

        // signal selection, but ldp, so MET>250, HT>400, >=1 b, mindelphi<4, 0L, nJets >= 3

      }
    }
  }
  
  printf("\n\n-----------------------------------------------------------------\n\n") ;

  //inFile << "Note: I've removed the b jet dependance of this result since it's always with zero b jets" << endl;
  // R_lsb  very low met sideband (50-100) ratio of mdp>4/mdp<4 (with zero b ratio)
  TString cutslsb = "MET>50&&MET<100&&nMu==0&&nEl==0&&nB==0&&";
  /////TString cutslsb = "MET>50&&MET<100&&nMu==0&&nEl==0&&";
  for (int j = 0 ; j < nBinsHT ; j++) {
    for (int k = 0 ; k < nBinsBjets ; k++) {

      TString Rlsb = "R_lsb" ;
      Rlsb = Rlsb+sHbins[j]+sBbins[k] ;
      
      TString cut = "HT>";  
      cut += Hbins[j];
      cut += "&&HT<";
      cut += Hbins[j+1];
      
      TString pass = "&&minDelPhiN>4";
      TString fail = "&&minDelPhiN<4";
      chain.Project("ht","HT",cutslsb+cut+pass);
      float npass = ht->GetSumOfWeights();
      TString allcutspass = cutslsb+cut+pass ;
      printf(" R_lsb -- HT,MET bins (%d,%d): npass=%10.1f, cuts=%s\n", j,k,npass, allcutspass.Data()) ;
        ht->Reset() ;
      chain.Project("ht","HT",cutslsb+cut+fail);
      float nfail = ht->GetSumOfWeights();
      TString allcutsfail = cutslsb+cut+fail ;
      printf(" R_lsb -- HT,MET bins (%d,%d): nfail=%10.1f, cuts=%s\n", j,k,nfail, allcutsfail.Data()) ;
        ht->Reset() ;
      
      inFile << Rlsb << "      \t" << npass/nfail << endl;
      
      float error = TMath::Sqrt( (1/npass) + (1/nfail) )*(npass/nfail);
      Rlsb = Rlsb+"_err" ;
      inFile << Rlsb << "  \t" << error << endl;

    }  
  }

  
  printf("\n\n-----------------------------------------------------------------\n\n") ;

  // Z -> ee observables 

  TString cutszee = "cat==2&&minDelPhiNee>4&&nVLB>=1&&";
  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {

      TString obs_Zee = "N_Zee" ;
      obs_Zee = obs_Zee+sMbins[i]+sHbins[j] ;
      
      TString cut = "HT>";  
      cut += Hbins[j];
      cut += "&&HT<";
      cut += Hbins[j+1];
      cut += "&&METee>";
      cut += Mbins[i];
      cut += "&&METee<";
      cut += Mbins[i+1];

      dyTree->Project("ht","HT",cutszee+cut);
      inFile << obs_Zee << "  \t" << (int)ht->GetSumOfWeights() << endl;
      TString allcuts = cutszee+cut ;
      printf(" N_Zee -- HT,MET bins (%d,%d): events=%7.1f, cuts=%s\n", j,i,ht->GetSumOfWeights(),allcuts.Data() ) ;
        ht->Reset() ;
      //Z->ee counts, with 1 VLb and sig selection, so so MET>250, HT>400, mindelphi>4, 2e, 0mu, nJets >= 3
      
    }
  }

  
  printf("\n\n-----------------------------------------------------------------\n\n") ;

  // Z -> mm observables

  TString cutszmm = "cat==1&&minDelPhiNmm>4&&nVLB>=1&&";
  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {

      TString obs_Zmm = "N_Zmm" ;
      obs_Zmm = obs_Zmm+sMbins[i]+sHbins[j] ;
      
      TString cut = "HT>";  
      cut += Hbins[j];
      cut += "&&HT<";
      cut += Hbins[j+1];
      cut += "&&METmm>";
      cut += Mbins[i];
      cut += "&&METmm<";
      cut += Mbins[i+1];
      
      dyTree->Project("ht","HT",cutszmm+cut);
      inFile << obs_Zmm << "  \t" << (int)ht->GetSumOfWeights() << endl;
      TString allcuts = cutszmm+cut ;
      printf(" N_Zmm -- HT,MET bins (%d,%d): events=%7.1f, cuts=%s\n", j,i,ht->GetSumOfWeights(),allcuts.Data() ) ;
        ht->Reset() ;
      //Z->mm counts, with 1 VLb and sig selection, so so MET>250, HT>400, mindelphi>4, 2mu, 0e, nJets >= 3

    }
  }

  //inFile << "Why are all three of these MC categories separate? I've just put them together (only QCD, Zinv, tt)" << endl;
  // Nttbarsingletopzjetsmc_ldp


  printf("\n\n-----------------------------------------------------------------\n\n") ;

  for (int i = 0 ; i < nBinsMET ; i++) {
    for (int j = 0 ; j < nBinsHT ; j++) {
      for (int k = 0 ; k < nBinsBjets ; k++) {

	TString obs_ttbarsingletopzjetsmc_ldp = "N_ttbarsingletopzjetsmc_ldp" ;
	obs_ttbarsingletopzjetsmc_ldp = obs_ttbarsingletopzjetsmc_ldp+sMbins[i]+sHbins[j]+sBbins[k] ;
	
	TString cut = "HT>";
	cut += Hbins[j];
	cut += "&&HT<";
	cut += Hbins[j+1];
	cut += "&&MET>";
	cut += Mbins[i];
	cut += "&&MET<";
	cut += Mbins[i+1];
	cut += "&&nB==";
	cut += k+1;

      chainTZ.Project("ht","HT",cutsldp+cut);
	inFile << obs_ttbarsingletopzjetsmc_ldp << "  \t" << ht->GetSumOfWeights() << endl;
        ht->Reset() ;
// signal selection, but ldp, so MET>250, HT>400, >=1 b, mindelphi<4, 0L, nJets >= 3
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
// signal selection, but ldp, so MET>250, HT>400, >=1 b, mindelphi<4, 0L, nJets >= 3
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
// signal selection, but ldp, so MET>250, HT>400, >=1 b, mindelphi<4, 0L, nJets >= 3
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

// use 2011 values for now.
  inFile << "Z_ee_eff  \t" << 0.6774 << endl;
  inFile << "Z_ee_eff_err  \t" << 0.0580 << endl;
  inFile << "Z_mm_eff  \t" << 0.7217 << endl;
  inFile << "Z_mm_eff_err  \t" << 0.0506 << endl;


  // Z -> ee VL to nominal scale factors
// would it make more sense to count events in only the loosest HT and/or MET bin and
// use the scale factors to translate between the different HT and MET bins?

/*  for (int k = 0 ; k < nBinsBjets ; k++) {
    TString knn_ee = "knn_ee" ;
    knn_ee = knn_ee+sBbins[k] ;
    inFile << knn_ee << "  \t" << dummyOne << endl;
    knn_ee = knn_ee+"_err" ;
    inFile << knn_ee << "  \t" << dummyErr << endl;
  } */
  

  // use 2011 values for now.
  inFile << "knn_1b     \t" << 0.401 << endl;
  inFile << "knn_1b_err \t" << 0.018 << endl;
  inFile << "knn_2b     \t" << 0.067 << endl;
  inFile << "knn_2b_err \t" << 0.009 << endl;
  inFile << "knn_3b     \t" << 0.009 << endl;
  inFile << "knn_3b_err \t" << 0.003 << endl;


  // Z -> ll purity

  // use 2011 values for now.
  inFile << "Z_ee_pur  \t" << 0.911 << endl;
  inFile << "Z_ee_pur_err  \t" << 0.104 << endl;
  inFile << "Z_mm_pur  \t" << 0.866 << endl;
  inFile << "Z_mm_pur_err  \t" << 0.079 << endl;


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
