#include <iostream>

void GenerateSusyFile() {

//gluino cross sections in pb.

  TFile prospino("referenceXSecs.root");
  TH1F *gluinoxsec = gluino;

// bin 1 = gluino mass = 100 GeV, 2 = 125, 3 = 150, 4 = 175, ...
// so gluino mass = 75+nbin*25; or nbin = (gluinomass-75)/25.

  TChain chainT1bbbb("tree");
  chainT1bbbb.Add("files5fb/T1bbbb.root");

  gROOT->Reset();

  const int nBinsMET   = 2 ;
  const int nBinsHT    = 2 ;
  const int nBinsBjets = 3 ;   // this must always be 3

  float Mbins[nBinsMET+1] = {150.,250.,99999.};
  float Hbins[nBinsHT+1] = {400.,600.,99999.};

  TString sMbins[nBinsMET];
  TString sHbins[nBinsHT];
  TString sBbins[3] = {"_1b","_2b","_3b"};

  TString cMbins[nBinsMET];
  TString cHbins[nBinsHT];
  TString cBbins[3] = {"&&nB==1","&&nB==2","&&nB>=3"};


  // dummy masses
  int minGlMass = 900 ;
  int maxGlMass = 910 ;

  double dummyYield = 9.9 ;
  double dummyCorr = 1. ;
  double dummyErr = 10. ; // error is in % 
  long dummyEvts = 10000 ;
  
  ofstream inFile;
  inFile.open("Susy.dat");

  // loop over gluino masses

  TH1F* ht = new TH1F("ht","ht",10,0,10000);
  TString cutsSig = "minDelPhiN>4&&nMu==0&&nEl==0&&";
  TString cutsSB = "minDelPhiN>4&&(nMu==1||nEl==1)&&";
  TString cutsLSB = "minDelPhiN>4&&nMu==0&&nEl==0&&MET>50&&MET<100&&";

  float xsec = -1.;
  for ( int mGl = minGlMass ; mGl < maxGlMass ; mGl = mGl + 25 ) {
  xsec = gluinoxsec->GetBinContent((mGl-75)/25);
//    for ( int mLsp = 50 ; mLsp < ( mGl - 25 ) ; mLsp = mLsp + 25 ) {
    for ( int mLsp = 300 ; mLsp < 310 ; mLsp = mLsp + 25 ) {

      inFile << mGl << " " << mLsp << " " << dummyEvts << " " ;

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBjets ; k++) {
// the three entries are the signal region, 1L region and very low met sideband
	
            TString cutSMS = "mgluino>";
	    cutSMS += mGl-1;
	    cutSMS += "&&mgluino<";
	    cutSMS += mGl+1;
	    cutSMS += "&&mlsp>";
	    cutSMS += mLsp-1;
	    cutSMS += "&&mlsp<";
	    cutSMS += mLsp+1;
	    cutSMS += "&&";	    
	
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

	    TString cutNoMET = "HT>";
	    cutNoMET += Hbins[j];
	    cutNoMET += "&&HT<";
	    cutNoMET += Hbins[j+1];
	    cutNoMET += "&&nB==";
	    cutNoMET += k+1;

// cross section in pb = gluinoxsec->GetBinContent((mGl-75)/25);
//cout << " bin = " << (mGl-75)/25 << " for mGl = " << mGl << " and xsec = " << gluinoxsec->GetBinContent((mGl-75)/25) << endl;

// each point has 10k events generated. The sigma is in pb and I want to normalized to 5 fb-1. 
// so multiple cross section by 0.5 to get events in 5 fb-1

	    chainT1bbbb.Project("ht","HT",cutSMS+cutsSig+cut);
            inFile << 0.5*xsec*ht->GetSumOfWeights() << " ";
            ht->Reset() ;
	    chainT1bbbb.Project("ht","HT",cutSMS+cutsSB+cut);
            inFile << 0.5*xsec*ht->GetSumOfWeights() << " ";
            ht->Reset() ;
	    chainT1bbbb.Project("ht","HT",cutSMS+cutsLSB+cutNoMET);
            inFile << 0.5*xsec*ht->GetSumOfWeights() << " ";
            ht->Reset() ;

	  }
	}
      }

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBjets ; k++) {
	    inFile << dummyCorr << " " << dummyCorr << " " << dummyCorr << " " ;
	  }
	}
      }

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBjets ; k++) {
	    inFile << dummyErr << " " << dummyErr << " " << dummyErr << " " ;
	  }
	}
      }

      inFile << endl ;

    }
  }

  
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


  return;

}
